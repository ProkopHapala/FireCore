

// No need to explicitely include the OpenCL headers 
//#include <clFFT.h>
#include "libOCLfft.h"

#include  "OCL.h"

#include  "approximation.h"

#include "Grid.h"

#include "FireCoreAPI.h"


bool loadWf_(const char* fname, float* out){
    const int nbuff = 1024;
    char buff[nbuff]; char* line;
    printf( "loadWf %s \n", fname );
    FILE *pFile = fopen(fname, "r" );
    if(pFile==0)return false;
    //printf( "loadWf_ 0 \n" );
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    //printf( "loadWf_ 1 \n" );
    //double xs[4];
    int n=0;
    while(true){
        line=fgets(buff,nbuff,pFile);
        //printf( "loadWf_ >>%s<< \n", line );
        for(int i=0; i<nbuff; i++){ if(line[i]=='D')line[i]='e'; }
        //int i = sscanf (line, "%lf %lf %lf %lf\n", &out[0], &out[1], &out[2], &out[3] );
        int i = sscanf (line, "%f %f %f %f\n", &out[0], &out[1], &out[2], &out[3] );
        if(i!=4) break;
        //printf( " %g %g %g %g \n", out[0], out[1], out[2], out[3] );
        out+=4;
        n+=4;
    }
    fclose(pFile);
    out-=n;
    //for(int i=0; i<n; i++){ printf( "DEBUG[%i] %g \n", i, out[i] ); }
    return true;
}

void resample1D(int nout, float x0to, float x0from, float dxTo, float dxFrom, float* from, float* to, int pitch, int off ){
    float invdx = 1/dxFrom;
    for(int i=0; i<nout; i++){
        float  x = ((x0to + i*dxTo)-x0from)*invdx;
        int   ix = (int)x;
        float  f = x-ix;
        int iout = i*pitch+off;
        to[iout] =  from[ix]*(1-f)-from[ix+1]*f;  // lerp
        //printf("resample1D %i : %g \n", i, to[iout] );
    }
}


OCLsystem ocl;
Approx::AutoApprox aaprox;
FireCore::Lib fireCore;



typedef struct{
    cl_uint  dim;
    size_t global[3];
    size_t local[3];
    //cl_int blocksize;
} KernelDims;

class OCLfft : public OCLsystem { public:


    clfftPlanHandle planHandle;
    clfftDim fft_dim = CLFFT_3D;

    //const size_t N0 = 4, N1 = 4, N2 = 4;
    int ndim=0;
    size_t Ns[4]; // = {N0, N1, N2};
    size_t Ntot;
    int4   Nvec;
    //size_t buffer_size;

    //static cl_mem data_cl;
    //static float *data;

    //int N=0;
    //int size = N * N;
    //static cl_mem d_a, d_b, d_c;

    int iKernell_mull=-1;
    int iKernell_project=-1;
    int iKernell_project_tex=-1;
    int iKernell_projectPos_tex=1;
    OCLtask cltask_mul;
    OCLtask cltask_project;
    OCLtask cltask_project_tex;
    OCLtask cltask_projectPos_tex;
    int itex_basis=-1;

    int    nAtoms;
    int    nPos;
    float4 pos0, dA, dB, dC;
    GridShape grid;
    int ibuff_atoms,ibuff_coefs;
    //int ibuff_poss;

    void updateNtot(){
        Ntot=1; for(int i=0; i<ndim; i++){ Ntot=Ns[i]; };
    }

    void planFFT(){
        // Create a default plan for a complex FFT. 
        err = clfftCreateDefaultPlan(&planHandle, context, fft_dim, Ns );
        err = clfftSetPlanPrecision (planHandle, CLFFT_SINGLE);
        err = clfftSetLayout        (planHandle, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED);
        err = clfftSetResultLocation(planHandle, CLFFT_INPLACE);
        err = clfftBakePlan         (planHandle, 1, &commands, NULL, NULL);
    }

    void initFFT( int ndim, size_t* Ns_ ){
        //printf("DEBUG initFFT() ndim %i Ns[%li,%li,%li]\n", ndim, Ns[0], Ns[1], Ns[2] );
        if     (ndim==1){ fft_dim=CLFFT_1D; }
        else if(ndim==2){ fft_dim=CLFFT_2D; }
        else if(ndim==3){ fft_dim=CLFFT_3D; };
        //buffer_size  = sizeof(float2);
        Ntot=1; for(int i=0; i<ndim;i++){  Ns[i]=Ns_[i];  Ntot*= Ns[i]; }
        printf( "initFFT ndim %i Ntot %li [%li,%li,%li]\n", ndim, Ntot, Ns[0],Ns[1],Ns[2] );
        clfftSetupData fftSetup;                printf("initFFT 1 \n");
        err  = clfftInitSetupData(&fftSetup);   printf("initFFT 2 \n");
        err  = clfftSetup        (&fftSetup);   printf("initFFT 3 \n");
        //data_cl = clCreateBuffer( context, CL_MEM_READ_WRITE, buffer_size, NULL, &err );
        planFFT(  );                            printf("initFFT 4 \n");
    }

    int newFFTbuffer( char* name ){
        return newBuffer( name, Ntot, sizeof(float2), 0, CL_MEM_READ_WRITE );
    }

    int initAtoms( int nAtoms_ ){
        nAtoms=nAtoms_;
        ibuff_atoms=newBuffer( "atoms", nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY );
        ibuff_coefs=newBuffer( "coefs", nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY );
        return ibuff_atoms;
    };

    int initBasisTable( int nx, int ny, float* data ){
        printf( "initBasisTable %i %i \n", nx, ny );
        itex_basis = newBufferImage2D( "BasisTable", ny, nx,   sizeof(float)*4,  data, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR , {CL_RGBA, CL_FLOAT} );
        //itex_basis = newBufferImage2D( "BasisTable", nx, ny,   sizeof(float),  data, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR , {CL_R, CL_FLOAT} ); // THIS WORKS FOR FLOAT TEXTURE
        //itex_basis = newBufferImage2D( "BasisTable", nx/4, ny,   sizeof(float)*4,  data, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR , {CL_RGBA, CL_FLOAT} );
        return itex_basis;
    }

    //void loadData( float* data_ ){
        //printf("DEBUG loadData() \n");
        //data = data_;
        //makeData();
        //printData( data ); 
        //err = clEnqueueWriteBuffer  ( commands, data_cl, CL_TRUE, 0, buffer_size, data, 0, NULL, NULL );
    //}

    void runFFT( int ibuff, bool fwd, float* data=0 ){
        //err = clEnqueueWriteBuffer ( queue, data_cl, CL_TRUE, 0, buffer_size, data, 0, NULL, NULL );
        cl_mem data_cl =  buffers[ibuff].p_gpu;
        if(fwd){
            err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &data_cl, NULL, NULL);  // Execute the plan. -> Forward Transform
        }else{
            err = clfftEnqueueTransform( planHandle, CLFFT_BACKWARD, 1, &commands, 0, NULL, NULL, &data_cl, NULL, NULL);  // Execute the plan.  -> Backward Transform      
        }                                                  
        if(data){
            err = clEnqueueReadBuffer  ( commands, data_cl, CL_TRUE, 0, buffers[ibuff].byteSize(), data, 0, NULL, NULL ); // Fetch results of calculations. 
            err = clFinish(commands);    // Wait for calculations to be finished. 
        }
        //printData( data ); 
    }

    void mul_buffs( int ibuffA, int ibuffB, int ibuff_result ){
        KernelDims kdim;
        kdim.dim        = 1;
        kdim.global[0]  = Ntot*2;
        kdim.local [0]  = 16;
        cl_kernel kernel = kernels[iKernell_mull]; 
        //int N2 = Ntot*2;
        err =  clSetKernelArg(kernel, 0, sizeof(int),    &kdim.global        );
        err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &(buffers[0].p_gpu) );
        err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &(buffers[1].p_gpu) );
        err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &(buffers[2].p_gpu) );
        //checkError(err, "Setting kernel args");
        //double start_time = wtime();
        printf( "mul_buffs kdim: dim %i global %li local %li \n", kdim.dim, kdim.global[0], kdim.local[0] ); 
        err = clEnqueueNDRangeKernel(  commands, kernel,   kdim.dim, NULL, kdim.global, kdim.local, 0, NULL, NULL);    
        OCL_checkError(err, "Enqueueing kernel");
        //err = clFinish(commands);
        //OCL_checkError(err, "Waiting for kernel to finish");
        //double run_time = wtime() - start_time;
    }


    void projectAtomPosTex( int ibuffAtoms, int ibuffCoefs, int nPos, float* poss, float* out ){
        KernelDims kdim;
        kdim.dim        = 1;
        kdim.global[0]  = Ntot*2;
        kdim.local [0]  = 16;
        cl_kernel kernel = kernels[iKernell_projectPos_tex]; 
        /*
            INTarg (nAtoms),        //1
            BUFFarg(ibuffAtoms),    //2
            BUFFarg(ibuffCoefs),    //3
            INTarg (nPos),          //4
            BUFFarg(ibuff_poss),    //5
            BUFFarg(ibuff_out),     //6
            BUFFarg(itex_basis),    //7
        */
        int ibuff_poss = newBuffer( "poss", nPos, sizeof(float4), poss, CL_MEM_READ_WRITE );
        int ibuff_out  = newBuffer( "out" , nPos, sizeof(float2), out , CL_MEM_READ_WRITE );
        buffers[ibuff_poss].toGPU(commands);
        //buffers[ibuff_out ].toGPU(commands);
        //err = clEnqueueWriteBuffer ( commands, buffers[ibuff_poss].p_gpu, CL_TRUE, 0, sizeof(float4)*nPos, poss, 0, NULL, NULL );
        //err = clEnqueueWriteBuffer ( commands, buffers[ibuff_out ].p_gpu, CL_TRUE, 0, sizeof(float2)*nPos, out,  0, NULL, NULL );
        err =  clSetKernelArg(kernel, 0, sizeof(int),    &nAtoms        );
        err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &(buffers[ibuffAtoms].p_gpu) );
        err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &(buffers[ibuffCoefs].p_gpu) );
        err =  clSetKernelArg(kernel, 3, sizeof(int),    &nPos          );
        err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &(buffers[ibuff_poss].p_gpu) );
        err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &(buffers[ibuff_out].p_gpu) );
        err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), &(buffers[itex_basis].p_gpu) );
        //checkError(err, "Setting kernel args");
        //double start_time = wtime();
        printf( "mul_buffs kdim: dim %i global %li local %li \n", kdim.dim, kdim.global[0], kdim.local[0] ); 
        err = clEnqueueNDRangeKernel(  commands, kernel,   kdim.dim, NULL, kdim.global, kdim.local, 0, NULL, NULL);    
        OCL_checkError(err, "Enqueueing kernel");
        //buffers[ibuff_poss].fromGPU();
        buffers[ibuff_out ].fromGPU(commands);
        clFinish(commands);
        buffers[ibuff_poss].release();
        buffers[ibuff_out ].release();
        //err = clFinish(commands);
        //OCL_checkError(err, "Waiting for kernel to finish");
        //double run_time = wtime() - start_time;
    }

    void initTask_mul( int ibuffA, int ibuffB, int ibuff_result ){
        //cltask_mul.setup( this, iKernell_mull, 1, Ntot*2, 16 );
        cltask_mul.setup( this, iKernell_mull, 1, 8, 1 ); printf( "WARRNING : initTask_mul() IS WRONG !!!! %i %i \n", ibuffA, ibuffB );
        cltask_mul.args = { 
            INTarg (cltask_mul.global[0]),
            BUFFarg(ibuffA),
            BUFFarg(ibuffB),
            BUFFarg(ibuff_result),
        };
    }

    void initTask_project( int ibuffAtoms, int ibuffCoefs, int ibuff_result ){
        Nvec  =(int4){(int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3]};
        cltask_project.setup( this, iKernell_project, 1, Ntot*2, 16 );
        cltask_project.args = { 
            INTarg (nAtoms),        //1
            BUFFarg(ibuffAtoms),    //2
            BUFFarg(ibuffCoefs),    //3
            BUFFarg(ibuff_result),  //4
            REFarg(Nvec),           //5
            REFarg(pos0),           //6
            REFarg(dA),           //7
            REFarg(dB),           //8
            REFarg(dC)            //9
        };
        //cltask_project.print_arg_list();
    }

    void initTask_project_tex( int ibuffAtoms, int ibuffCoefs, int ibuff_result ){
        printf("DEBUG initTask_project_tex() \n");
        Nvec  =(int4){(int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3]};
        cltask_project_tex.setup( this, iKernell_project_tex, 1, Ntot*2, 16 );
        cltask_project_tex.args = { 
            INTarg (nAtoms),        //1
            BUFFarg(ibuffAtoms),    //2
            BUFFarg(ibuffCoefs),    //3
            BUFFarg(ibuff_result),  //4
            BUFFarg(itex_basis),    //5
            REFarg(Nvec),           //6
            REFarg(pos0),           //7
            REFarg(dA),           //8
            REFarg(dB),           //9
            REFarg(dC)            //10
        };
        cltask_project_tex.print_arg_list();
    }

/*
    void initTask_projectPos_tex( int ibuffAtoms, int ibuffCoefs, int ibuff_poss, int ibuff_out ){
        cltask_projectPos_tex.setup( this, iKernell_project_tex, 1, Ntot*2, 16 );
        cltask_projectPos_tex.args = { 
            INTarg (nAtoms),        //1
            BUFFarg(ibuffAtoms),    //2
            BUFFarg(ibuffCoefs),    //3
            INTarg (nPos),        //4
            BUFFarg(ibuff_poss),  //5
            BUFFarg(ibuff_out),  //6
            BUFFarg(itex_basis),    //7
        };
        cltask_projectPos_tex.print_arg_list();
    }
*/

    void projectAtoms( float4* atoms, float4* coefs, int ibuff_result ){
        //for(int i=0; i<nAtoms;i++){printf( "atom[%i] xyz|e(%g,%g,%g|%g) coefs(%g,%g,%g|%g)\n", i, atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,  coefs[i].x,coefs[i].y,coefs[i].z,coefs[i].w  );}
        upload(ibuff_atoms,atoms);
        upload(ibuff_coefs,coefs);
        //ibuff_atoms=0;
        //ibuff_coefs=1;
        //printf("ibuff_atoms %i ibuff_coefs %i \n", ibuff_atoms, ibuff_coefs);
        //err = clEnqueueWriteBuffer  ( queue, data_cl,                   CL_TRUE, 0, buffer_size,           data,  0, NULL, NULL );
        //err = clEnqueueWriteBuffer ( commands, buffers[ibuff_atoms].p_gpu, CL_TRUE, 0, sizeof(float4)*nAtoms, atoms, 0, NULL, NULL ); OCL_checkError(err, "Creating ibuff_atoms");
        //err = clEnqueueWriteBuffer ( commands, buffers[ibuff_coefs].p_gpu, CL_TRUE, 0, sizeof(float4)*nAtoms, coefs, 0, NULL, NULL ); OCL_checkError(err, "Creating ibuff_coefs");
        clFinish(commands); 
        //initTask_project( ibuff_atoms, ibuff_coefs, ibuff_result );
        //cltask_project.enque( );
        initTask_project_tex( ibuff_atoms, ibuff_coefs, ibuff_result );
        cltask_project_tex.enque( );
        //initTask_mul( ibuff_atoms, ibuff_coefs,   ibuff_result   );
        //cltask_mul.enque( );
        clFinish(commands); 
    }
    
    void convolution( int ibuffA, int ibuffB, int ibuff_result ){
        err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffA].p_gpu, NULL, NULL);
        err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffB].p_gpu, NULL, NULL);  
        //mul_buffs( ibuffA, ibuffB, ibuff_result );
        initTask_mul( ibuffA, ibuffB, ibuff_result );  //cltask_mul.print_arg_list();
        cltask_mul.enque( );
        err = clfftEnqueueTransform( planHandle, CLFFT_BACKWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuff_result].p_gpu, NULL, NULL);  
    }

    void cleanup(){
        //clReleaseMemObject( data_cl );
        //free( data );
        err = clfftDestroyPlan( &planHandle );
        clfftTeardown( );
        //clReleaseCommandQueue( commands );
        //clReleaseContext( context );
        release_OCL();
    }


    void makeMyKernels(){
        printf( "DEBUG makeMyKernels \n" );
        buildProgram( "myprog.cl" );
        iKernell_mull    = newKernel( "mul" );
        iKernell_project = newKernel( "projectAtomsToGrid" );
        iKernell_project_tex    = newKernel( "projectAtomsToGrid_texture" );
        iKernell_projectPos_tex = newKernel( "projectWfAtPoints_tex" );
        printf( "DEBUG makeMyKernels END \n" );
        //exit(0);
    };


/*
    void runAll( ){
        //printf("DEBUG Ns[%li,%li,%li]\n", clLengths[0], clLengths[1], clLengths[2] );
        makeData();  printData( data );    printf( "DEBUG 1 \n" );
        initOCL();                         printf( "DEBUG 2 \n" );
        //printf("DEBUG Ns[%li,%li,%li]\n", clLengths[0], clLengths[1], clLengths[2] );
        initFFT( 3, clLengths );           printf( "DEBUG 3 \n" );
        loadData( data );
        runFFT();    printData(data );     printf( "DEBUG 5 \n" );
        cleanup();                         printf( "DEBUG 6 \n" );
    }
*/

    void loadWfBasis( const char* path, float Rcut, int nsamp, int ntmp, int nZ, int* iZs ){
        float* data     = new float[nsamp*2*nZ];
        float* data_tmp = new float[ntmp      ];
        char fname[64];
        float dxTmp =Rcut/ntmp;
        float dxSamp=Rcut/nsamp;
        printf( "loadWfBasis nsamp %i ntmp %i nZ %i Rcut %g \n", nsamp, ntmp, nZ, Rcut );
        for(int i=0; i<nZ; i++){
            int iz=iZs[i];
            sprintf( fname, "%s%03i_%03i.wf%i", path, iz, (int)(Rcut*100), 1 );
            loadWf_(fname, data_tmp );
            resample1D( nsamp, 0, 0, dxSamp, dxTmp, data_tmp, data+nsamp*(i*2), 2,0 );
            sprintf( fname, "%s%03i_%03i.wf%i", path, iz, (int)(Rcut*100), 2 );
            if( loadWf_(fname, data_tmp ) ){
                resample1D( nsamp, 0, 0, dxSamp, dxTmp, data_tmp, data+nsamp*(i*2), 2,1 );
            }else{
                for(int j=0; j<nsamp; j++){ data[nsamp*(i*2)+2*j+1]=data[nsamp*(i*2)+2*j]; }
            }
        }
        //for(int i=0; i<nsamp; i++){ printf("[%i] (%f,%f)   (%f,%f) \n", i, data[i*2],data[i*2+1],data[i*2+2*nsamp],data[i*2+1+2*nsamp] );  }
        delete [] data_tmp;
        itex_basis = newBufferImage2D( "BasisTable", nsamp, nZ,  sizeof(float)*2,  data, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR , {CL_RG, CL_FLOAT} );
        delete [] data;
    }

    void update_GridShape(){
        printf("saveToXsf() %g %g %g \n", dC.x, dC.y, dC.z );
        grid.cell.a=(Vec3d)(*(Vec3f*)&dA)*Ns[0];
        grid.cell.b=(Vec3d)(*(Vec3f*)&dB)*Ns[1];
        grid.cell.c=(Vec3d)(*(Vec3f*)&dC)*Ns[2];
        grid.pos0  =(Vec3d)(*(Vec3f*)&pos0);
        grid.n = (Vec3i){(int)Ns[0],(int)Ns[1],(int)Ns[2]}; //Ntot;
        grid.updateCell();
        grid.printCell();
    }

    void saveToXsf(const char* fname, int ibuff){
        update_GridShape();
        printf("saveToXsf() 1 \n");
        float* cpu_data = new float[Ntot*2]; // complex 2*float
        printf("saveToXsf() 2 %i \n", ibuff);
        download( ibuff,cpu_data);
        printf("saveToXsf() 2.5 %i \n", ibuff);
        finishRaw();
        printf("saveToXsf() 3 \n");
        grid.saveXSF( fname, cpu_data, 2, 0 );
        printf("saveToXsf() 4 \n");
        delete [] cpu_data;
    }

};

OCLfft oclfft;

extern "C" {

    void  init(){
        oclfft.init();
        oclfft.makeMyKernels();
    }

    int   upload(int i, const float* cpu_data ){ return oclfft  .upload(i,cpu_data);                        };
    int download(int i,       float* cpu_data ){ return oclfft.download(i,cpu_data);  oclfft.finishRaw();   };

    void initFFT( int ndim, size_t* Ns_ ){
        oclfft.initFFT( ndim, Ns_ );      printf( "C initFFT 1 \n" );
        oclfft.newFFTbuffer( "inputA" );  printf( "C initFFT 2 \n" );
        oclfft.newFFTbuffer( "inputB" );  printf( "C initFFT 3 \n" );
        oclfft.newFFTbuffer( "outputC" ); printf( "C initFFT 4 \n" );
        //oclfft.initTask_mul( 0, 1, 2 );    // If we know arguments in front, we may define it right now
    }

    int initAtoms( int nAtoms          ){  return oclfft.initAtoms( nAtoms ); };
    void runfft( int ibuff, bool fwd   ){ oclfft.runFFT( ibuff,fwd,0);     };
    //void runfft( int ibuff, bool fwd, float* data ){ oclfft.runFFT( ibuff,fwd, data); };
    void convolve( int ibuffA, int ibuffB, int ibuff_result ){
        oclfft.convolution( ibuffA, ibuffB, ibuff_result  );
        //oclfft.mul_buffs( ibuffA, ibuffB, ibuff_result );
    }
    void projectAtoms( float* atoms, float* coefs, int ibuff_result ){ oclfft.projectAtoms( (float4*)atoms, (float4*)coefs, ibuff_result ); }
    void cleanup(){ oclfft.cleanup(); }

    void setGridShape( float* pos0, float* dA, float* dB, float* dC ){
        oclfft.pos0=*(float4*)pos0;
        oclfft.dA  =*(float4*)dA;
        oclfft.dB  =*(float4*)dB;
        oclfft.dC  =*(float4*)dC;
        //printf( "setGridShape dA %g %g %g \n", oclfft.dA.x, oclfft.dA.y, oclfft.dA.z ); 
        //printf( "setGridShape dB %g %g %g \n", oclfft.dB.x, oclfft.dB.y, oclfft.dB.z ); 
        //printf( "setGridShape dC %g %g %g \n", oclfft.dC.x, oclfft.dC.y, oclfft.dC.z ); 
    }

    int initBasisTable( int nx, int ny, float* data ){  return oclfft.initBasisTable(nx,ny,data ); };



    void approx( int npoints, int npolys, double* xs, double* ys, double* ws ){
        //int npoints = 100;
        const int npows   = 4;
        //int npolys  = 15;
        double pows[npows] {1,2,4,8};
        //double pows [npows ]{-1,-2,-4,-8};
        //int    polys[npolys]{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        //double pows [];
        int    polys[npolys]; for(int i=0;i<npolys; i++)polys[i]=i;
        aaprox.bindOrRealloc( npolys, npows, npoints, polys, pows );
        aaprox.ws.alloc(aaprox.np);
        //aaprox.ws = aaprox.ys_ref;
        for(int i=0; i<npoints; i++){
            aaprox.xs    [i] = xs[i];
            aaprox.ys_ref[i] = ys[i];
            aaprox.ws    [i] = ws[i];
        }
        aaprox.reallocAux();
        aaprox.preparePowers();
        for(int i=0; i<aaprox.npows; i++){
            //int order = aaprox.tryVariant(10, 50, i ); // order of polynominal required for that accuracy?
            int order = aaprox.tryVariant(5, npoints, i );
            //int order = aaprox.ascendingPolyFit_(10, 50, i );
            if(order<0){ order=aaprox.npoly; printf("(not converged)"); }
            printf("[%i] pow %g err %g coefs[%i]: ", i, aaprox.pows[i], aaprox.err, order );
            //for(int j=0; j<order; j++ ) printf( "%i %g \n", i, aaprox.coefs[j] );
            for(int j=0; j<order; j++ ) printf( " %g ", aaprox.coefs[j] );
            printf("\n");
        }

    }

    void loadWf(const char* fname, float* out){ loadWf_(fname, out); };

    void loadWfBasis( const char* path, float Rcut, int nsamp, int ntmp, int nZ, int* iZs ){ oclfft.loadWfBasis(path, Rcut,nsamp,ntmp,nZ,iZs ); }

    void saveToXsf(const char* fname, int ibuff){ return oclfft.saveToXsf(fname, ibuff); }

    void initFireBall( int natoms, int* atypes, double* apos ){

        GridShape& g = oclfft.grid;

        // ======= Init Fireball
        fireCore.loadLib( "/home/prokop/git/FireCore/build/libFireCore.so" );
        fireCore.set_lvs( (double*)&( oclfft.grid.cell) );
        fireCore.preinit();
        fireCore.init   ( natoms, atypes, apos );
        //exit(0);

        // ======= Calculate Molecular Orbitals
        fireCore.assembleH( 0, 1, apos );
        double k0[3]{0.,0.,0.};
        fireCore.solveH( k0, 1  );
        double* pwfcoef; 
        fireCore.getPointer_wfcoef( &pwfcoef );
        for(int i=0; i<64; i++){ printf( "pwfcoef[%i] %g \n", i, pwfcoef[i] ); };

        int iMO=0;
        int Norb = 8;
        double* pwf = pwfcoef+Norb*iMO;

        double tmp[3]{0.,0.,0.};
        // Ecut_, ifixg0_, g0_,  |   ngrid, dCell ) 
        fireCore.setupGrid( 100.0, 0, tmp, (int*)&g.n, (double*)&g.dCell );
        //printf( "setupGrid N (%i,%i,%i) dA.x (%g,%g,%g) dB (%g,%g,%g) dC (%g,%g,%g) \n", oclfft.grid.n.x,oclfft.grid.n.x,oclfft.grid.n.x );
        g.updateCell_2();
        g.printCell();
        int ntot = g.getNtot();
        double* ewfaux = new double[ ntot ];
        fireCore.getGridMO( iMO+1, ewfaux );
        //for(int i=0; i<ntot; i+=100){ printf("%g \n", ewfaux[i] ); }
        //fireCore.getGridMO();
        g.saveXSF( "ref.xsf", ewfaux );

        //exit(0);

        // ==== Init OpenCL FFT
        oclfft.init();
        oclfft.makeMyKernels();
        size_t Ns[3]{ (size_t)g.n.x, (size_t)g.n.y, (size_t)g.n.z };
        int iZs[2]{1,6}; 
        initFFT( 3, Ns );
        oclfft.loadWfBasis( "Fdata/basis/", 4.50,100,1000, 2,iZs );
        initAtoms( natoms );
        // ==== Convert Wave-Function coefs and project using OpenCL 
        float pos0[4]{ (float)g.cell.a.x*-0.5f, (float)g.cell.b.y*-0.5f, (float)g.cell.c.z*-0.5f, 0.0};
        //printf( "pos0 (%g,%g,%g) \n", pos0[0],pos0[1],pos0[2] ); exit(0);
        float dA  [4]{ (float)g.dCell.a.x, (float)g.dCell.a.y, (float)g.dCell.a.z, 0.0};
        float dB  [4]{ (float)g.dCell.b.x, (float)g.dCell.b.y, (float)g.dCell.b.z, 0.0};
        float dC  [4]{ (float)g.dCell.c.x, (float)g.dCell.c.y, (float)g.dCell.c.z, 0.0};
        setGridShape( pos0, dA, dB, dC );
        oclfft.update_GridShape();
        float4* coefs  = new float4[natoms];
        float4* apos_  = new float4[natoms];
        Vec3d*  apos__ = (Vec3d*)apos; 
        int j=0;
        for(int i=0; i<natoms; i++){
            apos_[i] = (float4){  (float)apos__[i].x, (float)apos__[i].y, (float)apos__[i].z, atypes[i]-0.5f };
            if( atypes[i]==1 ){
                coefs[i]=(float4){0.f,0.f,0.f,(float)pwf[j]};  j++;
            }else{
                coefs[i]=(float4){(float)pwf[j+1],(float)pwf[j+2],(float)pwf[j+3],(float)pwf[j]};  j+=4;
            }
            printf( "coefs[%i] t %i | %g %g %g p|s %g \n", i, atypes[i],  coefs[i].x, coefs[i].y, coefs[i].z, coefs[i].w );
        }
        projectAtoms( (float*)apos_, (float*)coefs, 0 );
        //oclfft.saveToXsf( "test.xsf", 0 );


        oclfft.update_GridShape();
        float* cpu_data = new float[oclfft.Ntot*2]; // complex 2*float
        oclfft.download( 0, cpu_data);
        oclfft.finishRaw();

        int i0 = oclfft.Ns[0]*oclfft.Ns[1]*oclfft.Ns[3]/2 + oclfft.Ns[0]*oclfft.Ns[1]/2;
        printf( "Ntot %li i0 %i Ns (%li,%li,%li) \n", oclfft.Ntot*2, i0*2, oclfft.Ns[0],oclfft.Ns[1],oclfft.Ns[2]  );

        for(int i=0; i<44; i++){
            printf( "[%i] %g \n", i, cpu_data[ (i0+i)*2 ] );
        }




        //exit(0);
        //firecore_assembleH( iforce_, Kscf_, positions_ )
        //firecore_solveH( k_temp, ikpoint ) 

        delete [] coefs;
        delete [] apos;
        delete [] ewfaux;
        delete [] cpu_data;
    }

};