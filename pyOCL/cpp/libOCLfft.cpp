

// No need to explicitely include the OpenCL headers 
//#include <clFFT.h>
#include "OCLfft_errors.h"
#include "libOCLfft.h"
#include  "OCL.h"
#include  "approximation.h"
#include "Grid.h"
#include "FireCoreAPI.h"
#include "quaternion.h"

double const_Bohr_Radius = 0.529177210903;

inline static double dist2_PointBox( const Vec3d& p, const Vec3d& a, const Vec3d& b ){
    // from here : http://stackoverflow.com/questions/4578967/cube-sphere-intersection-test
    // assume C1 and C2 are element-wise sorted, if not, do that now
    double dist2 = 0.0;
    if (p.x < a.x){ dist2 += sq(p.x - a.x); }else if(p.x > b.x){ dist2 += sq(p.x - b.x); };
    if (p.y < a.y){ dist2 += sq(p.y - a.y); }else if(p.y > b.y){ dist2 += sq(p.y - b.y); };
    if (p.z < a.z){ dist2 += sq(p.z - a.z); }else if(p.z > b.z){ dist2 += sq(p.z - b.z); };
    return dist2;
}

int loadWf_(const char* fname, float* out){
    const int nbuff = 1024;
    char buff[nbuff]; char* line;
    //printf( "loadWf %s \n", fname );
    FILE *pFile = fopen(fname, "r" );
    if(pFile==0)return -1;
    //printf( "loadWf_ 0 \n" );
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    printf( "loadWf_ 1 \n" );
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
    return n;
}

void resample1D( int nin, float x0in, float x1in, float* from,   int nout, float x0out, float x1out, float* to,   int pitch, int off ){
    float dx_in    = (x1in -x0in )/nin;
    float dx_out   = (x1out-x0out)/nout;
    float invdx_in = 1/dx_in;
    //printf("resample1D  dx_in %f invdx_in %f x1in %f \n", dx_in, invdx_in, x1in );
    //for(int i=0; i<nin;  i++){ printf( "[%i] %g \n", i, from[i] );  }
    for(int i=0; i<nout; i++){
        float  x = x0out + i*dx_out;
        int iout = i*pitch+off;
        if( x>x1in ){
            to[iout] =0;
        }else{
            float  u = (x-x0in)*invdx_in;
            int   iu = (int)u;
            float  f = u-iu;
            to[iout] =  from[iu]*(1-f) + from[iu+1]*f;  // lerp
            //printf("resample1D [%i] x %g y %g xu %g \n", i, x, to[iout], u );
        }
        //printf("resample1D [%i] x %g y %g x1out %g \n", i, x, to[iout], x1out );
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
    void fitlocal( ){  for(cl_uint i=0; i<dim; i++){ global[i]=((int)(global[i]/(1.0*local[i]))+1)*local[i]; };  }
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
    int iKernell_project_dens_tex=-1;
    int iKernell_projectPos_tex=-1;
    int iKernell_poissonW=-1;

    OCLtask cltask_mul;
    OCLtask cltask_poissonW;
    OCLtask cltask_project;
    OCLtask cltask_project_tex;
    OCLtask cltask_project_den_tex;
    OCLtask cltask_projectPos_tex;
    int itex_basis=-1;

    int    nAtoms=0;
    int    nOrbs=0;
    int    nPos=0;
    float4 dcell_poisson{1.f,1.f,1.f,1.f};
    float4 pos0, dA, dB, dC;
    GridShape grid;
    int ibuff_atoms=0,ibuff_coefs=0;
    //int ibuff_poss;

    void updateNtot(){
        Ntot=1; for(int i=0; i<ndim; i++){ Ntot*=Ns[i]; };
    }

    void planFFT(){
        // Create a default plan for a complex FFT. 
        // https://github.com/clMathLibraries/clFFT/issues/148
        // The dimensions have to be powers of 2,3,5,7,11,13 or any combination of those.
        printf(  "planFFT fft_dim %i N(%li,%li,%li,%li)  \n", fft_dim, Ns[0], Ns[1], Ns[2], Ns[3] );
        err = clfftCreateDefaultPlan(&planHandle, context, fft_dim, Ns );        OCLfft_checkError(err, "clfftCreateDefaultPlan" );
        err = clfftSetPlanPrecision (planHandle, CLFFT_SINGLE);                  OCLfft_checkError(err, "clfftSetPlanPrecision" );
        err = clfftSetLayout        (planHandle, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED);   OCLfft_checkError(err, "clfftSetLayout" );
        err = clfftSetResultLocation(planHandle, CLFFT_INPLACE);                 OCLfft_checkError(err, "clfftSetResultLocation" );
        err = clfftBakePlan         (planHandle, 1, &commands, NULL, NULL);      OCLfft_checkError(err, "clfftBakePlan" );
    }

    void initFFT( int ndim, size_t* Ns_ ){
        //printf("DEBUG initFFT() ndim %i Ns[%li,%li,%li]\n", ndim, Ns[0], Ns[1], Ns[2] );
        if     (ndim==1){ fft_dim=CLFFT_1D; }
        else if(ndim==2){ fft_dim=CLFFT_2D; }
        else if(ndim==3){ fft_dim=CLFFT_3D; };
        //buffer_size  = sizeof(float2);
        Ntot=1; for(int i=0; i<ndim;i++){  Ns[i]=Ns_[i];  Ntot*= Ns[i]; }
        printf( "initFFT ndim %i Ntot %li Ns[%li,%li,%li]\n", ndim, Ntot, Ns[0],Ns[1],Ns[2] );
        clfftSetupData fftSetup;                //printf("initFFT 1 \n");
        err  = clfftInitSetupData(&fftSetup);   OCLfft_checkError(err, "clfftInitSetupData");
        err  = clfftSetup        (&fftSetup);   OCLfft_checkError(err, "clfftSetup" );
        //data_cl = clCreateBuffer( context, CL_MEM_READ_WRITE, buffer_size, NULL, &err );
        planFFT(  );                            //printf("initFFT 4 \n");
    }

    int newFFTbuffer( char* name ){
        return newBuffer( name, Ntot, sizeof(float2), 0, CL_MEM_READ_WRITE );
    }

    int initAtoms( int nAtoms_, int nOrbs_ ){
        nAtoms=nAtoms_;
        nOrbs  =nOrbs_;
        //printf("DEBUG initAtoms nAtoms %i nOrbs %i \n", nAtoms, nOrbs );
        ibuff_atoms   =newBuffer( "atoms",    nAtoms,       sizeof(float4), 0, CL_MEM_READ_ONLY );
        ibuff_coefs   =newBuffer( "coefs",    nAtoms*nOrbs, sizeof(float4), 0, CL_MEM_READ_ONLY );
        //ibuff_coefsAll=newBuffer( "coefsAll", nAtoms*nOrbs, sizeof(float4), 0, CL_MEM_READ_ONLY );
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
        OCLfft_checkError(err, " clfftEnqueueTransform " );                                                  
        if(data){
            err = clEnqueueReadBuffer  ( commands, data_cl, CL_TRUE, 0, buffers[ibuff].byteSize(), data, 0, NULL, NULL ); // Fetch results of calculations. 
            OCLfft_checkError(err, " clEnqueueReadBuffer " );
            err = clFinish(commands);    // Wait for calculations to be finished. 
            OCLfft_checkError(err, " clFinish " );
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


    void projectAtomPosTex(  float4* atoms, float4* coefs, int nPos, float4* poss, float2* out ){
        KernelDims kdim;
        kdim.dim        = 1;
        kdim.global[0]  = nPos;
        kdim.local [0]  = 16;
        kdim.fitlocal( ); printf( "projectAtomPosTex %li \n", kdim.global[0] );
        cl_kernel kernel = kernels[iKernell_projectPos_tex]; 
        //for(int i=0; i<nPos; i++){ printf("projectAtomPosTex %i (%g,%g,%g)\n", i, poss[i].x, poss[i].y, poss[i].z  );  };
        //printf("DEBUG projectAtomPosTex() 0 \n");
        int ibuff_poss = newBuffer( "poss", nPos, sizeof(float4), (float*)poss, CL_MEM_READ_WRITE );
        int ibuff_out  = newBuffer( "out" , nPos, sizeof(float2), (float*)out , CL_MEM_READ_WRITE );
        //printf("DEBUG projectAtomPosTex() 1 \n");
        buffers[ibuff_poss].toGPU(commands);
        upload(ibuff_atoms,atoms);
        upload(ibuff_coefs,coefs);
        //printf("DEBUG projectAtomPosTex() 2 \n");
        //buffers[ibuff_out ].toGPU(commands);
        //err = clEnqueueWriteBuffer ( commands, buffers[ibuff_poss].p_gpu, CL_TRUE, 0, sizeof(float4)*nPos, poss, 0, NULL, NULL );
        //err = clEnqueueWriteBuffer ( commands, buffers[ibuff_out ].p_gpu, CL_TRUE, 0, sizeof(float2)*nPos, out,  0, NULL, NULL );
        err =  clSetKernelArg(kernel, 0, sizeof(int),    &nAtoms        );
        err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &(buffers[ibuff_atoms].p_gpu) );
        err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &(buffers[ibuff_coefs].p_gpu) );
        err =  clSetKernelArg(kernel, 3, sizeof(int),    &nPos          );
        err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &(buffers[ibuff_poss].p_gpu) );
        err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &(buffers[ibuff_out ].p_gpu) );
        err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), &(buffers[itex_basis].p_gpu) );
        //printf("DEBUG projectAtomPosTex() 3 \n");
        //checkError(err, "Setting kernel args");
        //double start_time = wtime();
        //printf( "mul_buffs kdim: dim %i global %li local %li \n", kdim.dim, kdim.global[0], kdim.local[0] ); 
        err = clEnqueueNDRangeKernel(  commands, kernel,   kdim.dim, NULL, kdim.global, kdim.local, 0, NULL, NULL);    
        OCL_checkError(err, "Enqueueing kernel");
        //buffers[ibuff_poss].fromGPU();
        buffers[ibuff_out ].fromGPU(commands);
        //printf("DEBUG projectAtomPosTex() 4 \n");
        clFinish(commands);
        buffers[ibuff_poss].release();
        buffers[ibuff_out ].release();
        //printf("DEBUG projectAtomPosTex() 5 \n");
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

    void initTask_poissonW( int ibuffA, int ibuff_result ){
        printf( "BEGIN initTask_poissonW \n" );
        //cltask_poissonW.setup( this, iKernell_poissonW, 1, Ntot, 1 );
        cltask_poissonW.setup4( this, iKernell_poissonW, 3, *(size_t4*)Ns, (size_t4){1,1,1,1} );
        cltask_poissonW.args = { 
            INTarg ((int)Ntot),
            BUFFarg(ibuffA),
            BUFFarg(ibuff_result),
            REFarg(dcell_poisson)           
        };
        printf( "END initTask_poissonW \n" );
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

    void initTask_project_dens_tex( int ibuffAtoms, int ibuffCoefs, int ibuff_result ){
        printf("DEBUG initTask_project_dens_tex() \n");
        Nvec  =(int4){(int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3]};
        cltask_project_den_tex.setup( this, iKernell_project_dens_tex, 1, Ntot*2, 16 );
        cltask_project_den_tex.args = { 
            INTarg (nAtoms),        //1
            INTarg (0),             //2
            INTarg (0),             //3
            BUFFarg(ibuffAtoms),    //4
            BUFFarg(ibuffCoefs),    //5
            BUFFarg(ibuff_result),  //6
            //BUFFarg(ibuffCoefs),  //6
            BUFFarg(itex_basis),    //7
            REFarg(Nvec),           //8
            REFarg(pos0),           //9
            REFarg(dA),           //10
            REFarg(dB),           //11
            REFarg(dC)            //12
        };
        printf("DEBUG cltask_project_den_tex.args.size() %li  \n", cltask_project_den_tex.args.size() );
        printf("DEBUG initTask_project_dens_tex() END \n");
    }


    /*
    void initTask_projectPos_tex( int ibuffAtoms, int ibuffCoefs, int ibuff_poss, int ibuff_out ){
        cltask_projectPos_tex.setup( this, iKernell_project_tex, 1, Ntot*2, 16 );
        cltask_projectPos_tex.args = { 
            INTarg (nAtoms),        //1
            BUFFarg(ibuffAtoms),    //2
            BUFFarg(ibuffCoefs),    //3
            INTarg (nPos),          //4
            BUFFarg(ibuff_poss),    //5
            BUFFarg(ibuff_out),     //6
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

    void projectAtomsDens( float4* atoms, float4* coefs, int ibuff_result, int iorb1, int iorb2 ){
        printf( "DEBUG projectAtomsDens(%i,%i) | atoms* %li long* %li \n", iorb1, iorb2,  (long)atoms, (long)coefs );
        if( atoms ) upload(ibuff_atoms,atoms); 
        if( coefs ) upload(ibuff_coefs,coefs);
        clFinish(commands); 
        //initTask_project_dens_tex
        printf( "DEBUG projectAtomsDens 1 \n" );
        initTask_project_dens_tex( ibuff_atoms, ibuff_coefs, ibuff_result );
        //cltask_project_den_tex
        printf( "DEBUG projectAtomsDens 2 \n" );
        cltask_project_den_tex.args[1].i=iorb1;
        printf( "DEBUG projectAtomsDens 3 \n" );
        cltask_project_den_tex.args[2].i=iorb2;
        printf( "DEBUG projectAtomsDens 4 \n" );
        cltask_project_den_tex.print_arg_list();
        printf( "DEBUG projectAtomsDens 5 \n" );
        cltask_project_den_tex.enque( );
        printf( "cltask_project_den_tex.enque( ) END \n" );
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

    void poisson( int ibuffA, int ibuff_result, float4* dcell=0 ){
        printf( "BEGIN poisson %i -> %i ( %s -> %s ) \n", ibuffA, ibuff_result, buffers[ibuffA].name, buffers[ibuff_result].name );
        //err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffA].p_gpu, NULL, NULL);
        //OCLfft_checkError(err, " poisson::clfftEnqueueTransform " );
        if( dcell ){ dcell_poisson = *dcell; }
        initTask_poissonW( ibuffA, ibuff_result );      printf( "DEBUG 1 ");
        finishRaw(); saveToXsf( "rho.xsf", ibuffA );    printf( "DEBUG 2 ");
        err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffA].p_gpu, NULL, NULL);
        OCLfft_checkError(err, " poisson::clfftEnqueueTransform " ); printf( "DEBUG 4 ");
        finishRaw(); saveToXsf( "rho_w.xsf", ibuffA );               printf( "DEBUG 5 ");
        cltask_poissonW.enque( );                                    printf( "DEBUG 6 ");
        finishRaw(); saveToXsf( "Vw.xsf", ibuff_result );            printf( "DEBUG 7 ");
        err = clfftEnqueueTransform( planHandle, CLFFT_BACKWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuff_result].p_gpu, NULL, NULL);  
        finishRaw(); saveToXsf( "V.xsf", ibuff_result );             printf( "DEBUG 8 ");
        printf( "END poisson \n" );
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


    void makeMyKernels( const char*  cl_src_dir ){
        char srcpath[1024];
        sprintf( srcpath, "%s/myprog.cl", cl_src_dir );
        printf( "DEBUG makeMyKernels %s \n", srcpath );
        buildProgram( srcpath );
        iKernell_mull             = newKernel( "mul" );
        iKernell_poissonW         = newKernel( "poissonW" );
        iKernell_project          = newKernel( "projectAtomsToGrid" );
        iKernell_project_tex      = newKernel( "projectAtomsToGrid_texture"  );
        iKernell_project_dens_tex = newKernel( "projectOrbDenToGrid_texture" );
        iKernell_projectPos_tex   = newKernel( "projectWfAtPoints_tex" );
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

    void loadWfBasis( const char* path, float RcutSamp, int nsamp, int ntmp, int nZ, int* iZs, float* Rcuts ){
        float* data_tmp = new float[ntmp      ];
        float* data     = new float[nsamp*2*nZ];
        char fname[64];
        //float dxTmp =(Rcut*const_Bohr_Radius)/ntmp;
        //float dxSamp=Rcut/nsamp;
        //RcutSamp*=0.529177210903f;
        printf( "loadWfBasis nsamp %i ntmp %i nZ %i RcutSamp %g [A]\n", nsamp, ntmp, nZ, RcutSamp );
        for(int i=0; i<nZ; i++){
            int iz=iZs[i];
            float Ri = Rcuts[i];
            sprintf( fname, "%s%03i_%03i.wf%i", path, iz, (int)(Ri*100), 1 );
            int nin = loadWf_(fname, data_tmp );
            //resample1D( nsamp, 0, 0, dxSamp, dxTmp, data_tmp, data+nsamp*(i*2), 2,0 );
            resample1D( nin, 0.0, Ri*const_Bohr_Radius, data_tmp,   nsamp,0.0, RcutSamp, data+nsamp*(i*2),   2,0 );
            sprintf( fname, "%s%03i_%03i.wf%i", path, iz, (int)(Ri*100), 2 );
            if( loadWf_(fname, data_tmp ) ){
                resample1D( nin, 0.0, Ri*const_Bohr_Radius, data_tmp,   nsamp,0.0, RcutSamp, data+nsamp*(i*2),   2,1 );
            }else{
                for(int j=0; j<nsamp; j++){ data[nsamp*(i*2)+2*j+1]=data[nsamp*(i*2)+2*j]; }
            }
        }
        //for(int i=0; i<nsamp; i++){ printf("basis [%i] (%f,%f)   (%f,%f) \n", i, data[i*2],data[i*2+1],data[i*2+2*nsamp],data[i*2+1+2*nsamp] );  }
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

    void saveToXsf(const char* fname, int ibuff, int natoms=0, int* atypes=0, Vec3d* apos=0 ){
        printf("DEBUG saveToXsf 1 \n");
        update_GridShape();
        printf("DEBUG saveToXsf 2 \n");
        float* cpu_data = new float[Ntot*2]; // complex 2*float
        download( ibuff,cpu_data);
        printf("DEBUG saveToXsf 3 \n");
        finishRaw();
        printf("DEBUG saveToXsf 4 \n");
        //Vec3d pos0=grid.pos0;
        //grid.pos0=Vec3dZero;
        grid.saveXSF( fname, cpu_data, 2, 0,   natoms,atypes,apos );
        //grid.pos0=pos0;
        delete [] cpu_data;
    }

    void convOrbCoefs( int natoms, int* iZs, double* ocoefs, float4* coefs ){
        int io=0;
        for(int ia=0; ia<natoms; ia++){
            if(iZs[ia]==1){ 
                coefs[ia]=(float4){0.f,0.f,0.f, (float)ocoefs[io] };
                io+=1; 
            }else{ 
                coefs[ia]=(float4){ (float)ocoefs[io+3],(float)ocoefs[io+1],(float)ocoefs[io+2], (float)ocoefs[io] };   //  Fireball order:   s,py,pz,px   see https://nanosurf.fzu.cz/wiki/doku.php?id=fireball
                io+=4; 
            }
            //printf( "CPU [%i] coef(%g,%g,%g,%g)\n", ia, coefs[ia].x, coefs[ia].y, coefs[ia].z, coefs[ia].w );
        }
    }

    void convCoefs( int natoms, int* iZs, int* ityps, double* ocoefs, double* oatoms, int iorb0, int iorb1, bool bInit=false ){
        printf( "DEBUG convCoefs 1 \n" );
        // --- Count orbitals
        int norb=0;
        printf( "DEBUG convCoefs 2 \n" );
        for(int ia=0; ia<natoms; ia++){ if(iZs[ia]==1){ norb+=1; }else{ norb+=4; }; }
        //printf( "DEBUG convCoefs 2.5 norb %i \n", norb );
        //int ncoef=natoms*(iorb1-iorb0+1);
        int ncoef=natoms*norb;
        float4* atoms = new float4[ natoms ];
        float4* coefs = new float4[ ncoef  ];
        // --- trasnform positions
        for(int ia=0; ia<natoms; ia++){
            int i3=ia*3;
            atoms[ia]=(float4){ (float)oatoms[i3],(float)oatoms[i3+1],(float)oatoms[i3+2], ityps[ia]+0.1f };
        }
        printf( "DEBUG convCoefs 3 norb %i ncoef %i \n", norb, ncoef );
        // --- trasnform coefs
        //int io=norb*iorb0;
        //for(int iorb=iorb0; iorb<=iorb1; iorb++){
        for(int iorb=0; iorb<norb; iorb++){
            convOrbCoefs( natoms, iZs, ocoefs+iorb*norb, coefs+iorb*natoms );
        } // iorb
        
        /*
        printf( "DEBUG  print CPU ocoefs | norb %i  \n", norb );
        for(int iorb=0; iorb<norb; iorb++){
            printf( "ORB[%i] ", iorb );
            for (int j=0; j<norb; j++ ){
                int io = iorb*norb + j;
                printf( " %g ", ocoefs[io] );
            }
            printf( "\n");
        }
        printf( "DEBUG  print CPU coefs iorb0,iorb1,nAtoms %i, %i, %i \n", iorb0,iorb1,natoms );
        for(int iorb=iorb0; iorb<=iorb1; iorb++){
            for (int ia=0; ia<natoms; ia++ ){
                int io = iorb*natoms + ia;
                printf( "CPU [%i,%i] atom(%g,%g,%g,,%g) coef(%g,%g,%g,%g)\n", iorb, ia,  atoms[ia].x, atoms[ia].y, atoms[ia].z, atoms[ia].w,  coefs[io].x, coefs[io].y, coefs[io].z, coefs[io].w );
            }
        }
        */
        printf( "DEBUG convCoefs 4 \n" );
        if( bInit ){
            //printf( "DEBUG convCoefs 5 \n" );
            //nAtoms = natoms;
            //nOrbs  =  norb;
            initAtoms( natoms, norb );
        }
        printf( "DEBUG convCoefs 6 \n" );
        upload(ibuff_atoms,atoms, natoms );
        upload(ibuff_coefs,coefs, ncoef  );
        printf( "DEBUG convCoefs 7 \n" );
        delete [] atoms;
        delete [] coefs;
        printf( "DEBUG convCoefs END \n" );
    };

    void countOrbs( int natoms, int* iZs, int* offsets ){    
        int io=0;
        for(int i=0; i<natoms; i++){
            if(iZs[i]==1){ io++; }else{ io+=4; }
            offsets[i]=io;
        }
    }

    void buildDenmat( int nsel, int* sel, int* iZs, int* i0Cs, double* ocoefs, int iorb0, int iorb1, Quat4f* rho ){
        for(int i=0; i<(nsel*nsel*4); i++){ rho[i] = Quat4fZero; }
        Quat4f lcoefs[nsel];
        for(int iorb=iorb0; iorb<iorb1; iorb++ ){
            for(int i=0; i<nsel; i++){
                int ia = sel [i];
                int io = i0Cs[ia];
                if(iZs[ia]==1){ lcoefs[ia]=(Quat4f){0.f,0.f,0.f, (float)ocoefs[io] };}
                else          { lcoefs[ia]=(Quat4f){ (float)ocoefs[io+3],(float)ocoefs[io+1],(float)ocoefs[io+2], (float)ocoefs[io] }; }
            }
            for(int i=0; i<nsel; i++){
                Quat4f qi = lcoefs[i];
                for(int j=0; j<nsel; j++){
                    Quat4f qj = lcoefs[j];
                    rho[i  ].add_mul( qj, qi.x);
                    rho[i+1].add_mul( qj, qi.y);
                    rho[i+2].add_mul( qj, qi.z);
                    rho[i+4].add_mul( qj, qi.w);
                }
            } 
        }
    }

    int atoms2box( Vec3d p0, Vec3d p1, double Rcut, int natoms, Vec3d* apos, int* sel ){
        double R2=Rcut*Rcut;
        //int sel[natoms];
        int nsel=0;
        for(int i=0; i<natoms; i++){
            double r2 = dist2_PointBox( apos[i], p0,p1 );
            if( r2<R2 ){ sel[nsel]=i; nsel++; } 
        }
        return nsel;
    }

    void projectDenmat( int natoms, int* iZs, int* ityps, double* ocoefs, double* apos, int iorb0, int iorb1, double Rcut, bool bInit=false ){
        int sel [natoms];
        int i0Cs[natoms];
        countOrbs( natoms, iZs, i0Cs );

        Vec3d lbox  = (Vec3d){3.0,3.0,3.0};
        Vec3d cell  = (Vec3d){ grid.cell.a.x, grid.cell.b.y, grid.cell.c.z };
        Vec3i nbox  = (Vec3i){ (int)(1+cell.x/lbox.x), (int)(1+cell.y/lbox.y), (int)(1+cell.z/lbox.z)  };
        Vec3d dcell = (Vec3d){ cell.x/nbox.x, cell.y/nbox.y, cell.z/nbox.z, };
        for(int ix=0; ix<nbox.x; ix++){
            for(int iy=0; iy<nbox.y; iy++){
                for(int iz=0; iz<nbox.z; iz++){
                    Vec3d p0 = grid.pos0 + dcell*((Vec3d){(double)ix,(double)iy,(double)iz});
                    Vec3d p1 = p0 + dcell;
                    int nsel = atoms2box( p0, p1, Rcut, natoms, (Vec3d*)apos, sel );
                    float* rho = new float[4*4*nsel*nsel];
                    buildDenmat( nsel, sel, iZs, i0Cs, ocoefs, iorb0, iorb1, (Quat4f*)rho );
                    delete [] rho;
                }
            }
        }
    }

};

OCLfft oclfft;

extern "C" {

    void setErrorCheck(int ierr){
        bOCLCheckError = ierr>0;
    }

    void  init( const char* cl_src_dir ){
        oclfft.init();
        oclfft.makeMyKernels( cl_src_dir );
    }

    int   upload(int i, const float* cpu_data ){ return oclfft  .upload(i,cpu_data);                        };
    int download(int i,       float* cpu_data ){ return oclfft.download(i,cpu_data);  oclfft.finishRaw();   };

    int   upload_d(int ibuf, const double* data, bool bComplex ){ 
        int n=oclfft.buffers[ibuf].n; 
        printf( "DEBUG upload_d ibuf %i bComplex %i  n %i \n", ibuf, bComplex, n );
        float2* tmp = new float2[ n ];
        if  (bComplex){ for(int i=0; i<n; i++){ tmp[i]=(float2){(float)data[i*2],(float)data[i*2+1]};  } }
        else          { for(int i=0; i<n; i++){ tmp[i]=(float2){(float)data[i  ], 0.0f             };  } }
        //int nxy = oclfft.Ns[0]*oclfft.Ns[1];
        //if( (i%nxy)==0 ) printf( "CPU iz %i i %i data %g A(%g,%g) \n", i/nxy, i, data[i], tmp[i].x, tmp[i].y );
        return oclfft.upload(ibuf,tmp); 
        delete [] tmp;
    }

    void initFFT( int ndim, size_t* Ns_ ){
        oclfft.initFFT( ndim, Ns_ );      printf( "C initFFT 1 \n" );
        oclfft.newFFTbuffer( "inputA" );  printf( "C initFFT 2 \n" );
        oclfft.newFFTbuffer( "inputB" );  printf( "C initFFT 3 \n" );
        oclfft.newFFTbuffer( "outputC" ); printf( "C initFFT 4 \n" );
        //oclfft.initTask_mul( 0, 1, 2 );    // If we know arguments in front, we may define it right now
    }

    int initAtoms( int nAtoms, int nOrbs ){  return oclfft.initAtoms( nAtoms, nOrbs ); };
    void runfft( int ibuff, bool fwd     ){ oclfft.runFFT( ibuff,fwd,0);     };
    //void runfft( int ibuff, bool fwd, float* data ){ oclfft.runFFT( ibuff,fwd, data); };
    void convolve( int ibuffA, int ibuffB, int ibuff_result ){  oclfft.convolution( ibuffA, ibuffB, ibuff_result  );}
    void poisson ( int ibuffA, int ibuff_result, float* dcell ){  oclfft.poisson( ibuffA, ibuff_result, (float4*)dcell );}
    void projectAtoms    ( float* atoms, float* coefs, int ibuff_result                       ){ oclfft.projectAtoms    ( (float4*)atoms, (float4*)coefs, ibuff_result ); }
    void projectAtomsDens( float* atoms, float* coefs, int ibuff_result, int iorb1, int iorb2 ){ oclfft.projectAtomsDens( (float4*)atoms, (float4*)coefs, ibuff_result, iorb1, iorb2 ); }

    // void projectAtomPosTex(  float4* atoms, float4* coefs, int nPos, float4* poss, float2* out ){
    void projectAtomPosTex( float* atoms, float* coefs, int nPos, float* poss, float* out ){ oclfft.projectAtomPosTex( (float4*)atoms, (float4*)coefs,  nPos, (float4*)poss, (float2*)out ); }

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

    void convCoefs( int natoms, int* iZs, int* ityps, double* ocoefs, double* oatoms, int iorb0, int iorb1, bool bInit ){  oclfft.convCoefs( natoms, iZs, ityps, ocoefs, oatoms, iorb0, iorb1, bInit ); }

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

    void loadWfBasis( const char* path, float RcutSamp, int nsamp, int ntmp, int nZ, int* iZs, float* Rcuts ){ oclfft.loadWfBasis(path, RcutSamp,nsamp,ntmp,nZ,iZs,Rcuts ); }

    void saveToXsf     (const char* fname, int ibuff){ return oclfft.saveToXsf(fname, ibuff,0,0,0); }
    void saveToXsfAtoms(const char* fname, int ibuff, int natoms, int* atypes, double* apos ){ return oclfft.saveToXsf(fname, ibuff, natoms,atypes,(Vec3d*)apos); }

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
        oclfft.makeMyKernels( "../cl" );
        size_t Ns[3]{ (size_t)g.n.x, (size_t)g.n.y, (size_t)g.n.z };
        int   iZs[2]{1,6}; 
        float Rcuts[2]{4.50,4.50}; 
        initFFT( 3, Ns );
        oclfft.loadWfBasis( "Fdata/basis/", 4.50,100,1000, 2,iZs, Rcuts );
        initAtoms( natoms, 1 );
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