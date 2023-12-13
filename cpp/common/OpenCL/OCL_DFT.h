#ifndef  OCL_DFT_h
#define  OCL_DFT_h

#include "OCLfft_errors.h"
#include <clFFT.h>
#include "OCL.h"
#include "Grid.h"
#include "IO_utils.h"
#include "quaternion.h"

#include "VecN.h"

double const_Bohr_Radius = 0.529177210903;

void v2f4( const Vec3d& v, float4& f4 ){ f4.x=(float)v.x; f4.y=(float)v.y; f4.z=(float)v.z; };
//void v2f4( const Vec3d& v, cl_float4& f4 ){ f4.s[0]=(cl_float)v.x; f4.s[1]=(cl_float)v.y; f4.s[2]=(cl_float)v.z; };
cl_float4 cl_f4( const Vec3d& v ){ return (cl_float4){(cl_float)v.x,(cl_float)v.y,(cl_float)v.z,0.f}; };

//void print_(const Vec3d& v){ printf("(%g,%g,%g)\n", v.x,v.y,v.z ); };

inline static double dist2_PointBox( const Vec3d& p, const Vec3d& a, const Vec3d& b ){
    // from here : http://stackoverflow.com/questions/4578967/cube-sphere-intersection-test
    // assume C1 and C2 are element-wise sorted, if not, do that now
    double dist2 = 0.0;
    if (p.x < a.x){ dist2 += sq(p.x - a.x); }else if(p.x > b.x){ dist2 += sq(p.x - b.x); };
    if (p.y < a.y){ dist2 += sq(p.y - a.y); }else if(p.y > b.y){ dist2 += sq(p.y - b.y); };
    if (p.z < a.z){ dist2 += sq(p.z - a.z); }else if(p.z > b.z){ dist2 += sq(p.z - b.z); };
    return dist2;
}

/**
 * @brief Loads basis-function from a file.
 * 
 * The waveform data is expected to be in a specific format, with each line containing four floating-point values.
 * The function replaces any 'D' characters in the lines with 'e' characters before parsing the values.
 * The function continues reading lines until it encounters a line that does not contain four values.
 * 
 * @param fname file name of the wavefunction data file
 * @param out   The array to store the wavefunction data
 * @return      The total number of wavefunction values read from the file
 */
int loadWf_(const char* fname, float* out){
    const int nbuff = 1024;
    char buff[nbuff]; char* line;
    //printf( "loadWf %s \n", fname );
    FILE *pFile = fopen(fname, "r" );
    if(pFile==0)return -1;
    //printf( "loadWf_ 0 \n" );
    // skip header 
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    line=fgets(buff,nbuff,pFile);
    //double xs[4];
    int n=0;
    while(true){ // 
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

/**
 * Resamples a 1D array from the input range [x0in, x1in] to the output range [x0out, x1out].
 * 
 * @param nin   number of elements in the input array
 * @param x0in  starting value of the input range
 * @param x1in  ending value of the input range
 * @param from  input array
 * @param nout  number of elements in the output array
 * @param x0out starting value of the output range
 * @param x1out ending value of the output range
 * @param to    output array
 * @param pitch pitch of the output array
 * @param off   offset of the output array
 */
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

struct KernelDims{
    cl_uint  dim;
    size_t global[3];
    size_t local [3];
    //cl_int blocksize;
    void fitlocal( ){  for(cl_uint i=0; i<dim; i++){ global[i]=((int)(global[i]/(1.0*local[i]))+1)*local[i]; };  }
};

//=======================================================================
//=======================================================================

class OCL_DFT: public OCLsystem { public:

    clfftPlanHandle planHandle;
    clfftDim fft_dim = CLFFT_3D;

    int ndim=0;
    size_t Ns[4]; // = {N0, N1, N2};
    size_t Ntot;
    int4   Nvec;

    int iKernell_mull=-1;
    int iKernell_roll=-1;
    int iKernell_grad=-1;
    int iKernell_lincomb=-1;
    int iKernell_project=-1;
    int iKernell_project_tex=-1;
    int iKernell_project_dens_tex=-1;
    int iKernell_project_atom_dens_tex=-1;
    int iKernell_projectPos_tex=-1;
    int iKernell_poissonW=-1;
    int iKernell_gradient=-1;

    OCLtask cltask_mul;
    OCLtask cltask_lincomb;
    OCLtask cltask_poissonW;
    OCLtask cltask_gradient;
    OCLtask cltask_project;
    OCLtask cltask_project_tex;
    OCLtask cltask_project_den_tex;
    OCLtask cltask_project_atom_dens_tex;
    OCLtask cltask_projectPos_tex;
    int itex_basis=-1;

    int    nAtoms=0;
    int    nOrbs=0;
    int    nPos=0;
    float4 dcell_poisson{1.f,1.f,1.f,1.f};
    float4 dcell_gradient{1.f,1.f,1.f,1.f};
    float4 pos0, dA, dB, dC;
    float2 acumCoef = (float2){0.0,1.0};
    GridShape grid;
    int ibuff_atoms=-1,ibuff_coefs=-1,ibuff_aforces=-1,ibuff_neighs=-1,ibuff_neighCell=-1;

    int     nAtype=0;
    int*    atype_nOrb  =0; // number of orbitals per atomic type (hydrogen=1(s), carbon=4(s,px,py,pz))
    float2* atype_Qconfs=0; //

    void updateNtot(){
        Ntot=1; for(int i=0; i<ndim; i++){ Ntot*=Ns[i]; };
    }

    void planFFT(){
        // Create a default plan for a complex FFT. 
        // https://github.com/clMathLibraries/clFFT/issues/148
        // The dimensions have to be powers of 2,3,5,7,11,13 or any combination of those.
        if(verbosity>0)printf(  "planFFT fft_dim %i N(%li,%li,%li,%li)  \n", fft_dim, Ns[0], Ns[1], Ns[2], Ns[3] );
        int err=0;
        err = clfftCreateDefaultPlan(&planHandle, context, fft_dim, Ns );        OCLfft_checkError(err, "clfftCreateDefaultPlan" );
        err = clfftSetPlanPrecision (planHandle, CLFFT_SINGLE);                  OCLfft_checkError(err, "clfftSetPlanPrecision" );
        err = clfftSetLayout        (planHandle, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED);   OCLfft_checkError(err, "clfftSetLayout" );
        err = clfftSetResultLocation(planHandle, CLFFT_INPLACE);                 OCLfft_checkError(err, "clfftSetResultLocation" );
        err = clfftBakePlan         (planHandle, 1, &commands, NULL, NULL);      OCLfft_checkError(err, "clfftBakePlan" );
    }

    int setNs(int ndim, int* Ns_ ){
        if     (ndim==1){ fft_dim=CLFFT_1D; }
        else if(ndim==2){ fft_dim=CLFFT_2D; }
        else if(ndim==3){ fft_dim=CLFFT_3D; };
        Ntot=1; for(int i=0; i<ndim;i++){  Ns[i]=Ns_[i];  Ntot*= Ns[i]; }
        return Ntot;
    }

    void initFFT( int ndim, size_t* Ns_ ){
        //printf("DEBUG initFFT() ndim %i Ns[%li,%li,%li]\n", ndim, Ns[0], Ns[1], Ns[2] );
        if     (ndim==1){ fft_dim=CLFFT_1D; }
        else if(ndim==2){ fft_dim=CLFFT_2D; }
        else if(ndim==3){ fft_dim=CLFFT_3D; };
        //buffer_size  = sizeof(float2);
        Ntot=1; for(int i=0; i<ndim;i++){  Ns[i]=Ns_[i];  Ntot*= Ns[i]; }
        //printf( "initFFT ndim %i Ntot %li Ns[%li,%li,%li]\n", ndim, Ntot, Ns[0],Ns[1],Ns[2] );
        clfftSetupData fftSetup;                //printf("initFFT 1 \n");
        int err=0;
        err  = clfftInitSetupData(&fftSetup);   OCLfft_checkError(err, "clfftInitSetupData");
        err  = clfftSetup        (&fftSetup);   OCLfft_checkError(err, "clfftSetup" );
        //data_cl = clCreateBuffer( context, CL_MEM_READ_WRITE, buffer_size, NULL, &err );
        planFFT(  );                            //printf("initFFT 4 \n");
    }

    int newFFTbuffer( char* name, int nfloat=2, int ntot=-1 ){
        if(ntot<0)ntot=Ntot;
        return newBuffer( name, ntot, sizeof(float)*nfloat, 0, CL_MEM_READ_WRITE );
    }

    int newFFTimage( char* name, void* data=0, cl_int flags=CL_MEM_READ_ONLY ){
        if(data) flags |= CL_MEM_COPY_HOST_PTR;
        return newBufferImage3D( name, Ns[0], Ns[1], Ns[2], sizeof(float)*4, data, flags, {CL_RGBA, CL_FLOAT} );
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
        //printf( "DEBUG initBasisTable %i %i \n", nx, ny );
        itex_basis = newBufferImage2D( "BasisTable", ny, nx,   sizeof(float)*4,  data, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR , {CL_RGBA, CL_FLOAT} );
        //itex_basis = newBufferImage2D( "BasisTable", nx, ny,   sizeof(float),  data, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR , {CL_R, CL_FLOAT} ); // THIS WORKS FOR FLOAT TEXTURE
        //itex_basis = newBufferImage2D( "BasisTable", nx/4, ny,   sizeof(float)*4,  data, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR , {CL_RGBA, CL_FLOAT} );
        return itex_basis;
    }

    void runFFT( int ibuff, bool fwd, float* data=0 ){
        int err=0;
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
        int err=0;
        err =  clSetKernelArg(kernel, 0, sizeof(int),    &kdim.global        );
        err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &(buffers[0].p_gpu) );
        err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &(buffers[1].p_gpu) );
        err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &(buffers[2].p_gpu) );
        //checkError(err, "Setting kernel args");
        //double start_time = wtime();
        //printf( "mul_buffs kdim: dim %i global %li local %li \n", kdim.dim, kdim.global[0], kdim.local[0] ); 
        err = clEnqueueNDRangeKernel(  commands, kernel,   kdim.dim, NULL, kdim.global, kdim.local, 0, NULL, NULL);    
        OCL_checkError(err, "Enqueueing kernel");
        //err = clFinish(commands);
        //OCL_checkError(err, "Waiting for kernel to finish");
        //double run_time = wtime() - start_time;
    }

    void roll_buf( int ibuffA, int ibuffB, int4 shift ){
        int4 ngrid{ (int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3] };
        //printf( "DEBUG roll_buf iKernell_roll %i ibuffA %i ibuffB %i \n", iKernell_roll, ibuffA, ibuffB );
        useKernel( iKernell_roll );
        int err=0;
        err |= useArgBuff( ibuffA );
        err |= useArgBuff( ibuffB );
        err |= _useArg( shift );
        err |= _useArg( ngrid );
        OCL_checkError(err, "roll_bufs_1 ");
        printf( "DEBUG roll_buf 2 []\n" );
        //err = enque( 3, Ns, 0 ); 
        err = enque( 3, *(size_t4*)&Ns, (size_t4){1,1,1,1} );
        OCL_checkError(err, "roll_bufs_1 ");  
    }

    void gradient( int ibuffA, int ibuffB, float4 mask ){
        int4 ngrid{ (int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3] };
        //printf( "DEBUG roll_buf iKernell_roll %i ibuffA %i ibuffB %i \n", iKernell_roll, ibuffA, ibuffB );
        useKernel( iKernell_grad );
        int err=0;
        err |= useArgBuff( ibuffA );
        err |= useArgBuff( ibuffB );
        err |= _useArg( mask );
        err |= _useArg( ngrid );
        OCL_checkError(err, "gradient 1 ");
        err = enque( 3, *(size_t4*)&Ns, (size_t4){1,1,1,1} );
        OCL_checkError(err, "gradient 2 ");  
    }

    void projectAtomPosTex(  float4* atoms, float4* coefs, int nPos, float4* poss, float2* out ){
        KernelDims kdim;
        kdim.dim        = 1;
        kdim.global[0]  = nPos;
        kdim.local [0]  = 16;
        kdim.fitlocal( ); //printf( "projectAtomPosTex %li \n", kdim.global[0] );
        cl_kernel kernel = kernels[iKernell_projectPos_tex]; 
        //for(int i=0; i<nPos; i++){ printf("projectAtomPosTex %i (%g,%g,%g)\n", i, poss[i].x, poss[i].y, poss[i].z  );  };
        //printf("DEBUG projectAtomPosTex() 0 \n");
        int ibuff_poss = newBuffer( "poss", nPos, sizeof(float4), (float*)poss, CL_MEM_READ_WRITE );
        int ibuff_out  = newBuffer( "out" , nPos, sizeof(float2), (float*)out , CL_MEM_READ_WRITE );
        //printf("DEBUG projectAtomPosTex() 1 \n");
        buffers[ibuff_poss].toGPU(commands);
        upload(ibuff_atoms,atoms);
        upload(ibuff_coefs,coefs);
        int err=0;
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
        cltask_mul.setup( this, iKernell_mull, 1, {Ntot*2,0,0,0}, {16,0,0,0} );
        //cltask_mul.setup( this, iKernell_mull, 1, 8, 1 ); printf( "WARRNING : initTask_mul() IS WRONG !!!! %i %i \n", ibuffA, ibuffB );
        cltask_mul.args = { 
            INTarg (cltask_mul.global.x),
            BUFFarg(ibuffA),
            BUFFarg(ibuffB),
            BUFFarg(ibuff_result)         
        };
    }

    void initTask_lincomb( int ibuffA, int ibuffB, int ibuff_result ){
        Vec2f coefs;
        cltask_lincomb.setup( this, iKernell_lincomb, 1, {Ntot*2,0,0,0}, {16,0,0,0} );
        cltask_lincomb.args = { 
            INTarg (cltask_lincomb.global.x),
            BUFFarg(ibuffA),
            BUFFarg(ibuffB),
            BUFFarg(ibuff_result),
            REFarg(coefs)
        };
    }

    void initTask_poissonW( int ibuffA, int ibuff_result ){
        //printf( "BEGIN initTask_poissonW \n" );
        //cltask_poissonW.setup( this, iKernell_poissonW, 1, Ntot, 1 );
        cltask_poissonW.setup( this, iKernell_poissonW, 3, *(size_t4*)Ns, (size_t4){1,1,1,1} );
        cltask_poissonW.args = { 
            INTarg ((int)Ntot),
            BUFFarg(ibuffA),
            BUFFarg(ibuff_result),
            REFarg(dcell_poisson)           
        };
        //printf( "END initTask_poissonW \n" );
    }

    void initTask_gradient( int ibuffA, int ibuff_result ){
        //printf( "BEGIN initTask_gradient \n" );
        //cltask_poissonW.setup( this, iKernell_poissonW, 1, Ntot, 1 );
        cltask_gradient.setup( this, iKernell_gradient, 3, *(size_t4*)Ns, (size_t4){1,1,1,1} );
        Nvec  =(int4){(int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3]};
        cltask_gradient.args = { 
            REFarg(Nvec),           //5
            BUFFarg(ibuffA),
            BUFFarg(ibuff_result),
            REFarg (dcell_gradient)        
        };
        //printf( "END initTask_gradient \n" );
    }

    void initTask_project( int ibuffAtoms, int ibuffCoefs, int ibuff_result ){
        Nvec  =(int4){(int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3]};
        cltask_project.setup( this, iKernell_project, 1, {Ntot*2,0,0,0}, {16,0,0,0} );
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
        //printf("DEBUG initTask_project_tex() \n");
        Nvec  =(int4){(int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3]};
        cltask_project_tex.setup( this, iKernell_project_tex, 1, {Ntot*2,0,0,0}, {16,0,0,0} );
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
        //cltask_project_tex.print_arg_list();
    }

    void initTask_project_dens_tex( int ibuffAtoms, int ibuffCoefs, int ibuff_result ){
        printf("DEBUG initTask_project_dens_tex() \n");
        Nvec  =(int4){(int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3]};
        cltask_project_den_tex.setup( this, iKernell_project_dens_tex, 1, {Ntot*2,0,0,0}, {16,0,0,0} );
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
            REFarg(dC),           //12
            REFarg(acumCoef)      //13
        };
        //printf("DEBUG cltask_project_den_tex.args.size() %li  \n", cltask_project_den_tex.args.size() );
        //printf("DEBUG initTask_project_dens_tex() END \n");
    }

    void initTask_project_atom_dens_tex( int ibuffAtoms, int ibuffCoefs, int ibuff_result ){
        Nvec  =(int4){(int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3]};
        cltask_project_atom_dens_tex.setup( this, iKernell_project_atom_dens_tex, 1, {Ntot*2,0,0,0}, {16,0,0,0} );
        cltask_project_atom_dens_tex.args = { 
            INTarg (nAtoms),        //1
            BUFFarg(ibuffAtoms),    //4
            BUFFarg(ibuffCoefs),    //5
            BUFFarg(ibuff_result),  //6
            BUFFarg(itex_basis),    //7
            REFarg(Nvec),           //8
            REFarg(pos0),           //9
            REFarg(dA),           //10
            REFarg(dB),           //11
            REFarg(dC),           //12
            REFarg(acumCoef)      //13
        };
    }

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

    void projectAtomsDens( float4* atoms, float4* coefs, int ibuff_result, int iorb1, int iorb2, float2 acumCoef_ ){
        int ierr=0;
        printf( "projectAtomsDens() iorb=[%i .. %i] iresult=%i acumCoef(%g,%g) \n", iorb1, iorb2, ibuff_result, acumCoef_.x, acumCoef_.y  );
        //printf( "DEBUG projectAtomsDens acumCoef_ (%g,%g) \n", acumCoef_.x, acumCoef_.y  );
        //printf( "DEBUG projectAtomsDens(%i,%i) | atoms* %li long* %li \n", iorb1, iorb2,  (long)atoms, (long)coefs );
        if( atoms ) upload(ibuff_atoms,atoms); 
        if( coefs ) upload(ibuff_coefs,coefs);
        //clFinish(commands); 
        ierr = finishRaw();                               OCL_checkError( ierr, "upload");
        //initTask_project_dens_tex
        initTask_project_dens_tex( ibuff_atoms, ibuff_coefs, ibuff_result );   
        //cltask_project_den_tex
        cltask_project_den_tex.args[1].i=iorb1;
        cltask_project_den_tex.args[2].i=iorb2;
        acumCoef=acumCoef_;
        //cltask_project_den_tex.print_arg_list();
        ierr = cltask_project_den_tex.enque( );           OCL_checkError( ierr, "enque"  );
        finishRaw();                                      OCL_checkError( ierr, "finish" );
        //download( ibuff_result, (float2*)coefs );
        //printf( "DEBUG projectAtomsDens() END \n");
    }
    
    void projectAtomsDens0( int ibuff_result, float2 acumCoef_, int natoms_=0, int* ityps=0, Vec3d* oatoms=0 ){
        //printf( "DEBUG projectAtomsDens0 acumCoef_ (%g,%g) natoms_ %i \n", acumCoef_.x, acumCoef_.y, natoms_ );
        if(natoms_>0){
            makeAtomDensCoefs( natoms_, ityps, oatoms, false );
            finishRaw();
        }
        initTask_project_atom_dens_tex( ibuff_atoms, ibuff_coefs, ibuff_result );
        acumCoef=acumCoef_;
        cltask_project_atom_dens_tex.enque( );
        finishRaw();
        //printf( "DEBUG projectAtomsDens0() END \n");
    }


    void convolution( int ibuffA, int ibuffB, int ibuff_result ){
        //printf( "DEBUG convolution ibuffA,ibuffB,ibuff_result %i,%i,%i ", ibuffA,ibuffB,ibuff_result  );
        //saveToXsf( "DEBUG_buffA_in.xsf", ibuffA );
        //saveToXsf( "DEBUG_buffB_in.xsf", ibuffB );
        int err=0;
        err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffA].p_gpu, NULL, NULL);
        err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffB].p_gpu, NULL, NULL);  
        //saveToXsf( "DEBUG_buffA_w.xsf", ibuffA );
        //saveToXsf( "DEBUG_buffB_w.xsf", ibuffB );
        //mul_buffs( ibuffA, ibuffB, ibuff_result );
        initTask_mul( ibuffA, ibuffB, ibuff_result );  //cltask_mul.print_arg_list();
        cltask_mul.enque( );
        //saveToXsf( "DEBUG_buffOut_w.xsf", ibuff_result );
        err = clfftEnqueueTransform( planHandle, CLFFT_BACKWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuff_result].p_gpu, NULL, NULL);  
        //saveToXsf( "DEBUG_buffOut_out.xsf", ibuff_result );
    }

    void poisson( int ibuffA, int ibuff_result, float4* dcell=0 ){
        int err=0;
        //printf( "BEGIN poisson %i -> %i ( %s -> %s ) \n", ibuffA, ibuff_result, buffers[ibuffA].name, buffers[ibuff_result].name );
        //err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffA].p_gpu, NULL, NULL);
        //OCLfft_checkError(err, " poisson::clfftEnqueueTransform " );
        if( dcell ){ dcell_poisson = *dcell; }
        initTask_poissonW( ibuffA, ibuff_result );                   //printf( "DEBUG 1 ");
        //finishRaw(); saveToXsf( "rho.xsf", ibuffA );                 //printf( "DEBUG 2 ");
        err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffA].p_gpu, NULL, NULL);
        OCLfft_checkError(err, " poisson::clfftEnqueueTransform " ); //printf( "DEBUG 4 ");
        //finishRaw(); saveToXsf( "rho_w.xsf", ibuffA );               //printf( "DEBUG 5 ");
        cltask_poissonW.enque( );                                    //printf( "DEBUG 6 ");
        //finishRaw(); saveToXsf( "Vw.xsf", ibuff_result );            //printf( "DEBUG 7 ");
        err = clfftEnqueueTransform( planHandle, CLFFT_BACKWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuff_result].p_gpu, NULL, NULL);  
        //finishRaw(); saveToXsf( "V.xsf", ibuff_result );             //printf( "DEBUG 8 ");
        //printf( "END poisson \n" );
    }

    void gradient( int ibuffA, int ibuff_result, float4* dcell=0 ){
        //printf( "BEGIN gradient %i -> %i ( %s -> %s ) \n", ibuffA, ibuff_result, buffers[ibuffA].name.c_str(), buffers[ibuff_result].name.c_str() );
        if( dcell ){ dcell_gradient = *dcell; }
        initTask_gradient( ibuffA, ibuff_result );
        cltask_gradient.enque( );
        finishRaw();
    }

    void cleanup(){
        /*
        //clReleaseMemObject( data_cl );
        //free( data );
        err = clfftDestroyPlan( &planHandle );
        clfftTeardown( );
        //clReleaseCommandQueue( commands );
        //clReleaseContext( context );
        release_OCL();
        */
    }

    void makeMyKernels( const char*  cl_src_dir ){
        char srcpath[1024];
        sprintf( srcpath, "%s/myprog.cl", cl_src_dir );
        printf( "DEBUG makeKrenels() %s \n", srcpath );
        buildProgram( srcpath );  //printf( "DEBUG makeMyKernels 1 program %li \n", (long)program );

        iKernell_mull             = newKernel( "mul" );
        iKernell_roll             = newKernel( "roll" );
        iKernell_grad             = newKernel( "makeForceField" );
        iKernell_poissonW         = newKernel( "poissonW" );
        iKernell_gradient         = newKernel( "gradient" );
        iKernell_project          = newKernel( "projectAtomsToGrid" );
        iKernell_project_tex      = newKernel( "projectAtomsToGrid_texture"  );
        iKernell_project_dens_tex = newKernel( "projectOrbDenToGrid_texture" );
        iKernell_project_atom_dens_tex = newKernel( "projectAtomDenToGrid_texture" );
        iKernell_projectPos_tex   = newKernel( "projectWfAtPoints_tex" );
        //printf( "DEBUG makeMyKernels END \n" );
        //exit(0);
    };

    float* loadWfBasis( const char* path, float RcutSamp, int nsamp, int ntmp, int nZ, int* iZs, float* Rcuts, bool bDelete=true ){
        //printf( "loadWfBasis(%s) nsamp %i ntmp %i nZ %i RcutSamp %g [A] verbosity %i \n", path, nsamp, ntmp, nZ, RcutSamp, verbosity  );
        float* data_tmp = new float[ntmp      ];
        float* data     = new float[nsamp*2*nZ];
        char fname[64];
        //float dxTmp =(Rcut*const_Bohr_Radius)/ntmp;
        //float dxSamp=Rcut/nsamp;
        //RcutSamp*=0.529177210903f;
        if(verbosity>0)printf( "loadWfBasis(%s) nsamp %i ntmp %i nZ %i RcutSamp %g [A]\n", path, nsamp, ntmp, nZ, RcutSamp );
        for(int i=0; i<nZ; i++){
            int iz=iZs[i];
            float Ri = Rcuts[i];
            // --- wf1
            sprintf( fname, "%s%03i_%03i.wf%i", path, iz, (int)(Ri*100), 1 );
            int nin = loadWf_(fname, data_tmp );
            //resample1D( nsamp, 0, 0, dxSamp, dxTmp, data_tmp, data+nsamp*(i*2), 2,0 );
            resample1D( nin, 0.0, Ri*const_Bohr_Radius, data_tmp,   nsamp,0.0, RcutSamp, data+nsamp*(i*2),   2,0 );
            // --- wf2
            sprintf( fname, "%s%03i_%03i.wf%i", path, iz, (int)(Ri*100), 2 );
            if(verbosity>0)printf( "loadWfBasis[%i] (iZ=%2i,nin=%3i,Ri=%5.3f) %s \n", i, iz, Ri, nin, fname );
            if( loadWf_(fname, data_tmp ) ){
                printf( "resample1D \n" );
                resample1D( nin, 0.0, Ri*const_Bohr_Radius, data_tmp,   nsamp,0.0, RcutSamp, data+nsamp*(i*2),   2,1 );
            }else{
                printf( "copy from wf1 nsamp=%i \n", nsamp );
                for(int j=0; j<nsamp; j++){ int j0=nsamp*(i*2)+2*j; data[j0+1]=data[j0]; }
            }
        }
        //for(int i=0; i<nsamp; i++){ printf("basis [%i] (%f,%f) (%f,%f) \n", i, data[i*2],data[i*2+1],data[i*2+2*nsamp],data[i*2+1+2*nsamp] );  }
        delete [] data_tmp;
        //for(int i=0; i<nZ; i++){  printf( "wf[%2i]:", i );  for(int j=0; j<10; j++){ printf( "%g ", data[i*nsamp+j] ); }; printf( "\n" ); };
        itex_basis = newBufferImage2D( "BasisTable", nsamp, nZ,  sizeof(float)*2,  data, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR , {CL_RG, CL_FLOAT} );
        if(bDelete){ delete [] data; data=0; }
        return data;
    }

    void update_GridShape(){
        //printf("update_GridShape() %g %g %g \n", dC.x, dC.y, dC.z );
        grid.cell.a=(Vec3d)(*(Vec3f*)&dA)*Ns[0];
        grid.cell.b=(Vec3d)(*(Vec3f*)&dB)*Ns[1];
        grid.cell.c=(Vec3d)(*(Vec3f*)&dC)*Ns[2];
        grid.pos0  =(Vec3d)(*(Vec3f*)&pos0);
        grid.n = (Vec3i){(int)Ns[0],(int)Ns[1],(int)Ns[2]}; //Ntot;
        grid.updateCell();
        //grid.printCell();
    }

    void saveToXsfData(const char* fname, Vec3i ngrid, double* data, int natoms=0, int* atypes=0, Vec3d* apos=0 ){
        Ns[0]=ngrid.x; Ns[1]=ngrid.y; Ns[2]=ngrid.z;
        update_GridShape();
        grid.saveXSF( fname, data, 1, 0, natoms, atypes,apos );
    }

    void saveToXsf(const char* fname, int ibuff, int stride=2, int offset=0, int natoms=0, int* atypes=0, Vec3d* apos=0 ){
        update_GridShape();                   
        float* cpu_data = new float[Ntot*stride]; // complex 2*float
        download( ibuff,cpu_data);
        finishRaw();
        float amin,amax; VecN::bounds<float>( Ntot*stride, cpu_data, amin, amax ); printf( "saveToXsf() amin,amax %g %g \n" ); // DEBUG
        grid.saveXSF( fname, cpu_data, stride, offset,   natoms, atypes,apos );
        delete [] cpu_data;
    }

    void saveToBin(const char* fname, int ibuff){
        if(verbosity>0)printf( "saveToBin( %i, %s ) \n", ibuff, fname );
        //update_GridShape();
        float* cpu_data = new float[Ntot*2]; // complex 2*float
        download( ibuff,cpu_data);
        finishRaw();
        saveBin( fname, (Ntot*2)*sizeof(float), (char*)cpu_data );
        delete [] cpu_data;
    }

    void loadFromBin(const char* fname, int ibuff){
        if(verbosity>0)printf( "loadFromBin( %i, %s ) \n", ibuff, fname );
        //update_GridShape();
        float* cpu_data = new float[Ntot*2]; // complex 2*float
        loadBin( fname, (Ntot*2)*sizeof(float), (char*)cpu_data );
        upload( ibuff,cpu_data);
        finishRaw();
        delete [] cpu_data;
    }

    void prepareAtomCoords(  int natoms, int* ityps, Vec3d* oatoms ){
        float4* atoms = new float4[ natoms ];
        for(int ia=0; ia<natoms; ia++){
            atoms[ia]=(float4){ (float)oatoms[ia].x,(float)oatoms[ia].y,(float)oatoms[ia].z, ityps[ia]+0.1f };
        }
        upload(ibuff_atoms,atoms, natoms );
        finishRaw();
        delete [] atoms;
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

    void assignAtomDensCoefs( int natoms, int* ityps, float4* coefs ){
        //printf( "DEBUG assignAtomDensCoefs() \n" );
        for(int ia=0; ia<natoms; ia++){
            int ityp     = ityps[ia]-1;
            float Qs = (float)(atype_Qconfs[ityp].x);
            float Qp = (float)(atype_Qconfs[ityp].y/3.0);
            coefs[ia]=(float4){ Qp,Qp,Qp, Qs };
            //printf( "atom[%i] itip %i Qcoefs (%g,%g,%g,%g)\n", ia,ityp, coefs[ia].x,coefs[ia].y,coefs[ia].z,coefs[ia].w );
        }
        //for(int ia=0; ia<natoms; ia++){  printf(  "AtomQs[%i](%g|%g,%g,%g)\n", ia, coefs[ia].w, coefs[ia].x,coefs[ia].y,coefs[ia].z );  }
        //printf( "DEBUG assignAtomDensCoefs() DONE\n" );
    }

    void makeAtomDensCoefs( int natoms, int* ityps, Vec3d* oatoms, bool bInit=false ){
        //printf( "DEBUG makeAtomDensCoefs() \n" );
        if( bInit ){ initAtoms( natoms, natoms ); }
        prepareAtomCoords( natoms, ityps, oatoms );
        float4* coefs = new float4[ natoms  ];
        assignAtomDensCoefs( natoms, ityps, coefs );
        upload(ibuff_coefs,coefs, natoms  );
        delete [] coefs;
        //printf( "DEBUG makeAtomDensCoefs() DONE \n" );
    };

    int assignDiagonalOrbCoefs( int norb, int natoms, int* ityps, float4* coefs ){
        //printf( "DEBUG assignDiagonalOrbCoefs()\n" );
        int ia_orb=0;
        int io=0;
        int orbCount=0;
        for(int iorb=0; iorb<norb; iorb++){
            float4* cs = coefs+iorb*natoms;
            //printf( "iorb[%i] \n", iorb );
            int ityp = ityps[ia_orb];
            int nOrbAtom = atype_nOrb[ityp];
            if(io>=nOrbAtom){ io=0; ia_orb++; if(ia_orb>=natoms) break; }
            //printf( "iorb[%i|%i,%i] norb %i ityp %i \n", iorb, ia_orb,io,  norb, ityp );
            for(int ia=0; ia<natoms; ia++){ cs[ia]=(float4){0.f,0.f,0.f,0.f}; };
            if(io>0){ 
                double Q = atype_Qconfs[ityp].y/3.0; 
                //printf( "p Q %g \n", Q );
                ((float*)(cs+ia_orb))[io-1] = sqrt(Q);
            }else{
                double Q = atype_Qconfs[ityp].x;
                //printf( "s Q %g \n", Q ); 
                cs[ia_orb].w = sqrt(Q);
            }
            io++;
            orbCount++;
        }
        //printf( "DEBUG assignDiagonalOrbCoefs() PRINT OUT orbCount %i \n", orbCount );
        //for(int iorb=0; iorb<norb; iorb++){
        //    printf( "DEBUG ORB[%i]\n", iorb );
        //    float4* cs = coefs+iorb*natoms;
        //    for(int ia=0; ia<natoms; ia++){ 
        //        printf(  "(%g|%g,%g,%g)\n", cs[ia].w, cs[ia].x,cs[ia].y,cs[ia].z ); 
        //    }
        //}
        //printf( "DEBUG assignDiagonalOrbCoefs() DONE\n" );
        return orbCount;
    }


    int convCoefs( int natoms, int* iZs, int* ityps, double* ocoefs, double* oatoms, bool bInit=false, bool bDiagonal=false ){
        //printf( "DEBUG convCoefs() ocoefs %li \n", (long)ocoefs );
        int norb=0;
        for(int ia=0; ia<natoms; ia++){ if(iZs[ia]==1){ norb+=1; }else{ norb+=4; }; }   // --- Count orbitals
        int ncoef=natoms*norb;
        float4* coefs = new float4[ ncoef  ];
        if(bDiagonal){
            int norb_ = assignDiagonalOrbCoefs( norb, natoms, ityps, coefs );
            //printf( " convCoefs(): norb %i norb_ %i ", norb, norb_ );
            norb=norb_;
        }else for(int iorb=0; iorb<norb; iorb++){
            convOrbCoefs( natoms, iZs, ocoefs+iorb*norb, coefs+iorb*natoms );
        }
        if( bInit ){
            initAtoms( natoms, norb );
        }
        prepareAtomCoords( natoms, ityps, (Vec3d*)oatoms );
        upload(ibuff_coefs,coefs, ncoef  );
        delete [] coefs;
        return norb;
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
                if(iZs[ia]==1){ lcoefs[ia]=Quat4f{0.f,0.f,0.f, (float)ocoefs[io] };}
                else          { lcoefs[ia]=Quat4f{ (float)ocoefs[io+3],(float)ocoefs[io+1],(float)ocoefs[io+2], (float)ocoefs[io] }; }
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

        Vec3d lbox  = Vec3d{3.0,3.0,3.0};
        Vec3d cell  = Vec3d{ grid.cell.a.x, grid.cell.b.y, grid.cell.c.z };
        Vec3i nbox  = Vec3i{ (int)(1+cell.x/lbox.x), (int)(1+cell.y/lbox.y), (int)(1+cell.z/lbox.z)  };
        Vec3d dcell = Vec3d{ cell.x/nbox.x, cell.y/nbox.y, cell.z/nbox.z, };
        for(int ix=0; ix<nbox.x; ix++){
            for(int iy=0; iy<nbox.y; iy++){
                for(int iz=0; iz<nbox.z; iz++){
                    Vec3d p0 = grid.pos0 + dcell*(Vec3d{(double)ix,(double)iy,(double)iz});
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

#endif
