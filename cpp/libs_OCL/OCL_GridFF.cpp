

// No need to explicitely include the OpenCL headers 
//#include <clFFT.h>

#define CL_TARGET_OPENCL_VERSION 200

#include <CL/cl.h>
#include <clFFT.h>
#include "OCLfft_errors.h"
#include <clFFT.h>

#include "testUtils.h"
#include "OCL.h"
#include "Grid.h"
#include "IO_utils.h"
#include "quaternion.h"

int verbosity = 0;
#include "OCL_DFT.h"
#include "OCL_PP.h"

OCL_PP oclfft;

//#include "FireCoreAPI.h"
//#include "approximation.h"

//OCLsystem ocl;
//Approx::AutoApprox aaprox;
//FireCore::Lib      fireCore;

extern "C" {

    void setVerbosity(int verbosity_ ){ verbosity=verbosity_; }

    void setErrorCheck(int ierr){ bOCLCheckError=ierr>0; }

    void  init( const char* cl_src_dir ){
        oclfft.init();
        oclfft.makeMyKernels( cl_src_dir );
    }
    void printDeviceInfo( bool bDetails ){ oclfft.printDeviceInfo( bDetails ); }

    int   upload(int i, const float* cpu_data ){ return oclfft  .upload(i,cpu_data);                                   };
    int download(int i,       float* cpu_data ){ int iret=oclfft.download(i,cpu_data);  oclfft.finishRaw();  OCL_checkError_(iret,"OCL_GridFF.cpp::download()",i);  return iret; };

    int copy           ( int iBufFrom, int iBufTo, int nbytes, int  src_offset, int  dst_offset ){ int iret=oclfft.copy           ( iBufFrom, iBufTo, nbytes, src_offset, dst_offset );            OCL_checkError_(iret,"OCL_GridFF.cpp::copy()",iBufFrom);          return iret; };
    int copyBuffToImage( int iBuff, int itex, int nx,int ny,int nz                              ){ int iret=oclfft.copyBuffToImage( iBuff,      itex, size_t4{(size_t)nx,(size_t)ny,(size_t)nz} ); OCL_checkError_(iret,"OCL_GridFF.cpp::copyBuffToImage()",iBuff);  return iret; };
 
    void roll_buf( int ibuffA, int ibuffB, int* shift ){ return oclfft.roll_buf( ibuffA, ibuffB, *(int4*)shift ); }

    int   upload_d(int ibuf, const double* data, bool bComplex ){ 
        int n=oclfft.buffers[ibuf].n; 
        //printf( "DEBUG upload_d ibuf %i bComplex %i  n %i \n", ibuf, bComplex, n );
        float2* tmp = new float2[ n ];
        if  (bComplex){ for(int i=0; i<n; i++){ tmp[i]=(float2){(float)data[i*2],(float)data[i*2+1]};  } }
        else          { for(int i=0; i<n; i++){ tmp[i]=(float2){(float)data[i  ], 0.0f             };  } }
        //int nxy = oclfft.Ns[0]*oclfft.Ns[1];
        //if( (i%nxy)==0 ) printf( "CPU iz %i i %i data %g A(%g,%g) \n", i/nxy, i, data[i], tmp[i].x, tmp[i].y );
        return oclfft.upload(ibuf,tmp); 
        delete [] tmp;
    }

    void initFFT( int ndim, size_t* Ns_ ){
        oclfft.initFFT( ndim, Ns_ );      //printf( "C initFFT 1 \n" );
        oclfft.newFFTbuffer( "inputA" );  //printf( "C initFFT 2 \n" );
        oclfft.newFFTbuffer( "inputB" );  //printf( "C initFFT 3 \n" );
        oclfft.newFFTbuffer( "outputC" ); //printf( "C initFFT 4 \n" );
        //oclfft.initTask_mul( 0, 1, 2 );    // If we know arguments in front, we may define it right now
    }


    void release( bool bReleaseOCL, bool bReleaseOCLfft ){ oclfft.release_OCL_DFT( bReleaseOCL, bReleaseOCLfft); }

    // ================ PP

    int initPP( const char* cl_src_dir, size_t* Ns_ ){
        oclfft.init();
        oclfft.makeMyKernels ( cl_src_dir );
        oclfft.makeKrenels_PP( cl_src_dir );
        oclfft.initFFT( 3, Ns_ );
        oclfft.itex_FF = oclfft.newFFTimage( "FF" );  // make writable texture (save some memory by not requiring temporary buffer)
        oclfft.itex_FF = oclfft.newFFTimage( "FF", 0, CL_MEM_READ_WRITE );
        //float4* data = oclfft.debug_gen_FE();      ;printf("C DEBUG 5 \n");
        //oclfft.itex_FF = oclfft.newFFTimage( "FF", data ); ;printf("C DEBUG 6 \n");
        //delete [] data; ;printf("C DEBUG 7 \n");
        //oclfft.newFFTbuffer( "inputA" );  //printf( "C initFFT 2 \n" );
        //oclfft.newFFTbuffer( "inputB" );  //printf( "C initFFT 3 \n" );
        //oclfft.newFFTbuffer( "outputC" ); //printf( "C initFFT 4 \n" );
        return oclfft.itex_FF;
    }

    void makeStartPointGrid( int nx, int ny, double* p0, double* da, double* db ){ oclfft.makeStartPointGrid( (Vec2i){nx,ny}, *(Vec3d*)p0, *(Vec3d*)da, *(Vec3d*)db ); }
    void setGridShapePP    ( double* p0, double* dCell                          ){ if(p0)v2f4(*(Vec3d*)p0,oclfft.pos0); oclfft.setGridShape( *(Mat3d*)dCell ); }

    void relaxStrokesTilted( int ibuff_out, int nz, float dtip, int np=0, float* points=0 ){ oclfft.relaxStrokesTilted( ibuff_out, nz, dtip, np, (float4*)points ); };
    void getFEinStrokes    ( int ibuff_out, int nz, double* dTip, int np=0, float* points=0 ){ oclfft.getFEinStrokes    ( ibuff_out, nz, *(Vec3d*)dTip, np, (float4*)points ); };
    void evalLJC_QZs      ( int ibuff_out, int na=0, float* atoms=0, float* coefs=0 ){ oclfft.evalLJC_QZs( ibuff_out, na, (float4*)atoms, (float4*)coefs ); }
    void evalLJC_QZs_toImg(                int na=0, float* atoms=0, float* coefs=0 ){ oclfft.evalLJC_QZs_toImg(      na, (float4*)atoms, (float4*)coefs ); }


    // ================ END PP

    void newFFTbuffer( char* name, int nfloat, int ntot ){ oclfft.newFFTbuffer( name, nfloat, ntot ); }

    int initAtoms( int nAtoms, int nOrbs ){  return oclfft.initAtoms( nAtoms, nOrbs ); };
    void runfft( int ibuff, bool fwd     ){ oclfft.runFFT( ibuff,fwd,0);     };
    //void runfft( int ibuff, bool fwd, float* data ){ oclfft.runFFT( ibuff,fwd, data); };
    void convolve( int ibuffA, int ibuffB, int ibuff_result ){  oclfft.convolution( ibuffA, ibuffB, ibuff_result  );}
    void poisson ( int ibuffA, int ibuff_result, float* dcell ){  oclfft.poisson ( ibuffA, ibuff_result, (float4*)dcell );}
    void gradient( int ibuffA, int ibuff_result, float* dcell ){  oclfft.gradient( ibuffA, ibuff_result, (float4*)dcell );}
    void projectAtoms    ( float* atoms, float* coefs, int ibuff_result                       ){ oclfft.projectAtoms    ( (float4*)atoms, (float4*)coefs, ibuff_result ); }
    void projectAtomsDens( float* atoms, float* coefs, int ibuff_result, int iorb1, int iorb2, float* acumCoef ){  oclfft.projectAtomsDens( (float4*)atoms, (float4*)coefs, ibuff_result, iorb1, iorb2, *(float2*)acumCoef ); }
    void projectAtomsDens0( int ibuff_result, float* acumCoef, int natoms=0, int* ityps=0, Vec3d* oatoms=0 ){ oclfft.projectAtomsDens0( ibuff_result, *(float2*)acumCoef, natoms, ityps, (Vec3d*)oatoms ); }

    void projectDenmat( float* atoms, float* coefs, int ibuff_result, int iorb1, int iorb2, float* acumCoef ){  oclfft.projectDenmat( (float4*)atoms, (float4*)coefs, ibuff_result, iorb1, iorb2, *(float2*)acumCoef ); }
    
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

    void setTypes( int nAtype, int* atype_nOrb_, float* atype_Qconfs_, bool bInternal ){
        oclfft.nAtype=nAtype;
        if(bInternal){
            _realloc(oclfft.atype_nOrb,   nAtype );
            _realloc(oclfft.atype_Qconfs, nAtype );
            for(int i=0; i<nAtype;i++){
                oclfft.atype_nOrb[i]   =           atype_nOrb_   [i];
                oclfft.atype_Qconfs[i] = ((float2*)atype_Qconfs_)[i];
                //printf( "setTypes()[%i] atype_nOrb %i atype_Qconfs (%g,%g) \n", i, oclfft.atype_nOrb[i], oclfft.atype_Qconfs[i].x, oclfft.atype_Qconfs[i].y );
            }
        }else{

            oclfft.atype_nOrb  =         atype_nOrb_; // number of orbitals per atomic type (hydrogen=1(s), carbon=4(s,px,py,pz))
            oclfft.atype_Qconfs=(float2*)atype_Qconfs_;
        }
    }

    int initBasisTable( int nx, int ny, float* data ){  return oclfft.initBasisTable(nx,ny,data ); };

    int convCoefs( int natoms, int* iZs, int* ityps, double* ocoefs, double* oatoms, bool bInit, bool bDiagonal ){  return oclfft.convCoefs( natoms, iZs, ityps, ocoefs, oatoms, bInit, bDiagonal ); }

    void loadWf(const char* fname, float* out){ loadWf_(fname, out); };

    float* loadWfBasis( const char* path, float RcutSamp, int nsamp, int ntmp, int nZ, int* iZs, float* Rcuts, bool bDelete ){ return oclfft.loadWfBasis(path, RcutSamp,nsamp,ntmp,nZ,iZs,Rcuts, bDelete ); }

    void saveToBin(const char* fname, int ibuff){ oclfft.saveToBin(fname, ibuff); }
    void loadFromBin(const char* fname, int ibuff){ oclfft.loadFromBin(fname,ibuff); }

    void saveToXsf     (const char* fname, int ibuff, int stride, int offset ){ return oclfft.saveToXsf(fname, ibuff,stride,offset,0,0,0); }
    void saveToXsfAtoms(const char* fname, int ibuff, int stride, int offset, int natoms, int* atypes, double* apos ){ return oclfft.saveToXsf(fname, ibuff, stride, offset, natoms,atypes,(Vec3d*)apos); }
    void saveToXsfAtomsData(const char* fname, int* ngrid, double* data, int natoms, int* atypes, double* apos ){ return oclfft.saveToXsfData(fname, *(Vec3i*)ngrid, data, natoms,atypes,(Vec3d*)apos); }


    /*

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
        //for(int i=0; i<64; i++){ printf( "pwfcoef[%i] %g \n", i, pwfcoef[i] ); };
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
            //printf( "coefs[%i] t %i | %g %g %g p|s %g \n", i, atypes[i],  coefs[i].x, coefs[i].y, coefs[i].z, coefs[i].w );
        }
        projectAtoms( (float*)apos_, (float*)coefs, 0 );
        //oclfft.saveToXsf( "test.xsf", 0 );
        oclfft.update_GridShape();
        float* cpu_data = new float[oclfft.Ntot*2]; // complex 2*float
        oclfft.download( 0, cpu_data);
        oclfft.finishRaw();

        int i0 = oclfft.Ns[0]*oclfft.Ns[1]*oclfft.Ns[3]/2 + oclfft.Ns[0]*oclfft.Ns[1]/2;
        //printf( "Ntot %li i0 %i Ns (%li,%li,%li) \n", oclfft.Ntot*2, i0*2, oclfft.Ns[0],oclfft.Ns[1],oclfft.Ns[2]  );

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
    */

};