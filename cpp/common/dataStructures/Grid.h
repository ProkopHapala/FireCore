
#ifndef Grid_h
#define Grid_h

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "Vec3.h"
#include "Mat3.h"
//#include <string.h>
#include "quaternion.h"
#include "VecN.h"

#include <vector>


#include "Fourier.h"

// ================= MACROS

#define fast_floor_offset  1000
#define fast_floor( x )    ( ( (int)( (x) + fast_floor_offset ) ) - fast_floor_offset )
#define i3D( ix, iy, iz )  ( (iz)*nxy + (iy)*nx + (ix)  )

// ================= CONSTANTS


int nearPow2(int i){ return ceil( log(i)/log(2) ); }

// ================================
// ======== class GridShape
// ================================

// Force-Field namespace
class GridShape { public:
	Vec3d   pos0;
	Mat3d   cell;       // lattice vector
	Mat3d   dCell;      // basis vector of each voxel ( lattice vectors divided by number of points )
	Mat3d   diCell;     // inversion of voxel basis vector
	Vec3i   n;          // number of pixels along each basis vector
    //bool bCellSet=false;

    void center_cell( Vec3d c ){ cell.dot_to_T( c, pos0 ); }

	//inline Vec3d * allocateArray_Vec3d(){ return new Vec3d[n.x*n.y*n.z); }
	inline int getNtot() const {return n.x*n.y*n.z ; }

	inline double voxelVolume()const{ return dCell.determinant(); }

    inline int ip2i(const Vec3i& ip){ return ip.a + ( n.a*( ip.b + n.b*ip.c) );  }

    inline void autoN( double step){
        n.a=(int)(cell.a.norm()/step)+1;
        n.b=(int)(cell.b.norm()/step)+1;
        n.c=(int)(cell.c.norm()/step)+1;
    }

	inline void updateCell( double step=-1.0 ){
        if(step>0){ autoN(step); }
        dCell.a.set_mul( cell.a, 1.0/n.a );
		dCell.b.set_mul( cell.b, 1.0/n.b );
		dCell.c.set_mul( cell.c, 1.0/n.c );
		dCell.invert_T_to( diCell );
	}

    inline void updateCell_2(){
        cell.a.set_mul( dCell.a, n.a );
		cell.b.set_mul( dCell.b, n.b );
		cell.c.set_mul( dCell.c, n.c );
		dCell.invert_T_to( diCell );
	}

	inline void setCell( const Mat3d& cell_ ){
		//n.set( n_ );
		cell.set( cell_ );
        updateCell();
	}

    int cut1D( Vec3i iv0,  Vec3i di, double* data, double* line ){
        Vec3i iv = iv0;
        int i=0;
        //printf( "iv0(%i,%i,%i) di(%i,%i,%i) n(%i,%i,%i) \n", iv.a,iv.b,iv.c,  di.a,di.b,di.c,    n.a,n.b,n.c );
        while( iv.isBetween( Vec3iZero, n) ){
            line[i] = data[ ip2i(iv) ];
            //printf( "[%i] (%i,%i,%i) -> %g \n", i, iv.a,iv.b,iv.c, line[i] );
            iv.add(di);
            i++;
        }
        return i;
    }

	int init(double R, double step, bool bPow2=false){
        cell = Mat3d{ (2*R),0.0,0.0,  0.0,(2*R),0.0,  0.0,0.0,(2*R) };
        int ngx   = (2*R)/step;
        int mpow = -1;
        if(bPow2){
            printf( "n %i ", ngx );
            mpow = nearPow2( ngx );
            ngx  = 1<<mpow;
            printf( "-> n %i | mpow %i \n", ngx, mpow );
        }
        n    = {ngx,ngx,ngx};
        pos0 = Vec3d{-R,-R,-R};
        updateCell();
        return mpow;
	}

	//inline void set( int * n_, double * cell_ ){ set( *(Vec3d*) n_, *(Mat3d*)cell_ ); };

	inline void grid2cartesian( const Vec3d& gpos, Vec3d& cpos ) const {
		cpos.set_mul( dCell.a, gpos.x );
		cpos.add_mul( dCell.b, gpos.y );
		cpos.add_mul( dCell.c, gpos.z );
		cpos.add(pos0);
	}

	inline void cartesian2grid( Vec3d cpos, Vec3d& gpos ) const {
        cpos.sub( pos0 );
		gpos.a = cpos.dot( diCell.a );
		gpos.b = cpos.dot( diCell.b );
		gpos.c = cpos.dot( diCell.c );
	}

    int loadCell(const char * fname, double step=-1.0 ){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        fscanf( pFile, "%lf %lf %lf", &cell.a.x, &cell.a.y, &cell.a.z );
        fscanf( pFile, "%lf %lf %lf", &cell.b.x, &cell.b.y, &cell.b.z );
        fscanf( pFile, "%lf %lf %lf", &cell.c.x, &cell.c.y, &cell.c.z );
        fclose(pFile);
        updateCell(step);
        return 0;
    }

	void printCell() const {
	    printf( " n      %i %i %i \n", n.x,        n.y,       n.z        );
        printf( " pos0   %f %f %f \n", pos0.x,     pos0.y,     pos0.z    );
	    printf( " a      %f %f %f \n", cell.a.x,   cell.a.y,   cell.a.z  );
	    printf( " b      %f %f %f \n", cell.b.x,   cell.b.y,   cell.b.z  );
	    printf( " c      %f %f %f \n", cell.c.x,   cell.c.y,   cell.c.z  );
	    printf( " da     %f %f %f \n", dCell.a.x,  dCell.a.y,  dCell.a.z  );
	    printf( " db     %f %f %f \n", dCell.b.x,  dCell.b.y,  dCell.b.z  );
	    printf( " dc     %f %f %f \n", dCell.c.x,  dCell.c.y,  dCell.c.z  );
	    printf( " inv_da %f %f %f \n", diCell.a.x, diCell.a.y, diCell.a.z );
	    printf( " inv_db %f %f %f \n", diCell.b.x, diCell.b.y, diCell.b.z );
	    printf( " inv_dc %f %f %f \n", diCell.c.x, diCell.c.y, diCell.c.z );
    }

    void headerToXsf( FILE* fout )const{
        fprintf( fout, "BEGIN_BLOCK_DATAGRID_3D\n" );
        fprintf( fout, "   some_datagrid\n" );
        //fprintf( fout, "   BEGIN_DATAGRID_3D_whatever\n" );
        fprintf( fout, "   DATAGRID_3D\n" );
        fprintf( fout, "%i %i %i\n", n.x, n.y, n.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", pos0.x,   pos0.y,   pos0.z   );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.a.x, cell.a.y, cell.a.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.b.x, cell.b.y, cell.b.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.c.x, cell.c.y, cell.c.z );
    }

    void primcoordToXSF( FILE* fout,  int natoms, int* atyps, Vec3d* apos )const{
        //printf("primcoordToXSF() 1 natoms %i apos %li atyps %li \n", natoms, (long)apos, (long)atyps );
        fprintf( fout, "CRYSTAL\n" );
        fprintf( fout, "PRIMVEC\n" );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.a.x, cell.a.y, cell.a.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.b.x, cell.b.y, cell.b.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.c.x, cell.c.y, cell.c.z );
        fprintf( fout, "PRIMCOORD\n" );
        fprintf( fout, "    %i     1 \n", natoms );
        for(int i=0; i<natoms; i++){
            Vec3d p = apos[i] - pos0;
            fprintf( fout, "%3i %9.6f %9.6f %9.6f \n", atyps[i], p.x,p.y,p.z );
        }
        fprintf( fout, "\n" );
    }

    void atomsToXSF( FILE* fout,  int natoms, int* atyps, Vec3d* apos )const{
        fprintf( fout, "ATOMS\n" );
        for(int i=0; i<natoms; i++){
            //Vec3d p = apos[i] - pos0;
            Vec3d p = apos[i];
            fprintf( fout, "%3i %9.6f %9.6f %9.6f \n", atyps[i], p.x,p.y,p.z );
        }
        fprintf( fout, "\n" );
    }

    template<typename T>
    void toXSF( FILE* fout, const T* FF, int pitch, int offset ) const {
        //printf( "DEBUG GridShale::toXSF() stride %i offset %i \n", pitch, offset );
        headerToXsf( fout );
        int nx  = n.x; 	int ny  = n.y; 	int nz  = n.z; int nxy = ny * nx;
        for ( int ic=0; ic<nz; ic++ ){
            //printf("toXSF[%i] pitch,offset(%i,%i)\n", ic, pitch, offset );
            for ( int ib=0; ib<ny; ib++ ){
                for ( int ia=0; ia<nx; ia++ ){
                   int i = i3D( ia, ib, ic );
                   fprintf( fout, "%6.5e\n", FF[i*pitch+offset] );
                }
            }
        }
        fprintf( fout, "   END_DATAGRID_3D\n" );
        fprintf( fout, "END_BLOCK_DATAGRID_3D\n" );
    }

    template<typename T>
    void saveXSF( const char * fname,const T* FF, int pitch=1, int offset=0, int natoms=0, int* atypes=0, Vec3d* apos=0, bool bPrimCoord=true )const {
        printf( "saving %s\n", fname );
        FILE *fout;
        fout = fopen(fname,"w");
        if( fout==0 ){ printf( "ERROR saveXSF(%s) : Cannot open file for writing \n", fname ); return; }
        if(natoms>0){
            if (bPrimCoord){ primcoordToXSF( fout,  natoms, atypes, apos );  }
            else           { atomsToXSF    ( fout,  natoms, atypes, apos );  }
        } 
        toXSF( fout, FF, pitch, offset );
        fclose(fout);
    }

    double Laplace( const double* f, double* out )const{
        int nx  = n.x; 	int ny  = n.y; 	int nz  = n.z; int nxy = ny * nx;
        double idx2 = diCell.a.norm2();
        double idy2 = diCell.b.norm2();
        double idz2 = diCell.c.norm2();
        double Ek = 0;
        for ( int ic=0; ic<nz; ic++ ){
            for ( int ib=0; ib<ny; ib++ ){
                for ( int ia=0; ia<nx; ia++ ){
                    //printf(" %i %i %i \n", ic, ib, ia );
                    if((ia==0)||(ia==(nx-1))||(ib==0)||(ib==(ny-1))||(ic==0)||(ic==(nz-1))){ if(out)out[i3D(ia,ib,ic)]=0.0; continue; }
                    double f00 = f[i3D(ia,ib,ic)]*2;
                    double ddf =
                        (f[i3D(ia+1,ib,ic)]+f[i3D(ia-1,ib,ic)]-f00)*idx2
                      + (f[i3D(ia,ib+1,ic)]+f[i3D(ia,ib-1,ic)]-f00)*idy2
                      + (f[i3D(ia,ib,ic+1)]+f[i3D(ia,ib,ic-1)]-f00)*idz2;
                    Ek += ddf*f00;
                    if(out)out[i3D(ia,ib,ic)] = ddf;
                }
            }
        }
        return Ek;
    }

    double integrate( const double* f )const{
        int nx  = n.x; 	int ny  = n.y; 	int nz  = n.z; int nxy = ny * nx;
        double sum = 0;
        //int itot=0;
        //double vmax = -1e+300;
        for ( int ic=0; ic<nz; ic++ ){
           for ( int ib=0; ib<ny; ib++ ){
                for ( int ia=0; ia<nx; ia++ ){
                   int i = i3D( ia, ib, ic );
                   double fi = f[i];
                   //vmax=fmax(vmax,fi);
                   sum += fi;
                   //itot++;
                }
            }
        }
        //printf( "DEBUG itot %i(%i) fmax %g sum %g \n", itot, n.totprod(), vmax, sum*voxelVolume() );
        return sum * voxelVolume();
    }

    double integrateProductShifted( Vec3i ishift, const double* f1, const double* f2 )const{
		double sum = 0;
        int nx  = n.x; 	int ny  = n.y; 	int nz  = n.z; int nxy = ny * nx;
        for ( int ic=0; ic<nz-ishift.z; ic++ ){
           for ( int ib=0; ib<ny-ishift.y; ib++ ){
                for ( int ia=0; ia<nx-ishift.x; ia++ ){
                   int i1 = i3D( ia         , ib         , ic          );
                   int i2 = i3D( ia+ishift.x, ib+ishift.y, ic+ishift.z );
                   sum += f1[i1] * f2[i2];
                }
            }
            //printf( "[%u] sum %g \n,", ic, sum );
        }
        return sum * voxelVolume();
    }

    double integrateProductShifted( Vec3d shift, const double* f1, const double* f2 )const{
        Vec3i ishift;
		ishift.a = round(shift.dot( diCell.a )); // Should we use (int) instead ?
		ishift.b = round(shift.dot( diCell.b ));
		ishift.c = round(shift.dot( diCell.c ));
        return integrateProductShifted( ishift, f1, f2 ) * voxelVolume();
    }

    void gridIntegral( const double* f1, const double* f2, int nint, double x0, double dx, double* ys )const{
        int ng    = n.totprod();
        double dV = voxelVolume();
        for(int i=0; i<nint; i++){
            double x = x0+dx*i;
            double Q = dV * integrateProductShifted( Vec3d{x,0,0}, f1, f2 );
            ys[i] = Q;
        }
    }

};

// interpolation of vector force-field Vec3d[ix,iy,iz] in periodic boundary condition
inline double interpolate3DWrap( double * grid, const Vec3i& n, const Vec3d& r ){
	int xoff = n.x<<3; int imx = r.x +xoff;	double tx = r.x - imx +xoff;	double mx = 1 - tx;		int itx = (imx+1)%n.x;  imx=imx%n.x;
	int yoff = n.y<<3; int imy = r.y +yoff;	double ty = r.y - imy +yoff;	double my = 1 - ty;		int ity = (imy+1)%n.y;  imy=imy%n.y;
	int zoff = n.z<<3; int imz = r.z +zoff;	double tz = r.z - imz +zoff;	double mz = 1 - tz;		int itz = (imz+1)%n.z;  imz=imz%n.z;
	int nxy = n.x * n.y; int nx = n.x;
	//double out = grid[ i3D( imx, imy, imz ) ];

	//double out = mz * my * (  ( mx * grid[ i3D( imx, imy, imz ) ] ) +  ( tx * grid[ i3D( itx, imy, imz ) ] );
	//ty * ( mx * grid[ i3D( imx, ity, imz ) ] ) +  ( tx * grid[ i3D( itx, ity, imz ) ] ) );

	double out = mz * (
	my * ( ( mx * grid[ i3D( imx, imy, imz ) ] ) +  ( tx * grid[ i3D( itx, imy, imz ) ] ) ) +
	ty * ( ( mx * grid[ i3D( imx, ity, imz ) ] ) +  ( tx * grid[ i3D( itx, ity, imz ) ] ) ) )
               + tz * (
	my * ( ( mx * grid[ i3D( imx, imy, itz ) ] ) +  ( tx * grid[ i3D( itx, imy, itz ) ] ) ) +
	ty * ( ( mx * grid[ i3D( imx, ity, itz ) ] ) +  ( tx * grid[ i3D( itx, ity, itz ) ] ) ) );
	return out;
}

// interpolation of vector force-field Vec3d[ix,iy,iz] in periodic boundary condition
inline Vec3d interpolate3DvecWrap( Vec3d * grid, const Vec3i& n, const Vec3d& r ){
	int xoff = n.x<<3; int imx = r.x +xoff;	double tx = r.x - imx +xoff;	double mx = 1 - tx;		int itx = (imx+1)%n.x;  imx=imx%n.x;
	int yoff = n.y<<3; int imy = r.y +yoff;	double ty = r.y - imy +yoff;	double my = 1 - ty;		int ity = (imy+1)%n.y;  imy=imy%n.y;
	int zoff = n.z<<3; int imz = r.z +zoff;	double tz = r.z - imz +zoff;	double mz = 1 - tz;		int itz = (imz+1)%n.z;  imz=imz%n.z;
	int nxy = n.x * n.y; int nx = n.x;
	//printf( " %f %f %f   %i %i %i \n", r.x, r.y, r.z, imx, imy, imz );
	double mymx = my*mx; double mytx = my*tx; double tymx = ty*mx; double tytx = ty*tx;
	Vec3d out;
	out.set_mul( grid[ i3D( imx, imy, imz ) ], mz*mymx );   out.add_mul( grid[ i3D( itx, imy, imz ) ], mz*mytx );
	out.add_mul( grid[ i3D( imx, ity, imz ) ], mz*tymx );   out.add_mul( grid[ i3D( itx, ity, imz ) ], mz*tytx );
	out.add_mul( grid[ i3D( imx, ity, itz ) ], tz*tymx );   out.add_mul( grid[ i3D( itx, ity, itz ) ], tz*tytx );
	out.add_mul( grid[ i3D( imx, imy, itz ) ], tz*mymx );   out.add_mul( grid[ i3D( itx, imy, itz ) ], tz*mytx );
	return out;
}

inline Quat4f interpolate3DvecWrap( Quat4f * grid, const Vec3i& n, const Vec3d& r ){
	int xoff = n.x<<3; int imx = r.x +xoff;	double tx = r.x - imx +xoff;	double mx = 1 - tx;		int itx = (imx+1)%n.x;  imx=imx%n.x;
	int yoff = n.y<<3; int imy = r.y +yoff;	double ty = r.y - imy +yoff;	double my = 1 - ty;		int ity = (imy+1)%n.y;  imy=imy%n.y;
	int zoff = n.z<<3; int imz = r.z +zoff;	double tz = r.z - imz +zoff;	double mz = 1 - tz;		int itz = (imz+1)%n.z;  imz=imz%n.z;
	int nxy = n.x * n.y; int nx = n.x;
	//printf( " %f %f %f   %i %i %i \n", r.x, r.y, r.z, imx, imy, imz );
	double mymx = my*mx; double mytx = my*tx; double tymx = ty*mx; double tytx = ty*tx;
	Quat4f out;
	out.set_mul( grid[ i3D( imx, imy, imz ) ], mz*mymx );   out.add_mul( grid[ i3D( itx, imy, imz ) ], mz*mytx );
	out.add_mul( grid[ i3D( imx, ity, imz ) ], mz*tymx );   out.add_mul( grid[ i3D( itx, ity, imz ) ], mz*tytx );
	out.add_mul( grid[ i3D( imx, ity, itz ) ], tz*tymx );   out.add_mul( grid[ i3D( itx, ity, itz ) ], tz*tytx );
	out.add_mul( grid[ i3D( imx, imy, itz ) ], tz*mymx );   out.add_mul( grid[ i3D( itx, imy, itz ) ], tz*mytx );
	return out;
}

template<typename Func>
double evalOnGrid( const GridShape& grid, Func func ){
    int nx  = grid.n.x;
    int ny  = grid.n.y;
    int nz  = grid.n.z;
    int nxy = ny * nx;
    int ii = 0;
    double res=0.0;
    for ( int ic=0; ic<nz; ic++ ){
        for ( int ib=0; ib<ny; ib++ ){
            for ( int ia=0; ia<nx; ia++ ){
                Vec3d pos;
                grid.grid2cartesian( Vec3d{ia,ib,ic}, pos );
                func( ii, pos, res );
                ii++;
            }
        }
    }
    return res;
}

// iterate over field
//template< void FUNC( int ibuff, const Vec3d& pos_, void * args ) >
template<typename FUNC>
void interateGrid3D( const GridShape& grid, FUNC func ){
    printf( "interateGrid3D() pos0(%g,%g,%g) a(%g,%g,%g) b(%g,%g,%g) c(%g,%g,%g)\n", grid.pos0.x,grid.pos0.y,grid.pos0.z,  grid.dCell.a.x,grid.dCell.a.y,grid.dCell.a.z, grid.dCell.b.x,grid.dCell.b.y,grid.dCell.b.z, grid.dCell.c.x,grid.dCell.c.y,grid.dCell.c.z );
	int nx  = grid.n.x; 	int ny  = grid.n.y; 	int nz  = grid.n.z;
	//int nx  = n.z; 	int ny  = n.y; 	int nz  = n.x;
	int nxy = ny * nx;
	//printf( "interateGrid3D nx,y,z (%i,%i,%i) nxy %i\n", nx,ny,nz, nxy );
	Vec3d pos; // pos.set( grid.pos0 );
	//printf(" interateGrid3D : args %i \n", args );
	for ( int ic=0; ic<nz; ic++ ){
        //printf("ic %i \n", ic );
		for ( int ib=0; ib<ny; ib++ ){
	        for ( int ia=0; ia<nx; ia++ ){
			    int ibuff = i3D( ia, ib, ic );
                //FUNC( ibuff, {ia,ib,ic}, pos );
                pos = grid.pos0 + grid.dCell.c*ic + grid.dCell.b*ib + grid.dCell.a*ia;
                func( ibuff, pos );
                //printf("(%i,%i,%i)(%3.3f,%3.3f,%3.3f)\n",ia,ib,ic,pos.x,pos.y,pos.z);
				//pos.add( grid.dCell.a );
			}
			//pos.add_mul( grid.dCell.a, -nx );
			//pos.add( grid.dCell.b );
		}
		//exit(0);
		//pos.add_mul( grid.dCell.b, -ny );
		//pos.add( grid.dCell.c );
	}
    printf("!!!! DEBUG interateGrid3D pmin(%g,%g,%g) pmax(%g,%g,%g) \n",  grid.pos0.x, grid.pos0.y, grid.pos0.z,    pos.x, pos.y, pos.z );
    //printf ("\n");
}


// ================================
// ======== Free Functions
// ================================


char* DEBUG_saveFile1="temp/f1.xsf";
char* DEBUG_saveFile2="temp/f2.xsf";
char* DEBUG_saveFile12="temp/prod_f1_f2.xsf";
double* DEBUG_f1=0;
double* DEBUG_f2=0;

template<typename Func1, typename Func2>
void gridNumIntegral( int nint, double gStep, double Rmax, double Lmax, double* ys, Func1 func1, Func2 func2, bool bDebugXsf = false ){
    GridShape grid;
    //grid.cell = Mat3d{ (Rmax+Lmax),0.0,0.0,  0.0,Rmax,0.0,  0.0,0.0,Rmax };
    //grid.n    = {(int)((2*Rmax+Lmax)/gStep)+1,(int)(2*Rmax/gStep)+1,(int)(2*Rmax/gStep)+1};
    grid.cell = Mat3d{ (2*Rmax+Lmax),0.0,0.0,  0.0,(2*Rmax),0.0,  0.0,0.0,(2*Rmax) };
    int ngx = (2*Rmax)/gStep;
    grid.n    = {(2*Rmax+Lmax)/gStep,ngx,ngx};  //printf( "gridNumIntegral grid.n.x %i \n", grid.n.x);
    grid.pos0 = Vec3d{-Rmax,-Rmax,-Rmax};
    grid.updateCell();
    int ng = grid.n.totprod();
    double  dV = grid.voxelVolume(); //printf( "dV %e \n", dV );
    double* f1 = new double[ ng ];
    double* f2 = new double[ ng ];
    //func1( grid, f1, 0.0 );
    //printf( "DEBUG integral{f1} %g |f^2| %g \n", grid.integrate( f1 ), VecN::dot(ng, f1, f1 )*dV );
    double dx = Lmax/nint;
    //if(bDebugXsf)grid. saveXSF( DEBUG_saveFile1, f1, -1 );
    //int iplot = nint-1;
    int iplot = 0;
    for(int i=0; i<nint; i++){
        double x = dx*i;
        //printf( "DEBUG 0 \n" );
        func1( grid, f1, x );
        func2( grid, f2, x );
        //printf( " |f1| %g \n" , VecN::sum2( ng, f1 )*dV*dV );
        //printf( " |f2| %g \n" , VecN::sum2( ng, f2 )*dV*dV );
        //printf( " |f2| %g \n" , VecN::sum( ng, f2 )*dV );
        //printf( "DEBUG 1 \n" );
        if(i==iplot){
            Vec3i iv0 = (Vec3i){0,grid.n.b/2,grid.n.c/2};
            Vec3i di  = (Vec3i){1,0,0};
            //printf( "DEBUG_f1 %li DEBUG_f2 %li \n", (long)DEBUG_f1, (long)DEBUG_f2 );
            if(DEBUG_f1) grid.cut1D( iv0, di, f1, DEBUG_f1 );
            if(DEBUG_f2) grid.cut1D( iv0, di, f2, DEBUG_f2 );
            if(bDebugXsf){ 
                grid.saveXSF( DEBUG_saveFile1, f1, 1,0 );
                grid.saveXSF( DEBUG_saveFile2, f2, 1,0 );
            }
        }
        //printf( "DEBUG 2 \n" );

        //double Q = dV * VecN::dot( ng, f1, f2 );
        double Q11=0,Q22=0,Q=0;
        for(int j=0; j<ng; j++){ Q11+=sq(f1[j]); Q22+=sq(f2[j]); f2[j]*=f1[j]; Q+=f2[j]; };
        Q*=dV;
        //printf( "[%i]  x %g Q %g Q11 %g Q22 %g dV %e \n", i, x, Q, Q11*dV, Q22*dV, dV );
        ys[i] = Q;
        //printf( "DEBUG 3 \n" );

        if(bDebugXsf&&(i==iplot)){
            //for(int j=0; j<ng; j++) f2[j]*=f1[j];
            //for(int j=0; j<ng; j++) f2[j]*=0;
            grid. saveXSF( DEBUG_saveFile12, f2, 1,0 );
        }
        //printf( "DEBUG 4 \n" );
        //printf( "[i] Q %g |  CPUtime %g [Mticks]\n", i, Q, (getCPUticks()-timeStart)*1e-6  );
    }
    delete [] f1;
    delete [] f2;
}

void getIsovalPoints_a( const GridShape& grid, double isoval, Vec3d  *FF, std::vector<Vec3d>& points ){
    int nx  = grid.n.x; 	int ny  = grid.n.y; 	int nz  = grid.n.z; int nxy = ny * nx;
    int ii = 0;
    for ( int ic=0; ic<(nz-1); ic++ ){
        for ( int ib=0; ib<ny; ib++ ){
            for ( int ia=0; ia<nx; ia++ ){
                int     i1 = i3D( ia, ib,  ic    );
                int     i2 = i3D( ia, ib, (ic+1) );
                double df1 = FF[i1].z-isoval;
                double df2 = FF[i2].z-isoval;
                //if(ii<1000)printf( " %i (%i,%i,%i) (%g,%g)\n", ii, ia,ib,ic, df1, df2 );
                if( (df1*df2)<0 ){
                    double fc = df1/(df1-df2);
                    points.push_back( grid.dCell.a*ia + grid.dCell.b*ib + grid.dCell.c*(ic+fc) );

                    int ip = points.size()-1;
                    if( ip < 1000 ) printf( " %i (%i,%i,%i) (%g,%g,%g)\n", ip, ia,ib,ic, points[ip].x, points[ip].y, points[ip].z );
                }
                ii++;
            }
        }
    }
}

void getIsoSurfZ( const GridShape& grid, double isoval, bool sign, Quat4f  *FF, Vec3d *pos, Vec3d * normal ){
    int nx  = grid.n.x; 	int ny  = grid.n.y; 	int nz  = grid.n.z; int nxy = ny * nx;
    int ii = 0;
    //printf("%i %i %i \n", nx,ny,nxy );
    for ( int ib=0; ib<ny; ib++ ){
        for ( int ia=0; ia<nx; ia++ ){
            int ibuff  = i3D( ia, ib,  0 );
            double ofz = FF[ibuff].x;
            double fz  = 0;
            int ic;
            //printf( "iba (%i,%i)\n", ib,ia );
            for ( ic=nz-1; ic>1; ic-- ){
                int ibuff_ = ibuff + nxy*ic;
                fz = FF[ibuff_].z;
                if( (fz>isoval)==sign ){
                    ibuff = ibuff_;
                    break;
                }
                ofz = fz;
            }
            //double fc = (ofz-isoval)/(ofz-fz);
            double fc = 1-((ofz-isoval)/(ofz-fz));
            //double fc = 0;
            int ibxy  = ib*nx + ia;
            //printf( "ibxy %i %i \n", ibxy, ibuff );
            pos   [ibxy] = grid.dCell.a*ia + grid.dCell.b*ib + grid.dCell.c*(ic+fc);
            //normal[ibxy] = FF[ibuff-1];
            normal[ibxy] = (Vec3d)FF[ibuff].f;
            normal[ibxy].normalize();
            //normal[ibxy] = interpolate3DvecWrap( FF, grid.n, {ia,ib, ic+fc } );
        }
    }
}

#endif









