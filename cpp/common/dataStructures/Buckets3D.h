﻿#ifndef  Buckets3D_h
#define  Buckets3D_h
/// @file Buckets.h @brief contains Buckets class, for accelarating neighbourhood search for particle-in-cell on rectangular 3D grids.
/// @ingroup Grids
/// @ingroup Neighbours

#include "integerOps.h"
#include "fastmath.h"
#include "Vec3.h"
#include "grids3D.h"
#include "Buckets.h"

/// @brief Class for accelarating neighbourhood search for particle-in-cell on rectangular 3D grids.
class Buckets3D : public Buckets, public CubeGridRuler { public:

    bool bResizing=false;
    int* neighs_in=0;
    int  neighs_in_size=0; 

    void setup_Buckets3D( Vec3d pmin, Vec3d pmax, double step ){
        CubeGridRuler::setup( pmin, pmax, step  );
        printf( "Buckets3D::setup_Buckets3D() ntot=%i n(%i,%i,%i) step=%g pmin(%g,%g,%g) pmax(%g,%g,%g)\n", ntot, n.x,n.y,n.z, step, pmin.x,pmin.y,pmin.z, pmax.x,pmax.y,pmax.z );
        Buckets::resizeCells( ntot );
    }

    void updateNeighsBufferSize(){
        int new_N = maxInBucket*14;
        if(neighs_in_size<new_N ){  printf("Buckets3D::updateNeighsBufferSize() new_N %i \n", new_N);  neighs_in_size=new_N;  _realloc(neighs_in,neighs_in_size); }
    }

    void pointsToCells( int np, Vec3d* ps, bool* ignore=0 ){
        //printf( "Buckets3D::pointsToCells(np=%i) @ps=%li @ignore=%li @obj2cell=%li \n", np, (long)ps, (long)ignore, (long)obj2cell );
        if( bResizing || (obj2cell==0) ){ Buckets::resizeObjs( np, true ); }
        //printf( "Buckets3D::pointsToCells(np=%i) @ps=%li @ignore=%li @obj2cell=%li \n", np, (long)ps, (long)ignore, (long)obj2cell );
        if( np>nobjSize ){ printf( "Buckets3D::pointsToCells(bResizing=%i) ERROR np(%i)>nobjSize(%i) => Exit()\n", bResizing, np, nobjSize ); exit(0); }
        //printf( "Buckets3D::pointsToCells() obj2cell %li \n", (long)obj2cell  );
        if  (ignore){ 
            for( int i=0; i<np; i++ ){ 
                //printf( "Buckets3D::pointsToCells()[i=%i]  A\n", i  );
                //bool bi = ignore[i];
                //const int ic = icell( ps[i] );
                int ic; 
                if( ignore[i] ){ ic = -1; }else{ ic=icell_bound( ps[i] ); }
                //if((ic<0)||(ic>=ncell  )){ printf( "Buckets3D::pointsToCells() ic=%i ncell=%i \n", ic, ncell ); }      // Debug
                //if((i <0)||(i>=nobjSize)){ printf( "Buckets3D::pointsToCells() i=%i nobjSize=%i \n", i, nobjSize ); } // Debug
                obj2cell[i] = ic; 
                //if(!ignore[i])obj2cell[i] = icell( ps[i] );}
            } 
        }
        else { for( int i=0; i<np; i++ ){
                //printf( "Buckets3D::pointsToCells()[i=%i]  B\n", i  );
                //if(i>=nobjSize){ printf( "Buckets3D::pointsToCells() i=%i nobjSize=%i \n", i, nobjSize ); }   // Debug
                //const int ic = icell( ps[i] );
                const int ic = icell_bound( ps[i] );
                //if((ic<0)||(ic>=ncell  )){ printf( "Buckets3D::pointsToCells() ic=%i ncell=%i \n", ic, ncell ); }      // Debug
                obj2cell[i] = ic;
            } 
        }
        //printf( "Buckets3D::pointsToCells() np %i \n", np  );
        updateCells( np, obj2cell );
        //printf( "Buckets3D::pointsToCells() DONE\n"  );
    }

    int getNeighbors( Vec3i ip, int* neighs ){
        int nfound = 0;
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z+1} ), neighs+nfound );  // 1
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z+1} ), neighs+nfound );  // 2
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z+1} ), neighs+nfound );  // 3
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z+1} ), neighs+nfound );  // 4
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z+1} ), neighs+nfound );  // 5
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z+1} ), neighs+nfound );  // 6
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z+1} ), neighs+nfound );  // 7
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z+1} ), neighs+nfound );  // 8
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z+1} ), neighs+nfound );  // 9

        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z  } ), neighs+nfound );  // 10
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z  } ), neighs+nfound );  // 11
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z  } ), neighs+nfound );  // 12
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z  } ), neighs+nfound );  // 13
        //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z  } ), neighs+nfound );  // 14 Allready inside
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z  } ), neighs+nfound ); // 15
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z  } ), neighs+nfound ); // 16
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z  } ), neighs+nfound ); // 17
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z  } ), neighs+nfound ); // 18

        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z-1} ), neighs+nfound ); // 19
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z-1} ), neighs+nfound ); // 20
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z-1} ), neighs+nfound ); // 21
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z-1} ), neighs+nfound ); // 22
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z-1} ), neighs+nfound ); // 23
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z-1} ), neighs+nfound ); // 24
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z-1} ), neighs+nfound ); // 25
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z-1} ), neighs+nfound ); // 26
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z-1} ), neighs+nfound ); // 27
        return nfound;
    }


    int getForwardNeighbors( Vec3i ip, int* neighs ){
        //printf("Buckets3D::getForwardNeighbors() neighs %li \n", (long)neighs );
        int nfound = 0;
      //nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z+1} ), neighs+nfound );  // 1
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z-1} ), neighs+nfound ); // 27
      //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z+1} ), neighs+nfound );  // 2
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z-1} ), neighs+nfound ); // 26
      //nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z+1} ), neighs+nfound );  // 3
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z-1} ), neighs+nfound ); // 25
        //printf("Buckets3D::getForwardNeighbors() 3 \n");
      //nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z+1} ), neighs+nfound );  // 4
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z-1} ), neighs+nfound ); // 24
      //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z+1} ), neighs+nfound );  // 5
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z-1} ), neighs+nfound ); // 23
      //nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z+1} ), neighs+nfound );  // 6
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z-1} ), neighs+nfound ); // 22
        //printf("Buckets3D::getForwardNeighbors() 6 \n");
      //nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z+1} ), neighs+nfound );  // 7
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z-1} ), neighs+nfound ); // 21
      //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z+1} ), neighs+nfound );  // 8
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z-1} ), neighs+nfound ); // 20
      //nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z+1} ), neighs+nfound );  // 9
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z-1} ), neighs+nfound ); // 19
        //printf("Buckets3D::getForwardNeighbors() 9 \n");
      //nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z  } ), neighs+nfound );  // 10
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z  } ), neighs+nfound ); // 18
      //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z  } ), neighs+nfound );  // 11
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z  } ), neighs+nfound ); // 17
      //nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z  } ), neighs+nfound );  // 12
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z  } ), neighs+nfound ); // 16
        //printf("Buckets3D::getForwardNeighbors() 12 \n");
      //nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z  } ), neighs+nfound );  // 13
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z  } ), neighs+nfound ); // 15
        //printf("Buckets3D::getForwardNeighbors() 13 \n");
        //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z  } ), neighs+nfound );  // 14 Allready inside
        return nfound;
    }

};


#endif
