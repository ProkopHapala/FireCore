#ifndef  BucketsND_h
#define  BucketsND_h

#include "integerOps.h"
#include "fastmath.h"
#include "VecN.h"
#include "Buckets.h"

template<typename T>
struct BoxND{
    int n=0;
    T* pmin=0;
    T* pmax=0;

    BoxND()=default;
    BoxND(int n_,T* buff=0, bool bInit=true){ if(buff){ bind(n,buff,bInit); }else{ alloc(n,bInit); } };

    void alloc(int n_,bool bInit=true){
        n=n_;
        pmin=new T[n*2];
        pmax=pmin+n;
        if(bInit)init(int n);
    }

    void bind(int n_,T* buff, bool bInit=true){
        n=n_;
        pmin=buff;
        pmax=buff+n;
        if(bInit)init(int n);
    }

    void init(){
        for(int i=0; i<n; i++){
            pmin[i]= 1e+300;
            pmax[i]=-1e+300;
        }
    }

    bool inside( const T* p ){
        for(int i=0; i<n; i++){
            if( (p[i]<pmin[i])||(p[i]>pmax[i]) ) return false;
        }
        return true;
    }

    void enclose(T* p, T R ){
        for(int i=0; i<n; i++){
            pmin[i] = fmin( pmin[i], p[i]-R );
            pmax[i] = fmax( pmax[i], p[i]+R );
        }
    }

    void copy(const BoxND& b ){ for(int i=0; i<n; i++){ pmin[i] = b.pmin[i]; pmax[i] = b.pmax[i]; } }

    void unionWith( const BoxND& b ){
        for(int i=0; i<ndim; i++){
            pmin[i] = fmin( pmin[i], b.pmin[i] );
            pmax[i] = fmax( pmax[i], b.pmax[i] );
        }
    }

    void intersection( const BoxND& b ){
        for(int i=0; i<ndim; i++){
            pmin[i] = fmax( pmin[i], b.pmin[i] );
            pmax[i] = fmin( pmax[i], b.pmax[i] );
        }
    }

};


template<typename T>
class BucketsND : public Buckets{ public:

    // void pointsToCells( int np, ** ps, bool* ignore=0 ){
    //     //printf( "Buckets3D::pointsToCells(np=%i) @ps=%li @ignore=%li @obj2cell=%li \n", np, (long)ps, (long)ignore, (long)obj2cell );
    //     if( bResizing || (obj2cell==0) ){ Buckets::resizeObjs( np, true ); }
    //     //printf( "Buckets3D::pointsToCells(np=%i) @ps=%li @ignore=%li @obj2cell=%li \n", np, (long)ps, (long)ignore, (long)obj2cell );
    //     if( np>nobjSize ){ printf( "Buckets3D::pointsToCells(bResizing=%i) ERROR np(%i)>nobjSize(%i) => Exit()\n", bResizing, np, nobjSize ); exit(0); }
    //     //printf( "Buckets3D::pointsToCells() obj2cell %li \n", (long)obj2cell  );
    //     if  (ignore){ 
    //         for( int i=0; i<np; i++ ){ 
    //             //printf( "Buckets3D::pointsToCells()[i=%i]  A\n", i  );
    //             //bool bi = ignore[i];
    //             //const int ic = icell( ps[i] );
    //             int ic; 
    //             if( ignore[i] ){ ic = -1; }else{ ic=icell_bound( ps[i] ); }
    //             //if((ic<0)||(ic>=ncell  )){ printf( "Buckets3D::pointsToCells() ic=%i ncell=%i \n", ic, ncell ); }      // Debug
    //             //if((i <0)||(i>=nobjSize)){ printf( "Buckets3D::pointsToCells() i=%i nobjSize=%i \n", i, nobjSize ); } // Debug
    //             obj2cell[i] = ic; 
    //             //if(!ignore[i])obj2cell[i] = icell( ps[i] );}
    //         } 
    //     }
    //     else { for( int i=0; i<np; i++ ){
    //             //printf( "Buckets3D::pointsToCells()[i=%i]  B\n", i  );
    //             //if(i>=nobjSize){ printf( "Buckets3D::pointsToCells() i=%i nobjSize=%i \n", i, nobjSize ); }   // Debug
    //             //const int ic = icell( ps[i] );
    //             const int ic = icell_bound( ps[i] );
    //             //if((ic<0)||(ic>=ncell  )){ printf( "Buckets3D::pointsToCells() ic=%i ncell=%i \n", ic, ncell ); }      // Debug
    //             obj2cell[i] = ic;
    //         } 
    //     }
    //     //printf( "Buckets3D::pointsToCells() np %i \n", np  );
    //     updateCells( np, obj2cell );
    //     //printf( "Buckets3D::pointsToCells() DONE\n"  );
    // }


};


#endif
