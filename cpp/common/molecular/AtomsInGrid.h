
#ifndef AtomsInGrid_h
#define AtomsInGrid_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "Mat3.h"

#include <unordered_map>

#include "grids3D.h"
#include "Atoms.h"

struct AtomInBox{
    int type; // 1=corner, 2=edge, 3=face, 4=center
    int neighs[4]; // indexes of the nearest grid points
    // if type==1(corner) neighs={ corner, -1,-1,-1 }
    // if type==2(edge)   neighs={ ,   edge, -1,-1 }   
};

struct GridPointDynamics{
    int ic = -1;
    int ia = -1;
    bool  fixed = false;
    Vec3d pos   = Vec3dZero;
    Vec3d force = Vec3dZero;
};

/**
 * @brief Class representing quadrature mesh created by embeding spherical atoms in 3D rectangular grid. 
 * It is used for integration of electron density and other functions.
 * 
 * The atomic spheres are iseerted into the grid using the following algorithm:
 *   -# we map each atom into grid box (ix,iy,iz) by rounding its center position
 *   -# we decide if the atom is closer to (i) the center of the box (ii) or to the edge or (iii)  or to the face (iv) to the corner vertex of the rectangular grid box
 *    - if the atom is closer to the corer we move (snap) the corner vertex position to the atom center
 *   -# we conect the remainging nearest atom center by edges to the nearby vertices (corners) of the rectangular grid
 *   -# we eventually elasticaly deform the grid by moving the vertices (corners) to the distribute the sample points more evenly
 */
class AtomsInGrid : public CubeGridRuler { public:

    // --- from CubeGridRuler
    // double step;
    // double invStep;
    // Vec3d  pos0;
    // Vec3d  pmax;
    // Vec3d  span;
    // Vec3i  n;
    // int    ntot,nxy;

    double Korig = 0.2; // stiffness of the original grid points
    double KL = 0.5; // stiffness along the edge
    double KT = 1.0; // stiffness perpendicular to the edge ( strightening the edge )

    int* ginds = 0;           // grid point indexes
    std::vector<GridPointDynamics> gpoints;   // grid points
    //std::vector<Vec2i> edges; 

    // --- for global optimization

    bool make_neigh( Vec3i ipos ){
        //if( !isIndexValid( ipos ) ) return false;
        //int ic = ixyz2i( ipos );
        int ic = ixyz2i_wrap( ipos );
        if( ginds[ic]<0 ){ // if grid point is not yet occupied
            ginds[ic] = gpoints.size();
            Vec3d pos = box2pos( ipos );
            gpoints.push_back( GridPointDynamics{ ic, -1, false, pos, Vec3dZero } );  
        }
        return true;
    } 

    Vec3d get_gpos( Vec3i ipos ){
        //int ic = ixyz2i( ipos );
        int ic = ixyz2i_wrap( ipos );
        int ig = ginds[ic];
        if( ig>=0 ){ // if grid point is not yet occupied
            return gpoints[ig].pos;
        }else{
            return box2pos( ipos );
        }
    }

    bool get_gpos2( Vec3i ipos, Vec3d& gpos ){
        //int ic = ixyz2i( ipos );
        int ic = ixyz2i_wrap( ipos );
        int ig = ginds[ic];
        if( ig>=0 ){ // if grid point is not yet occupied
            gpos = gpoints[ig].pos;
            return gpoints[ig].fixed;
        }else{
            gpos = box2pos( ipos );
        }
        return false;
    }

    void snap_corners( const Atoms& atoms, bool bAutoSetup=true, double margin=2.0, double step=0.5 ){
        //printf( "snap_corners(): atoms.n=%i bAutoSteup=%i, margin=%g, step=%g \n", atoms.natoms, bAutoSetup, margin, step );
        if( bAutoSetup ){
            atoms.getAABB( pos0, pmax );
            pos0.add( -margin );
            pmax.add(  margin );
            setup( pos0, pmax, step );
            //printf( "snap_corners(): n(%i,%i,%i) pmin(%g,%g,%g) pmax(%g,%g,%g) \n", n.x,n.y,n.z, pmin.x,pmin.y,pmin.z, pmax.x,pmax.y,pmax.z );
        }
        printf( "snap_corners(): n(%i|%i,%i,%i) pmin(%g,%g,%g) pmax(%g,%g,%g) \n", ntot, n.x,n.y,n.z, pos0.x,pos0.y,pos0.z, pmax.x,pmax.y,pmax.z );
        _realloc( ginds, ntot ); for(int i=0; i<ntot; i++){ ginds[i]=-1; }
        gpoints.clear();
        const Vec3d p0{0.5*step,0.5*step,0.5*step}; 
        for(int ia=0; ia<atoms.natoms; ia++){
            //printf( "snap_corners[%i] \n", ia );
            Vec3d p = atoms.apos[ia];
            Vec3i ipos;
            Vec3d dpos;
            pos2box( p+p0, ipos, dpos );
            int ic = ixyz2i( ipos );
            //printf( "snap_corners[%i] pos(%6.3f,%6.3f,%6.3f) ipos(%i,%i,%i) ic=%i \n", ia, p.x,p.y,p.z, ipos.x,ipos.y,ipos.z, ic );
            // snap grid point to atom center (i.e. set position and fix it )
            int ig = ginds[ic];
            if( ig<0 ){ // if grid point is not yet occupied
                ginds[ic] = gpoints.size();
                gpoints.push_back( GridPointDynamics{ ic, ia, true, p, Vec3dZero } );  
            }else if( !gpoints[ig].fixed ){     // if grid point is already occupied
                gpoints[ig].pos.set( p  );
                gpoints[ig].ia   =  ia ;
                gpoints[ig].fixed = true; 
            }else{
                printf( "ERROR in snap_corners[ia=%i] grid point %i is already occupied by atom #%i \n", ia, ic, gpoints[ig].ia );
                exit(0);
            }

            // create edges to the nearest grid points, and apply forces to them
            make_neigh({ipos.x+1, ipos.y,   ipos.z   });
            make_neigh({ipos.x-1, ipos.y,   ipos.z   });
            make_neigh({ipos.x,   ipos.y+1, ipos.z   });
            make_neigh({ipos.x,   ipos.y-1, ipos.z   });
            make_neigh({ipos.x,   ipos.y,   ipos.z+1 });
            make_neigh({ipos.x,   ipos.y,   ipos.z-1 } );
        }
        printf( "snap_corners() DONE: gpoints.size()=%i \n", gpoints.size() );
    }

    void eval_forces(){
        for( GridPointDynamics& gp : gpoints ){ gp.force.set(0.0); }
        for( GridPointDynamics& gp : gpoints ){
            if( gp.fixed ) continue;
            Vec3i ipos; i2ixyz( gp.ic, ipos );
            Vec3d p1,p2,f=Vec3dZero,p=gp.pos;
            f.add( (box2pos( ipos )-p)*Korig );  // force to snap to the original grid point
            p1 = get_gpos( {ipos.x+1,   ipos.y,   ipos.z } ); 
            p2 = get_gpos( {ipos.x-1,   ipos.y,   ipos.z } ); 
            f.add( ((p1+p2)*0.5 - p)*KT );                     // force to strighten the edge along x
            p1.sub(p); f.add( p1*(p1.norm()*invStep-1.0)*KL ); // force to keep the edge length
            p2.sub(p); f.add( p2*(p2.norm()*invStep-1.0)*KL ); // force to keep the edge length
            p1 = get_gpos( {ipos.x,   ipos.y+1,   ipos.z } );  
            p2 = get_gpos( {ipos.x,   ipos.y-1,   ipos.z } );  
            f.add( ((p1+p2)*0.5 - p)*KT  );                    // force to strighten the edge along y
            p1.sub(p); f.add( p1*(p1.norm()*invStep-1.0)*KL );
            p2.sub(p); f.add( p2*(p2.norm()*invStep-1.0)*KL );
            p1 = get_gpos( {ipos.x,   ipos.y,   ipos.z+1 } ); 
            p2 = get_gpos( {ipos.x,   ipos.y,   ipos.z-1 } );
            f.add( ((p1+p2)*0.5 - p)*KT  );                   // force to strighten the edge along z
            p1.sub(p); f.add( p1*(p1.norm()*invStep-1.0)*KL );
            p2.sub(p); f.add( p2*(p2.norm()*invStep-1.0)*KL );
            
            gp.force.add( f );
        }

    }

    double move_points( double dt=0.1 ){
        double f2sum = 0;
        for( GridPointDynamics& gp : gpoints ){
            if( gp.fixed ) continue;
            f2sum += gp.force.norm2();
            gp.pos.add_mul( gp.force, dt );
        }
        return f2sum;
    }




/*
    void generate_mesh_points(){
        mesh_poits.clear();
        mesh_edges.clear();
        for(int ix=0; ix<n.x; ix++){
            for(int iy=0; iy<n.y; iy++){
                for(int iz=0; iz<n.z; iz++){
                    int i   = ixyz2i( {ix,iy,iz} );
                    Vec3d p = box2pos( {ix,iy,iz} );
                    mesh_poits.push_back( p );
                    //printf( "ixyz(%i,%i,%i) i=%i p(%3.3f,%3.3f,%3.3f) \n", ix,iy,iz, i, p.x,p.y,p.z );
                    if( ix<n.x-1 ){ mesh_edges.push_back( {i,ixyz2i( {ix+1,iy,iz} )} ); }
                    if( iy<n.y-1 ){ mesh_edges.push_back( {i,ixyz2i( {ix,iy+1,iz} )} ); }
                    if( iz<n.z-1 ){ mesh_edges.push_back( {i,ixyz2i( {ix,iy,iz+1} )} ); }
                }
            }

        }

    }
*/

/*
    // THIS IS TOO COMPLICATED
    void atoms2grid( const Atoms& atoms, double* Rs=0 ){
        corner_snap.clear();
        face_snap.clear();
        edge_snap.clear();
        new_edges.clear();
        a2edge0  .clear();
        for(int ia=0; ia<atoms.natoms; ia++){
            Vec3d p = atoms.apos[ia];
            Vec3i ipos;
            Vec3d dpos;
            pos2box( p, ipos, dpos );
            //printf( "(%3.3f,%3.3f,%3.3f) (%i,%i,%i)\n", p.x, p.y, p.z, ipos.x,ipos.y,ipos.z);
            // NOTE: for indexes, every second box is edge, every 1st is center
            bool bx = ipos.x>>1;
            bool by = ipos.y>>1;
            bool bz = ipos.z>>1;
            int nc = bx + by + bz;
            //int btyp = bx + (by<<1) + (bz<<2);
            a2edge0.push_back( new_edges.size() );
            if( nc == 1 ){ // edge
                if      ( bx ){

                }else if( by ){  

                }else         {  

                }
            }else if (nc == 2 ){ //  face
                if      ( !bx ){

                }else if( !by ){

                }else          {  

                }
            }else if( nc == 0 ){ // center
                new_edges.push_back( ixyz2i( {ipos.x,   ipos.y,   ipos.z  } ) );
                new_edges.push_back( ixyz2i( {ipos.x,   ipos.y,   ipos.z+1} ) );
                new_edges.push_back( ixyz2i( {ipos.x,   ipos.y+1, ipos.z  } ) );
                new_edges.push_back( ixyz2i( {ipos.x,   ipos.y+1, ipos.z+1} ) );
                new_edges.push_back( ixyz2i( {ipos.x+1, ipos.y,   ipos.z  } ) );
                new_edges.push_back( ixyz2i( {ipos.x+1, ipos.y,   ipos.z+1} ) );
                new_edges.push_back( ixyz2i( {ipos.x+1, ipos.y+1, ipos.z  } ) );
                new_edges.push_back( ixyz2i( {ipos.x+1, ipos.y+1, ipos.z+1} ) );
            }else{ // corner
                int ic = ixyz2i( ipos );
                corner_snap.push_back( Vec2i{ic,ia} );
            }
        }
    }
    */


};

#endif
