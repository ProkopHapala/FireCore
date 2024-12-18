
#ifndef CMesh_h
#define CMesh_h
/// @file CMesh.h @brief 3D mesh class
/// @ingroup Geometry
/// @ingroup Topology

#include "Vec2.h"
#include "Vec3.h"

/// @brief minimalistic C-like class for basic 3D mesh comprising 3D vertices, edges, triangles, and polygons (faces).
class CMesh{ public:
    int nvert ;
    int nedge ;
    int ntri  ;
    int nfaces;
    Vec3d * verts;
    Vec2i * edges;
    Vec3i * tris ;  // later we may do polygon faces ?
    int   * ngons;
    int   * faces;
};

#endif






