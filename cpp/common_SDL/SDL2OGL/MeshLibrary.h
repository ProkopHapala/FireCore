#ifndef MESH_LIBRARY_H
#define MESH_LIBRARY_H

#include "GLMesh.h"

namespace MeshLibrary {
    extern GLMesh<MPOS> point;
    extern GLMesh<MPOS> pointCross;
    
    extern GLMesh<MPOS> line;
    extern GLMesh<MPOS> line2D;

    extern GLMesh<MPOS> wireCube;
    extern GLMesh<MPOS,MNORMAL> cubeWithNormals;

    extern GLMeshBase<MPOS> sphere;
    extern GLMesh<MPOS> rect;
    extern GLMesh<MPOS> circle;
}

#endif // MESH_LIBRARY_H
