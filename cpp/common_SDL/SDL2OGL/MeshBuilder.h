#ifndef _MESH_BUILDER_H_
#define _MESH_BUILDER_H_

#include "GLvbo.h"

namespace MeshBuilder{

    // ========================
    // ===== GL_TRIANGLES =====
    // ========================

    template<attrib...As>
    inline void addTriangle(GLvbo<As...>& vbo, Vec3f p1, Vec3f p2, Vec3f p3, Vec3f color=COLOR_WHITE){
        //     p3
        //    / |
        //  p1--p2
        Vec3f normal = cross((p2-p1), (p3-p1));
        vbo.template push_back_t<MPOS,MCOLOR,MNORMAL>(p1, color, normal);
        vbo.template push_back_t<MPOS,MCOLOR,MNORMAL>(p2, color, normal);
        vbo.template push_back_t<MPOS,MCOLOR,MNORMAL>(p3, color, normal);
    }

    template<attrib...As>
    inline void addTriangle(GLvbo<As...>& vbo, Vec3f p1, Vec3f p2, Vec3f p3, Vec3f normal, Vec3f color){
        vbo.template push_back_t<MPOS,MCOLOR,MNORMAL>(p1, color, normal);
        vbo.template push_back_t<MPOS,MCOLOR,MNORMAL>(p2, color, normal);
        vbo.template push_back_t<MPOS,MCOLOR,MNORMAL>(p3, color, normal);
    }

    template<attrib...As>
    inline void addQuad(GLvbo<As...>& vbo, Vec3f p1, Vec3f p2, Vec3f p3, Vec3f p4, Vec3f color=COLOR_WHITE){
        //  p4--p3
        //  |  / |
        //  | /  |
        //  p1--p2
        addTriangle(vbo, p1, p2, p3, color);
        addTriangle(vbo, p1, p3, p4, color);
    }

    template<attrib...As>
    inline void addCube(GLvbo<As...>& vbo, Vec3f p1, Vec3f p2, Vec3f p3, Vec3f p4, Vec3f p5, Vec3f p6, Vec3f p7, Vec3f p8, Vec3f color=COLOR_WHITE){
        //     .8------7
        //   .' |    .'|
        //  4---+--3'  |
        //  |   |  |   |
        //  |  .5--+---6
        //  |.'    | .' 
        //  1------2'   

        addQuad(vbo, p1, p2, p3, p4, color);
        addQuad(vbo, p2, p6, p7, p3, color);
        addQuad(vbo, p6, p5, p8, p7, color);
        addQuad(vbo, p5, p1, p4, p8, color);
        addQuad(vbo, p4, p3, p7, p8, color);
        addQuad(vbo, p5, p6, p2, p1, color);
    }

    template<attrib...As>
    inline void addCube(GLvbo<As...>&vbo, Vec3f pMin, Vec3f pMax, Vec3f color=COLOR_WHITE){
        Vec3f p1 = pMin;
        Vec3f p2 = (Vec3f){pMax.x, pMin.y, pMin.z};
        Vec3f p3 = (Vec3f){pMax.x, pMax.y, pMin.z};
        Vec3f p4 = (Vec3f){pMin.x, pMax.y, pMin.z};
        Vec3f p5 = (Vec3f){pMin.x, pMin.y, pMax.z};
        Vec3f p6 = (Vec3f){pMax.x, pMin.y, pMax.z};
        Vec3f p7 = pMax;
        Vec3f p8 = (Vec3f){pMin.x, pMax.y, pMax.z};
        addCube(vbo, p1, p2, p3, p4, p5, p6, p7, p8, color);
    }


    // ========================
    // ======= GL_LINES =======
    // ========================

    template<attrib...As>
    inline void addLine(GLvbo<As...>& vbo, Vec3f p1, Vec3f p2, Vec3f color=COLOR_WHITE){
        vbo.template push_back_t<MPOS,MCOLOR>(p1, color);
        vbo.template push_back_t<MPOS,MCOLOR>(p2, color);
    }

    template<attrib...As>
    inline void addLine(GLvbo<As...>& vbo, std::vector<Vec3f> ps, Vec3f color=COLOR_WHITE){
        for(int i=0; i<ps.size()-1; i++){
            addLine(vbo, ps[i], ps[i+1], color);
        }
    }

    template<attrib...As>
    inline void addLineLoop(GLvbo<As...>& vbo, std::vector<Vec3f> ps, Vec3f color=COLOR_WHITE){
        addLine(vbo, ps, color);
        addLine(vbo, ps[ps.size()-1], ps[0], color);
    }

    template<attrib...As>
    inline void addCircleAxis( GLvbo<As...>& vbo, int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R, float dca, float dsa, Vec3f color=COLOR_WHITE){
        Vec3f v; v.set(v0);
        std::vector<Vec3f> ps;
        for( int i=0; i<n; i++ ){
            ps.push_back(pos+v*R);
            v.rotate_csa( dca, dsa, uaxis );
        }
        addLineLoop(vbo, ps, color);
    }

    template<attrib...As>
    inline void addCircleAxis(GLvbo<As...>& vbo, int n, Vec3f pos, Vec3f v0, Vec3f normal, float r, Vec3f color=COLOR_WHITE){
        // add a polygon made of n vertices
        // v0 is the the direction of the first vertex
        float dphi = 2*M_PI/n;
        float dca  = cos( dphi );
        float dsa  = sin( dphi );
        addCircleAxis(vbo, n, pos, v0, normal, r, dca, dsa, color);
    }

    template<attrib...As>
    inline void addSphereOct(GLvbo<As...>& vbo, int n, float r, Vec3f pos, Vec3f color=COLOR_WHITE){
        addCircleAxis(vbo, n, pos, {0, 1, 0}, {1, 0, 0}, r, color);
        addCircleAxis(vbo, n, pos, {0, 0, 1}, {0, 1, 0}, r, color);
        addCircleAxis(vbo, n, pos, {1, 0, 0}, {0, 0, 1}, r, color);
    }
};

#endif // _MESH_BUILDER_H_
