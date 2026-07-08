/// @file TopologyNpz.h
/// @brief Topology NPZ surface-group AABB overlay (group_bbox_min/max, icolGroup).
/// Parity reference: web/common_js/nanocrystalSvg.js enrichViewerJsonWithTopology
#ifndef TopologyNpz_h
#define TopologyNpz_h

#include "NpzIO.h"
#include "CrystalNpz.h"
#include <vector>
#include <math.h>

struct TopologyBboxOverlay { public:
    bool loaded = false;
    int n_groups = 0;
    std::vector<double> group_bbox_min; // G*3
    std::vector<double> group_bbox_max; // G*3
    std::vector<int> icolGroup;         // N (optional)
    void clear(){ loaded=false; n_groups=0; group_bbox_min.clear(); group_bbox_max.clear(); icolGroup.clear(); }
};

inline TopologyBboxOverlay topology_bboxes_from_npz(const NpzFile& npz){
    const NpyArray* pmin = npz.get("group_bbox_min");
    const NpyArray* pmax = npz.get("group_bbox_max");
    if(!pmin || !pmax){ return TopologyBboxOverlay(); }
    const double* mn = npy_as_f8(*pmin);
    const double* mx = npy_as_f8(*pmax);
    int G = 0;
    if(const NpyArray* pg = npz.get("n_groups")) G = npy_as_i4(*pg)[0];
    else if(pmin->ndim==2) G = pmin->shape[0];
    else G = (int)(pmin->count()/3);
    if(G <= 0){ printf("TopologyNpz ERROR: n_groups=%i\n", G); exit(1); }
    if((int)pmin->count() != G*3 || (int)pmax->count() != G*3){
        printf("TopologyNpz ERROR: bbox array size mismatch G=%i\n", G); exit(1);
    }
    TopologyBboxOverlay o;
    o.loaded = true;
    o.n_groups = G;
    o.group_bbox_min.assign(mn, mn + G*3);
    o.group_bbox_max.assign(mx, mx + G*3);
    if(const NpyArray* pic = npz.get("icolGroup")){
        const int32_t* ig = npy_as_i4(*pic);
        o.icolGroup.assign(ig, ig + (int)pic->count());
    }
    return o;
}

inline bool topology_has_geometry_keys(const NpzFile& npz){
    return npz.get("pos") && npz.get("Z");
}

inline bool topology_has_bbox_keys(const NpzFile& npz){
    return npz.get("group_bbox_min") && npz.get("group_bbox_max");
}

/// Single-read topology/crystal NPZ load.
inline void loadTopologyNpzOnce(NpzFile& npz, Molecule& mol, MMFFparams* params, TopologyBboxOverlay& overlay){
    if(topology_has_geometry_keys(npz)){
        CrystalData c = crystal_from_npz(npz);
        crystal_fill_molecule(mol, c, params);
    }
    overlay = topology_bboxes_from_npz(npz);
    if(!overlay.loaded && !topology_has_geometry_keys(npz)){
        printf("TopologyNpz ERROR: npz has neither pos/Z nor group_bbox keys\n");
        exit(1);
    }
}

/// Wireframe color per group — parity with nanocrystalViewer.html wgColor(gid).
inline void topology_wgColor(int gid, float& r, float& g, float& b){
    double h = fmod(gid * 0.61803398875, 1.0);
    const double pi2 = 6.28318530718;
    r = (float)((100.0 + 155.0 * fabs(sin(h * pi2))) / 255.0);
    g = (float)((100.0 + 155.0 * fabs(sin((h + 0.3333) * pi2))) / 255.0);
    b = (float)((100.0 + 155.0 * fabs(sin((h + 0.6666) * pi2))) / 255.0);
}

#endif
