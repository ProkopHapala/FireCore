// Linearized Force-Field (LFF) OpenCL kernel
// Jacobi iterations with Chebyshev/heavy-ball mixing performed entirely on local data.

#define LFF_MAX_VERTS   512
#define MAX_NEIGHBORS   8

__kernel void lff_projective_jacobi(
    __global const int*    mols,        // [nworkgroups] {i0} starting index of first atom in each molecule in the bbuffer
    __global       float4* pos,            // [nSystems * nverts] positions  (xyz | charge )
    __global       float4* vel,            // [nSystems * nverts] velocities (xyz | mass   )
    __global const int*    neighs,         // [nSystems * nverts * LFF_MAX_NEIGHBORS] neighbor indices
    __global const float*  kngs,           // [nSystems * nverts * LFF_MAX_NEIGHBORS] spring stiffness per neighbor
    __global const float*  l0ngs,          // [nSystems * nverts * LFF_MAX_NEIGHBORS] rest lengths per neighbor
    __global const int*    fixed_mask,     // optional 0/1 mask
    const float3           Efield,         // global external electric field
    const float            dt,
    const int              nOuter,
    const int              nInner,
    const float            bMix           // momentum mixing parameter
){
    const int isys = get_group_id(0);
    const int lid  = get_local_id(0);
    const int lsz  = get_local_size(0);

    const int ia0    = mols[isys];
    const int ia1    = mols[isys + 1];
    const int nverts = ia1 - ia0;

    if (lid >= nverts) return;

    const int idx        = ia0 + lid;          // vertex index in global buffers
    const int neigh_base = idx * MAX_NEIGHBORS;

    __local float4 lpos[LFF_MAX_VERTS];  // xyz: y, w: diag mass term
    float3 pi = pos[idx].xyz;
    float3 vi = vel[idx].xyz;
    float  mi = vel[idx].w;
    float  Qi = pos[idx].w;

    lpos[lid] = (float4){pi,mi};
    const float  inv_mass = 1.0f/mi;
    const float  inv_dt   = 1.0f/dt;
    const float  inv_dt2  = inv_dt*inv_dt;
    const float  Ii       = mi*inv_dt2;

    const int is_fixed = fixed_mask[idx] != 0;

    barrier(CLK_LOCAL_MEM_FENCE);

    // --- private copy of neighbor data
    int   ng_idx[MAX_NEIGHBORS];
    float ng_k  [MAX_NEIGHBORS];
    float ng_l0 [MAX_NEIGHBORS];
    int nneigh = 0;
    
    for (int jj = 0; jj < MAX_NEIGHBORS; ++jj) {
        const int j = neighs[neigh_base + jj];
        ng_idx[jj]  = j; // make sure j is local (in-molecule, in-group index)
        ng_k  [jj]  = kngs [neigh_base + jj];
        ng_l0 [jj]  = l0ngs[neigh_base + jj];
        if (j < 0) break;
        ++nneigh;
    }

    // --- Projective Dynamics MD-loop
    for (int outer = 0; outer < nOuter; ++outer) {

        if(is_fixed) continue;

        // --- external force
        const float3 fex    = Efield * Qi;
        // TODO: non-covalent forces can be accumulated between outer iterations.

        // --- Leapfrog predictor
        float3 pi_old = pi;
        vi   += fex *  (dt * inv_mass);
        pi   += vi *dt;
        lpos[lid] = (float4){pi,mi};

        // ---- Jacobi-fly Linear Constraint Solver loop
        float3 mom_vec  = (float3)(0.0f, 0.0f, 0.0f);        
        for (int iter = 0; iter < nInner; ++iter) {
            barrier(CLK_LOCAL_MEM_FENCE);
            //const float  mi = lpos[lid].w;
            // inertial contribution: M_i/dt^2 * p_i
            //const float  Ii = mi*inv_dt2;
            float3  bi = pi*Ii;
            float  Aii = Ii;
            // ---- Loop over constraints (neighbors of given vertex)
            for (int jj=0; jj<nneigh; ++jj){
                int j      = ng_idx[jj];
                if (j < 0) break;
                // load parameters (float)
                float  k   = ng_k[jj];
                float  l0  = ng_l0[jj];
                float3 pj  = lpos[j].xyz; // make sure j is local (in-molecule, in-group index)
                float3 dij = pi - pj;   
                float len  = length(dij);
                float c    = k * ( l0/len - 1.0f );
                bi  += dij * c;      //  b_i    +=  \sum_j ( K_{ij} d_{ij} )   
                bi  += pj  * k;      //  s_j    +=  \sum_j ( K_{ij} p_j    )
                Aii += k;            //  A_{ii} +=  \sum_j   K_{ij} 
            }
            const float   invA   = 1.0f / Aii;
            const float3  pi_new = bi * invA;

            // momentum acceleration
            const float3  pi_    = pi_new + mom_vec * bMix;
            mom_vec              = pi_ - pi;
            pi                   = pi_;
            
            barrier(CLK_LOCAL_MEM_FENCE);
            lpos[lid]            = (float4)(pi, mi);
        }

        vi = (pi - pi_old)*inv_dt;
        
    }

    pos[idx] = (float4)(pi, Qi);
    vel[idx] = (float4)(vi, mi);
}
