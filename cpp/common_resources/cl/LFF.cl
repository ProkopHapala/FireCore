// Linearized Force-Field (LFF) OpenCL kernel
// Jacobi iterations with Chebyshev/heavy-ball mixing performed entirely on local data.

//#define LFF_MAX_VERTS   512
//#define LFF_MAX_VERTS   512

#define LFF_WG_SIZE    32
#define MAX_NEIGHBORS   8

__kernel void lff_projective_jacobi(
    __global const int*    mols,       // 1  [nworkgroups] {i0} starting index of first atom in each molecule in the bbuffer
    __global       float4* pos,        // 2  [nSystems * nverts] positions  (xyz | charge )
    __global       float4* vel,        // 3  [nSystems * nverts] velocities (xyz | mass   )
    __global const int*    neighs,     // 4  [nSystems * nverts * LFF_MAX_NEIGHBORS] neighbor indices
    __global const float2* KLs,        // 5  [nSystems * nverts * LFF_MAX_NEIGHBORS] {K,l0} bond params per neighbor, K is stiffness, l0 is rest length
    __global const int*    fixed_mask, // 6  optional 0/1 mask of fixed vertices
    const float3           Efield,     // 7  global external electric field
    const float            dt,         // 8  global time step
    const int              nOuter,     // 9  number of outer iterations
    const int              nInner,     // 10 number of inner iterations
    const float            bMix        // 11 momentum mixing parameter
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

    __local float4 lpos[LFF_WG_SIZE];  // xyz: y, w: diag mass term
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

    // --- debug if neighbor params are correct
    if(lid == 0){
        printf("GPU: lff_projective_jacobi() dt=%f bMix=%f nOuter=%i nInner=%i Efield=(%f,%f,%f)\n", dt, bMix, nOuter, nInner, Efield.x, Efield.y, Efield.z);
        for(int i=0; i<nverts; i++){
            int ing0 = i * MAX_NEIGHBORS;
            printf("GPU: atom %3i (%6.4f,%6.4f,%6.4f|%6.4f) neigs:", i, pos[i].x, pos[i].y, pos[i].z, pos[i].w) ;
            for (int jj = 0; jj < MAX_NEIGHBORS; jj++) {
                int j     = neighs[ing0 + jj];
                float2 kl = KLs   [ing0 + jj];
                if (j < 0) break;
                printf(" (%3i,%6.4f,%6.4f)", j, kl.y, kl.x);
            }
            printf("\n");
        }
    }

    // --- private copy of neighbor data
    int    ng_idx[MAX_NEIGHBORS];
    float2 ng_KLs[MAX_NEIGHBORS];
    int nneigh = 0;
    
    for (int jj = 0; jj < MAX_NEIGHBORS; ++jj) {
        const int j = neighs[neigh_base + jj];
        ng_idx[jj]  = j; // make sure j is local (in-molecule, in-group index)
        ng_KLs[jj]  = KLs [neigh_base + jj];
        if (j < 0) break;
        ++nneigh;
    }

    // --- Projective Dynamics MD-loop
    for (int outer = 0; outer < nOuter; ++outer) {

        if(is_fixed) continue;

        // --- external force
        float3 fex    = Efield * Qi;
        // TODO: non-covalent forces can be accumulated between outer iterations.


        // DEBUG - turn of dynamics
        //fex    = (float3)(0.0f, 0.0f, 0.0f);
        //vi     = (float3)(0.0f, 0.0f, 0.0f);

        // --- Leapfrog predictor
        float3 pi_old = pi;
        vi   += fex *  (dt * inv_mass);
        pi   += vi *dt;
        lpos[lid] = (float4){pi,mi};

        // ---- Jacobi-fly Linear Constraint Solver loop
        float3 mom_vec  = (float3)(0.0f, 0.0f, 0.0f);        
        for (int iter = 0; iter < nInner; ++iter) {
            if(lid==0)printf( "GPU iter %i\n", iter );
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
                float2  kl = ng_KLs[jj];
                float3 pj  = lpos[j].xyz; // make sure j is local (in-molecule, in-group index)
                float3 dij = pi - pj;   
                float len  = length(dij);
                if(lid==0)printf( "GPU bond (%3i,%3i) l: %6.4f l0: %6.4f K: %6.4f\n", lid, j, len, kl.y, kl.x );
                float inv_len = len > 1e-8f ? 1.0f/len : 0.0f;
                float3 rest_pos = pj + dij * (kl.y * inv_len);
                bi  += rest_pos * kl.x;   // accumulate projected target
                Aii += kl.x;              // accumulate diagonal weight
            }
            const float   invA   = 1.0f / Aii;
            const float3  pi_new = bi * invA;

            // momentum acceleration
            // const float3  pi_    = pi_new + mom_vec * bMix;
            // mom_vec              = pi_ - pi;
            // pi                   = pi_;
            pi = pi_new;
            
            barrier(CLK_LOCAL_MEM_FENCE);
            lpos[lid]            = (float4)(pi, mi);
        }

        vi = (pi - pi_old)*inv_dt;
        
    }

    pos[idx] = (float4)(pi, Qi);
    vel[idx] = (float4)(vi, mi);
}
