
// OpenCL kernel for projecting sparse density matrix to a real-space grid.
// Data layout in memory (originally Fortran column-major): 
// rho(imu, inu, ineigh, iatom) -> rho[iatom][ineigh][inu][imu] in C-order indexing

typedef struct {
    float4 origin;
    float4 dA;
    float4 dB;
    float4 dC;
    int4 ngrid;
} GridSpec;

typedef struct {
    float4 pos_rcut; // x, y, z, Rcut
    int type;        // index into basis data
    int i0orb;       // start index in global orbital list
    int norb;        // number of orbitals
    int pad;
} AtomData;

typedef struct {
    int4 block_idx;  // ix, iy, iz, pad
    int n_atoms;     // number of overlapping atoms
    //int nj;          // start of jatom block ( if off-diagonal block )
    int atom_start;  // index into active_atoms list
    int pair_start;  // index into active_pairs list
    int n_pairs;     // number of active pairs (i,j)
} TaskData;

// Cubic B-spline interpolation for radial part
float evaluate_radial(
    float r, 
    int ityp, int ish, 
    __global const float* basis_data,
    int n_nodes, float dr, int max_shells
) {
    if (r >= (n_nodes - 1) * dr) return 0.0f;
    float x = r / dr;
    int i = (int)x;
    float t = x - i;
    
    // index: [ityp][ish][node]
    int base = (ityp * max_shells + ish) * n_nodes;
    
    int i0 = max(0, i - 1);
    int i1 = i;
    int i2 = min(n_nodes - 1, i + 1);
    int i3 = min(n_nodes - 1, i + 2);
    
    float y0 = basis_data[base + i0];
    float y1 = basis_data[base + i1];
    float y2 = basis_data[base + i2];
    float y3 = basis_data[base + i3];
    
    // Catmull-Rom spline
    float a = -0.5f * y0 + 1.5f * y1 - 1.5f * y2 + 0.5f * y3;
    float b = y0 - 2.5f * y1 + 2.0f * y2 - 0.5f * y3;
    float c = -0.5f * y0 + 0.5f * y2;
    float d = y1;
    
    return ((a * t + b) * t + c) * t + d;
}

__kernel void project_density_sparse(
    __global const GridSpec* grid,
    const int n_tasks,
    __global const TaskData* tasks,
    __global const AtomData* atoms,
    __global const int* active_atoms,
    __global const int4* active_pairs, // (iatom, jatom, ineigh, pad)
    __global const float* rho,        // [natoms][neigh_max][numorb_max][numorb_max]
    __global const int* neigh_j,      // [natoms][neigh_max]
    __global const float* basis_data, // [n_species][max_shells][n_nodes]
    __global const int4* species_info, // [n_species] -> (nssh, l0, l1, l2)
    const int n_nodes, 
    const float dr_basis,
    const int max_shells,
    const int neigh_max,
    const int numorb_max,
    __global float* out_grid          // [nx][ny][nz]
) {
    int gid = get_global_id(0);
    // DEBUG: emit one-line params to verify kernel entry
    if (get_global_id(0) == 0) {
        printf("GPU: kernel entry: n_tasks=%d vox_per_task=%d n_nodes=%d max_shells=%d neigh_max=%d numorb_max=%d ngrid=(%d,%d,%d)\\n",
               n_tasks, 512, n_nodes, max_shells, neigh_max, numorb_max,
               grid->ngrid.x, grid->ngrid.y, grid->ngrid.z);
    }

    int n_pairs;
    int pair_start;
    float3 r_vox;
    { // sanitize local memory space
        const int vox_per_task = 512; // 8*8*8
        const int i_task = gid / vox_per_task;
        if (i_task >= n_tasks) return;

        const int v = gid - i_task * vox_per_task; // 0..511
        const int lx = v & 7;          // v % 8
        const int ly = (v >> 3) & 7;   // (v/8) % 8
        const int lz = (v >> 6) & 7;   // (v/64) % 8
        const int3 b_idx     = tasks[i_task].block_idx.xyz;
        const int n_pairs    = tasks[i_task].n_pairs;
        const int pair_start = tasks[i_task].pair_start;
        const int gx = b_idx.x * 8 + lx;
        const int gy = b_idx.y * 8 + ly;
        const int gz = b_idx.z * 8 + lz;
        const int3 ngrid_dim    = grid->ngrid.xyz;
        if (gx >= ngrid_dim.x || gy >= ngrid_dim.y || gz >= ngrid_dim.z) return;
        r_vox = grid->origin.xyz + (float)gx * grid->dA.xyz + (float)gy * grid->dB.xyz + (float)gz * grid->dC.xyz;
        if(v==0) {
            printf("GPU task[%3i] b_idx=(%i,%i,%i) n_atoms=%i n_pairs=%i \n", i_task, b_idx.x, b_idx.y, b_idx.z, tasks[i_task].n_atoms, n_pairs);
        }
    }

    float den = 0.0f;
    // Loop over active pairs in this block
    for (int ip = 0; ip < n_pairs; ip++) {
        int4 pair_info = active_pairs[pair_start + ip];
        int i_atom = pair_info.x;
        int j_atom = pair_info.y;
        int ineigh_ij = pair_info.z;

        AtomData ad_i = atoms[i_atom];
        AtomData ad_j = atoms[j_atom];

        float4 dri,drj;
        dri.xyz = r_vox - ad_i.pos_rcut.xyz;
        dri.w   = dot(dri, dri);
        drj.xyz = r_vox - ad_j.pos_rcut.xyz;
        drj.w   = dot(drj, drj);
        
        float rcut_i = ad_i.pos_rcut.w;
        float rcut_j = ad_j.pos_rcut.w;

        if (dri.w < rcut_i*rcut_i && drj.w < rcut_j*rcut_j) {

            int4 sp_i = species_info[ad_i.type];
            int4 sp_j = species_info[ad_j.type];
            int rho_base = i_atom * neigh_max * numorb_max * numorb_max + ineigh_ij * numorb_max * numorb_max;

            const __global float4* rho_ij = (const __global  float4*)(rho + rho_base);

            dri.w = sqrt(dri.w); dri.xyz /= dri.w; // Yi = dri_hat
            drj.w = sqrt(drj.w); drj.xyz /= drj.w; // Yj = drj_hat

            dri*=exp(-dri.w);
            drj*=exp(-drj.w);

            // from dr -> WF
            // dri.w    =  evaluate_radial(dri.w, ad_i.type, 0, basis_data, n_nodes, dr_basis, max_shells); // si-orb
            // dri.xyz *=  evaluate_radial(dri.w, ad_i.type, 1, basis_data, n_nodes, dr_basis, max_shells); // pi-orb
            // drj.w    =  evaluate_radial(drj.w, ad_j.type, 0, basis_data, n_nodes, dr_basis, max_shells); // sj-orb
            // drj.xyz *=  evaluate_radial(drj.w, ad_j.type, 1, basis_data, n_nodes, dr_basis, max_shells); // pj-orb

            // 4x4 block (px,py,pz,s)_i * (px,py,pz,s)_j 
            den += dot( dri,  (
            rho_ij[0]  * drj.x + 
            rho_ij[1]  * drj.y + 
            rho_ij[2]  * drj.z + 
            rho_ij[3]  * drj.w   ) );  
            
        }
    }

    //int g_idx = (gx * ngrid_dim.y + gy) * ngrid_dim.z + gz;
    //out_grid[g_idx] = den;
    out_grid[gid] = den;

    //out_grid[gid] = gid;
   
    
}
