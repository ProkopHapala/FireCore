// OpenCL kernel for projecting sparse density matrix to a real-space grid.
// Data layout in memory (originally Fortran column-major):
// rho(imu, inu, ineigh, iatom) -> rho[iatom][ineigh][inu][imu] in C-order indexing

#ifndef DEBUG_EARLY_EXIT
#define DEBUG_EARLY_EXIT 0
#endif

#ifndef DEBUG_CLEAR_ONLY
#define DEBUG_CLEAR_ONLY 0
#endif

#ifndef DEBUG_RETURN0
#define DEBUG_RETURN0 0
#endif

#ifndef DEBUG_READ_TASK
#define DEBUG_READ_TASK 0
#endif

#ifndef DEBUG_READ_GRID
#define DEBUG_READ_GRID 0
#endif

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
    int x, y, z, w;  // block_idx. x,y,z are coordinates, w is padding
    int na;          // number of overlapping atoms
    int nj;          // start of jatom block ( if off-diagonal block )
    int pad1;
    int pad2;
} TaskData;

// Real spherical harmonic normalization prefactors (same as myprog.cl)
// Y_00 = pref_s, Y_1m = pref_p * (x,y,z)/r
#define PREF_S 0.28209479f   // 1/sqrt(4*pi)
#define PREF_P 0.48860251f   // sqrt(3/(4*pi))
#define PREF_D 1.09254843f   // sqrt(15/(4*pi))
#define PREF_D_Z2 0.31539157f // for 3z^2-r^2
#define PREF_D_X2Y2 0.54627422f // for x^2-y^2

// Cubic B-spline interpolation for radial part
float evaluate_radial(
    float r, 
    int ityp, int ish, 
    __global const float* basis_data,
    int n_nodes, float dr, int max_shells
) {
    if (ityp < 0) return 0.0f;
    if (ish < 0) return 0.0f;
    if (ish >= max_shells) return 0.0f;
    if (r >= (n_nodes - 1) * dr) return 0.0f;
    // NOTE: basis_data is packed as float2 per node: (wf, wf_spline_second_derivative)
    const __global float2* basis2_ptr = (const __global float2*)basis_data;

    const float x = r / dr;
    int i = (int)floor(x);
    if (i < 0) i = 0;
    if (i > (n_nodes - 2)) i = (n_nodes - 2);
    const float t = x - (float)i;

    // interval [i, i+1]
    const int base = (ityp * max_shells + ish) * n_nodes;
    const float2 lo = basis2_ptr[base + i];
    const float2 hi = basis2_ptr[base + i + 1];

    const float a = 1.0f - t;
    const float b = t;
    // Fortran getpsi(): psi = a*ylo + b*yhi + ((a^3-a)*d2lo + (b^3-b)*d2hi)*(h^2)/6
    const float h2_6 = (dr * dr) * (1.0f/6.0f);
    const float corr = ((a*a*a - a) * lo.y + (b*b*b - b) * hi.y) * h2_6;
    return a * lo.x + b * hi.x + corr;
}

// ============================================================================
// Dense-matrix orbital/density projection helpers with d-orbital support
// ============================================================================

// Evaluate real spherical harmonic (tesseral) for given (l, mm) and unit vector rhat.
// mm ranges from -l to +l, following Fortran realTessY ordering.
// For l=2 (d orbitals): mm=-2(xy), -1(yz), 0(z2), 1(xz), 2(x2-y2)
inline float eval_angular_dense(int l, int mm, float3 rhat) {
    switch(l) {
    case 0:
        return PREF_S;
    case 1:
        switch(mm) {
        case -1: return rhat.y * PREF_P;  // py
        case 0:  return rhat.z * PREF_P;  // pz
        case 1:  return rhat.x * PREF_P;  // px
        }
        break;
    case 2:
        switch(mm) {
        case -2: return rhat.x * rhat.y * PREF_D;                                    // d_xy
        case -1: return rhat.y * rhat.z * PREF_D;                                    // d_yz
        case 0:  return (2.0f*rhat.z*rhat.z - rhat.x*rhat.x - rhat.y*rhat.y) * PREF_D_Z2; // d_z2
        case 1:  return rhat.x * rhat.z * PREF_D;                                    // d_xz
        case 2:  return (rhat.x*rhat.x - rhat.y*rhat.y) * PREF_D_X2Y2;                // d_x2y2
        }
        break;
    }
    return 0.0f;
}

// Evaluate all basis functions for one atom at a point, returning psi values into out_buf.
// out_buf must have at least norb slots. Returns number of orbitals written.
inline int eval_atom_orbitals(
    float3 r_vox,
    AtomData ad,
    __global const float* basis_data,
    int n_nodes, float dr_basis, int max_shells,
    float* out_buf
) {
    float3 d = r_vox - ad.pos_rcut.xyz;
    float r2 = dot(d, d);
    float rcut = ad.pos_rcut.w;
    if (r2 > rcut*rcut) return 0;
    float r = sqrt(r2);
    float3 rhat = d / (r + 1e-12f);

    int norb = ad.norb;
    int i0orb = ad.i0orb;
    int iorb = 0;
    for (int ish = 0; ish < max_shells && iorb < norb; ++ish) {
        float R = evaluate_radial(r, ad.type, ish, basis_data, n_nodes, dr_basis, max_shells);
        int l = ish;  // shell index = angular momentum for STO basis
        int n_m = 2*l + 1;
        for (int mm = -l; mm <= l && iorb < norb; ++mm, ++iorb) {
            out_buf[iorb] = R * eval_angular_dense(l, mm, rhat);
        }
    }
    return iorb;
}

__kernel void project_density_sparse(
    __global const GridSpec* grid,
    const int n_tasks,
    __global const TaskData* tasks,
    __global const AtomData* atoms,
    __global const int* task_atoms,   // [n_tasks][nMaxAtom]
    __global const float* rho,        // [natoms][neigh_max][numorb_max][numorb_max]
    __global const int* neigh_j,      // [natoms][neigh_max]
    __global const float* basis_data, // [n_species][max_shells][n_nodes]
    __global const int4* species_info, // [n_species] -> (nssh, l0, l1, l2)
    const int n_nodes, 
    const float dr_basis,
    const int max_shells,
    const int neigh_max,
    const int numorb_max,
    const int nMaxAtom,
    __global float* out_grid          // [nx][ny][nz]
) {
    const int gid = get_global_id(0);
    const int threads_per_task = get_local_size(0);
    // DEBUG: emit one-line params to verify kernel entry
    // if (0 && get_global_id(0) == 0) {
    //     printf("GPU: kernel entry: n_tasks=%d vox_per_task=%d n_nodes=%d max_shells=%d neigh_max=%d numorb_max=%d ngrid=(%d,%d,%d)\n",
    //            n_tasks, 512, n_nodes, max_shells, neigh_max, numorb_max,
    //            grid->ngrid.x, grid->ngrid.y, grid->ngrid.z);
    // }
    const int i_task = get_group_id(0);
    const int t_idx  = get_local_id(0);

    //{ // sanitize local memory space
    if (i_task >= n_tasks) return;
    if (DEBUG_RETURN0) return;
    __global const int* my_atoms = task_atoms + i_task * nMaxAtom;
    const TaskData task = tasks[i_task];
    const int  na    = task.na;
    const int  nj    = task.nj;
    if (DEBUG_READ_TASK) return;
    if (DEBUG_READ_GRID) { (void)grid->ngrid.x; return; }
    //const int3 b_idx = task.block_idx.xyz;
    
    //}

    // Each thread processes 32 voxels
    for (int v = t_idx; v < 512; v += threads_per_task) {
        float3 r_vox;
        int    g_idx;
        const int lx =  v       & 7;    // v     % 8
        const int ly = (v >> 3) & 7;   // (v/8 ) % 8
        const int lz = (v >> 6) & 7;   // (v/64) % 8
        { 
            //const int lx =  v       & 7;    // v     % 8
            //const int ly = (v >> 3) & 7;   // (v/8 ) % 8
            //const int lz = (v >> 6) & 7;   // (v/64) % 8
            const int gx = task.x * 8 + lx;
            const int gy = task.y * 8 + ly;
            const int gz = task.z * 8 + lz;
            const int3 ngrid_dim = grid->ngrid.xyz;
            if (gx >= ngrid_dim.x || gy >= ngrid_dim.y || gz >= ngrid_dim.z) continue;
            g_idx = (gx * ngrid_dim.y + gy) * ngrid_dim.z + gz;
            r_vox = grid->origin.xyz + (float)gx * grid->dA.xyz + (float)gy * grid->dB.xyz + (float)gz * grid->dC.xyz;
            if(v==0) { 
            //    printf("GPU task[%3i] b_idx=(%i,%i,%i) na=%i nj=%i g_idx=%i <? nxyz=%i \n", i_task, task.x, task.y, task.z, na, nj, g_idx, ngrid_dim.x*ngrid_dim.y*ngrid_dim.z ); 
            }
        }

        // if( t_idx==0 ){ 
        //     printf("GPU task[%3i] b_idx=(%i,%i,%i) na=%i nj=%i \n", i_task, task.block_idx.x, task.block_idx.y, task.block_idx.z, na, nj); 
        // }

        if (DEBUG_CLEAR_ONLY) {
            out_grid[g_idx] = 0.0f;
            continue;
        }
        float den = 0.0f;
        if (DEBUG_EARLY_EXIT) {
            out_grid[g_idx] = 0.0f;
            continue;
        }
        // Loop over active pairs in this block
        for (int i = 0; i < na; i++) {
            
            int j_start, j_end;
            if (nj < 0) { // Diagonal block: i interacting with j >= i
                j_start = i;  j_end= na;
            } else {      // Off-diagonal block: i in [0, nj) interacting with j in [nj, na)
                if (i>=nj) break; 
                j_start=nj; j_end=na;
            }

            const int i_atom = my_atoms[i]; // <-- GLOBAL READ
            AtomData ad_i    = atoms[i_atom];  // <-- GLOBAL READ
            float    rcut_i2 = ad_i.pos_rcut.w; rcut_i2*=rcut_i2;
            float4 dri;
            dri.xyz       = r_vox - ad_i.pos_rcut.xyz;
            const float ri2 = dot(dri.xyz, dri.xyz);
            if (ri2 > rcut_i2) continue;
            const float ri = sqrt(ri2);
            dri.xyz /= (ri + 1e-12f);
            dri.w =  evaluate_radial(ri, ad_i.type, 0, basis_data, n_nodes, dr_basis, max_shells) * PREF_S;
            dri.xyz *= evaluate_radial(ri, ad_i.type, 1, basis_data, n_nodes, dr_basis, max_shells) * PREF_P;
            
            for (int j = j_start; j < j_end; j++) {
                const int j_atom = my_atoms[j]; // <-- GLOBAL READ
                // Find neighbor index ineigh_ij such that neigh_j[i_atom * neigh_max + k] == j_atom + 1
                int ineigh_ij = -1;
                for (int k = 0; k < neigh_max; k++) {
                    if (neigh_j[i_atom * neigh_max + k] == j_atom + 1) {
                        ineigh_ij = k;
                        break;
                    }
                }
                if (ineigh_ij < 0) continue;

                AtomData ad_j    = atoms[j_atom];  // <-- GLOBAL READ
                float    rcut_j2 = ad_j.pos_rcut.w; rcut_j2*=rcut_j2;
                float4 drj;
                drj.xyz = r_vox - ad_j.pos_rcut.xyz;
                const float rj2 = dot(drj.xyz, drj.xyz);

                // if( t_idx==0 ){
                //     float3 rij = ad_j.pos_rcut.xyz - ad_i.pos_rcut.xyz;
                //     int rho_base = i_atom * neigh_max * numorb_max * numorb_max + ineigh_ij * numorb_max * numorb_max;
                //     if( dot(rij,rij)<(2*rcut_i2) ) printf("GPU task[%3i] rho[%i,%i] %f \n", i_task, i, j, rho[rho_base+0] ); 
                // } 
                
                if (rj2 <= rcut_j2) {
                    //int4 sp_i = species_info[ad_i.type];
                    //int4 sp_j = species_info[ad_j.type];

                    int rho_base = i_atom * neigh_max * numorb_max * numorb_max + ineigh_ij * numorb_max * numorb_max;
                    const __global float4* rho_ij_ptr = (const __global float4*)(rho + rho_base); // <-- GLOBAL READ
                    float4 rho_ij_0 = rho_ij_ptr[0];
                    float4 rho_ij_1 = rho_ij_ptr[1];
                    float4 rho_ij_2 = rho_ij_ptr[2];
                    float4 rho_ij_3 = rho_ij_ptr[3];
                    const float rj = sqrt(rj2);
                    drj.xyz /= (rj + 1e-12f);
                    drj.w =  evaluate_radial(rj, ad_j.type, 0, basis_data, n_nodes, dr_basis, max_shells) * PREF_S;
                    drj.xyz *= evaluate_radial(rj, ad_j.type, 1, basis_data, n_nodes, dr_basis, max_shells) * PREF_P;
                    // 4x4 block (px,py,pz,s)_i * (px,py,pz,s)_j 
                    // den += dot( dri,  (
                    // rho_ij[0]  * drj.x +     // <-- GLOBAL READ
                    // rho_ij[1]  * drj.y + 
                    // rho_ij[2]  * drj.z + 
                    // rho_ij[3]  * drj.w   ) );  
                    

                    // Correct formula: Σ_αβ ρ_ij[α,β] φ_i[α] φ_j[β]
                    // Compute full 4x4 matrix multiplication
                    float4 rho_i0 = rho_ij_0;  // [ρ_sx, ρ_sy, ρ_sz, ρ_ss] or similar
                    float4 rho_i1 = rho_ij_1;
                    float4 rho_i2 = rho_ij_2;
                    float4 rho_i3 = rho_ij_3;
                    
                    // den = dri · (ρ_ij · drj)
                    // where ρ_ij is 4x4 block, dri and drj are 4-vectors
                    // den = Σ_α (Σ_β ρ_ij[α,β] * drj[β]) * dri[α]
                    
                    float4 rho_dot_drj;
                    rho_dot_drj.x = rho_i0.x * drj.x + rho_i0.y * drj.y + rho_i0.z * drj.z + rho_i0.w * drj.w;
                    rho_dot_drj.y = rho_i1.x * drj.x + rho_i1.y * drj.y + rho_i1.z * drj.z + rho_i1.w * drj.w;
                    rho_dot_drj.z = rho_i2.x * drj.x + rho_i2.y * drj.y + rho_i2.z * drj.z + rho_i2.w * drj.w;
                    rho_dot_drj.w = rho_i3.x * drj.x + rho_i3.y * drj.y + rho_i3.z * drj.z + rho_i3.w * drj.w;
                    
                    den += dot(dri.wxyz, rho_dot_drj);  
                }
            }
        }
        //if( v==0 ){  // <---works
        //if( lx == 0 ){ // <---works
        // if( ly == 0 ){ // <---works
        // //if( lz == 0 ){ // <--- crashs pyopencl._cl.LogicError: clFinish failed: INVALID_COMMAND_QUEUE
        // //if( 0  == 0 ){ // <--- crashs pyopencl._cl.LogicError: clFinish failed: INVALID_COMMAND_QUEUE
            out_grid[g_idx] = den;
        //}
    }
}

typedef struct {
    int x, y, z, w;
    int na;
    int nj;
    int pad1, pad2;
} TaskData_local;

__kernel void count_atoms_per_block(
    __global const GridSpec* grid,
    const int natoms,
    __global const AtomData* atoms,
    const int block_res,
    const int n_blocks_x,
    const int n_blocks_y,
    const int n_blocks_z,
    __global int* block_counts
) {
    const int ia = get_global_id(0);
    if (ia >= natoms) return;

    AtomData ad = atoms[ia];
    float3 pos = ad.pos_rcut.xyz;
    float rcut = ad.pos_rcut.w;

    // Find range of blocks this atom can touch
    float3 r_min = (pos - rcut - grid->origin.xyz);
    float3 r_max = (pos + rcut - grid->origin.xyz);
    
    // Convert to block indices using floor (since origin is grid zero)
    // NOTE: dCell is dA.x, dB.y, dC.z assuming orthogonal grid for simplicity in indexing
    float3 block_size = (float)block_res * (float3)(grid->dA.x, grid->dB.y, grid->dC.z);
    int3 b0 = convert_int3(floor(r_min / block_size));
    int3 b1 = convert_int3(floor(r_max / block_size));

    b0 = clamp(b0, (int3)0, (int3)(n_blocks_x-1, n_blocks_y-1, n_blocks_z-1));
    b1 = clamp(b1, (int3)0, (int3)(n_blocks_x-1, n_blocks_y-1, n_blocks_z-1));

    for (int ix = b0.x; ix <= b1.x; ix++) {
        for (int iy = b0.y; iy <= b1.y; iy++) {
            for (int iz = b0.z; iz <= b1.z; iz++) {
                // Sphere-AABB check for each candidate block
                float3 b_min = grid->origin.xyz + (float)ix * block_res * grid->dA.xyz + (float)iy * block_res * grid->dB.xyz + (float)iz * block_res * grid->dC.xyz;
                float3 b_max = b_min + (float)block_res * (grid->dA.xyz + grid->dB.xyz + grid->dC.xyz);
                
                float3 closest_p = clamp(pos, b_min, b_max);
                float3 diff = pos - closest_p;
                if (dot(diff, diff) < rcut * rcut) {
                    int b_idx = (ix * n_blocks_y + iy) * n_blocks_z + iz;
                    atomic_inc(&block_counts[b_idx]);
                }
            }
        }
    }
}

__kernel void fill_task_atoms(
    __global const GridSpec* grid,
    const int natoms,
    __global const AtomData* atoms,
    const int block_res,
    const int n_blocks_x,
    const int n_blocks_y,
    const int n_blocks_z,
    __global int* block_offsets, // used for atomic fetch-add to write atom ids
    __global int* task_atoms,    // [n_blocks][nMaxAtom]
    const int nMaxAtom
) {
    const int ia = get_global_id(0);
    if (ia >= natoms) return;

    AtomData ad = atoms[ia];
    float3 pos = ad.pos_rcut.xyz;
    float rcut = ad.pos_rcut.w;

    float3 r_min = (pos - rcut - grid->origin.xyz);
    float3 r_max = (pos + rcut - grid->origin.xyz);
    float3 block_size = (float)block_res * (float3)(grid->dA.x, grid->dB.y, grid->dC.z);
    int3 b0 = convert_int3(floor(r_min / block_size));
    int3 b1 = convert_int3(floor(r_max / block_size));
    b0 = clamp(b0, (int3)0, (int3)(n_blocks_x-1, n_blocks_y-1, n_blocks_z-1));
    b1 = clamp(b1, (int3)0, (int3)(n_blocks_x-1, n_blocks_y-1, n_blocks_z-1));

    for (int ix = b0.x; ix <= b1.x; ix++) {
        for (int iy = b0.y; iy <= b1.y; iy++) {
            for (int iz = b0.z; iz <= b1.z; iz++) {
                float3 b_min = grid->origin.xyz + (float)ix * block_res * grid->dA.xyz + (float)iy * block_res * grid->dB.xyz + (float)iz * block_res * grid->dC.xyz;
                float3 b_max = b_min + (float)block_res * (grid->dA.xyz + grid->dB.xyz + grid->dC.xyz);
                float3 closest_p = clamp(pos, b_min, b_max);
                float3 diff = pos - closest_p;
                if (dot(diff, diff) < rcut * rcut) {
                    int b_idx = (ix * n_blocks_y + iy) * n_blocks_z + iz;
                    int slot = atomic_inc(&block_offsets[b_idx]);
                    if (slot < nMaxAtom) {
                        task_atoms[b_idx * nMaxAtom + slot] = ia;
                    }
                }
            }
        }
    }
}

__kernel void compact_tasks(
    const int n_blocks_x,
    const int n_blocks_y,
    const int n_blocks_z,
    __global const int* block_counts,
    __global const int* task_offsets, // prefix sum of (block_counts > 0)
    __global const int* task_atoms_raw, // [n_blocks][nMaxAtom]
    __global TaskData_local* tasks_out,
    __global int* task_atoms_out,
    const int nMaxAtom
) {
    const int ix = get_global_id(0);
    const int iy = get_global_id(1);
    const int iz = get_global_id(2);

    if (ix >= n_blocks_x || iy >= n_blocks_y || iz >= n_blocks_z) return;

    int b_idx = (ix * n_blocks_y + iy) * n_blocks_z + iz;
    int na = block_counts[b_idx];
    if (na > 0) {
        int t_idx = task_offsets[b_idx];
        TaskData_local task;
        task.x = ix; task.y = iy; task.z = iz; task.w = 0;
        task.na = na;
        task.nj = -1;
        task.pad1 = 0; task.pad2 = 0;
        tasks_out[t_idx] = task;
        
        for (int k = 0; k < nMaxAtom; k++) {
            task_atoms_out[t_idx * nMaxAtom + k] = task_atoms_raw[b_idx * nMaxAtom + k];
        }
    }
}

// Tiled kernel to avoid atomic adds and minimize global reads
#ifndef TILE_ATOMS
#define TILE_ATOMS 8
#endif
__kernel void project_density_sparse_tiled(
    __global const GridSpec* grid,
    const int n_tasks,
    __global const TaskData* tasks,
    __global const AtomData* atoms,
    __global const int* task_atoms,   // [n_tasks][nMaxAtom]
    __global const float* rho,        // [natoms][neigh_max][numorb_max][numorb_max]
    __global const int* neigh_j,      // [natoms][neigh_max]
    __global const float* basis_data, // [n_species][max_shells][n_nodes]
    __global const int4* species_info, // [n_species] -> (nssh, l0, l1, l2)
    const int n_nodes, 
    const float dr_basis,
    const int max_shells,
    const int neigh_max,
    const int numorb_max,
    const int nMaxAtom,
    __global float* out_grid          // [nx][ny][nz]
) {
    const int i_task = get_group_id(0);
    const int t_idx  = get_local_id(0);
    const int threads_per_task = get_local_size(0);
    
    if (i_task >= n_tasks) return;
    if (DEBUG_RETURN0) return;

    TaskData task = tasks[i_task];
    if (DEBUG_READ_TASK) return;
    if (DEBUG_READ_GRID) { (void)grid->ngrid.x; return; }

    // Clean up output grid so we can accumulate into it   --- to avoid need for  l_den[512];
    for (int v = t_idx; v < 512; v += threads_per_task) {
        const int lx =  v       & 7;
        const int ly = (v >> 3) & 7;
        const int lz = (v >> 6) & 7;
        const int gx = task.x * 8 + lx;
        const int gy = task.y * 8 + ly;
        const int gz = task.z * 8 + lz;
        const int3 ngrid_dim = grid->ngrid.xyz;
        if (gx < ngrid_dim.x && gy < ngrid_dim.y && gz < ngrid_dim.z) {
            int g_idx = (gx * ngrid_dim.y + gy) * ngrid_dim.z + gz;
            out_grid[g_idx] = 0.0f;
        }
    }

    const int na = task.na;
    __local AtomData l_atom_i[TILE_ATOMS];
    __local AtomData l_atom_j[TILE_ATOMS];
    __local float4   l_rho[TILE_ATOMS*TILE_ATOMS*4]; // TILE_ATOMS x TILE_ATOMS atoms * 4 orbitals
    //__local float  l_den[512];  // Store accumulated density for each voxel in the block

    // Initialize local density buffer
    //for (int v = t_idx; v < 512; v += threads_per_task) {    l_den[v] = 0.0f;}

    barrier(CLK_LOCAL_MEM_FENCE);

    if (DEBUG_CLEAR_ONLY) return;

    if (DEBUG_EARLY_EXIT) return;

    // Tiled interaction loop: TILE_ATOMS x TILE_ATOMS atoms at a time
    for (int it = 0; it < na; it += TILE_ATOMS) {
        // Load atom_i block to local memory
        if (t_idx < TILE_ATOMS && (it + t_idx) < na) {
            int i_atom = task_atoms[i_task * nMaxAtom + it + t_idx];
            l_atom_i[t_idx] = atoms[i_atom];
        }
        
        for (int jt = 0; jt < na; jt += TILE_ATOMS) {
            if (task.nj < 0 && jt < it) continue; 

            // Load atom_j block to local memory
            if (t_idx < TILE_ATOMS && (jt + t_idx) < na) {
                int j_atom = task_atoms[i_task * nMaxAtom + jt + t_idx];
                l_atom_j[t_idx] = atoms[j_atom];
            }

            // Load rho_ij blocks for the tile
            // Each thread can help load TILE_ATOMS*TILE_ATOMS*4 float4s
            for (int k = t_idx; k < TILE_ATOMS*TILE_ATOMS*4; k += threads_per_task) {
                int pair_idx = k / 4;
                int orb_idx  = k % 4;
                int i_in_tile = pair_idx / TILE_ATOMS;
                int j_in_tile = pair_idx % TILE_ATOMS;
                
                int i = it + i_in_tile;
                int j = jt + j_in_tile;
                
                bool active = (i < na && j < na);
                if (active && task.nj < 0 && j < i) active = false;
                if (active && task.nj >= 0 && (i >= task.nj || j < task.nj)) active = false;

                if (active) {
                    int i_atom = task_atoms[i_task * nMaxAtom + i];
                    int j_atom = task_atoms[i_task * nMaxAtom + j];
                    
                    int ineigh_ij = -1;
                    for (int n = 0; n < neigh_max; n++) {
                        if (neigh_j[i_atom * neigh_max + n] == j_atom + 1) {
                            ineigh_ij = n;
                            break;
                        }
                    }
                    
                    if (ineigh_ij >= 0) {
                        int rho_base = i_atom * neigh_max * numorb_max * numorb_max + ineigh_ij * numorb_max * numorb_max;
                        float4 rho_val = ((__global float4*)(rho + rho_base))[orb_idx];
                        l_rho[k] = rho_val;
                    } else {
                        l_rho[k] = (float4)(0.0f);
                    }
                } else {
                    l_rho[k] = (float4)(0.0f);
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);

            // Loop over voxels in the 8x8x8 block and accumulate density from this tile
            for (int v = t_idx; v < 512; v += threads_per_task) {
                const int lx =  v       & 7;
                const int ly = (v >> 3) & 7;
                const int lz = (v >> 6) & 7;
                const int gx = task.x * 8 + lx;
                const int gy = task.y * 8 + ly;
                const int gz = task.z * 8 + lz;
                const int3 ngrid_dim = grid->ngrid.xyz;
                
                if (gx >= ngrid_dim.x || gy >= ngrid_dim.y || gz >= ngrid_dim.z) continue;
                
                float3 r_vox = grid->origin.xyz + (float)gx * grid->dA.xyz + (float)gy * grid->dB.xyz + (float)gz * grid->dC.xyz;
                float den_tile = 0.0f;

                for (int i_in_tile = 0; i_in_tile < TILE_ATOMS && (it + i_in_tile) < na; i_in_tile++) {
                    const AtomData ad_i = l_atom_i[i_in_tile];
                    const float3 dpos_i = r_vox - ad_i.pos_rcut.xyz;
                    const float  r2_i   = dot(dpos_i, dpos_i);
                    const float  rcut_i = ad_i.pos_rcut.w;
                    if (r2_i > rcut_i*rcut_i) continue;
                    const float r_i      = sqrt(r2_i);
                    const float si = evaluate_radial(r_i, ad_i.type, 0, basis_data, n_nodes, dr_basis, max_shells) * PREF_S;
                    const float pi = evaluate_radial(r_i, ad_i.type, 1, basis_data, n_nodes, dr_basis, max_shells) * PREF_P;
                    const float sc_i = (pi / (r_i + 1e-12f));
                    const float4 psi_i = (float4)(dpos_i.x*sc_i, dpos_i.y*sc_i, dpos_i.z*sc_i, si);

                    const int j_start_tile = (task.nj < 0 && jt == it) ? i_in_tile : 0;
                    for (int j_in_tile = j_start_tile; j_in_tile < TILE_ATOMS && (jt + j_in_tile) < na; j_in_tile++) {
                        const int i = it + i_in_tile;
                        const int j = jt + j_in_tile;
                        if (task.nj >= 0 && (i >= task.nj || j < task.nj)) continue;

                        const AtomData ad_j = l_atom_j[j_in_tile];
                        const float3 dpos_j = r_vox - ad_j.pos_rcut.xyz;
                        const float r2_j    = dot(dpos_j, dpos_j);
                        const float rcut_j  = ad_j.pos_rcut.w;
                        if (r2_j > rcut_j*rcut_j) continue;
                        const float r_j  = sqrt(r2_j);
                        const float sj = evaluate_radial(r_j, ad_j.type, 0, basis_data, n_nodes, dr_basis, max_shells) * PREF_S;
                        const float pj = evaluate_radial(r_j, ad_j.type, 1, basis_data, n_nodes, dr_basis, max_shells) * PREF_P;
                        const float sc_j = (pj / (r_j + 1e-12f));
                        const float4 psi_j = (float4)(dpos_j.x*sc_j, dpos_j.y*sc_j, dpos_j.z*sc_j, sj);

                        const int tile_rho_base = (i_in_tile * TILE_ATOMS + j_in_tile) * 4;
                        const float pairsym = (task_atoms[i_task * nMaxAtom + i] == task_atoms[i_task * nMaxAtom + j]) ? 1.0f : 2.0f;
                        
                        // In Ortega convention (s, py, pz, px): s-s is [0,0], s-pz is [0,2], pz-pz is [2,2]
                        den_tile += pairsym * dot( psi_i.wyzx, (
                            l_rho[tile_rho_base + 0] * psi_j.w +
                            l_rho[tile_rho_base + 1] * psi_j.y +
                            l_rho[tile_rho_base + 2] * psi_j.z +
                            l_rho[tile_rho_base + 3] * psi_j.x ) );
                    }
                }
                //l_den[v] += den_tile;

                const int g_idx = (gx * ngrid_dim.y + gy) * ngrid_dim.z + gz;
                out_grid[g_idx] += den_tile;

            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }

/*
    // Final write back to global memory
    for (int v = t_idx; v < 512; v += threads_per_task) {
        const int lx =  v       & 7;
        const int ly = (v >> 3) & 7;
        const int lz = (v >> 6) & 7;
        const int gx = task.x * 8 + lx;
        const int gy = task.y * 8 + ly;
        const int gz = task.z * 8 + lz;
        const int3 ngrid_dim = grid->ngrid.xyz;

        if (gx < ngrid_dim.x && gy < ngrid_dim.y && gz < ngrid_dim.z) {
            int g_idx = (gx * ngrid_dim.y + gy) * ngrid_dim.z + gz;
            out_grid[g_idx] = l_den[v];
        }
    }
*/
}

// ============================================================================
// Orbital projection kernel (simpler than density projection)
// Computes ψ(r) = Σ_i C_i φ_i(r) where C_i are MO coefficients
// This is for projecting a single orbital, not density
// 
// IMPORTANT: Coeffs must be in OpenCL order [s, px, py, pz] (remapped from Fortran [s, py, pz, px])
// Angular dependence: psi_s = R_s(r) * Y_00
//                     psi_px = R_p(r) * (x/r) * PREF_P  [Y_1,+1]
//                     psi_py = R_p(r) * (y/r) * PREF_P  [Y_1,-1]  
//                     psi_pz = R_p(r) * (z/r) * PREF_P  [Y_1,0]
// ============================================================================
__kernel void project_orbital(
    __global const GridSpec* grid,
    const int n_tasks,
    __global const TaskData* tasks,
    __global const AtomData* atoms,
    __global const int* task_atoms,
    __global const float* coeffs,     // [natoms][numorb_max] - MO coefficients in OpenCL order [s, px, py, pz]
    __global const float* basis_data,
    const int n_nodes,
    const float dr_basis,
    const int max_shells,
    const int numorb_max,
    const int nMaxAtom,
    __global float* out_grid
) {
    const int gid = get_global_id(0);
    const int threads_per_task = get_local_size(0);
    const int i_task = get_group_id(0);
    const int t_idx = get_local_id(0);
    if (i_task >= n_tasks) return;

    const TaskData task = tasks[i_task];
    const int na = task.na;

    for (int v = t_idx; v < 512; v += threads_per_task) {
        float3 r_vox;
        int g_idx;
        const int lx = v & 7;
        const int ly = (v >> 3) & 7;
        const int lz = (v >> 6) & 7;
        {
            const int gx = task.x * 8 + lx;
            const int gy = task.y * 8 + ly;
            const int gz = task.z * 8 + lz;
            const int3 ngrid_dim = grid->ngrid.xyz;
            if (gx >= ngrid_dim.x || gy >= ngrid_dim.y || gz >= ngrid_dim.z) continue;
            g_idx = (gx * ngrid_dim.y + gy) * ngrid_dim.z + gz;
            r_vox = grid->origin.xyz + (float)gx * grid->dA.xyz + (float)gy * grid->dB.xyz + (float)gz * grid->dC.xyz;
        }

        float psi = 0.0f;

        for (int i = 0; i < na; i++) {
            const int i_atom = task_atoms[i_task * nMaxAtom + i];
            AtomData ad_i = atoms[i_atom];
            float rcut_i2 = ad_i.pos_rcut.w;
            rcut_i2 *= rcut_i2;

            float3 dri;
            dri = r_vox - ad_i.pos_rcut.xyz;
            const float ri2 = dot(dri, dri);
            if (ri2 > rcut_i2) continue;
            const float ri = sqrt(ri2);

            // Evaluate radial parts
            const float rs = evaluate_radial(ri, ad_i.type, 0, basis_data, n_nodes, dr_basis, max_shells);
            const float rp = evaluate_radial(ri, ad_i.type, 1, basis_data, n_nodes, dr_basis, max_shells);
            
            // Angular unit vector (x/r, y/r, z/r) for p-orbital angular dependence
            const float3 rhat = dri / (ri + 1e-12f);

            // Sum over orbitals: ψ += C_i * φ_i(r)
            // Coeffs are in order [px, py, pz, s] matching existing FireballOCL convention
            // See DFT/utils.py convCoefs() which produces [px, py, pz, s] order
            const int coeff_base = i_atom * numorb_max;
            const int norb = ad_i.norb;

            // Compute basis function values: [px*rhat.x, py*rhat.y, pz*rhat.z, s]
            // This matches sp3_tex() in myprog.cl: return (float4)(dp*ir*pref_p, pref_s)*fr
            const float4 basis_val = (float4)(
                rp * rhat.x * PREF_P,   // px component
                rp * rhat.y * PREF_P,   // py component
                rp * rhat.z * PREF_P,   // pz component
                rs * PREF_S             // s component
            );

            // Dot product with coefficients [px, py, pz, s]
            const __global float4* coeffs_ptr = (const __global float4*)(coeffs);
            float4 coeff_val = coeffs_ptr[coeff_base / 4];
            psi += dot(coeff_val, basis_val);
    }

        out_grid[g_idx] = psi;
    }
}

// ============================================================================
// Orbital projection at arbitrary points with exponential (vacuum) radial decay
// Computes ψ(p_k) = Σ_atoms Σ_orb C_{atom,orb} φ_{atom,orb}(p_k)
// Radial part is replaced by exp(-beta*(r-r0)) instead of Fireball basis tables.
// Coeff convention: [px, py, pz, s] per atom (float4)
// Angular part: cartesian p: rhat.x/y/z (same as project_orbital_points)
// ============================================================================
__kernel void project_orbital_points_exp(
    const int n_points,
    __global const float4* points,    // [n_points] xyz
    __global const AtomData* atoms,   // [natoms]
    const int natoms,
    __global const float* coeffs,     // [natoms*4] packed as float4
    const float beta,
    const float r0,
    __global float* out_psi           // [n_points]
) {
    const int ip = get_global_id(0);
    if (ip >= n_points) return;

    const float3 p = points[ip].xyz;
    float psi = 0.0f;

    for (int ia = 0; ia < natoms; ia++) {
        const AtomData ad = atoms[ia];
        float3 d = p - ad.pos_rcut.xyz;
        const float r2 = dot(d, d);
        const float rcut2 = ad.pos_rcut.w * ad.pos_rcut.w;
        if (r2 > rcut2) continue;
        const float r = sqrt(r2);

        const float invr = 1.0f / (r + 1e-12f);
        const float3 rhat = d * invr;

        const float f = exp(-beta * (r - r0));
        const float rs = f;
        const float rp = f;

        const float4 basis_val = (float4)(
            rp * rhat.x * PREF_P,
            rp * rhat.y * PREF_P,
            rp * rhat.z * PREF_P,
            rs * PREF_S
        );

        const int coeff_base = ad.i0orb;
        const __global float4* coeffs_ptr = (const __global float4*)(coeffs);
        float4 coeff_val = coeffs_ptr[coeff_base / 4];
        psi += dot(coeff_val, basis_val);
    }

    out_psi[ip] = psi;
}

// ============================================================================
// MO overlap for molecular tip vs molecular sample using exponential radial + SK angular
// Each work-item corresponds to one tip-center position (scan pixel).
//
// For given tip center R_tip, tip atoms are at R_tip + tip_pos_rel[ia].
// We compute overlap amplitude:
//   t(R_tip) = Σ_{ia∈tip} Σ_{ja∈smp}  cT(ia)^T  S_ij(R_ij)  cS(ja)
// where S_ij is a 4x4 (s,px,py,pz) SK-like block with radial f=exp(-beta*(r-r0)).
// Coeff convention on input: float4 per atom in order [px,py,pz,s].
// Output:
//   out_t[ip]  = signed amplitude (float)
//   out_I[ip]  = |t|^2 (float)
// ============================================================================
inline float sk_contract_sp(
    const float4 cT_pxpy_pz_s,
    const float4 cS_pxpy_pz_s,
    const float l, const float m, const float n,
    const float Vss, const float Vsp, const float Vps,
    const float Vpp_sig, const float Vpp_pi
){
    // reorder [px,py,pz,s] -> (s,px,py,pz)
    const float sT  = cT_pxpy_pz_s.w;
    const float pxT = cT_pxpy_pz_s.x;
    const float pyT = cT_pxpy_pz_s.y;
    const float pzT = cT_pxpy_pz_s.z;
    const float sS  = cS_pxpy_pz_s.w;
    const float pxS = cS_pxpy_pz_s.x;
    const float pyS = cS_pxpy_pz_s.y;
    const float pzS = cS_pxpy_pz_s.z;

    const float pT_dot_pS = pxT*pxS + pyT*pyS + pzT*pzS;
    const float pT_dot_u  = pxT*l  + pyT*m  + pzT*n;
    const float pS_dot_u  = pxS*l  + pyS*m  + pzS*n;

    const float t_ss = sT * Vss * sS;
    const float t_sp = sT * Vsp * (l*pxS + m*pyS + n*pzS);
    const float t_ps = (l*pxT + m*pyT + n*pzT) * Vps * sS;
    const float d = Vpp_sig - Vpp_pi;
    const float t_pp = Vpp_pi * pT_dot_pS + d * pT_dot_u * pS_dot_u;
    return t_ss + t_sp + t_ps + t_pp;
}

inline float3 quat_rotate3( const float4 q, const float3 v ){
    // q = (x,y,z,w), unit quaternion; rotate v by q
    const float3 qv = (float3)(q.x,q.y,q.z);
    const float3 t  = 2.0f*cross(qv, v);
    return v + q.w*t + cross(qv, t);
}

__kernel void mo_overlap_points_exp_sk(
    // Scan grid: each work-item computes overlap for one tip-center position
    const int n_points,                       // number of scan points (pixels)
    __global const float4* tip_centers,      // [n_points] tip-center positions (x,y,z,NA) in Å; : These are the lateral shifts of the rigid tip relative to sample
    __global const float4* tip_quat,         // [n_points] tip rotation quaternion (x,y,z,w), unit; : Applied per pixel/work-item (can vary across the scan)
    __global const float4* tip_pos_rel,      // [ntip_atoms] tip atom positions relative to tip center (x,y,z,NA) in Å; :  This geometry is already rotated on the host side (if rotation is desired)
    __global const float4* smp_pos,          // [nsmp_atoms] sample atom absolute positions (x,y,z, NA) in Å
    const int ntip_atoms,                    // number of tip atoms
    const int nsmp_atoms,                    // number of sample atoms
    __global const float4* coeffs_tip,       // [ntip_atoms] MO coefficients for tip, per atom as float4 [px,py,pz,s]: Order: .x=px, .y=py, .z=pz, .w=s (cartesian p orbitals)
    __global const float4* coeffs_smp,       // [nsmp_atoms] MO coefficients for sample, per atom as float4 [px,py,pz,s]
    // Exponential radial decay parameters: f(r) = exp(-beta*(r - r0))
    const float beta,                        // decay constant (Å^-1)
    const float r0,                          // reference distance (Å) where f=1
    const float rcut,                        // distance cutoff (Å); atom pairs beyond rcut are skipped
    __global float* out_t,                   // [n_points] output signed overlap amplitude t = c_tip^T S_ts c_sample
    __global float* out_I                    // [n_points] output intensity |t|^2
){
    const int ip = get_global_id(0);
    if(ip >= n_points) return;
    const float3 cen = tip_centers[ip].xyz;
    const float4 q   = tip_quat[ip];
    const float rcut2 = rcut*rcut;
    float t = 0.0f;

    for(int ia=0; ia<ntip_atoms; ia++){
        const float3 pT = cen + quat_rotate3(q, tip_pos_rel[ia].xyz);
        const float4 cT0 = coeffs_tip[ia];
        const float3 pcoefT = quat_rotate3(q, (float3)(cT0.x,cT0.y,cT0.z));
        const float4 cT = (float4)(pcoefT.x,pcoefT.y,pcoefT.z,cT0.w);
        for(int ja=0; ja<nsmp_atoms; ja++){
            const float3 d = pT - smp_pos[ja].xyz;
            const float r2 = dot(d,d);
            if(r2 > rcut2) continue;
            const float r = sqrt(r2);
            const float invr = 1.0f/(r + 1e-12f);
            const float l = d.x*invr;
            const float m = d.y*invr;
            const float n = d.z*invr;
            const float f = exp(-beta*(r - r0));
            // Pure exp radial + SK angular (no extra per-channel amplitudes):
            // Use fixed sign convention to keep p-p sigma/pi anisotropy.
            const float Vss = -f;
            const float Vsp = -f;
            const float Vps = -f;
            const float Vpp_sig = -f;
            const float Vpp_pi  = +f;
            t += sk_contract_sp(cT, coeffs_smp[ja], l,m,n, Vss, Vsp, Vps, Vpp_sig, Vpp_pi);
        }
    }
    out_t[ip] = t;
    out_I[ip] = t*t;
}

__kernel void mo_overlap_points_exp_sk_2mol(
    // Explicit two-molecule entrypoint (tip and sample may be different molecules)
    // NOTE: Implementation is identical to mo_overlap_points_exp_sk; we keep it separate
    //       to avoid breaking existing workflows and to make call sites self-documenting.
    const int n_points,
    __global const float4* tip_centers,
    __global const float4* tip_quat,
    __global const float4* tip_pos_rel,
    __global const float4* smp_pos,
    const int ntip_atoms,
    const int nsmp_atoms,
    __global const float4* coeffs_tip,
    __global const float4* coeffs_smp,
    const float beta,
    const float r0,
    const float rcut,
    __global float* out_t,
    __global float* out_I
){
    const int ip = get_global_id(0);
    if(ip >= n_points) return;
    const float3 cen = tip_centers[ip].xyz;
    const float4 q   = tip_quat[ip];
    const float rcut2 = rcut*rcut;
    float t = 0.0f;

    for(int ia=0; ia<ntip_atoms; ia++){
        const float3 pT = cen + quat_rotate3(q, tip_pos_rel[ia].xyz);
        const float4 cT0 = coeffs_tip[ia];
        const float3 pcoefT = quat_rotate3(q, (float3)(cT0.x,cT0.y,cT0.z));
        const float4 cT = (float4)(pcoefT.x,pcoefT.y,pcoefT.z,cT0.w);
        for(int ja=0; ja<nsmp_atoms; ja++){
            const float3 d = pT - smp_pos[ja].xyz;
            const float r2 = dot(d,d);
            if(r2 > rcut2) continue;
            const float r = sqrt(r2);
            const float invr = 1.0f/(r + 1e-12f);
            const float l = d.x*invr;
            const float m = d.y*invr;
            const float n = d.z*invr;
            const float f = exp(-beta*(r - r0));
            const float Vss = -f;
            const float Vsp = -f;
            const float Vps = -f;
            const float Vpp_sig = -f;
            const float Vpp_pi  = +f;
            t += sk_contract_sp(cT, coeffs_smp[ja], l,m,n, Vss, Vsp, Vps, Vpp_sig, Vpp_pi);
        }
    }
    out_t[ip] = t;
    out_I[ip] = t*t;
}
// Coeff convention: [px, py, pz, s] per atom (float4)
// basis_data: packed as float2 per node (wf, wf_spline second derivative)
// ============================================================================
__kernel void project_orbital_points(
    const int n_points,
    __global const float4* points,    // [n_points] xyz
    __global const AtomData* atoms,   // [natoms]
    const int natoms,
    __global const float* coeffs,     // [natoms*4] packed as float4
    __global const float* basis_data,
    const int n_nodes,
    const float dr_basis,
    const int max_shells,
    __global float* out_psi           // [n_points]
) {
    const int ip = get_global_id(0);
    if (ip >= n_points) return;

    const float3 p = points[ip].xyz;
    float psi = 0.0f;

    for (int ia = 0; ia < natoms; ia++) {
        const AtomData ad = atoms[ia];
        float3 d = p - ad.pos_rcut.xyz;
        const float r2 = dot(d, d);
        const float rcut2 = ad.pos_rcut.w * ad.pos_rcut.w;
        if (r2 > rcut2) continue;
        const float r = sqrt(r2);

        const float rs = evaluate_radial(r, ad.type, 0, basis_data, n_nodes, dr_basis, max_shells);
        const float rp = evaluate_radial(r, ad.type, 1, basis_data, n_nodes, dr_basis, max_shells);
        const float invr = 1.0f / (r + 1e-12f);
        const float3 rhat = d * invr;

        const float4 basis_val = (float4)(
            rp * rhat.x * PREF_P,
            rp * rhat.y * PREF_P,
            rp * rhat.z * PREF_P,
            rs * PREF_S
        );

        const int coeff_base = ia * 4;
        const __global float4* coeffs_ptr = (const __global float4*)(coeffs);
        const float4 c = coeffs_ptr[coeff_base / 4];
        psi += dot(c, basis_val);
    }

    out_psi[ip] = psi;
}



// ============================================================================
// Orbital projection at arbitrary points using dense MO coefficient vector.
// Computes psi(p) = sum_{atoms} sum_{mu} coeffs[i0orb+mu] * phi_mu(p)
// Supports arbitrary number of orbitals per atom (s, p, d) via i0orb/norb.
// ============================================================================
__kernel void project_orbital_dense_points(
    const int n_points,
    __global const float4* points,
    __global const AtomData* atoms,
    const int natoms,
    __global const float* coeffs,       // [norb_total] dense coefficient vector
    __global const float* basis_data,
    const int n_nodes,
    const float dr_basis,
    const int max_shells,
    __global float* out_psi
) {
    const int ip = get_global_id(0);
    if (ip >= n_points) return;

    const float3 p = points[ip].xyz;
    float psi = 0.0f;

    for (int ia = 0; ia < natoms; ++ia) {
        AtomData ad = atoms[ia];
        float3 d = p - ad.pos_rcut.xyz;
        float r2 = dot(d, d);
        float rcut2 = ad.pos_rcut.w * ad.pos_rcut.w;
        if (r2 > rcut2) continue;

        float r = sqrt(r2);
        float3 rhat = d / (r + 1e-12f);

        int i0 = ad.i0orb;
        int norb = ad.norb;
        int iorb = 0;
        for (int ish = 0; ish < max_shells && iorb < norb; ++ish) {
            float R = evaluate_radial(r, ad.type, ish, basis_data, n_nodes, dr_basis, max_shells);
            int l = ish;
            for (int mm = -l; mm <= l && iorb < norb; ++mm, ++iorb) {
                float ang = eval_angular_dense(l, mm, rhat);
                psi += coeffs[i0 + iorb] * R * ang;
            }
        }
    }

    out_psi[ip] = psi;
}

// ============================================================================
// Orbital projection at arbitrary points using dense MO coefficient vector
// with exponential radial decay (for STM at large distances).
// Computes psi(p) = sum_{atoms} sum_{mu} coeffs[i0orb+mu] * phi_mu(p)
// Uses f(r) = exp(-beta*(r - r0)) instead of spline basis.
// NO CUTOFF - truly long-range for STM simulation.
// Supports arbitrary number of orbitals per atom (s, p, d) via i0orb/norb.
// ============================================================================
__kernel void project_orbital_dense_points_exp(
    const int n_points,
    __global const float4* points,
    __global const AtomData* atoms,
    const int natoms,
    __global const float* coeffs,       // [norb_total] dense coefficient vector
    const float beta,                   // exponential decay constant (Å^-1), must be > 0
    const float r0,                     // reference distance (Å) where f=1
    const int max_shells,
    __global float* out_psi
) {
    const int ip = get_global_id(0);
    if (ip >= n_points) return;

    const float3 p = points[ip].xyz;
    float psi = 0.0f;

    for (int ia = 0; ia < natoms; ++ia) {
        AtomData ad = atoms[ia];
        float3 d = p - ad.pos_rcut.xyz;
        float r = sqrt(dot(d, d));
        float3 rhat = d / (r + 1e-12f);

        // Exponential radial decay: f(r) = exp(-beta*(r - r0))
        // beta > 0 ensures decaying function
        float R = exp(-beta * (r - r0));

        int i0 = ad.i0orb;
        int norb = ad.norb;
        int iorb = 0;
        for (int ish = 0; ish < max_shells && iorb < norb; ++ish) {
            int l = ish;
            for (int mm = -l; mm <= l && iorb < norb; ++mm, ++iorb) {
                float ang = eval_angular_dense(l, mm, rhat);
                psi += coeffs[i0 + iorb] * R * ang;
            }
        }
    }

    out_psi[ip] = psi;
}

// ============================================================================
// Density projection at arbitrary points using dense density matrix.
// Computes rho(p) = sum_{i,j} sum_{mu,nu} coeffs[i0_i+mu] * coeffs[i0_j+nu] * phi_mu(p) * phi_nu(p)
// For diagonal density matrix (MO occupancy): dm is a 1D vector of diagonal elements
// For full density matrix: dm is a 2D matrix flattened as [norb_total][norb_total]
// ============================================================================
__kernel void project_density_dense_points(
    const int n_points,
    __global const float4* points,
    __global const AtomData* atoms,
    const int natoms,
    __global const float* dm,          // dense density matrix [norb_total*norb_total]
    const int norb_total,
    __global const float* basis_data,
    const int n_nodes,
    const float dr_basis,
    const int max_shells,
    __global float* out_rho
) {
    const int ip = get_global_id(0);
    if (ip >= n_points) return;

    const float3 p = points[ip].xyz;

    // Evaluate all orbitals at this point into private arrays
    // Max orbitals per atom = 9 (spd), max atoms in a typical task = 64
    // We use a fixed-size buffer for simplicity; for very large systems,
    // this kernel should be called with smaller point batches.
    float phi_i[9];  // max orbitals per atom (s=1, p=3, d=5)
    float phi_j[9];

    float density = 0.0f;

    for (int ia = 0; ia < natoms; ++ia) {
        AtomData ad_i = atoms[ia];
        int i0_i = ad_i.i0orb;
        int norb_i = eval_atom_orbitals(p, ad_i, basis_data, n_nodes, dr_basis, max_shells, phi_i);
        if (norb_i == 0) continue;

        for (int ja = ia; ja < natoms; ++ja) {
            AtomData ad_j = atoms[ja];
            int i0_j = ad_j.i0orb;
            int norb_j = eval_atom_orbitals(p, ad_j, basis_data, n_nodes, dr_basis, max_shells, phi_j);
            if (norb_j == 0) continue;

            float pairsym = (ia == ja) ? 1.0f : 2.0f;

            for (int mu = 0; mu < norb_i; ++mu) {
                for (int nu = 0; nu < norb_j; ++nu) {
                    float rho_munu = dm[(i0_i + mu) * norb_total + (i0_j + nu)];
                    density += pairsym * rho_munu * phi_i[mu] * phi_j[nu];
                }
            }
        }
    }

    out_rho[ip] = density;
}

// ============================================================================
// Orbital projection onto 3D grid using dense MO coefficient vector.
// Same execution model as project_orbital: one workgroup per 8x8x8 task block.
// ============================================================================
__kernel void project_orbital_dense(
    __global const GridSpec* grid,
    const int n_tasks,
    __global const TaskData* tasks,
    __global const AtomData* atoms,
    __global const int* task_atoms,
    __global const float* coeffs,       // [norb_total] dense coefficient vector
    __global const float* basis_data,
    const int n_nodes,
    const float dr_basis,
    const int max_shells,
    const int nMaxAtom,
    __global float* out_grid
) {
    const int i_task = get_group_id(0);
    const int t_idx  = get_local_id(0);
    const int threads_per_task = get_local_size(0);

    if (i_task >= n_tasks) return;

    const TaskData task = tasks[i_task];
    const int na = task.na;

    for (int v = t_idx; v < 512; v += threads_per_task) {
        const int lx = v & 7;
        const int ly = (v >> 3) & 7;
        const int lz = (v >> 6) & 7;
        const int gx = task.x * 8 + lx;
        const int gy = task.y * 8 + ly;
        const int gz = task.z * 8 + lz;
        const int3 ngrid_dim = grid->ngrid.xyz;
        if (gx >= ngrid_dim.x || gy >= ngrid_dim.y || gz >= ngrid_dim.z) continue;
        const int g_idx = (gx * ngrid_dim.y + gy) * ngrid_dim.z + gz;
        float3 r_vox = grid->origin.xyz + (float)gx * grid->dA.xyz + (float)gy * grid->dB.xyz + (float)gz * grid->dC.xyz;

        float psi = 0.0f;
        for (int i = 0; i < na; ++i) {
            const int i_atom = task_atoms[i_task * nMaxAtom + i];
            AtomData ad = atoms[i_atom];
            float3 d = r_vox - ad.pos_rcut.xyz;
            float r2 = dot(d, d);
            float rcut2 = ad.pos_rcut.w * ad.pos_rcut.w;
            if (r2 > rcut2) continue;
            float r = sqrt(r2);
            float3 rhat = d / (r + 1e-12f);

            int i0 = ad.i0orb;
            int norb = ad.norb;
            int iorb = 0;
            for (int ish = 0; ish < max_shells && iorb < norb; ++ish) {
                float R = evaluate_radial(r, ad.type, ish, basis_data, n_nodes, dr_basis, max_shells);
                int l = ish;
                for (int mm = -l; mm <= l && iorb < norb; ++mm, ++iorb) {
                    float ang = eval_angular_dense(l, mm, rhat);
                    psi += coeffs[i0 + iorb] * R * ang;
                }
            }
        }
        out_grid[g_idx] = psi;
    }
}

// ============================================================================
// Density projection onto 3D grid using dense density matrix.
// Same execution model as project_density_sparse: one workgroup per 8x8x8 task block.
// ============================================================================
__kernel void project_density_dense(
    __global const GridSpec* grid,
    const int n_tasks,
    __global const TaskData* tasks,
    __global const AtomData* atoms,
    __global const int* task_atoms,
    __global const float* dm,          // [norb_total*norb_total] dense density matrix
    const int norb_total,
    __global const float* basis_data,
    const int n_nodes,
    const float dr_basis,
    const int max_shells,
    const int nMaxAtom,
    __global float* out_grid
) {
    const int i_task = get_group_id(0);
    const int t_idx  = get_local_id(0);
    const int threads_per_task = get_local_size(0);

    if (i_task >= n_tasks) return;

    const TaskData task = tasks[i_task];
    const int na = task.na;

    for (int v = t_idx; v < 512; v += threads_per_task) {
        const int lx = v & 7;
        const int ly = (v >> 3) & 7;
        const int lz = (v >> 6) & 7;
        const int gx = task.x * 8 + lx;
        const int gy = task.y * 8 + ly;
        const int gz = task.z * 8 + lz;
        const int3 ngrid_dim = grid->ngrid.xyz;
        if (gx >= ngrid_dim.x || gy >= ngrid_dim.y || gz >= ngrid_dim.z) continue;
        const int g_idx = (gx * ngrid_dim.y + gy) * ngrid_dim.z + gz;
        float3 r_vox = grid->origin.xyz + (float)gx * grid->dA.xyz + (float)gy * grid->dB.xyz + (float)gz * grid->dC.xyz;

        float phi_i[9];  // max orbitals per atom (s=1, p=3, d=5)
        float phi_j[9];
        float density = 0.0f;

        for (int i = 0; i < na; ++i) {
            const int i_atom = task_atoms[i_task * nMaxAtom + i];
            AtomData ad_i = atoms[i_atom];
            int i0_i = ad_i.i0orb;
            int norb_i = eval_atom_orbitals(r_vox, ad_i, basis_data, n_nodes, dr_basis, max_shells, phi_i);
            if (norb_i == 0) continue;

            for (int j = i; j < na; ++j) {
                const int j_atom = task_atoms[i_task * nMaxAtom + j];
                AtomData ad_j = atoms[j_atom];
                int i0_j = ad_j.i0orb;
                int norb_j = eval_atom_orbitals(r_vox, ad_j, basis_data, n_nodes, dr_basis, max_shells, phi_j);
                if (norb_j == 0) continue;

                float pairsym = (i_atom == j_atom) ? 1.0f : 2.0f;
                for (int mu = 0; mu < norb_i; ++mu) {
                    for (int nu = 0; nu < norb_j; ++nu) {
                        float rho_munu = dm[(i0_i + mu) * norb_total + (i0_j + nu)];
                        density += pairsym * rho_munu * phi_i[mu] * phi_j[nu];
                    }
                }
            }
        }
        out_grid[g_idx] = density;
    }
}






// ============================================================================
// STM Response Amplitude — exponential-basis coupling, single s-tip orbital
// ============================================================================
// Precompute on CPU:  G0 = inv((E+iη)S_s - H_s),  v = C^T G0
// GPU kernel builds a_st = (E+iη)S_ts - H_ts per grid point and computes:
//   resp = |v·a_st^H|^2 / |(E+iη-E_tip) - a_st·G0·a_st^H|^2
//
// Buffers:
//   points      [n_points]   float4 xyz (tip positions)
//   atoms_s     [natoms_s]   AtomData for sample atoms
//   starts_s    [natoms_s+1] int   orbital offsets
//   v_re, v_im  [ns]         float  precomputed v = C^T G0
//   G0_re, G0_im[ns*ns]      float  precomputed sample Green's function
//   params: E_re, E_im, E_tip, beta, r0, A_ss, A_sp, rcut
//
// out_resp      [n_points]   float  response amplitude
// ----------------------------------------------------------------------------
__kernel void response_amplitude_exp(
    const int n_points,
    __global const float4* points,
    const int natoms_s,
    __global const AtomData* atoms_s,
    __global const int* starts_s,
    const int ns,
    __global const float* v_re,
    __global const float* v_im,
    __global const float* G0_re,
    __global const float* G0_im,
    const float E_re,
    const float E_im,
    const float E_tip,
    const float beta,
    const float r0,
    const float A_ss,
    const float A_sp,
    const float rcut,
    __global float* out_resp
) {
    const int ip = get_global_id(0);
    if (ip >= n_points) return;

    const float3 p = points[ip].xyz;
    const float rcut2 = rcut * rcut;

    // Build a_st in private memory (max 256 orbitals)
    float2 a_st[256];
    for (int i = 0; i < ns; i++) { a_st[i] = (float2)(0.0f, 0.0f); }

    for (int ia = 0; ia < natoms_s; ia++) {
        const AtomData ad = atoms_s[ia];
        float3 d = p - ad.pos_rcut.xyz;
        const float r2 = dot(d, d);
        if (r2 > rcut2 || r2 < 1e-16f) continue;
        const float r = sqrt(r2);
        const float invr = 1.0f / r;
        const float3 rhat = d * invr;
        const float l = rhat.x, m = rhat.y, n = rhat.z;
        const float f = exp(-beta * (r - r0));

        const float Vss = A_ss * f;
        const float Vsp = A_sp * f;

        const int i0 = starts_s[ia];
        const int nj = starts_s[ia+1] - i0;

        // a = z*S - H; with S=0 for tunneling => a = -H
        // s-tip row: [Vss, l*Vsp, m*Vsp, n*Vsp]
        a_st[i0] = (float2)(-Vss, 0.0f);
        if (nj > 1) {
            a_st[i0+1] = (float2)(-l * Vsp, 0.0f);
            a_st[i0+2] = (float2)(-m * Vsp, 0.0f);
            a_st[i0+3] = (float2)(-n * Vsp, 0.0f);
        }
    }

    // s1 = sum_i v_i * conj(a_i) = sum_i (v_re_i * a_re_i) + i * sum_i (-v_im_i * a_re_i)
    float s1_re = 0.0f, s1_im = 0.0f;
    for (int i = 0; i < ns; i++) {
        s1_re += v_re[i] * a_st[i].x;
        s1_im += -v_im[i] * a_st[i].x;
    }

    // s2 = sum_i a_i * sum_j G0_ij * conj(a_j)
    // conj(a_j) = a_j since a_im = 0
    // b_i = sum_j (G0_re_ij + i*G0_im_ij) * a_j
    float s2_re = 0.0f, s2_im = 0.0f;
    for (int i = 0; i < ns; i++) {
        float b_re = 0.0f, b_im = 0.0f;
        for (int j = 0; j < ns; j++) {
            const int ij = i * ns + j;
            b_re += G0_re[ij] * a_st[j].x;
            b_im += G0_im[ij] * a_st[j].x;
        }
        s2_re += a_st[i].x * b_re;
        s2_im += a_st[i].x * b_im;
    }

    const float d_re = (E_re - E_tip) - s2_re;
    const float d_im = E_im - s2_im;
    const float d_norm2 = d_re * d_re + d_im * d_im;
    const float s1_norm2 = s1_re * s1_re + s1_im * s1_im;

    out_resp[ip] = (d_norm2 > 1e-30f) ? (s1_norm2 / d_norm2) : 0.0f;
}







// ====================


/*

### 1. Analysis of the Constraints & The "Dyson Subspace" Trick

**The Memory Constraint:** 
A dense $400 \times 400$ complex matrix (`float2`) takes about **1.28 MB**. The absolute maximum shared memory (`__local`) per workgroup on most GPUs is **32 KB to 64 KB**.
*Conclusion:* You **cannot** load the full matrices into shared memory, and running a $400 \times 400$ Gauss-Jordan elimination in global memory for every single pixel would be far too slow.

**The Physical Optimization (The "Active Subspace"):**
The hopping matrix $H_{TS}(R)$ between the tip and the sample is mostly zeros! It only has non-zero entries for the atoms physically close to each other (e.g., the tip apex and the $\approx 1-4$ sample atoms directly beneath it). 
If we use a cutoff radius $R_{cut}$, the number of "active" orbitals $N_{act}$ is very small. For example, 4 tip atoms and 4 sample atoms = $16 \times 16$ active orbitals. 

Instead of solving the full $400 \times 400$ system $Ax=b$, we use the **Dyson Equation projection**:
1.  **Host/Global:** Precalculate the full isolated Green's functions $G_T$ (Tip) and $G_S$ (Sample). Store them in `__global` memory.
2.  **Kernel/Local:** Identify the small active subsets.
3.  **Kernel/Local:** Extract only the small active blocks of $G_T$ and $G_S$ into `__local` memory (e.g., $16 \times 16$ complex = **2 KB**, easily fitting in shared memory!).
4.  **Kernel/Local:** Solve the multiple-scattering transmission matrix strictly inside this active subspace using a local **Gauss-Jordan solver**.

---

### 2. The Algorithm per Workgroup (Pixel)

Each workgroup handles **1 pixel (1 tip position)**. The threads (e.g., 64 threads) collaborate:
1.  **Distance Filter:** Threads loop over atoms, measure distances, and build a list of active tip and sample atom indices.
2.  **Preload (Tiling):** Threads collaboratively copy the required sub-blocks of $G_S$ and $G_T$ from global to local memory.
3.  **Build $V_{TS}$:** Threads compute the $4 \times 4$ Slater-Koster overlap blocks for the active pairs.
4.  **Matrix Multiplications:** Compute $W = I - G_S V_{TS}^\dagger G_T V_{TS}$ in local memory.
5.  **Linear Solver:** Run an in-place parallel Gauss-Jordan elimination on $W$ to solve $W x_{S} = b_S$.
6.  **Current Integration:** Calculate the trace/dot product representing the transmission.

*/

// Complex arithmetic helpers
inline float2 c_add(float2 a, float2 b) { return (float2)(a.x+b.x, a.y+b.y); }
inline float2 c_sub(float2 a, float2 b) { return (float2)(a.x-b.x, a.y-b.y); }
inline float2 c_mul(float2 a, float2 b) { return (float2)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }
inline float2 c_div(float2 a, float2 b) {
    float den = b.x*b.x + b.y*b.y + 1e-30f;
    return (float2)((a.x*b.x + a.y*b.y)/den, (a.y*b.x - a.x*b.y)/den);
}

// Define the maximum active subspace size. 
// 32 orbitals = 8 atoms (4 orbitals per atom: px, py, pz, s)
// Fits perfectly in shared memory (32x32 complex matrix = 8 KB)
#define MAX_ACT_ORB 32 

__kernel void solve_stm_dyson_wg(
    const int n_pixels,
    __global const float4* tip_centers,
    __global const float4* tip_pos_rel,
    __global const float4* smp_pos,
    const int ntip_atoms,
    const int nsmp_atoms,
    // Precalculated full Green's functions from Host (Global Memory)
    __global const float2* GT_global, // [4*ntip_atoms * 4*ntip_atoms]
    __global const float2* GS_global, // [4*nsmp_atoms * 4*nsmp_atoms]
    // Incident wave vector injected from the source lead into the tip
    __global const float2* uT_source, // [4*ntip_atoms]
    // Hopping parameters
    const float beta, const float r0, const float rcut,
    // Output
    __global float* out_current
) {
    // 1 Workgroup = 1 Pixel
    const int pixel_id = get_group_id(0);
    const int t_idx    = get_local_id(0); // Thread ID (e.g., 0 to 63)
    const int threads  = get_local_size(0);
    
    if (pixel_id >= n_pixels) return;

    // --- SHARED MEMORY ALLOCATIONS ---
    __local int active_T_atoms[8]; // Max 8 active tip atoms
    __local int active_S_atoms[8]; // Max 8 active sample atoms
    __local int num_act_T, num_act_S;

    // Local Matrices (32x32 floats2 = 8KB each)
    __local float2 GS_loc[MAX_ACT_ORB][MAX_ACT_ORB];
    __local float2 GT_loc[MAX_ACT_ORB][MAX_ACT_ORB];
    __local float2 V_ts[MAX_ACT_ORB][MAX_ACT_ORB]; // Hopping
    __local float2 W[MAX_ACT_ORB][MAX_ACT_ORB]; // The Dyson Matrix to invert
    
    // Local Vectors
    __local float2 uT_loc[MAX_ACT_ORB];
    __local float2 bS_loc[MAX_ACT_ORB]; // Right-hand side

    if (t_idx == 0) { num_act_T = 0; num_act_S = 0; }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 1. DYNAMICALLY IDENTIFY ACTIVE ATOMS (Distance < rcut)
    const float3 cen   = tip_centers[pixel_id].xyz;
    const float rcut2  = rcut * rcut;

    // Tip: take first up to 8 atoms (this is a simple, deterministic choice)
    if (t_idx == 0) {
        const int nt = min(ntip_atoms, 8);
        num_act_T = nt;
        for (int i = 0; i < nt; i++) active_T_atoms[i] = i;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Sample: collect up to 8 atoms within rcut from the tip center
    for (int ja = t_idx; ja < nsmp_atoms; ja += threads) {
        const float3 ps = smp_pos[ja].xyz;
        const float3 d  = ps - cen;
        const float  r2 = dot(d, d);
        if (r2 < rcut2) {
            const int slot = atomic_inc((volatile __local int*)&num_act_S);
            if (slot < 8) active_S_atoms[slot] = ja;
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (t_idx == 0) {
        if (num_act_S > 8) num_act_S = 8;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    const int N_T = min(num_act_T * 4, MAX_ACT_ORB);
    const int N_S = min(num_act_S * 4, MAX_ACT_ORB);
    if ((N_T <= 0) || (N_S <= 0)) {
        if (t_idx == 0) out_current[pixel_id] = 0.0f;
        return;
    }

    // 2. PRELOAD ACTIVE BLOCKS FROM GLOBAL TO SHARED MEMORY
    // Threads loop through the 2D local arrays and fetch from global memory
    for(int i = t_idx; i < N_S * N_S; i += threads) {
        int r = i / N_S; int c = i % N_S;
        int glob_r = active_S_atoms[r/4]*4 + (r%4);
        int glob_c = active_S_atoms[c/4]*4 + (c%4);
        GS_loc[r][c] = GS_global[glob_r * (4*nsmp_atoms) + glob_c];
    }
    for(int i = t_idx; i < N_T * N_T; i += threads) {
        int r = i / N_T; int c = i % N_T;
        int glob_r = active_T_atoms[r/4]*4 + (r%4);
        int glob_c = active_T_atoms[c/4]*4 + (c%4);
        GT_loc[r][c] = GT_global[glob_r * (4*ntip_atoms) + glob_c];
    }
    for(int i = t_idx; i < N_T; i += threads) {
        int glob_i = active_T_atoms[i/4]*4 + (i%4);
        uT_loc[i] = uT_source[glob_i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 3. COMPUTE SLATER-KOSTER HOPPING MATRIX V_ts (real, stored as complex with imag=0)
    for (int i = t_idx; i < N_T * N_S; i += threads) {
        const int it = i / N_S;
        const int is = i - it * N_S;
        const int ia = active_T_atoms[it / 4];
        const int ja = active_S_atoms[is / 4];
        const int ot = it & 3;
        const int os = is & 3;
        const float3 pt = cen + tip_pos_rel[ia].xyz;
        const float3 ps = smp_pos[ja].xyz;
        float3 d = ps - pt;
        const float r2 = dot(d, d);
        if (r2 > rcut2 || r2 < 1e-16f) {
            V_ts[it][is] = (float2)(0.0f, 0.0f);
            continue;
        }
        const float r = sqrt(r2);
        const float invr = 1.0f / r;
        const float l = d.x * invr;
        const float m = d.y * invr;
        const float n = d.z * invr;
        const float f = exp(-beta * (r - r0));

        // Use simple isotropic exponential SK parameters (units are arbitrary scaling for now)
        const float Vss     = 1.0f * f;
        const float Vsp     = 1.0f * f;
        const float Vps     = 1.0f * f;
        const float Vpp_sig = 1.0f * f;
        const float Vpp_pi  = 0.2f * f;

        // Orbital order in this kernel: (px,py,pz,s) index 0..3
        float val = 0.0f;
        if (ot == 3 && os == 3) {
            val = Vss;
        } else if (ot == 3 && os != 3) {
            val = Vsp * ((os == 0) ? l : (os == 1) ? m : n);
        } else if (ot != 3 && os == 3) {
            val = Vps * ((ot == 0) ? l : (ot == 1) ? m : n);
        } else {
            const float ut = (ot == 0) ? l : (ot == 1) ? m : n;
            const float us = (os == 0) ? l : (os == 1) ? m : n;
            const float dV = Vpp_sig - Vpp_pi;
            const float delta = (ot == os) ? 1.0f : 0.0f;
            val = Vpp_pi * delta + dV * ut * us;
        }
        V_ts[it][is] = (float2)(val, 0.0f);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 3b) M1 = GT * V_ts   (N_T x N_S)
    __local float2 M1[MAX_ACT_ORB][MAX_ACT_ORB];
    for (int i = t_idx; i < N_T * N_S; i += threads) {
        const int it = i / N_S;
        const int is = i - it * N_S;
        float2 acc = (float2)(0.0f, 0.0f);
        for (int kt = 0; kt < N_T; kt++) {
            acc = c_add(acc, c_mul(GT_loc[it][kt], V_ts[kt][is]));
        }
        M1[it][is] = acc;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 4. Construct W = I - GS * (V^H * (GT*V)) = I - GS * M2
    for (int i = t_idx; i < N_S * N_S; i += threads) {
        const int is = i / N_S;
        const int js = i - is * N_S;

        // M2[ks,js] = sum_t conj(V[t,ks]) * M1[t,js]
        float2 m3 = (float2)(0.0f, 0.0f);
        for (int ks = 0; ks < N_S; ks++) {
            float2 m2 = (float2)(0.0f, 0.0f);
            for (int it = 0; it < N_T; it++) {
                const float2 v = V_ts[it][ks];
                const float2 vH = (float2)(v.x, -v.y);
                m2 = c_add(m2, c_mul(vH, M1[it][js]));
            }
            m3 = c_add(m3, c_mul(GS_loc[is][ks], m2));
        }

        const float2 Iij = (is == js) ? (float2)(1.0f, 0.0f) : (float2)(0.0f, 0.0f);
        W[is][js] = c_sub(Iij, m3);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // RHS: bS = GS * (V^H * (GT * uT))
    __local float2 tvec[MAX_ACT_ORB];
    __local float2 svec[MAX_ACT_ORB];
    for (int it = t_idx; it < N_T; it += threads) {
        float2 acc = (float2)(0.0f, 0.0f);
        for (int kt = 0; kt < N_T; kt++) {
            acc = c_add(acc, c_mul(GT_loc[it][kt], uT_loc[kt]));
        }
        tvec[it] = acc;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int is = t_idx; is < N_S; is += threads) {
        float2 acc = (float2)(0.0f, 0.0f);
        for (int it = 0; it < N_T; it++) {
            const float2 v = V_ts[it][is];
            const float2 vH = (float2)(v.x, -v.y);
            acc = c_add(acc, c_mul(vH, tvec[it]));
        }
        svec[is] = acc;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int is = t_idx; is < N_S; is += threads) {
        float2 acc = (float2)(0.0f, 0.0f);
        for (int ks = 0; ks < N_S; ks++) {
            acc = c_add(acc, c_mul(GS_loc[is][ks], svec[ks]));
        }
        bS_loc[is] = acc;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 5. IN-PLACE PARALLEL GAUSS-JORDAN SOLVER (W * xS = bS)
    // We solve the linear system for the active sample subspace
    for (int k = 0; k < N_S; k++) {
        // Thread 0 finds pivot (partial pivoting optional for N<32, but good for stability)
        if (t_idx == 0) {
            float2 pivot = W[k][k];
            // Normalize pivot row
            for(int j=k; j < N_S; j++) W[k][j] = c_div(W[k][j], pivot);
            bS_loc[k] = c_div(bS_loc[k], pivot);
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // All threads eliminate the other rows
        for (int i = t_idx; i < N_S; i += threads) {
            if (i != k) {
                float2 factor = W[i][k];
                for (int j = k; j < N_S; j++) {
                    W[i][j] = c_sub(W[i][j], c_mul(factor, W[k][j]));
                }
                bS_loc[i] = c_sub(bS_loc[i], c_mul(factor, bS_loc[k]));
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // Now bS_loc contains the exact response wave x_S.

    // 6. CALCULATE FINAL TRANSMISSION CURRENT
    // I = x_S^H * Gamma_R * x_S
    float current = 0.0f;
    if (t_idx == 0) {
        for(int i=0; i<N_S; i++) {
            // Assuming Gamma_R is a simple density-of-states weighting for the drain
            current += (bS_loc[i].x*bS_loc[i].x + bS_loc[i].y*bS_loc[i].y);
        }
        out_current[pixel_id] = current;
    }
}

// ======================================================================
// STM GF Dyson 2-molecule MO scan — full-matrix GF approach (work-item per pixel)
//
// Math: amp = c_tip^H · GT · M_ts · GS · c_smp
// Precompute on CPU:
//   v_S = GS @ c_smp        (smp_norb_ocl vector, remapped to OCL [px,py,pz,s] order)
//   u_T = c_tip^H @ GT      (tip_norb_ocl vector, remapped to OCL [px,py,pz,s] order)
//
// On GPU per pixel:
//   amp = Σ_{it∈tip, is∈smp} u_T[it] · H_hop(it,is) · v_S[is]
//   out = |amp|²
//
// H_hop is computed via simplified exponential Slater-Koster:
//   f = exp(-beta*(r-r0))
//   SK with Vss=1, Vsp=1, Vps=1, Vpp_sig=1, Vpp_pi=0.2 all multiplied by f
//
// Orbital convention in this kernel: OpenCL Grid order [px,py,pz,s] per atom.
// All vectors (u_T, v_S) and orb2atom arrays are already remapped on CPU.
// ======================================================================
__kernel void stm_gf_dyson_2mol_mo_scan(
    const int n_points,
    __global const float4* tip_centers,
    __global const float4* tip_pos_rel,
    __global const float4* smp_pos,
    __global const int* tip_orb2atom,     // [tip_norb_ocl] 0-based atom index per orbital
    __global const int* smp_orb2atom,     // [smp_norb_ocl]
    __global const float2* u_T,            // [tip_norb_ocl] c_tip^H @ GT  in OCL order
    __global const float2* v_S,            // [smp_norb_ocl] GS @ c_smp   in OCL order
    const int ntip_atoms,
    const int nsmp_atoms,
    const int tip_norb_ocl,                // = ntip_atoms * 4
    const int smp_norb_ocl,                // = nsmp_atoms * 4
    const float beta,
    const float r0,
    const float rcut,
    __global float* out_current
) {
    const int ip = get_global_id(0);
    if (ip >= n_points) return;

    const float3 cen = tip_centers[ip].xyz;
    const float rcut2 = rcut * rcut;
    float2 amp = (float2)(0.0f, 0.0f);

    for (int it = 0; it < tip_norb_ocl; it++) {
        const int ia = tip_orb2atom[it];
        if (ia < 0 || ia >= ntip_atoms) continue;
        const float2 uTit = u_T[it];
        if (uTit.x == 0.0f && uTit.y == 0.0f) continue;  // skip zero-padded orbitals (e.g. H px,py,pz)

        const float3 pt = cen + tip_pos_rel[ia].xyz;
        const int ot = it & 3;  // 0=px, 1=py, 2=pz, 3=s in OCL convention

        for (int is = 0; is < smp_norb_ocl; is++) {
            const int ja = smp_orb2atom[is];
            if (ja < 0 || ja >= nsmp_atoms) continue;
            const float2 vSis = v_S[is];
            if (vSis.x == 0.0f && vSis.y == 0.0f) continue;

            const float3 ps = smp_pos[ja].xyz;
            const float3 d = pt - ps;
            const float r2 = dot(d, d);
            if (r2 > rcut2 || r2 < 1e-16f) continue;

            const float r = sqrt(r2);
            const float invr = 1.0f / r;
            const float l = d.x * invr;
            const float m = d.y * invr;
            const float n = d.z * invr;
            const float f = exp(-beta * (r - r0));
            const int os = is & 3;

            // Simplified exponential SK hopping (real, symmetric)
            float V = 0.0f;
            if (ot == 3 && os == 3) {
                V = f;                                   // Vss
            } else if (ot == 3 && os != 3) {
                V = f * ((os == 0) ? l : (os == 1) ? m : n);  // Vsp * dir
            } else if (ot != 3 && os == 3) {
                V = f * ((ot == 0) ? l : (ot == 1) ? m : n);  // Vps * dir
            } else {
                const float ut = (ot == 0) ? l : (ot == 1) ? m : n;
                const float us = (os == 0) ? l : (os == 1) ? m : n;
                const float Vpp_pi = 0.2f * f;
                const float dV = f - Vpp_pi;              // Vpp_sig - Vpp_pi
                const float delta = (ot == os) ? 1.0f : 0.0f;
                V = Vpp_pi * delta + dV * ut * us;
            }

            // amp += u_T[it] * V * v_S[is]   (V is real, stored as scalar)
            amp.x += uTit.x * V * vSis.x - uTit.y * V * vSis.y;
            amp.y += uTit.x * V * vSis.y + uTit.y * V * vSis.x;
        }
    }
    out_current[ip] = amp.x * amp.x + amp.y * amp.y;
}