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
    int x, y, z, w;  // block_idx. x,y,z are coordinates, w is padding
    int na;          // number of overlapping atoms
    int nj;          // start of jatom block ( if off-diagonal block )
    int pad1;
    int pad2;
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
    if (get_global_id(0) == 0) {
        printf("GPU: kernel entry: n_tasks=%d vox_per_task=%d n_nodes=%d max_shells=%d neigh_max=%d numorb_max=%d ngrid=(%d,%d,%d)\n",
               n_tasks, 512, n_nodes, max_shells, neigh_max, numorb_max,
               grid->ngrid.x, grid->ngrid.y, grid->ngrid.z);
    }
    const int i_task = get_group_id(0);
    const int t_idx  = get_local_id(0);

    //{ // sanitize local memory space
    if (i_task >= n_tasks) return;
    __global const int* my_atoms = task_atoms + i_task * nMaxAtom;
    const TaskData task = tasks[i_task];
    const int  na    = task.na;
    const int  nj    = task.nj;
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

        float den = 0.0f;
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
            dri.w         = dot(dri.xyz, dri.xyz);
            dri.w = sqrt(dri.w); dri.xyz /= (dri.w + 1e-12f);
            dri*=exp(-dri.w);
            
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
                drj.w   = dot(drj.xyz, drj.xyz);

                // if( t_idx==0 ){
                //     float3 rij = ad_j.pos_rcut.xyz - ad_i.pos_rcut.xyz;
                //     int rho_base = i_atom * neigh_max * numorb_max * numorb_max + ineigh_ij * numorb_max * numorb_max;
                //     if( dot(rij,rij)<(2*rcut_i2) ) printf("GPU task[%3i] rho[%i,%i] %f \n", i_task, i, j, rho[rho_base+0] ); 
                // } 
                
                if (dri.w < rcut_i2 && drj.w < rcut_j2) {
                    //int4 sp_i = species_info[ad_i.type];
                    //int4 sp_j = species_info[ad_j.type];

                    int rho_base = i_atom * neigh_max * numorb_max * numorb_max + ineigh_ij * numorb_max * numorb_max;
                    const __global float4* rho_ij = (const __global  float4*)(rho + rho_base); // <-- GLOBAL READ
                    drj.w = sqrt(drj.w); drj.xyz /= (drj.w + 1e-12f);
                    // -- for now we test with STO basis, so we just use exp(-r), later we will use spline radial basis functions
                    drj*=exp(-drj.w);
                    // from dr -> WF
                    // dri.w    =  evaluate_radial(dri.w, ad_i.type, 0, basis_data, n_nodes, dr_basis, max_shells); // si-orb
                    // dri.xyz *=  evaluate_radial(dri.w, ad_i.type, 1, basis_data, n_nodes, dr_basis, max_shells); // pi-orb
                    // drj.w    =  evaluate_radial(drj.w, ad_j.type, 0, basis_data, n_nodes, dr_basis, max_shells); // sj-orb
                    // drj.xyz *=  evaluate_radial(drj.w, ad_j.type, 1, basis_data, n_nodes, dr_basis, max_shells); // pj-orb
                    // 4x4 block (px,py,pz,s)_i * (px,py,pz,s)_j 
                    // den += dot( dri,  (
                    // rho_ij[0]  * drj.x +     // <-- GLOBAL READ
                    // rho_ij[1]  * drj.y + 
                    // rho_ij[2]  * drj.z + 
                    // rho_ij[3]  * drj.w   ) );  
                    

                    den += dot( dri.wxyz,  (
                    rho_ij[0]  * drj.w +     // <-- GLOBAL READ
                    rho_ij[1]  * drj.x + 
                    rho_ij[2]  * drj.y + 
                    rho_ij[3]  * drj.z   ) );  
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

    TaskData task = tasks[i_task];

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
                        l_rho[k] = ((__global float4*)(rho + rho_base))[orb_idx];
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
                    const float radial_i = exp(-r_i);
                    const float4 psi_i = (float4){dpos_i * (radial_i / (r_i + 1e-12f)), radial_i};

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
                        const float radial_j = exp(-r_j);
                        const float4 psi_j = (float4){dpos_j * (radial_j / (r_j + 1e-12f)), radial_j};

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
