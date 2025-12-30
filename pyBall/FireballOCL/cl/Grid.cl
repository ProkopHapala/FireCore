
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
    int atom_start;  // index into active_atoms list
    int n_atoms;     // number of overlapping atoms
    int pair_start;  // index into active_pairs list
    int n_pairs;     // number of active pairs (i,j)
} TaskData;

// Real Spherical Harmonics (s, p, d)
// Ortega order: s, py, pz, px, dxy, dyz, dz2, dxz, dx2-y2
float getYlm(int l, int m, float3 rhat) {
    if (l == 0) return 0.28209479177387814f; // 1/sqrt(4*pi)
    if (l == 1) {
        // Fireball order: py, pz, px (m=-1, 0, 1 mapped to py, pz, px)
        if (m == 0) return 0.4886025119029199f * rhat.y; // py
        if (m == 1) return 0.4886025119029199f * rhat.z; // pz
        if (m == 2) return 0.4886025119029199f * rhat.x; // px
    }
    return 0.0f;
}

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
    
    // Cubic B-spline coefficients (a, b, c, d)
    // For now we only have raw values. Fireball's .wf files contain radial values.
    // If we want cubic interpolation, we need to precompute splines or pass 4 nodes.
    // Let's use 4-point cubic Hermite interpolation for better accuracy than linear.
    
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

typedef struct {
    int nssh;
    int lssh[4]; // max 4 shells for now
} SpeciesInfo;

__kernel void project_density_sparse(
    const GridSpec grid,
    const int n_tasks,
    __global const TaskData* tasks,
    __global const AtomData* atoms,
    __global const int* active_atoms,
    __global const int4* active_pairs, // (iatom, jatom, ineigh, pad)
    __global const float* rho,        // [natoms][neigh_max][numorb_max][numorb_max]
    __global const int* neigh_j,      // [natoms][neigh_max]
    __global const float* basis_data, // [n_species][max_shells][n_nodes]
    __global const SpeciesInfo* species_info,
    const int n_nodes, 
    const float dr_basis,
    const int max_shells,
    const int neigh_max,
    const int numorb_max,
    __global float* out_grid          // [nx][ny][nz]
) {
    int i_task = get_group_id(0);
    if (i_task >= n_tasks) return;

    TaskData task = tasks[i_task];
    int3 b_idx = task.block_idx.xyz;
    
    int3 ngrid_dim = grid.ngrid.xyz;

    int lx = get_local_id(0);
    int ly = get_local_id(1);
    int lz = get_local_id(2);

    int gx = b_idx.x * 8 + lx;
    int gy = b_idx.y * 8 + ly;
    int gz = b_idx.z * 8 + lz;

    if (gx >= ngrid_dim.x || gy >= ngrid_dim.y || gz >= ngrid_dim.z) return;

    float3 r_vox = grid.origin.xyz + (float)gx * grid.dA.xyz + (float)gy * grid.dB.xyz + (float)gz * grid.dC.xyz;

    float den = 0.0f;

    // Loop over active pairs in this block
    for (int ip = 0; ip < task.n_pairs; ip++) {
        int4 pair_info = active_pairs[task.pair_start + ip];
        int i_atom = pair_info.x;
        int j_atom = pair_info.y;
        int ineigh_ij = pair_info.z;

        AtomData ad_i = atoms[i_atom];
        AtomData ad_j = atoms[j_atom];

        float3 dri = r_vox - ad_i.pos_rcut.xyz;
        float3 drj = r_vox - ad_j.pos_rcut.xyz;
        float ri = length(dri);
        float rj = length(drj);

        if (ri < ad_i.pos_rcut.w && rj < ad_j.pos_rcut.w) {
            float3 rhat_i = (ri > 1e-12f) ? dri / ri : (float3)(0.0f, 0.0f, 1.0f);
            float3 rhat_j = (rj > 1e-12f) ? drj / rj : (float3)(0.0f, 0.0f, 1.0f);

            SpeciesInfo sp_i = species_info[ad_i.type];
            SpeciesInfo sp_j = species_info[ad_j.type];

            // rho layout (Fortran order): rho(imu, inu, ineigh, iatom)
            // Flat index: imu + inu*numorb_max + ineigh*numorb_max*numorb_max + iatom*numorb_max*numorb_max*neigh_max
            int rho_base = i_atom * neigh_max * numorb_max * numorb_max + ineigh_ij * numorb_max * numorb_max;
            
            int imu_idx = 0;
            for (int ish_i = 0; ish_i < sp_i.nssh; ish_i++) {
                int l_i = sp_i.lssh[ish_i];
                float R_i = evaluate_radial(ri, ad_i.type, ish_i, basis_data, n_nodes, dr_basis, max_shells);
                
                for (int m_i = 0; m_i < 2 * l_i + 1; m_i++) {
                    float phi_i = R_i * getYlm(l_i, m_i, rhat_i);
                    
                    int inu_idx = 0;
                    for (int ish_j = 0; ish_j < sp_j.nssh; ish_j++) {
                        int l_j = sp_j.lssh[ish_j];
                        float R_j = evaluate_radial(rj, ad_j.type, ish_j, basis_data, n_nodes, dr_basis, max_shells);
                        
                        for (int m_j = 0; m_j < 2 * l_j + 1; m_j++) {
                            float phi_j = R_j * getYlm(l_j, m_j, rhat_j);
                            
                            // rho(imu, inu, ineigh, iatom)
                            float r_val = rho[rho_base + inu_idx * numorb_max + imu_idx];
                            den += r_val * phi_i * phi_j;
                            inu_idx++;
                        }
                    }
                    imu_idx++;
                }
            }
        }
    }

    int g_idx = (gx * ngrid_dim.y + gy) * ngrid_dim.z + gz;
    out_grid[g_idx] = den;
}
