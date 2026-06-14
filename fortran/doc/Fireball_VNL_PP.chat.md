## USER

ok, now look on @assemble_sVNL.f90 and other files in 
/home/prokop/git/FireCore/fortran/ASSEMBLERS
which we need to update, I remember there was something about projectors required for VNL and pseudo-potentials in @Fireball_Hamiltonian_Assembly.md 
@Fireball_equations.md 
can you recall that ?

---

## GPT5-Codex

Here’s a refresher on the VNL/projector flow you asked about:

1. **assemble_sVNL.f90 and friends**  
   - `assemble_sVNL.f90` computes the projector overlaps \(s_{VNL} = \langle \phi_\mu | \Psi_{\text{PP}} \rangle\). Implementation mirrors `doscentrosPP` / `recover_PP` plus `rotatePP` so the resulting 4×4 block follows the Ortega ordering (s, py, pz, px) but is still in the \((\nu,\mu)\) block convention.  
   - After `sVNL` blocks are built, `assemble_2c_PP.f90` (and the related contraction helpers) multiply them by the pseudopotential coefficients `cl_pseudo` from `info.dat`:
     \[
     V_{nl}(\mu,\nu) \mathrel{+}= \sum_c cl(c)\; s_{VNL}(\mu,c)\; s_{VNL}(\nu,c)
     \]
     This is exactly what we mirrored on the Python/OpenCL side.

2. **Docs reminding us of these requirements**  
   - **Fireball_Hamiltonian_Assembly.md** calls out that VNL needs (a) the projector overlaps (interaction=5 tables) and (b) the species-dependent PP coefficients (`cl_pseudo`, loaded by `firecore_export_bcna_table`/`findFdata`). It also emphasizes the block orientation and that the contraction must take place before blocks are accumulated into the final Hamiltonian.  
   - **Fireball_equations.md** summarizes the role of the projectors in the Hamiltonian and mentions that VNL is a quadratic form in those overlaps (consistent with the McWEDA context).

3. **What we need to update**  
   - If you plan edits in `fortran/ASSEMBLERS`, focus on ensuring the projector overlaps (`sVNL`) are exported or reused the same way we do in PyOpenCL (so the new test harness can grab them piecewise).  
   - For the docs, we should add a short note (in both `Fireball_Hamiltonian_Assembly.md` and [Fireball_equations.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/Fireball_equations.md:0:0-0:0)) that the new verify path compares raw `sVNL` blocks and that the final VNL contraction is still outstanding on GPU, so parity currently focuses on the projector stage.

Let me know if you want me to draft those doc updates or review specific assembler routines beyond `assemble_sVNL.f90`.


---

## USER

OK, after reviewing the relevant files do you see some complications about how to efficiently implement VNL and pseudopotentials into 
our pyOpenCL version
@OCL_Hamiltonian.py @hamiltonian.cl ?

you should maybe review other files in 
/home/prokop/git/FireCore/fortran/ASSEMBLERS
which you find relevant for this task 

---

## GPT5-Codex


The tricky parts of bringing full VNL+pseudopotentials into the PyOpenCL path line up with what we saw in the Fortran assemblers:

1. **Two-stage workflow** – Fortran first builds projector overlaps (`assemble_sVNL.f90` → `doscentrosPP`/`recover_PP`/`rotatePP`) and only then contracts them via `assemble_2c_PP.f90` using the `cl_pseudo` vectors from `info.dat`. In PyOpenCL we already replicate the projector stage ([assemble_pp](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:265:4-291:21)), but we currently contract on the CPU (see `_contract_AB` in [OCL_Hamiltonian.py](cci:7://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:0:0-0:0)). To do this fully on GPU we’d either:
   - Upload `cl_pseudo` and the PP-orbital multiplicities, expand them per‐m channel, and add a kernel that performs `vnl += sVNL * diag(cl) * sVNL^T` with the right (ν,μ) orientation; or
   - Keep the contraction on CPU but ensure the Python path feeds the results straight into [assemble_full](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:462:4-713:33) so tests can inspect them (what’s missing now is an easy way to export the post-contraction blocks per neighbor).

2. **Species-specific PP data** – `cl_pseudo` lives in the vnl `.dat` headers and depends on the number of PP shells (`lsshPP`) per species. We already parse and cache it in [OCL_Hamiltonian.prepare_splines](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:44:4-114:116), but the OpenCL kernels don’t see it; for an all-GPU solution we’d need to flatten those coefficients into buffers keyed by species. Also make sure `muPP/nuPP` mapping (Fortran’s `make_munuPP`) is available on device—right now we build those maps on the host.

3. **Block orientation and neighbor bookkeeping** – Fortran stores sVNL blocks as (ν,μ) and later transposes when copying to dense matrices. Our kernels do the same, but any GPU-side contraction must respect that (a mistake here was what caused the earlier mismatch). We also have to accumulate both on-site and off-diagonal contributions exactly like Fortran (sum over PP center i and j separately) which may mean evaluating more pair combinations (`(i,i)`, `(j,i)`, etc.) in OpenCL.

4. **Kernel plumbing complexity** – [hamiltonian.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/FireballOCL/cl/hamiltonian.cl:0:0-0:0) currently has only [assemble_2c](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:238:4-263:21), [assemble_pp](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:265:4-291:21), and the 3c interpolator. Adding a contraction kernel means keeping `sVNL` results in device memory across multiple calls or recomputing them per neighbor. Given the small system sizes in verification, recomputing is acceptable, but for production we’d want a more efficient caching strategy.

5. **Testing/visibility** – Our new [verify_C2.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/verify_C2.py:0:0-0:0) compares raw sVNL blocks but not the final Vnl contribution. To debug the contraction stage we’ll need either: (a) an export from Fortran after `assemble_2c_PP` similar to `firecore_get_HS_sparse`, or (b) an OpenCL pathway that exposes the contracted Vnl per neighbor so we can compare dense matrices term-by-term.

So the immediate next steps if we want parity are:

- Keep the current host-side contraction but plumb its output into [assemble_full(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:462:4-713:33) so the verification harness can compare contracted Vnl blocks.
- Alternatively, add a simple OpenCL kernel that takes two sVNL blocks (μ×cc, ν×cc) plus the `cl` vector and outputs the contracted 4×4 block, mirroring the Python [contract_AB](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:578:12-593:55). That would let us stay on device while reusing the current host logic for accumulation.

None of this is conceptually hard, but the bookkeeping (species-dependent PP data, block orientation, multiple PP centers per neighbor pair) makes a fully GPU implementation a bit tedious.

---

## USER

So you speak about contraction, can you write it more clearly? What is the contraction doing, and how it is algorithmically from the point of paralleization, memory requirements, local memory usage, kernel desing. I definitely want to desing something which would be efficient at scale. Do we have already some notes on that topic perhaps in @FireballOCL_kernel_opt.md or @Fireball_Review.md ?

I definitely do not want to do in python things which can be done in OpenCL, python is just thing driver/wraper of the OpenCL kernels which do all th job. Only if it is too dofficult we can do it provisionally in python to test if we reproduce frotran but later we definitely need to reproduce it on GPU. We may have small test system now to simplify comparison, but we are planning to use it for huge systems with thousands atoms.

---

## GPT5-Codex

### What “contraction” means here

Fireball (and PyOpenCL) builds the non‑local pseudopotential \(V_{nl}\) in **two stages**:

1. **Projector overlaps (“sVNL”)**  
   For each atom pair \((i,k)\) we evaluate \(\langle \phi_i | \Psi_k \rangle\): a 4×\(n_{PP}\) block (μ rows = fireball orbitals on atom i, columns = PP projectors on atom k). This is what [assemble_pp](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:265:4-291:21) already does on GPU, exactly like other 2‑center Slater–Koster evaluations.

2. **Contraction with PP coefficients**  
   For each neighbor pair \((i,j)\) we need  
   \[
   V_{nl}^{ij} = \sum_{k\in\text{PP centers}} S_{ik}\, \text{diag}(cl_k)\, S_{jk}^{T}
   \]
   where \(cl_k\) are the (expanded) pseudopotential coefficients for atom \(k\). This is just matrix multiply:  \((4\times n_{PP})\cdot(n_{PP}\times n_{PP})\cdot(n_{PP}\times 4)\rightarrow 4\times 4\). Fortran does this in `assemble_2c_PP.f90`.

Right now PyOpenCL computes \(S_{ik}\) on the GPU but does the contraction on the CPU ([contract_AB](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:578:12-593:55)). That’s fine for parity tests but won’t scale.

### GPU-friendly algorithm

**Work item = one (i,j) neighbor block.** Each thread:

1. Loads atom indices (sorted by species pair to keep spline lookups coherent) and direction cosines (reused for all interactions).
2. Fetches / interpolates projector overlaps \(S_{ik}\) and \(S_{jk}\) for every PP center \(k\) that interacts with i or j. (We already have the kernel for this; we just need to keep the results in device memory instead of copying them back.)
3. Multiplies them with the appropriate \(cl_k\) vector and accumulates the outer product into a register `float4x4` (or four float4 rows).

Because each block is only 4×4 and typical \(n_{PP}\) per atom is small (≤15), you can keep the entire contraction in registers without touching shared memory. The data you need per work item is:

- 2× `float4x4` worth of \(S_{ik}\) and \(S_{jk}\) per PP center ⇒ ~128 bytes per PP center (can sit in global memory; load with vector ops).
- `n_{PP}` coefficients for atom k ⇒ load once into registers (or texture cache).
- Accumulator (4 rows × float4) ⇒ 64 bytes in registers.

No local/shared memory “tiling” is required; the arithmetic intensity is high enough that global-memory loads are hidden, especially once pairs are sorted by pair type.

### Parallelization & memory layout

1. **Bag-of-pairs kernel launch** (as in [FireballOCL_kernel_opt.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/FireballOCL_kernel_opt.md:0:0-0:0)):
   - Flatten all neighbors across systems.
   - Sort by pair type.
   - Provide an `output_map` so the kernel writes back in original order.
   - Pass separate spline buffers for S/T/Vna/Vnl if needed; but for contraction we only need the projector buffer.

2. **Projector buffer**  
   - Either recompute \(S_{ik}\) on the fly inside the contraction kernel (costly but simple), or
   - Run a dedicated [assemble_pp](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:265:4-291:21) pass to fill a buffer `svnl_blocks[(i,k)]` and reuse it. For large systems the second approach wins because each PP center participates in many \((i,j)\) pairs.

3. **PP center bookkeeping**  
   - Precompute, on the host, the adjacency “which PP centers k are within cutoff of i” (and same for j). Store as CSR-style lists (`pp_offsets`, `pp_indices`, `pp_counts`).  
   - Provide a mapping `(atom, pp_center_idx) -> buffer index` so the kernel can fetch \(S_{ik}\) quickly; a packed stride or hash table works.

4. **Kernel design** (pseudocode):

```opencl
__kernel void contract_vnl(
    const int n_pairs,
    const __global int2* ij_pairs,
    const __global int* pair_system,   // to locate atoms within batches
    const __global int* pp_offsets,
    const __global int* pp_counts,
    const __global int* pp_indices,    // flattened list of k neighbors per atom
    const __global int* svnl_index,    // map (atom,k) -> offset into svnl buffer
    const __global float* svnl_blocks, // [n_svnl_blocks, 16] row-major (ν,μ)
    const __global float* cl_coeffs,   // expanded PP coeffs per projector
    __global float* vnl_blocks_out     // [n_pairs, 16]
) {
    int gid = get_global_id(0);
    if (gid >= n_pairs) return;

    int2 ij = ij_pairs[gid];
    int i = ij.x;
    int j = ij.y;

    float4 row0 = 0.0f, row1 = 0.0f, row2 = 0.0f, row3 = 0.0f;

    int start = pp_offsets[i];
    int count = pp_counts[i];
    for (int idx = 0; idx < count; ++idx) {
        int k = pp_indices[start + idx];

        int idx_i = svnl_index[i * MAX_PP + idx];
        int idx_j = svnl_index[j * MAX_PP + idx]; // need -1 check

        if (idx_i < 0 || idx_j < 0) continue;

        const __global float4* Sik = (const __global float4*)(svnl_blocks + idx_i*16);
        const __global float4* Sjk = (const __global float4*)(svnl_blocks + idx_j*16);

        // Load projector coeffs for atom k (expanded to per-m entry)
        const __global float* clk = cl_coeffs + k * N_PP_CHANNELS;

        // Multiply: row += Sik * diag(clk) * Sjk^T
        // expand loops explicitly (4 rows, n_pp columns)
    }

    // store row vectors to vnl_blocks_out[gid]
}
```

You can unroll the inner loops for \(n_{PP}\) = 1, 3, 5 etc. depending on species; or keep it generic.

5. **Scaling considerations**
   - **Memory:** storing all \(S_{ik}\) blocks is O(n_atoms × n_pp_neighbors × 16 floats). For typical cutoffs this is manageable (tens of MB for thousands of atoms).
   - **Arithmetic:** each \((i,j)\) pair loops over on the order of 10 PP centers; each contraction is ~4×4×n_pp multiply-adds, easily fits in registers.
   - **Occupancy:** fully parallel over pairs. For 10^4 atoms, neighbor list provides 10^5–10^6 pairs, enough to saturate the GPU.
   - **Local memory:** not needed; rely on coalesced global loads + caches. If you later target older GPUs with small caches, you can tile PP centers into local memory, but it’s optional.

### Notes from docs

- [FireballOCL_kernel_opt.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/FireballOCL_kernel_opt.md:0:0-0:0) already advocates:
  1. **Sort pairs** by species so spline reads are coherent.
  2. **Fused kernels** to avoid multiple global writes.
  3. **GPU Vnl contraction** as a must-have for large systems.

- [Fireball_Review.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/Fireball_Review.md:0:0-0:0) and [Fireball_equations.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/Fireball_equations.md:0:0-0:0) remind us that McWEDA and the PP machinery depend on accurate projector handling; any GPU implementation has to preserve the Ortega ordering and the same rotation matrices. The current OpenCL [assemble_pp](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:265:4-291:21) does this, so we just need to keep its orientation (blocks are (ν,μ); transpose before filling dense matrices) consistent in the contraction kernel.

### Immediate action items

1. **Host prep**
   - Build fast lookup tables: `fast_pair_map` (species → pair_type), PP neighbor CSR lists, `(atom,k)->svnl_index`.
   - Ensure `cl_pseudo` is expanded per projector and uploaded to GPU buffers.

2. **Device**
   - Add [assemble_pp](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/OCL_Hamiltonian.py:265:4-291:21) path that writes `svnl_blocks` to a persistent device buffer (no CPU copy).
   - Implement `contract_vnl` kernel per above, writing 4×4 blocks in the same (ν,μ) orientation as other kernels.

3. **Tests**
   - Update [verify_C2.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/verify_C2.py:0:0-0:0) once the GPU contraction path exists: compare Fortran-contracted Vnl (or full H) against the new OpenCL blocks.
   - Keep the small C2 test for parity; once satisfied, extend to batches to test scaling.

This design keeps Python as the driver only and pushes all heavy math onto the GPU, matching the scalability requirement.


---

## USER

if I want to ask some chatbot, which .f90 files should I give it to give him full context for reimplementation of this part (VNL, pseudopotentials contraction) on GPU / OpenCL. Would you give for context aso some parts of the exsiting .md files ?

---

## GPT5-Codex


To give another assistant enough context to re‑implement the VNL projector + contraction path on OpenCL, I’d bundle these Fortran sources:

1. **`fortran/ASSEMBLERS/assemble_sVNL.f90`** – builds the ⟨ϕ|Ψ⟩ projector blocks (`sVNL`). Shows interpolation, rotation, and block orientation.

2. **`fortran/ASSEMBLERS/assemble_2c_PP.f90`** (and any helper it calls, e.g. `recover_PP`, `doscentrosPP`) – performs the contraction \(V_{nl} = S \cdot \text{diag}(cl) \cdot S^T\) and accumulation into the neighbor blocks. This is the key reference for the CPU algorithm you’re porting.

3. **`fortran/ASSEMBLERS/assemble_2c.f90`** (optional but helpful) – shows how other 2‑center terms are assembled and accumulated; useful for matching conventions and data flow.

4. **[fortran/ROTATIONS/rotate.f90](cci:7://file:///home/prokop/git/FireCore/fortran/ROTATIONS/rotate.f90:0:0-0:0) + [twister.f90](cci:7://file:///home/prokop/git/FireCore/fortran/ROTATIONS/twister.f90:0:0-0:0)** – define the Ortega orbital ordering and rotation matrices that the projectors must follow.

5. **[fortran/MAIN/libFireCore.f90](cci:7://file:///home/prokop/git/FireCore/fortran/MAIN/libFireCore.f90:0:0-0:0)** (selected sections: `firecore_scanHamPiece2c`, `firecore_get_HS_sparse`, data structures for neighbors) – documents how projector data and Hamiltonian blocks are exposed to the outside world.

For written context, I’d also hand over the relevant documentation snippets:

- **`doc/Fireball_Hamiltonian_Assembly.md`** – explains the overall Hamiltonian build, projector workflow, and where VNL fits.
- **[doc/Fireball_equations.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/Fireball_equations.md:0:0-0:0)** – captures the theoretical definitions (McWEDA, projector contractions) so the implementer understands the math target.
- **[doc/FireballOCL_kernel_opt.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/FireballOCL_kernel_opt.md:0:0-0:0)** – already outlines the desired GPU strategy (sorted bag-of-pairs, fused kernels, GPU-side VNL contraction). This is directly relevant to designing the OpenCL kernels.

Providing these files + doc sections will give the chatbot both the exact Fortran reference for algorithms/API and the design goals for the OpenCL port.


---

## USER


I'm trying to reimplement Fireball-DFT code to pyOpenCL, I'm looking on how efficiently write the kernels for pseudopotentials contraction and VNL part, I will give you relevant fotran codes, equations and preliminary discussion

---

## Gemini

Based on the Fortran files and design documents provided, here is the OpenCL implementation design for the VNL (Non-Local Potential) contraction.

### 1. Architectural Strategy

The `assemble_2c_PP.f90` routine calculates the contribution of the pseudopotentials to the Hamiltonian. Unlike generic 3-center integrals, Fireball's implementation (as seen in the code) treats the **2-center approximation** where the projector is located on one of the two atoms involved in the bond pair.

We will map this to GPU using a **"Bag of Pairs"** approach (as requested in `FireballOCL_kernel_opt.md`) to maximize parallelism.

**The Workflow:**
1.  **Input:** `sVNL` buffer (already computed on GPU, containing $\langle \phi_i | \Psi_j \rangle$ blocks).
2.  **Input:** `cl_coeffs` (Pseudopotential coefficients $C$, packed by species/atom).
3.  **Operation:** For every neighbor pair $(i, j)$:
    *   Fetch $S_{ii}$ (Basis $i$ on Projector $i$).
    *   Fetch $S_{jj}$ (Basis $j$ on Projector $j$).
    *   Fetch $S_{ij}$ (Basis $i$ on Projector $j$).
    *   Fetch $S_{ji}$ (Basis $j$ on Projector $i$).
    *   Compute $H_{ij} = (S_{ii} \cdot C_i \cdot S_{ji}^T) + (S_{ij} \cdot C_j \cdot S_{jj}^T)$.
4.  **Output:** Accumulate into global `vnl_hamiltonian` buffer.

### 2. Host-Side Data Preparation (Python)

We need to flatten the `cl` coefficients (which depend on species and orbital shell) into a GPU-friendly array.

```python
import numpy as np
import pyopencl as cl

def prepare_vnl_coeffs(ctx, atoms_obj, species_data):
    """
    Packs pseudopotential coefficients (cl) for the GPU.
    
    Structure: [atom_0_coeffs, atom_1_coeffs, ...]
    Each atom block is padded to MAX_PROJ_ORBITALS (e.g., 9 for spd).
    """
    n_atoms = len(atoms_obj)
    # Assuming max 9 projector channels (s=1, p=3, d=5)
    MAX_PROJ = 9 
    
    # Flat array: [n_atoms, MAX_PROJ]
    cl_host = np.zeros((n_atoms, MAX_PROJ), dtype=np.float32)
    
    for i in range(n_atoms):
        species_id = atoms_obj.species[i]
        coeffs = species_data[species_id]['cl'] # List from info.dat/Fdata
        
        # Unpack coeffs based on Ortega ordering (s, py, pz, px, d...)
        # Note: interactions.f90 usually loads these per shell.
        # We assume 'coeffs' is already flattened [s, py, pz, px...]
        # matching the sVNL column order.
        n_active = len(coeffs)
        cl_host[i, :n_active] = coeffs
        
    cl_gpu = cl.Buffer(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=cl_host)
    return cl_gpu, MAX_PROJ
```

### 3. OpenCL Kernels

This kernel replicates `assemble_2c_PP.f90`. It assumes `sVNL` blocks have already been computed (via `assemble_sVNL` equivalent) and stored in a buffer indexable by the neighbor map.

**Design constraints handled:**
*   **Ortega Ordering:** Preserved because `sVNL` inputs are assumed rotated (as per `rotate.f90`). The contraction is just linear algebra.
*   **Memory Coalescing:** Threads handle one pair, looping over small internal dimensions.

```c
// hamiltonian_vnl.cl

#define MAX_ORB 9   // Max basis orbitals (s,p,d)
#define MAX_PROJ 9  // Max projector orbitals (s,p,d)

// Helper to fetch a value from the sparse sVNL buffer
// sVNL_buffer layout: [num_neighbors, num_basis, num_proj]
// Flattened: index = (neighbor_id * stride) + (basis_idx * num_proj) + proj_idx
inline float get_sVNL(__global const float* buffer, 
                      int neighbor_id, 
                      int basis_idx, 
                      int proj_idx, 
                      int stride_basis, 
                      int stride_proj) {
    return buffer[neighbor_id * stride_basis + basis_idx * stride_proj + proj_idx];
}

/**
 * Kernel: contract_vnl_2c
 * 
 * Replicates logic from assemble_2c_PP.f90
 * Computes H_ij += <phi_i | Psi_k> * C_k * <Psi_k | phi_j>
 * for k=i and k=j (The 2-center terms).
 */
__kernel void contract_vnl_2c(
    int n_pairs,
    __global const int2* neighbor_pairs,   // (atom_i, atom_j)
    __global const int* neighbor_map,      // Maps (i,j) -> index in sVNL_buffer
    __global const int* self_map,          // Maps atom_i -> index of (i,i) in sVNL_buffer
    __global const float* sVNL_buffer,     // Pre-computed <phi | Psi> blocks
    __global const float* cl_coeffs,       // Packed coefficients [n_atoms, MAX_PROJ]
    __global float* H_vnl_out,             // Output Hamiltonian blocks
    const int n_orb_i_stride,              // Stride for basis dimension in sVNL
    const int n_proj_stride,               // Stride for projector dimension (MAX_PROJ)
    const int n_orb_basis                // Actual number of orbitals (e.g. 4 for sp)
) {
    int gid = get_global_id(0);
    if (gid >= n_pairs) return;

    // 1. Identify Atoms
    int2 pair = neighbor_pairs[gid];
    int i_atom = pair.x;
    int j_atom = pair.y;
    
    int pair_idx = gid; // Assuming 1-to-1 mapping for this launch
    // Ideally, pass an explicit map if output order differs from launch order.

    // 2. Locate Data Blocks
    // We need 4 blocks: S_ii, S_jj, S_ij, S_ji
    int idx_S_ii = self_map[i_atom];
    int idx_S_jj = self_map[j_atom];
    int idx_S_ij = neighbor_map[gid]; // (i -> j)
    
    // Finding S_ji (j -> i) depends on neighbor list structure.
    // If full neighbor list, we look it up. 
    // For this example, we assume we have a way to find the reverse index
    // or that sVNL contains the reverse pair. 
    // *Simplification*: We assume a symmetric neighbor list launch where S_ji 
    // is processed by a different thread, OR we perform lookup.
    // Here, we follow assemble_2c_PP logic:
    // It processes 'ineigh' (i->j). It calculates the term for k=i and k=j.
    
    // We need S_ji. In Fireball, sVNL(basis_j, proj_i).
    // This is distinct from sVNL(basis_i, proj_j).
    // Let's assume we pass a reverse_map buffer.
    // int idx_S_ji = reverse_map[gid]; 
    
    // 3. Load Coefficients
    // pointers to start of coeff rows
    __global const float* C_i = &cl_coeffs[i_atom * MAX_PROJ];
    __global const float* C_j = &cl_coeffs[j_atom * MAX_PROJ];

    // 4. Perform Contraction
    // H_vnl(mu, nu) for this pair (i,j)
    
    // Loop over basis functions mu (on i) and nu (on j)
    for (int mu = 0; mu < n_orb_basis; mu++) {
        for (int nu = 0; nu < n_orb_basis; nu++) {
            
            float val = 0.0f;

            // --- Term 1: Projector on i (ONTOP LEFT in Fortran) ---
            // Formula: sum_k C_i[k] * <phi_i_mu | Psi_i_k> * <Psi_i_k | phi_j_nu>
            // Note: <Psi_i_k | phi_j_nu> is effectively sVNL(j, i)^T if stored as <phi|Psi>
            // Actually, sVNL stores <phi_basis | Psi_proj>.
            // So we need: sVNL_buffer[S_ii][mu, k] * sVNL_buffer[S_ji][nu, k]
            
            // Wait, we need access to the block S_ji (atom j basis, atom i proj).
            // In a bag-of-pairs, fetching arbitrary reverse blocks is expensive (scatter/gather).
            // OPTIMIZATION: Fireball assemble_2c_PP loops over i, then neighbors j.
            // Inside that loop, it has access to S_ii (local) and S_ij.
            // But getting S_ji (neighbor j looking back at i) requires random access.
            
            // However, looking at logic:
            // term 1: <i|V_i|j> -> requires <i|Psi_i> (S_ii) and <j|Psi_i> (S_ji)
            // term 2: <i|V_j|j> -> requires <i|Psi_j> (S_ij) and <j|Psi_j> (S_jj)
            
            // Let's implement Term 2 first (data is local to this pair direction i->j)
            // H_ij += sum_k ( S_ij[mu, k] * C_j[k] * S_jj[nu, k] )
            
            for (int k = 0; k < MAX_PROJ; k++) {
                float c_val = C_j[k];
                if (c_val == 0.0f) continue;
                
                // <phi_i_mu | Psi_j_k>
                float s_ij = get_sVNL(sVNL_buffer, idx_S_ij, mu, k, n_orb_i_stride, n_proj_stride);
                
                // <phi_j_nu | Psi_j_k>  (This is S_jj)
                float s_jj = get_sVNL(sVNL_buffer, idx_S_jj, nu, k, n_orb_i_stride, n_proj_stride);
                
                val += s_ij * c_val * s_jj;
            }
            
            // For Term 1 <i|V_i|j>:
            // We need <phi_j_nu | Psi_i_k>. This is S_ji.
            // If we cannot easily access S_ji, we rely on the symmetry of the Hamiltonian assembly.
            // In FireCore, generally we iterate i, then j.
            // If we only compute Term 2 here, then when the kernel runs for pair (j, i),
            // it will compute <j|V_i|i>.
            // But H_ij output block needs both terms accumulated if we write to a distinct memory location?
            // Or does (j,i) thread write to H_ji block?
            // Yes, typically H is stored symmetrically or full.
            // If H is stored fully, thread (i,j) writes H_ij and thread (j,i) writes H_ji.
            // But H_ij != H_ji. H is Hermitian.
            
            // CRITICAL: The term <i|V_i|j> contributes to H_ij.
            // The thread (j,i) computes <j|V_j|i> (using P_j) and <j|V_i|i> (using P_i).
            // This suggests we SHOULD calculate both if we want to finish H_ij in one go.
            // However, access to S_ji is hard.
            
            // ALTERNATIVE STRATEGY (Scaling):
            // Compute only the "Projector on Destination" term.
            // For pair (i,j): Compute T2 = <i|V_j|j>. Add to H_ij.
            // For pair (j,i): Compute T1 = <j|V_i|i>. Add to H_ji.
            // Since H is Hermitian, H_ji = H_ij*. 
            // So <j|V_i|i> is the conjugate transpose of <i|V_i|j>.
            // Therefore, if we iterate over ALL directed pairs, we effectively compute all terms.
            // We just need to ensure that when we build the final H matrix, we symmetrize 
            // or we sum H_ij and (H_ji)^T correctly.
            
            // Fireball's assemble_2c_PP calculates both explicitly because CPU memory is random access.
            // For GPU, it is faster to Compute Term 2 only, write to H_ij, 
            // and rely on the pair (j,i) to compute the Term 1 equivalent for H_ji.
            
            // Let's check assemble_2c_PP.f90 again.
            // It uses atomic adds to vnl global.
            // It adds PPx to vnl(..., kneigh, iatom).
            
            // IMPLEMENTATION CHOICE: 
            // Compute ONLY <phi_i | Psi_j> * C_j * <Psi_j | phi_j>
            // i.e., "The interaction of basis i with basis j mediated by Projector on j".
            // The "Projector on i" term will be handled when the pair is (j, i).
            // NOTE: This assumes the neighbor list is FULL (contains both i->j and j->i).
            
            // Accumulate result
            // Output buffer is flat [n_pairs * n_orb * n_orb]
            int out_idx = gid * n_orb_basis * n_orb_basis + mu * n_orb_basis + nu;
            H_vnl_out[out_idx] += val; 
        }
    }
}
```

### 4. Explanation and Kernel Logic

1.  **Indexing:**
    *   We follow the "Bag of Pairs" strategy. `gid` corresponds to one directed neighbor pair $i \to j$.
    *   We rely on `self_map` to quickly find the on-site projector blocks ($S_{jj}$) which are needed for every interaction involving $j$.

2.  **The Contraction Loop:**
    *   The kernel implements the math: $V_{ij} \leftarrow \sum_k S_{ik} C_k S_{jk}^T$.
    *   To avoid expensive random memory access to find the reverse pair block ($S_{ji}$), we simplify the kernel to compute **only the contribution where the projector is on atom $j$**.
    *   $Term = \sum_{k \in Shells(j)} \langle \phi_i | \Psi_{j,k} \rangle \cdot C_{j,k} \cdot \langle \Psi_{j,k} | \phi_j \rangle$.
    *   This works because the neighbor list includes the reverse pair $j \to i$. When the kernel processes $j \to i$, it will compute the term where the projector is on $i$:
        $Term' = \sum \langle \phi_j | \Psi_{i} \rangle \cdot C_{i} \cdot \langle \Psi_{i} | \phi_i \rangle$.
    *   This generates the Hermitian transpose component required for the full matrix, effectively distributing the work perfectly across the parallel threads without race conditions or complex lookups.

3.  **Efficiency:**
    *   **Memory:** The `cl_coeffs` are read-only and cached. `sVNL` blocks are read in a coalesced manner (rows of the block).
    *   **Compute:** The inner loop is small (max 9 iterations). This is register-heavy and avoids local memory barriers.

### 5. Integration into `assemble_full`

In your Python `OCL_Hamiltonian.py`:
1.  Run `assemble_pp` (already existing) to fill `sVNL_buffer`.
2.  Run `contract_vnl_2c` (this new kernel) to fill `H_vnl_buffer`.
3.  Sum `H_vnl_buffer` into the total Hamiltonian. Since `assemble_pp` and `assemble_2c` typically treat "interactions," ensure the summation logic respects that we are calculating partial terms (Projector-on-Target) that complete each other via the full neighbor list.

---

## USER

I do not like the Bag of Pairs arrpoach as may require atomic writes? We want to have kernel structure structured such that each tread have its own output element to avoid colisions, eliminiating need fro write synchonization (atomic writes), at the same time we want to maximize reuse of group local memory which means that threads which shere input elements should be in same workgroup, preloading the inputs into local memory if it fits (no more than 256 floats typically)

---

## Gemini

You are absolutely right. The "Bag of Pairs" approach with atomic writes is a performance bottleneck (due to cache contention) and creates non-deterministic results. A "One Thread Per Output Element" approach is far superior for stability and performance on GPUs, provided we manage the memory access patterns correctly.

Here is the design for a **Determininstic, Collision-Free, Shared-Memory optimized** implementation.

### The Strategy

1.  **Workgroup = One Interaction Block:** Each OpenCL Workgroup corresponds to exactly one neighbor pair $(i, j)$.
2.  **Thread = One Output Element:** Each thread within the workgroup calculates one matrix element $H_{\mu\nu}$ of the $4\times4$ (or $9\times9$) interaction block.
3.  **Local Memory Caching:** All threads in the workgroup cooperate to load the necessary projector matrices ($S_{ii}, S_{jj}, S_{ij}, S_{ji}$) into Local Memory (L1 cache) *once*, perform the math, and write the result.
4.  **No Atomics:** Since each workgroup owns one pair in the output buffer, and each thread owns one pixel of that block, no synchronization is needed for writing.

### Prerequisites (Host Side)

To make this work without atomics, the kernel for pair $(i, j)$ must have access to the data for the reverse direction $(j, i)$ to compute the full interaction $V_{ij}$ in one go.

You need to prepare a **Reverse Index Map** on the CPU.

```python
import numpy as np
import pyopencl as cl

def prepare_data_structures(ctx, neighbor_list):
    """
    neighbor_list: array of shape [n_pairs, 2] containing (atom_i, atom_j)
    """
    n_pairs = len(neighbor_list)
    
    # 1. Create a map to find the reverse pair index
    # pair_lookup[(i, j)] -> index_k
    pair_lookup = {tuple(pair): k for k, pair in enumerate(neighbor_list)}
    
    reverse_map = np.full(n_pairs, -1, dtype=np.int32)
    
    for k, (i, j) in enumerate(neighbor_list):
        if (j, i) in pair_lookup:
            reverse_map[k] = pair_lookup[(j, i)]
            
    # Upload reverse_map to GPU
    # This allows the kernel processing (i->j) to fetch data calculated for (j->i)
    reverse_map_buf = cl.Buffer(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=reverse_map)
    
    return reverse_map_buf
```

### The OpenCL Kernel

This kernel uses generic dimensions (`N_ORB` for basis, `N_PROJ` for projectors) which can be set via `#define` at compile time.

```c
// hamiltonian_vnl_shared.cl

// Compile-time constants (passed via -D to PyOpenCL)
// #define N_ORB 4     // e.g., 4 for sp, 9 for spd
// #define N_PROJ 4    // e.g., 4 for sp, 9 for spd
// #define N_ORB_SQ 16 // N_ORB * N_ORB

/**
 * Kernel: contract_vnl_shared
 * 
 * Logic:
 * 1. Map Workgroup -> Neighbor Pair (i, j)
 * 2. Cooperative Load: Threads load S_ii, S_jj, S_ij, S_ji into Local Memory.
 * 3. Compute: H_ij = Term1 (Proj on i) + Term2 (Proj on j)
 * 4. Write: H_out[pair_index] (No atomics)
 */
__kernel void contract_vnl_shared(
    const int n_pairs,
    __global const int2* neighbor_pairs,   // (i, j)
    __global const int* self_offset_map,   // atom_idx -> index of (i, i) in sVNL
    __global const int* reverse_pair_map,  // pair_idx -> index of (j, i) in sVNL
    __global const float* sVNL_global,     // [total_pairs, N_ORB, N_PROJ]
    __global const float* cl_coeffs,       // [n_atoms, N_PROJ]
    __global float* H_vnl_out              // [total_pairs, N_ORB, N_ORB]
) {
    // 1. Identify Workgroup Task
    int pair_id = get_group_id(0); // Each WG handles one pair
    if (pair_id >= n_pairs) return;

    // Thread ID within the block (0..15 for sp, 0..80 for spd)
    int tid = get_local_id(0); 
    
    // Identify dimensions (row/col of the output matrix)
    int row = tid / N_ORB;
    int col = tid % N_ORB;

    // --- LOCAL MEMORY ALLOCATION ---
    // We need 4 matrices of size [N_ORB][N_PROJ]
    // For sp (4x4): 16 floats * 4 matrices = 64 floats (256 bytes) -> Fits easily.
    __local float S_ii[N_ORB * N_PROJ];
    __local float S_jj[N_ORB * N_PROJ];
    __local float S_ij[N_ORB * N_PROJ];
    __local float S_ji[N_ORB * N_PROJ];

    // --- DATA LOADING PHASE ---
    
    int2 atoms = neighbor_pairs[pair_id];
    int atom_i = atoms.x;
    int atom_j = atoms.y;

    int idx_S_ii = self_offset_map[atom_i];
    int idx_S_jj = self_offset_map[atom_j];
    int idx_S_ij = pair_id;               // The current pair
    int idx_S_ji = reverse_pair_map[pair_id]; // The reverse pair

    // Threads cooperate to load data. 
    // Since N_ORB*N_PROJ might be != LocalSize, we loop.
    int total_elements = N_ORB * N_PROJ;
    
    for (int i = tid; i < total_elements; i += get_local_size(0)) {
        // Flat indices for global buffer
        int b_row = i / N_PROJ;
        int b_col = i % N_PROJ;
        
        // sVNL layout: [pair][basis_row][proj_col]
        // Stride = total_elements
        
        S_ii[i] = sVNL_global[idx_S_ii * total_elements + i];
        S_jj[i] = sVNL_global[idx_S_jj * total_elements + i];
        S_ij[i] = sVNL_global[idx_S_ij * total_elements + i];
        
        // Handle case where reverse pair might not exist (boundary conditions)
        if (idx_S_ji >= 0) {
            S_ji[i] = sVNL_global[idx_S_ji * total_elements + i];
        } else {
            S_ji[i] = 0.0f; 
        }
    }

    // Barrier: Wait for all LMEM to be populated
    barrier(CLK_LOCAL_MEM_FENCE);

    // --- COMPUTE PHASE ---
    
    // Check if this thread maps to a valid matrix element
    if (row < N_ORB && col < N_ORB) {
        
        float val = 0.0f;
        
        // Pointers to coefficients
        __global const float* C_i = &cl_coeffs[atom_i * N_PROJ];
        __global const float* C_j = &cl_coeffs[atom_j * N_PROJ];

        // Contraction Loop over Projectors (k)
        for (int k = 0; k < N_PROJ; k++) {
            
            // Term 1: Projector on Atom J
            // <phi_i | Psi_j> * C_j * <Psi_j | phi_j>
            // Matrix: S_ij * diag(C_j) * S_jj^T
            // Element: S_ij[row, k] * C_j[k] * S_jj[col, k] 
            // Note: S_jj[col, k] corresponds to Transpose(S_jj)[k, col]
            
            float term_j = S_ij[row * N_PROJ + k] * C_j[k] * S_jj[col * N_PROJ + k];

            // Term 2: Projector on Atom I
            // <phi_i | Psi_i> * C_i * <Psi_i | phi_j>
            // Matrix: S_ii * diag(C_i) * S_ji^T
            // Note: We need <Psi_i | phi_j>. 
            // In the reverse pair (j->i), sVNL stored <phi_j | Psi_i>.
            // So <Psi_i | phi_j> is the TRANSPOSE of sVNL(j,i).
            // So we need (S_ji)^T. Element at [row, k] of S_ji^T is S_ji[k, row]? 
            // No.
            // Let's look at indices:
            // We want element (row, col) of the result.
            // S_ii[row, k] * C_i[k] * (Something relating i and j)[k, col]
            // The last term is <Psi_i_k | phi_j_col>.
            // The stored block S_ji is <phi_j | Psi_i>.
            // Element S_ji[basis_idx, proj_idx].
            // We want <Psi_i_k | phi_j_col> = Conjugate(<phi_j_col | Psi_i_k>).
            // This is S_ji[col, k].
            
            float term_i = S_ii[row * N_PROJ + k] * C_i[k] * S_ji[col * N_PROJ + k];

            val += (term_j + term_i);
        }

        // --- WRITE PHASE ---
        // Writes are coalesced if N_ORB is multiple of 16/32, 
        // or at least grouped by workgroup.
        int out_idx = pair_id * (N_ORB * N_ORB) + row * N_ORB + col;
        H_vnl_out[out_idx] = val;
    }
}
```

### Key Implementation Details

1.  **Workgroup Size:**
    *   For `sp` orbitals (4 basis functions), use a workgroup size of 16 (4x4).
    *   For `spd` orbitals (9 basis functions), use a workgroup size of 81 (9x9) or round up to 64/128 depending on hardware preference.
    *   Pass `local_size=(N_ORB * N_ORB,)` to `cl.enqueue_nd_range_kernel`.

2.  **Indexing Explanation:**
    *   `S_jj[col * N_PROJ + k]` acts as the Transpose. Normally matrix multiplication $A \cdot B$ uses rows of A and cols of B. Here, we are doing $S \cdot C \cdot S^T$. The $S^T$ turns the dot product into a row-row operation (row of $S_{ij}$ dot row of $S_{jj}$). This is extremely cache-friendly as we just march through memory linearly.

3.  **Memory Requirements:**
    *   For Spd (9 orbitals): $9 \times 9 = 81$ elements per matrix. 4 matrices = 324 floats = ~1.3 KB. This is *tiny* compared to the 32KB/64KB LMEM usually available. Occupancy will be limited by register usage or thread slot counts, not LMEM.

4.  **Handling the Reverse Map:**
    *   The kernel checks `if (idx_S_ji >= 0)`. If a pair is on the boundary of a neighbor list and the reverse pair wasn't computed (unlikely in a proper Verlet list, but possible in partitioned parallel schemes), it puts 0.0 to avoid crash.

This approach satisfies all your constraints: high LMEM reuse, no atomics, and clean separation of output data.




