# Neighbor Exclusions for Non-Covalent OpenCL Kernels

## Context and Goals
Computing non-covalent interactions on the GPU requires excluding 1-2 (bonded) and 1-3 (angle) pairs while preserving high throughput. Each atom can have at most 16 exclusions. Rather than performing up to 16 comparisons for every candidate pair `(i,j)`, we exploit the monotonic traversal of `j` and maintain an exclusion pointer that advances at most once per iteration, delivering amortized *O(1)* exclusion checks per pair.@doc/Topics/ForceFields/NeighborExclusion_ChatGPT.md#608-715

see:
* https://chatgpt.com/share/69021f05-8290-8003-aa5e-3fffa8afa154
* https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221Z3oDXimOUCh8gy4cmErDN_vcrIiFpMki%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing

## Assumptions and Preconditions
- Exclusion lists are generated on the host, limited to ≤16 entries per atom, and sorted strictly by atom index (lower 24 bits).@doc/Topics/ForceFields/NeighborExclusion_ChatGPT.md#639-715
- No atom appears multiple times in an exclusion list, so each `(atom, cell)` pair is unique. This enables single-step pointer advancement without missing matches.@doc/Topics/ForceFields/NeighborExclusion_ChatGPT.md#592-610
- Atom indices fit into 24 bits (≤16,777,215) and cell identifiers into 8 bits (≤255), allowing them to be packed into a single `uint32_t` value `cell<<24 | atom`.@doc/Topics/ForceFields/NeighborExclusion_Gemini.md#681-909
- Host preprocessing verifies the above guarantees (duplicate detection, bounds checks) before uploading buffers to the GPU.

## Evaluated Alternatives
1. **Hash/Bloom prefilters**: A bitmask pre-check can suppress many scans, but still needs a follow-up linear scan on hash hits, introducing extra state and host-side mask maintenance.@doc/Topics/ForceFields/NeighborExclusion_Gemini.md#214-547
2. **Binary search per pair**: Although `O(log N)` in theory, branching causes warp divergence and irregular memory access, making it slower than linear scans on short lists.@doc/Topics/ForceFields/NeighborExclusion_Gemini.md#372-442
3. **Pointer-based merge scan (final choice)**: Sorting exclusions and keeping a per-thread pointer yields deterministic, branch-light control flow with just one comparison per `(i,j)` step while respecting GPU execution constraints.@doc/Topics/ForceFields/NeighborExclusion_ChatGPT.md#637-715

## Final Data Layout and Algorithm
1. **Host packing**: Each exclusion entry becomes a packed `uint32_t` with cell index in the top 8 bits and atom index in the lower 24 bits. Entries are sorted by atom index only, preserving stable grouping of identical atoms (if any), and stored in a flat array with per-atom offsets.@doc/Topics/ForceFields/NeighborExclusion_Gemini.md#723-909
2. **Kernel traversal**:
   - Maintain pointer `iex` into the exclusion slice for atom `i`.
   - Because `j` increments monotonically and exclusions are unique, a single conditional increment advances `iex` when the pointed atom is below `j`.
   - When `UNPACK_ATOM_IDX(excl[iex]) == j`, compare the packed cell against the current PBC image index, skipping only that specific image. Otherwise the interaction is computed.@doc/Topics/ForceFields/NeighborExclusion_Gemini.md#805-900

This architecture ensures at most one exclusion comparison per candidate pair, lifting the dominant cost compared with per-pair scans, while keeping the kernel branch structure simple and GPU-friendly.@doc/Topics/ForceFields/NeighborExclusion_ChatGPT.md#691-705

## Host-Side Preparation (C++)
```cpp
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#define PACK_EXCLUSION(atom_idx, cell_idx) \
    (static_cast<uint32_t>(cell_idx) << 24 | static_cast<uint32_t>(atom_idx))
#define UNPACK_ATOM_IDX(packed_val) ((packed_val) & 0x00FFFFFFu)

struct ExclusionPair { uint32_t atom; uint32_t cell; };

void prepare_packed_exclusions(
    const std::vector<std::vector<ExclusionPair>>& per_atom_raw,
    std::vector<uint32_t>& packed_flat,
    std::vector<uint32_t>& offsets
){
    const size_t nAtoms = per_atom_raw.size();
    packed_flat.clear();
    offsets.assign(nAtoms + 1u, 0u);

    uint32_t cursor = 0u;
    std::vector<uint32_t> tmp;
    tmp.reserve(16);

    for(size_t i=0; i<nAtoms; ++i){
        offsets[i] = cursor;
        tmp.clear();
        tmp.reserve(per_atom_raw[i].size());

        for(const auto& pr : per_atom_raw[i]){
            if(pr.atom >= (1u<<24) || pr.cell >= (1u<<8))
                throw std::runtime_error("exclusion overflow");
            tmp.push_back(PACK_EXCLUSION(pr.atom, pr.cell));
        }

        std::sort(tmp.begin(), tmp.end(), [](uint32_t a, uint32_t b){
            return UNPACK_ATOM_IDX(a) < UNPACK_ATOM_IDX(b);
        });

        // optional duplicate check to guarantee uniqueness
        for(size_t k=1; k<tmp.size(); ++k){
            if(UNPACK_ATOM_IDX(tmp[k-1]) == UNPACK_ATOM_IDX(tmp[k]))
                throw std::runtime_error("duplicate exclusion atom");
        }

        packed_flat.insert(packed_flat.end(), tmp.begin(), tmp.end());
        cursor += tmp.size();
    }

    offsets[nAtoms] = cursor;
}
```
This host utility enforces packing constraints, sorts by atom index, and guards against duplicates, ensuring the kernel’s single-step pointer advance remains correct.@doc/Topics/ForceFields/NeighborExclusion_Gemini.md#723-909

## Kernel Fragment (OpenCL C)
```c
#define UNPACK_ATOM_IDX(v) ((v) & 0x00FFFFFFu)
#define UNPACK_CELL_IDX(v) ((v) >> 24)

__kernel void getNonBond(
    const int4 nDOFs,
    __global const float4* apos,
    __global float4* aforce,
    __global const float4* REQs,
    __global const uint* packed_excl,
    __global const uint* excl_offsets,
    __global const cl_Mat3* lvecs,
    const int4 nPBC,
    const float4 GFFParams
){
    const int iG  = get_global_id(0);
    const int nat = nDOFs.x;
    if(iG >= nat) return;

    const int iS  = get_global_id(1);
    const int nnode = nDOFs.y;
    const int nvec  = nat + nnode;
    const int i0a = iS*nat;
    const int i0v = iS*nvec;
    const int iaa = i0a + iG;
    const int iav = i0v + iG;

    const uint excl_begin = excl_offsets[iaa];
    const uint excl_end   = excl_offsets[iaa+1];
    uint iex = excl_begin;

    // ... load per-atom data, local memory tiles, etc. ...

    for(int j=0; j<nat; ++j){
        if(j==iG) continue;

        if(iex < excl_end && UNPACK_ATOM_IDX(packed_excl[iex]) < (uint)j)
            ++iex;

        bool skip_atom = (iex < excl_end && UNPACK_ATOM_IDX(packed_excl[iex]) == (uint)j);

        if(skip_atom && (nPBC.x+nPBC.y+nPBC.z)==0){
            continue; // excluded in non-PBC case
        }

        // compute displacement, mix parameters ...

        if((nPBC.x+nPBC.y+nPBC.z)>0){
            int ipbc = 0;
            // loop over 27 images (omitted for brevity)
            // ...
            bool skip_img = skip_atom && (UNPACK_CELL_IDX(packed_excl[iex]) == (uint)ipbc);
            if(skip_img) continue;
            // accumulate force contribution
        }else{
            // accumulate force contribution for non-excluded pair
        }
    }

    // write back accumulated force
}
```
The kernel keeps `iex` monotonic, so each exclusion entry is evaluated at most once per `(i,j)` traversal. Packed indices enable comparing atom and cell in a single register load, ensuring the PBC case only skips the matching image while computing all others.@doc/Topics/ForceFields/NeighborExclusion_Gemini.md#805-900

## Key Takeaways
- Sorting per-atom exclusions and guaranteeing uniqueness enable single-step pointer advancement, cutting exclusion checks to one comparison per pair.
- Packing `(cell, atom)` into `uint32_t` reduces memory bandwidth and keeps comparisons simple.
- Host validation is essential: it enforces the guarantees the kernel relies on, preventing divergence and missed exclusions.
- This pattern generalizes to any small, monotonic exclusion sets and is well-suited to OpenCL SIMT execution.
