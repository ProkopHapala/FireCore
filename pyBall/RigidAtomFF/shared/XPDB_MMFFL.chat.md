## USER

look we have this position based dynamics solver of molecular dynamics XPDB in pyOpenCL, @XPDB_new.cl @XPDB_new.py @test_TiledJacobi_molecules.py 

we need to couple it to proper  type assignment using this @MMFFL.py @MMFF.py @MMparams.py @MolecularDynamics.py @AtomicSystem.py 

@AtomTypes.dat @ElementTypes.dat 

identify which exact function we need to use and where, the goal is for each atom in molecule emit and print proper type name, number of bonds , pi-orbitals and number of electron paits (note pi-orbitals are proxy for bond order sp3 has npi=0, sp2 npi=1, sp1 npi=2) electron pairs are created and counted by this function add_electron_pairs in @AtomicSystem.py 

generally for main period atoms like C,N,O, we use octed rule so nbond+npi+nepair = 4 (two electron per each sigma bond, pi-orbital or electron pair)

look carefully on @MMFFL.py this is exactly implementing the angle-derived bonds which we need in @test_TiledJacobi_molecules.py , pehrpas we can eliminate/replace build_bond_arrays_with_angles but I do not remember if @MMFFL.py is properly tested, we need to test it

now make detailed notes and a plan

---

## GPT 5.2

Got it. I’ll first map the existing type-assignment + bonding logic in `pyBall/OCL/MMFFL.py`, `MMFF.py`, `MMparams.py`, and `AtomicSystem.py`, then map how XPDB builds/uses bonds in `test_TiledJacobi_molecules.py`. After that I’ll identify the *exact* functions you should call (and where) to: (1) derive bonds (including angle-derived), (2) assign atom types, and (3) compute/print `nbond`, `npi`, and `nepair` (via `add_electron_pairs`). I’ll also propose a small test harness to validate `MMFFL`’s bond/angle builder against the current XPDB builder.

## Detailed notes (what exists today, and what to reuse)

### 1) XPDB solver inputs (what [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0) currently builds)
- XPDB ([XPDB_new.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB_new.cl:0:0-0:0)) only knows:
  - `bond_indices_local[ i*N_MAX_BONDED + slot ]`
  - `bond_lengths[...]`
  - `bond_stiffness[...]`
- In [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0), bonds are built by:
  - [AtomicSystem.AtomicSystem(fname=...)](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78)
  - [mol.findBonds()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:135:4-139:29) if missing
  - [mol.neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23)
  - then `build_bond_arrays_with_angles(...)` which:
    - takes **first neighbors** from `mol.bonds`
    - adds **second neighbors** by a *single global* relaxed angle `alpha0_deg` and global stiffness `k_angle`

So today: angle-derived “bonds” are **purely geometric** and **type-blind**.

### 2) Where “proper type assignment” lives (and the exact function)
There is no robust type assignment in [AtomicSystem.py](cci:7://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:0:0-0:0) itself (it mostly stores element symbols like `"C"`, `"H"` in `enames` from XYZ/MOL loaders).

The **actual UFF/MMFF-like type assignment** (e.g. `C_3`, `C_2`, `C_R`, `O_2`, `N_3`, …) is implemented in:

- **[pyBall/OCL/UFFbuilder.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/UFFbuilder.py:0:0-0:0)**
  - **Class:** `UFF_Builder`
  - **Key method (the one you want):** `UFF_Builder.assign_uff_types()`
    - It computes per-atom UFF types (and also bond orders internally) based on bonding topology.
    - It uses the canonical parameter tables:
      - [MMparams.MMFFparams(data_path)](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMparams.py:244:0-261:51) which loads:
        - [cpp/common_resources/ElementTypes.dat](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/ElementTypes.dat:0:0-0:0)
        - [cpp/common_resources/AtomTypes.dat](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat:0:0-0:0)

Important detail: `UFF_Builder` stores assigned types in its internal `builder_atoms[i].type` (an integer index into `params.atypes`). It does **not** automatically rewrite `mol.enames` for you.

So the **exact function you need to call** for “proper type name per atom” is:

- `UFF_Builder(mol, ...).assign_uff_types()`

### 3) Where angle-derived “linear bonds” are implemented (MMFFL)
In [pyBall/OCL/MMFFL.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:0:0-0:0):

- **[MMFFL.build_linearized(mol, ...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:40:4-71:19)** is the high-level entry point.
- It calls:
  - `super().toMMFFsp3_loc(...)` to build MMFF-like neighbor arrays and per-bond equilibrium parameters.
  - then [_build_angle_bonds(mol, atom_types)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:80:4-95:78) which does:
    - [angs, triplets = mol.findAngles(ngs=mol.ngs)](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:178:4-181:68)
    - for each angle triplet `(ja, ib, jc)`:
      - uses the **central atom type’s** `atype.Ass` (equilibrium angle)
      - uses equilibrium bond lengths `r_ab`, `r_bc` (from `self.bLs` if possible)
      - computes the implied “angle bond” length `l0` via law-of-cosines
      - stiffness `k` from either `self.K_ang` override or `atype.Kss`
      - stores it as `self.linear_bonds.append((ja, jc, l0, k, ("angle", (...))))`

This is **exactly** the “angle-derived bonds” you said XPDB needs (but now type-aware).

Crucial caveat: [MMFFL](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:12:0-210:76) expects `mol.enames[i]` to be **atom type names that exist in `atom_types`**, not plain element symbols. Otherwise [_build_angle_bonds](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:80:4-95:78) will do:
```py
atype = atom_types.get(mol.enames[ib])
```
and fail to find e.g. `"C"` vs `"C_3"`.

So we need a small bridge step: **after UFF typing, temporarily set `mol.enames` to the assigned UFF type names** (e.g. `"C_3"`, `"C_R"`, ...), at least for the MMFFL build.

### 4) Electron pairs ([add_electron_pairs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52)) and what it really does
[AtomicSystem.add_electron_pairs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52):
- Only acts for elements in `VALENCE_DICT` = `{ "O":(2,2), "N":(3,1) }`
- Uses current neighbor list `self.ngs` (so bonds must exist and [neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23) must be run).
- Computes:
  - `nsigma = len(self.ngs[i])`
  - `npi = nb - nsigma` (with `nb` from `VALENCE_DICT`)
- Then calls [make_epair_geom(i, npi, nsigma)](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:903:4-947:52) to **actually append “E” atoms** and bonds.

So if your goal is “print number of electron pairs per atom”:
- The most reliable approach is **not** to mutate the molecule just to count.
- Since you already have [AtomTypes.dat](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat:0:0-0:0) where each type has `nepair` and `valence`, you can compute counts directly from assigned type:
  - `nepair = atype.nepair`
  - `npi = atype.valence - nbond`  (matches the algebra used in [MMFF.toMMFFsp3_loc](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFF.py:244:4-539:31), and matches your “octet-like” rule)

Then you can still optionally call [add_electron_pairs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52) later if/when you actually want explicit E-pair particles.

## What to print per atom (and how to compute it robustly)

Once UFF typing is done, for each atom `i`:
- **Type name:** `tname = uff.params.atypes[uff.builder_atoms[i].type].name`
- **Number of sigma bonds (`nbond`):** `nbond = len(mol.ngs[i])` (after [mol.neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23))
- **Electron pairs (`nepair`):** `nepair = uff.params.atypes[type_idx].nepair`
- **Pi orbitals proxy (`npi`):**
  - `npi = uff.params.atypes[type_idx].valence - nbond`
  - clamp if needed to `0..2` for main group (if you want to be strict)
- **Octet check for main period C/N/O (your invariant):**
  - assert: `nbond + npi + nepair == 4` for `tname[0] in {"C","N","O"}` (or element from `params.atypes[type_idx].element_name`)

This uses the authoritative data in [AtomTypes.dat](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat:0:0-0:0) and the actual neighbor topology.

## Integration point (where to put the calls in XPDB test)
The natural place is inside `load_molecule_bonds(...)` in [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0), right after:
- [mol.findBonds()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:135:4-139:29)
- [mol.neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23)

Add a “typing + printout” step there:
1) Run UFF typing (`UFF_Builder`)
2) Print per atom: `type name`, `nbond`, `npi`, `nepair`, and assert the octet rule for C/N/O.
3) Build XPDB bond arrays:
   - Replace (or gate) `build_bond_arrays_with_angles(...)`
   - Instead build “angle bonds” using [MMFFL](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:12:0-210:76):
     - Temporarily set `mol.enames = assigned_uff_type_names`
     - Call [MMFFL(...).build_linearized(mol, include_linear=True, bUFF=True/False depending on what you want for bond lengths)](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:12:0-210:76)
     - Extract only angle bonds from `mmffl.linear_bonds` with tag `"angle"`
     - Merge these with the primary bonds into XPDB fixed-size arrays.

I would keep the existing `build_bond_arrays_with_angles` behind a flag for comparison testing rather than deleting it (per your “never delete” rule).

## Risks / things to test immediately
- **[MMFFL](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:12:0-210:76) depends on type names in `mol.enames`.** Without rewriting `enames`, [_build_angle_bonds()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:80:4-95:78) will silently skip everything (`atype is None`).
- **[MMFF.toMMFFsp3_loc](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFF.py:244:4-539:31) currently assumes `mol.ngs` provides both neighbor atoms and bond indices (`ngbs`)**, but in your [AtomicSystem.neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23) you get a dict of neighbor->bond_index (via `au.neigh_bonds`). This is probably fine, but it’s a high-risk area; if it breaks, [MMFFL](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:12:0-210:76) won’t get correct equilibrium bond lengths (`bLs`) and will fall back to geometric distances.
- **[AtomicSystem.add_electron_pairs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52) mutates topology** (adds atoms/bonds). Don’t call it in the XPDB build path unless you really want to simulate explicit electron pairs in XPDB.

# Plan (concrete milestones)

## 1) Add a type-assignment + reporting step (UFF-based) in `load_molecule_bonds`
- **Goal**
  - For each atom print:
    - `i`, element symbol (original), assigned `UFF type name`, `nbond`, `npi`, `nepair`
  - Assert for C/N/O:
    - `nbond + npi + nepair == 4`
- **Exact function to use**
  - `pyBall/OCL/UFFbuilder.py : UFF_Builder.assign_uff_types()`

## 2) Bridge UFF types into MMFFL and generate type-aware angle-derived bonds
- **Goal**
  - Use [MMFFL._build_angle_bonds()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:80:4-95:78) via [MMFFL.build_linearized()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:40:4-71:19) to generate angle-derived neighbor bonds (second neighbors) using:
    - per-central-atom `Ass`
    - per-bond equilibrium lengths
- **Integration**
  - In `load_molecule_bonds`, after typing:
    - create `mol.enames_typed = [...]`
    - swap `mol.enames = mol.enames_typed` temporarily for MMFFL call

## 3) Build XPDB fixed-size bond arrays from:
- primary bonds: `mol.bonds`
- plus angle bonds: `MMFFL.linear_bonds` filtered by tag `"angle"`
- **Goal**
  - Replace/augment current `build_bond_arrays_with_angles` with a new builder (or a selectable path) that uses MMFFL angle bonds.

## 4) Test that MMFFL path is sane (before relying on it)
- **Minimal test set**
  - A small sp3 molecule (e.g. methane/ethane)
  - A planar sp2 system (benzene/pentacene)
  - A carbonyl-containing molecule (formaldehyde / acetone fragment)
- **Checks**
  - Angle-bond count per atom matches expectation:
    - For atom with degree `d`: number of angle pairs = `d*(d-1)/2`
    - Each pair creates a j–k “angle bond” (so per central atom expect that many)
  - Compare the *connectivity* of second neighbors:
    - existing `build_bond_arrays_with_angles` vs MMFFL-derived
  - Compare rest lengths `Ljk`:
    - global-alpha version vs type-specific `Ass` version (should differ for sp3 vs sp2)

## Status
- **Completed**
  - Identified the exact typing function: `UFF_Builder.assign_uff_types()`
  - Identified the exact angle-derived bond function: [MMFFL.build_linearized()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:40:4-71:19) → [_build_angle_bonds()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:80:4-95:78)
  - Located [add_electron_pairs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52) and clarified it mutates the system (best used optionally)
- **Next**
  - You apply the integration changes in [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0) (or switch to **Code mode** and I can propose an exact patch).

If you want, tell me whether you prefer:
- **Option A (recommended):** UFF typing only for printing + MMFFL only for angle-bond generation
- **Option B:** fully convert the molecule to typed `enames` and keep it that way throughout the XPDB test

---

## USER

OK, now please implement it in one step. Our goal should be to maximally reuse existing functionality in @MMFFL.py , eventually finish it if it is unfinished. We want to replace ad-hoc functions in @test_TiledJacobi_molecules.py , this should be only lighwitht scrip using @MMFFL.py . I'm not sure if we want to use UFF and its types exactly but I'm not purely agains it. But we definitely should have poption to use types as defined here
@AtomTypes.dat 

Note that we should add electron pairs as optinal dummy atoms (they may and may not be added). Similarly for pi-orbitals as additinal explicit dummy atoms (there should be flag to swithcit on), and we should have choice if these pi-orbital dummy atoms should be one sided or both sided (above and below atomatic ring). The purpose of these pi-orbitals dummies is to enhace rigidity of aromatic systems because with only angle-derived bonds between real atoms (next nearest neighbors) we will never reproduce out-off plane stiffness which is solved by dihedrals (torsions) and inversions in normal forcefield like UFF. In place of these we want to make pi-derived bonds between pi-dummy and neighbors of host atom, and epair-derived bonds between epair and neighbor of host atom. (this is similar to angle-derived bonds which are between two neibhor atoms of a central atom)

Be carefull that these angle-derived, pi-derived, and epair-derived bonds must be put into different second_neighbor buffer not into original neighbor list, otherwise we will suffer from order-dependnce inifinte recusion and other problems. only just before we export the oribitals to GPU buffers we add both - 1st-order neighbors, and then 2nd order neibhors.

Note that in the final forcefield evaluation, the pi-dummies are involved only in boding interactions (to increase stiffness) but not in collision interactions. So we need to have some counter (natom, ndumise) to distinguish this, as both atoms and dummies are involved in the groups of atoms and reindexing kernel. ( see build_local_topology in @XPDB_new.cl )


Now our phases of tst driven development is
1) we make script just to test type-assignemnt and topology kreation (angle/pi/epair derived bonds), we produce visalization using mapltlib with explicitly showing the bonds as sticks, and eventually aslo exporting it as mol2 file.
2) then we run the kernel and check if the atoms are properly re-indexed to local, we may write some CPU test to map them back from local to global.

now analyse the input files once agaon, and make updated, more specific plan

---

## GPT5.1-Codex

**Updated plan (specific to MMFFL-driven topology builder + XPDB script):**

1. **Type assignment options & reporting (Phase 1 foundations)**
   1.1 Wrap `UFF_Builder.assign_uff_types()` so we can populate per-atom type metadata without mutating the core [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) (store `{index, type_idx, type_name}` plus `valence/nepair` from `MMparams`). @pyBall/OCL/UFFbuilder.py#1-186  
   1.2 Add an alternate path that loads types directly from [AtomTypes.dat](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat:0:0-0:0) (no UFF inference) so users can force specific labels when desired. @cpp/common_resources/AtomTypes.dat#1-71  
   1.3 After either path, compute `(nbond, npi, nepair)` using neighbor counts from [AtomicSystem.neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23) and the type’s valence data, then print + assert the octet rule for C/N/O before we touch XPDB arrays. @pyBall/AtomicSystem.py#135-155, @pyBall/OCL/MMFF.py#300-315  

2. **Topology builder (new helper module)**
   2.1 Use [MMFFL.build_linearized()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:40:4-71:19) to generate angle-derived bonds; ensure `mol.enames` is temporarily remapped to the assigned type names so [_build_angle_bonds](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:80:4-95:78) finds `atype.Ass`/`Kss`. @pyBall/OCL/MMFFL.py#40-97  
   2.2 Extend the helper to toggle pi-dummy creation (`two_pi_dummies`, `lone_pairs_pi`, etc.) and lone-pair dummy creation. Reuse existing MMFFL logic for inserting dummies and recording `linear_bonds`, but keep separate lists for:
        - primary bonds (`mol.bonds`)
        - angle-derived bonds
        - pi-dummy bonds (host↔dummy and dummy↔neighbors)
        - epair-dummy bonds (host↔dummy and dummy↔neighbors)  
      @pyBall/OCL/MMFFL.py#96-139, @pyBall/OCL/MMFF.py#493-513
   2.3 Maintain bookkeeping for counts: `n_real_atoms`, `n_pi_dummies`, `n_epair_dummies`. Tag each dummy so downstream XPDB code can exclude them from collision logic.
   2.4 Provide an aggregation API that returns **two buffers**: first-order neighbors (real atoms only) and second-order/dummy neighbors (angle/pi/epair). Only right before GPU upload do we concatenate them into the fixed `bond_indices_*` arrays.

3. **XPDB test script refactor (Phase 1 driver)**
   3.1 Replace `build_bond_arrays_with_angles` in [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0) with calls into the new topology builder. Keep a flag to fall back to the legacy builder for A/B comparison. @pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py#31-162  
   3.2 Add CLI flags to control: type-source (`--type-source=uff|table`), enable pi dummies, one-sided vs two-sided pi dummies, enable epair dummies, include second-neighbor bonds, etc.
   3.3 After building the topology, emit:
        - per-atom report (type, nbond, npi, nepair)
        - matplotlib visualization showing primary bonds vs angle/pi/epair bonds with distinct colors
        - optional MOL2 export for inspection
   3.4 When preparing buffers for XPDB, append real atoms and dummies to the position/mass arrays, but keep a bitmask so the OpenCL solver can skip dummy atoms in the collision loop while still processing their bond constraints. (XPDB kernels already separate bond vs collision handling, so we just need a `num_real_atoms` threshold when looping in `solve_cluster_jacobi` / `build_local_topology`.) @pyBall/XPDB_AVBD/XPDB_new.cl#250-350, @pyBall/XPDB_AVBD/XPDB_new.py#44-156

4. **Testing roadmap**
   4.1 **Phase 1 (CPU/visual):** run the refactored script on representative molecules (sp³, sp², aromatic, heteroatoms) to verify:
        - type assignment correctness (octet rule)
        - counts of angle/pi/epair bonds per central atom
        - visualization/MOL2 output matches expectations
   4.2 **Phase 2 (GPU reindex):** once topology looks correct, run the existing XPDB kernels, then dump the local-index buffers to confirm:
        - real atoms and dummies maintain consistent global↔local mapping
        - second-order bonds survive the `build_local_topology` reindexing
        - collision kernel ignores dummy indices (based on the mask or count).

---

## USER

add
1) Add "an alternate path that loads types directly from AtomTypes.dat" 
 - we do not need to avoid @UFF.py if it is useful, I only siad I'm not sure we need it. But definitely use @MMparams.py for loading any paremeters, and that should be our only source of parameters, if we need any other params we should implement it there.
2) any topology builder related to angle-derived bonds, epari-derived or pi-derive should be implemented in @MMFFL.py , that is the class dedicated to this (replacement of angular forces by distance constrain, that is why we call it MMFFL - linearized MMFF)

clarify and specify in more detail that

---

## GPT5.1-Codex

Here’s how we’ll address your two requirements:

1. **Alternate type source (AtomTypes.dat via `MMparams`)**
   - [MMparams.MMFFparams](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMparams.py:244:0-261:51) already exposes [read_AtomAndElementTypes](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMparams.py:202:0-208:36) which loads [ElementTypes.dat](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/ElementTypes.dat:0:0-0:0) and [AtomTypes.dat](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat:0:0-0:0) into `self.atom_types_map` / `self.element_types_map`. We can add a lightweight helper (e.g., `load_atom_types_from_table(path)`), but under the hood it just calls [MMparams.MMFFparams(data_path)](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMparams.py:244:0-261:51) and keeps the `atom_types_map` as-is. That gives us direct access to every entry from @/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat without running the UFF topology inference in @/home/prokop/git/FireCore/pyBall/OCL/UFF.py. In practice we’ll offer a switch in the XPDB loader:  
     - `type_source="uff"` → run `UFF_Builder.assign_uff_types()` (same as before).  
     - `type_source="table"` → fetch the [AtomType](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMparams.py:40:0-82:88) records straight from [MMparams.MMFFparams](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMparams.py:244:0-261:51) and assign them to atoms by name (falling back to parent types as needed).  
     Both paths still rely exclusively on `MMparams` for parameter values (Req, Kss, Ass, valence, nepair, etc.), per your requirement that this module is the single source of truth.

2. **All topology builders live inside [MMFFL](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:12:0-210:76)**
   - We will expand @/home/prokop/git/FireCore/pyBall/OCL/MMFFL.py so that [_build_angle_bonds](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:80:4-95:78), [_build_pi_dummies](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:97:4-137:61), and the planned electron-pair helpers cover every “linearized” constraint:
     * **Angle-derived bonds**: already handled by [_build_angle_bonds](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:80:4-95:78) using Ass/Kss from the active atom type. We’ll route all XPDB angle bonds through this path and remove the ad-hoc builder in [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0).
     * **Pi-derived bonds**: [_build_pi_dummies](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:97:4-137:61) already creates dummy nodes and records host↔dummy plus dummy↔neighbor bonds. We’ll gate one-side vs two-side creation via `two_pi_dummies`, and expose a flag to enable/disable the entire pi block when we create the topology descriptor.
     * **Electron-pair (epair) dummies**: we’ll add a sibling helper (e.g., `_build_epair_dummies`) to MMFFL that mirrors what [_build_pi_dummies](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:97:4-137:61) does, but driven by each atom’s `nepair` count (again from `MMparams`). These epair dummies will generate host↔epair and epair↔neighbor linear bonds, just like your description.
   - MMFFL will maintain separate lists for primary bonds (from `mol.bonds`) vs. each derived class (`angle`, `pi`, `epair`). The XPDB loader will simply ask MMFFL for the “first-order” set and for the “second-order/dummy” set, so we only merge them right before uploading to GPU buffers—preventing any recursion or re-use inside MMFFL itself.

With those clarifications, we can keep all parameter IO centralized in `MMparams`, keep all linearization logic in [MMFFL](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:12:0-210:76), and let the XPDB script act purely as a lightweight client that selects the type source, toggles pi/epair options, and formats the resulting buffers/visualizations.

---

## USER

Now implement it, focus on steps 1,2,3) resulging in the 3D view of plot .mol2 file for selected molecule (notice .mol2 export in save_mol2 @AtomicSystem.py @atomicUtils.py .
when you finish coding run it for test and check the ouputed atoms and bonds (including the dumise) if ti make sense. They sime simple molecule for starter like 
@H2O.xyz @OCH2.xyz 

only then 
if you are able you should also try to implement 4, good

---

## GPT 5.2

## What’s implemented (Steps 1–3)

### 1) **Alternate type assignment path using [AtomTypes.dat](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat:0:0-0:0) (via `MMparams` only)**
Implemented in [pyBall/OCL/MMFFL.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:0:0-0:0):

- [MMFFL.assign_type_names(..., type_source="table")](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:71:4-93:18)
  - Maps `mol.enames` like `"O","C","H"` to a valid key in [AtomTypes.dat](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat:0:0-0:0) using:
    - direct match (`"O"`)
    - fallback `"H_"` if needed
    - fallback by matching `element_name` in the table
- [MMFFL.assign_type_names(..., type_source="uff")](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:71:4-93:18)
  - Optional: runs UFF inference via `UFF_Builder.assign_uff_types()`, but still uses `MMparams` as the authoritative parameter DB for the types.

Also added:
- [MMFFL.report_types_and_valence()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:103:4-121:125) prints per atom:
  - `type`, `nbond`, `npi`, `nepair`
  - asserts octet rule for C/N/O (`nbond + npi + nepair == 4`) and fails loudly if violated.

### 2) **All topology building moved into [MMFFL.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:0:0-0:0)**
Implemented in [pyBall/OCL/MMFFL.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:0:0-0:0):

- [MMFFL.build_topology(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:123:4-204:9) returns a topology dict containing:
  - `apos` (including dummy atoms)
  - `type_names` (including `E` and `Pi` where relevant)
  - `is_dummy` mask
  - bond lists separated by class:
    - `bonds_primary`
    - `bonds_angle`
    - `bonds_pi`
    - `bonds_epair`

Derived-bond generation:
- **Angle bonds**: already in [_build_angle_bonds()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:249:4-264:78)
- **Pi dummy bonds**: [_build_pi_dummies()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:266:4-306:61) (now correctly gated by `add_pi`)
- **Epair dummy bonds**: newly added [_build_epair_dummies()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:307:4-346:64) + [_epair_dirs()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:349:4-429:21)  
  (mirrors [AtomicSystem.make_epair_geom](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:903:4-947:52) *in spirit*, but produces dummy particles + “orthogonal stiffening” bonds to neighbors)

Important fix:
- [MMFFL.build_linearized()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:206:4-237:19) now sets `mol.natoms` when missing (AtomicSystem didn’t have it), preventing the earlier crash.

### 3) **XPDB test script becomes a lightweight driver + 3D plot + MOL2 export**
Implemented in [pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0):

- Adds robust repo-root `sys.path` setup (fixes `import pyBall`)
- Adds [load_molecule_topology_mmffl(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:178:0-257:96) which:
  - builds topology through [MMFFL.build_topology()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:123:4-204:9)
  - merges selected bond classes into XPDB fixed-size arrays
  - optionally exports MOL2 via [atomicUtils.save_mol2](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:1710:0-1769:95)
  - prints full topology report
- Adds [plot_topology_3d(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:260:0-296:14) (matplotlib 3D view, colors bonds by class)
- Adds CLI flags:
  - `--type_source table|uff`
  - `--add_epair 0/1`
  - `--add_pi 0/1`
  - `--two_pi 0/1`
  - `--topology_only 0/1`
  - `--mol2_out path`

## Test results (H2O and OCH2)

### H2O
Command used:
```bash
python -u pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py \
  --molecule cpp/common_resources/xyz/H2O.xyz \
  --use_mmffl 1 --type_source table \
  --enable_angles 1 --add_epair 1 --add_pi 0 \
  --topology_only 1 --mol2_out /tmp/H2O_topo.mmffl.mol2 --print_buffers 1
```

Output summary:
- Types:
  - O: `nbond=2 npi=0 nepair=2` (passes octet rule)
- Topology:
  - `n_real=3`, `n_all=5` (2 epair dummies)
  - `bonds_primary`: O–H x2
  - `bonds_angle`: H–H (one)
  - `bonds_epair`: O–E x2 plus E–H “orth” bonds

MOL2 [/tmp/H2O_topo.mmffl.mol2](cci:7://file:///tmp/H2O_topo.mmffl.mol2:0:0-0:0) shows 5 atoms and 9 bonds (including dummy bonds) and looks consistent.

### OCH2
Command used:
```bash
python -u pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py \
  --molecule cpp/common_resources/xyz/OCH2.xyz \
  --use_mmffl 1 --type_source table \
  --enable_angles 1 --add_epair 1 --add_pi 0 \
  --topology_only 1 --mol2_out /tmp/OCH2_topo.mmffl.mol2 --print_buffers 1
```

Output summary:
- Types:
  - O: `nbond=1 npi=1 nepair=2` (passes)
  - C: `nbond=3 npi=1 nepair=0` (passes)
- Topology:
  - `n_real=4`, `n_all=6` (2 epair dummies)
  - `bonds_angle` includes (0,2), (0,3), (2,3) consistent with angle-pair expansion around central atoms
  - `bonds_epair`: O–E and E–neighbor stiffeners

MOL2 [/tmp/OCH2_topo.mmffl.mol2](cci:7://file:///tmp/OCH2_topo.mmffl.mol2:0:0-0:0) contains the expected dummy atoms and bonds.

## Notes / remaining gap (Step 4)
You asked to also try Phase 2 (GPU reindex / local topology mapping). I haven’t implemented that yet. The key missing pieces are:

- **Collision exclusion for dummies** in [XPDB_new.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB_new.cl:0:0-0:0) (and/or via a mask buffer) so epair/pi dummies only participate in bond constraints, not collision constraints.
- **Download and validate** `bond_indices_local`, `ghost_indices_flat`, `ghost_counts` after `build_local_topology` to CPU-check global↔local mapping.

If you want me to proceed, I’ll:
- Add `num_real_atoms` argument to the OpenCL kernels and skip collision loops for `gid >= num_real_atoms`.
- Add a debug download function in [XPDB_new.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB_new.py:0:0-0:0) to read the local topology buffers and a CPU checker in the test script.

## Status
- **Done:** Steps **1–3** (MMparams table typing, MMFFL-based topology builder with angle/pi/epair options, 3D plot + MOL2 export, tested on H2O/OCH2).
- **Pending (optional):** Step **4** (GPU reindex validation + dummy collision exclusion).