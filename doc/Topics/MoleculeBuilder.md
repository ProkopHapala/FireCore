## USER

I want to make simple graphical user interface for placing molecules on ionic substrate. My usin substrate fro starter is NaCl crystal, but have problem that traditional crystalopgraphy scemes does not generale really stable surfaces of ionic crystal. Basicall for me NaCl is simple cubic with laternating Na and Cl ions where the nearest neighbor Na-Cl bonds are 2.82A long and are aligned along cartesian axis x,y,z this is something different than what I ge from ASE or CIF when using lattice constant of NaCl 5.62A (I think it is some diagonal not sure if face diagonal or body diagonal of the cube), I have script to geenrate it, I would like to generalize it.

Then I want to generate step edge which is this surface but tilted so each is it offset acetly by one layer along the superlattice, analyse how it is done in this script
@NaCl_step_.py @NaCl_step_2.py @NaCl_step.py 
let me know what are the differences between these scripts?

I would like to generalize this. making it possible to make multiple layers, 

then I want to place row of molecules along this substrate under certain angle. One think is the angle of the row (two lattice vectros), the other is the rotation of the molecule. Let me control both independnetly. 

Make sure I have interactive visualization and it is fast, using matplotlib 3D view and CLIP (fast redraw)

as a reference take this molecular placer @MolecularPlacer.py , this is very good xample of similar program just made for trinagular lattices, I want sometthing similar for cubic lattices as stepd edges. Think deeply how to make this really user friendly and feature full but also lightweight, performant, and resuing existing functionality.

Plese use this AtomicSystem as underlaying representation of the molecule and surface  (molecule and surface should be two independent instances), @AtomicSystem.py 

use AtomicSystem also for saving and loading of the molecule topology

now your job is investigate provided probramgs, make detailed notes about usefull functions properties and features, and formulate more specific and systematics goals and also a plan.


---

## GPT-5.1-Codex

### Differences between NaCl step scripts
- [tests/pyutils/NaCl_step.py](cci:7://file:///home/prokop/git/FireCore/tests/pyutils/NaCl_step.py:0:0-0:0) (@tests/pyutils/NaCl_step.py#1-58): ASE-based, builds alternating Na/Cl step by skewing x–z (`ax=-1/nx`, shift z by `ax*x`, shift x by `-ax*z`, add a +a height jump for x>Lx/2). Writes an ASE `Atoms`, sets cell, centers with vacuum, saves XYZ.
- [tests/pyutils/NaCl_step_.py](cci:7://file:///home/prokop/git/FireCore/tests/pyutils/NaCl_step_.py:0:0-0:0) (@tests/pyutils/NaCl_step_.py#1-65): Pure writer to XYZ with charges, same geometry logic, different grid (nx=15, ny=8), writes lvec header and per-atom charge ±Q0.
- [tests/pyutils/NaCl_step_2.py](cci:7://file:///home/prokop/git/FireCore/tests/pyutils/NaCl_step_2.py:0:0-0:0) (@tests/pyutils/NaCl_step_2.py#1-64): Same as `_` but grid nx=13, ny=12 (matches first script), otherwise identical to `_` (XYZ with charges, same header).
Core algorithm identical: cubic Na/Cl alternating on simple cubic, step made by linear skew + half-slab height offset past midline in x.

### Useful pieces in AtomicSystem
- File loading/saving: [AtomicSystem(fname=...)](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) reads xyz/mol/mol2 and parses `lvs` comment; [saveXYZ](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:80:4-86:144), [save_mol2](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:91:4-92:75), [save_mol](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:88:4-89:80).
- Topology helpers: [findBonds](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:135:4-139:29), [findHBonds](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:141:4-142:149), [_bond_adjacency](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:296:4-306:18) + selection grow/shrink/connected.
- PBC replication: [clonePBC(nPBC)](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:217:4-259:87) builds tiled system; [symmetrized](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:261:4-265:104), [selectSubset](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:267:4-286:87).
- Orientation utilities: [makeRotMat](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:369:4-372:38), [orient_vs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:381:4-387:48), [orient](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:389:4-395:69), [orientPCA](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:397:4-398:43).
- Geometry math comes from `pyBall.atomicUtils`; use AtomicSystem as the in-memory representation for both substrate and molecule, and for I/O.

### Useful pieces in MolecularPlacer (template for GUI features)
- Matplotlib interactive GUI with sliders/buttons; efficient redraw via cached colors/sizes.
- Rotation sliders Rx/Ry/Rz and application ([apply_rotation](cci:1://file:///home/prokop/git/FireCore/pyBall/MolecularPlacer.py:381:4-393:42)).
- Lattice-based positioning: Hex lattice helper (triangular), cell indices, offsets; dual-flower duplication; PBC replication slider.
- Export: save xyz, mol2 (with bond detection), save poses/state.
- 2-panel view (XY + XZ), grid toggles; status text.

### Proposed goals (systematic)
1. **NaCl step substrate generator (cubic)**  
   - Parameterize `a`, `nx, ny, nz`, step offset direction/height, optional charge assignment.  
   - Produce an [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) for substrate with `lvec` set; allow multiple step edges and layer count.
2. **Generalized step builder**  
   - Allow multiple terraces/steps (integer number of layers per terrace).  
   - Control tilt/offset vector (default along +x) and step height (integer layers).  
   - Optionally add vacuum padding in z and trimming window.
3. **Molecule placement on substrate**  
   - Two independent [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) instances (substrate, molecule).  
   - Controls: row direction (two lattice vectors or angle in surface plane), row spacing, per-molecule rotation, per-row tilt/offset.  
   - Independent molecule orientation (Rx/Ry/Rz) vs. row direction.
4. **Interactive GUI (matplotlib 3D fast redraw)**  
   - Reuse MolecularPlacer widget patterns: sliders for rotations, lattice params, offsets; buttons for load/save substrate/molecule; toggles for grid/axes.  
   - Add 3D view with quick scatter redraw (clip/fast path, avoid heavy artists).  
   - Live update of both systems, with optional PBC replication preview for substrate.
5. **I/O and persistence**  
   - Load/save substrate/molecule via [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) (xyz with `lvs`, mol2 for bonds).  
   - Export combined placement and rigid-body poses.
6. **Stability/validation**  
   - Assertions on lattice sizes >0, correct alternation, step height integral; loud errors on missing data.

### High-level plan
1. Extract the NaCl step generator into a reusable function returning [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) (cubic simple alternating, step via skew + offset), with charge option and lvec.
2. Extend generator to multi-step/terrace: parameters for number of steps, terrace width, step height (layers), direction.
3. Design a cubic lattice helper (analogous to HexLattice) exposing lattice vectors, cell origins, and step geometry.
4. Build a GUI (matplotlib) that:
   - Loads molecule/substrate as separate [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78).
   - Controls molecule rotations and row direction/spacing independently.
   - Places a row of molecules on the surface at chosen angle; optional PBC preview.
   - Fast redraw in 3D and side views.
5. Add save/export paths (xyz/mol2 for combined structure, plus poses), reuse `atomicUtils` for bonds if needed.
6. Keep code lightweight (numpy + matplotlib), fail loudly on invalid inputs.