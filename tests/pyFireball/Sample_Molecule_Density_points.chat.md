## USER

I want to sample prepare samples of electron density from Fireball DFT program at selected importaint points (like centers of atoms, centers of bonds, electron pairs, sigma holes (apexes of capping atoms) etc. This taks has several sub tasks.

1) generate the geometry of samples from molecular geometry. We hsoul use @AtomicSystem.py with make_epair_geom, similar to as is illustrated in @run_add_epairs.py , but we should some more - extension of apexes etc. we should make nice CLI where one can desing how exactly to gnerate these samples.
2) we should test it for few molecules. Some in this dectionary
/home/prokop/git/FireCore/tests/pyFireball/molecules
we will start with CH4, H2O, NH3, later C2H4, CH2O, CH2NH
3) we should sample the density at those points this is done using python interface to fireball fotran @libFireCore.f90 @FireCore.py 
firecore_dens2points, see example in @plot_h2o_density_plane.py 
4) when all of this is tested and we are fine, we can try scale it up to larger training sets where we have distorted geometries of those molecules. see @distort_molecule.py 

now please read all the provided files and make detailed notes how to use them, that means (1) how to develope the sampling framework using these tool, what excat function and classes to use, note all relevant technical details (2) how to use it, how it should work in the end from user perspective (3) plan phases of devlompment in detailed steps.


---

## GPT5-Codex

Here’s a concise technical brief based on the provided files and rules, plus a proposed development plan for the sampling framework.

## 1) Geometry generation (AtomicSystem & epairs)
- [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) core behaviors:
  - Load geometries from `.xyz/.mol/.mol2` with debug prints; see initializer @pyBall/AtomicSystem.py#19-80.
  - [findBonds](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:135:4-139:29) → fills `self.bonds` using covalent radii; [neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23) builds neighbor map `self.ngs` (dict per atom) @pyBall/AtomicSystem.py#136-175.
  - [preinitialize_atomic_properties](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:438:4-469:72) sets `qs` (valence e⁻), `Rs` (vdW radii), `aux_labels`, and ensures neighbors exist @pyBall/AtomicSystem.py#439-471.
- Electron pair placement:
  - `VALENCE_DICT` defines expected bond count and lone pairs per element (`O: (2,2)`, `N: (3,1)`).
  - [add_electron_pairs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52) computes `nb`, `npi=nb-nsigma`, calls [make_epair_geom](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:903:4-947:52) @pyBall/AtomicSystem.py#875-889.
  - [make_epair_geom(i, npi, nb)](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:903:4-947:52) builds directions from neighbor vectors:
    - NH₃-like (nb=3, npi=0): place epair along flipped plane normal.
    - H₂O-like (nb=2, npi=0): two tetrahedral directions from bisector/cross with coefficients cc=√(1/3), cb=√(2/3).
    - =N⁻-like (nb=2, npi=1): one pair opposite bisector.
    - =O-like (nb=1, npi=1): uses [get_atomi_pi_direction](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:889:4-901:32) of neighbor to create two pairs around pi-plane @pyBall/AtomicSystem.py#904-945.
  - [place_electron_pair(i, direction, distance=0.5, ename='E', atype=200, qs=0.0, Rs=1.0)](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:949:4-975:78) appends a new “atom” for the pair, updates `apos`, `atypes`, `enames`, `qs`, `Rs`, and (if present) `bonds` and `ngs` @pyBall/AtomicSystem.py#939-976.
- Example usage (adds epairs, saves XYZ): @tests/pyFireball/run_add_epairs.py#21-28.
- Other geometry utilities relevant for sampling locations:
  - Neighbor queries, angle/dihedral finders, selection helpers, rotations, translations, cloning with PBC; see @pyBall/AtomicSystem.py#151-428.
  - Distortion/scanning examples (distance/angle/rotation/scale sweeps) show how to iterate geometries and optionally add epairs or run Fireball: @tests/pyFireball/distort_molecule.py#23-255.

## 2) Density sampling via Fireball Fortran (Python wrapper)
- Python interface is in [pyBall/FireCore.py](cci:7://file:///home/prokop/git/FireCore/pyBall/FireCore.py:0:0-0:0) (ctypes).
- Basic flow (as in [plot_h2o_density_plane.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/plot_h2o_density_plane.py:0:0-0:0)):
  1. [fc.setVerbosity(0)](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:122:0-123:61)
  2. [fc.initialize(atomType=mol.atypes, atomPos=mol.apos, verbosity=…)](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:782:0-789:15) → calls [preinit()](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:204:0-205:33) + [init](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:129:0-136:57).
  3. Run SCF: [forces, energies = fc.evalForce(mol.apos, nmax_scf=100)](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:142:0-146:21)
  4. Sample density: [densities = fc.dens2points(points_3d, f_den=1.0, f_den0=0.0)](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:344:0-349:21) @tests/pyFireball/plot_h2o_density_plane.py#127-133.
- [dens2points](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:344:0-349:21) is bound in [FireCore.py](cci:7://file:///home/prokop/git/FireCore/pyBall/FireCore.py:0:0-0:0) (via generated argDict; not shown in snippet but present).
- Grid example for H₂O:
  - Build plane from atom positions, define orthonormal `x_axis`, `y_axis`, create mesh, map to 3D points, call [dens2points](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:344:0-349:21), reshape to 2D, plot.
- Other wrappers:
  - [evalForce](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:142:0-146:21), [relax](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:152:0-166:17), [set_export_mode](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:667:0-668:57), [get_HS_dims](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:455:0-484:5), etc., for more advanced data export. Initialize once per geometry; reuse for multiple point batches.

## 3) Test molecules and distortion scaffolding
- Molecule files live in `tests/pyFireball/molecules/`.
- Distortion script ([distort_molecule.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/distort_molecule.py:0:0-0:0)) shows:
  - How to scan distances/angles/rotations, optionally add epairs and/or call Fireball to compute energies.
  - Saving successive XYZ frames for training data.

## 4) User-facing goal
- Provide a CLI to design sampling points derived from molecular geometry:
  - Sources: atom centers, bond midpoints, epair positions (lone pairs), sigma-hole apexes (extend along bond/coordination axes), optional custom offsets/extensions.
  - Output: augmented point set (with labels/roles) and an XYZ (or CSV/NPZ) containing sample locations.
  - Then run Fireball to compute density at those points via [dens2points](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:344:0-349:21), storing densities alongside coordinates and labels.

## 5) Relevant technical considerations (from code)
- Always call [findBonds()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:135:4-139:29) then [neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23) (implicit in [preinitialize_atomic_properties](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:438:4-469:72)) before [add_electron_pairs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52).
- Epairs are appended as new atoms with ename `'E'`, atype `200`, default distance `0.5 Å`. We can parameterize distance/atype/name in CLI.
- [get_atomi_pi_direction](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:889:4-901:32) uses neighbors’ cross products; requires at least two neighbors.
- [dens2points](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:344:0-349:21) expects a contiguous (n,3) float64 array; returns flat densities length n.
- FireCore initialization uses integer `atypes` (atomic numbers) and float64 positions; ensure contiguous arrays.

## 6) Proposed phased development plan
1) **Geometry sampling CLI (no density yet)**
   - Build a script (e.g., `sample_points.py`) to load XYZ, find bonds, build neighbors, optionally add epairs, and generate point sets:
     - Atom centers.
     - Bond centers (with optional offsets along bond).
     - Epairs (reuse [add_electron_pairs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52), but expose distance/atype/ename as options).
     - Sigma-hole/apex points: extend beyond heavy atoms along vector(s) derived from neighbors or bond directions (parameterize extension length).
   - Output: augmented XYZ for visualization + tabular point list with labels/roles.

2) **Density sampling integration**
   - Add options to initialize FireCore, run SCF once, and call [dens2points](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:344:0-349:21) on the assembled point set.
   - Store results in NPZ/CSV with coordinates, labels, and density values. Keep XYZ for inspection.

3) **Validation on small set**
   - Test CH4, H2O, NH3 using molecules in `tests/pyFireball/molecules/`. Visual sanity: epairs appear correctly; densities finite.
   - Add CH2O, CH2NH, C2H4 after initial pass.

4) **Extension/advanced sampling**
   - Add “apex extension” controls (length scaling, directional choice: bisectors, opposite lone pair, etc.).
   - Allow multiple templates (e.g., along each bond, along opposite of bond, along angle bisector).

5) **Scaling to distorted sets**
   - Reuse [distort_molecule.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/distort_molecule.py:0:0-0:0) patterns to generate batches of distorted geometries; loop through files, sampling points + densities into datasets.

6) **Optional visualization**
   - Minimal plotting (matplotlib) for 2D slices or histograms of densities; keep plotting separate from core compute.

## 7) How a user should run it (intended UX sketch)
- Minimal:
  ```
  python sample_points.py --xyz tests/pyFireball/molecules/H2O.xyz \
      --atoms --bonds --epairs --apex \
      --apex-length 0.8 --epair-distance 0.5 \
      --out-points out/h2o_points.npz --out-xyz out/h2o_augmented.xyz
  ```
- With density:
  ```
  python sample_points.py --xyz tests/pyFireball/molecules/H2O.xyz \
      --atoms --bonds --epairs --apex \
      --fireball --nscf 100 \
      --out-points out/h2o_points.npz --out-dens out/h2o_dens.csv
  ```

## 8) Key functions/classes to reuse
- Geometry:
  - [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) constructor, [findBonds](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:135:4-139:29), [neighs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23), [add_electron_pairs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52), [place_electron_pair](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:949:4-975:78) (parameterize distance), selection utilities.
  - `atomicUtils` helpers (`normalize`, bond finding, rotation, etc.).
- Density:
  - [FireCore.initialize](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:782:0-789:15), [evalForce](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:142:0-146:21), [dens2points](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:344:0-349:21).

Let me know if you want me to draft the CLI structure and argument set, or outline the data schema for saved samples.

---

## USER

Now lets investigate in details the idea of how to generate the geometry which is more questionable.
1) you do not need to care about number of pi-orbitals, or free electron pairs on the atom, this should be able @AtomicSystem.py to figure out automatically - it builds the bonding topology, from that assigns types sp,sp2,sp3 ... 
2) making bond centers is perhaps not implemented, we should supplement it.
3) we need to augument the sampling by some samples outside molecule in reasonable distance where density should go already to zero, to make sure density which we fit in the end is confined to proximity of atoms. It is questionable where exactly do these augumentation points but perahs it is bets to do it by prolongation of epairs to some distance from molecule (e.g. 6 Angstroems? should be selectable by CLI)
4) then it perhaps make sence to add as augumentation inversions of all bonds and electron pairs (like vectro directed to oposite side)
6) still may be usefull to add something like centers of all angles, where by angles we mean pairs of bonds, and by bonds we mean not only bonds between atoms, but also bonds between atom and its electron pair (electron pair conidered as dummy atom)

all these augumentation layers should be selectable independncetly over the CLI, and probably it would be good to provided certain radius (as absolute number of multiple of atomic radius?)

Try to think it through more and make detailed notes


---

## GPT5-Codex

Here are detailed notes on geometry generation and augmentation options, with how to implement and expose them via CLI. No code changes yet.

### Baseline capabilities (AtomicSystem)
- Topology & hybridization: [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) already builds neighbors ([findBonds](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:135:4-139:29), [neighs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23)) and derives epair directions automatically via [add_electron_pairs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52)/[make_epair_geom](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:903:4-947:52); you do not need to manually supply pi counts. See @pyBall/AtomicSystem.py#136-175 and epair block @pyBall/AtomicSystem.py#874-976.
- Epairs are added as dummy atoms (`ename='E'`, `atype=200`, default distance=0.5 Å). Directions depend on neighbor geometry (sp/sp2/sp3-like).

### New sampling layers to add
Each should be toggleable independently via CLI flags, with distance controls either absolute (Å) or scaled by atomic radius.

1) **Bond centers (midpoints)**
   - For each bond `(i,j)` in `self.bonds`, add point `p = 0.5*(r_i + r_j)`; optionally allow offset along bond axis: `p = r_i + t*(r_j - r_i)` with `t in [0,1]` (default 0.5).
   - Option: multiple samples per bond (e.g., at 0.25, 0.5, 0.75).
   - CLI: `--bond-centers` and `--bond-fracs 0.25,0.5,0.75` (or `--bond-offset 0.5`).

2) **Outer “far” samples to capture decay**
   - Idea: place points outside molecule where density should be ~0.
   - Practical implementation: prolong existing epair directions and/or bond directions to a user-defined distance (e.g., 6 Å from reference atom).
   - For each epair direction `d` on atom i, add `r_i + L_far * d`; default `L_far=6.0 Å`, CLI `--far-epairs --far-dist 6.0`.
   - Also consider atom-centered radial samples: `r_i + L_far * u_k` for a small set of quasi-uniform directions (e.g., tetrahedral set) to ensure coverage; CLI `--far-shell --far-dist 6.0 --far-normals tetra`.

3) **Inversions (opposite direction of bonds/epairs)**
   - For each bond `(i,j)`, also add mirrored point: if you have midpoint or along-bond point at `p = r_i + t*(r_j - r_i)`, add inverted `p_inv = r_i - t*(r_j - r_i)` (or vice versa relative to j). Good for symmetry and decay sampling.
   - For each epair at `r_ep = r_i + d`, add `r_i - d` as inversion.
   - CLI: `--invert-bonds` (works with chosen bond fractions), `--invert-epairs`.

4) **Angle centers (bond-pair bisectors)**
   - For each angle (i-j-k), define bisector at central atom j: `v1 = norm(r_i - r_j)`, `v2 = norm(r_k - r_j)`, bisector `b = norm(v1+v2)` (guard small norms). Place point `r_j + L_ang * b`, with L_ang as absolute or scaled by avg radius of (i,k) or j.
   - Include angles where one leg is an epair (treat epair as neighbor) once epairs are added.
   - CLI: `--angle-centers --angle-dist 0.7` (or `--angle-scale 0.8 --angle-scale-mode=avgR|Rj`).

5) **Sigma-hole/apex extensions**
   - For atoms with multiple neighbors, extend along “anti-bond” directions:
     - Option A: opposite of bond vector: for each neighbor j, add `r_i - L_apex * norm(r_j - r_i)`.
     - Option B: opposite of average neighbor direction: `d = norm(sum_j norm(r_j - r_i)); p = r_i - L_apex * d`.
   - CLI: `--apex` with `--apex-mode=per-bond|avg` and `--apex-dist 0.8` (or scale by Rs).

6) **Atom-centered samples**
   - Always include atom centers (baseline).
   - Optional shells at scaled radii for robustness: `r_i + s*Rs_i*u_k` (u_k from a small direction set). CLI: `--atom-shells --shell-scales 0.5,1.0`.

### Distance / scaling choices
- Support both absolute Å and scaled by atomic radius `Rs_i` (from [preinitialize_atomic_properties](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:438:4-469:72)):
  - CLI pattern: `--dist-mode=abs|scale` plus `--dist-value X`.
  - For bond fractions, keep dimensionless 0–1 along bond; for external/inversion, use abs or scaled.

### Order of operations (to ensure neighbors include epairs)
1) Load geometry.
2) [findBonds](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:135:4-139:29), [neighs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23), [preinitialize_atomic_properties](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:438:4-469:72).
3) If `--epairs` enabled: [add_electron_pairs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52) (appends epair “atoms”, updates bonds/neighs).
4) Build sampling sets per selected layers (atoms, bonds, angles, epairs, apex, far shells, inversions).
5) De-duplicate optionally (within tolerance).
6) Optionally export augmented XYZ (with labels for visualization) and point table for density.

### Labeling
- Keep a parallel list of labels/roles: e.g., `ATOM:C1`, `BOND:C1-H2@0.5`, `BOND_INV:C1-H2@-0.5`, `EPAIR:O3@0.5`, `EPAIR_INV:O3@-0.5`, `ANGLE:O-H-H@0.7`, `APEX:O@avg`, `FAR_EP:O@6.0`.
- These labels travel with the point coordinates and later with density values.

### CLI option sketch (independent toggles)
- `--atoms`
- `--bonds --bond-fracs 0.5` (comma list)
- `--invert-bonds`
- `--epairs --epair-dist 0.5`
- `--invert-epairs`
- `--angle-centers --angle-dist 0.7 --angle-scale-mode abs|Rs`
- `--apex --apex-mode per-bond|avg --apex-dist 0.8 --apex-dist-mode abs|Rs`
- `--far-epairs --far-dist 6.0`
- `--far-shell --far-dist 6.0 --far-normals tetra|octa`
- `--dedup-tol 1e-4`
- `--out-xyz ... --out-points ...` (NPZ/CSV)
- `--with-density --nscf 100` (if sampling density; reuse FireCore.init + dens2points).

### Implementation hooks in AtomicSystem
- Bonds: `self.bonds`, `self.apos`; fractions along `rj - ri`.
- Angles: use [findAngles](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:178:4-181:68) or `ngs` adjacency to enumerate (i-j-k).
- Epairs: after [add_electron_pairs](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52), epairs appear as atoms with ename 'E'; bonds updated if `self.bonds` existed before; you can also track epair origins if needed by storing metadata when placing.
- Radii: `self.Rs` available after [preinitialize_atomic_properties](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:438:4-469:72).

### Data flow to density
- Build point array (n,3) float64 contiguous.
- Init FireCore once per geometry; call [evalForce](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:142:0-146:21) (SCF) then [dens2points(points, f_den=1.0, f_den0=0.0)](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:344:0-349:21).
- Save densities + labels.

If you want, I can draft the exact CLI spec and the sampling builder function signatures next, or propose how to store point metadata (e.g., NPZ with coords, labels, roles, source atoms).

---

## USER


OK, please implement this sampling into new script. For now we are debugging just the sampling, so make script with CLI where user chose molecule (H2O by default), and it will generate the points, and view them using matplotlib (can be switched off), also write them into terminal so you can check the generation, and also save .xyz file for extarnal view (actually you can read the dumped .xyz directly, but maybe it would be good to have labels which type of sample point it is like angle centre, bonds centers, epairs, bond inversions etc. should be distinguishable tyles and classes)...

when you finish coding run it for the test and check output if it looks correct. Then I will test it visualy using the matplolib 3D render. 

maybe it would be good to store not just point but also their connections (bonds) to the host atom. Maybe in that case we should save .mol or .mol2 not just .xyz. Note that subroutines for these files are in @AtomicSystem.py yu shold not reimplement it. Try to maximally reuse existing functionality of @AtomicSystem.py avoiding code duplication.