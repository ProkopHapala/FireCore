# Heterocycle / Kekule Tools

This directory contains two geometry generators and a Kekulé pi-bond optimizer for planar conjugated molecules. They are designed to let you sketch a molecule, obtain a 2-D SVG rendering, and optionally solve/verify a Kekulé bond-order pattern with hard atom-sum constraints.

The geometry is built on a graphene-like hexagonal lattice and stored as a `pyBall.AtomicSystem`, so it can be passed directly to other pyBall / plotUtils workflows.

## Files

- `heterocycle_generator.py` – programatic geometry builder
- `ascii_art_heterocycle.py` – ASCII-art geometry builder
- `heterocycle_generator.py` also contains plotting helpers (`plot_system`, `plot_kekule_phases`)
- `pyBall/KekulePure.py` – pure-Python NumPy Kekulé pi-bond optimizer

---

## 1. The ASCII art interpreter (`ascii_art_heterocycle.py`)

This is the most flexible entry point: you draw the molecule as text and the script turns it into an `AtomicSystem`, optionally runs the Kekulé solver, and renders the result.

### Two input formats

The parser auto-detects the format:

1. **Dimer format** – rows alternate between atom rows and bond rows.
   - Atom rows contain element symbols (always treated as carbon in this format).
   - Bond rows contain `|` (vertical dimer), `-` (horizontal dimer), or `.` (single atom).
   - Example:
     ```
       C C
      | | |
       C C
     ```

2. **Single-atom format** – every non-space character is either an atom or an explicit bond mark.
   - Each row is an atom layer.
   - Adjacent rows are connected by honeycomb topology.
   - Example:
     ```
       C C
       C C C
       C C C
        C C
     ```

### Element symbols and `n_pi` (sp2 vs sp3)

Case is preserved and used by the Kekulé solver to decide how many pi electrons each atom contributes:

- **Uppercase** `C`, `N`, `O` → sp2, `n_pi = 1`
- **Lowercase** `c`, `n`, `o` → sp3, `n_pi = 0`

This is the mechanism used by `make_n_pi()` in `KekulePure.py`. For example, a saturated bridgehead carbon can be written as `c`, while a conjugated one is `C`.

### Explicit bond markers for non-hexagonal rings

The single-atom format also recognizes three special characters:

- `_` – horizontal bond between the nearest atom on the left and the nearest atom on the right in the same row.
- `/` – “shortcut” bond connecting an atom in the row above to an atom in the row below-left. This is used to make pentagons (skipping a hexagon vertex).
- `\` – analogous shortcut, connecting an atom in the row above to an atom in the row below-right.

These marks are **not** atoms; they only add extra bonds. The ASCII examples `purin`, `karbazol`, `biphenylene`, `guanin` show how to use them.

### Hydrogen bonds (`:`)

Both formats recognize `:` as a hydrogen-bond marker between the nearest atom in the row above and the nearest atom in the row below:

```
  O n O
   : :          <-- two H-bonds connecting the two "n" / "O" pairs
  O n O
```

The vertical spacing of the rows is adjusted so that the donor–acceptor distance matches the CLI parameter `--hbond_length` (default 3.0 Å). After the heavy-atom geometry and hydrogen passivation are complete, the `:` pairs are resolved into `(H_index, acceptor_index)` segments and drawn as dashed magenta lines in the SVG.

H-bond donors are identified automatically:
- `N … O` → N is the donor (must carry H)
- `n … N` → the sp3 side (`n`, lowercase) is the donor
- If only one side already has an attached H after passivation, that side is chosen regardless of element.

### Automatic hydrogen passivation (default on)

By default `--hydrogens 1` adds capping H atoms based on pure electron counting, *independent* of H-bonds. The algorithm is:

```
nsigma = nval - n_pi - n_epair
```

For 2nd-period atoms `nval = 4`, `n_epair = 1` for N, `2` for O, `0` for C, and `n_pi` is determined from case:

| Symbol | Hybridization | `n_pi` | `n_epair` | `nsigma` | Result |
|--------|---------------|--------|-----------|----------|--------|
| `N` (uppercase) | sp2 | 1 | 1 | **2** | two sigma neighbors, no H |
| `n` (lowercase) | sp3 | 0 | 1 | **3** | three sigma neighbors, one H if only two heavy neighbors |
| `O` (uppercase) | sp2 (=O) | 1 | 2 | **1** | one double bond (=O), no H |
| `o` (lowercase) | sp3 (-OH) | 0 | 2 | **2** | two sigma neighbors, one H (like -OH) |
| `C` / `c` | sp2 / sp3 | 1 / 0 | 0 | **3 / 4** | standard valence |

This means:
- Uppercase `N` in a ring (two heavy neighbors) stays as NH
- Lowercase `n` with two heavy neighbors automatically gets one H → NH
- Lowercase `o` automatically gets one H → OH
- Uppercase `O` with one heavy neighbor stays as =O, no H added

If the heavy-atom valence already equals the target, **no H is added**.

Turn passivation off with `--hydrogens 0` if you want a bare heavy-atom skeleton.

### Coordinate system

ASCII rows are read from top to bottom, and the plot places the top row at the **largest** Y coordinate. This keeps the visual orientation the same as the text. Rows of the same parity are spaced by `aCC`, while rows of opposite parity are spaced by `aCC/2`, giving the correct zig-zag geometry.

### One-step bond-length relaxation

Because the `_`/`/`/`\` shortcuts can create bonds that are longer or shorter than `aCC`, a fast Jacobi bond-length relaxation is available:

```bash
python3 ascii_art_heterocycle.py -e purin --relax_bonds 10 --relax_bmix 0.5 ...
```

- `--relax_bonds` – number of Jacobi steps.
- `--relax_bmix` – momentum factor from the second step onward.

The implementation is order-independent: it first accumulates all per-atom displacements from the frozen snapshot and only then updates the positions.

---

## 2. The Kekulé pi-bond optimizer (`pyBall/KekulePure.py`)

### Problem statement

For every atom `i` the solver enforces a hard equality:

```
Σ_j x_ij = n_pi[i]
```

where `x_ij` is the **pi** bond order on bond `(i, j)` and `n_pi[i]` is the number of pi electrons on atom `i` (typically 1 for sp2, 0 for sp3). Each bond variable is bounded by `0 ≤ x_ij ≤ 1`.

### Energy model

The objective is a sum of three terms:

1. **Atom valence (constraint)** – `Kval * Σ_i (Σ_j x_ij - n_pi[i])²`. This is the dominant term and is treated as a hard equality in the quadratic solver.
2. **Aromatic stabilization** – `Karo * Σ_bonds (x - 0.5)²`. Favors the aromatic value `x = 0.5`.
3. **Localization** – `Kloc * Σ_bonds (x - target)²`. The `target` is chosen from `{0, 0.5, 1}` depending on the current value of `x` and on whether aromatic bonds are allowed.

### Constrained solver

The default `--solver linsolve` solves the constrained quadratic problem exactly by forming the KKT matrix:

```
[  H   A^T ] [ x ] = [ rhs_x ]
[  A   0   ] [ λ ]   [ n_pi  ]
```

with an active-set method for the `[0, 1]` bounds. After the solve, the code checks `max |A @ x - n_pi|` and fails loudly if it is not on the order of machine precision.

### Two-stage solution

`ascii_art_heterocycle.py` runs the solver in two phases:

1. **Phase 1** (`Kloc = 0.0`) – finds the constrained aromatic/delocalized solution.
2. **Phase 2** (`Kloc = Kloc_final`) – adds the localization bias and drives the bonds toward the nearest discrete minimum (`0`, `0.5`, or `1`).

The right panel in the SVG is the **Phase 2 result**, not a naive rounding. Naive rounding was removed because it breaks the hard constraints.

### Aromatic vs discrete

- `--aromatic 1` (default): allows the 0.5 basin. Use this for benzene, naphthalene, etc.
- `--aromatic 0`: forces localization to 0/1. Use this for molecules like biphenylene where a valid aromatic pattern cannot exist for all atoms.

### Symmetry breaking

For highly symmetric graphs, the deterministic phase-2 solver can remain locked in a symmetric fractional fixed point (e.g. the central four-membered ring of biphenylene). Add a small random perturbation before phase 2:

```bash
python3 ascii_art_heterocycle.py -e biphenylene --aromatic 0 --sym_break 0.2 --seed 1
```

- `--sym_break` – amplitude of uniform noise added to `x` before phase 2.
- `--seed` – reproducible random seed (0 = no seeding).

### Visualization of the KKT matrix

For debugging the solver you can save the KKT matrix:

```bash
python3 ascii_art_heterocycle.py -e naphthalene --kkt 1 --kkt_mode signed --out mol.svg
```

- `--kkt_mode signed` – `seismic` colormap of matrix values.
- `--kkt_mode logabs` – reversed `magma`/`inferno` colormap of log absolute values.

---

## 3. The programatic geometry builder (`heterocycle_generator.py`)

This is the older, more explicit builder. It still works and uses the same Kekulé integration.

### Geometry model

- **x-axis**: zig-zag direction
- **y-axis**: armchair direction
- Rows are either **E (edge)** layers or **D (dimer)** layers
  - First and last rows are **E** layers (single atom per unit cell)
  - Internal rows are **D** layers (one or two atoms per unit cell)

### D-layer modes

| Mode | Description | Effect |
|------|-------------|--------|
| `\|` (default) | vertical dimer | normal hexagonal lattice |
| `.` | single atom | defect / under-coordinated site |
| `-` | horizontal dimer | creates pentagons / heptagons |

### Run

```bash
python3 heterocycle_generator.py --out mol.svg --xyz mol.xyz
python3 heterocycle_generator.py -e substituted --out mol.svg
python3 heterocycle_generator.py -i my_system.py --out mol.svg
```

---

## 4. CLI reference

### Common options

| Option | Default | Description |
|--------|---------|-------------|
| `-e`, `--example` | `naphthalene` | Built-in example (see list below) |
| `-o`, `--out` | `/tmp/kekule/heterocycle.svg` | Output SVG file |
| `--xyz` | `None` | Optional XYZ output |
| `-t`, `--title` | `None` | Plot title |
| `--size` | `120` | Atom marker size |
| `--aCC` | `1.42` | C–C bond length (Å) |
| `--plot` | `1` | Show matplotlib figure (1) or suppress (0) |
| `--kekule` | `1` | Run Kekulé solver and draw bond orders |
| `--kekule_single` | `0` | Draw a single Kekulé panel instead of two-phase |
| `--aromatic` | `1` | Allow aromatic 0.5 bonds (0 = force 0/1) |
| `--solver` | `linsolve` | `linsolve` (KKT) or `gd` (gradient descent) |
| `--Kval` | `50.0` | Atom-sum stiffness (largest) |
| `--Kloc` | `5.0` | Localization stiffness (must be < Kval) |
| `--Karo` | `0.5` | Aromatization stiffness (must be < Kloc) |
| `--kkt` | `0` | Save KKT matrix plot |
| `--kkt_mode` | `signed` | `signed` or `logabs` |
| `--kkt_cmap` | `None` | Override colormap |

### ASCII-art specific options

| Option | Default | Description |
|--------|---------|-------------|
| `--relax_bonds` | `0` | Number of Jacobi bond-length relaxation steps |
| `--relax_bmix` | `0.5` | Momentum factor for relaxation |
| `--sym_break` | `0.0` | Symmetry-breaking noise before phase 2 |
| `--seed` | `0` | Random seed for symmetry breaking |
| `--dump_bonds` | `0` | Print atom list and bond list to stdout |
| `--hydrogens` | `1` | Add passivation H based on electron counting (0 = off) |
| `--hbond_length` | `3.0` | Donor–acceptor distance for `:` H-bonds (Å) |
| `--mol` | `auto` | Save MOL file (`auto` / `off` / explicit path) |

---

## 5. Examples

### Naphthalene (aromatic, two panels)

```bash
python3 ascii_art_heterocycle.py -e naphthalene --out /tmp/kekule/naphthalene.svg
```

### Purine (with shortcut `/` bond)

```bash
python3 ascii_art_heterocycle.py -e purin --out /tmp/kekule/purin.svg
```

### Biphenylene (discrete, symmetry-broken)

```bash
python3 ascii_art_heterocycle.py -e biphenylene --aromatic 0 --sym_break 0.2 --seed 1 \
    --out /tmp/kekule/biphenylene.svg
```

### Relax bond lengths

```bash
python3 ascii_art_heterocycle.py -e purin --relax_bonds 10 --relax_bmix 0.5 \
    --out /tmp/kekule/purin_relaxed.svg
```

### Single-panel Kekulé plot

```bash
python3 ascii_art_heterocycle.py -e cytosin --kekule 1 --kekule_single 1 \
    --out /tmp/kekule/cytosin.svg
```

### Plain geometry with n_pi labels

```bash
python3 ascii_art_heterocycle.py -e cytosin --kekule 0 --out /tmp/kekule/cytosin_plain.svg
```

### Hydrogen bonds with passivation (2NCI dimer format)

```bash
python3 ascii_art_heterocycle.py -e 2NCI --out /tmp/kekule/2nci.svg
```

Lowercase `n` → sp3 NH, uppercase `O` → =O, `:` markers are resolved into dashed magenta H-bond lines.

### Adjust H-bond distance

```bash
python3 ascii_art_heterocycle.py -e 2Quinolone --hbond_length 3.5 \
    --out /tmp/kekule/2quinolone.svg
```

### No hydrogens (bare heavy-atom skeleton)

```bash
python3 ascii_art_heterocycle.py -e NTCDI --hydrogens 0 --out /tmp/kekule/ntcdi.svg
```

### MOL export with correct bond orders

```bash
python3 ascii_art_heterocycle.py -e 2Quinolone --mol /tmp/kekule/2quinolone.mol
```

The MOL file contains all atoms (including H) and bond types: `1` = single, `2` = double, `4` = aromatic. Use `--mol auto` (default) to derive the name from `--out`.

---

## 6. Built-in ASCII examples

`naphthalene`, `naphthalene2`, `fulvalene`, `pentacross`, `biphenyl`, `biphenyl2`, `phenanthrene`, `phenanthrene2`, `perylene`, `perylene2`, `purin`, `purin_x`, `purin_y`, `7azaindol`, `karbazol`, `biphenylene`, `uracil`, `cytosin`, `guanin`, `NTCDA`, `NTCDI`, `TAP`, `2NCI`, `2Quinolone`, `2Quinolinone`, `2purin`.

---

## 7. Output

- `*.svg` – 2-D rendering using `pyBall.plotUtils.plotSystem` (atom colors from `pyBall.elements`)
- `*.xyz` (optional) – Cartesian coordinates via `AtomicSystem.saveXYZ`
- `*.mol` (optional) – MOL V2000 with bond types (`1` = single, `2` = double, `4` = aromatic) and all atoms including H

Internally the geometry is stored as a `pyBall.AtomicSystem` with positions, element names, bonds, and a neighbor list, so it can be passed directly to other pyBall / plotUtils workflows.

The `atoms.hbonds_ascii` attribute holds the resolved H-bond segments as `(H_idx, acceptor_idx)` pairs for downstream analysis or custom plotting.

---

## 8. Tips for building on top of this

- Use `parse_ascii_art(text)` in your own scripts to get an `AtomicSystem`.
- Use `make_n_pi(system)` to obtain the `n_pi` vector from element case.
- Use `KekulePure(system, n_pi=npi, ...)` directly if you need custom control over the solver.
- The solver is intentionally pure NumPy and self-contained: only `pyBall.AtomicSystem` and `pyBall.elements` are required for the geometry part.

### Reusable helpers in `ascii_art_heterocycle.py`

These are the functions the CLI script uses internally; you can call them directly:

- `_build_target_valence(atoms, n_pi)` – returns a `{atom_idx: target_sigma}` dict using `KekuleBackend._target_sigma` for every C/N/O. Works with uppercase and lowercase symbols.
- `resolve_hbond_pairs(atoms)` – resolves `atoms._hbonds_pairs` into `atoms.hbonds_ascii = [(H_idx, acceptor_idx), ...]`. Donor selection is data-driven (prefers the side that actually has an attached H).
- `run_kekule_solver(atoms, ...)` – runs the two-phase solver and returns `{'bo_raw', 'bo_snap', 'n_pi', 'k', 'err', 'report'}`.
- `mol_bond_types(atoms, bo_snap, ...)` – maps Kekulé pi bond orders to MOL bond types (`1`, `2`, `4`).
- `export_mol(atoms, mol_opt, out_path, title, bond_types)` – handles the `auto`/`off`/path logic and calls `save_mol`.
