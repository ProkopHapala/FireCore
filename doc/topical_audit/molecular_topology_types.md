# Molecular Topology — Type Assignment

> Cross-language audit of atom type assignment: from bonding topology to C_sp2/C_sp3/O_hydroxyl etc. See [base topology](molecular_topology.md) for graph representations and [interactive codemap](https://windsurf.com/codemaps/692593e6-1efe-495f-bbf6-2ad291a285c9-fe86ab10a43f3d18).

## Quick Navigation

| I want to... | Go to |
|--------------|-------|
| Understand bond finding or graph representations | [molecular_topology.md](molecular_topology.md) |
| Build/edit molecules interactively | [molecular_topology_editors.md](molecular_topology_editors.md) |
| Understand the type assignment pipeline | [The Type Assignment Pipeline](#the-type-assignment-pipeline) below |
| Compare type assignment across languages | [Cross-Language Parity](#cross-language-parity) |

---

## The Type Assignment Pipeline

All three languages implement the same pipeline:

```
Bond Topology → Valence Count → Hybridization (sp1/sp2/sp3) → Subtype → Parameters
     ↑                ↑                  ↑                      ↑          ↑
  autoBonds()    nSigma, nPi       npi/nep              C_sp2, C_sp3   l0, k, Ass
```

1. **Bond Topology**: Distance-based bond detection (see [base doc](molecular_topology.md))
2. **Valence Count**: Count sigma bonds, pi orbitals, electron pairs
3. **Hybridization**: sp1 (2 domains), sp2 (3 domains), sp3 (4 domains)
4. **Subtype**: Map (element, npi, nep) → specific type (C_sp2, O_hydroxyl, etc.)
5. **Parameters**: Look up bond lengths, angles, force constants from tables

---

## Parameter File Formats

All implementations share the same `.dat` file formats. The **canonical default parameter files** are located in:

```
tests/tUFF/data_UFF/
├── ElementTypes.dat   # Element symbols, atomic numbers, radii, vdW params
├── AtomTypes.dat      # Hybridized atom types (C_sp2, O_hydroxyl, etc.)
├── BondTypes.dat      # Bond length & stiffness by atom type pair
├── AngleTypes.dat     # Equilibrium angle & stiffness by atom type triple
└── DihedralTypes.dat  # Torsion parameters by atom type quadruple
```

These files are loaded by C++ (`MMFFparams.h`), Python (`OCL/MMFF.py`, `atomicUtils.py`), and JavaScript (`MMParams.js`). When modifying parameters, update this central directory — copies in other locations should be considered derived or legacy.

### ElementTypes.dat
| Field | Description |
|-------|-------------|
| `name` | Element symbol (C, H, O, etc.) |
| `iZ` | Atomic number |
| `neval` | Valence electrons |
| `valence` | Valence |
| `piMax` | Maximum pi orbitals |
| `color` | RGB color for rendering |
| `Rcov` | Covalent radius (used for bond detection) |
| `RvdW` | van der Waals radius |
| `EvdW` | vdW energy parameter |
| `Quff`, `Uuff`, `Vuff` | UFF parameters |
| `Eaff`, `Ehard`, `Ra`, `eta` | QEq parameters (optional, cols 12-15) |

### AtomTypes.dat
| Field | Description |
|-------|-------------|
| `name` | Type name (C_sp2, O_hydroxyl, etc.) |
| `parent_name` | Parent type for fallback |
| `element_name` | Element symbol |
| `epair_name` | Electron pair type |
| `valence` | Valence |
| `nepair` | Number of electron pairs |
| `npi` | Number of pi orbitals |
| `sym` | Symmetry |
| `Ruff`, `RvdW`, `EvdW`, `Qbase`, `Hb` | General params |
| `Ass`, `Asp`, `Kss`, `Ksp`, `Kep`, `Kpp` | MMFF params (cols 13-18) |

### BondTypes.dat
`typeA typeB order l0 k` — bond length and stiffness for atom type pairs.

### AngleTypes.dat
`typeA typeB typeC ang0 k` — equilibrium angle and stiffness.

---

## Parameter Loading

Three separate parsers for the same file format:

| Language | Class | Load Method | File |
|----------|-------|-------------|------|
| **C++** | `MMFFparams` | `loadElementTypes()`, `loadAtomTypes()` | `MMFFparams.h:280` |
| **Python** | `MMparams` | `read_element_types()`, `read_atom_types()` | `MMFF.py:114` |
| **JavaScript** | `MMParams` | `parseElementTypes()`, `parseAtomTypes()` | `MMParams.js:243` |

**C++**: Custom `sscanf` parsing, stores in `std::vector<ElementType>`, `std::vector<AtomType>`.

**Python**: Python `open()` + line splitting, stores in dicts.

**JavaScript**: `fetch()` + `split('\n')`, stores in objects. WebGPU version has additional `logger` integration and `VERBOSITY_LEVEL` gating.

**Redundancy**: Three separate implementations of identical parsing logic. Consider JSON intermediate format or keep separate for language-specific needs.

---

## Type Resolution

Resolving an element symbol to an atom type:

| Language | Function | Algorithm |
|----------|----------|-----------|
| **C++** | `MMFFparams::assignSubTypes()` | Map (iZ, npi, nep) → subtype index |
| **Python** | `MMFF.initAtomProperties()` | Calculate npi/nep from atom_types dict |
| **JavaScript** | `MMParams.resolveTypeNameTable()` | Symbol → exact key → symbol+`_` → first matching element_name |

**C++ subtype assignment** (`MMFFparams.h:533`):
```cpp
// Example: Carbon (iZ=6)
npi=0, nep=0 → C_sp3
npi=1, nep=0 → C_sp2
npi=2, nep=0 → C_sp1
```

**JavaScript resolution priority** (`MMParams.js:155`):
1. Exact symbol key (e.g., `'C'`)
2. Symbol + `'_'` (e.g., `'C_'`)
3. First atom-type whose `element_name` matches

---

## Hybridization Detection (npi/nep)

### C++

`MM::Builder::assignSp3Type()` (`MMFFBuilder.h:445`):
- Iterates atom confs
- Counts `npi` from pi-orbital flags
- Counts `nep` from `nEPair`
- Uses `MMFFparams::assignSubTypes()` to map → subtype

`AtomConf.fillIn_npi_sp3()`: Automatically fills in npi for sp3 atoms based on neighbor count.

### Python

`MMFF.initAtomProperties()` (`pyBall/OCL/MMFF.py:51`):
```python
def initAtomProperties(mol, atom_types, capping_atoms={'H'}):
    # Calculate npi from atom_types dict
    npi_list = [...]
    # Calculate nep from valence
    nep_list = [...]
    # Mark isNode (non-H, non-capping)
    isNode = [...]
```

Key detail: Python separates **node atoms** (heavy atoms with pi orbitals) from **capping atoms** (H, etc.). This separation is crucial for the `toMMFFsp3_loc()` conversion.

### JavaScript

`MMFFLTopology.computePiOrientations()` (`MMFFLTopology.js:248`):
- Computes pi-orbital directions from neighbor geometry
- Uses `computeAtomiPiDirectionFromNeighs()` — sum of normalized cross products
- Propagates low-norm pi dirs from neighbors (min_norm=0.7, max_iter=4)
- Falls back to axis direction for single-neighbor atoms

---

## Pi-Fragment Detection

Identifying conjugated/aromatic systems for pi-pi interactions.

| Language | Function | File |
|----------|----------|------|
| **C++** | `MM::Builder::assignPiFragments()` | `MMFFBuilder.h:398` |
| **Python** | Not explicit (part of type assignment) | — |
| **JavaScript** | `computePiOrientations()` (direction only) | `MMFFLTopology.js:248` |

**C++**: `assignPiFragments()` walks the graph to find connected pi-orbital networks. Sets `piFragment` IDs on atoms.

---

## Bond Parameter Assignment

Assigning bond length (`l0`) and stiffness (`k`) from atom types and bond order.

| Language | Function | File |
|----------|----------|------|
| **C++** | `MM::Builder::assignBondParams()` | `MMFFBuilder.h:554` |
| **Python** | Implicit in `toMMFFsp3_loc()` | `MMFF.py` |
| **JavaScript** | `buildMMFFLTopology()` → bonds in `topo.bonds_linear` | `MMFFLTopology.js` |

**C++**: `params->getBondParams(ai.type, aj.type, order, b.l0, b.k)` — looks up (typeA, typeB, order) in bond parameter table. Considers pi-bonds and lone pairs for pi-pi alignment.

**JavaScript**: Angle constraints are generated as linearized bonds via `buildAngleBonds()` using law of cosines:
```javascript
function lawOfCosines(rab, rbc, cosTheta) {
    const t = rab * rab + rbc * rbc - 2.0 * rab * rbc * cosTheta;
    return Math.sqrt(Math.max(0.0, t));
}
```

---

## Cross-Language Parity

### Type Assignment Equivalents

| Step | C++ | Python | JavaScript |
|------|-----|--------|------------|
| Load params | `MMFFparams::loadAtomTypes()` | `MMparams.read_atom_types()` | `MMParams.parseAtomTypes()` |
| Count npi/nep | `AtomConf.npi/nEPair` | `initAtomProperties()` | `computePiOrientations()` |
| Map → subtype | `assignSubTypes()` | Implicit in type dict | `resolveTypeNameTable()` |
| Assign bond params | `assignBondParams()` | Implicit in `toMMFFsp3_loc()` | `buildMMFFLTopology()` |

### Known Differences

1. **C++** has explicit `assignPiFragments()` for conjugated systems; Python/JS handle this implicitly
2. **JavaScript** `computePiOrientations()` is more geometry-aware (uses actual neighbor positions) than C++/Python (topology-only)
3. **Python** explicitly separates node vs capping atoms; C++ and JS handle this differently
4. **C++** `AtomConf` supports multiple confs per atom; Python/JS do not

---

## File Index

### C++
- `cpp/common/molecular/MMFFBuilder.h` (837 lines) — `assignSp3Type`, `assignPiFragments`, `assignBondParams`
- `cpp/common/molecular/MMFFparams.h` (615 lines) — parameter loading, `assignSubTypes`

### Python
- `pyBall/OCL/MMFF.py` (1107 lines) — `initAtomProperties`, `toMMFFsp3_loc`
- `pyBall/atomicUtils.py` (2186 lines) — bond/angle utilities

### JavaScript
- `web/molgui_webgpu/MMFFLTopology.js` (826 lines) — `buildMMFFLTopology`, `computePiOrientations`, `buildAngleBonds`
- `web/molgui_webgpu/MMParams.js` (524 lines) — `parseAtomTypes`, `resolveTypeNameTable`, `bondLengthEstimate`

---

## See Also

- [Topical Audit Index](topical_audit.md) — priority ranking, dependency graph, missing topics
- [Base Topology](molecular_topology.md) — graph representations, bond finding, rings, bridges
- [Editors & Editing](molecular_topology_editors.md) — GUIs, crystal building gap analysis, consolidation roadmap
- [GUI Feature Audit](gui_audit.md) — visualization & editor feature matrices, VisPy consolidation plan
- [Intramolecular Forcefields](intramolecular_forcefields.md) — FF evaluation using assigned types
- [Interactive Codemap](https://windsurf.com/codemaps/692593e6-1efe-495f-bbf6-2ad291a285c9-fe86ab10a43f3d18) — visual navigation
