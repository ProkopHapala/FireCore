---
name: molecular-structure-sync
description: Use when touching anything related to molecular topology, bond orders, rendering, export, or syncing structure across subsystems — ensures AtomicGraph stays authoritative
trigger:
  glob:
    - "**/AtomicGraph.py"
    - "**/MoleculeEditorBackend.py"
    - "**/KekuleExtension.py"
    - "**/SPAMMM_GUI.py"
    - "**/VispyUtils.py"
    - "**/ascii_art_heterocycle.py"
    - "**/KekulePure.py"
    - "**/atomicUtils.py"
    - "**/AtomicSystem.py"
    - "**/topology/**"
    - "**/GUI/**"
---

## AtomicGraph is the Authoritative Molecular Structure

`AtomicGraph` (`spammm/topology/AtomicGraph.py`) is the **single source of truth** for molecular topology. All other representations are **derived** and must sync from it, never the reverse.

### Representations and Their Roles

| Layer | Class | Data format | Role | Direction |
|-------|-------|-------------|------|-----------|
| **Authoritative** | `AtomicGraph` | Object graph (dicts of `Atom`/`Bond`/`Ring` with stable `_id`) | Editable topology, bond orders, hybridization | Read/write here |
| **Bridge** | `MoleculeEditorBackend` | Dense numpy arrays + ephemeral ID mappings | Syncs graph → arrays for rendering/FF/export | Graph → arrays only |
| **Rendering** | `AtomScene` (Vispy) | Dense arrays, picks emit `Atom._id` | 3D visualization, mouse interaction | Read-only consumer |
| **I/O** | `AtomicSystem` | Dense arrays (`apos`, `enames`, `bonds`) | File export/import, FF computation | Read-only consumer |
| **Solver** | `KekulePure` | Pi-bond order arrays | Kekule bond order optimization | Writes results back to graph |

### The Sync Flow

```
AtomicGraph (authoritative)
    │
    ├── MoleculeEditorBackend._sync_sys()  →  AtomicSystem (dense arrays for rendering/export)
    │                                    _atom_ids, _bond_ids (ephemeral mappings)
    │
    ├── MoleculeEditorBackend.get_graph_bond_orders()  →  Vispy rendering (reads Bond.order)
    │
    ├── MoleculeEditorBackend._graph_bond_types_mol()  →  MOL/MOL2 export (reads Bond.order)
    │
    └── KekuleExtension._store_bond_orders_on_graph()  →  writes Bond.order after solving
```

### Critical Rules

1. **Write to graph, read from graph.** If you compute something (bond orders, atom types, positions), store it on `AtomicGraph` objects (`Bond.order`, `Atom.npi`, `Atom.pos`). Never store it in a parallel array on `window` or `backend` that duplicates the graph.

2. **No parallel arrays.** Do not create `window.bond_orders`, `window.bond_order_bonds`, or similar. These are duplicates of `Bond.order` on the graph and will desync. The old `window.bond_orders` / `window.bond_order_bonds` were removed for this reason.

3. **`_sync_sys()` is the export gate.** Before any consumer (rendering, export, FF) reads dense arrays, `_sync_sys()` must have been called. It rebuilds `AtomicSystem` arrays + 4 ephemeral mapping arrays from `graph.to_arrays()`. The mappings are:
   - `_atom_ids[i]` = `Atom._id` at dense index `i` (dense→graph)
   - `_atom_idx_map[id]` = dense index (graph→dense)
   - `_bond_ids[i]` = `Bond._id` at dense index `i` (dense→graph)
   - `_bond_idx_map[id]` = dense index (graph→dense)

4. **`Bond.order` is total bond order** (1.0=single, 1.5=aromatic, 2.0=double). Kekule pi bond order is `Bond.order - 1.0`. When Kekule solver finishes, write `bond.order = 1.0 + pi_order` for pi bonds, `bond.order = 1.0` for sigma-only.

5. **Export must read from graph.** `save_structure()` calls `_graph_bond_types_mol()` which reads `Bond.order` from `_bond_list`. Both `.mol` and `.mol2` paths use this. Never hardcode bond types in export functions.

6. **Rendering must read from graph.** `refresh_view()` calls `backend.get_graph_bond_orders()` which reads `Bond.order` from `_bond_list`. Never read bond orders from a `window.*` attribute.

7. **Stable identity is `Atom._id` / `Bond._id`, not array index.** Dense indices are ephemeral — they change on every `_sync_sys()`. All signals, selections, and callbacks use `_id`. Mappings bridge the two worlds but are rebuilt on every sync.

### Common Mistakes to Avoid

- **Storing solver results on `window`** instead of on `Bond.order` → desync between rendering and export
- **Hardcoding bond types in export** (e.g. `bond_type = 1` in `save_mol2`) → ignores Kekule results
- **Adding export buttons to plugin panels** → export should go through the main `save_structure()` which reads from the graph
- **Creating parallel bond-order arrays** → always use `Bond.order` on the graph
- **Forgetting `_sync_sys()` before reading** `_bond_list` or `_atom_ids` → stale mappings

### Key Files

| File | Role |
|------|------|
| `spammm/topology/AtomicGraph.py` | `Atom`, `Bond`, `Ring` classes, `to_arrays()`, `add_bond(order=...)` |
| `spammm/topology/MoleculeEditorBackend.py` | `_sync_sys()`, `get_graph_bond_orders()`, `_graph_bond_types_mol()`, `save_structure()` |
| `spammm/GUI/SPAMMM_GUI.py` | `refresh_view()` — reads bond orders from graph via backend |
| `spammm/GUI/KekuleExtension.py` | `_store_bond_orders_on_graph()` — writes solver results to `Bond.order` |
| `spammm/AtomicSystem.py` | `save_mol(bond_types=...)`, `save_mol2(bond_types=...)` — pass-through to atomicUtils |
| `spammm/atomicUtils.py` | `save_mol()`, `save_mol2()` — file writers that accept `bond_types` parameter |

### When Adding a New Subsystem or Plugin

1. **Write to graph:** If your plugin computes molecular properties (bond orders, charges, types), store them on `Atom` or `Bond` objects in `AtomicGraph`.
2. **Read from graph:** If your plugin needs molecular structure, call `backend._sync_sys()` then read from `backend.sys` or `backend._bond_list`.
3. **No duplicate storage:** Do not cache molecular data on `window` or in plugin-local variables beyond the minimal needed for UI display.
4. **Export through main save:** Do not add export buttons to your plugin panel. Bond orders/charges stored on the graph will automatically be used by `save_structure()`.
