# `cpp/common/io` — NumPy NPZ/NPY readers for C++ viewers

Header-only parsers for nanocrystal **viewer stages** (`01`/`02`/`03` NPZ) and optional standalone `.npy`. Deflate parity with [`web/common_js/npzIO.js`](../../../web/common_js/npzIO.js).

**Schema contract (v1.2):** [`doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md`](../../../doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md)  
**Pipeline stages + consumers:** [`doc/Topics/FTIR_Nanocrystals/Nanocrystal_NPZ_Pipeline.guide.md`](../../../doc/Topics/FTIR_Nanocrystals/Nanocrystal_NPZ_Pipeline.guide.md)

- **NpyIO.h** — `.npy` v1 header parse, optional mmap; decodes `<f8`, `<f4`, `<i8`, `<i4`, `|u1`, `|b1`; skips unknown dtypes (unicode metadata in `04`/`05` analysis NPZs).
- **NpzIO.h** — `.npz` zip local-header walk, raw deflate; loads decodable arrays, warns on skip.
- **CrystalNpz.h** — `01_crystal` / `02_relaxed`: `pos`, `Z`, `bonds_ij` → `Molecule`; `pos.npy` + `Z.npy` merge path.
- **TopologyNpz.h** — `03_topology`: `group_bbox_min/max`, `icolGroup`; geometry when `pos`/`Z` present. MMFFL `neigh_idx`/`bond_k` not yet drawn as bond sticks in C++ (distance fallback).

**Viewer listing:** MolecularBrowser shows `01_*`, `02_*`, `03_*` only — hides `04_hessian`, `05_spectrum` (`BrowserView::filterNonViewerNpzFiles`).
