# `cpp/common/io` — NumPy NPZ/NPY readers for C++ viewers

Header-only parsers for nanocrystal **viewer stages** (`01`/`02`/`03` NPZ) and optional standalone `.npy`. Deflate parity with [`web/common_js/npzIO.js`](../../../web/common_js/npzIO.js).

**Schema contract (v1.2):** [`doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md`](../../../doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md)  
**Pipeline stages + consumers:** [`doc/Topics/FTIR_Nanocrystals/Nanocrystal_NPZ_Pipeline.guide.md`](../../../doc/Topics/FTIR_Nanocrystals/Nanocrystal_NPZ_Pipeline.guide.md)

- **NpyIO.h** — `.npy` v1 header parse, optional mmap; decodes `<f8`, `<f4`, `<i8`, `<i4`, `|u1`, `|b1`; skips unknown dtypes (unicode metadata in `04`/`05` analysis NPZs).
- **NpzIO.h** — `.npz` zip local-header walk (ZIP64 extra 0x0001; NumPy `np.savez` always uses `force_zip64=True`), stored or raw deflate; loads decodable arrays, warns on skip. `npz_try_read_file` returns false with a loud error (browser skips the tile); `npz_read_file` still `exit(1)` for `--verify-npz`.
- **CrystalNpz.h** — `01_crystal` / `02_relaxed`: `pos`, `Z`, `bonds_ij` → `Molecule`; `pos.npy` + `Z.npy` merge path.
- **TopologyNpz.h** — `03_topology`: `group_bbox_min/max`, `icolGroup`; geometry when `pos`/`Z` present. MMFFL `neigh_idx`/`bond_k` not yet drawn as bond sticks in C++ (distance fallback).

**Viewer listing:** MolecularBrowser shows `01_*`, `02_*`, `03_*` only — hides `04_hessian`, `05_spectrum` (`BrowserView::filterNonViewerNpzFiles`).
