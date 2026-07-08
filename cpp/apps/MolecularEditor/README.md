# `cpp/apps/MolecularEditor` — SDL molecular browser & editor apps

Lightweight OpenGL viewers and editors built on `SDL2OGL` + `MMFFparams`.

- **MolecularBrowser.cpp** — ACDsee-style file browser: BROWSE (thumbnail grid) ↔ VIEW (3D). Loads `.xyz`, `.mol`, `.npz`, `.npy`; `--verify-npz` headless parse. See [`doc/Topics/FTIR_Nanocrystals/CPP_MolecularBrowser_NPZ.md`](../../../doc/Topics/FTIR_Nanocrystals/CPP_MolecularBrowser_NPZ.md).
- **BrowserView.h** — Directory grid, folder tiles, scroll preservation, NPZ filtering (`*hessian*`, `*spectrum*` hidden from grid).
- **MolView.h** — Single-molecule 3D view: atoms/bonds, topology AABB wireframes, MMFF bond-length color map + per-type legend.
- **MolGUI** (if present in tree) — Full MMFF editor; heavier than MolecularBrowser.

**Build:** `cpp/Build/apps/MolecularEditor/MolecularBrowser` (CMake target `MolecularBrowser`, links `z` for NPZ).

**Run:** `./tests/tSiNCs/run_cpp_mol_browser.sh` (defaults to `tests/tSiNCs/fixtures/`).
