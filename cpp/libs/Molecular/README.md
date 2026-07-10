# cpp/libs/Molecular — C-wrapper shared libraries for Python ctypes binding

Each `*_lib.cpp` file here wraps a header-only C++ module from `cpp/common/molecular/`
into `extern "C"` functions callable via Python `ctypes`. The corresponding Python
bindings live in `pyBall/`.

- **FFfit_lib.cpp** — FFfit class handle + standalone graph/topology functions (CSR bond-graph BFS, 1-4 pairs, dihedral enumeration, batch dihedral sensitivity). Wraps `FFfit.h`.
- **MMFF_lib.cpp** — MMFFsp3 force field evaluation, Hessian, neighbor lists. Wraps `MMFFsp3.h`.
- **MMFFsp3_lib.cpp** — Local MMFF variant with OpenCL buffers. Wraps `MMFFsp3_loc.h`.
- **Forces_lib.cpp** — Generic force calculation interface. Wraps `ForceField.h`.
- **Kekule_lib.cpp** — Kekulé bond-order optimization. Wraps `Kekule.h`.
- **FitREQ_lib.cpp** — Charge equilibrium parameter fitting. Wraps `FitREQ.h`.
- **eFF_lib.cpp** — Electron force field. Wraps `eFF.h`.
- **FF2D_lib.cpp** — 2D force field. Wraps `FF2D.h`.
- **CLCFGO_lib.cpp** — Crystal lattice CG optimizer. Wraps `CLCFGO.h`.
- **Molecular.cpp** — Molecular system utilities. Wraps `molecular.h`.
- **ReactiveFF.cpp** — Reactive force field. Wraps `FlexibleAtomReactiveFF.h`.
- **RigidMol.cpp** — Rigid body molecular dynamics. Wraps `RigidBodyFF.h`.
- **SchroedingerGreen*_lib.cpp** — 1D/2D Green's function solvers.
- **utils.cpp** — Shared utility functions.
- **CMakeLists.txt** — Builds all libs as shared `.so` files for ctypes loading.
