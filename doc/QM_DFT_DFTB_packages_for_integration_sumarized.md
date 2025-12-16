# QM / DFT / DFTB / Semiempirical codes for Python-driven QM/MM integration

This document summarizes LLM-collected notes about **open-source (and a few reference closed-source) electronic-structure codes** that may be integrated into our **Pythonic QM/MM package and scripting workflows**.

Focus:
- **Fast / lower-accuracy solvers first** (DFTB / semiempirical), then DFT, then higher-level reference methods.
- **Integration & hackability** (Python APIs, ability to extract density matrix / Hamiltonian / density on grid or points).
- **Acceleration features** (GPU, linear-scaling / sparse / density-matrix methods, periodic support).

## Selection criteria emphasized in the notes

- **Integration / hackability (Python-first)**  
  Priority given to tools that expose **density matrix**, **Fock/Hcore**, **MO coefficients**, and allow **density evaluation on arbitrary grids/points** inside Python.

- **Efficiency features**
  - **GPU acceleration** (how “native”, what parts accelerated)
  - **Linear-scaling / fast algorithms** (SP2 purification, XL-BOMD, sparse methods)
  - **Periodic support** (PBC, k-points) where mentioned

- **Basis / representation variety** (as stated)
  - **Gaussian** (PySCF, Psi4, NWChem, GATOR, ChronusQ-GPU)
  - **Numerical atomic orbitals / LCAO** (GPAW LCAO; DFTB+ STO + Slater–Koster tables)
  - **Gaussian + numerical / GPW** (CP2K)
  - **Multi-wavelets** (MADNESS-DFT)
  - **Orbital-free DFT** (DFTpy; grid-based, no orbitals)


## 1) DFTB / Semiempirical (fast)

| Code | Methods | Basis / representation (as stated) | GPU | Linear-scaling / fast techniques | Python integration / hackability | Density access (as stated) |
|---|---|---|---|---|---|---|
| [PYSEQM](https://github.com/lanl/PYSEQM) | MNDO, AM1, PM3 (NDDO) | PyTorch tensor implementation (basis not detailed here) | Yes (native PyTorch/CUDA) | XL-BOMD; SP2 density-matrix purification | High (pure PyTorch tensors, hack SCF) | `get_density()` tensor; evaluate on points/grids via tensor ops |
| [DFTB+](https://github.com/dftbplus/dftbplus) | DFTB2/3 etc. | STO basis + Slater–Koster integral tables | **External eigensolver GPU only** (ELPA/MAGMA/cuSolverMP via ELSI is claimed); no CUDA kernels inside DFTB+ itself | linear-scaling SP2 solver CPU-only (claimed) | Medium-High (Python bindings via ASE/sockets claimed; internals Fortran-heavy) | charges easy; density grids/points mainly via outputs/postprocessing |
| [DFTpy](https://github.com/pyscf/dftpy) | orbital-free DFT (OFDFT) | real-space grid / FFT (orbital-less) | **No native CUDA kernels**; speedups depend on backend dispatch (PyTorch/CuPy mentioned); GitHub is claimed to be a mirror, upstream on GitLab | O(N log N) via grids/FFT (claimed) | High (pure Python; integrates with PySCF/ASE claimed) | density on grids native; eval at points trivial (claimed) |
| [SCINE Sparrow](https://github.com/qcscine/sparrow) | NDDO + DFTB variants | not specified | No | no explicit linear scaling claimed | High (Python API); C++ core | density matrix available; density-on-grid not shipped (needs contraction) |
| [PySCF semiempirical](https://github.com/pyscf/semiempirical) | MINDO/3 | Gaussian toolchain via PySCF | Partial mention (via GPU4PySCF for PySCF parts; semi-emp limited) | no explicit linear scaling claimed | Top-tier (same PySCF API) | DM via `make_rdm1()`, density via `numint.eval_rho` |
| PyQuante / PyQuante2 | MINDO/3, MNDO, AM1, PM3 | minimal AO (implied by NDDO); pure Python emphasized | No | linear scaling claimed for large molecules (in source table) | High (pure Python, hackable SCF) | DM accessible post-SCF (claimed) |
| [tblite](https://github.com/tblite/tblite) | GFN0/1/2-xTB etc. | TB-based | CPU-only (OpenMP); no GPU (claimed) | linear-scaling SP2 solver included (claimed) | High (Python bindings via `pip install tblite`) | density matrix exposed via API (`calc.get("density-matrix")` for tblite ≥ 0.4.0 is claimed); density-on-grid not shipped (reconstruct from AO values) |
| MOPAC + pyMOPAC | PM6/PM7/PM3 etc. | NDDO family (implied) | papers about GPU acceleration linked | not detailed | Low for internals (wraps binary) | no DM inside Python; reconstruct from output is clumsy (claimed) |
| [MDTorch](https://github.com/bytedance/MDTorch) | semiempirical QM on GPUs | PyTorch implementation | GPU (implied by title/description) | not detailed | PyTorch-based | not detailed |

 Per-package notes / links (compact):
 - **PYSEQM**
   - Code: https://github.com/lanl/PYSEQM
   - Install: `pip install git+https://github.com/lanl/PYSEQM` (claimed)
   - Notes: pure PyTorch tensors; can backprop through SCF (claimed). Limitations mentioned: no PBC yet; only neutral molecules.
 - **DFTB+**
   - Paper: https://aip.scitation.org/doi/10.1063/1.5143190
   - Code: https://github.com/dftbplus/dftbplus
   - Notes: periodic/k-points/forces emphasized in the source; density on grid via external tooling is mentioned (e.g. *slateratom* utility or *sisl* in one claim-set).
 - **DFTpy**
   - Code (GitHub): https://github.com/pyscf/dftpy
   - Code (GitLab, claimed upstream): https://gitlab.com/dftpy/dftpy
   - Paper: https://doi.org/10.1063/5.0097487
   - Notes: GPU is described as backend-dispatch (PyTorch/CuPy mentioned); no native CUDA kernels.
 - **SCINE Sparrow**
   - Code: https://github.com/qcscine/sparrow
   - Paper: https://aip.scitation.org/doi/10.1063/5.0136404
   - Tools: https://github.com/qcscine/readuct (workflow/driver CLI; structure opt, TS search etc. listed in source)
   - Notes: DM access via `calc.results['density_matrix']` is claimed; density-on-grid is not shipped (contract yourself).
 - **PySCF semiempirical**
   - Code: https://github.com/pyscf/semiempirical
   - Notes: MINDO/3 only is emphasized; same PySCF API, so DM via `mf.make_rdm1()` and ρ via `dft.numint.eval_rho`.
 - **PyQuante / PyQuante2**
   - Homepage mentioned: pyquante.sourceforge.net (no explicit URL in source beyond this text)
   - Notes: emphasized as a pure-Python “sandbox” where SCF is short and hackable.
 - **xTB / tblite**
   - Code (tblite upstream): https://github.com/tblite/tblite
   - Paper: https://doi.org/10.1063/5.0048841
   - Notes: DM API example in source uses `calc.set("save-density-matrix", True)` then `calc.get("density-matrix")`; density-on-grid/points requires AO evaluation + contraction.
 - **MOPAC + pyMOPAC**
   - GPU papers: https://pubs.acs.org/doi/10.1021/ct3001798 ; https://pubs.accs.org/doi/10.1021/ct3004645 ; https://pubs.acs.org/doi/10.1021/acs.jctc.0c00243
   - Notes: source emphasizes wrapper-based integration (binary), no DM inside Python; GPU builds exist but not via the packaged route.
 - **MDTorch**
   - Paper: https://pubs.acs.org/doi/10.1021/acs.jctc.1c00443
   - Code: https://github.com/bytedance/MDTorch

## 2) DFT (LDA/GGA/Hybrid/meta-hybrid)

| Code | Basis (as stated) | Hybrids / B3LYP (as stated) | GPU (as stated) | Periodic (as stated) | Python integration / hackability | Density on grid / points (as stated) |
|---|---|---|---|---|---|---|
| [PySCF](https://pyscf.org/user/gpu.html) | Gaussian | Yes (B3LYP etc.) | GPU via GPU4PySCF | PBC DFT/HF with k-points; hybrids are challenging | Excellent (pure Python, modular SCF) | `numint.eval_rho` at points/grids; DM exposed |
| Psi4 | Gaussian | Yes (B3LYP via LibXC) | via BrianQC plugin (2–10x claimed) | not stated | Good (Python API) | DM via `wfn.Da()`, grid eval via low-level API |
| GPAW (LCAO mode) | numerical atomic (LCAO) | Yes (B3LYP via LibXC) | GPU only in PW mode (not LCAO) | implied via ASE, not detailed | High (ASE/Python) | density via `calc.get_density()` etc. |
| CP2K | Gaussian + numerical (GPW) | Yes (B3LYP full) | GPU via DBCSR (CUDA/HIP/OpenCL) | Yes (explicitly highlighted) | Medium (pycp2k; core Fortran) | density via outputs; less direct in Python |
| NWChem | Gaussian | Yes (B3LYP) | some GPU via MAGMA/cuBLAS (claimed) | Yes (claimed) | Medium-Low (PySCF bridge) | density via directives / postprocessing (claimed) |
| Siesta | NAO (numerical atomic orbitals) | hybrids like HSE/PBE0 mentioned; no full B3LYP (claimed) | no native GPU (claimed) | excellent for periodic (claimed) | medium-high via ASE; Fortran internals | outputs .RHO grids; sisl suggested for postprocessing |
| [QUICK](https://github.com/merzlab/QUICK) | HF/DFT (LibXC implied) | hybrids via LibXC (claimed) | multi-GPU CUDA/HIP (claimed) | not stated | not detailed | not detailed |
| [gpu4pyscf](https://github.com/pyscf/gpu4pyscf) | PySCF extension | hybrids yes (implied) | GPU | not stated | high (PySCF ecosystem) | implied |
| ByteQC | HF/DFT (claimed) | hybrids claimed | GPU; multi-GPU strong scaling claimed | not stated | not detailed | not detailed |
| [Quantum ESPRESSO (QE-GPU)](https://github.com/QEF/q-e-gpu) | plane-wave | hybrids including B3LYP are claimed (with EXX support) | CUDA Fortran kernels for FFT/DGEMM/XC/EXX are claimed | periodic; k-points claimed | `qe_tools` Python pkg mentioned for reading density (points via interpolation claimed) | pp.x -> cube; qe_tools -> interpolate points (claimed) |
| [Octopus](https://octopus-code.org/documentation/main/) | real-space | not stated | GPU noted | not stated | not stated | not stated |
| [JDFTx](https://github.com/jdftx/JDFTx) | plane-wave / pseudopotential (claimed) | hybrids incl. HSE/PBE0/B3LYP claimed | **DISPUTED**: Kimi says no GPU (`Performance.html` says no GPU); Grok/DeepSeek claim CUDA | k-points and periodic claimed | Python `jdftx.io` mentioned in one claim-set; other claim-set not specific | `dump ElecDensity` cube; `dump ElecDensityAccum` grid (claimed) |

 Per-package links / notes (compact):
 - **PySCF / GPU4PySCF**
   - https://pyscf.org/user/gpu.html (GPU docs page)
   - https://rowansci.substack.com/p/gpu-accelerated-dft (blog-style overview/benchmark discussion)
   - https://github.com/pyscf/gpu4pyscf (code)
- **QUICK**
  - https://github.com/merzlab/QUICK (code)
  - https://pubs.acs.org/doi/10.1021/acs.jctc.1c00943 (paper)
  - https://mattermodeling.stackexchange.com/a/10033/175 (thread-linked overview with pubs/benchmarks)
- **ByteQC**
  - https://docs.google.com/presentation/d/1bR6sZGq7lL3p1xQ6Kyk8vFPjf7F4L4PvXz5uGvl8X0E (developer slide deck linked in thread)
- **Octopus**
  - https://octopus-code.org/documentation/main/ (documentation)
- **Quantum ESPRESSO (QE)**
  - https://github.com/QEF/q-e (official repo)
  - https://github.com/QEF/q-e-gpu (GPU branch / workflow, claimed merged into develop)
- **JDFTx**
  - https://github.com/jdftx/JDFTx (official repo)
  - https://jdftx.sourceforge.net (docs)
  - https://jdftx.sourceforge.net/Performance.html (performance page; used in the GPU dispute)
  - https://doi.org/10.1016/j.softx.2017.09.004 (paper)

## 3) HF / post-HF (reference / high-accuracy, GPU-hybrid focused)

| Code | Methods (as stated) | GPU (as stated) | Python integration (as stated) | Density access (as stated) |
|---|---|---|---|---|
| [ChronusQ-GPU fork](https://github.com/chronusq/chronusq/tree/gpu_dft) | hybrids; EOM-CCSD; SF-TDDFT | C++/CUDA | Python binding for density is claimed (`cq.SCF.getDensity()`) | DM exposed; no built-in grid evaluator (claimed) |
| [GATOR](https://github.com/uu-chem-students/GATOR) | HF/DFT; B3LYP/PBE0/CAM-B3LYP; TD-DFT(TDA); CIS | CUDA kernels | PyBind11 layer (claimed) | DM and `eval_rho(grid)` claimed |
| MADNESS-DFT (GPU branch) | B3LYP with exact exchange | GPU branch | PyBind interface (claimed) | density on arbitrary grid (claimed) |
| [TeraChem](http://www.petachem.com/products.html) (closed source reference) | GPU DFT + more | native multi-GPU CUDA | not stated | not stated |

 Per-package links / notes (compact):
 - **ChronusQ**
   - https://github.com/chronusq/chronusq/tree/gpu_dft (GPU branch)
   - https://doi.org/10.1021/acs.jctc.1c01108 (paper)
 - **GATOR**
   - https://github.com/uu-chem-students/GATOR (repo mentioned explicitly in the source)
   - https://pubs.acs.org/doi/10.1021/acs.jctc.3c00352 (poster abstract / citation)
   - **Notes**: described as an open-source “TeraChem child” with CUDA kernels and a Python layer exposing tensors needed for embedding (DM/Fock/MO coeffs are claimed elsewhere in the compilation).
 - **MADNESS-DFT**
   - https://pubs.acs.org/doi/10.1145/2833157 (MADNESS environment paper)
   - https://doi.org/10.26434/chemrxiv-2022-k2hhf (GPU acceleration preprint)
 - **TeraChem**
   - http://www.petachem.com/products.html (product page)
   - https://pubs.acs.org/doi/10.1021/ct9003004 (paper)
   - https://pubs.acs.org/doi/10.1021/ct3007046 (paper)

## To investigate next (missing info / tick-the-box gaps)

Consensus vs disputed (from USER8 additions):
- **Consensus**
  - **DFTB+ GPU**: no native GPU kernels; GPU acceleration is via external eigensolver libraries (ELPA/MAGMA/etc.) only.
  - **xTB/tblite**: no GPU; density matrix API exists; density-on-grid/points require reconstruction ("hack").
  - **ChronusQ-GPU**: GPU branch `gpu_dft` exists; density matrix accessible via Python binding; no built-in grid evaluator.
  - **QE-GPU**: GPU coverage includes FFT/XC/EXX in the claimed QE-GPU workflow; density-at-points requires helper tooling (`qe_tools`).
- **Disputed / inconsistent across the appended LLM answers**
  - **DFTpy canonical repo URL**: one reconciliation claims GitLab upstream with GitHub mirror; the exact GitLab URL differs across sources.
  - **JDFTx GPU**: Kimi claims CPU-only (and links `Performance.html`); Grok/DeepSeek claim native CUDA.

Uniform checklist (copied from the compilation; note that JDFTx GPU remains disputed elsewhere in the same file):

| Code | DM | H/Fock | ρ(grid) | ρ(points) | forces | stress | k-points | GPU cov. |
|---|---|---|---|---|---|---|---|---|
| DFTpy | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | FFT/XC |
| JDFTx | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | none *(DISPUTED in other appended text)* |
| QE-GPU | ✓ | ✓ | ✓ | ✓* | ✓ | ✓ | ✓ | FFT/XC/EXX |
| xTB/tblite | ✓ | — | hack | hack | ✓ | — | ✓ | none |
| ChronusQ | ✓ | ✓ | — | — | ✓ | — | Γ only | XC/integrals |
| DFTB+ | — | — | — | — | ✓ | ✓ | ✓ | diag only |
| PySCF-GPU | ✓ | ✓ | ✓ | ✓ | ✓ | — | Γ/k | XC/grad |
| GATOR | ✓ | ✓ | ✓ | ✓ | ✓ | — | Γ only | all |

*through post-processing helper (`qe_tools`)