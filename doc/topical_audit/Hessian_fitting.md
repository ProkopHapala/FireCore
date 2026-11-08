---
type: TopicalAudit
title: Hessian / vibration force-field fitting
description: Inventory of bonded-FF stiffness fitting to QM Hessians (Wilson/hybrid FFfit, own-min k-scans, planned restricted-subspace work). Not FitREQ, not GridFF.
tags: [fffit, hessian, wilson, vibration, mmff, nanocrystal]
timestamp: 2026-09-03
---

# Hessian / vibration force-field fitting

**This file is the inventory + method checklist + idea catalog.** Theory SSOT: [`../Topics/FFfit/HessianFitting_Theory.md`](../Topics/FFfit/HessianFitting_Theory.md). Numbers: [`../../tests/tSiNCs/SiNCs_FFfit_summary.md`](../../tests/tSiNCs/SiNCs_FFfit_summary.md). Own-min C L1: [`../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md`](../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md). GPT 5.6 design notes (2026-09-03): [`../Topics/FTIR_Nanocrystals/Hessian_fitting.chat.md`](../Topics/FTIR_Nanocrystals/Hessian_fitting.chat.md) from line ~1005.

## Summary

At **fixed geometry**, a bonded FF \(E=\sum_p k_p f_p(q_p)\) has Hessian **linear in \(k\)**: \(H=\sum_p k_p A_p\). That is ordinary (bounded, regularized) least squares. Random annealing is the wrong primary solver for this problem. The current production fitter (`fit_hybrid_hessian`) already does that linear solve with three objectives: DFT **modes**, **graph-local** Cartesian blocks, and **Wilson row-space**.

**Do not abandon Wilson.** The old mistake was treating the least-norm \(F=C^{+T}DC^+\) diagonal as per-bond DFT springs (central Si–Si looked soft: gauge, not physics). The current row-space objective is mathematically correct. GPT 5.6: once bonds+angles span \(3N-6\), that “internal” term is almost the **full vibrational Hessian**, not a local bond/angle fit. For FTIR Si–H (~2000 cm⁻¹) restrict \(B\) to chemically chosen subspaces (Si–H stretches first).

Two protocols that must not be mixed:

| Protocol | Geometry | Valid product | Code |
|----------|----------|---------------|------|
| **FFfit** | Shared QM/DFTB min | \(k\) matching \(H\) at that \(q\) | `test_FFfit.py` |
| **Own-min spectrum** | Each FF’s FIRE min | Frequencies / DOS | `fit_mmff_kss_pyscf.py` |

Spectrum rule: [`Hessian_at_own_minimum.md`](../Topics/FTIR_Nanocrystals/Hessian_at_own_minimum.md).

**Not this topic:** FitREQ (H-bond REQ), GridFF B-spline fit, FDBM. `pyBall/FFFit.py` (capital F) is Fireball linear-scan, not Hessian FFfit.

**JS / OpenCL:** no Hessian fitter. JS is visualization only (`FFfit_plots.build_stiffness_html` p5.js). OpenCL is GridFF/MD, not \(A_p\) assembly.

**HARD (2026-09-03, GPT 5.6 + user):** **Do not rigid-mode-project \(H\) for fitting.** \(P=I-RR^T\) (and Löwdin \((CC^T)^{-1/2}\)) are dense and smear local Wilson directions. Restricted Wilson already has \(B R_{\mathrm{rigid}}=0\). Spectrum: diagonalize raw \(D=M^{-1/2}HM^{-1/2}\), **identify** the six rigid eigenvectors by overlap with analytic \(T,R\), ignore them — do not delete/project before storing or fitting. Response \(Hx=f\) may use a *matrix-free* \(v\leftarrow v-R(R^T v)\), never a stored dense \(P\). Current `reference_vibrational_modes` does \(Q^T D Q\) — analysis only, **not** the new local-fit path. Target is **qualitative/relative** (SiH vs SiH₂ order, in-phase vs out-of-phase, facet shifts), not 1% cm⁻¹. Prefer coupling ratios \(\eta_{ij}=K_{ij}/\sqrt{K_{ii}K_{jj}}\) over global RMSD.

Systematic protocol: [§ Protocol](#systematic-protocol-small--large).

---

## Method checklist

Marks: `[x]` done and usable · `[*]` tried, rejected or diagnostic-only · `[~]` exists, needs retune · `[ ]` planned, not coded.

### A. What we compare (objectives)

- [x] **Linear \(H=\sum k_p A_p\)** at frozen \(q\). Core of `FFfit.h` / `fit_hybrid_hessian`. Condition after scaling ~1.3–1.7. **Use this**, not annealing, for stiffness.
- [x] **Mode Rayleigh** \(v_s^T D(k) v_s=\lambda_s^{\mathrm{DFT}}\). Default `mode_balance=frequency`, `mode_mixing=0.1`. Good for spectra; Hungarian not needed inside the LS solve.
- [x] **Graph-local Cartesian** \(d_{\mathrm{graph}}\le 2\) (`--local-graph-distance`). Does **not** punish FF for missing long-range QM blocks. \(d=3\) available for 1–4 tests.
- [x] **Gauge-invariant Wilson row-space** \(\|Q^T(D^{\mathrm{FF}}-D^{\mathrm{DFT}})Q\|_F\). Correct replacement for fitting \(F_{ii}\). **Caveat (GPT 5.6):** full bond+angle \(B\) often has \(\mathrm{rank}=3N-6\), so this is nearly a **global** vibrational Hessian, redundant with all-mode Frobenius. Default CLI still weights mode/local/internal all 1.0 — too opaque for a surface-FTIR fit.
- [*] **Least-norm Wilson \(F=C^{+T}DC^+\)** (`internal_hessian_projection`, `--plot-stiffness`). Diagnostic only. Individual \(F_{ii}\) are **not** unique springs. Central-bond softness was this artifact. [`HessianFit.chat.md`](../Topics/FFfit/HessianFit.chat.md).
- [*] **Full Cartesian Frobenius** \(\|H^{\mathrm{FF}}-H^{\mathrm{DFT}}\|_F\). Punishes physics the local FF cannot represent. Graph-local + projected \(K_Q\) are the replacements.
- [*] **Hungarian frequency RMSE as the optimizer loss.** Used for reporting and for own-min anneal. Nonlinear, mode-crossing. Fine as a metric; not the frozen-geom solver.
- [x] **Restricted Wilson hydride \(K_{ab}\)** — `wilson_kab_design` / `--ladder si`: raw \(g_a\), no \(C^+\), no Löwdin. Default `wilson_rows='hydride'` = XH + any angle with ≥1 H. Stretch-only via `wilson_rows='bonds'`.
- [ ] **Restricted \(B_{\mathrm{surf}}\)** vs \(B_{\mathrm{core}}\) as *separate weighted blocks* (current hydride rows mix HXH+XXH).
- [ ] **Si–H window mode fit** \(u_m^T H^{\mathrm{FF}} u_n\) in the **unprojected** DFT eigenbasis (modes identified after eigh of raw \(D\); rigid ones dropped by overlap, not by \(PDP\)). Still linear.
- [ ] **Directional Hessian** \(Hg_a\) from \(\pm\delta g_a\) DFT forces (\(2m\) gradients, not \(6N\)). Needed for large NCs. Same object as local Wilson response.
- [x] **No Löwdin / no \(PDP\) on the fitting Hessian** (ladder path: `mode_weight=0`, `internal_weight=0`). Classify rigid modes after eigh. [`Hessian_fitting.chat.md`](../Topics/FTIR_Nanocrystals/Hessian_fitting.chat.md) ~2893.

### B. Parameterizations (what \(k\) we vary)

- [x] **Elemental 5-param Si** (\(k_{\mathrm{SiSi}},k_{\mathrm{SiH}}\), three angles). Mean freq RMSE **41 cm⁻¹** on SiH₄ + Si spheres. Minimal honest basis.
- [x] **Hierarchical Si / SiH / SiH₂ / SiH₃** (`--si-subtypes`, `--subtype-shrinkage 0.001`). RMSE **32 cm⁻¹**. Main spectral win. Fitted Si–Si \(k\sim 9.2\)–\(9.8\) eV/Å² (not the Wilson-diagonal fake softness).
- [x] **Stretch–stretch \(K_{rr}\)** (`--stretch-stretch`). Almost negligible on Si network. `[*]` as a Si-backbone extra.
- [x] **Stretch–bend \(K_{r\theta}\)** (`--stretch-bend`). Hessian relFrob 6.9% → **5.9%**; spectrum barely better. Largest useful extras were **H–Si–H**. Requires `--equilibrium local` (type-average prestress Hessian not implemented).
- [*] **UFF/Prokop torsions** (`--dihedrals`). Implemented + C++/Python parity. **Not production** for tetrahedral Si/C: indefinite prestress, correlated with angles. [`Dihedral_Torsion.md`](../Topics/FFfit/Dihedral_Torsion.md).
- [*] **1–4 distance springs** (`--third-bonds`). Exists. GPT 5.6: postpone; helps cage/soft modes more than the ~2000 cm⁻¹ manifold. Prefer explicit H···H or \(J q_i q_j\).
- [~] **Type-average \(r_0,\theta_0\)** (`--equilibrium type-average`). Transferable geometry; prestress \(f'C\) kept. Cross terms **disabled** on this path.
- [ ] **SiH–SiH stretch coupling \(J q_i q_j\)** — same-Si (SiH₂) and neighboring-Si. Direct control of \(\lambda_\pm=k\pm J\). Not in the hybrid CLI yet. GPT 5.6 preferred extra over generic 1–4.
- [ ] **H···H Pauli** \(A e^{-br}\) (fix \(b\), fit \(A\)). **Stage 2** after harmonic \(k\): nonzero force at DFT \(q\), not linear in the same way. Own-min C L1 already showed surface NB helps cubes and hurts octahedra ([`MMFF_C_CH_vs_CH2_kfit.md`](../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md)) — that was nonlinear FIRE, not projected \(K_H\).
- [ ] **Typed CH vs CH₂ \(k\)** (C FTIR). Shared \(k_{\mathrm{XH}}\) failed (~72 cm⁻¹ SA). FFfit hierarchy analog of Si subtypes, on **stationary** C Hessians.

### C. Solvers

- [x] **`scipy.optimize.lsq_linear`** + column scaling + Tikhonov + \(k\ge 0\) (`solve_regularized_lsq`). Production hybrid solver.
- [*] **Normal equations** `G k = y` (`FFfit.solve_normal_equations`). Squares the condition number. Kept in C++; hybrid path does **not** use it.
- [*] **Gradient descent** on Hessian residual (`fit_gradient_descent_multi`). Legacy; linear LS supersedes it for \(k\).
- [*] **Simulated annealing / random walk on \(k\)** at frozen \(q\). Obsolete for harmonic stiffness. Still used in `fit_mmff_kss_pyscf.py --ref anneal` because Morse α + FIRE geometry make the **own-min spectrum** nonlinear. That SA did **not** unify cube/octa C.
- [x] **SVD / condition / unobservable-column checks** in hybrid diagnostics. Use before adding a new term (GPT 5.6: identifiability before accepting \(H_p\)).

### D. Locality / DFT-Hessian diagnostics (GPT 5.6 — do before more FF terms)

A local spring network already has **global** phonons (\(H\) sparse, \(H^{-1}\) dense). DFT \(H\) is the Born–Oppenheimer (electrons eliminated) Hessian, so it can look dense even when chemically local. Wide-gap C should be nearsighted; Si less so; surface X–H has weak dipoles that can still split a near-degenerate stretch band (\(J/k\sim 0.01\) → ~10 cm⁻¹ at 2000 cm⁻¹).

- [ ] **Plot \(\|H_{ij}\|_F\) vs graph distance on the RAW stored Hessian** (before any \(P\)). First check **acoustic sum rule** \(\sum_j H_{ij}\approx 0\) and **six near-null eigenvalues**. If far blocks sit on a \(1/N\) floor, the file was already projected — do not use it for locality. PySCF `mf.Hessian().kernel()` and DFTB+ `hessian.out` are saved unprojected; `frequencies_cm1.npy` / `modes.npy` may come from `thermo.harmonic_analysis` (analysis only). **Never reconstruct \(H\) from projected eigenpairs.**
- [ ] **Truncate DFT \(H\) at \(d=1,2,3,\ldots\)**. For the *spectrum* of the truncated matrix, restoring ASR is a **diagnostic**, not a fit-time densify. Report full RMSE **and** hydride-window RMSE / \(\eta_{ij}\) signs.
- [ ] **Build \(K_H=g_a^T H g_b\)** (raw \(g\), no \(C^+\)) vs H···H, same-Si vs neighbor, facet. Also \(\eta_{ij}\). This decides \(J\) vs \(Ae^{-br}\) vs \(K_{r\theta}\).
- [ ] Same \(K_H,\eta\) from current FF vs DFT residual.

### E. Own-min / nonlinear (separate protocol)

- [x] Frozen-geom \(k_{\mathrm{XH}}\) grids on cages (`fit_mmff_kss_pyscf.py`, CH₄/adamantane/…). [`MMFF_stiffness_scaling.md`](../Topics/FTIR_Nanocrystals/MMFF_stiffness_scaling.md).
- [x] Own-min stretch RMSE vs PBE L1 (`--ref joint`). C \(k_{\mathrm{XH}}=1.775\), mean 83 cm⁻¹. Si numbers in `OUT_mmff_kss_fit/pyscf_pbe_ncs/*_Si` were vs **off-min** Batch-1 refs — do not cite as FTIR; re-run vs `pyscf_vib_results` (5 Si now 0-imag).
- [*] Surface LJ H-EvdW scan; Morse α × \(k_{\mathrm{XH}}\); SA on \(k_{\mathrm{XH}},\alpha,K_{\mathrm{XXX}},K_{\mathrm{XXH}},K_{\mathrm{HXH}}\). Families see-saw. **Closed (fail)** as a shared C pack.
- [ ] Apply **fitted FFfit \(k\)** (hierarchy + optional \(J\)) then FIRE+Hessian on motif atlas. That is the FTIR product; hybrid LS is the \(k\) factory.

---

## Idea catalog (second pass)

Four goals, often in tension: **Acc**uracy of the spectrum we care about · **Rob**ustness (identifiable \(k\), not overfit, not fragile to geometry/gauge) · **Int**erpretability (one knob ↔ one physical effect) · **Cheap** DFT (fewer than \(6N\) force evals, smaller molecules first).

The current hybrid LS already scores well on **Rob** for frozen-\(q\) stiffness. The remaining FTIR gap is mostly **Acc** on the hydride *manifold* (splitting/dispersion), not mean stretch, plus **Cheap** for Wulff-sized NCs. Adding more generic UFF terms usually hurts **Rob** and **Int** without helping the 2000 cm⁻¹ window.

### Master table

Marks in *Status* match the checklist. Goal columns: `+` helps, `0` neutral, `−` hurts or risks.

| Idea | Acc | Rob | Int | Cheap | Status | Note |
|------|:---:|:---:|:---:|:-----:|--------|------|
| Linear \(H=\sum k A_p\) + bounded LS | + | + | + | + | [x] | No anneal, no \(\partial\omega/\partial k\). SVD tells identifiability. |
| Don’t fit least-norm \(F_{ii}\) as springs | 0 | + | + | 0 | [*] | Gauge; central Si–Si softness. Keep as diagnostic plot only. |
| Don’t fit full dense \(\\|H\\|_F\) | + | + | + | 0 | [*] | Punishes physics a local FF cannot have. |
| Graph-local Cartesian \(d\le 2\) | + | + | 0 | 0 | [x] | Matches FF interaction range. |
| Full Wilson row-space (all bonds+angles) | ~ | − | − | 0 | [~] | If \(\mathrm{rank}(B)=3N-6\), ≈ global Hessian; redundant with all-mode term. Default weight 1.0 too high for FTIR. |
| Restricted \(B_H\) (XH stretches only) | + | + | + | 0 | [ ] | **Main next Acc+Int lever.** \(K_H=G^T D G\) is the hydride Hamiltonian. |
| \(B_{\mathrm{surf}}\) vs \(B_{\mathrm{core}}\) split weights | + | + | + | 0 | [ ] | Don’t let cage modes dominate the Si–H fit. |
| All-mode Rayleigh, frequency-weighted | + | 0 | 0 | 0 | [x] | Good global spectra; swamps hydride fine structure if weight=1. |
| Windowed modes \(U\) in 1800–2200 (C–H analog) | + | + | + | 0 | [ ] | Still linear; off-diagonals = wrong mixing. FTIR-shaped loss. |
| Hungarian \(\omega\)-RMSE as *optimizer* | − | − | − | 0 | [*] | Nonlinear, crossings. Metric only. |
| Elemental 5-param Si | 0 | + | + | 0 | [x] | Honest baseline (41 cm⁻¹). Start here, add one term. |
| Hierarchy Si/SiH/SiH₂/SiH₃ + shrinkage | + | + | + | 0 | [x] | Best Acc win so far (→32). Shrinkage stops sparse bulk from overfit. |
| Freeze \(k_{\mathrm{XX}}=1\), fit only \(k_{\mathrm{XH}}\) | − | + | + | 0 | [x] | Own-min C; too few knobs for CH vs CH₂. |
| Stretch–stretch \(K_{rr}\) on Si | 0 | − | 0 | 0 | [*] | Extra params, almost no spectral gain. |
| Stretch–bend \(K_{r\theta}\) | ~ | − | + | 0 | [x] | Helps Hessian more than FTIR; HSiH is the useful one. Local-\(q\) only. |
| UFF torsions / 1–4 springs | ~ | − | − | 0 | [*] | Cage/soft modes; indefinite / correlated. Not the 2000 cm⁻¹ term. |
| \(J q_i q_j\) same-Si and neigh-Si | + | + | + | 0 | [ ] | Direct splitting \(\lambda_\pm=k\pm J\). Two numbers, SVD-checkable. |
| H···H \(Ae^{-br}\) (fix \(b\), fit \(A\)) | + | 0 | + | 0 | [ ] | Geometry-dependent coupling. **Stage 2** (forces ≠ 0 at DFT \(q\)). |
| Split \(K_{\mathrm{XXX}},K_{\mathrm{XXH}},K_{\mathrm{HXH}}\) | ~ | − | + | 0 | [*] | Own-min SA: right extra, still see-saw CH vs CH₂. |
| Typed CH vs CH₂ \(k\) | + | + | + | 0 | [ ] | C analog of Si subtypes. Shared \(k_{\mathrm{XH}}\) failed. |
| Type-average \(r_0,\theta_0\) + prestress | 0 | + | + | 0 | [~] | Transferable; needed for a real FF. Cross terms not on this path yet. |
| Local \(r_0=\!r^{\mathrm{DFT}}\) (zero prestress) | + | + | + | 0 | [x] | Cleanest linear Hessian fit. Not a transferable FF by itself. |
| Tikhonov toward UFF / previous \(k\) | 0 | + | 0 | 0 | [x] | `--reg`. Merge or freeze if \(\sigma_{\min}\) small. |
| SVD / drop unobservable columns | 0 | + | + | 0 | [x] | Do **before** accepting a new \(H_p\). |
| Multi-system joint fit (SiH₄ + spheres) | + | + | 0 | 0 | [x] | Shares \(k\); held-out crystal still missing as a routine check. |
| Leave-one-crystal-out / min-type-count | 0 | + | 0 | 0 | [~] | `--min-type-count` flags rare types; no automated LOO. |
| Two-stage: harmonic \(k\) then steric \(A\) | + | + | + | 0 | [ ] | Don’t let H···H wreck a well-conditioned valence LS. |
| \(K_H\) residual vs \(r_{\mathrm{HH}}\), same-Si, facet | + | 0 | + | + | [ ] | Decides \(J\) vs exponential vs \(K_{r\theta}\) from data, not guessing. |
| Truncate DFT \(H\) at \(d\), repair ASR, re-eigh | 0 | + | + | + | [ ] | Uses *existing* Hessian. If \(H^{(2)}\) already matches Si–H, stop adding range. |
| Plot \(\\|H_{ij}\\|_F\) vs graph distance | 0 | + | + | + | [ ] | “Dense” vs “tiny tails”. Cheap on Si₁₀H₁₆. |
| Directional \(Hg_a\) (\(2m\) DFT forces) | 0 | 0 | + | + | [ ] | **Cheap path for large NCs.** \(g_a^T H g_b\) for all \(b\) from one pair. Not the same as dropping distant \(H_{ij}\) by hand. |
| Cartesian NN blocks only (skip far \(H_{ij}\)) | − | 0 | 0 | + | [*] | GPT: worse than \(Hg_a\). Misses “stretch this bond, force on whom?”. |
| Small motifs first (XH₄, X₂H₆, adamantane) | 0 | + | + | + | [x] | Conditioning lab. Transfer, then check NC. |
| Full \(6N\) Hessian on every Wulff L1 | + | 0 | 0 | − | [~] | Have it for 12 L1 jobs; **10/12 now 0-imag** (`pyscf_vib_results`). Two saddles pending mode-follow. |
| Own-min FIRE+Hungarian after LS | + | 0 | 0 | 0 | [ ] | FTIR product. LS is the \(k\) factory. |
| Anneal Morse+angles on own-min stretch RMSE | − | − | − | 0 | [*] | Nonlinear, see-saw. Not the harmonic fitter. |
| Don’t rigid-project \(H\) / don’t Löwdin \(g_a\) for the fit | + | + | + | 0 | [x] | HARD. Ladder path skips `reference_vibrational_modes`. |
| Fit \(Hg_a\) or \(K_{ab}=g_a^THg_b\) (raw local \(g\)) | + | + | + | + | [x] | `wilson_kab_design` + graph-local Cartesian. \(Hg_a\)-only still open. |
| Coupling ratio \(\eta_{ij}=K_{ij}/\sqrt{K_{ii}K_{jj}}\) | + | + | + | 0 | [ ] | Sign = in-phase vs out-of-phase order. Qualitative FTIR target. |
| One global scale \(s\) + fit \(\eta\) | + | + | + | 0 | [ ] | Tolerate 5–10% overall \(\omega\) error if order/splitting is right. |
| Hierarchical lock: XH₄ → cage → NC | + | + | + | + | [x] | `--ladder si` freezes A hydride \(k\). Si 2026-09-03: \(k_{\mathrm{XX}}\) D vs E within 2%; \(k_{\mathrm{XXX}}\) see-saw 19%. |
| Equal weight mode=local=internal=1 | − | − | − | 0 | [~] | Opaque. Ladder path uses local=1, \(K_{ab}=1\), mode=internal=0. |
| Cosθ angle form, \(k_\theta\) in eV/rad² | 0 | + | + | 0 | [x] | Avoid \(1/\sin\theta\). |
| \(k\ge 0\) bounds; signed \(J\) unbounded-but-capped | 0 | + | + | 0 | [x] | Negative springs = unstable Hessian. |

### How the goals pull in different directions

**Accuracy of mean X–H stretch** is already mostly a single \(k_{\mathrm{XH}}\) (linear, cheap). Hierarchy and windowed \(B_H\) are for **Acc of the band shape** (SiH vs SiH₂, {100} vs {111}).

**Robustness** is “can the data uniquely determine this extra \(k\)?” More types, torsions, and 1–4 make \(A\) columns nearly parallel. Hierarchy shrinkage, SVD, one-extra-at-a-time, and *not* fitting gauge-dependent \(F_{ii}\) are the robustness tools we already trust. Own-min annealing failed robustness: two crystals wanted opposite corners of the same 5-D box.

**Interpretability** wants a short list of chemically named coordinates (this Si–H, this H–Si–H, this H···H pair), not a 3N−6 SVD soup. Restricted Wilson + \(K_H\) plots are that. Full Wilson row-space is the opposite, even though it is mathematically clean.

**Cheap DFT:** a full analytic Hessian is \(6N\) displacements (or one coupled-perturbed job). For \(m\) hydrides, \(2m\) directional forces give the whole \(K_H\). Locality truncations and \(K_H\) scatter plots reuse Hessians we already have (SiH₄, Si₁₀H₁₆, spheres) — do those *before* new PySCF. Do not pay \(6N\) on L2 to hunt a term the truncated-L0 Hessian would have ruled out.

### Decision rule (when tempted to add a term)

1. Show it on a \(K_{ab}\) / \(\eta_{ij}\) or truncated-\(H^{(d)}\) plot (Cheap + Int).
2. Build \(H_p\), inspect \(\sigma_{\min}\) of the stacked \(A\) (Rob). If the new column is parallel to an old one, merge or drop.
3. Fit with the extra on; require a drop in **window** RMSE **or** correct \(\eta\) signs, not only global relFrob. Pinned hydride \(k\) must not jump.
4. Own-min FIRE last (Acc of the actual FTIR protocol).

Skip 1–4 and torsions unless step 1 shows the residual living in the Si backbone / soft cage, not in X–H.

---

## Systematic protocol (small → large)

**Do not start over on each molecule.** SiH₄ (CH₄) owns hydride springs; the crystal is not allowed to forget them. Regularize new \(k\) toward the previous stage; freeze or tightly pin what is already determined.

C and Si are **separate tracks**. PySCF and DFTB+ are **separate references** (do not LS-mix Hessians from different Hamiltonians). DFTB+ is the cheap locality / transfer check; PySCF is the stiffness source when both exist.

### 0. Assemble the reference bank

Inventory what we actually have (do not recompute until this table is filled and locality-checked). Paths are off-repo unless noted.

| Stage | Molecules | PySCF Hessian | DFTB+ Hessian | Notes |
|-------|-----------|---------------|---------------|-------|
| L0 mol | CH₄, SiH₄ | cages / `jobs_pyscf_vib_OUT_small_nc/SiH4`; CompChemUtils `run_vib_spectra` | 3ob / pbc / mio as used in `fit_mmff_kss_pyscf.py` | **Fit \(k_{\mathrm{XH}}\), \(k_{\mathrm{HXH}}\) only** |
| L0 dimer | C₂H₆, Si₂H₆ | same | same | Unlock \(k_{\mathrm{XX}}\), \(k_{\mathrm{XXH}}\); pin hydride \(k\) |
| L0 cage | adamantane C₁₀H₁₆, Si₁₀H₁₆ | `jobs_pyscf_vib_chem_atlas_L0`; small_nc | L1 DFTB tree + cage SK jobs | First XX/XXX in a Td network. Mixed XH/XH₂ on one molecule. |
| L1 Wulff | cube, octa, ±5-ring, ±7-ring (C and Si) | **`/home/prokop/SIMULATIONS/SiNCs/pyscf_vib_results`** (canonical 2026-09-03). Legacy Batch-1: `.../pySCF/jobs/results` (Si off-min). Dirs may be `luna_*` / `tight_*`; loaders resolve via `resolve_pyscf_vib_case_dir`. | `/home/prokop/SIMULATIONS/SiNCs/DFTB/L1` | **10/12 at min.** Cite: five C (not `cube_C`) + five Si (not `octahedron_7ring_Si`). Pending mode-follow: `cube_C` (−82.5), `octahedron_7ring_Si` (−25.8). |
| Spheres (legacy) | Si_R3p8 … R6p0, SiH₄ | `jobs_pyscf_vib_OUT_small_nc` — this is the 2026-07 FFfit report set | — | Keep as regression, not the motif ladder. |
| L2 | octa / trunc / rhombic | — | 3ob C own-min; MMFF-at-DFTB-geom is **quarantined** | FTIR overlays; not the first LS target. |

For each case store/keep: `relaxed.xyz`, **raw** `hessian.npy` (or `hessian.out`), `masses.npy`, `status.json` (`n_imag`). Do **not** use `modes.npy` to rebuild \(H\). PySCF `frequencies_cm1.npy` is **3N−6** (rigid dropped in analysis); `n_{\mathrm{rigid}}\approx 0` in that array is expected. The stored Hessian is \((N,N,3,3)\) Hartree/Bohr² — still the full Cartesian kernel, not \(PDP\).

**L1 PBE bank (2026-09-03, `pyscf_vib_results`):** 10 of 12 are minima. Loaders: `pyBall.FFfit_utils.resolve_pyscf_vib_case_dir` (picks `modefollow` > `tight` > `luna`).

| Case | Job dir | \(n_{\mathrm{imag}}\) | FTIR? | Note |
|------|---------|----------------------:|:-----:|------|
| `cube_5ring_C` … `octahedron_7ring_C` (5) | `luna_*` | 0 | PDOS yes | Same as Batch-1 0-imag C |
| `cube_C` | `tight_cube_C` | 1 (−82.5) | **No** | Td saddle; mode-follow queued |
| `cube_Si`, `cube_5ring_Si`, `cube_7ring_Si`, `octahedron_Si`, `octahedron_5ring_Si` | `tight_*` | 0 | **yes** (own-min) | FMAX 0.01→0.001 |
| `octahedron_7ring_Si` | `tight_octahedron_7ring_Si` | 1 (−25.8) | **No** | Tip-local saddle; mode-follow queued |

Old `pySCF/jobs/results` Si (5–11 imag) are **superseded**. Do not mix Batch-1 Si freqs with tight geometries.

### 1. Sanity: is this Hessian still sparse / unprojected?

On **every** file before fitting:

1. Acoustic sum rule: \(\max_i\|\sum_j H_{ij}\|\) (translation). Near machine/SCF noise.
2. Rotations: \(H R_\alpha\) small **only if** the geometry is a stationary point of that Hamiltonian. If not, do not “fix” by projecting.
3. Six eigenvalues of \(D\): classify by overlap \(w^{\mathrm{rigid}}=\sum_{\alpha}|R_\alpha^T v|^2\). Expect ~6 with \(w\approx 1\) and tiny \(\omega\). Fail loud if a stretch-like mode is mixed into “rigid”.
4. Locality histogram: \(\|H_{ij}\|_F\) vs graph distance **on this raw \(H\)**. A flat far-block floor \(\sim 1/N\) ⇒ likely already \(PDP\)’d. Compare to a known-raw MMFF FD Hessian of the same graph (should be exactly sparse beyond \(d=2\)).

L1 PBE `hessian.npy` is \((N,N,3,3)\) with ASR \(\max|\sum_j H_{ij}|\sim 10^{-3}\) Ha/Bohr² (SCF noise, not a \(1/N\) projector floor). `frequencies_cm1.npy` is already 3N−6 — that is analysis, not evidence that \(H\) was projected.

If (4) looks projected, **stop** — refit/re-export from `Hessian().kernel()` / `hessian.out`, do not proceed.

### 2. Fit ladder (one Hamiltonian at a time)

At each stage: frozen QM geometry, \(r_0,\theta_0\) = that geometry for the linear \(A_p\) (prestress off). Solver: bounded LS + SVD. **Loss:** local \(K_{ab}\) / \(Hg_a\) (and optional windowed \(u^T H u\)). **Not** Hungarian \(\omega\)-RMSE. **Not** full-Wilson \(Q\) with \(\mathrm{rank}=3N-6\). **Not** `internal_weight=1` default.

| Stage | Unlock | Pin / prior | Success = |
|-------|---------|-------------|-----------|
| **A** XH₄ | \(k_{\mathrm{XH}}\), \(k_{\mathrm{HXH}}\) | — | \(K_{\mathrm{XH}}\), \(K_{\mathrm{HXH}}\) match; 4 stretch freqs; \(\eta\) N/A |
| **B** X₂H₆ | \(k_{\mathrm{XX}}\), \(k_{\mathrm{XXH}}\) | Tikhonov hydride \(k\) to A (strong) | XX stretch + XXH; hydride springs barely move |
| **C** X₁₀H₁₆ | \(k_{\mathrm{XXX}}\) refine \(k_{\mathrm{XX}}\) | hydride still pinned to A; XX prior from B | cage fingerprint; XH vs XH₂ **order** if both present |
| **D** L1 octa (XH-rich) | maybe \(J_{\mathrm{neigh}}\), tiny subtype offsets | all of A–C | octa stretch window vs PySCF; \(\eta\) vs \(r_{\mathrm{HH}}\) |
| **E** L1 cube (XH₂-rich) | \(J_{\mathrm{sameSi}}\), typed XH vs XH₂ if SVD allows | same | cube window; **must not** destroy D |
| **F** 5/7-ring | ring-tagged extras only if residual lives there | all previous | ΔS of ring vs parent; no see-saw |
| **G** own-min FIRE | none (use locked pack) | — | FTIR protocol on C L1 0-imag; Si only after re-relax |

If a later stage wants to move a pinned \(k\) by more than a set relative tolerance (e.g. 5%), **fail loud** and inspect — that is overfitting, not transfer.

`--subtype-shrinkage` and `--reg` already implement the math for D–F; the missing piece is **explicit freeze/pin groups** (XH₄-determined vs crystal-determined) in the CLI.

### 3. What to compute at each stage (same plots, growing molecules)

After (1) passes:

1. \(K_{aa}\) for each chemically typed local \(g_a\) (XH, XX, HXH, XXH, XXX).
2. Nearby \(K_{ab}\) and \(\eta_{ij}\) (same-Si hydrides, vicinal, 1–4). Scatter vs \(r_{\mathrm{HH}}\).
3. \(H^{\mathrm{FF}}g_a\) vs \(H^{\mathrm{DFT}}g_a\) with optional \(d\le 2\) atom mask (local → local).
4. DFT eigenmodes of **unprojected** \(D\): use hydride-character modes as extra linear targets \(u^T H^{\mathrm{FF}} u\) (user point 4). Off-diagonal \(u_m^T H^{\mathrm{FF}} u_n\) penalizes wrong mixing. Rigid modes excluded by overlap, not by projecting \(H\).
5. SVD of the stacked \(A\) **before** unlocking a new parameter. Drop or merge if \(\sigma_{\min}\) collapses.
6. Only if 2–4 say the residual is in the backbone/soft cage: consider 1–4 / torsion (historically rejected).

Restricted Wilson vs SVD, in one sentence: **Wilson supplies sparse directions \(g_a\)** (bonds/angles on 2–3 atoms). **SVD of \(A=\partial^2 E/\partial k\partial r^2\)** tells whether those parameters are independently seen by the data. Do **not** SVD-pseudoinverse the redundant \(B\) into an \(F\) matrix (the old artifact). Optional SVD of \(B\) is only to *report* \(\mathrm{rank}(B)\) — if it is \(3N-6\), you accidentally used a global basis; restrict \(B\).

### 4. Gate to own-min / FTIR

LS at DFT \(q\) is the \(k\) factory. FTIR plots require FIRE with **those locked \(k\)** then Hessian (same switches). C L1 first (except `cube_C` until mode-follow). Si L1: five tight minima in `pyscf_vib_results` are now valid own-min references; `octahedron_7ring_Si` is not. Do not anneal Morse+angles as a substitute for this ladder. Do not reuse Batch-1 Si RMSE.

### 5. Implementation order (code)

1. Hessian audit script: sum rule, rigid overlap, locality vs \(d\), dump PNG/JSON per case in the bank.
2. [x] `K_ab` from existing `build_wilson_matrix` rows (no \(C^+\), no Löwdin). \(Hg_a\) residual plots still open.
3. [x] Pin/freeze + stage CLI: `python3 tests/tSiNCs/test_FFfit.py --ladder si`. Default hybrid CLI is still `--legacy-hybrid` weights.
4. [x] Ladder path: `internal_weight=0`, `mode_weight=0`. Default `test_FFfit.py` without `--ladder` still uses mode=local=internal=1.
5. Then \(J\) / typed XH₂ only if \(\eta\) plots demand it (cube vs octa \(k_{\mathrm{XXX}}\) already see-saw).

---

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/molecular/FFfit.h` | **active** | \(A_p\), Wilson B, CSR graph, local mask, 1–4, dihedral FD batch, LSQ/GD |
| C++ | `cpp/libs/Molecular/FFfit_lib.cpp` | **active** | ctypes `libFFfit_lib.so` |
| C++ | `cpp/libs/Molecular/MMFF_lib.cpp` `getHessianContext` / `assembleHessianFromParams` | **active** | UFF H₂O-scale Hessian LS; **not** the hybrid Wilson pipeline |
| Python | `pyBall/FFfit.py` | **active** | ctypes + `assemble_hybrid_hessian_system` / `fit_hybrid_hessian` / `internal_hessian_projection` |
| Python | `pyBall/FFfit_utils.py` | **active** | types, topology, `load_pyscf_hessian_case`, `run_si_pbe_pin_ladder` |
| Python | `pyBall/FFfit_plots.py` | **active** | spectra, Wilson-diagonal plots (diagnostic), p5.js HTML maps |
| Python CLI | `tests/tSiNCs/test_FFfit.py` | **active** | `--ladder si` pin-ladder; default still spheres hybrid |
| Python | `tests/tSiNCs/fit_mmff_kss_pyscf.py` | **active** | own-min / frozen stretch RMSE; **not** hybrid LS |
| Python | `tests/tMMFF/test_hessian_fitting.py` | **experimental** | UFF `assembleHessianFromParams` on H₂O |
| Python tests | `tests/tSiNCs/test_fffit_hybrid.py` | **active** | hybrid/Wilson/\(K_{ab}\)/freeze |
| Python tests | `tests/tSiNCs/test_parity_py_cpp.py`, `test_parity_graph_cpp.py` | **active** | sensitivity + 14 graph/dihedral tests |
| JS | `pyBall/FFfit_plots.py` → p5.js HTML | **viz only** | `--stiffness-html`; not a fitter |
| OpenCL / WebGPU | — | **none** | no \(A_p\) kernel |
| backups | `pyBall/FFfit copy.py`, `FFfit_lib copy.cpp`, `FFfit copy.h`, `test_FFfit copy.py` | **deprecated** | do not edit |

Related (compute Hessian, not fit \(k\)): `pyBall/MMFF.py` `getHessian3Nx3N`; `pyBall/FTIR.py`; `pyBall/nanocrystal_pipeline.py`.

---

## Documents

| File | Kind | Use |
|------|-------|-----|
| **This file** | topical audit | inventory + checklist |
| [`../Topics/FFfit/HessianFitting_Theory.md`](../Topics/FFfit/HessianFitting_Theory.md) | theory SSOT | Wilson vs row-space, prestress, hybrid, subtypes |
| [`../Topics/FFfit/README.md`](../Topics/FFfit/README.md) | folder index | FFfit topic folder |
| [`../Topics/FFfit/Dihedral_Torsion.md`](../Topics/FFfit/Dihedral_Torsion.md) | topical audit | torsions tried / rejected |
| [`../Topics/FFfit/HessianFit.chat.md`](../Topics/FFfit/HessianFit.chat.md) | chat | Wilson-diagonal softness; may predate theory |
| [`../Topics/FTIR_Nanocrystals/Hessian_fitting.chat.md`](../Topics/FTIR_Nanocrystals/Hessian_fitting.chat.md) | chat | GPT 5.6: restricted Wilson, \(Hg_a\), **no rigid-\(P\)**, \(\eta_{ij}\); ~2893 |
| [`../../tests/tSiNCs/SiNCs_FFfit_summary.md`](../../tests/tSiNCs/SiNCs_FFfit_summary.md) | results | 2026-07-10 all-Si spheres 41→32 cm⁻¹ |
| [`../Topics/FTIR_Nanocrystals/MMFF_stiffness_scaling.md`](../Topics/FTIR_Nanocrystals/MMFF_stiffness_scaling.md) | report | cage \(k\) grids; angle×0.30 bug |
| [`../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md`](../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md) | report | own-min C L1; Morse/SA fail |
| [`../Topics/FTIR_Nanocrystals/Hessian_at_own_minimum.md`](../Topics/FTIR_Nanocrystals/Hessian_at_own_minimum.md) | HARD | FFfit vs spectrum |
| [`../../tests/tSiNCs/FFfit_python_to_cpp_port.plan.md`](../../tests/tSiNCs/FFfit_python_to_cpp_port.plan.md) | plan | graph/dihedral C++ port (mostly done) |
| [`molecular_topology.md`](molecular_topology.md) § Bond-Graph | audit | BFS / Wilson fill / parity tols |
| [`SiNCs.md`](SiNCs.md) §3.4 | briefing | one-paragraph pointer |
| [`Nanocrystal_Vibrations.md`](Nanocrystal_Vibrations.md) | audit | Hessian APIs, not the fitter checklist |

Generated plots (gitignored): `tests/tSiNCs/OUT_FFfit_plots/`, `OUT_mmff_kss_fit/`. References: `/home/prokop/SIMULATIONS/jobs_pyscf_vib_OUT_small_nc/results` (spheres); L1 PBE **canonical** `/home/prokop/SIMULATIONS/SiNCs/pyscf_vib_results`.

---

## Parity status

| Pair | Tolerance | Test |
|------|-----------|------|
| Python vs C++ sensitivity / Wilson / graph / dihedral FD | `atol=1e-8` graph/Wilson; `1e-6` FD Hessian | `test_parity_graph_cpp.py`, `test_parity_py_cpp.py` |
| Hybrid assembly / Wilson scaling / hierarchy rows | exact small systems | `test_fffit_hybrid.py` |

No GPU parity (no kernel).

---

## Open issues

- Full-Wilson default weight 1.0 on the *legacy* CLI is chemically misleading when \(\mathrm{rank}(B)=3N-6\). The pin-ladder (`--ladder si`) uses raw \(K_{ab}\) + graph-local Cartesian instead.
- `--ladder si` run 2026-09-03: hydride freeze holds; \(k_{\mathrm{SiSi}}\) D vs E **1.9%**; \(k_{\mathrm{XXX}}\) **19%** cube vs octa (see-saw, same family as C CH vs CH₂). Stretch RMSE at frozen \(q\) grows A→E because SiH₄ \(k_{\mathrm{XH}}\) is stiffer than L1 hydride stretches (2237 vs ~2080–2110). **Not FTIR** until FIRE+Hessian.
- C pin-ladder not run (would mix B3LYP small-mol with PBE L1 if stacked; pin-transfer only).
- Cross terms lack type-average prestress Hessian.
- `test_FFfit.py` still has leftover `fit_hessian` / GD imports; production path is `fit_hybrid_hessian` / `--ladder`.
- Directional \(Hg_a\) DFT driver not wired to CompChemUtils/PySCF jobs.
- Copy-files (`FFfit copy.*`) should stay unused or be deleted in a later cleanup (not this edit).





