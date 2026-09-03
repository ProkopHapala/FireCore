# MMFF vibration-fit scales — where 1.995 / 0.30 came from, and why angle×0.30 collapsed L2

**Intent (valid):** fit MMFF **stiffness** \(k\) so harmonic frequencies match QM. Equilibrium geometry \(l_0\), \(\theta_0\) must not change.  
**What actually happened:** C–H \(k\) was scaled correctly; “angle ×0.30” multiplied **\(\sin(\theta_0/2)\)**, not \(K\). FIRE then followed a broken tetrahedral target.

---

## Detective chain (June 2026 → 2026-09-02)

### 1. Small-molecule Hessian fit (the real idea)

Canonical report: [`CompChemUtils/examples/tSiNCs/MMFF_VIBRATION_FITTING_REPORT.md`](../../../../CompChemUtils/examples/tSiNCs/MMFF_VIBRATION_FITTING_REPORT.md) (**2026-06-11**).

| Molecule | Reference | What was fitted | Scripts |
|----------|-----------|-----------------|---------|
| CH₄ | PySCF B3LYP/cc-pVDZ | bond scale **1.995**, “angle scale” **0.30** | `fit_mmff_ch4.py` |
| C₂H₆ | same | C–H frozen 1.995; C–C **1.0**; angle **0.30** | `fit_mmff_c2h6.py` |
| Adamantane C₁₀H₁₆ | DFTB+ mio-1-1 | transferred those two numbers | `analyze_adamantane_modes.py` |

Grid: `sess.set_scales(sb, sa)` then `compute_frequencies` — **no FIRE in the grid**. Geometry was a frozen `relaxed.xyz` (CH₄) or ASE molecule. Criterion: RMSE of sorted frequencies vs QM (stretch window for C₂H₆/adamantane).

That is legitimate **FFfit** (shared geometry, vary \(k\)). Related phonon grid: `CompChemUtils/examples/phonons/` (`scale_bond=1.62`, `scale_angle=0.93` on **diamond**, same `apars[:,1]` bug; 0.93 is close to 1 so it did not implode a lattice).

### 2. What the fit actually multiplied

Live kernel (`cpp/common/molecular/MMFFsp3_loc.h`):

```text
cs0_ss = (apar.x, apar.y)   # unit pair: (cos(θ0/2), sin(θ0/2)) for evalAngleCosHalf
ssK    = apar.z             # angle stiffness Kss
```

Measured on CH₄ (this tree, `libMMFF_lib.so`):

| col | value | meaning |
|----:|------:|---------|
| 0 | 0.577 | \(\cos(54.74^\circ)=\cos(\theta_0/2)\), \(\theta_0=109.5^\circ\) |
| 1 | 0.816 | \(\sin(\theta_0/2)\), \(\sqrt{c^2+s^2}=1\) |
| 2 | 25.28 | **Kss** (builder `Kss*4` from `AngleTypes.dat` \(k=5.1\)) |
| 3 | −0.33 | π-related |

Bonds (correct): `bKs = 16.457` (C–H \(k\) from `BondTypes.dat`), `bLs = 1.095` (\(l_0\), **never scaled**).

Python / fitting code copied a **stale comment**:

```text
# pyBall/MMFF.py getBuffs:
apars  (c0, Kss, Ksp, c0_e)     # WRONG vs live kernel
# MMFFsp3_loc.h header default_NeighParams comment: same stale [c0, Kss, …]
# MMFFMolecularSession._angle_col = 1
# MMFF_VIBRATION_FITTING_REPORT.md: "apars[:,1] for angles"
# vib_utils.MMFFCalc: apars[:,1] *= scale_angle
# phonon_backends.MMFFPhononSession._angle_col = 1
# setAngleParamsByType(..., k=): apars[i, 1] = k
```

So the June fit’s “optimal angle 0.30” is **not** \(K \times 0.30\). It is \(\sin(\theta_0/2)\times 0.30\) at **frozen** QM/MMFF geometry. That changes the Hessian (prestress / wrong \(cs0\)) and can look like “softer bends” without moving atoms. C–H **1.995 on `bKs`** is a real stiffness scale and is why stretches jumped ~2240 → ~3100 cm⁻¹.

### 3. Transfer onto L2 nanocrystals (2026-09-02)

This chat reused CH×1.995 / angle×0.30 from that report as “adamantane-scaled MMFF” on L2 diamonds.

1. First bake: Hessian at **DFTB geometry**, no MMFF FIRE → not a spectrum ([`Hessian_at_own_minimum.md`](Hessian_at_own_minimum.md)).
2. Second bake: FIRE **with** those scales. Bond \(k\)×1.995 is harmless to the minimum. Angle ×0.30 on **col 1** denormalizes \((c,s)\) and moves \(\theta_0\). FIRE followed that → Rg 6.8→1.2 Å, C–C ~0.2 Å. **Not** “\(k\) collapsed diamond.”

---

## Correct protocol

Harmonic \(E=\tfrac12 k(q-q_0)^2\): the minimum is \(q_0\), independent of \(k\).

| Buffer | Scale this for frequencies | Never scale for a \(k\)-fit |
|--------|----------------------------|-----------------------------|
| `bKs` | C–H / C–C stiffness | — |
| `bLs` | — | \(l_0\) |
| `apars[:,2]` | \(K_{\mathrm{ss}}\) | — |
| `apars[:,0:2]` | — | \(\cos(\theta_0/2),\sin(\theta_0/2)\) |

After any “angle stiffness” scale: `hypot(apars[:,0], apars[:,1])` must stay 1 and match the unscaled pair. Fail loud otherwise.

**Do not reuse June 2026 angle×0.30** — that multiplied \(\sin(\theta_0/2)\), not \(K_{\mathrm{ss}}\).

Redo (2026-09-02): `tests/tSiNCs/fit_mmff_kss_pyscf.py`, frozen QM `relaxed.xyz`, scale `bKs` and `apars[:,2]` only. Output: `tests/tSiNCs/OUT_mmff_kss_fit/`.

| Fit | Reference | \(k_{\mathrm{XH}}\) | \(K_{\mathrm{ss}}\) | \(k_{\mathrm{XX}}\) | RMSE |
|-----|-----------|--------------------:|--------------------:|--------------------:|------|
| CH₄ | PySCF B3LYP/cc-pVDZ | **2.038** | **0.488** | — | 42 cm⁻¹ (all 9) |
| CH₄ | DFTB+ **3ob-3-1** | **1.812** | **0.488** | — | 45 cm⁻¹ (all 9) |
| C₂H₆ | PySCF; \(k_{\mathrm{CH}}\) frozen | 2.038 | 0.514 | 2.079 | 96 cm⁻¹ (ν>500) |
| SiH₄ | PySCF B3LYP/cc-pVDZ | **1.138** | **0.400** | — | 28 cm⁻¹ (all 9) |
| SiH₄ | DFTB+ **pbc-0-3** | **1.025** | **0.300** | — | 44 cm⁻¹ (all 9) |
| Si₂H₆ | PySCF; \(k_{\mathrm{SiH}}\) frozen 1.138 | 1.138 | 0.314 | **0.95** | 75 cm⁻¹ (ν>200) |
| Si₂H₆ | DFTB pbc; \(k_{\mathrm{SiH}}\) frozen 1.025 | 1.025 | 0.314 | **0.95** | 70 cm⁻¹ |
| adamantane X–H×K, CC=1 | DFTB 3ob-3-1 cage job | **1.76** | **0.68** | 1 | 129 all / 29 stretch |
| adamantane X–X×K, CH=1.76 | same | 1.76 | 0.68 | **2.20** | 87 all / 28 stretch |
| Si₁₀H₁₆ X–H×K | DFTB pbc-0-3 | 1.12 | 0.20† | 1 | 54 all / **79** stretch |

† Si₁₀ \(K=0.20\) is the **grid floor**. All-mode RMSE 348→54 by killing bends, but **Si–H stretch RMSE got worse** (default 44 → 79). Do not transfer that \(K\). Default MMFF Si–H stretches on Si₁₀ are already close (νmax 2102 vs pbc 2157).

Si `AtomTypes.dat` already has \(K_{\mathrm{ss}}=6.32\) (same as C → `apars.z=25.28`); \(\mathrm{hypot}(c,s)=1\). Default SiH₄ bends are too stiff (1330 vs PySCF 908), same pattern as CH₄.

**Transfer (FTIR / NC spectra), stiffness only, \(k_{\mathrm{XX}}=1\):**

| Crystal | SK | \(k_{\mathrm{XH}}\) | \(K_{\mathrm{ss}}\) |
|---------|----|--------------------:|--------------------:|
| C diamond L2 | 3ob-3-1 | 1.812 | 0.488 |
| Si NC | pbc-0-3 | 1.025 | 0.400 (SiH₄ PySCF) or 0.300 (SiH₄ pbc) |

Adamantane wants \(k_{\mathrm{CH}}\approx1.76\) (matches CH₄-3ob) and a slightly larger \(K\) (0.68 vs 0.49). Si₂H₆ wants \(k_{\mathrm{SiSi}}\approx1\) — unlike ethane’s C–C ridge.

Default MMFF vs PySCF CH₄ was RMSE **710** cm⁻¹ (bends ~1812/2062, stretches ~2107/2265). June wrong-column point (1.995, 0.30) is RMSE **244** on the new grid — not the minimum.

CH₄ vs PySCF, best vs default:

| | T₂ bend | E bend | A₁ stretch | T₂ stretch |
|--|--------:|-------:|-----------:|-----------:|
| PySCF | 1309 | 1531 | 3030 | 3152 |
| MMFF default | 1812 | 2062 | 2107 | 2265 |
| MMFF \(k\)-fit | 1300 | 1449 | 3008 | 3178 |

The leftover mismatch is the E/T₂ bend splitting (MMFF cannot open that gap with one \(K\)). Stretches are now in the 3000 cm⁻¹ window.

C–C on ethane is a **ridge**: at \(K=0.51\), \(k_{\mathrm{CC}}=1.0\) is RMSE 112 vs 96 at 2.08. Do **not** transfer \(k_{\mathrm{CC}}\times2\) onto diamond NCs. Diamond phonon grid (same session, old col-1 bug) liked bond ~1.62. For L2 FTIR, use CH₄ \(k_{\mathrm{CH}}\) + \(K_{\mathrm{ss}}\) and **leave C–C at 1**.

L2 nanocrystals were computed with DFTB+ **3ob-3-1**, so the transfer that belongs on those spectra is \(k_{\mathrm{CH}}=1.812\), \(K_{\mathrm{ss}}=0.488\), \(k_{\mathrm{CC}}=1\). PySCF \(k_{\mathrm{CH}}=2.04\) overshoots 3ob C–H (adamantane frozen-geom mean stretch 3137 vs 3ob 2931 cm⁻¹).

Spectrum of a scaled-\(k\) MMFF: relax with **those \(k\)** (geometry will barely move if only \(k\) changed), then Hessian. FFfit diagnostic: freeze geometry, vary \(k\), compare Hessians/frequencies — label it FFfit, not a new crystal.

### PySCF PBE/cc-pVDZ L1 NCs (2026-09-02)

Jobs: `/home/prokop/SIMULATIONS/SiNCs/pyscf_vib_results` (canonical PBE / ccECP-cc-pVDZ, 2026-09-03). Fitter: `tests/tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref joint`. Output: `tests/tSiNCs/OUT_mmff_kss_fit/pyscf_pbe_ncs/`. The Si pack below was fitted to **Batch-1 off-min** Hessians — do not cite as FTIR; re-run vs tight Si.

**General FF (one pack per element, not per crystal).** Criterion: minimax own-min stretch RMSE across all 6 C / 6 Si L1 jobs. FIRE with those \(k\), then Hessian. \(k_{\mathrm{XX}}=1\). \(K_{\mathrm{ss}}\) frozen at XH₄ all-mode (C 0.488, Si 0.400). 7-ring MMFF sometimes stops at a saddle (\(|F|\sim10^{-5}\), one imag ~−20…−55 cm⁻¹); kick 0.02 Å and re-FIRE.

| Element | \(k_{\mathrm{XH}}\) | \(K_{\mathrm{ss}}\) | worst stretch RMSE | mean stretch RMSE |
|---------|--------------------:|--------------------:|-------------------:|------------------:|
| C (all 6 L1) | **1.775** | 0.488 | 98 (cube_C) | 83 |
| Si (all 6 L1) | **1.097** | 0.400 | 71 (octahedron_5ring_Si) | 55 |

Cube vs octahedron still want different \(k\) (C cube ~1.88, octa ~1.68; Si 5-ring softer hydrides). Minimax is the compromise that keeps every crystal’s stretch RMSE under ~100 cm⁻¹. Si 5-ring vanilla was already close (RMSE 29–48); the general \(k\) slightly **hurts** those two (→57–71) in exchange for cube/7-ring/octa (155–160 → 40–65).

Whole-spectrum overlays (vanilla / fitted / PySCF): `pyscf_pbe_ncs/<case>/ownmin/spectra_overlay.png` and `gallery_ownmin_{C,Si}.png`. Fitted MMFF moves C–H from ~2200 to ~2900 and Si–H from ~2100 to ~2230; it also clears the fake MMFF band in the 1500–2300 (C) / 1000–1600 (Si) gap. Fingerprint (≲1500) is improved but not a spectroscopic match — \(K_{\mathrm{ss}}\) was not jointly fitted.

Do **not** put B3LYP-CH₄ \(k=2.038\) on these PBE NCs. DFTB-3ob \(k=1.812\) remains the right transfer for **L2 3ob** spectra. For this PBE L1 set use **C \(k_{\mathrm{CH}}=1.775\)** as the bonded-only starting pack — **not** as FTIR-ready. Si \(k_{\mathrm{SiH}}=1.097\) is from off-min Batch-1 refs; re-fit after the tight bank.

Cube vs octahedron still disagree (~100 cm⁻¹ shared minimax). Surface H···H + Morse α + split XXX/XXH/HXH **do not** close that gap (SA mean 72, see-saw). Full numbers, CLI (`nbscan` / `morse2d` / `anneal`), and plots: [`MMFF_C_CH_vs_CH2_kfit.md`](MMFF_C_CH_vs_CH2_kfit.md). Si FTIR: five tight minima exist; do not reuse Batch-1 Si RMSE.

Spectrum of a scaled-\(k\) MMFF: relax with **those \(k\)** (geometry will barely move if only \(k\) changed), then Hessian. FFfit diagnostic: freeze geometry, vary \(k\), compare Hessians/frequencies — label it FFfit, not a new crystal.

---

## Code that had to change

- `CompChemUtils/examples/tSiNCs/mmff_molecular_session.py` — `_angle_col = 2`
- `pyBall/MMFF.py` — buffer comment + `setAngleParamsByType`
- `CompChemUtils/examples/tSiNCs/vib_utils.py`
- `tests/tSiNCs/fit_mmff_kss_pyscf.py` — CH₄/C₂H₆/SiH₄/Si₂H₆/adamantane/Si₁₀H₁₆ \(k\)-grids on `apars[:,2]`; `--mol pyscf_nc` vs PBE L1 (`joint` / `nbscan` / `morse2d` / `anneal`)
