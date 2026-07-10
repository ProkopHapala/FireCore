---
type: ResultsReport
title: Si Nanocrystal Hessian Force-Field Fit
description: Current all-Si PySCF benchmark for hierarchical bond-angle fitting and optional local valence couplings.
tags: [fffit, silicon, hessian, vibration, forcefield]
timestamp: 2026-07-10
---

# Si Nanocrystal FFfit Results

This report summarizes the current six-system PySCF benchmark: `SiH4`, `Si_R3p8`, `Si_R4p5`, `Si_R5p0`, `Si_R5p5`, and `Si_R6p0`. The model is fitted to the DFT Hessians with the hybrid mode, graph-local Hessian, and orthonormal Wilson row-space objective described in [`HessianFitting_Theory.md`](../../doc/Topics/FFfit/HessianFitting_Theory.md).

## Model Definitions

| Model | Parameters | Physical constraint |
|---|---:|---|
| Elemental bond + angle | 5 | One parameter per element bond/angle family |
| Hierarchical Si subtypes | 16 | `Si`, `SiH`, `SiH2`, `SiH3` subtype offsets shrunk to a data-supported elemental-family mean |
| Hierarchical + stretch--stretch | 23 | Adds one signed local $K_{rr}$ per angle type |
| Hierarchical + both cross terms | 30 | Adds signed local $K_{rr}$ and symmetric $K_{r\theta}$ per angle type |

The least-norm Wilson matrix $F=C^{+T}DC^+$ is a diagnostic only. Its diagonal is not a unique DFT bond stiffness in a redundant bond-angle basis and is not used as the fitted parameter target.

## All-Si Comparison

| Model | Mean RMSE cm$^{-1}$ | Mean MAE cm$^{-1}$ | Mean Hessian relFrob % | Condition |
|---|---:|---:|---:|---:|
| Elemental bond + angle | 40.96 | 27.15 | 6.86 | 1.32 |
| Hierarchical Si subtypes | 31.59 | 20.44 | 6.86 | 1.61 |
| Hierarchical + stretch--stretch | 30.46 | 19.64 | 6.74 | 1.62 |
| Hierarchical + stretch--stretch + stretch--bend | 31.09 | 19.66 | **5.86** | 1.67 |

The subtype hierarchy is the main spectrum improvement: it lowers mean frequency RMSE by about 23% without worsening the global Hessian error. Stretch--bend coupling is the useful optional Hessian refinement; stretch--stretch alone is nearly negligible for the Si network.

## Transferable Parameters

For the hierarchical bond-angle fit, Si-Si subtype springs remain physically coherent despite the sparse bulk class:

| Si-Si subtype | $k$ eV/Å$^2$ | Observations |
|---|---:|---:|
| Si-Si | 9.17 | 8 |
| Si-SiH | 9.42 | 42 |
| Si-SiH2 | 9.57 | 18 |
| SiH-SiH | 9.25 | 24 |
| SiH-SiH2 | 9.79 | 66 |
| SiH2-SiH2 | 9.54 | 6 |

The corresponding elemental fit gives $k_{\rm SiSi}=8.34$ eV/Å$^2$. The constrained subtype values are therefore interpreted as modest environment corrections, not evidence for anomalously soft bulk bonds.

## Optional Cross Terms

The optional local terms are

$$E_{rr}=K_{rr}q_1q_2, \qquad E_{r\theta}=K_{r\theta}(q_1+q_2)\Delta\theta/\sqrt2,$$

with $q=\Delta r/r_0$. They are signed coefficients in eV, zero-centered regularized, bounded, and valid only with local equilibrium coordinates. In the all-Si fit, Si--Si--Si cross coefficients are small ($K_{rr}\approx0.012$ eV and $K_{r\theta}\approx4\times10^{-4}$ eV); the largest fitted couplings are H--Si--H stretch--bend terms of about 0.43--0.52 eV.

## Results and Reproduction

- [Hierarchical subtype outputs](OUT_FFfit_plots/typing_comparison_all_Si_hierarchical_joint_v2/) — parameters, per-system frequency tables, and six spectrum overlays.
- [Stretch--stretch outputs](OUT_FFfit_plots/typing_all_Si_cross_stretch_stretch/) — optional $K_{rr}$ comparison.
- [Both cross-term outputs](OUT_FFfit_plots/typing_all_Si_cross_both/) — best global Hessian comparison.

```bash
CPP_BUILD_PATH=$PWD/cpp/Build-opt/libs python3 tests/tSiNCs/test_FFfit.py \
  --cases all_Si --compare-typing --subtype-shrinkage 0.001 \
  --stretch-stretch --stretch-bend \
  --plot-dir tests/tSiNCs/OUT_FFfit_plots/typing_all_Si_cross_both
```

## Scope and Remaining Limits

- UFF torsions and 1-4 distance springs remain optional diagnostics, not recommended Si-network production terms.
- Cross terms currently require `--equilibrium local`; their type-average prestress Hessian is intentionally not approximated.
- A lower fit residual is accepted only after checking model stability, condition number, and held-out/independent spectra.
