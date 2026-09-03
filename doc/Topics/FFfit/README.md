# FFfit / Force-Field Fitting

This folder documents the PySCF-to-classical Hessian fitting workflow centred on `tests/tSiNCs/test_FFfit.py`. Pin-ladder: `python3 tests/tSiNCs/test_FFfit.py --ladder si` (local \(K_{ab}\), hydride freeze). The default CLI still uses the older hybrid mode/local/Wilson-row-space objective.

- [HessianFitting_Theory.md](HessianFitting_Theory.md) — authoritative derivation: Wilson diagnostic versus gauge-invariant fit, hierarchical regularization, and cross terms.
- [`../../topical_audit/Hessian_fitting.md`](../../topical_audit/Hessian_fitting.md) — **inventory + method checklist + idea catalog** (accuracy / robustness / interpretability / cheap DFT).
- [Dihedral_Torsion.md](Dihedral_Torsion.md) — torsion implementation/status audit; torsions are currently diagnostic rather than recommended Si-network terms.
- [HessianFit.chat.md](HessianFit.chat.md) — chronological development log; historical statements may predate the current theory document.
- [`../FTIR_Nanocrystals/Hessian_fitting.chat.md`](../FTIR_Nanocrystals/Hessian_fitting.chat.md) — GPT 5.6: restricted Wilson, \(Hg_a\), no rigid-mode projector, \(\eta_{ij}\).

Own-min MMFF stretch RMSE vs PySCF PBE L1 (FIRE then Hessian, not this hybrid Hessian objective): [`../FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md`](../FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md). One C pack did not unify cube CH₂ vs octahedron CH.
