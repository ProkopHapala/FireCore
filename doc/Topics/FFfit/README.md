# FFfit / Force-Field Fitting

This folder documents the PySCF-to-classical Hessian fitting workflow centred on `tests/tSiNCs/test_FFfit.py`. The fitted model uses a hybrid mode/local/Wilson-row-space objective, transferable Si environment types with hierarchy regularization, and optional local valence cross terms.

- [HessianFitting_Theory.md](HessianFitting_Theory.md) — authoritative derivation: Wilson diagnostic versus gauge-invariant fit, hierarchical regularization, and cross terms.
- [Dihedral_Torsion.md](Dihedral_Torsion.md) — torsion implementation/status audit; torsions are currently diagnostic rather than recommended Si-network terms.
- [HessianFit.chat.md](HessianFit.chat.md) — chronological development log; historical statements may predate the current theory document.
