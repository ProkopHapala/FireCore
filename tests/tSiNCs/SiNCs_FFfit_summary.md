
# FFfit Hessian Fitting — Results Report

## 1. Systems Studied

6 Si nanocrystals passivated with H, from PySCF DFT reference Hessians:

| System | Atoms | Bonds | Angles | Si Types Present |
|--------|------:|------:|-------:|------------------|
| SiH4 | 5 | 4 | 6 | SiH3 |
| Si_R3p8 | 26 | 28 | 60 | SiH, SiH2 |
| Si_R4p5 | 47 | 52 | 114 | Si, SiH, SiH2 |
| Si_R5p0 | 50 | 58 | 132 | Si, SiH, SiH2 |
| Si_R5p5 | 68 | 82 | 192 | Si, SiH, SiH2 |
| Si_R6p0 | 76 | 92 | 216 | Si, SiH, SiH2 |

## 2. Three Typing Strategies Compared

### Strategy A: Elemental (5 params)
- Bond types: `H-Si`, `Si-Si` (2)
- Angle types: `H-Si-H`, `H-Si-Si`, `Si-Si-Si` (3)

### Strategy B: Central-only Si subtyping (16 params)
- Bond types: `H-SiH`, `H-SiH2`, `H-SiH3`, `Si-Si`, `Si-SiH`, `Si-SiH2`, `SiH-SiH`, `SiH-SiH2`, `SiH2-SiH2` (9)
- Angle types: keyed by **central atom type only** — `H-SiH3-H`, `H-SiH2-H`, `H-SiH2-Si`, `H-SiH-Si`, `Si-SiH-Si`, `Si-SiH2-Si`, `Si-Si-Si` (7)

### Strategy C: Full Si subtyping (32 params)
- Same 9 bond types
- Angle types: keyed by **all three atom types** — 23 distinct angle types

## 3. Frequency Fit Results

### Headline: **Central-only subtyping (16 params) is the winner**

| Typing | Params | Mean RMSE cm⁻¹ | Mean MAE cm⁻¹ | Mean relFrob % |
|--------|-------:|---------------:|--------------:|---------------:|
| Elemental | 5 | 37.12 | 24.44 | 6.56 |
| **Si subtyping (central)** | **16** | **29.42** | **19.95** | 7.29 |
| Si subtyping (full) | 32 | 29.76 | 20.32 | 7.29 |

### Per-system breakdown (central-only subtyping):

| System | Elemental RMSE | Subtyped RMSE | Improvement |
|--------|---------------:|--------------:|------------:|
| SiH4 | 54.10 | **29.98** | **45%** |
| Si_R3p8 | 33.58 | **27.95** | 17% |
| Si_R4p5 | 39.26 | 38.98 | 1% |
| Si_R5p0 | 28.50 | **23.02** | 19% |
| Si_R5p5 | 33.84 | **28.84** | 15% |
| Si_R6p0 | 33.46 | **27.74** | 17% |

**Mean improvement: 21% in RMSE, 18% in MAE.**

### Key observations:

- **Full subtyping (32 params) ≈ central-only (16 params)** — doubling parameter count gives negligible improvement. The extra angle types are underdetermined (many have only 3-12 observations).
- **relFrob slightly worsens** with subtyping (6.56% → 7.29%). This is expected: the hybrid fit optimizes vibrational modes, not raw Hessian elements. Better mode fitting can introduce local Hessian deviations.
- **SiH4 improves most** (54→30 cm⁻¹, 45%) — the small molecule benefits enormously from distinguishing `H-SiH3-H` angle stiffness (2.51 eV/rad²) from the generic `H-Si-H` (2.21 eV/rad²).
- **Si_R4p5 barely changes** — this system has 3 bulk Si atoms and a more complex surface, suggesting the model needs additional coupling terms to improve further.

## 4. Fitted Parameters (central-only, 16 params)

### Bond stiffnesses (eV/Å²)

| Bond Type | k (eV/Å²) | Count | r₀ (Å) |
|-----------|----------:|------:|-------:|
| H-SiH3 | 18.26 | 4 | 1.480 |
| H-SiH2 | 18.59 | 96 | 1.473 |
| H-SiH | 18.13 | 52 | 1.478 |
| Si-SiH | 10.36 | 42 | 2.354 |
| Si-SiH2 | 9.88 | 18 | 2.353 |
| SiH-SiH2 | 9.74 | 66 | 2.350 |
| SiH2-SiH2 | 9.68 | 6 | 2.348 |
| SiH-SiH | 9.63 | 24 | 2.352 |
| Si-Si (bulk) | 9.91 | 8 | 2.355 |

**Physical interpretation**: All bonds positive and physical. H-Si bonds are ~18 eV/Å² (stiff), Si-Si bonds ~9.6-10.4 eV/Å² (softer). The Si-SiH bond (surface Si bonded to 1 H) is slightly stiffer than bulk Si-Si (10.36 vs 9.91), consistent with surface reconstruction effects. Bond lengths vary by <0.3% within Si-Si family.

### Angle stiffnesses (eV/rad²)

| Angle Type | k (eV/rad²) | Count | θ₀ (deg) |
|------------|------------:|------:|----------:|
| H-SiH3-H | 2.51 | 6 | 109.47 |
| H-SiH2-H | 2.37 | 48 | 108.86 |
| H-SiH2-Si | 1.35 | 192 | 109.61 |
| H-SiH-Si | 1.64 | 156 | 109.41 |
| Si-SiH2-Si | 1.04 | 48 | 109.54 |
| Si-SiH-Si | 0.97 | 156 | 109.53 |
| Si-Si-Si (bulk) | 0.99 | 114 | 109.47 |

**Physical interpretation**: All angles positive and physical. Clear hierarchy:
- **H-Si-H angles stiffest** (~2.4-2.5 eV/rad²) — H-Si-H bending is strongly constrained
- **H-Si-Si intermediate** (~1.4-1.6 eV/rad²) — mixed surface angle
- **Si-Si-Si softest** (~1.0 eV/rad²) — bulk tetrahedral angle is floppier

The `H-SiH2-H` angle (108.86°, k=2.37) is slightly more compressed and softer than `H-SiH3-H` (109.47°, k=2.51), consistent with the additional Si neighbors providing geometric constraint.

## 5. Model Comparison: 1-4 and Torsions

From `@/home/prokop/git/FireCore/tests/tSiNCs/OUT_FFfit_plots/model_comparison_Si_R3p8/model_comparison.md`:

| Model | Params | RMSE cm⁻¹ | relFrob % |
|-------|-------:|----------:|----------:|
| Bond + angle | 5 | 39.01 | 6.29 |
| + 1-4 springs | 6 | 38.37 | 6.20 |
| + UFF torsion | 8 | 39.01 | 6.29 |
| + both | 9 | 38.37 | 6.20 |

**Findings**:
- **1-4 springs give marginal improvement**: 0.6 cm⁻¹ RMSE reduction (1.6%). The fitted 1-4 stiffnesses are tiny (0.07-0.10 eV/Å² for Si-Si, ~0.005 for H-Si, ~0.07 for H-H).
- **UFF torsions contribute exactly zero**: fitted V values are ~1e-15 to ~1e-25 eV — the optimizer drives them to zero. The torsion sensitivity matrices are not positive semidefinite (prestress term `f'·C` produces negative eigenvalues), so any nonzero V degrades the fit.
- **Torsion + 1-4 = 1-4 alone**: torsions add nothing on top of 1-4 springs.

## 6. Spectrum Quality

From the spectrum plots:

- **SiH4**: Subtyping dramatically improves the spectrum. The elemental model misses the bending mode region (~950 cm⁻¹) by ~50 cm⁻¹; subtyping brings it to ~30 cm⁻¹ error. The Si-H stretch region (~2200 cm⁻¹) is well-matched in both cases.
- **Si_R3p8 through Si_R6p0**: The spectral envelope is well reproduced. Main deviations are in the 800-1000 cm⁻¹ region (Si-Si-Si bending modes) where the reduced model lacks sufficient coupling terms. The Si-H stretch peak (~2100 cm⁻¹) is accurately matched.
- **Larger nanocrystals** (R5.5, R6.0) show more low-frequency modes that the model captures well, with the subtyped model showing visibly tighter peak alignment.

## 7. What Works and Why

### Wilson B Matrix / GF Method
The internal-coordinate projection (`F = C^{+T} D C^+` where `C = B·M^{-1/2}`) provides a physically meaningful decomposition of the Hessian into valence force constants. This:
- Automatically removes rigid translations/rotations (Wilson matrix annihilates them)
- Gives dimensionless internal coordinates (δr/r₀ for bonds, radians for angles)
- Allows the fit to target physically meaningful force constants rather than raw Cartesian matrix elements

### Hybrid Objective
The three-component hybrid (mode-basis + local-Hessian + Wilson internal) is the key to the good fits:
1. **Mode objective** ensures vibrational frequencies are correct — this is what matters physically
2. **Local Hessian** prevents the model from ignoring short-range matrix elements
3. **Wilson internal** ensures the fitted force constants are physically meaningful in valence space

The mode objective alone would give good frequencies but unphysical force constants; the local and internal components regularize this.

### Si Type Splitting
The central-only subtyping works because:
- **Bond stiffness varies with coordination**: H-SiH3 (fully hydrogenated) vs H-SiH (one H, three Si neighbors) differ by ~0.5 eV/Å²
- **Angle stiffness varies strongly with central atom type**: H-SiH3-H (2.51) vs Si-SiH-Si (0.97) differ by 2.5×
- The central atom determines the local electronic environment, which controls bending stiffness
- Full subtyping (all three atoms) overfits: many angle types have <10 observations, and the extra parameters don't improve frequencies

## 8. Limitations and Open Issues

- **relFrob ~6-8%**: The reduced bond+angle model captures ~93% of the Hessian spectral content. The remaining ~7% requires coupling terms (stretch-stretch, stretch-bend, bend-bend).
- **Unmatched modes**: 6-12 modes per system can't be matched one-to-one between DFT and model — these are likely low-frequency collective modes requiring longer-range interactions.
- **UFF torsions are broken**: The prestress term `f'·C` makes the torsion sensitivity indefinite. The proposed fix (pure curvature `f''·b⊗b` only) is documented but not yet implemented.
- **1-4 springs are negligible**: The fitted values are ~0.1 eV/Å² — too small to matter. This confirms that for these Si nanocrystals, the physics is dominated by nearest-neighbor bond and angle terms.
- **Si_R4p5 is an outlier**: Only 1% improvement from subtyping, possibly because it has 3 bulk Si atoms creating modes that the surface-typed model can't capture well.

## 9. Summary

| Aspect | Status | Quality |
|--------|--------|---------|
| Bond+angle Hessian fitting | Working | relFrob ~6-8%, freq RMSE ~30-40 cm⁻¹ |
| Si type splitting (central) | **Working** | **21% RMSE improvement over elemental** |
| Si type splitting (full) | Working but overfits | Same as central with 2× params |
| Wilson B matrix / GF method | Working | Provides physical internal-coordinate projection |
| Hybrid mode+local+internal fit | Working | Best method, all three components contribute |
| 1-4 interactions | Negligible | <2% improvement, k ~0.1 eV/Å² |
| UFF torsions | **Broken** | V→0, prestress makes sensitivity indefinite |
| All fitted parameters | Physical | All positive, clear physical trends |