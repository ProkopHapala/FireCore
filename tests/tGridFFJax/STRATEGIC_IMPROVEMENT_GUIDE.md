# Strategic Improvement Guide: Metal GridFF Framework

## Purpose
This document records the **scientific reasoning** behind each strategic decision in
improving GridFF accuracy for molecule-metal surface interactions. It serves as a
template for extending the framework to other MLIPs, other metals, and other systems.

---

## 1. The Fundamental Physics: Why Metals Are Different from NaCl

### NaCl (ionic surface) — solved case
The NaCl interaction is fundamentally **pairwise**:
```
E_total = Σ_ij [ Morse(r_ij) + q_i*q_j/r_ij ]
```
- **Pauli (P):** Morse repulsion from each Na⁺/Cl⁻ ion. Error vs LAMMPS: ~1e-8 eV
- **London (L):** Morse attraction from each ion. Error vs LAMMPS: ~1e-8 eV
- **Coulomb (Q):** Point charges on ions + Ewald sum. Error vs LAMMPS: ~1e-5 eV
- The PLQ grid faithfully tabulates these pairwise sums on a 3D grid
- B-spline interpolation introduces < 1e-5 eV error at 0.1 A grid spacing

### Metals (Ag, Cu, Au) — the challenge
The metal-molecule interaction is fundamentally **NOT pairwise**:

1. **Delocalized electrons:** No point charges. The metallic electron gas creates a
   smooth, collective potential that cannot be decomposed into atom-atom pairs.

2. **Pauli repulsion:** Arises from overlap of adsorbate electron cloud with the
   metallic electron density (Jellium model). Scales roughly as the local electron
   density → our `V_P = rho^p` approach is physically motivated.

3. **van der Waals dispersion:** Has THREE distinct regimes:
   - **Short range (z < 3 A):** Density-overlap attraction. Our `V_L = -rho^0.5`
     captures this regime. Physically: electron-electron correlation in the
     overlap region.
   - **Medium range (3-8 A):** Screened pairwise C₆/R⁶ dispersion. Each metal
     atom contributes C₆/R⁶ attraction, but metallic screening reduces C₆ by
     30-50% vs free atoms (Tkatchenko, DiStasio, Car, Scheffler, PRL 2012).
     This regime has **site dependence** (hollow > bridge > top) because hollow
     sites have more Ag neighbors contributing.
   - **Long range (z > 8 A):** Non-retarded vdW: C₃/z³ asymptote from the
     Lifshitz-Zaremba-Kohn (LZK) theory. For Ag(111), C₃ ≈ 1.5 eV·A³.

4. **Image charge interaction:** A charge q at height z above a metal surface
   induces an image charge -q at depth z, giving E_image = -q²/(16πε₀z).
   For CHONH2 partial charges at z=3 A, this is ~50-200 meV — significant!

5. **Pillow effect / Pauli pushback:** The adsorbate compresses the surface
   electron spillover, reducing the surface dipole. This is captured partially
   by the LOCPOT (if computed with adsorbate present) but our clean-slab LOCPOT
   misses this.

6. **d-band hybridization:** At chemisorption distances (z < 2.5 A), the
   adsorbate frontier orbitals hybridize with the metal d-band. The d-band
   center varies by site (hollow vs top), creating site-specific bonding.
   This cannot be captured by any smooth-density-based approach.

---

## 2. Strategic Improvement Timeline

### Phase 1: Fix the Machinery (completed, 81.8 → 64.5 meV)

**Goal:** Ensure the existing code works correctly before adding physics.

**Diagnosis:** The B-spline coefficient overwrite bug meant we were evaluating
raw grid values with B-spline kernels, causing systematic smoothing errors.

**Strategy:** Always fix bugs before adding features. A broken implementation
can mask the true error floor of the functional form.

**Lesson:** Interpolation bugs are silent — they don't crash, they just give
wrong (but plausible-looking) results. Comparing B-spline vs trilinear output
is a good sanity check.

### Phase 2: Systematic Parameter Exploration (completed, 64.5 meV plateau)

**Goal:** Find the optimal existing parameters before changing the functional form.

**What we swept:**
- Density powers: pauli_power = {1.0, 1.5, 2.0}, london_power = {0.5, 0.65, 0.75, 0.8, 1.0}
- Direct PLQ scalars vs REQ Morse coupling
- Additional degrees of freedom: z-shift, static charges
- Optimizer: scipy L-BFGS-B vs JAX optax Adam

**Key findings:**
- Default powers (pauli=1.0, london=0.5) are optimal — changing them worsens fit
- REQ Morse coupling is essential (direct PLQ: 190 meV)
- Adding parameters (z-shift, charges) doesn't help — the error is in the SHAPE, not the scale
- JAX optax is 79x faster than scipy (21.5s vs 1700s) with identical results

**Strategy:** Exhaust the parameter space of the existing model before modifying
the functional form. This establishes the baseline ceiling and identifies WHAT
is wrong (shape vs magnitude vs position).

**Lesson:** When adding parameters doesn't help, the problem is in the functional
form, not the parameterization. This is a crucial diagnostic.

### Phase 3: Physics-Based Functional Form Fix (completed, 64.5 → 36.7 meV)

**Goal:** Fix the London tail decay problem identified in Phase 2.

**Diagnosis from Phase 2:**
- The z-region analysis showed -90 meV bias at z=3.5-5.0 A at ALL sites
- This is systematic, site-independent → it's in the London field shape
- Physics: `V_L = -rho^0.5` decays as `exp(-kappa*z)`, but true vdW decays
  as `C₆/R⁶` which is faster at medium range

**Solution:** Fermi damping function `f(z) = 1/(1+exp((z-d₀)/w))` applied
to the London grid during construction. Physically motivated by:
- Tkatchenko-Scheffler vdW framework (PRL 2009)
- DFT-D3 Becke-Johnson damping (JCC 2011)
- The density-overlap London is valid only where electron clouds overlap;
  beyond that, the true interaction is pairwise C₆/R⁶

**Optimization approach:** Three-stage sweep:
1. Coarse: 10 d₀ values, 1 width → identify optimum region (d₀ ≈ 3.5-4.0)
2. Fine: 32 combos (8 d₀ × 4 widths) → narrow to d₀ = 3.5-3.8, w = 0.3-0.5
3. Ultra-fine: 48 combos → optimal d₀ = 3.70, w = 0.35 (robust, not fragile)

**Key: Run MACE teacher ONCE, sweep grid parameters.** The damping changes
only the grid builder, not the pose evaluations. This made the 48-combo
ultra-fine sweep take ~8 minutes instead of ~48 × 3 minutes.

**Strategy:** When you identify a physics-based correction, implement it as
a continuous parameter that can be swept. Use hierarchical sweeps (coarse →
fine → ultra-fine) to find the optimum efficiently.

**Lesson:** The optimal d₀ = 3.70 A is consistent with the sum of vdW radii
for C-Ag (~3.5-3.6 A from TS parameters), confirming the physical basis.

### Phase 4: Remaining Error Diagnosis (current)

**Current error breakdown (z ≥ 2.0 A, 36.7 meV):**

| Source | RMSE | % of total | Fixable? |
|--------|------|-----------|----------|
| Hollow well depth | ~102 meV at well | ~60% | Needs pairwise C₆ |
| Short-range shape (z < 2.5) | ~120 meV | ~25% | Needs better P form |
| Bridge/top residual | ~10 meV | ~10% | Near floor |
| Tail (z > 3.5) | < 3 meV | < 5% | Solved |

**Root cause of hollow-site error:**
The PLQ fields are derived from the metallic electron density, which is
laterally smooth (metallic bonding → delocalized electrons). The PLQ field
values at top/bridge/hollow sites at the same z height differ by only ~5-10%,
but MACE gives hollow 30% deeper well than top.

The physical reason: at hollow sites, the molecule sits above a 3-fold coordinated
pocket with MORE Ag neighbors at close range. A pairwise sum over Ag atoms would
naturally give stronger interaction at hollow. The smooth density doesn't capture
this coordination-number effect.

---

## 3. The Path Forward: What Physics Is Missing

### 3.1 Pairwise C₆/R⁶ Dispersion (Medium Priority — fixes hollow site)

The density-based London works at short range but fails at medium range.
The pairwise C₆/R⁶ sum works at medium range but double-counts at short range.
The correct approach: **blend them.**

```
V_vdW(r) = V_density_London(r) * f_short(z) + V_pairwise_C6(r) * f_long(z)
```

where `f_short + f_long ≈ 1` is enforced by complementary Fermi/Tang-Toennies functions.

The pairwise term `V_pairwise_C6(r) = -Σ_j C₆(Ag-X) / |r - r_j|⁶` naturally
counts Ag neighbors → differentiates hollow from top.

**Implementation path:**
1. The pairwise Morse infrastructure already exists (`accumulate_pairwise_fields`)
2. Need to replace Morse with C₆/R⁶ functional form
3. Need screened C₆ values for Ag-{C,H,N,O} from vdW-surf (Ruiz et al. 2012)
4. Blend with density-London using complementary damping functions
5. The blend point is our existing d₀ parameter

**Expected improvement:** 20-50 meV reduction in hollow error → overall
RMSE potentially below 25 meV.

### 3.2 Image Charge Interaction (Medium Priority — adds ~50-200 meV physics)

For polar molecules like CHONH2 with partial charges (N: -0.8e, H: +0.3e, etc.),
the classical image charge interaction with a metal surface is:

```
E_image = -Σ_i q_i² / (16π ε₀ · 2 · (z_i - z_image))
```

This is already implemented in the framework (`use_image_charge_fixed`) but
currently disabled. The key parameters:
- `image_plane`: z-position of the image plane (typically 0.5-1.0 A above top layer)
- `image_damping`: regularization at short distances
- `image_scale`: overall scaling

**Why it matters:** For CHONH2 at z=3 A, with charges ~0.3-0.8 e, the
image interaction is ~50-200 meV. This is comparable to our current error!

**Implementation path:**
1. Enable `use_image_charge_fixed = True`
2. Fit `image_plane`, `image_scale` as free parameters alongside REQ
3. The LOCPOT already captures some image-charge-like physics, so `image_scale`
   may end up less than 1.0 (avoiding double-counting)

**Expected improvement:** 10-30 meV if not double-counting with LOCPOT.

### 3.3 Element-Specific Density Powers (Low Priority)

Currently all elements share the same Pauli and London powers. Different
adsorbate atoms (H vs C vs O vs N) may have different electron density
profile shapes, leading to different optimal powers.

This adds 2×N_elements parameters but may help the short-range shape.

### 3.4 Anisotropic Corrections (Future)

Some adsorbate atoms (especially O, N with lone pairs) have anisotropic
charge distributions. The current PLQ treats each atom as a point sampler
of the grid fields. An angular-dependent correction could help for
oriented molecules.

---

## 4. Generic Framework Design Principles

### 4.1 Hierarchy of Interactions (ordered by importance for metals)

```
E_total = E_Pauli(rho) + E_London_short(rho) + E_London_long(C₆/R⁶)
        + E_Coulomb(LOCPOT) + E_image(q²/z) + E_reactive(site-specific)
```

Each term has a physical basis, is independently testable, and adds
predictive power. No term is a "patch" — each addresses a distinct
physical mechanism.

### 4.2 Validation Protocol

For each new term:
1. Show the isolated contribution (turn off other terms, check sign/magnitude)
2. Show it improves the RIGHT z-region (don't improve tail while worsening well)
3. Check parameter sensitivity (is the optimum robust or fragile?)
4. Check transferability (does it work for multiple molecules/sites?)

### 4.3 When to Add vs Subtract Complexity

**Add a term when:**
- There's a systematic bias in a specific z-region
- The bias has a known physical origin
- The existing parameters can't compensate (Phase 2 test)

**Don't add a term when:**
- The error is random/unsystematic
- The improvement requires fine-tuned, fragile parameters
- The same effect is already partially captured (double-counting risk)

---

## 5. Metal-Specific Damping Parameters

The London damping d₀ should scale with the metal's electron density decay length.
For the Fermi damping `f(z) = 1/(1+exp((z-d₀)/w))`:

| Metal | Work function (eV) | Decay κ (A⁻¹) | Estimated d₀ (A) | Status |
|-------|-------------------|---------------|-------------------|--------|
| Ag(111) | 4.74 | 1.12 | 3.70 | Optimized |
| Cu(111) | 4.94 | 1.14 | 3.55 | To test |
| Au(111) | 5.31 | 1.18 | 3.40 | To test |

The trend: higher work function → faster electron density decay → damping
starts closer to the surface.

The width parameter w ≈ 0.35 A should be similar across metals (it reflects
the transition width between density-overlap and pairwise-C₆ regimes, which
is a property of the adsorbate electron cloud, not the metal).

---

## 6. Component Error Budget

To understand where to focus improvement effort, we need to know the
error contribution from each PLQ channel. The method:

1. Fix optimal parameters from the best fit
2. Evaluate P, L, Q contributions at each z-point
3. Compare total GridFF against MACE teacher
4. The difference = what's MISSING from PLQ

This analysis answers: "Is the Pauli shape wrong, or the London magnitude
wrong, or is there a missing interaction (image, C₆, site-specific)?"

**Key insight:** If P and L shapes match well but the total is wrong,
then a missing interaction (image charge, site correction) is the cause.
If individual channel shapes are wrong, then the functional form needs
improvement.

---

## 7. Data Files and Reproducibility

Every result is traceable:
- CHGCAR/LOCPOT: `/home/niel/git/ORR_HER_Ag_Colab/.../slab_clean/...`
- MACE model: `tests/tGridFFJax/mad-surf_data/models/.../MACE_model.model`
- Sweep results: `tests/tGridFFJax/damping_sweep_ultra/`
- Optimal run: `tests/tGridFFJax/gridff_proof_run_damped_optimal/`
- Benchmark doc: `tests/tGridFFJax/METAL_GRIDFF_BENCHMARK.md`

Scripts print full paths at startup and save configs as JSON for reproducibility.

---

## 8. Summary: Error Reduction Roadmap

```
81.8 meV — Initial (B-spline bug)
  ↓ -21%  Bug fix (correct B-spline prefiltering)
64.5 meV — Baseline (density-only PLQ)
  ↓ -43%  London Fermi damping (d₀=3.70, w=0.35)
36.7 meV — Current best ← WE ARE HERE
  ↓ ~-30%  Pairwise C₆/R⁶ blending (fixes hollow site)
~25 meV  — Projected with C₆ blending
  ↓ ~-20%  Image charge interaction
~20 meV  — Projected with image charges
  ↓ ~-25%  Site-specific reactive corrections
~15 meV  — Projected theoretical floor for PLQ+corrections
```

The theoretical floor of ~15 meV is set by:
- d-band hybridization effects (not capturable by classical FF)
- Charge transfer / covalent bonding at short range
- Many-body dispersion beyond pairwise C₆

For a grid-based force field that enables real-time molecular scanning,
15-20 meV RMSE would be excellent — comparable to DFT-D3 accuracy.
