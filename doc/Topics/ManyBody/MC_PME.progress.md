# MC + Pauli Master Equation Solver: Implementation Progress Report

## 1. General Formulation of the Problem

### 1.1 Physical System
The system consists of **N quantum dots** (sites) that can be occupied by electrons. The state of the system is defined by an occupancy vector:
```
|n⟩ = |n₁, n₂, ..., n_N⟩  where nᵢ ∈ {0, 1}
```

### 1.2 Hamiltonian
The many-body Hamiltonian is:
```
H_sys(n) = Σᵢ εᵢ nᵢ + Σᵢ<ⱼ Wᵢⱼ nᵢ nⱼ
```
- **εᵢ**: On-site energy (includes intrinsic level + electrostatic potential from STM Tip)
- **Wᵢⱼ**: Coulomb interaction energy between electrons at sites i and j

### 1.3 Tunneling & Rates (Fermi's Golden Rule)
Transitions occur via **single-electron tunneling** between a site and one of the leads (Tip or Substrate).

**For adding an electron** (Tunneling IN):
```
W_{a→b}^L = (1/ħ) · Γ_{L,i} · f(ΔE - μ_L)
```
where f(E) = 1/(1 + e^{E/k_B T}) is the Fermi-Dirac distribution.

**For removing an electron** (Tunneling OUT):
```
W_{a→b}^L = (1/ħ) · Γ_{L,i} · [1 - f(ΔE - μ_L)]
```

**Chemical Potentials:**
- μ_substrate = 0 (ground reference)
- μ_tip = e · V_bias (shifted by bias voltage)

**Tunneling Matrix Elements:**
- Substrate: Assumed constant coupling Γ_sub
- Tip: Depends exponentially on tip position:
  ```
  Γ_tip,i(r_tip) = C_pre · e^{-2β|r_tip - r_i|}
  ```
  (In code, T_factor = e^{-βr}, so rate ∝ T_factor²)

### 1.4 Pauli Master Equation (PME)
In steady state, the probability distribution P solves:
```
dP/dt = W · P = 0
```
where W is the rate matrix (generator matrix). The solution gives stationary probabilities from which the **current** is computed:
```
I = Σ_{a,b} (P_a · W_tip(a→b) - P_b · W_tip(b→a))
```

### 1.5 The Challenge: Exponential Scaling
- **Full PME**: Number of many-body states scales as 2^N (128 sites → 2^128 states, impossible)
- **Monte Carlo**: Can handle 128 sites easily (ground state search)
- **Goal**: Combine MC ground state with reduced-basis PME solver

### 1.6 Subspace Selection Strategy
The key insight is to select only **kinetically relevant states** based on:
1. **Single-flip excitations** from ground state (Hamming distance = 1)
2. **Kinetic score**: P_i^est ≈ Γ_in / (Γ_in + Γ_out) — captures non-equilibrium physics better than energy alone
3. **Double excitations**: Formed from pairs of highly active single excitations

## 2. Implementation Architecture

### 2.1 Two-Solver Approach

**Full Solver (PME.cl)**:
- Enumerates ALL 2^N states (for N ≤ 6-8)
- Solves full rate matrix via Gauss-Jordan elimination
- Used as reference for validation

**Hubbard Solver (hubbard.cl)**:
- Monte Carlo for ground state
- Kinetic basis scanner for subspace selection
- Dense PME solver on reduced basis (max 64 states)
- Scalable to large N

### 2.2 Pipeline
```
1. MC Ground State Search (solve_mc_2phase)
       ↓
2. Precalc Esite, Tsite (precalc_Esite_Thop kernel)
       ↓
3. Kinetic Basis Scanner (kinetic_basis_scanner kernel)
       ↓
4. Dense PME Solver (solve_pme_dense_batch kernel)
       ↓
5. Current Calculation
```

## 3. Challenges Faced and Solutions

### Challenge 1: NaN in Esite Calculation
**Symptom**: `max_dI=nan` in test results, Esite values were NaN

**Root Cause**: In `HubbardSolver.py:145`, the `realloc_precalc_buffers` function allocated:
```python
sz_f*4*nSingle*nMulti  # WRONG: float4 size
```
for multipole coefficients that are plain floats. This oversized buffer caused subsequent buffer reuse issues that zeroed out the rotation matrix.

**Solution**: Fixed to:
```python
sz_f*nSingle*nMulti  # CORRECT: plain float size
```

### Challenge 2: Singular Rate Matrix for Embedded Systems
**Symptom**: `curr=[nan]`, `probs=[0,0,...,0]` from dense solver

**Root Cause**: For 2-site embedded in 4-site system, generating all 2^4=16 basis states included spectator sites (E=100 eV), producing a singular rate matrix.

**Solution**: Modified `make_basis_u4` to generate only states where active sites vary:
```python
def make_basis_u4(nTips, nSite, active_sites=None):
    # Only active sites vary; spectator sites stay 0
```

### Challenge 3: NaN at Low Bias Voltages
**Symptom**: `max_dI=nan` for 4-site V-scan at very low bias (V ≈ 0 to 0.04)

**Root Cause**: Gauss-Jordan solver in `hubbard.cl` lacked **partial pivoting**, causing numerical instability for nearly singular matrices at equilibrium.

**Solution**: Added partial pivoting matching PME.cl implementation:
```c
// Find best pivot row
int pivot_row = k;
float max_val = fabs(Mat[k][k]);
for (int i = k + 1; i < N; i++) {
    float val = fabs(Mat[i][k]);
    if (val > max_val) { max_val = val; pivot_row = i; }
}
// Swap rows if needed
```

### Challenge 4: Rate Matrix dE Sign Bug (CRITICAL)
**Symptom**: `max_dP ≈ 1.0`, `max_dI ≈ 3e-6` (~100% error) despite perfect Esite/Thop parity

**Root Cause**: In rate matrix construction, for **removing transitions** (c→r where r has fewer electrons):
- `dE = E_r - E_c = -ε_add` (where ε_add = E(more) - E(fewer))
- Passing this to `calc_total_rate(..., false)` gave:
  ```
  Γ·(1-f(-ε_add)) = Γ·f(ε_add)  // WRONG: this is the ADDING rate!
  ```
- Should have been: `Γ·(1-f(ε_add))`

**Solution**: One-line fix:
```c
float dE_add = adding ? dE : -dE;  // Always use ε_add = E(more) - E(fewer)
```

## 4. Comprehensive Results Table

| Test | max\|dE\| | max\|dT\| | max\|dP\| | max\|dI\| | I_spread |
|------|-----------|-----------|-----------|-----------|----------|
| 2-site x-scan | 0 | 4.5e-8 | 1.2e-6 | **1.6e-12** | [2.8e-7, 3.0e-6] |
| 2-site V-scan | 0 | 0 | 3.4e-6 | **5.7e-14** | [0, 5.6e-7] |
| 2-site 2D x-V | — | — | — | **2.1e-11** | — |
| 4-site x-scan | 0 | 4.5e-8 | 2.0e-6 | **2.3e-12** | [1.3e-6, 3.2e-6] |
| 4-site V-scan | 0 | 0 | 3.4e-6 | **1.9e-12** | [~0, 1.8e-6] |
| 4-site 2D x-V | — | — | — | **5.5e-11** | — |

**Interpretation**: All errors are at **float32 rounding level** (1e-6 to 1e-12). Current spreads are **nontrivial** (order 1e-6), confirming real physics is being computed.

## 5. Input Parameters

### 5.1 Physical Parameters Used in Tests

| Parameter | Value | Description |
|-----------|-------|-------------|
| Rtip | 3.0 | Characteristic tip distance |
| zV0 (zMirror) | -0.5 | Mirror plane offset |
| zVd (zRampOffset) | 0.0 | Ramp offset |
| beta (Thop_decay) | 0.5 | Hopping decay constant |
| Gamma | 1.0 | Lead coupling |
| bMirror | 1.0 | Mirror term enabled |
| bRamp | 1.0 | Ramp term enabled |
| W_scalar | 0.0 | Inter-site Coulomb (off for parity) |
| T_K | 4.0 | Temperature (Kelvin) |
| zTip | 3.0 | Tip height |

### 5.2 Multipole Coefficients
- **Monopole (Q0)**: 1.0
- **Quadrupole (Qzz)**: 0.0 (for parity tests)
- **Order**: 0 (monopole only)

### 5.3 Scan Ranges

**X-scan (tip position)**:
- Range: x ∈ [-2.0, 7.0] (2-site), [-10.0, 10.0] (4-site, original)
- Fixed bias: V = 0.50 V
- Points: 91 (simplified) or 201 (original)

**V-scan (bias voltage)**:
- Range: V ∈ [0.0, 1.0] (simplified) or [0.0, 1.5] (original)
- Fixed x: 2.5 (2-site), 0.0 (4-site)
- Points: 101 (simplified) or 151 (original)

**2D x-V scan**:
- X: 101 points
- V: 76 points

## 6. Output Directories

### 6.1 Original Parity Tests (with Ruslan geometry)
```
/home/prokop/git/FireCore/figs/parity_2site_ruslan/
├── xscan/
│   ├── scan_2site_xscan.csv
│   ├── plot_2site_xscan.png
│   ├── summary_2site_xscan.txt
│   └── errors_2site_xscan.json
├── Vscan/
│   ├── scan_2site_Vscan.csv
│   ├── plot_2site_Vscan.png
│   ├── summary_2site_Vscan.txt
│   └── errors_2site_Vscan.json
└── xVmap/
    ├── I_full_2site_xV.csv
    ├── I_hub_2site_xV.csv
    ├── plot2d_2site_xV.png
    └── summary_2site_xV.txt

/home/prokop/git/FireCore/figs/parity_4site_ruslan/
├── xscan/
├── Vscan/
└── xVmap/
```

### 6.2 Simplified Parity Tests (newer version)
```
/home/prokop/git/FireCore/figs/parity_2site/
├── xscan/
└── Vscan/

/home/prokop/git/FireCore/figs/parity_4site/
├── xscan/
└── Vscan/
```

## 7. Relevant Files

### 7.1 OpenCL Kernels
| File | Description |
|------|-------------|
| `/home/prokop/git/FireCore/pyBall/OCL/cl/hubbard.cl` | Hubbard/Ising model kernels: MC, precalc, kinetic scanner, dense PME solver |
| `/home/prokop/git/FireCore/pyBall/OCL/cl/PME.cl` | Full PME solver (reference, for small N) |

### 7.2 Python Modules
| File | Description |
|------|-------------|
| `/home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py` | Python harness for Hubbard solver |
| `/home/prokop/git/FireCore/pyBall/OCL/pauli_ocl.py` | Python harness for full PME solver |
| `/home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_runall.py` | Main parity test script (original) |
| `/home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_single.py` | Single-point parity test |
| `/home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_sweep.py` | 1D sweep parity test |

### 7.3 Geometry Files
| File | Description |
|------|-------------|
| `/home/prokop/git/FireCore/pyBall/OCL/Ruslan_long.txt` | 2-site geometry |
| `/home/prokop/git/FireCore/pyBall/OCL/Ruslan_kite.txt` | 4-site kite geometry |

## 8. Summary

### What Was Achieved
1. **Implemented** modular MC + PME pipeline with kinetic subspace selection
2. **Achieved near-machine-precision parity** (max|dI| ≈ 1e-12) between full and Hubbard solvers
3. **Fixed 4 critical bugs**:
   - Buffer size in precalc
   - Active-sites basis generation
   - Partial pivoting in Gauss-Jordan
   - Rate matrix dE sign convention

### Current Status
- ✅ Esite/Thop precalc: **Perfect parity** (max|dE|=0)
- ✅ Probabilities: **Excellent parity** (max|dP| ≈ 3e-6)
- ✅ Currents: **Excellent parity** (max|dI| ≈ 1e-12)
- ✅ All 6 test configurations pass

### Next Steps (Future Work)
1. Implement analytic "Star" solver for cases with no lateral hopping
2. Optimize kinetic basis scanner (parallel sort instead of serial)
3. Test with finite W (inter-site Coulomb)
4. Scale to larger systems (N > 8) with kinetic subspace selection
5. Add per-site current decomposition

*Report generated: February 2026*
*Project: FireCore - Monte Carlo + Pauli Master Equation Solver*