# Fast Particle Collisions by Projective dynamics

This document outlines the detailed implementation plan for integrating fast particle collisions using projective dynamics. The key innovation is splitting particle interactions into:
1. Linear (harmonic) part - handled by iterative linear solver (e.g. Jacobi, Gauss-Seidel)
2. Non-linear smoothening part - handled as external force

## Background

[Position based dynamics](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf) and [projective dynamics](https://www.projectivedynamics.org/Projective_Dynamics/index.html) are very efficient methods for solving hard-degrees of freedom in molecules  like bonds ( i.e. stiff ineraction between particles with fixed topology.

Advantage of these schemes is that they are unconditionaly stable, which allows to significantly increase time-step of molecule (or other dynamical) simulations, since the hard dregrees of freedom are treated by implicit solver, rather than explicitly propagated by dynamics through force. We avoid limitations of traditional explicit force-based dynamics simulation where the time step, is generally inversly proportional to highest vibration frequency $f \approx \sqrt{k/m}$ by stiffnes of these hard degrees for freedom. The limitation of implicit schemes is, on the other hand, that they are applicable only to linear constrains. Therefore additional non-linear (anharmonic) forces, or non-conservative forces (e.g. air drag) must be added externlly (as external force). These forces than determine the required time step.  

Here we try to extend the method to solve collisions between the particles (i.e. stiff short range interactions bettween atoms).

## Problems:

1. The interactions are not-linear (or harmonic) unlike bonds, but are anharmonic and dissociative (i.e. like Leanrd-Jones). This means we can use linear solver only approximatively arround the minimum.

2. Due to rapid changes of the interaction pairs, it does not make much sense to explicitly construct matrixes (like in Choolesky factorization) as these matrixes would rapidly change with rapidly changing interaction topology. Therefore it makes much more sense to use iterative methods, like Jacobi or Gauss-Seidel.

## Split Quasi Linear Poential

To efficiently use linear solver like Jacobi or Gauss-Seidel we will split our anharmonic potential (like Lenard-Jones) to two parts.

1. Linear (harmonic) part
2. non-linear dissociative part

Within the framework of projective dynamics we can than implement the harmonic part in the linear solver and the non-linear part as external force. The importaint assumption is that the non-linear part is small and soft (much smaller stiffness) than the harmonic part.

## The Linear poential cutoff

The linear solver tries to find such position of each particle which minimize the potential energy form all linear (harmonic) constrains. In order to function properly the linear constrains must be defined in some reasonable range around the minimum (where force from this constrain is zero). although for collision is most importaint just for the repulsive part, for numerical stability and convenience we typically extend the linear region slightly beyond the minimum to some $R_{cut} > R_{min}$. Therefore the particles-particle force can be expressed as $F=k(r_{ij}-R_{min})$ with constant stiffness $k$ along the whole interval of distance $r<R_{cut}$.

This Linear harmonic potential is implemented into the linear solver (Jacobi or Gauss seidel).
 - At every time step, we need to update the neighbor list. The neighbor list for contrain solver should contain (i) all bonded neighbors (which are constant over the simulation, as the bond topology is fixed) (ii) all collision pairs of particles within $r<R_{cut}$.  
 - If the particles are beyond $r<R_{cut}$ linear solver ignores them.

## The smoothening potential

The linear poential has serious problem: I is swithing abruptly between zero force for $r>R_{cut}$ and non-zero force $f=-k(R_{cut}-R_{min})$ when the particles cross $R_{cut}$ boundary. **This would lead to unstable simulations**, dificult convergence and non-physical behaviour.

For this reason we need to implement additional anharmonic disociative correction poential which we call **smoothening potential**. This poential should be defined for $R_{min}<r<R_{cut2}$ and which transition smoothly between the harmonic one, and zero  force and energy (E=0,F=0) at distance $r>R_{cut2}$. Beyond $R_{cut2}$ the particles should not interact by any means. That means that we can use fast spatial-partitioning algorithms for solving non-covalent interactions betwee the particles.

There are many possible functional forms for the smoothening potential. Here we will use perhaps the simplest one, which is antoher quadratic potential but with negative stiffness. The functional form is as follows:

```c++
inline double getSR_x2_smooth(const Vec3d& dp, Vec3d& f, double E_min, double R_min, double R_cut, double R_cut2) {
    double r2 = dp.norm2();
    if( r2 > (R_cut2*R_cut2) ){
        f = Vec3dZero;
        return 0.0;
    }
    double r  = sqrt(r2);
    /// TODO: we may precompute these constants  k1,k2 (and d1,d2)
    double d1 = R_cut  - R_min;  // harmonic region:  [R_min, R_cut]
    double d2 = R_cut2 - R_cut;  // smoothing region: [R_cut, R_cut2]  
    double k1 = -2.0 * E_min / ( d1 * (d1 + d2) );
    double k2 = - k1 * (d1 / d2);
    double E, F;            
    if (r < R_cut) {
        double x = r - R_min;
        E = 0.5 * k1 * x * x + E_min;
        F =     - k1 * x;
    } else {
        double x  = r - R_cut2;
        E = 0.5 * k2 * x * x;
        F =     - k2 * x;
    }
    f.set_mul(dp, F/r );
    return E;
}
```

```c++
Vec2d computeCollisionParameters(double E_min, double R_min, double R_cut, double R_cut2) {
    double d1 = R_cut - R_min;
    double d2 = R_cut2 - R_cut;
    double k1 = -2.0 * E_min / (d1 * (d1 + d2));
    double k2 = -k1 * (d1 / d2);
    return Vec2d{k1, k2};
}
```


# Implementation plan

## Key Components and Files

### Core Files:
- `MolWorld_sp3.h`: Main simulation loop and integration
- `ProjectiveDynamics_d.h` and `ProjectiveDynamics_d.cpp`: Linear solver implementation in `run_LinSolve()`
   - the iterative linear solvers for constrains like `updateIterativeJacobi()` and `updateIterativeMomentum()` and `updatePD_dRHS()`use variables like 
     - `ProjectiveDynamics_d::neighs` - [natom*nNeighMax] array of neighbor atom indices
     - `ProjectiveDynamics_d::params` - [natom*nNeighMax] params=0; neighbor distance constrain parameters (l0,kP,kT,damp), we only need l0 and kP
- `NBFF.h`: Neighbor list management and non-linear force evaluation by modified `evalSortRange_BBs()`
   - `NBFF:neighs` is neighborlist for fixed bond topology, these should be copied to `ProjectiveDynamics_d::neighs`
- `Forces.h`: 
   -  `getSR_x2_smooth()` - the short-range potential and forces for particle collisions using quadratic potential for both harmmonic part and smoothing part.
   - `_mixREQ` is used for mixing non-covalent parameters `REQ = {RvdW,EvdW,Q,Hbond}`.
- `MMFFsp3_loc.h`: Force field implementation, defines variables `bLs` and `bKs` storing bond lengths and bond stiffnesses. These parameters should be copied to `ProjectiveDynamics_d::params` for fixed bond topology.

### Key Functions:
- `getSR_x2_smooth()`: Calculates short-range potential and forces for particle collisions
- `NBFF::evalSortRange_BBs()`: Use axis-aligned bounding boxes to efficiently find neighbors and evaluate collision interactions between particles when $r_{ij}<R_{cut2}$. We will modify this function to update neighbor list for Projective Dynamics for $r_{ij}<R_{cut}$ and add non-linear external force for $R_{cut}<r_{ij}<R_{cut2}$.
   - we need to make sure that for added neighbors originating from the collisions, we set the params l0,kP properly based on non-covalent potential parameter l0=R_min, E_min, kP = k1 from `getSR_x2_smooth()`
   - R_min and E_min for `getSR_x2_smooth()` are obtained from `NBFF::REQs` as `REQ_ij=_mixREQ(REQ[i],REQ[j]); R_min=REQ_ij.x; E_min=REQ_ij.y;`
- `ProjectiveDynamics_d::run_LinSolve()`: solves the dynamics of the system by splitting the interaction into linear constrain solver and non-linear external force part
  - within its body `run_LinSolve()` calls function pointer `user_update`. That is where we can call `NBFF::evalSortRange_BBs()` to update neighbor list and add non-linear external force.
- Alternatively, we may define a new function `MMFFsp3_loc::run_projective_dynamics()` which will call `NBFF::evalSortRange_BBs()` and re-implements relevant parts of `ProjectiveDynamics_d::run_LinSolve()` usin lower-level functions from `ProjectiveDynamics_d` class such as `ProjectiveDynamics_d::updateIterativeMomentum()` or `ProjectiveDynamics_d::updateIterativeJacobi()`

## Implementation Steps

This plan details integrating fast particle collisions into a molecular dynamics simulation using projective dynamics. The approach modifies the `NBFF` class to handle collision detection and force calculations efficiently, integrating with the existing `ProjectiveDynamics_d` class, specifically the `updatePD_dRHS` function,  without altering `prepare_LinearSystem` or `rhs_ProjectiveDynamics`.  `MMFFsp3_loc` will be adapted to use the new collision handling.

**I. Data Structures and Initialization:**

1.  **Modify `MolWorld_sp3`:** Add member variables to `MolWorld_sp3` to store collision parameters:
    *   `double R_cut2`: Outer cutoff distance for the smoothing potential.  This parameter will control the range of the collision detection.
2.  **Modify `NBFF`:**  No changes are needed to the data structures in `NBFF`.  We will work with the existing member variables.
3.  **Modify `ProjectiveDynamics_d`:** No changes to the data structures are required. We will leverage the existing `neighs` and `params` arrays.

**II. Collision Detection and Force Calculation (`NBFF.h`):**

1.  **write `evalSortRange_BBs_PD()` based on `NBFF::evalSortRange_BBs()`:** 
    * This function will be heavily modified to handle collision detection and force calculation efficiently.
    *   The function should take as output arguments  pointers  `pd.neighs` and `pd.params` from  `ProjectiveDynamics_d& pd`
    *   Use the existing bounding box structure (`BBs`, `pointBBs`) to efficiently identify atom pairs within `R_cut2`.
    *   For each atom `i` list all atoms j within `r_{ij} < R_cut2`:
        *  Obtain `R_min` and `E_min` from `_mixREQ(REQs[i], REQs[j])`
        *  Calculate `k1` and `k2` parameters from `getSR_x2_smooth` using `computeCollisionParameters()`
        *  If `r < R_cut`:  Add `(i, j)` to `pd.neighs`. Set the corresponding entry in `pd.params` to `l0 = R_min` and `kP = k1`.
        *  If `R_cut ≤ r < R_cut2`: Apply the non-linear smoothing force `f` directly to `pd.forces` 
           * Note: this is basically the external force which is used inside `ProjectiveDynamics_d::run_LinSolve()`, 
              * see line : `ps_pred[i] = points[i].f + vel[i].f*dt + forces[i].f*(dt2/points[i].w); `


**III. Projective Dynamics Solver Integration (`ProjectiveDynamics_d.cpp`):**

The key to integrating the collision handling efficiently is to modify the existing `updatePD_dRHS` function.  We will avoid modifying `rhs_ProjectiveDynamics` and will not remove `prepare_LinearSystem` as it may be used for other purposes.

*  ** `updatePD_dRHS()`:** This function is the core of the iterative solver. It updates the right-hand side (RHS) vector (`bvec`) based on the current positions (`pnew`) at each time step before we start the iterative solver loop. This function **DOES NOT NEED TO BE MODIFIED**, but we should check if we properly updated all variables it uses, namely `pd.neighs` and `pd.params`. These should contain both the fixed bonded neighbors (obtained from `NBFF:neights` and `bLs`, `bKs`) as well as dynamically updated collision neighbors (obtained from `NBFF:evalSortRange_BBs_PD`).

**IV.  `MMFFsp3_loc` Integration:**

1.  **Call `NBFF::evalSortRange_BBs_PD()`:**  Modify the `eval_atoms()` function in `MMFFsp3_loc` to call the updated `NBFF::evalSortRange_BBs_PD()` function *before* any other force calculations within the main simulation loop. This will handle both bonded and non-bonded forces, including the pre-computation of non-linear collision forces and populating the `ProjectiveDynamics_d` neighbor list.  The non-linear forces will be applied directly to atom positions, so you should remove any redundant collision force calculations.

**V. `MolWorld_sp3` Integration:**

1.  **Add Parameters:** Add `R_cut2` and `alphaMorse` (if needed by `getSR_x2_smooth`) as member variables to `MolWorld_sp3` to store collision parameters. These parameters should be configurable (command-line arguments, configuration file).

2.  **Create `run_projective_dynamics()`:** Create a new function `MolWorld_sp3::run_projective_dynamics()` (or modify existing MD loop functions) to integrate the new collision handling. This function should:
    *   Initialize collision parameters (`R_cut2`, `alphaMorse`).
    *   Call `ffl.evalSortRange_BBs_PD(pd, R_cut2, alphaMorse)` within the main loop *before* each call to `pd.run_LinSolve()`.  This ensures the neighbor list and parameters are up-to-date for the projective dynamics solver.

**VI. Code Snippets (Illustrative):**

**A.  `NBFF::evalSortRange_BBs_PD()` (Partial):**

```c++
double NBFF::evalSortRange_BBs_PD(ProjectiveDynamics_d& pd, double R_cut2, double alphaMorse) {
    // ... existing code for bonded interactions ...

    for (int ib = 0; ib < nBBs; ib++) {
        // ... existing code ...

        int n_collisions = selectFromOtherBucketsInBox(ib, R_cut2, ...);
        for (int k = 0; k < n_collisions; k++) {
            int i = ...; // Index of atom i
            int j = ...; // Index of atom j
            Vec3d dp = apos[j] - apos[i];
            double r = dp.norm();
            double R_min_ij = _mixREQ(REQs[i], REQs[j]).x;
            Vec3d f;
            double E = getSR_x2_smooth(dp, f, R_min_ij, /*R_cut calculated inside*/, R_cut2, alphaMorse);
            if (r < R_cut2) {
                if (r < R_cut) {
                    pd.addNeighbor(i, j, R_min_ij, k1); // Add to projective dynamics neighbor list
                } else {
                    apos[i].add(f);
                    apos[j].sub(f);
                }
            }
            E += E_ij;
        }
    }
    // ... existing code ...
    return E;
}
