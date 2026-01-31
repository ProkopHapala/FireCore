## Rigid local Procrustes (ports → neighbors)

Goal: for each atom (center) with local ports $r_j$ (body frame) and world neighbors $p_j$, find position $t$ and rotation $R$ (or quaternion $q$) minimizing

$$E = \sum_j k_j \| t + R r_j - p_j \|^2$$

This is the ARAP local step and the XPBD port-matching constraint. Translation and rotation decouple: translation is determined by centroids once $R$ is known.

### Physical motivation
- Each port-neighbor pair is a stiff spring (weight/stiffness $k_j$) enforcing relative orientation/spacing of a rigid motif.
- Minimizing $E$ is equivalent to zeroing net torque/force from those springs in the local frame. The optimal $R$ aligns the weighted covariance of $(r_j, p_j)$; the optimal $t$ aligns weighted centroids.

### Problem formulation (3D)
1) Weighted centroids
$$c_r = \frac{\sum k_j r_j}{\sum k_j}, \quad c_p = \frac{\sum k_j p_j}{\sum k_j}$$

2) Centered vectors
$$r'_j = r_j - c_r, \quad p'_j = p_j - c_p$$

3) Covariance (3×3)
$$H = \sum_j k_j\, p'_j (r'_j)^T$$

4) Optimal rotation is the closest orthogonal matrix to $H$ (polar factor). Position follows:
$$t = c_p - R c_r$$

### Rotation solvers

#### A) Newton–Schulz polar (matrix, used in solve_shape_matching_3d)
- Initialize $R_0 = H$.
- Iterate (3–5 steps):
$$R_{k+1} = \tfrac{1}{2} R_k (3I - R_k^T R_k)$$
- Converges to the orthogonal factor of $H$ when $H$ is nonsingular and close to rotation. Pure matmul/add; GPU friendly. Risk: if $H$ is near singular/degenerate, iterations can drift or produce reflection; column renorm/det fix may be needed in edge cases.

#### B) Horn quaternion via K matrix + power iteration (used in solve_optimal_transform_3d)
- Build scalar covariance components $s_{ab}$ of $H$ and construct symmetric 4×4 K matrix:
  - $\text{tr} = s_{xx}+s_{yy}+s_{zz}$
  - diagonal: $(k_{00},k_{11},k_{22},k_{33}) = (\text{tr}, s_{xx}-s_{yy}-s_{zz}, s_{yy}-s_{xx}-s_{zz}, s_{zz}-s_{xx}-s_{yy})$
  - off-diagonal uses $(s_{yz}-s_{zy}, s_{zx}-s_{xz}, s_{xy}-s_{yx})$ and symmetric sums $(s_{xy}+s_{yx}, s_{zx}+s_{xz}, s_{yz}+s_{zy})$.
- Dominant eigenvector of K is the optimal unit quaternion $q$ (maximizes alignment). Solve via a few power iterations with normalize each step, warm-started from previous frame.
- Translation: $t = c_p - R(q) c_r$.
- Pros: inherently enforces det=+1 and numerically robust for small neighbor counts; warm-start exploits temporal coherence. Cons: still sensitive if all weights ~0 or points colinear/coplanar.

### Edge cases / stability
- If $\sum k_j \approx 0$: no neighbors → leave state unchanged or write zeros; kernels currently early-return with minimal guard.
- Degenerate geometry (colinear/coplanar) makes $H$ near singular: Newton–Schulz may diverge; quaternion method tends to be more stable but eigenvector can flip sign; downstream should normalize and, if needed, enforce det=+1 or reflection fix.
- Neighbor slots: unused slots are marked j<0 with k=0; ensure stiffness matches the mask to avoid garbage torque.
- Warm starts: quaternions should be normalized before iteration to avoid blow-up.

### Summary of GPU kernels (XPDB_new.cl)
- `solve_shape_matching_3d`: builds centroids, covariance H, runs 5 Newton–Schulz iterations, outputs translation and rotation matrix rows.
- `solve_optimal_transform_3d`: builds centroids and covariance scalars, constructs Horn K, runs 4 power iterations on quaternion (warm start), outputs translation and updated quaternion.

## Test plan (single-node, deterministic)
Objective: verify recovered translation/rotation from both kernels on a minimal rigid motif before using in full MD/XPBD.

1) Geometry: load CH4 (@cpp/common_resources/xyz/CH4.xyz) via `load_xyz`; center C is the node, H atoms are neighbors (4 ports). Ports = H_i - C in rest frame. Stiffness all 1.
2) Ground truth transform: choose known rotation (e.g., axis-angle about z) and translation; generate target neighbor positions: $p_j^{true} = t_{true} + R_{true} r_j$.
3) Perturb initial state: use `perturb_state(pos, quat, perturb_pos, perturb_rot, rng)` to offset both center and neighbors plus initial quaternion.
4) Device data: build `neighs` (node→H indices, fill -1 after 4), `port_local`, `stiffness`, `curr_quat` (identity or perturbed), `pos` with perturbed neighbor coords.
5) Run `solve_optimal_transform_3d` once; read back `target_pos` and `curr_quat`. Compute errors:
   - translation error $\|t_{est}-t_{true}\|$
   - rotation error via rotated ports RMS against $p_j^{true}-t_{true}$.
6) Run `solve_shape_matching_3d`; reconstruct R from rows; same error metrics. Check det(R)>0.
7) Tolerances: expect <1e-4–1e-3 for small perturbations. Fail loudly if non-finite or norms explode.
8) Edge-case note: if any slot weight=0 or neighbor masked, ensure arrays consistent to avoid false torque.

Once this passes, integrate as a new `--method` in `test_rigid_XPBD_molecules.py` that calls the chosen kernel inside the MD loop.
