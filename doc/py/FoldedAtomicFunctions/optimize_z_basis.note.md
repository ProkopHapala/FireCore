# Optimal Analytical Basis along Z-Direction

This document describes a method and its implementation in `optimize_z_basis.py` for finding an optimal, compact set of analytical basis functions to represent a collection of 1D potential profiles, particularly those varying along the z-direction (e.g., above a material surface).

## 1. Motivation and Goals

When dealing with numerous 1D potential profiles, such as those generated from varying system parameters (atomic charges, vdW parameters, lattice constants), it's often desirable to:

*   **Efficiently Describe Data:** Find a compact representation for these profiles rather than storing all raw data points.
*   **Reduce Dimensionality:** Identify the most significant features or shapes present across the dataset.
*   **Obtain Analytical Forms:** If the basis functions are analytical (e.g., polynomials or combinations thereof), the potential can be represented by a few coefficients, facilitating further analysis or use in other models.
*   **Physical Insight:** The shapes of the optimal basis functions might offer insights into the underlying physics governing the potential variations.

This method aims to derive such an optimal basis set using Singular Value Decomposition (SVD).

## 2. Theoretical Background

The core idea is to represent a given set of sample functions using a small number of basis functions that are themselves linear combinations of a larger, predefined "library" of analytical functions.

### 2.1. Problem Statement

We have a set of $M$ sample functions, $y_j(z)$, for $j=1, \dots, M$. Each function is sampled at $N_z$ discrete points along the z-axis, $z_i$, for $i=1, \dots, N_z$.
We can represent these sample functions as a matrix $\mathbf{Y}_T$ of shape $(M \times N_z)$, where each row is a sample function $y_j(z_i)$.

### 2.2. Library Basis Representation

First, we choose a "library" of $P$ analytical basis functions, $\phi_p(z)$, for $p=1, \dots, P$. These can be polynomials, exponentials, Gaussian functions, or any other suitable set. These library functions, evaluated at the $N_z$ points, form a matrix $\mathbf{\Phi}_T$ of shape $(P \times N_z)$, where each row is a library basis function $\phi_p(z_i)$.

Each sample function $y_j(z)$ can be approximated as a linear combination of these library basis functions:
$$ y_j(z) \approx \sum_{p=1}^{P} S_{pj} \phi_p(z) $$
In matrix notation, considering all sample functions and their evaluations at $z_i$:
Let $\mathbf{Y}$ be an $(N_z \times M)$ matrix where columns are sample functions $y_j(z_i)$.
Let $\mathbf{\Phi}$ be an $(N_z \times P)$ matrix where columns are library functions $\phi_p(z_i)$.
Let $\mathbf{S}$ be a $(P \times M)$ matrix of coefficients.
Then, the approximation is $\mathbf{Y} \approx \mathbf{\Phi} \mathbf{S}$.

The script uses a transposed convention for convenience with row-major NumPy arrays:
*   $\mathbf{Y}_T$: $(M \times N_z)$ - sample functions (input `Y_T`)
*   $\mathbf{\Phi}_T$: $(P \times N_z)$ - library basis functions (input `Phi_T`)
*   $\mathbf{S}$: $(P \times M)$ - coefficient matrix (output `S_coeffs`)

The relationship is then $\mathbf{Y}_T^T \approx \mathbf{\Phi}_T^T \mathbf{S}$.
The coefficient matrix $\mathbf{S}$ is found by solving this system, typically using (weighted) least squares:
$$ \mathbf{S} = (\mathbf{\Phi}_T^T)^\dagger \mathbf{Y}_T^T $$
where $\dagger$ denotes the Moore-Penrose pseudoinverse. If weights $w(z_i)$ are used, the problem becomes minimizing $\sum_i (y_j(z_i) - \sum_p S_{pj} \phi_p(z_i))^2 w_i^2$.

### 2.3. Optimal Basis via Singular Value Decomposition (SVD)

The coefficient matrix $\mathbf{S}$ ($P \times M$) contains the representation of each of the $M$ sample functions in the $P$-dimensional space spanned by the library basis. To find a more compact representation, we perform SVD on $\mathbf{S}$:
$$ \mathbf{S} = \mathbf{U}_{\text{svd}} \mathbf{\Sigma}_{\text{svd}} \mathbf{V}_{\text{svd}}^T $$
*   $\mathbf{U}_{\text{svd}}$ ($P \times P$ or $P \times \min(P,M)$) is an orthogonal matrix whose columns are the left singular vectors. These vectors represent principal directions (combinations of library basis functions) in the coefficient space.
*   $\mathbf{\Sigma}_{\text{svd}}$ is a diagonal matrix containing the singular values, which indicate the importance of each principal direction.
*   $\mathbf{V}_{\text{svd}}^T$ ($M \times M$ or $\min(P,M) \times M$) is an orthogonal matrix whose rows are the right singular vectors.

We select the top $K$ left singular vectors (columns of $\mathbf{U}_{\text{svd}}$ corresponding to the $K$ largest singular values) to form a matrix $\mathbf{U}_k$ of shape $(P \times K)$. These $K$ vectors are the coefficients that define the $K$ new optimal basis functions in terms of the original library basis $\mathbf{\Phi}_T$.

### 2.4. Constructing the Optimal Basis

The $k$-th new optimal basis function, $B_{\text{opt},k}(z)$, is a linear combination of the original library functions $\phi_p(z)$:
$$ B_{\text{opt},k}(z) = \sum_{p=1}^{P} (U_k)_{pk} \phi_p(z) $$
where $(U_k)_{pk}$ is the element in the $p$-th row and $k$-th column of $\mathbf{U}_k$.

Evaluated on the grid $z_i$, the set of $K$ optimal basis functions forms a matrix $\mathbf{B}_{\text{opt},T}$ of shape $(K \times N_z)$:
$$ \mathbf{B}_{\text{opt},T} = \mathbf{U}_k^T \mathbf{\Phi}_T $$

## 3. Implementation Details and Workflow (`optimize_z_basis.py`)

The script `optimize_z_basis.py` implements this methodology through a series of modular functions orchestrated by `run_optimization_pipeline`.

### 3.1. Core Script
`optimize_z_basis.py`

### 3.2. Key Steps and Functions (Algorithm Outline)

1.  **Define Z-Coordinates (`common_z_coords`, `zs`):**
    *   A 1D NumPy array specifying the $N_z$ points along the z-axis where all functions are defined. This should respect any desired $z_{min}$.

2.  **(Optional) Define Z-Weights (`ws`):**
    *   A 1D NumPy array of $N_z$ weights.
    *   Generated by: `generate_z_weights(zs, z0, z_decay_width, min_weight)`
    *   Purpose: To emphasize or de-emphasize certain regions of $z$. For example, to dampen the highly repulsive part of potentials, a weighting function like $w(z) = 1.0$ for $z > z_0$ (attractive region) and $w(z) = 1.0 - ((z - z_0) / d)^2$ for $z \le z_0$ (repulsive region, clamped at `min_weight`) can be used.

3.  **Generate or Provide Sample Functions $\mathbf{Y}_T$ (variable `_Y_T`, shape $M \times N_z$):**
    *   **Option A (Generate):**
        *   Define system parameters: `define_system_configurations(num_systems, param_ranges, ...)` returns `system_definitions` (a list of dictionaries, each specifying atom parameters like $R_i, E_i, q_i$, and lattice constant $L_x$).
        *   Generate 1D potential profiles: `generate_1d_potential_profiles(system_definitions, common_z_coords, x_slice_coords, ...)` uses `FoldedAtomicFunctions.PotentialCalculator` to compute 2D potentials and extracts 1D slices along $z$ at specified `x_slice_coords`. These are interpolated onto `common_z_coords`.
    *   **Option B (Provide):** The user can directly supply the `Y_T` NumPy array.

4.  **Generate or Provide Library Basis $\mathbf{\Phi}_T$ (variable `_Phi_T`, shape $P \times N_z$):**
    *   **Option A (Generate):**
        *   `create_library_basis(zs, config)` returns `Phi_T`, `z_scale_info`, `basis_labels`.
        *   The `config` dictionary specifies the type:
            *   `{'type': 'polynomial', 'degree': d, 'scale_z': True}`: Generates polynomials $(z')^p$ up to degree $d$. If `scale_z` is true, $z$ is scaled to $[0,1]$ as $z' = (z - z_{\text{min}})/(z_{\text{max}} - z_{\text{min}})$ for numerical stability. `z_scale_info` stores $(z_{\text{min}}, z_{\text{range}})$.
            *   `{'type': 'custom', ...}` or `{'type': 'list_of_funcs', ...}`: Allows user-defined functions.
    *   **Option B (Provide):** The user can directly supply the `Phi_T` NumPy array and optionally `basis_labels` and `z_scale_info_provided`.

5.  **(Optional) Orthogonalize Library Basis:**
    *   If `orthogonalize_library` is `True`:
        `_Phi_T = orthogonalize_basis_gram_schmidt(_Phi_T, ws)`
    *   This makes the rows of `_Phi_T` orthonormal with respect to the weights `ws` (if provided). The `_basis_labels` are updated, and `_z_scale_info` is typically invalidated for simple polynomial interpretation.

6.  **Compute Coefficients $\mathbf{S}$ (variable `S_coeffs`, shape $P \times M$):**
    *   `S_coeffs = compute_coefficients_in_library(_Y_T, _Phi_T, ws, is_phi_orthogonal_weighted)`
    *   This function solves $\mathbf{Y}_T^T \approx \mathbf{\Phi}_T^T \mathbf{S}$ for $\mathbf{S}$.
    *   If `is_phi_orthogonal_weighted` is `True` (e.g., after Gram-Schmidt), coefficients are found by direct projection: $S_{pj} = \sum_i Y_{ji} \Phi_{pi} w_i$.
    *   Otherwise, (weighted) `scipy.linalg.lstsq` is used. For weighted least squares, the problem effectively transforms to minimizing $\|\mathbf{W}(\mathbf{Y}_T^T - \mathbf{\Phi}_T^T \mathbf{S})\|_F^2$, where $\mathbf{W}$ is a diagonal matrix with $\sqrt{w_i}$ (or $w_i$ depending on formulation, script uses $\sqrt{w_i}$ for `lstsq`).

7.  **Find Optimal Basis via SVD:**
    *   `B_opt_T, U_k_coeffs, s_vals_svd = find_optimal_analytical_basis_svd(S_coeffs, _Phi_T, num_optimal_K)`
    *   Performs SVD on `S_coeffs`: $\mathbf{S}_{\text{coeffs}} = \mathbf{U}_{\text{svd}} \mathbf{\Sigma}_{\text{svd}} \mathbf{V}_{\text{svd}}^T$.
    *   `U_k_coeffs` ($P \times K$) are the first $K$ columns of $\mathbf{U}_{\text{svd}}$.
    *   `B_opt_T` ($K \times N_z$) is the set of $K$ optimal basis functions evaluated on `zs`: $\mathbf{B}_{\text{opt},T} = \mathbf{U}_k^T \mathbf{\Phi}_T$.
    *   `s_vals_svd` are the singular values from $\mathbf{\Sigma}_{\text{svd}}$.

8.  **Analysis and Plotting:**
    *   `plot_1d_profiles`: Visualize sample functions, library basis, or optimal basis.
    *   `plot_singular_values`: Plot `s_vals_svd` to help choose $K$.
    *   `print_analytical_form_polynomial`: If the library basis was polynomial and not subsequently orthogonalized in a way that obscures its form, this prints the analytical expression of the optimal basis functions using `U_k_coeffs` and `z_scale_info`.

## 4. Key Variables and Their Meanings

*   `common_z_coords` (often aliased as `zs`): 1D NumPy array ($N_z$) of z-coordinates.
*   `Y_T`: 2D NumPy array ($M \times N_z$) of $M$ sample functions.
*   `Phi_T` (or `_Phi_T` in pipeline): 2D NumPy array ($P \times N_z$) of $P$ library basis functions.
*   `ws`: 1D NumPy array ($N_z$) of weights for z-coordinates.
*   `S_coeffs`: 2D NumPy array ($P \times M$), coefficients of sample functions in the library basis.
*   `num_optimal_K` (or `K`): Integer, desired number of optimal basis functions.
*   `U_k_coeffs`: 2D NumPy array ($P \times K$), coefficients defining the optimal basis functions in terms of the library basis $\mathbf{\Phi}_T$.
*   `B_opt_T`: 2D NumPy array ($K \times N_z$), the $K$ optimal basis functions evaluated on `common_z_coords`.
*   `s_vals_svd`: 1D NumPy array of singular values from the SVD of `S_coeffs`.
*   `z_scale_info`: Tuple `(z_min_orig, z_range_orig)` used if polynomial basis functions were generated on scaled z-coordinates.
*   `basis_labels`: List of strings describing each function in `Phi_T`.

## 5. Considerations

*   **Region of Interest ($z_{min}$):**
    It's crucial to define `common_z_coords` to start from a relevant $z_{min}$. Potentials at very small $z$ (highly repulsive region) can have extremely large values, which might numerically dominate the fitting process and overshadow the more subtle features in the attractive or equilibrium regions.

*   **Weighting Function ($w(z_i)$):**
    A weighting function allows for emphasizing certain parts of the z-range during the fitting process (specifically, when calculating `S_coeffs`). The goal is to minimize a weighted mean square error:
    $$ \text{Error} = \sum_j \sum_i \left( (y_j(z_i) - \sum_p S_{pj} \phi_p(z_i)) w(z_i) \right)^2 $$
    This is implemented in `compute_coefficients_in_library` by effectively solving for $S_{pj}$ in $y_j(z_i)w(z_i) \approx \sum_p S_{pj} (\phi_p(z_i)w(z_i))$. If Gram-Schmidt orthogonalization is used, the weights `ws` are incorporated into the definition of the inner product.

*   **Performance for Large Datasets:**
    *   The SVD is performed on the `S_coeffs` matrix ($P \times M$). The complexity of SVD is roughly $O(\min(P^2 M, M^2 P))$. If $P$ (number of library basis functions) is kept moderate, the method can handle a large number of samples $M$.
    *   The primary computational cost before SVD is often the calculation of `S_coeffs` via `lstsq(Phi_T.T, Y_T.T)`.
    *   For extremely large $M$, if an orthonormal library basis is used (e.g., after `orthogonalize_basis_gram_schmidt`), the calculation of each column of `S_coeffs` (i.e., for each sample function) becomes independent and can be parallelized.

*   **Choice of Library Basis $\mathbf{\Phi}_T$:**
    *   **Polynomials:** A common default choice (e.g., $z^0, z^1, \dots, z^d$). They are general but might require a high degree $d$ for complex functions or wide z-ranges, potentially leading to numerical instability. Scaling $z$ to a range like $[0,1]$ (via `scale_z` in `create_library_basis`) is highly recommended for polynomials.
    *   **Custom Functions:** Users can provide their own functions (e.g., exponentials, Gaussians, or functions with known physical relevance). This can lead to a more compact and meaningful library basis $\mathbf{\Phi}_T$, potentially requiring a smaller $P$.

*   **Orthogonalization of Library Basis:**
    *   The `orthogonalize_basis_gram_schmidt` function can make the library basis $\mathbf{\Phi}_T$ orthonormal (with respect to the weights `ws`, if provided).
    *   **Pros:** If `is_phi_orthogonal_weighted` is then set to `True`, `compute_coefficients_in_library` becomes a simple projection, which can be faster. An orthonormal basis can also have desirable numerical properties.
    *   **Cons:** The orthogonalization process transforms the original library functions. If the original $\phi_p(z)$ had simple analytical forms (like $z^p$), the orthogonalized functions will be more complex linear combinations, making the interpretation of `U_k_coeffs` (and thus the final `B_opt_T`) in terms of the *original* simple forms less direct. The `print_analytical_form_polynomial` function is most straightforward if the library basis used for SVD is the original (scaled) polynomial set.

*   **Choosing the Number of Optimal Basis Functions ($K$):**
    *   The parameter `num_optimal_K` is crucial. A small $K$ leads to a very compact basis but may not capture all features of the sample functions, resulting in higher reconstruction errors. A large $K$ improves accuracy but reduces compactness.
    *   The plot of singular values (`s_vals_svd` from `plot_singular_values`) is a key diagnostic. A sharp drop in singular values suggests that subsequent components contribute less significantly. $K$ is often chosen where this "elbow" occurs.
    *   Alternatively, $K$ can be chosen to capture a certain percentage of the "variance" (sum of squared singular values) or by evaluating the reconstruction error of the sample functions for different $K$.

*   **Interpretation of Optimal Basis Functions:**
    *   If the library basis $\mathbf{\Phi}_T$ consists of simple analytical functions (e.g., scaled polynomials) and is *not* orthogonalized before SVD, the `U_k_coeffs` provide a direct analytical expression for each optimal basis function $B_{\text{opt},k}(z)$ in terms of these simple functions.
    *   If $\mathbf{\Phi}_T$ is complex or has been orthogonalized, the primary interpretation of $B_{\text{opt},k}(z)$ comes from visualizing their shapes (the rows of `B_opt_T` as plotted by `plot_1d_profiles`).

*   **Numerical Stability:**
    *   As mentioned, scaling z-coordinates (e.g., to $[0,1]$ via `scale_z=True` in `library_basis_config` when using polynomials) is important for the numerical stability of the library basis generation and subsequent least-squares fitting.
    *   If the chosen library basis functions $\mathbf{\Phi}_T$ are nearly linearly dependent (ill-conditioned), the least-squares problem for `S_coeffs` can be sensitive. Orthogonalizing the basis via `orthogonalize_basis_gram_schmidt` can mitigate this, though it changes the interpretation as noted above. Regularization in `compute_coefficients_in_library` (if `is_phi_orthogonal_weighted` is False) also helps stabilize the solution.

## 6. Main Pipeline Function

The primary entry point for using this methodology is the `run_optimization_pipeline` function. It takes numerous arguments to control each step of the process, from data generation/loading to the final SVD and plotting. This allows for flexible usage, whether generating all data from scratch or providing pre-computed sample functions or library bases.
