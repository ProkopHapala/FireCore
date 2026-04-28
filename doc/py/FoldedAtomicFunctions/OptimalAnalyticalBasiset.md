## Problem: Finding an Optimal Analytical Basis Set

We are given:
1.  A set of $M$ sample functions, $\{y_j(x)\}_{j=1}^M$. Each $y_j(x)$ is typically known at $N_x$ discrete points $x_i$, so we have data $(x_i, y_{ji})$.
2.  A predefined library of $P$ simple, linearly independent analytical functions, $L = \{\phi_p(x)\}_{p=1}^P$. For example, polynomials $\phi_p(x) = x^{p-1}$ for $p=1, \dots, P$.
3.  A desired number of new, optimal basis functions, $K$, where $K < P$.

**Objective:**
Find a set of $K$ new basis functions, $B = \{b_k(x)\}_{k=1}^K$, such that:
*   Each new basis function $b_k(x)$ is an analytical function, specifically a linear combination of the functions from the library $L$:
    $$ b_k(x) = \sum_{p=1}^{P} C_{kp} \phi_p(x) \quad (*)$$
    where $C_{kp}$ are coefficients we need to determine.
*   This new basis set $B$ can represent the original sample functions $\{y_j(x)\}$ with minimum reconstruction error, in a least-squares sense.

Let $\mathbf{Y}$ be an $N_x \times M$ matrix where each column $j$ contains the values $y_j(x_i)$ of a sample function.
Let $\mathbf{\Phi}$ be an $N_x \times P$ matrix where each column $p$ contains the values $\phi_p(x_i)$ of a library function.

We want to approximate each $y_j(x)$ as:
$$ y_j(x) \approx \hat{y}_j(x) = \sum_{k=1}^{K} A_{jk} b_k(x) $$
where $A_{jk}$ are the coefficients of $y_j(x)$ in the new basis $B$.

The total error to minimize is:
$$ \mathcal{E} = \sum_{j=1}^{M} \sum_{i=1}^{N_x} (y_{ji} - \hat{y}_{ji})^2 = \sum_{j=1}^{M} \left\| \mathbf{y}_j - \sum_{k=1}^{K} A_{jk} \mathbf{b}_k \right\|^2_2 $$
where $\mathbf{y}_j$ and $\mathbf{b}_k$ are column vectors of function values at points $x_i$.
This is a non-trivial optimization problem as we need to find the matrix $C$ (defining $b_k$) and $A$ (coefficients for $y_j$) simultaneously.

## Proposed Solution: Two-Step Approach using SVD on Coefficients

This approach decouples the problem into more manageable steps:

**Step 1: Represent Sample Functions in the Full Analytical Library Basis**

First, we approximate each sample function $y_j(x)$ using the *full* library of $P$ analytical functions $\{\phi_p(x)\}$.
For each sample function $y_j(x)$, we find coefficients $S_{jp}$ such that:
$$ y_j(x) \approx \sum_{p=1}^{P} S_{jp} \phi_p(x) $$
In matrix form, for all sample functions simultaneously:
$$ \mathbf{Y} \approx \mathbf{\Phi} \mathbf{S}_{\text{coeffs}} $$
where $\mathbf{S}_{\text{coeffs}}$ is a $P \times M$ matrix. Each column $j$ of $\mathbf{S}_{\text{coeffs}}$ contains the $P$ coefficients $(S_{1j}, S_{2j}, \dots, S_{Pj})^T$ that represent the sample function $y_j(x)$ in terms of the library $\{\phi_p(x)\}$.
These coefficients can be found using a standard least-squares solution:
$$ \mathbf{S}_{\text{coeffs}} = (\mathbf{\Phi}^T \mathbf{\Phi})^{\dagger} \mathbf{\Phi}^T \mathbf{Y} $$
where $(\cdot)^{\dagger}$ denotes the Moore-Penrose pseudoinverse, which is equivalent to $(\mathbf{\Phi}^T \mathbf{\Phi})^{-1}$ if $\mathbf{\Phi}^T \mathbf{\Phi}$ is invertible.

**Step 2: Find Principal Components in the Coefficient Space**

The matrix $\mathbf{S}_{\text{coeffs}}$ ($P \times M$) now contains $M$ coefficient vectors (its columns), each of dimension $P$. We want to find the $K$ most significant "directions" or "patterns" within this set of coefficient vectors. This is precisely what Principal Component Analysis (PCA), often implemented via Singular Value Decomposition (SVD), achieves.

We perform SVD on the coefficient matrix $\mathbf{S}_{\text{coeffs}}$:
$$ \mathbf{S}_{\text{coeffs}} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T $$
where:
*   $\mathbf{U}$ is a $P \times P$ orthogonal matrix. Its columns are the left singular vectors (principal components or principal directions in the $P$-dimensional coefficient space).
*   $\mathbf{\Sigma}$ is a $P \times M$ rectangular diagonal matrix of singular values, sorted in descending order.
*   $\mathbf{V}^T$ is an $M \times M$ orthogonal matrix (its rows are right singular vectors).

The columns of $\mathbf{U}$ represent orthonormal basis vectors for the space of coefficients. The first $K$ columns of $\mathbf{U}$, corresponding to the $K$ largest singular values, capture the most variance in the coefficient data. Let this sub-matrix be $\mathbf{U}_K$ (a $P \times K$ matrix).

**Step 3: Construct the Optimal Analytical Basis Functions**

The columns of $\mathbf{U}_K$ provide the coefficients that define our new $K$ optimal analytical basis functions $b_k(x)$ in terms of the original library $\{\phi_p(x)\}$.
Specifically, the $k$-th new basis function $b_k(x)$ is defined by the $k$-th column of $\mathbf{U}_K$:
$$ b_k(x) = \sum_{p=1}^{P} (\mathbf{U}_K)_{pk} \phi_p(x) $$
Comparing this to equation $(*)$, we see that $C_{kp} = (\mathbf{U}_K)_{pk}$.
In matrix form, if $\vec{b}(x) = [b_1(x), \dots, b_K(x)]$ is a row vector of the new basis functions, and $\vec{\phi}(x) = [\phi_1(x), \dots, \phi_P(x)]$ is a row vector of the library functions, then:
$$ \vec{b}(x) = \vec{\phi}(x) \mathbf{U}_K $$
The matrix of evaluated new basis functions on the grid $x_i$ is $\mathbf{B}_{\text{new}} = \mathbf{\Phi} \mathbf{U}_K$. This $\mathbf{B}_{\text{new}}$ is an $N_x \times K$ matrix, where each column is one of the new optimal analytical basis functions.

**Properties of the New Basis:**
*   The functions $b_k(x)$ are analytical because they are linear combinations of the analytical library functions $\phi_p(x)$.
*   They are "optimal" in the sense that they capture the maximum possible variance of the original samples $y_j(x)$ *when those samples are first projected into the space spanned by the library $\{\phi_p(x)\}$*, for a given $K$.
*   If the library functions $\phi_p(x)$ are chosen to be orthonormal (e.g., Legendre polynomials over a specific interval, or orthonormalized monomials), and the discrete points $x_i$ and weights (if any) in the least-squares fit respect this orthogonality, then the resulting $b_k(x)$ will also be orthonormal. The columns of $\mathbf{U}_K$ are orthonormal, so $\int b_k(x) b_l(x) dx = \sum_{p,q} (\mathbf{U}_K)_{pk} (\mathbf{U}_K)_{ql} \int \phi_p(x) \phi_q(x) dx = \sum_p (\mathbf{U}_K)_{pk} (\mathbf{U}_K)_{pl} = \delta_{kl}$. If simple monomials $x^{p-1}$ are used (which are not orthogonal over a general interval), the $b_k(x)$ will not generally be orthogonal, but they still form the desired optimal basis in the sense described.

This method provides a constructive way to find a reduced set of analytical basis functions tailored to the specific characteristics of the sample data.
