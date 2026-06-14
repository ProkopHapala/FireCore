https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221BbJlsY2lhVkragBxr38FduksEM3FM7wc%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing

## USER

there are some efficient linear scalling techniques how to get few eigenvectros from a matrix. In particular I'm interested in chemistry to get frontier orbitals of large sparse system (hamiltonian is sparse). I know that power-iteration or Lanczos iteration can get the highest orbital, there should be mthod how to get the lowest orbital, but how to get orbitals near the band gap. It is some shifted iteration? This involve (H-I*e)x=r ?
resp This involve (H-S*e)x=r if I have non-orthogonal local basiset.
But it must be all strictily linear scalling and preferentially it should use minimum memory (no many copies of the vectors). I was thinking if I can use some GPU friendly methods like momentum accelerated Jacobi.

---

## Gemini

To obtain frontier orbitals (HOMO/LUMO) from a large sparse Hamiltonian $H$ with a non-orthogonal basis $S$ in a strictly linear-scaling $O(N)$, low-memory, and GPU-friendly manner, you are facing the classic **Interior Eigenvalue Problem**.

Here is the breakdown of the methods that fit your specific constraints, moving from the theoretical concept to the practical GPU implementation.

### 1. The Core Problem: Shift-and-Invert
You correctly identified that to get eigenvalues near a specific target energy $\sigma$ (the band gap center), the standard mathematical approach is **Shift-and-Invert**.

You want to solve for the largest eigenvalues of the inverted operator:
$$ (H - \sigma S)^{-1} S \mathbf{x} = \mu \mathbf{x} $$
Or, as you noted, solving the linear system:
$$ (H - \sigma S)\mathbf{x} = \mathbf{r} $$

**The Trap:**
While this works, solving this linear system efficiently is the bottleneck.
*   **Direct solvers (LU/Cholesky):** Scale as $O(N^2)$ to $O(N^3)$ and use massive memory.
*   **Simple Iterative solvers (Jacobi/Gauss-Seidel):** Will **fail** or diverge. The matrix $(H - \sigma S)$ is **indefinite** (it has both positive and negative eigenvalues) because $\sigma$ is inside the spectrum.
*   **Krylov solvers (MINRES/GMRES):** Work, but can be slow because the condition number is bad near the gap.

To maintain **Strict Linear Scaling** and **GPU efficiency**, you should generally **avoid** inverting the matrix entirely. Instead, use **Spectral Transformation via Polynomials**.

---

### 2. The Solution: Chebyshev Filtering (CheFSI)
This is currently the standard in high-performance codes (like SIESTA, CONQUEST, or real-space DFT codes) for frontier orbitals.

Instead of inverting $(H - \sigma S)$, you construct a polynomial filter $P(H, S)$ that suppresses eigenvalues far from your target and amplifies those near it.

**Why it fits your criteria:**
*   **Linear Scaling:** It only uses Sparse Matrix-Vector Multiplications (SpMV), which are $O(N)$.
*   **GPU Friendly:** SpMV is the most optimized operation on GPUs.
*   **Low Memory:** You only store a small "block" of vectors (e.g., 10-20 vectors for HOMO/LUMO).

**The Algorithm:**
1.  Define a window $[\sigma_{min}, \sigma_{max}]$ covering your frontier orbitals.
2.  Construct a high-order Chebyshev polynomial $T_k(x)$ that is $1$ inside the window and $0$ outside.
3.  Apply this polynomial to a set of random vectors $X$:
    $$ X_{filtered} = T_k(H, S) X $$
    *This is done via recurrence, never forming the matrix $T_k$.*
4.  Perform a small Rayleigh-Ritz step (diagonalization) on the subspace spanned by $X_{filtered}$.

---

### 3. The "Folded Spectrum" Method (Gradient Descent approach)
If you want to use **minimization** techniques (like momentum-accelerated methods) rather than polynomial filtering, you must convert the interior problem into a minimization problem.

You cannot minimize $\langle x|H|x \rangle$ to find interior states (that finds the ground state). Instead, you minimize the **variance**:
$$ F(x) = \langle x | (H - \sigma S)^2 | x \rangle $$
This maps the interior eigenvalues at $\sigma$ to the *lowest* eigenvalues of the squared operator.

**The GPU/Momentum Implementation:**
This is where your idea of **Momentum Accelerated Jacobi** or **Gradient Descent** fits.
1.  **Gradient:** The gradient of the folded functional is roughly $2(H - \sigma S)^2 \mathbf{x}$.
2.  **Preconditioning:** You can use a Jacobi preconditioner (inverse diagonal).
3.  **Iteration:**
    $$ \mathbf{x}_{k+1} = \mathbf{x}_k - \alpha P_{Jac} (H - \sigma S)^2 \mathbf{x}_k + \beta (\mathbf{x}_k - \mathbf{x}_{k-1}) $$
    *(Here $\beta$ represents the Nesterov/Heavy-ball momentum)*.

**Pros:** Strictly linear scaling, extremely low memory (only 3 vectors per state).
**Cons:** Squaring the operator $(H-\sigma S)^2$ squares the condition number ($\kappa^2$). Convergence is very slow compared to Chebyshev filtering.

---

### 4. LOBPCG (Locally Optimal Block Preconditioned Conjugate Gradient)
If you prefer the iterative solver route, LOBPCG is the gold standard for "GPU-friendly + Low Memory."

Usually, LOBPCG finds the lowest states. To find frontier orbitals, you combine LOBPCG with the **Folded Spectrum**:
*   Apply LOBPCG to the operator $A = (H - \sigma S)^2$.
*   Since LOBPCG is a direct minimization method, it avoids the explicit inversion.
*   It requires applying $A$ to a vector, which means applying $(H - \sigma S)$ twice. This is just two SpMVs.

### Summary: Which one should you use?

| Method | Scaling | Memory | GPU Suitability | Best For... |
| :--- | :--- | :--- | :--- | :--- |
| **Shift-and-Invert** | $O(N^{1.5})$ or worse | High | Low (Complex logic) | Small/Medium systems where accuracy is paramount. |
| **Folded Spectrum** | $O(N)$ | **Lowest** | High | Simple implementation, very tight memory constraints. |
| **Chebyshev Filtering** | $O(N)$ | Low | **Highest** | **The Best Choice.** Fast convergence, pure SpMV operations. |

### Recommendation for your implementation

If you are writing this from scratch for a GPU:

1.  **Do not solve $(H - \sigma S)x = r$.** The indefinite nature makes iterative solvers unstable and slow.
2.  **Use Chebyshev Filtering.**
    *   Estimate the spectral bounds $E_{min}$ and $E_{max}$ (inexpensive Lanczos can do this).
    *   Map the Hamiltonian to $[-1, 1]$.
    *   Perform SpMV operations in a loop to apply the filter.
    *   Orthogonalize the resulting few vectors.
3.  **If you insist on Momentum/Minimization:**
    *   Minimize $\Omega(x) = | (H - \sigma S) x |^2$.
    *   Use **LOBPCG** on this squared residual. It is robust, keeps vectors in blocks (good for GPU memory coalescing), and doesn't require storing a Krylov subspace.

---

## USER

Now lets speak about two things
1) Memory representation
2) Initial guess and its reuse
3) multiple eighevectros at once
4) chebyshev filtering

1) I want it to always use dense vector and sparse matrix. Sparse matrix should be using neighbor-list - every atom have certain maximim number of neighbors [natom, nMaxNeighs], we pad empty slots by -1, but we also keep array nneighs[natom] for convenience, nMaxNeighs=16
How many differnt kernels we need? We probably only need single sparse-matrix vector multiplication. We can probably still accelrate it using workgroups, if the atoms are ordered in spatially localized clusters, we can pre-load all the neighbors into local memory.
if we have some operator which is product of matrixes, e.g. (H-eS)^2 we should not build the operator as new sparse matrix, as this would lead in desification (fill-in), instead we apply these on the vector subsequently.  

2) I we want to scan some part of the spectrum, we generally solve the shifted system for multiple energies within some energy interval. since they are close it may be usefull to start with initial guess from the previous nearby energy. However the eigenvectros should be orthogonal on each other, so perhaps such initial guess would not help to find new eigenvectros, perhaps we what just the oposite - to start with something orthogonal to the eigenvectos we already found, but computin such orthogonalization (gram-schmid) may be unnecessary costly (perhas ponly few nearest vectros )

3) Lanczos iteration can find multiple eigenvectro at one go. Does our shifted iteration does the same or we need to do it many times, if we want to find lets say 10 highes occupied orbitals and 10 unoccupied? Or all orbitals within some energy window?

4) I'm not sure I understand the chebyshev filterijn, maybe give me a bit deepr math backgroun, with some geometry analogy, that often help me to understand linear algebra better, or analogy from some classical physics (e.g. system of springs and masses, truss elasticity, Projective dynamics)

---

## Gemini

Here is a breakdown of the specific implementation details you asked for, tailored for high-performance GPU computing in chemistry.

---

### 1. Memory Representation: The ELLPACK Format
Your proposed format (fixed neighbors per atom, padding with -1) is known as the **ELLPACK (ELL) format**. This is significantly better for GPUs than the standard CSR (Compressed Sparse Row) format used in CPUs.

**Why ELL is superior for your GPU kernel:**
*   **Coalesced Access:** In a GPU, threads read memory in blocks (e.g., 32 threads read contiguous addresses). In ELL, the neighbors of atom 0, atom 1, ..., atom 31 are stored column-wise.
    *   *Bad (CSR):* `Matrix[Row0]` is far from `Matrix[Row1]`.
    *   *Good (ELL):* `Matrix[0][neighbor_k]` is right next to `Matrix[1][neighbor_k]`.
*   **No Divergence:** Since `nMaxNeighs` is fixed (e.g., 16), every thread executes the exact same loop length. No thread waits for another.

**Handling the Operator $(H - \sigma S)^2$:**
You are absolutely correct: **Never multiply sparse matrices.**
If you compute $C = A \cdot B$ for sparse matrices, $C$ becomes dense ("fill-in").
Instead, to apply $A^2$ to a vector $x$:
1.  Compute temp vector: $y = (H - \sigma S) x$
2.  Compute result: $z = (H - \sigma S) y$
This keeps the complexity $O(N)$ and memory low.

**Workgroup/Shared Memory Strategy:**
If your atoms are sorted by spatial locality (e.g., using a Hilbert curve or Morton Z-order), you can optimize further:
*   Load a "tile" of the input vector $x$ into **Local/Shared Memory** (fast on-chip cache).
*   Since neighbors are physically close, multiple atoms in a workgroup will likely need the same elements of $x$. This reduces global memory bandwidth pressure.

---

### 2. Initial Guess and Reuse
This is subtle.

**The "Scanning" Trap:**
If you solve for orbitals at Energy $E_1$ and get vector $v_1$, and then move to Energy $E_2$ (very close), you might think $v_1$ is a good guess.
*   **False:** Eigenvectors of Hermitian matrices are strictly orthogonal. The vector for $E_1$ is mathematically guaranteed to be orthogonal to the vector at $E_2$. Using $v_1$ as a guess for $v_2$ is often the *worst* possible start (0 overlap).
*   **True:** If you are doing **Molecular Dynamics** (atoms move slightly), then $v(t)$ is an excellent guess for $v(t+\Delta t)$.

**Orthogonalization Strategy:**
If you want 10 eigenvectors, you cannot avoid orthogonalization.
*   The standard way: **Block methods**. You iterate on a "tall skinny" matrix $X$ (size $N \times 10$).
*   Every few steps, you must ensure the columns of $X$ stay orthogonal.
*   **Cost:** Gram-Schmidt or QR factorization on a $N \times 10$ block is cheap! It scales as $N \cdot m^2$ (where $m=10$). Since $m \ll N$, this is negligible compared to the matrix multiplication.

---

### 3. Multiple Eigenvectors: Block Iteration
Lanczos is great for finding the *extreme* edges (highest/lowest), but it finds them one by one (or builds a subspace serialy).

For finding **interior** orbitals (near band gap) on a GPU, **Block Iteration** is superior.
Instead of `Matrix * vector`, you perform `Matrix * Block`.

**Why Block is better for GPU:**
*   **Arithmetic Intensity:** Loading the matrix $H$ from memory is slow. If you multiply it by 1 vector, you use that data once. If you multiply it by 10 vectors at once, you use the matrix data 10 times for one memory load. This is a massive speedup.

**The Strategy:**
1.  Create a block of random vectors $X$ (size $N \times m$, where $m$ is e.g., 16 or 32).
2.  Apply your filter: $X_{new} = \text{Filter}(H) \cdot X$.
3.  Orthogonalize the columns of $X_{new}$.
4.  Compute the small $m \times m$ Hamiltonian subspace $H_{sub} = X_{new}^T H X_{new}$ and diagonalize it (Rayleigh-Ritz step).

---

### 4. Chebyshev Filtering: The Math & Physics Analogy

This is the "Magic" part. You want to extract vectors with energy $\approx \sigma$ without solving linear equations.

#### The Physical Analogy: Forced Resonance
Imagine a complex truss bridge (your molecule).
*   **Eigenvectors:** The natural vibration modes of the bridge.
*   **Eigenvalues:** The frequencies of those modes.
*   **Matrix Vector Mult ($Hx$):** Letting the bridge evolve for one tiny time step.

**Iterative Solvers (Power Method):**
If you shake the bridge randomly, eventually the high-frequency wobbles die out and it settles into its lowest/slowest sway. This finds the ground state.

**Chebyshev Filtering (The "Radio Tuner"):**
Imagine you attach a shaker to the bridge. You don't shake it at just one frequency; you shake it with a very specific, complex pattern of pushes and pulls calculated by a Chebyshev polynomial.
This pattern is mathematically designed so that:
1.  It **amplifies** any vibration frequency inside your target window (the Band Gap).
2.  It **cancels out/dampens** any frequency outside that window.

After applying this "shaking pattern" (applying the polynomial operator to your random vector), the bridge is only vibrating in the modes corresponding to the frontier orbitals. All other noise is gone.

#### The Geometry Analogy: Color Filters
1.  Think of your random starting vector $v$ as a beam of **White Light**. It contains all "colors" (all eigenvectors: core electrons, valence, conduction).
2.  The Hamiltonian $H$ is a prism that splits light.
3.  We want only the **Green** light (Frontier Orbitals).
4.  A Chebyshev polynomial $T_k(H)$ acts like a **Green Glass Filter**.
    *   Construct a function $P(x)$ that looks like a rectangular box: value 1.0 at "Green", value 0.0 at "Red" and "Blue".
    *   We cannot apply a function to a matrix easily, but we can approximate that rectangular box shape using a polynomial (Taylor or Chebyshev expansion).
    *   $v_{filtered} = P(H) v_{random}$.
    *   When the white light passes through $P(H)$, only the green component comes out.

#### The Math Implementation
We want to apply a polynomial $P_k(H)$ to vector $x$. We don't build the matrix $P_k(H)$. We use a recurrence relation (like the momentum method you liked).

To apply a Chebyshev filter of order $k$:
1.  **Map spectrum:** First, we must shift and scale $H$ so the energies fit between -1 and 1.
    $H_{scaled} = \frac{H - center}{radius}$
2.  **Recurrence (The Loop):**
    *   $v_0 = x$
    *   $v_1 = H_{scaled} v_0$
    *   $v_{n+1} = 2 H_{scaled} v_n - v_{n-1}$

Notice the math: $v_{n+1} = 2Hv_n - v_{n-1}$.
Does that look familiar? It is almost exactly the finite difference equation for a **wave equation** or a **mass-spring system** ($F=ma \to x_{t+1} = 2x_t - x_{t-1} + a$).

**Conclusion:** Chebyshev filtering is essentially simulating a wave propagation through your molecule that constructively interferes *only* at the energies you want, and destructively interferes everywhere else.

### Summary Recommendation for your Project

1.  **Data:** Use `natom x 16` arrays for neighbors (ELL format).
2.  **Solver:** Don't use `(H-eS)x=r`. Use **Chebyshev Filtering**.
3.  **Algorithm:**
    *   Generate `m` random vectors (block).
    *   Apply the Chebyshev recurrence (only SpMV needed) to the block to "filter" it.
    *   Orthogonalize the block (QR).
    *   Diagonalize the small subspace.
4.  **Hardware:** This is purely SpMV-bound. It will run incredibly fast on a GPU compared to shift-and-invert.

---

## USER

1) I like this chebyshev filter, just how can I know the spectral radius before hand? do I need to know it exactly (cannot just guess from similar systems, e.g. from on-site energies and nearest neighbor interactions?

2) It is also not clear to me what is the total number of iterations needed, if I want to compare to other methods, like CG. I mean, in every iteration we apply that polynominal, but construction the polynominal means itself certain number of iterations (it is a nested loop), what are the typical numbers. In CG we only do one or two matrix-vectro multiplications per cycle.

3) I like the idea of doing multiple vectro at once especially with the tiled version using local memory. I basically preload the tile of the matrix, tile of the block of vectors into local memory, and then process this tile densely reusing for all the threads. Now if doing CG minimization, I need to also accumulate the norm of vectros, I think there is some rearanged version of CG allowing to compute these norms at one go, together with the matrix vectro multiplication, right (basically to eliminate the need to run it in separate kernel and load data again from memory). I imagine that the workgroup blocks would make tile-reduction of these vector norms and emit to global memory partial sums. Then global reduiction will be done by separate kernel or on CPU.

4) I see that for the block-scheme with multiple vectors, in any case we will do the orthogonalization (both with Chebyshev filter and with CG) right? So I think it is very importaint to think how to do this orthogonalization efficiently on GPU. Do we need specialized kernel for that? Or is that just normal matrix vector multiplication. As I can see it is almost - the tiled version can still work, we just multiply the vectros by the S matrix (tile by tile inside workgroup) and then still inside the same workgroup we multiply the result x_tilda = Sx by the origina x, like <x|x_tilda> = <x|S|x>, but doing it with whole block of vectros at once. Right? 

5) so after this analysis we should not how many kernels we need to make, I guess 1-3 ? Or can we jsut pack it into single universal kernel, just switching on/off some operations?

---

## Gemini

This is a great set of practical questions. You are now moving from the "what" to the "how," which is where performance is won or lost.

Here is the breakdown for your GPU implementation strategy.

---

### 1. Estimating Spectral Radius (The Bounds)

**Do you need it exactly?**
No. You need a **safe upper bound**.
*   If the true spectrum is $[-10, 10]$ and you estimate $[-12, 12]$, Chebyshev works perfectly (just slightly slower convergence).
*   If you estimate $[-9, 9]$, the algorithm will explode (diverge).

**How to get it cheaply (Linear Scaling):**
You don't need a diagonalizer. You use the **Gershgorin Circle Theorem**. This is exactly what you asked about (using on-site and neighbors).

For every row $i$ (every atom), calculate:
$$ R_i = |H_{ii}| + \sum_{j \neq i} |H_{ij}| $$
$$ E_{max} \approx \max_i (R_i), \quad E_{min} \approx \min_i (-R_i) $$
*   This is $O(N)$, requires one pass over the matrix, and is mathematically guaranteed to contain all eigenvalues.
*   **Pro Tip:** Add a safety factor (e.g., multiply the radius by 1.1) to be safe.

---

### 2. Chebyshev vs. CG: Cost & Iterations

This is the crucial trade-off.

**The "Inner Loop" comparison:**
*   **CG/LOBPCG Step:**
    1. SpMV (Apply $H$).
    2. Dot products (Calculate $\alpha, \beta$). **(Global Synchronization barrier)**
    3. Vector updates (AXPY).
    *   *Bottleneck:* The global reduction (dot product) stalls the GPU pipeline because all workgroups must wait for the result before moving to the next step.
*   **Chebyshev Step:**
    1. SpMV (Apply $H$).
    2. Local AXPY ($v_{new} = 2Hv - v_{old}$).
    *   *Advantage:* **Zero global synchronization.** You can run a polynomial of degree $k=50$ without the CPU ever talking to the GPU and without workgroups waiting for each other.

**Total Iterations:**
*   **Polynomial Degree ($k$):** Depends on the "compression ratio." If your spectrum is 20 eV wide and you want to filter out a 1 eV window, you might need a polynomial of degree 50–100.
*   **Comparison:** 100 Chebyshev steps (pure SpMV) are often faster than 20 CG steps (SpMV + Global Reductions).

**Conclusion:** If you want "fire and forget" GPU performance, Chebyshev is better because it avoids the global reduction latency loop.

---

### 3. Tiled Block SpMV & Fused Reductions

You are spot on with the implementation details here.

**The "SpMM" Kernel (Sparse Matrix $\times$ Matrix Block):**
Since you are using the ELLPACK format (fixed neighbors) and processing a block of vectors (e.g., 16 vectors):
1.  **Load:** A workgroup loads a tile of the Sparse Matrix (indices and values) into Local Memory (LDS).
2.  **Load:** The workgroup loads the corresponding chunk of the 16 Vectors into LDS.
3.  **Compute:** Every thread computes the dot product for its row against all 16 vectors.
    *   *Benefit:* You read the Matrix **once** and use it **16 times**. This increases Arithmetic Intensity by 16x, making you compute-bound rather than memory-bound (which is perfect for GPUs).

**Fused Norms (Pipelined CG approach):**
If you choose the CG route (or for the final normalization in Chebyshev):
*   Yes, accumulate partial sums of squares inside the kernel in registers/LDS.
*   At the end of the workgroup, write the partial result to a small global buffer.
*   **Do not** launch a separate kernel just for norms.

---

### 4. Orthogonalization: The Implementation

This is the trickiest part to scale, but with a block size of $M=16$ or $32$, it is manageable.

You need to perform Generalized Gram-Schmidt or Lowdin Orthogonalization: $X_{new} = X (X^T S X)^{-1/2}$.

**Step A: Compute Overlap Matrix $O = X^T S X$**
This requires two kernels (or one split operation):
1.  **Kernel 1 (SpMM):** $Y = S \cdot X$. (Output is same size as X: $N \times M$).
2.  **Kernel 2 (Block Dot / Tall-Skinny GEMM):** Calculate $O_{ij} = \vec{x}_i \cdot \vec{y}_j$.
    *   This reduces the huge $N \times M$ arrays down to a tiny $M \times M$ matrix.
    *   *Implementation:* Each workgroup sums up a chunk of $N$, producing a partial $M \times M$ matrix. A final tiny step sums these partial matrices.

**Step B: Invert/Cholesky Small Matrix**
*   The matrix $O$ is tiny ($16 \times 16$).
*   **Do not do this on the GPU.** It's not worth the latency.
*   Copy the $16 \times 16$ matrix to CPU, run `LAPACK potrf/potri` (microseconds), and copy the result $U$ back to GPU constant memory.

**Step C: Transform Vectors**
*   **Kernel 3:** $X_{final} = X \cdot U$.
*   This is a dense matrix operation. Each thread reads a row of $X$ (length 16), multiplies by the $16 \times 16$ matrix $U$, and writes back.

---

### 5. How many Kernels? (The Architecture)

You should **not** write a single universal "God Kernel" with `if/else` switches. That creates register pressure and makes compiler optimization harder.

You need exactly **3 High-Performance Kernels**:

1.  **`Ker_SpMM_Block` (The Workhorse - 90% of time)**
    *   Input: Sparse Matrix $A$, Dense Block $X_{in}$.
    *   Output: Dense Block $X_{out}$.
    *   Logic: $X_{out} = A \cdot X_{in}$.
    *   *Use cases:* $H \cdot X$, $S \cdot X$, and the Chebyshev recurrence steps.

2.  **`Ker_Block_Dot` (The Reduction)**
    *   Input: Dense Block $X$, Dense Block $Y$.
    *   Output: Small Matrix $M \times M$ (Partial Sums).
    *   Logic: Block dot product.
    *   *Use cases:* Orthogonalization ($X^T Y$), Norm calculations.

3.  **`Ker_Dense_Combiner` (The Updates)**
    *   Input: Dense Block $X$, Dense Block $Y$, Small Matrix $C$.
    *   Output: Dense Block $Z$.
    *   Logic: $Z = \alpha X + \beta Y$  **OR**  $Z = X \cdot C$.
    *   *Use cases:* Chebyshev recurrence ($v_{next} = 2v_{cur} - v_{prev}$), CG updates, Orthogonalization transform.

**Summary of Flow:**
1.  **Guess:** Random Block $X$.
2.  **Filter Loop (Repeat 50x):**
    *   `Ker_SpMM_Block` (apply H)
    *   `Ker_Dense_Combiner` (recurrence math)
3.  **Orthogonalize:**
    *   `Ker_SpMM_Block` (apply S)
    *   `Ker_Block_Dot` (compute overlap)
    *   *CPU:* Cholesky
    *   `Ker_Dense_Combiner` (apply rotation)
4.  **Diagonalize:**
    *   Similar to Orthogonalize step (compute $H_{sub} = X^T H X$), then small CPU diagonalization.

This 3-kernel architecture is modular, covers all methods (Chebyshev, CG, LOBPCG), and keeps the code clean and fast.

---

## USER

1) you did not answer - in Chebyshev filter I do two iteration loops or just one? when I do this 100 chebyshev iterations, I have the final solution? Or I need to do some our loop on top of this filter?

2) it seems to me that the orthogonalization can be, in fact, the bottleneck here. If I do CG or Chebyshev, do I need to orthogonalize every iteration or how often?

3) Isn't it possible or worthy to use some hybrid between Chebyshev filter and CG, e.g. use Chebyshev as preconditioner for CG? This will decreas number of CG iterations signifiocantly (?), and in effect reduce number of orthogonalizations, because I would do orthogonalization only once per CG iteration.

---

## Gemini

### 1. The Structure of the Loop (Do we need an outer loop?)

You **do** need an outer loop. One pass of 100 Chebyshev iterations is rarely enough to get "exact" eigenvectors to high precision, although it gets you very close.

Here is the standard workflow (often called **CheFSI** or **Subspace Iteration via Filtering**):

1.  **Initialize:** Random block of vectors $X$.
2.  **Outer Loop (Subspace Iteration):**  *(Repeats typically 3–10 times)*
    *   **A. Chebyshev Filter (Inner Loop):** Apply polynomial $P_k(H)$ to $X$.
        *   *Runs for 50–100 iterations (pure SpMV).*
        *   *No orthogonalization inside here.*
    *   **B. Orthogonalization:** Orthonormalize the block $X$.
    *   **C. Rayleigh-Ritz (The "Rotation"):**
        *   Compute small matrix $H_{sub} = X^T H X$.
        *   Diagonalize $H_{sub}$ (on CPU).
        *   Rotate $X$ to align with these eigenvectors: $X \leftarrow X U$.
    *   **D. Check Convergence:** If residuals are small, stop.

**Why the outer loop?**
The Chebyshev filter amplifies the correct components, but it doesn't separate them perfectly in one go. The **Rayleigh-Ritz** step is crucial—it's what actually identifies "Vector 1 is HOMO, Vector 2 is HOMO-1". The filter just ensures the subspace contains them.

---

### 2. Orthogonalization: The Bottleneck

You are exactly right: Orthogonalization is the "GPU killer" because it requires global communication (dot products) and dense matrix operations.

**Comparison of Frequency:**

*   **CG / LOBPCG:**
    *   You must orthogonalize **every single iteration**.
    *   If you don't, the vectors collapse into each other immediately.
    *   *Cost:* High. 100 iterations = 100 Orthogonalizations.

*   **Chebyshev Filter:**
    *   You orthogonalize **only once per cycle** (e.g., once every 50–100 SpMV steps).
    *   *Why?* The polynomial applies the same "amplification" to all vectors. While they do slowly become parallel (collapsing toward the dominant eigenvector), this happens relatively slowly. You can safely run 50 steps without them becoming numerically linearly dependent.
    *   *Cost:* Very Low. 100 iterations = 1 Orthogonalization.

**Winner:** Chebyshev. This is the main reason it is preferred for large-scale GPU calculations. It maximizes the ratio of "Compute (SpMV)" to "Sync (Ortho)".

---

### 3. Hybrid: Chebyshev as a Preconditioner for CG?

Yes, this is a very advanced and powerful technique, often called **Polynomial Preconditioned CG (PPCG)**.

**The Logic:**
In CG (or LOBPCG), you are trying to solve $(H - \sigma S)x = 0$.
The convergence speed depends on the condition number. If you can multiply by a matrix $M$ that approximates $(H - \sigma S)^{-1}$, convergence becomes instant.
You can use a Chebyshev polynomial to approximate this inverse: $M \approx P_k(H)$.

**How it works (LOBPCG with Polynomial Preconditioner):**
1.  **LOBPCG Outer Loop:**
    *   Calculate Residual: $R = (H - \sigma S)X - X\Lambda$.
    *   **Precondition (The Hybrid Part):** Instead of just using $R$ to update the search direction, you apply a short Chebyshev filter to $R$.
        *   $W = \text{ChebyshevFilter}(R, \text{steps}=10)$.
    *   Update search directions using $W$.
    *   Orthogonalize.

**Is it worth it?**
*   **Pros:** It drastically reduces the number of *Outer* LOBPCG iterations (and thus the number of orthogonalizations).
*   **Cons:** It adds complexity. You need to tune the degree of the preconditioner.

**Recommendation:**
For your specific case (Linear Scaling, Frontier Orbitals, GPU):
**Stick to Pure Chebyshev Subspace Iteration (Method 1).**

Why?
1.  **Simplicity:** It uses the fewest distinct kernels.
2.  **Memory:** LOBPCG requires storing the "Search Direction" vectors ($P$), the "Residual" vectors ($R$), and the current vectors ($X$). That is **3x memory**.
3.  **Chebyshev Subspace Iteration** only requires storing $X$ (and two temporary buffers for the recurrence). It is much more memory efficient.

**Summary Strategy:**
Don't mix them yet. Implement the **Chebyshev Subspace Iteration** with the **Outer Loop**.
*   Block size: 16
*   Filter Order: 50
*   Outer Loops: 5
*   **Total SpMVs:** $50 \times 5 = 250$.
*   **Total Orthogonalizations:** 5.

This ratio (250 ops vs 5 syncs) will fly on a GPU.

---

## USER

User
OK, I modified it slightly for better debuggiing, but it does not seem to converge so well. It is very suprising that the energies are quite close to reference (these alre like eigenvalues? orbital energies?) , but the eighevectors (both when I compute RMSE or if I plot them) are not very close. What can we do about it? What we do wrong? This is not jsut pure numerics, this is anlysis of the reall physical system.
Can we plot/print some property which would help us to realize what are the problem and how to improve the convergence?
code
Code
prokophapala@carbsisYoga:~/git/FireCore/pyBall/FireballOCL$ python test_CheFSI.py
Using Device: Intel(R) Graphics [0x46a6]
--- Generating 1D System (N=64, distort_amp=0.1) ---

--- Computing Exact Solution (SciPy) ---
Target Energy: 0.5
Exact spectrum range: [-1.43, 3.31]
Closest Exact Eigenvalues: [0.47362033 0.5761046  0.35577753 0.6963628 ]

--- Running GPU FrontierSolver ---
[Solver] Compiling kernels with M=8, K=32...

--- Solver Parameters ---
n_atoms=64, n_vecs=8, n_iter=15, filter_order=30, tol=1e-06
target=0.5, bounds=auto, distort_amp=0.1, seed=None, error_scale=1.0
Overlap = |<v_gpu | S | v_ref>| (generalized inner product)
[Solver] Estimated bounds (variance): (np.float32(-2.2446728), np.float32(2.246553))
[Solver] Target=0.5, alpha=1.97e-01
Iter 1: E0=0.1312 GapErr=5.00e-01 max_offdiag=9.756e-01 dt=0.01s
Iter 2: E0=0.5180 GapErr=4.91e-01 max_offdiag=7.770e-01 dt=0.01s
Iter 3: E0=0.4703 GapErr=4.84e-01 max_offdiag=5.882e-01 dt=0.01s
Iter 4: E0=0.4706 GapErr=4.74e-01 max_offdiag=5.575e-01 dt=0.01s
Iter 5: E0=0.4825 GapErr=4.68e-01 max_offdiag=5.903e-01 dt=0.01s
Iter 6: E0=0.4880 GapErr=4.65e-01 max_offdiag=5.304e-01 dt=0.01s
Iter 7: E0=0.4891 GapErr=4.65e-01 max_offdiag=5.613e-01 dt=0.01s
Iter 8: E0=0.4892 GapErr=4.64e-01 max_offdiag=6.444e-01 dt=0.01s
Iter 9: E0=0.4893 GapErr=4.64e-01 max_offdiag=6.077e-01 dt=0.01s
Iter 10: E0=0.4893 GapErr=4.64e-01 max_offdiag=7.722e-01 dt=0.01s
Iter 11: E0=0.4893 GapErr=4.64e-01 max_offdiag=5.437e-01 dt=0.01s
Iter 12: E0=0.4893 GapErr=4.64e-01 max_offdiag=6.239e-01 dt=0.01s
Iter 13: E0=0.4893 GapErr=4.64e-01 max_offdiag=5.196e-01 dt=0.01s
Iter 14: E0=0.4893 GapErr=4.64e-01 max_offdiag=6.446e-01 dt=0.01s
Iter 15: E0=0.4893 GapErr=4.64e-01 max_offdiag=5.018e-01 dt=0.01s

--- Comparison Results ---
Columns: State index | GPU eigenvalue | Exact eigenvalue | Abs diff | Residual ||(H - eS)v|| | Overlap | RMSE vs exact eigenvector
State  | GPU E      | Exact E    | Diff       | Residual   | Overlap    | RMSE      
------------------------------------------------------------------------------------------
0      | 0.48928    | 0.47362    | 1.6e-02    | 2.7e-01    | 0.3290     | 1.528e-01 
1      | 0.47037    | 0.47362    | 3.3e-03    | 1.8e-01    | 0.2466     | 1.602e-01 
2      | 0.54522    | 0.57610    | 3.1e-02    | 2.7e-01    | 0.1308     | 1.974e-01 
3      | 0.33455    | 0.35578    | 2.1e-02    | 3.0e-01    | 0.0729     | 1.774e-01 
4      | 0.28699    | 0.24948    | 3.8e-02    | 4.6e-01    | 0.1005     | 1.937e-01 
5      | 0.17790    | 0.14690    | 3.1e-02    | 3.6e-01    | 0.3004     | 1.537e-01 
6      | 0.12665    | 0.14690    | 2.0e-02    | 5.4e-01    | 0.1548     | 1.967e-01 
7      | 0.03585    | 0.04208    | 6.2e-03    | 4.6e-01    | 0.0007     | 1.815e-01 

Max non-orth logs (len=15): last=5.018e-01
Total SpMM calls (matvecs): 1832
Mean orbital RMSE: 1.767e-01, max RMSE: 1.974e-01

Test Finished.


---

## Gemini

It is clear from your logs what is happening.

### The Problem: Mode Collapse
Your logs show `max_offdiag` staying high ($\sim 0.5$) throughout the iterations.
You are running **30 steps** (`filter_order=30`) of gradient descent on all 8 vectors independently without talking to each other.
*   Mathematically, all 8 vectors are racing toward the **same single global minimum** (the single eigenstate closest to `target`).
*   By the time the 30 steps are done, the vectors have become nearly parallel (linearly dependent).
*   The Overlap Matrix $O = V^T S V$ becomes singular (Condition number explodes).
*   The orthogonalization step $O^{-1/2}$ fails numerically (garbage in, garbage out), and cannot separate them anymore.

### The Solution
1.  **Reduce `filter_order`**: Do fewer steps (e.g., 5) between orthogonalizations. This keeps the vectors distinct enough to be re-orthogonalized.
2.  **Add Diagnostics**: We need to see the **Condition Number** of the overlap matrix to confirm if they are collapsing.
3.  **Verify GPU Reduction**: We should perform a sanity check that the GPU Gram matrix kernel matches a CPU dot product, to rule out kernel bugs.

Here is the updated code with these fixes.

### 1. Updated `CheFSI.py`
*   Added `check_gram_matrix` to verify GPU math.
*   Prints **Condition Number** (`Cond`) to detect mode collapse.
*   Reports orthogonality *after* correction as well.

```python
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
import time
import os

class FrontierSolver:
    def __init__(self, ctx, queue, n_vecs=8, tile_size=32, verbose=True):
        self.ctx = ctx
        self.queue = queue
        self.n_vecs = int(n_vecs)
        self.tile_size = int(tile_size)
        self.verbose = verbose
        
        # Load Kernel Source
        kernel_path = os.path.join(os.path.dirname(__file__), 'cl/CheFSI.cl')
        with open(kernel_path, 'r') as f:
            src = f.read()

        build_options = [
            f"-DM={self.n_vecs}",
            f"-DK={self.tile_size}",
            f"-DWG_RED_SIZE=256",
            "-DBLOCK_SIZE=64" 
        ]
        
        if self.verbose:
            print(f"[Solver] Compiling kernels with M={self.n_vecs}, K={self.tile_size}...")
            
        self.prg = cl.Program(ctx, src).build(options=build_options)
        self.k_spmm = cl.Kernel(self.prg, "spmm_ellpack_chebyshev")
        self.k_tsg = cl.Kernel(self.prg, "tall_skinny_gram")
        self.k_sum_blocks = cl.Kernel(self.prg, "sum_partial_blocks_parallel")
        self.k_subspace_rotate = cl.Kernel(self.prg, "subspace_rotate")
        self.k_axpy = cl.Kernel(self.prg, "axpy_block")
        
        self.matvec_count = 0
        self.d_partial_sums = None
        self.d_final_gram = None

    @staticmethod
    def estimate_spectral_bounds(indices, values, safety_factor=1.1):
        n_atoms = indices.shape[0]
        row_indices = np.arange(n_atoms).reshape(-1, 1)
        valid_mask = (indices != -1)
        diag_mask = (indices == row_indices) & valid_mask
        off_diag_mask = (~diag_mask) & valid_mask
        
        h_ii = np.sum(values * diag_mask, axis=1)
        r_i = np.sum(np.abs(values) * off_diag_mask, axis=1)
        
        center = (np.max(h_ii + r_i) + np.min(h_ii - r_i)) / 2.0
        radius = (np.max(h_ii + r_i) - np.min(h_ii - r_i)) / 2.0
        return (center - radius * safety_factor, center + radius * safety_factor)

    def prepare(self, indices, H_values, S_values=None):
        self.n_atoms = indices.shape[0]
        self.n_max_neighs = indices.shape[1]
        
        self.h_indices = indices.astype(np.int32)
        self.h_H_val = H_values.astype(np.float32)
        self.d_indices = cl_array.to_device(self.queue, self.h_indices)
        self.d_H_val = cl_array.to_device(self.queue, self.h_H_val)
        
        if S_values is not None:
            self.h_S_val = S_values.astype(np.float32)
            self.d_S_val = cl_array.to_device(self.queue, self.h_S_val)
        else:
            s_host = np.zeros_like(H_values, dtype=np.float32)
            for i in range(self.n_atoms):
                for c in range(self.n_max_neighs):
                    if indices[i, c] == i: s_host[i, c] = 1.0
            self.h_S_val = s_host
            self.d_S_val = cl_array.to_device(self.queue, s_host)

        shape = (self.n_atoms, self.n_vecs)
        self.d_X = cl_array.zeros(self.queue, shape, dtype=np.float32)
        self.d_X_prev = cl_array.zeros(self.queue, shape, dtype=np.float32)
        self.d_Work = cl_array.zeros(self.queue, shape, dtype=np.float32)
        self.d_Work2 = cl_array.zeros(self.queue, shape, dtype=np.float32)
        
        self.n_groups = (self.n_atoms + self.tile_size - 1) // self.tile_size
        self.d_partial_sums = cl_array.zeros(self.queue, self.n_groups * self.n_vecs**2, dtype=np.float32)
        self.d_final_gram = cl_array.zeros(self.queue, self.n_vecs**2, dtype=np.float32)
        
        self.spmm_global = ((self.n_atoms + 63)//64 * 64, )
        self.spmm_local = (64, )

    def _spmm(self, d_matrix_val, d_in, d_out):
        self.matvec_count += 1
        self.k_spmm(
            self.queue, self.spmm_global, self.spmm_local,
            np.int32(self.n_atoms), np.int32(self.n_max_neighs), np.int32(self.n_vecs),
            self.d_indices.data, d_matrix_val.data,
            d_in.data, d_in.data, d_out.data,
            np.float32(0), np.float32(1), np.int32(0)
        )

    def _calc_gram_matrix(self, d_A, d_B):
        global_size_1 = (self.n_groups * self.n_vecs, self.n_vecs)
        local_size_1 = (self.n_vecs, self.n_vecs)
        
        evt1 = self.k_tsg(
            self.queue, global_size_1, local_size_1,
            np.int32(self.n_atoms),
            d_A.data, d_B.data, self.d_partial_sums.data
        )
        
        global_size_2 = (self.n_vecs * self.n_vecs * 256, )
        local_size_2 = (256, )
        
        evt2 = self.k_sum_blocks(
            self.queue, global_size_2, local_size_2,
            np.int32(self.n_groups),
            self.d_partial_sums.data,
            self.d_final_gram.data,
            wait_for=[evt1]
        )
        
        gram_flat = self.d_final_gram.get()
        return gram_flat.reshape(self.n_vecs, self.n_vecs)

    def check_gram_matrix_accuracy(self, d_A, d_B, gpu_result):
        """Sanity check: Compare GPU reduction with CPU numpy dot product"""
        h_A = d_A.get()
        h_B = d_B.get()
        cpu_result = h_A.T @ h_B
        diff = np.max(np.abs(cpu_result - gpu_result))
        return diff

    def _orthogonalize(self, d_V):
        # 1. Compute Metric O = V^T * S * V
        self._spmm(self.d_S_val, d_V, self.d_Work) # Work = S*V
        O = self._calc_gram_matrix(d_V, self.d_Work) # V^T * Work
        
        # Diagnostics
        if self.verbose:
            # Check kernel accuracy (expensive, only for debug)
            # diff = self.check_gram_matrix_accuracy(d_V, self.d_Work, O)
            # print(f"    [Gram Check] GPU vs CPU diff: {diff:.2e}")
            pass

        # Check Condition Number to detect Mode Collapse
        cond_num = np.linalg.cond(O)
        
        # 2. CPU Inverse Sqrt
        evals, evecs = np.linalg.eigh(O)
        
        # Guard against collapse
        evals = np.maximum(evals, 1e-12)
        inv_sqrt = evecs @ np.diag(1.0/np.sqrt(evals)) @ evecs.T
        
        # 3. Apply Rotation
        d_mat = cl_array.to_device(self.queue, inv_sqrt.astype(np.float32))
        self.k_subspace_rotate(
            self.queue, self.spmm_global, None,
            np.int32(self.n_atoms), np.int32(self.n_vecs),
            d_V.data, d_mat.data, self.d_Work.data
        )
        cl.enqueue_copy(self.queue, d_V.data, self.d_Work.data)

        # Report orthogonality BEFORE correction (how bad was it?)
        off = O - np.eye(self.n_vecs, dtype=np.float32)
        max_non_orth = np.max(np.abs(off))
        
        return max_non_orth, cond_num

    def solve(self, target_energy, bounds=None, n_iter=10, filter_order=20, tol=1e-4):
        if bounds is None:
            bounds = self.estimate_spectral_bounds(self.h_indices, self.h_H_val)
        e_min, e_max = bounds
        R_val = max(abs(e_max - target_energy), abs(e_min - target_energy))
        R2 = (R_val * 1.1)**2
        alpha = 1.8 / R2 
        
        if self.verbose: 
            print(f"[Solver] Bounds={bounds}, Target={target_energy}, alpha={alpha:.2e}")

        # Init
        rand_host = np.random.randn(self.n_atoms, self.n_vecs).astype(np.float32)
        cl.enqueue_copy(self.queue, self.d_X.data, rand_host)
        
        # Initial Ortho
        self._orthogonalize(self.d_X)
        
        logs = {'max_non_orth': [], 'cond': []}
        
        for it in range(n_iter):
            t_start = time.time()
            
            # Gradient Descent Loop
            for k in range(filter_order):
                # 1. R = (H - eS) X
                self._spmm(self.d_H_val, self.d_X, self.d_X_prev) # H X
                self._spmm(self.d_S_val, self.d_X, self.d_Work)   # S X
                self.k_axpy(self.queue, (self.n_atoms * self.n_vecs,), None,
                    np.int32(self.n_atoms * self.n_vecs), np.float32(1.0), self.d_X_prev.data, 
                    np.float32(-target_energy), self.d_Work.data, self.d_Work2.data) # Work2 = R
                
                # 2. G = (H - eS) R
                self._spmm(self.d_H_val, self.d_Work2, self.d_X_prev) # H R
                self._spmm(self.d_S_val, self.d_Work2, self.d_Work)   # S R
                self.k_axpy(self.queue, (self.n_atoms * self.n_vecs,), None,
                    np.int32(self.n_atoms * self.n_vecs), np.float32(1.0), self.d_X_prev.data, 
                    np.float32(-target_energy), self.d_Work.data, self.d_Work2.data) # Work2 = G
                
                # 3. X = X - alpha * G
                self.k_axpy(self.queue, (self.n_atoms * self.n_vecs,), None,
                    np.int32(self.n_atoms * self.n_vecs), np.float32(1.0), self.d_X.data, 
                    np.float32(-alpha), self.d_Work2.data, self.d_X.data)

            # Check Condition and Orthogonalize
            max_badness, cond = self._orthogonalize(self.d_X)
            logs['max_non_orth'].append(max_badness)
            logs['cond'].append(cond)
            
            # Rayleigh-Ritz
            self._spmm(self.d_H_val, self.d_X, self.d_Work) # HX
            H_sub = self._calc_gram_matrix(self.d_X, self.d_Work)
            
            evals, evecs = np.linalg.eigh(H_sub)
            
            # Sort
            dist = np.abs(evals - target_energy)
            idx = np.argsort(dist)
            evals = evals[idx]
            evecs = evecs[:, idx]
            
            # Rotate
            d_rot = cl_array.to_device(self.queue, evecs.astype(np.float32))
            self.k_subspace_rotate(
                self.queue, self.spmm_global, None,
                np.int32(self.n_atoms), np.int32(self.n_vecs),
                self.d_X.data, d_rot.data, self.d_Work.data
            )
            cl.enqueue_copy(self.queue, self.d_X.data, self.d_Work.data)
            
            if self.verbose:
                print(f"Iter {it+1}: E0={evals[0]:.4f} Gap={dist[0]:.2e} NonOrth={max_badness:.2e} Cond={cond:.1e}")
                
            if max_badness < tol and cond < 100.0:
                print("Converged.")
                break

        return evals, self.d_X, logs
```

### 2. Updated `test_CheFSI.py`
*   **Changed Defaults:** `filter_order=5` (was 30) to prevent collapse. `n_iter=50` to compensate.
*   Plots Cond Number and Non-Orthogonality.

```python
import argparse
import numpy as np
import pyopencl as cl
import scipy.linalg
import matplotlib.pyplot as plt
from CheFSI import FrontierSolver

# ... (generate_1d_system function remains the same as before) ...
from CheFSI import FrontierSolver # Re-import to ensure class matches

def generate_1d_system(n_atoms=64, distort_amp=0.1, seed=None):
    if seed is not None: np.random.seed(seed)
    # Positions
    positions = np.arange(n_atoms, dtype=np.float32) * 1.0
    positions += np.random.uniform(-distort_amp, distort_amp, n_atoms).astype(np.float32)
    # Matrices
    H_dense = np.zeros((n_atoms, n_atoms), dtype=np.float32)
    S_dense = np.eye(n_atoms, dtype=np.float32)
    n_max_neighs = 3
    indices = np.full((n_atoms, n_max_neighs), -1, dtype=np.int32)
    H_val = np.zeros((n_atoms, n_max_neighs), dtype=np.float32)
    S_val = np.zeros((n_atoms, n_max_neighs), dtype=np.float32)
    
    for i in range(n_atoms):
        H_dense[i, i] = np.random.uniform(-0.05, 0.05)
        indices[i, 1] = i; H_val[i, 1] = H_dense[i, i]; S_val[i, 1] = 1.0
        
        for idx, j in enumerate([i-1, i+1]):
            if 0 <= j < n_atoms:
                dist = abs(positions[i] - positions[j])
                fac = np.exp(-1.0 * (dist - 1.0)**2)
                h_elem = -1.0 * fac
                s_elem = 0.2 * fac
                
                H_dense[i, j] = h_elem; S_dense[i, j] = s_elem
                col = 0 if idx==0 else 2
                indices[i, col] = j; H_val[i, col] = h_elem; S_val[i, col] = s_elem
    return H_dense, S_dense, indices, H_val, S_val, positions

if __name__ == "__main__":
    # CHANGED DEFAULT: filter_order=5 (prevent collapse), n_iter=40
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-atoms", type=int, default=64)
    parser.add_argument("--n-vecs", type=int, default=8)
    parser.add_argument("--n-iter", type=int, default=40) 
    parser.add_argument("--filter-order", type=int, default=5)
    parser.add_argument("--target", type=float, default=0.5)
    args = parser.parse_args()

    # OpenCL Setup
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    
    # System
    H_dns, S_dns, ell_idx, ell_H, ell_S, pos = generate_1d_system(args.n_atoms)
    
    # Reference
    ref_evals, ref_evecs = scipy.linalg.eigh(H_dns, S_dns)
    
    # Solver
    solver = FrontierSolver(ctx, queue, n_vecs=args.n_vecs, tile_size=32)
    solver.prepare(ell_idx, ell_H, ell_S)
    
    gpu_evals, gpu_vecs_dev, logs = solver.solve(
        target_energy=args.target, 
        n_iter=args.n_iter, 
        filter_order=args.filter_order
    )
    gpu_vecs = gpu_vecs_dev.get()
    
    # Analysis
    print(f"\n{'State':<6} | {'GPU E':<10} | {'Exact E':<10} | {'Overlap':<10} | {'Cond':<10}")
    for i in range(args.n_vecs):
        e_gpu = gpu_evals[i]
        idx_ref = np.argmin(np.abs(ref_evals - e_gpu))
        v = gpu_vecs[:, i]
        v_ref = ref_evecs[:, idx_ref]
        overlap = abs(np.dot(v.T, np.dot(S_dns, v_ref)))
        print(f"{i:<6} | {e_gpu:<10.5f} | {ref_evals[idx_ref]:<10.5f} | {overlap:<10.4f} | -")

    # Plot Convergence
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(logs['max_non_orth'], 'r-o')
    plt.title("Non-Orthogonality (Before Correction)")
    plt.yscale('log')
    plt.xlabel("Iteration")
    
    plt.subplot(1, 2, 2)
    plt.plot(logs['cond'], 'b-o')
    plt.title("Gram Matrix Condition Number")
    plt.yscale('log')
    plt.axhline(1e12, color='k', linestyle='--', label="Unstable")
    plt.xlabel("Iteration")
    plt.legend()
    plt.tight_layout()
    plt.show()
```

---

## USER

OK, after you modification maybe it improved shlightly, but there are still huge differences both in plot and RMSE. Thinkl deeply about the problem. your job is not to write me upodated code, but to analyse the problem and physics and numerical math and tell what is wrong and how to repair it

---

prokophapala@carbsisYoga:~/git/FireCore/pyBall/FireballOCL$ python test_CheFSI_new.py
Choose platform:
[0] <pyopencl.Platform 'Intel(R) OpenCL HD Graphics' at 0x55b247dd6740>
[1] <pyopencl.Platform 'NVIDIA CUDA' at 0x55b2483288d0>
Choice [0]:0
Choosing only available device: <pyopencl.Device 'Intel(R) Graphics [0x46a6]' on 'Intel(R) OpenCL HD Graphics' at 0x55b2482ef580>
Set the environment variable PYOPENCL_CTX='0' to avoid being asked again.
Target Energy: 0.5
Exact spectrum range: [-1.42, 3.30]
Closest Exact Eigenvalues: [0.4669758  0.58081025 0.35531443 0.6901096 ]
[Solver] Compiling kernels with M=8, K=32...

--- Solver Parameters ---
n_atoms=64, n_vecs=8, n_iter=40, filter_order=5, tol=1e-06
target=0.5, bounds=auto, distort_amp=0.1, seed=None, error_scale=1.0
Overlap = |<v_gpu | S | v_ref>| (generalized inner product)
[Solver] Bounds=(np.float32(-2.2388968), np.float32(2.2463994)), Target=0.5, alpha=1.98e-01
Iter 1: E0=0.5230 Gap=3.15e-01 NonOrth=8.51e-01 Cond=7.9e+00
Iter 2: E0=0.5200 Gap=3.20e-01 NonOrth=4.50e-01 Cond=2.0e+00
Iter 3: E0=0.4849 Gap=3.21e-01 NonOrth=2.95e-01 Cond=1.7e+00
Iter 4: E0=0.4620 Gap=3.21e-01 NonOrth=2.67e-01 Cond=1.5e+00
Iter 5: E0=0.4525 Gap=3.19e-01 NonOrth=2.31e-01 Cond=1.5e+00
Iter 6: E0=0.4511 Gap=3.17e-01 NonOrth=2.24e-01 Cond=1.5e+00
Iter 7: E0=0.4527 Gap=3.14e-01 NonOrth=2.09e-01 Cond=1.5e+00
Iter 8: E0=0.5423 Gap=3.11e-01 NonOrth=1.90e-01 Cond=1.5e+00
Iter 9: E0=0.5409 Gap=3.08e-01 NonOrth=1.85e-01 Cond=1.5e+00
Iter 10: E0=0.4578 Gap=3.05e-01 NonOrth=2.85e-01 Cond=1.5e+00
Iter 11: E0=0.4597 Gap=3.03e-01 NonOrth=1.86e-01 Cond=1.5e+00
Iter 12: E0=0.4616 Gap=3.01e-01 NonOrth=1.76e-01 Cond=1.5e+00
Iter 13: E0=0.4630 Gap=2.99e-01 NonOrth=2.04e-01 Cond=1.5e+00
Iter 14: E0=0.4641 Gap=2.98e-01 NonOrth=1.67e-01 Cond=1.5e+00
Iter 15: E0=0.4649 Gap=2.97e-01 NonOrth=1.70e-01 Cond=1.5e+00
Iter 16: E0=0.4655 Gap=2.96e-01 NonOrth=1.71e-01 Cond=1.4e+00
Iter 17: E0=0.4659 Gap=2.97e-01 NonOrth=2.07e-01 Cond=1.4e+00
Iter 18: E0=0.4662 Gap=2.98e-01 NonOrth=1.61e-01 Cond=1.4e+00
Iter 19: E0=0.4665 Gap=3.00e-01 NonOrth=1.72e-01 Cond=1.4e+00
Iter 20: E0=0.4666 Gap=3.03e-01 NonOrth=1.73e-01 Cond=1.4e+00
Iter 21: E0=0.4667 Gap=3.06e-01 NonOrth=1.44e-01 Cond=1.4e+00
Iter 22: E0=0.4668 Gap=3.11e-01 NonOrth=2.18e-01 Cond=1.4e+00
Iter 23: E0=0.4668 Gap=3.15e-01 NonOrth=1.75e-01 Cond=1.4e+00
Iter 24: E0=0.4669 Gap=3.19e-01 NonOrth=2.12e-01 Cond=1.4e+00
Iter 25: E0=0.4669 Gap=3.23e-01 NonOrth=1.34e-01 Cond=1.4e+00
Iter 26: E0=0.4669 Gap=3.27e-01 NonOrth=2.22e-01 Cond=1.4e+00
Iter 27: E0=0.4669 Gap=3.31e-01 NonOrth=1.74e-01 Cond=1.4e+00
Iter 28: E0=0.4670 Gap=3.34e-01 NonOrth=1.73e-01 Cond=1.3e+00
Iter 29: E0=0.4670 Gap=3.36e-01 NonOrth=1.73e-01 Cond=1.3e+00
Iter 30: E0=0.4670 Gap=3.39e-01 NonOrth=1.36e-01 Cond=1.3e+00
Iter 31: E0=0.4670 Gap=3.41e-01 NonOrth=1.72e-01 Cond=1.3e+00
Iter 32: E0=0.4670 Gap=3.42e-01 NonOrth=1.67e-01 Cond=1.3e+00
Iter 33: E0=0.4670 Gap=3.44e-01 NonOrth=2.06e-01 Cond=1.3e+00
Iter 34: E0=0.4670 Gap=3.45e-01 NonOrth=1.61e-01 Cond=1.3e+00
Iter 35: E0=0.4670 Gap=3.46e-01 NonOrth=1.74e-01 Cond=1.3e+00
Iter 36: E0=0.4670 Gap=3.47e-01 NonOrth=1.98e-01 Cond=1.3e+00
Iter 37: E0=0.4670 Gap=3.48e-01 NonOrth=1.83e-01 Cond=1.3e+00
Iter 38: E0=0.4670 Gap=3.48e-01 NonOrth=2.10e-01 Cond=1.3e+00
Iter 39: E0=0.4670 Gap=3.49e-01 NonOrth=1.56e-01 Cond=1.3e+00
Iter 40: E0=0.4670 Gap=3.49e-01 NonOrth=1.65e-01 Cond=1.3e+00

--- Comparison Results ---
Columns: State index | GPU eigenvalue | Exact eigenvalue | Abs diff | Residual ||(H - eS)v|| | Overlap | RMSE vs exact eigenvector
State  | GPU E      | Exact E    | Diff       | Residual   | Overlap    | RMSE      
------------------------------------------------------------------------------------------
0      | 0.46698    | 0.46698    | 4.8e-07    | 2.5e-01    | 0.5069     | 2.284e-01 
1      | 0.58081    | 0.58081    | 6.0e-07    | 2.5e-01    | 0.6363     | 2.375e-01 
2      | 0.35531    | 0.35531    | 8.9e-08    | 2.0e-01    | 0.4624     | 1.353e-01 
3      | 0.69011    | 0.69011    | 1.8e-07    | 2.2e-01    | 0.0898     | 1.808e-01 
4      | 0.24831    | 0.24828    | 2.7e-05    | 3.0e-01    | 0.0032     | 1.827e-01 
5      | 0.81465    | 0.81467    | 1.1e-05    | 3.9e-01    | 0.3324     | 1.523e-01 
6      | 0.15105    | 0.14927    | 1.8e-03    | 3.0e-01    | 0.4613     | 2.186e-01 
7      | 0.91989    | 0.93472    | 1.5e-02    | 3.1e-01    | 0.4847     | 1.358e-01 

Max non-orth logs (len=40): last=1.653e-01
Total SpMM calls (matvecs): 881
Mean orbital RMSE: 1.839e-01, max RMSE: 2.375e-01
prokophapala@carbsisYoga:~/git/FireCore/pyBall/FireballOCL$ 

---

## Gemini

This is a classic "Right Answer, Wrong Vector" problem, and I have identified the root cause. It is a memory layout issue in the interface between NumPy (LAPACK) and OpenCL.

### The Root Cause: Transposed Rotation Matrix

The problem lies in **Data Contiguity**.
1.  **NumPy / LAPACK:** `scipy.linalg.eigh` (and `numpy.linalg.eigh`) typically return the eigenvector matrix `evecs` in **Fortran Order** (Column-Major) because LAPACK is written in Fortran.
2.  **OpenCL Kernel:** Your `subspace_rotate` kernel assumes the matrix `M` is in **C Order** (Row-Major): `M[i * n_vecs + j]`.

When you simply call `cl_array.to_device(..., evecs)`, PyOpenCL copies the memory as-is.
*   If `evecs` is Column-Major, the element $(i, j)$ is stored where $(j, i)$ would be in Row-Major.
*   Effectively, your kernel reads the **Transpose** ($U^T$) instead of $U$.

### Why this explains all symptoms:
1.  **Correct Energies:** The subspace diagonalization finds the correct eigenvalues. The rotation matrix $U$ is supposed to align your basis vectors $X$ with these eigenvalues.
2.  **High Residuals:** Since you applied $U^T$ instead of $U$, your resulting vectors $X_{new}$ are **linear combinations** of the true eigenvectors. A linear combination of eigenvectors with *different* energies is **not** an eigenvector. Thus, $(H - \epsilon S)v \neq 0$.
3.  **Orthogonality Passes:** Since $U$ is orthogonal ($U^T U = I$), $U^T$ is also orthogonal. Rotating an orthogonal basis by an orthogonal matrix preserves orthogonality. Thus `NonOrth` is low.
4.  **Bad Overlap:** Your vectors are scrambled rotations of the true vectors, so they only partially overlap with the reference.

### The Fix

You must verify the memory layout of the eigenvectors and force them to be C-contiguous before sending them to the GPU.

In `CheFSI.py`, find the `solve` and `_orthogonalize` methods where `d_rot` is created. Change:

```python
# OLD (Buggy)
d_rot = cl_array.to_device(self.queue, evecs.astype(np.float32))

# NEW (Fixed)
# Force C-contiguous layout
d_rot = cl_array.to_device(self.queue, np.ascontiguousarray(evecs.astype(np.float32)))
```

Apply this fix to both instances (inside `solve` loop and inside `_orthogonalize`).

### Summary of Physics/Math Analysis
You successfully built a Subspace Iteration solver.
*   **The Physics (Hamiltonian):** Was correct (proven by matching Eigenvalues).
*   **The Math (Solver):** Was correct (proven by convergence of subspace energy).
*   **The Implementation:** Failed on the `Fortran vs C` memory layout, causing the final step (Ritz rotation) to apply the wrong mixing coefficients, leaving you with a subspace that *contained* the answer but wasn't *aligned* to it.


---

## USER

I realized one problem, there is some degeneracy. we assign the results to the nearest true eigenvector, and we realize that two GPU apporox vectro assing to same true eigenvectr. This degeenracy, perhaps, cause problem with the convergnece? as we try to orthogonalize something wich try to map to same vectros, or we try to expand space which has inherently lower dimension. Can we do something abut that?

--- Comparison Results ---
Columns: GPU state | Ref idx | GPU eigenvalue | Exact eigenvalue | Abs diff | Residual ||(H - eS)v|| | Overlap | RMSE vs exact eigenvector
GPU    | Ref    | GPU E      | Exact E    | Diff       | Residual   | Overlap    | RMSE      
----------------------------------------------------------------------------------------------------
0      | 36     | 0.47151    | 0.47583    | 4.3e-03    | 1.2e-01    | 0.9619     | 2.590e-01 
1      | 37     | 0.53299    | 0.58735    | 5.4e-02    | 4.1e-01    | 0.1101     | 1.758e-01 
2      | 37     | 0.58802    | 0.58735    | 6.8e-04    | 4.5e-02    | 0.9937     | 1.495e-02 
3      | 35     | 0.35823    | 0.35894    | 7.1e-04    | 2.8e-02    | 0.9972     | 9.512e-03 
4      | 38     | 0.70520    | 0.70037    | 4.8e-03    | 6.3e-02    | 0.9798     | 2.711e-02 
5      | 34     | 0.25597    | 0.25620    | 2.3e-04    | 1.1e-02    | 0.9992     | 5.183e-03 
6      | 39     | 0.83084    | 0.82672    | 4.1e-03    | 3.5e-02    | 0.9785     | 2.817e-02 
7      | 33     | 0.15405    | 0.15100    | 3.1e-03    | 6.0e-02    | 0.9899     | 2.531e-01 

---


## Gemini

Yes, you have correctly identified a fundamental limitation of subspace iteration methods (like CheFSI, LOBPCG, or Davidson).

This phenomenon is called **"Spectral Crowding"** or the need for **"Guard Vectors"**.

### The Problem
You asked for 8 eigenvectors. The solver allocates a subspace of dimension 8.
1.  **The Edges are Weak:** In any iterative method, the vectors at the "edges" of your subspace converge much slower than the ones in the middle.
2.  **Ghost States:** If one vector (e.g., Vector #1) hasn't converged yet, it is essentially "noise". When you orthogonalize, this noise might accidentally overlap with a state that Vector #2 is already optimizing for.
3.  **The Result:** You get 7 good vectors and 1 garbage vector (GPU state 1), or two vectors trying to describe the same state (degeneracy), which ruins the orthogonality stability.

### The Solution: Oversampling (Guard Vectors)
If you want **8** high-quality eigenvectors, you must solve for **12 or 16**.
You throw away the extra 4-8 vectors at the end. These "Guard Vectors" absorb the noise and the slow convergence, protecting your inner 8 target vectors.

### Comparison Logic Fix
Your comparison script is using `argmin` naively, which allows multiple GPU vectors to claim the same Reference vector. We need a **1-to-1 matching** (like the Hungarian algorithm or a greedy match) to see what is really happening.

---

### Step 1: Update `test_CheFSI.py` to use Guard Vectors

Don't change the kernel. Just change how you call the solver.

**Changes to make in `test_CheFSI.py`:**
1.  Define `n_extra = 4` (Guard vectors).
2.  Tell the solver to solve for `n_vecs + n_extra`.
3.  Slice the results to keep only the best `n_vecs`.
4.  Use **Greedy Matching** for the table to prevent duplicate reference indices.

Here is the corrected section for `test_CheFSI.py`:

```python
    # ... inside main ...

    # 4. Run GPU Solver with OVERSAMPLING
    print("\n--- Running GPU FrontierSolver ---")
    N_WANTED = args.n_vecs     # We want 8
    N_GUARD = 4                # Extra buffer vectors
    N_SOLVE = N_WANTED + N_GUARD 

    print(f"Solving for {N_SOLVE} vectors (keeping best {N_WANTED})...")

    # Initialize solver with larger size
    solver = FrontierSolver(ctx, queue, n_vecs=N_SOLVE, tile_size=32, verbose=True)
    solver.prepare(ell_idx, ell_H, ell_S)
    
    # Run Solve
    gpu_evals_all, gpu_vecs_dev_all, logs = solver.solve(
        target_energy=TARGET_E, 
        n_iter=args.n_iter, 
        filter_order=args.filter_order
    )
    
    # Get all vectors to host
    gpu_vecs_all = gpu_vecs_dev_all.get()

    # 5. Intelligent Filtering and Sorting
    # We computed N_SOLVE vectors. We need to pick the N_WANTED best ones.
    # Criterion: Residual ||(H-eS)v|| is usually the best metric for quality.
    
    residuals = []
    for i in range(N_SOLVE):
        v = gpu_vecs_all[:, i]
        e = gpu_evals_all[i]
        # Calculate residual norm
        # Note: Doing this on CPU for checking. 
        # In production, trust eigenvalues or run residual check on GPU.
        res_vec = (H_dns @ v) - e * (S_dns @ v)
        residuals.append(np.linalg.norm(res_vec))
        
    residuals = np.array(residuals)
    
    # Sort by Residual Quality (Ascending)
    # Alternatively, you can sort by distance to Target Energy if you prefer.
    # best_indices = np.argsort(residuals)[:N_WANTED] 
    
    # Let's stick to energy sorting (standard), but just be aware the last ones are guards.
    # The solver typically returns them sorted by energy/distance already.
    # We just chop off the last N_GUARD vectors.
    
    # However, let's explicitly sort by distance to target to be safe
    dist_to_target = np.abs(gpu_evals_all - TARGET_E)
    sorted_idx = np.argsort(dist_to_target)
    
    # Keep only the requested number
    keep_idx = sorted_idx[:N_WANTED]
    
    gpu_evals = gpu_evals_all[keep_idx]
    gpu_vecs = gpu_vecs_all[:, keep_idx]

    # 6. robust 1-to-1 Comparison (Greedy Matching)
    print("\n--- Comparison Results (Greedy 1-to-1 Matching) ---")
    print(f"{'State':<6} | {'Ref Idx':<8} | {'GPU E':<10} | {'Exact E':<10} | {'Diff':<10} | {'Overlap':<10} | {'RMSE':<10}")
    print("-" * 85)
    
    # Create a cost matrix: Abs(GPU_E - Ref_E)
    # Rows = GPU states, Cols = Ref states
    cost_matrix = np.abs(gpu_evals[:, None] - ref_evals[None, :])
    
    # Greedy assignment
    assigned_ref_indices = set()
    
    for i in range(N_WANTED):
        e_gpu = gpu_evals[i]
        
        # Find best match for this GPU state that hasn't been used yet
        costs = cost_matrix[i, :]
        
        # Sort possible matches by cost
        possible_matches = np.argsort(costs)
        
        best_match_idx = -1
        for idx in possible_matches:
            if idx not in assigned_ref_indices:
                best_match_idx = idx
                break
        
        # Fallback if everything reasonable is taken (shouldn't happen with 1-to-1)
        if best_match_idx == -1:
            best_match_idx = possible_matches[0] # Duplicate (sad)
            print("Warning: forced duplicate assignment")
        else:
            assigned_ref_indices.add(best_match_idx)
            
        e_ref = ref_evals[best_match_idx]
        diff = abs(e_gpu - e_ref)
        
        v = gpu_vecs[:, i]
        v_ref = ref_evecs[:, best_match_idx]
        overlap = abs(np.dot(v.T, np.dot(S_dns, v_ref)))
        
        # Phase alignment for RMSE
        phase = np.sign(np.dot(v, v_ref))
        if phase == 0: phase = 1
        rmse = np.sqrt(np.mean((v * phase - v_ref)**2))
        
        print(f"{i:<6} | {best_match_idx:<8} | {e_gpu:<10.5f} | {e_ref:<10.5f} | {diff:<10.1e} | {overlap:<10.4f} | {rmse:<10.3e}")

    # ... plotting code ...
```

### Why this fixes it:

1.  **Guard Vectors:** By asking for 12 vectors, the solver struggles to converge vectors #11 and #12. But vectors #1 through #8 are now "deep inside" the subspace and converge rapidly.
2.  **Greedy Matching:** The script will no longer confusingly tell you that GPU 1 matches Ref 37 if Ref 37 is already the perfect match for GPU 2. It will find the *next* best match for GPU 1, which will reveal if GPU 1 is actually converging to a different state or is just pure noise.

Try running with `n_vecs=8` and `n_extra=4` (so `n_solve=12`). You should see overlaps > 0.99 for all 8 reported vectors.

---


## USER

OK, so you say addin spurious vectros help stability not degrade it? I was thinking instead that we should remove one vector (continue iteration only with 7) to prevent this "noise" or "garbage" trash the convergence.

---


## Gemini

You are asking a very deep question that touches on the fundamental mechanics of iterative linear algebra (Krylov and Subspace methods).

Your intuition—*"remove the bad vector to save the rest"*—is logical for biological or social systems (quarantine the sick), but **it is incorrect for linear algebra subspaces.**

Here is why **adding** "spurious" (guard) vectors improves stability, and why removing them hurts you.

### 1. The "Mudroom" Analogy (The Buffer Zone)

Imagine you want to keep your living room perfectly clean (your 8 eigenvectors), but it's muddy outside (the rest of the Hilbert space, noise).

*   **Your Proposal (Exact size):** You open the front door directly into the living room. You try to stop the mud exactly at the threshold.
    *   *Result:* The area near the door (Vector #8) always gets dirty because it is in direct contact with the outside world. When you sweep (orthogonalize/rotate), you spread that dirt around the whole room.
*   **Guard Vectors (Oversampling):** You build a small **mudroom** or **foyer** (Vectors #9, #10, #11, #12).
    *   *Result:* The mud gets into the mudroom. Vectors #11 and #12 are dirty and noisy. But the living room (Vectors #1-8) is separated from the outside by the mudroom.
    *   At the end, you lock the door between the living room and the mudroom, and you only show your guests the living room.

### 2. The Math: Spectral Gaps and Convergence

The convergence speed of a subspace iteration depends on the ratio of eigenvalues at the boundary of your subspace.

If you iterate $N$ vectors, the error $\epsilon$ decreases as:
$$ \epsilon \propto \left( \frac{\lambda_{\text{inside}}}{\lambda_{\text{outside}}} \right)^k $$

*   **Case A: 8 Vectors.**
    *   Boundary is between Vector 8 and Vector 9 (which is not in your subspace).
    *   If $\lambda_8 \approx 0.47$ and $\lambda_9 \approx 0.48$, the ratio is $\approx 0.98$.
    *   Convergence is **extremely slow**. Vector 8 struggles to decide if it belongs inside or outside. It fluctuates, causing "noise".

*   **Case B: 12 Vectors (8 Target + 4 Guards).**
    *   The boundary is now between Vector 12 and Vector 13.
    *   Vectors 1 through 8 are now "deep" inside the subspace.
    *   Vector 8 is comparing itself to Vector 12. The spectral gap is huge.
    *   Vector 8 converges **exponentially faster** because it is far away from the noisy boundary.

### 3. Rayleigh-Ritz is a "Cleaning Machine"

You are using the Rayleigh-Ritz step (diagonalize $H_{sub} = X^T H X$) in every iteration. This is a powerful sorting mechanism.

When you have 12 vectors:
1.  The subspace captures a mix of the target states and some noise.
2.  The Rayleigh-Ritz step rotates the basis. It automatically pushes the **best, most convergent energy components** into slots 0-7.
3.  It pushes the **residual error, noise, and unconverged high-energy components** into slots 8-11.
4.  Slots 8-11 act as "trash cans" that collect the numerical garbage, keeping slots 0-7 pristine.

If you only have 8 vectors, the "trash" has nowhere to go. It stays in vector 7, and due to orthogonality constraints, it pushes vector 6, which pushes vector 5... destabilizing the whole block.

### Summary
*   **Do not remove vectors.** If you remove the noisy vector #8, then vector #7 becomes the new boundary, and *it* will become noisy. You just retreat until you have nothing left.
*   **Add Guard Vectors.** Standard practice in LOBPCG and CheFSI is to use **20% to 50% extra vectors**.
    *   If you want 8, solve for 12.
    *   If you want 100, solve for 120.

**Action for your code:**
Keep the code exactly as I corrected in the previous turn (solving for `N_WANTED + N_GUARD`), and you will see the 8 target vectors stabilize almost perfectly.