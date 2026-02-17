# USER

I have pyOpenCL Monte-Carlo solver of Hubard/Ising model searchin minimum energy occupancy of array of quantum dot depening of position of STM tip which gate the on site energies by bias voltage (these are not spins but charge states).
The problem is that this is just ground-state solver. But In experiments we see that shard ptansitions beweeen states are non-equalibirum and depend on hopping from tip, and we must solve it at least on leel of Pauli Master Quation. The problem is that number of many body states scales exponentially with number of sites ( while we can solve array of 128 sites easily by monte carlo in parallel for many pixels=positins of STM tip, we can solve max 6-8 sites by brute force Pauli-Master Equations ).
So I was thinking to connect the two but only in some sub-space of lowest energy many body states. We fisrts find ground-state by monte carlo, then construct some low energy excitations (and that is the main think we need to think deeply about, how to search these low energy excitations, what should be the criterium - energy difference or hamming distance = number of charge flips=bit-flips in occupancy vector? perhaps both combined)
Then once we have the sub space (only few states per each site e.g. max 64x64 matrix, maybe even less) we should solve theis Pauli Kernel Matrix by some suitabl e method. Not sure if best method is Gauss-Elimination, if we cannot find some better more sutiable method  - e.g. iterative to use on GPU? that each pixel (each position of tip) should be solved by one workgroup on GPU.

---

# Gemini-3-pro

This is a sophisticated problem connecting combinatorial optimization (Ground State search) with non-equilibrium statistical mechanics (Pauli Master Equation) on a GPU.

Here is the strategy breakdown regarding the Physics, the Sub-space selection, and the Solver mathematics.

### 1. Physics & Topology of the Sub-Space

The fundamental issue with simple "Low Energy" selection is **Connectivity**.

In the Pauli Master Equation (PME), transitions usually represent single-electron processes (tunneling from tip/substrate). This means transitions are valid only between states with a **Hamming distance of 1** (one charge flip).

*   **The Trap of Energy Selection:** If you simply pick the 64 lowest energy states from the $2^{128}$ Hilbert space, you might find the Ground State (GS) and a cluster of low-energy states that are, say, 4 bit-flips away (Hamming distance 4). If you do not include the intermediate states (dist 1, 2, and 3), the GS is **kinetically disconnected** from those excitations. The matrix will block-diagonalize, and the probability flow will be zero (or incorrect).
*   **The "Star" Topology:** The most physically relevant excitations for transport are those directly connected to the Ground State. These are the **Single-Flip Excitations**.
    *   There are exactly $N$ such states (where $N=128$).
    *   If you restrict your basis to $\{GS\} \cup \{ \text{Single Flips} \}$, the connectivity graph is a **Star Graph**. The GS is the hub; the excitations are the spokes. Excitations cannot talk to each other directly (that would be distance 2).
    *   **Analytical Shortcut:** A Star Graph PME can be solved **analytically** in $O(N)$ time. You do not need Gaussian elimination. The steady-state probability of an excitation $P_i$ is simply related to $P_{GS}$ by the ratio of rates $W_{GS \to i} / W_{i \to GS}$.

**Recommendation:**
Since you want to solve this for $N=128$ sites, but can only afford a matrix size of $\sim 64$, you cannot even include *all* single-flip excitations.
You should construct the basis by:
1.  Taking the **Ground State (GS)** found by Monte Carlo.
2.  Scanning all 128 possible **Single Bit Flips**.
3.  Calculating the Energy Cost $\Delta E$ for each flip.
4.  Selecting the top $K$ (e.g., 63) flips with the **lowest excitation energy**.
5.  (Optional) If you suspect strong second-order effects (cotunneling via virtual states), you might include "Double Flips", but this explodes the count. For a first pass, the "Lowest Energy Connected Neighbors" is the robust choice.

### 2. The Numerical Solver Strategy (GPU)

You asked about Iterative vs. Gaussian Elimination.

**Verdict: Direct Gaussian Elimination (Gauss-Jordan) in Shared Memory.**

**Why not Iterative (Arnoldi/Lanczos)?**
1.  **Stiffness:** Tunneling rates scale exponentially with distance ($e^{-\kappa r}$). This leads to a rate matrix where elements vary by orders of magnitude (e.g., $10^9$ Hz vs $10^3$ Hz). Such matrices are extremely "stiff." Iterative solvers converge very slowly on stiff matrices unless preconditioned perfectly.
2.  **Divergence:** Simple iterative maps ($P_{t+1} = P_t + \Delta t K P_t$) require extremely small time steps to remain stable for the fast rates, making them useless for finding the steady state of the slow rates.

**Why Gaussian Elimination?**
1.  **Size:** A $64 \times 64$ float matrix takes 16KB. This fits perfectly into the 32KB-64KB Local Data Share (LDS) / Shared Memory of a standard GPU workgroup.
2.  **Predictability:** The runtime is deterministic $O(M^3)$. For $M=64$, this is manageable per pixel.
3.  **Stability:** With Partial Pivoting (swapping rows), it is robust enough for single-precision floats, even with stiff rates.

### 3. Proposed Pipeline

Do not try to fuse everything into one kernel. Split it into logical steps.

#### Step A: The "Basis Scanner" Kernel
*   **Input:** The `occ_out` (GS bitmask) from your Monte Carlo solver.
*   **Operation:**
    *   Each workgroup (or thread) loads the GS for one tip position.
    *   Iterate through all sites $i \in [0, 128)$.
    *   Calculate $\Delta E_i$: The energy change if site $i$ is flipped.
        *   This is efficient: $\Delta E_i = \pm (E_{site, i} + \sum_{j \in GS} W_{ij})$.
        *   You only need the row $i$ of the sparse $W$ matrix.
    *   **Sort:** Store pairs `(index, delta_energy)` and sort them.
    *   **Select:** Keep the best $M-1$ indices (e.g., top 63).
*   **Output:** A buffer `basis_indices[nTips, M]` containing the bit-indices of the chosen subspace. Index 0 can be a sentinel for "No Flip (GS)".

#### Step B: The "Subspace Solver" Kernel
*   **Input:** `occ_out` (GS), `basis_indices`, Physics Params.
*   **Operation (One Workgroup per Tip):**
    1.  **Load Basis:** Load the $M$ states into Local Memory.
        *   State 0: `GS`
        *   State $k$: `GS` with bit `basis_indices[k]` flipped.
    2.  **Build Matrix (in Registers/LDS):**
        *   Loop $a = 0..M$, $b = 0..M$.
        *   Compute Hamming distance between Basis[$a$] and Basis[$b$].
        *   **Crucial Logic:**
            *   If $a=0$ (GS) and $b>0$ (Excitation): Dist is always 1. Calculate Rate $GS \to Exc$.
            *   If $a>0$ and $b=0$: Dist is always 1. Calculate Rate $Exc \to GS$.
            *   If $a>0$ and $b>0$: The distance is 2 (since they are different single flips from the same GS). **Rate is 0**.
        *   *Note:* This confirms the **Star Graph** topology. The matrix is exceedingly sparse (only first row and first column are non-zero, plus diagonal).
    3.  **Solve:**
        *   **If Star Graph:** You don't need a matrix solver.
            *   $P_{GS} \sum_{i} \Gamma_{GS \to i} = \sum_{i} P_i \Gamma_{i \to GS}$
            *   $P_i = P_{GS} \frac{\Gamma_{GS \to i}}{\Gamma_{i \to GS}}$
            *   normalization: $1 = P_{GS} (1 + \sum \frac{\Gamma_{GS \to i}}{\Gamma_{i \to GS}})$
            *   This is $O(M)$. Extremely fast.
        *   **If you include 2-flips (Dense Matrix):**
            *   Fill the full $M \times M$ matrix in LDS.
            *   Run the parallel Gauss-Jordan elimination (as in your existing code).

### 4. Mathematical Refinement for the Basis

To ensure the "Star Graph" assumption holds (and thus extreme speed), you must verify that direct transitions between excitations are negligible.
*   Excitation A: GS + electron at site $i$.
*   Excitation B: GS + electron at site $j$.
*   Transition A $\to$ B requires moving electron $i \to j$.
*   In your model (Hubbard with Tip), is there direct hopping $t_{ij}$?
    *   **If No (only Tip-Site tunneling):** Then A $\to$ B is impossible (requires 2 steps: $i \to Tip \to j$). The Star Graph is exact.
    *   **If Yes (inter-site hopping):** Then A and B are connected. The matrix is not a Star; it is a "Wheel" (Hub + rim connections). You **must** use the Gauss-Jordan solver.

Given you mentioned "hopping from tip", I assume the dominant transport is Tip-driven. However, if lateral hopping exists, use the Gauss-Jordan solver on the sub-space.

### 5. Summary of the Plan

1.  **Use MC** to get the global minimum (Ground State).
2.  **Use a Scanner Kernel** to find the $K$ lowest energy *single-charge excitations* relative to that Ground State.
3.  **Construct a reduced Hamiltonian/Rate Matrix** for just those $K+1$ states.
    *   State definitions are stored simply as "Index of flipped bit" relative to GS.
4.  **Solve** using:
    *   **Analytic Sum** if lateral hopping $t_{ij} \approx 0$.
    *   **LDS Gauss-Jordan** if lateral hopping is significant.

This approach scales linearly with $N$ (for the scan) and effectively solves the exponential bottleneck by physically truncating the Hilbert space to the "transport active" shell.

---

# USER

OOD,
I agree we should consider only single-flip excitations because we consider only single electron hopping (no cotuneling etc, only 1st order perturbative approach, that has highest amplitudes), this makes it quite easy. Architectural question is how to split this?
i) the MC groundstate search kernel should be definitely separate, this is the most costly part and we should not complromise it (it can take hudersts or thousands of iterations (kernel passes) to find groud state)
ii) then we can directly generate the single flip excitations, evaluate their energies, sort them, choose the top N, and directly build and solve the PME kernel matrix, or we can split that into 2-3 kernels. How? Lets think about it before we start coding. We should consider registre pressure. etc. Although I think that sorting and solving PME will be much less demanding than searching the ground state.
We should assume that number of single flip excited many body states for one pixel always fit to local memory (max 128 states?) and that number of selected low energy sub-space is smaller then number of threads in the workgroups (nManyBodySubSpace<workgroupsize), so we can sort them efficiently in local memory.
I agree, when doing Gauss-Jordan elimination in workgroup, worgroup size (number of threads) and the dimension of matrix is the same, I can do one line elimination in parallel right? Please write this elemination in parallel way within workgroup. But as we reduce it to trinagular shape the treads start to be idle (right)? Aren't there some other more efficient direct solver methods which can be efficiently paralleized so that all threads are utilized equaly
Note any other consideration, maybe I forgot about something. Either about some problem, or about some opportunity, some method I was not thinking about and can help us.
Now do the analysis we are still not coding.

---

# Gemini-3-pro

This is a very robust plan. Splitting the **"Heavy" MC Optimizer** from the **"Light" PME Analyzer** is the correct architectural choice. The MC kernel is stochastic and iterative; the PME kernel is deterministic and linear algebra-heavy. They have different convergence criteria and thread divergence behaviors.

Here is the breakdown of the strategies and the math for the Solver and the Architecture.

### 1. Architecture: The "Post-Processing" Kernel

We should treat the Pauli-Master Equation (PME) solver as a **Post-Processing Kernel**. It takes the result of the MC (the ground state configuration) and refines it.

**Kernel Split:**
1.  **MC Kernel:** (Existing) Outputs `occ_out` (best bitmask) and `E_out`.
2.  **PME Analyzer Kernel:** (New)
    *   Reads `occ_out`.
    *   Generates Excitations -> Sorts -> Builds Matrix -> Solves -> Writes Current.
    *   **Recommendation:** Do *not* split this into "Generate" and "Solve" kernels. The data passed between them (the sorted list of subspace indices and the matrix) is small enough to stay in Local Memory (LDS) but awkward to serialize to global VRAM. Doing it in one fused kernel keeps everything in the fastest memory tier.

**Memory & Constraints:**
*   **Subspace Size ($M$):** You mentioned solving a $64 \times 64$ matrix. This is the sweet spot.
    *   Matrix Storage: $64 \times 64 \times 4 \text{ bytes (float)} = 16 \text{ KB}$.
    *   This fits comfortably in the 32KB or 64KB Local Memory of modern GPUs.
    *   Workgroup Size: 64 threads. One thread per matrix row/column.
*   **Input Size:** 128 sites.
    *   The "Scanner" phase needs to iterate 128 times, but this is fast.

**Proposed PME Kernel Flow (per Workgroup/Pixel):**
1.  **Load:** Load the Ground State (GS) bitmask from Global to Local.
2.  **Scan & Assess:**
    *   Threads collaborate to evaluate the energy of flipping each of the 128 sites.
    *   Math: $\Delta E_i = \pm (\epsilon_i + V_{tip} + \sum_{j} W_{ij} n_j)$.
    *   Store `(Energy, SiteIndex)` tuples in a local buffer.
3.  **Parallel Sort:**
    *   Sort these 128 tuples by energy.
    *   Keep the top $M-1$ (e.g., 63) single-flip excitations.
    *   Index 0 is reserved for the original GS.
4.  **Build Matrix:**
    *   Populate the $M \times M$ matrix in LDS.
5.  **Solve:**
    *   In-place Gauss-Jordan Elimination.
6.  **Compute Current:**
    *   Dot product of probability vector and rate vector.

### 2. The Parallel Solver: Full Gauss-Jordan

You asked about the efficiency of elimination and idle threads.

**The "Triangular" Problem:**
In standard Gaussian Elimination (LU decomposition), you reduce the matrix to Upper Triangular.
*   Step 1: Eliminate col 0 from rows 1..63.
*   Step 50: Eliminate col 50 from rows 51..63. (Most threads idle).
*   Then you need Back-Substitution (serial bottleneck).

**The Solution: Gauss-Jordan Elimination (Diagonalization)**
Instead of stopping at Upper Triangular, we eliminate **above and below** the pivot at every step. This transforms the matrix directly into the Identity matrix.

*   **Algorithm:**
    For each pivot column $k$ from $0$ to $M-1$:
    1.  **Pivot Selection:** (Thread 0 or reduction) Find row with max value in column $k$, swap with row $k$.
    2.  **Normalization:** Thread $k$ (or all threads) divides Row $k$ by $A[k][k]$. Now $A[k][k] = 1$.
    3.  **Elimination (The Parallel Magic):**
        *   Every Thread $i$ (where $i \neq k$) is responsible for **Row $i$**.
        *   It reads the value $val = A[i][k]$.
        *   It subtracts $val \times \text{Row}[k]$ from $\text{Row}[i]$.
        *   **Result:** Column $k$ becomes 0 for everyone except the pivot row.

*   **Why this is better:**
    *   At step $k$, you update $M-1$ rows.
    *   Utilization is constant: $(M-1)/M \approx 100\%$ throughout the entire algorithm.
    *   No Back-Substitution phase is required. The RHS column *is* the solution vector at the end.

### 3. Missing Considerations & Opportunities

Here are three critical physics/math points that could break the solver if ignored.

#### A. The "Disconnected" Subspace Risk (Grand Canonical vs Canonical)
In your MC solver, you likely flip bits. This changes the total charge $N$.
In the PME, you are coupling to a lead (Tip/Substrate).
*   **The Check:** Does your model allow $N$ to change?
    *   YES: Tunneling events ($GS \leftrightarrow GS \pm 1e$) are valid.
    *   **The Trap:** If you select 63 excitations, ensure they are kinetically connected.
    *   If you pick a state that requires **two** simultaneous tunneling events to reach from the GS (a 2nd order process), and your Hamiltonian only has 1st order tunneling terms, that state will be an "Island" in the matrix.
    *   **Solution:** The "Single Flip" constraint guarantees connectivity to the GS. Do *not* include double-flips unless you add cotunneling terms to the matrix.

#### B. Lateral Hopping vs. Star Topology (The Huge Optimization)
We need to clarify the Hamiltonian term $\sum t_{ij} c^\dagger_i c_j$.
*   **Scenario A: No Lateral Hopping ($t_{ij} = 0$).**
    *   Electrons only tunnel Tip $\leftrightarrow$ Site.
    *   Excitations (Site $i$ flip) can only relax back to GS. They cannot jump to Site $j$.
    *   **Topology:** Star Graph.
    *   **Opportunity:** You do **not** need Gauss-Jordan. The solution is analytic.
        $$ P_{exc, i} = P_{GS} \frac{\Gamma_{GS \to i}}{\Gamma_{i \to GS}} $$
        $$ P_{GS} = \frac{1}{1 + \sum P_{exc, i}/P_{GS}} $$
    *   Speedup: 100x.

*   **Scenario B: Lateral Hopping Exists ($t_{ij} \neq 0$).**
    *   An electron can tunnel onto Site $i$ (State A), then hop to Site $j$ (State B).
    *   This connects Excitation A and Excitation B directly.
    *   **Topology:** Connected Graph.
    *   **Action:** You **must** use the Gauss-Jordan solver.

*   **Strategy:** Write the kernel with the Gauss-Jordan solver. It works for both. If you confirm $t_{ij}==0$ later, you can replace the solver function with the analytic one easily.

#### C. The Sorting Bottleneck & Register Pressure
Sorting 128 elements in a workgroup is fast, but it consumes registers if you aren't careful.
*   **Problem:** If each thread holds 2 elements and you try a bitonic sort, register usage spikes, potentially limiting occupancy.
*   **Optimized Selection:** You don't actually need a *full sort*. You only need the *top K*.
*   **Algorithm:**
    *   Each thread computes energies for 2 sites (since $128 \text{ sites} / 64 \text{ threads} = 2$).
    *   Thread keeps its "Best of 2".
    *   Perform a "Parallel Reduction" style collection to gather the top 64 into the matrix indices.
    *   *Actually, for K=64 and Threads=64, it's trivial:*
        1.  Calculate Energy for Site $i$ and Site $i+64$.
        2.  Store `(E, index)` in a local array `buffer[128]`.
        3.  Barrier.
        4.  One-pass parallel Bitonic Sort (or even 64 passes of Bubble sort, since data is in LDS and $N$ is tiny, overhead is negligible compared to matrix solve).

### Summary of Analysis

1.  **Workgroup:** 1 Workgroup per Pixel (Tip Position).
2.  **Thread Count:** 64 threads (ideal for Wave64 on AMD or 2x Warps on NV).
3.  **Matrix Size:** $64 \times 64$ (Ground State + 63 Single-Flip Excitations).
4.  **Solver:** Full Gauss-Jordan Elimination (diagonalization) to maximize parallel efficiency.
5.  **Data Flow:**
    *   Read `GS_occupancy` (global).
    *   Compute 128 $\Delta E$ values.
    *   Sort/Select top 63.
    *   Fill LDS Matrix $M_{64 \times 64}$.
    *   Solve $M \cdot P = 0$.
    *   Compute $I = \sum P_i \cdot \Gamma_{current}$.
    *   Write `Current` (global).

This looks like a solid plan. We are ready to approach the code structure.

---

# USER

Wow, you say that if I'm limited to single electron tuneling (ingle bit flip) then I can solve it analytically without kernel matrix, I did not noticed that before. That is very importaint if true. Elaborate on that. Please make didactic explanation - like for student, and also grigrous derivation why/how 

---

# Gemini-3-pro

### The "Star Graph" Solution

This is a massive optimization opportunity. If your physics allows it, you turn an $O(N^3)$ matrix inversion problem into an $O(N)$ summation problem.

Here is the explanation, starting with the intuitive picture and moving to the rigorous math.

### 1. The Didactic Explanation (The "Hub & Spoke" Analogy)

Imagine a train system (Probability Flow) with a central station and several dead-end suburbs.

*   **The Central Hub (GS):** This is your Ground State (say, empty dots).
*   **The Suburbs (Excitations):** These are the states with 1 electron added to site $A$, site $B$, site $C$, etc.

**The Rules of the Tracks:**
1.  **Hub $\leftrightarrow$ Suburb:** You can travel between the Hub and any Suburb (tunneling from Tip/Substrate to a Site).
2.  **No Inter-Suburb Tracks:** You **cannot** travel directly from Suburb A to Suburb B. To go from A to B, you must go back to the Hub first.
    *   *Physics:* This means an electron cannot hop directly from Site A to Site B. It must leave Site A (going back to the Lead/GS) and then a new electron enters Site B.

**The Flow:**
In a steady state (equilibrium or non-equilibrium), the number of trains entering a specific Suburb must equal the number of trains leaving it. Since there is only one track connected to Suburb A (the one to the Hub), the traffic is simple:

$$ \text{Flow}(Hub \to A) = \text{Flow}(A \to Hub) $$

If we know the transition rates (the speed of the trains), we can instantly find the population of Suburb A relative to the Hub without caring about Suburb B or C.

**The Conclusion:**
The suburbs are **independent** of each other. You don't need to solve a complex system of equations where A affects B. You just solve A vs Hub, B vs Hub, etc., and then normalize the total population to 100%.

### 2. Rigorous Derivation

Let $P_0$ be the probability of the Ground State.
Let $P_i$ be the probability of the $i$-th Single-Flip Excitation ($i = 1 \dots N$).

#### The Master Equation
The rate of change for probability $P_i$ is the sum of flows in minus flows out:
$$ \frac{dP_i}{dt} = \sum_{j} (W_{j \to i} P_j - W_{i \to j} P_i) $$

In the **Steady State**, $\frac{dP_i}{dt} = 0$.

#### The Topology Constraint (The "Star")
If we allow **only** single-electron tunneling from the leads, and **no** direct hopping between sites ($t_{ij} = 0$), the transition matrix $W$ becomes extremely sparse.

1.  Transitions $0 \leftrightarrow i$ exist (Hub to Spoke).
2.  Transitions $i \leftrightarrow j$ are **zero** (No Spoke to Spoke).

#### The Algebraic Solution
Look at the equation for an excited state $i$:
$$ 0 = \sum_{j} W_{j \to i} P_j - P_i \sum_{j} W_{i \to j} $$

Because there is no connection between $i$ and any other $j$ (where $j \neq 0$), the only state that talks to $i$ is the Ground State ($0$).
The sum collapses to a single term:

$$ 0 = W_{0 \to i} P_0 - W_{i \to 0} P_i $$

Rearranging this gives us the **Detailed Balance**-like relation (valid even in non-equilibrium):

$$ P_i = P_0 \frac{W_{0 \to i}}{W_{i \to 0}} $$

Let's call the ratio of rates $R_i = \frac{W_{0 \to i}}{W_{i \to 0}}$.
So, $P_i = P_0 \cdot R_i$.

#### Normalization
We know the total probability must be 1:
$$ \sum_{k=0}^{N} P_k = 1 $$
$$ P_0 + \sum_{i=1}^{N} P_i = 1 $$

Substitute $P_i = P_0 R_i$:
$$ P_0 + \sum_{i=1}^{N} (P_0 R_i) = 1 $$
$$ P_0 \left( 1 + \sum_{i=1}^{N} R_i \right) = 1 $$

#### Final Analytic Formula
1.  **Ground State Probability:**
    $$ P_0 = \frac{1}{1 + \sum_{i=1}^{N} \left( \frac{W_{0 \to i}}{W_{i \to 0}} \right)} $$

2.  **Excited State Probability:**
    $$ P_i = P_0 \cdot \frac{W_{0 \to i}}{W_{i \to 0}} $$

### 3. Current Calculation (The payoff)

The current $I$ through the tip is the sum of currents from all transitions.
$$ I = \sum_{\text{transitions } a \to b} (W_{a \to b}^{\text{Tip}} P_a - W_{b \to a}^{\text{Tip}} P_b) $$

For the Star Graph, we only have transitions $0 \leftrightarrow i$.
$$ I = \sum_{i=1}^{N} \left( W_{0 \to i}^{\text{Tip}} P_0 - W_{i \to 0}^{\text{Tip}} P_i \right) $$

Substitute $P_i$:
$$ I = \sum_{i=1}^{N} \left( W_{0 \to i}^{\text{Tip}} P_0 - W_{i \to 0}^{\text{Tip}} \left( P_0 \frac{W_{0 \to i}^{\text{Total}}}{W_{i \to 0}^{\text{Total}}} \right) \right) $$

$$ I = P_0 \sum_{i=1}^{N} \left( W_{0 \to i}^{\text{Tip}} - W_{i \to 0}^{\text{Tip}} \frac{W_{0 \to i}^{\text{Total}}}{W_{i \to 0}^{\text{Total}}} \right) $$

*(Note: $W^{Total} = W^{Tip} + W^{Substrate}$)*.

### 4. Implementation Algorithm (O(N))

This is how your kernel would look. No sorting required. No matrix inversion.

**Inputs:** `GS_config`, `E_site`, `W_interact`, `Tip_Params`...

**Kernel Logic (Per Pixel):**
1.  **Initialize:**
    *   `Sum_Rate_Ratios = 0.0`
    *   `Sum_Current_Terms = 0.0`
    *   Calculate `E_GS` (Energy of Ground State).

2.  **Loop over Sites ($i = 0 \dots 127$):**
    *   **Get Excitation:** Imagine flipping bit $i$.
    *   **Energy:** Calculate $\Delta E = E_{excited} - E_{GS}$.
        *   *Check:* Is $\Delta E$ too high? (Optional cutoff).
    *   **Rates:**
        *   Calculate $W_{0 \to i}$ (Rate to enter excited state).
            *   Depends on Fermi function $f(\Delta E - \mu_{tip})$.
        *   Calculate $W_{i \to 0}$ (Rate to return to GS).
            *   Depends on $1 - f(\Delta E - \mu_{tip})$.
    *   **Accumulate:**
        *   `Ratio = W_in / W_out`
        *   `Sum_Rate_Ratios += Ratio`
        *   `Current_Term = (W_in_Tip - W_out_Tip * Ratio)`
        *   `Sum_Current_Terms += Current_Term`

3.  **Finalize:**
    *   `P_GS = 1.0 / (1.0 + Sum_Rate_Ratios)`
    *   `Total_Current = P_GS * Sum_Current_Terms`

### 5. The Critical Condition

**You must check this:**
Does your Hamiltonian have a term $\sum t_{ij} c^\dagger_i c_j$ (hopping between quantum dots)?

*   **YES:** You have a "Wheel" graph. The Star solution is wrong. You must use Gauss-Jordan.
*   **NO:** You have a "Star" graph. This analytic method is exact and 100x faster.

*Self-Check:* You said "Hubbard/Ising model... these are not spins but charge states... hopping from tip".
Usually, if you call it "Ising-like", it implies purely diagonal interactions ($n_i n_j$) and no kinetic hopping term between sites ($c^\dagger_i c_j$). If that is true, **you have a Star Graph.**

---

# USER

OK, I understad, if we always must make transitions over the ground stata, then it is true. Now I need to think if this is the case in our physical system. Up to now I somulated only small system (3-4 sites) where I considered all many body state. 

1) It is definitely true that we must consider only single flip transitions in single event.
2) what is not so clear is if there is some non-equlibirum occupancy also the initial stats other than ground state can be populated, and then they can hop into other state (with hamming distance from ground state 2, ie >1), and these can transition back co singly excited state.

Paths like this

|...0...0...>   Ground state
|...1...0...>   singly excited state |i>
|...1...1...>   doubly excited state |i,j>
|...0...1...>   singly excited state |j>
|...0..01...>   back to ground state

all these are single electron transitions, but the path goes over multiply occupied states

the question is if we can necglect (from stat disallow) these paths. I don't know we should critically analyse this. 

But definitely it is worth to try make simple solver with this assumptions, it will save use a lot of trubles. Then we can see how it differs from full solver for some small systems.

analyse the problem from math/physics perspective, do not write code yet

---

# Gemini-3-pro

This is a critical physical insight. You have identified the distinction between the **Star Topology** (purely radial, independent excitations) and the **Diamond/Mesh Topology** (connected excitations allowing sequential transport loops).

Here is the rigorous analysis of that specific path:
$|00\rangle \to |10\rangle \to |11\rangle \to |01\rangle \to |00\rangle$.

### 1. Physics Analysis: Sequential vs. Virtual

Whether we can neglect the doubly excited state $|11\rangle$ depends entirely on **Energy** and **Temperature/Bias**.

#### A. The Coulomb Blockade Regime (Star is Valid)
If the onsite interaction energy $W_{ij}$ (or $U$) is large:
$$ E_{11} \gg E_{10}, E_{01} $$
Specifically, if $E_{11}$ is well above the chemical potential of the tip ($\mu_{tip}$):
1.  **Population:** The probability $P_{11}$ is exponentially suppressed ($e^{-(E_{11}-\mu)/kT} \approx 0$).
2.  **Transport:** The transition $|10\rangle \to |11\rangle$ requires an electron to enter from the tip. If $E_{11} - E_{10} > eV_{bias}$, this transition is **blocked**.
3.  **Result:** The path is cut. The loop is broken. The system effectively behaves as a Star Graph.

#### B. The Conducting/High-Bias Regime (Star Fails)
If the bias is large ($eV_{bias} > W_{ij}$) or interaction is attractive:
1.  **Population:** State $|11\rangle$ becomes energetically accessible.
2.  **Transport:** The system can realistically hop $|10\rangle \to |11\rangle$.
3.  **Result:** The "Star" approximation would explicitly forbid this step, forcing the electron at $|10\rangle$ to exit back to $|00\rangle$ before $|01\rangle$ can be populated. **This yields incorrect current.** It underestimates the total flow because it blocks a valid transport channel ($|10\rangle \to |11\rangle \to |01\rangle$).

### 2. The "Subspace" Solution (The Robust Approach)

You are correctly concerned that the "Star" topology is too restrictive for general cases. However, we also cannot afford the full exponential Hilbert space.

**The Solution is in your selection algorithm.**

If you select the **Top $N$ Lowest Energy States**, you do not need to force a topology. The physics will dictate the topology automatically.

**Scenario:**
Imagine a 2-site system.
*   GS: $|00\rangle$ (Energy 0)
*   Excitation A: $|10\rangle$ (Energy $\epsilon$)
*   Excitation B: $|01\rangle$ (Energy $\epsilon$)
*   Double: $|11\rangle$ (Energy $2\epsilon + U$)

**Algorithm:**
1.  **Sort:** You scan all states (GS, single flips, maybe double flips if you want, but let's stick to single flips from the current "basis").
    *   *Wait, here is a subtle trap in the "Single Flip from GS" definition.*
    *   If $|11\rangle$ is the global ground state, your MC finds it.
    *   If $|00\rangle$ is the global ground state, your "Single Flip" scanner checks $|10\rangle$ and $|01\rangle$. **It does NOT check $|11\rangle$ because that is 2 flips away.**

**CRITICAL FLAW REVEALED:**
If you strictly scan only "Hamming Distance = 1" from the MC Ground State, **you will never find $|11\rangle$** if the GS is $|00\rangle$.
*   This explicitly breaks the $|10\rangle \to |11\rangle$ link because $|11\rangle$ is simply not in your basis set.
*   **Consequence:** You are mathematically forcing the Star Topology. You are forbidding sequential tunneling through higher occupancy states.

### 3. Is this Acceptable? (The 1st Order Argument)

You stated: *"I consider only single electron hopping... 1st order perturbative approach."*

*   **Sequential Tunneling (Real):** $|00\rangle \xrightarrow{t} |10\rangle \xrightarrow{t} |11\rangle$.
    *   This is a sequence of two real, energy-conserving (within $kT$) 1st-order events.
    *   **Validity:** This happens only if $|11\rangle$ is populated.
*   **Cotunneling (Virtual):** $|10\rangle \xrightarrow{virtual |11\rangle} |01\rangle$.
    *   This is a 2nd order coherent process.
    *   **Validity:** This happens even if $|11\rangle$ is high energy (not populated).

**Verdict:**
By limiting your basis to $Hamming(GS, State) \le 1$, you are correctly modeling **1st Order Tunneling** near the Ground State.
You are neglecting **Sequential Tunneling through high-energy states** (which is correct, as they aren't populated) AND **Sequential Tunneling through low-energy states that are 2 flips away** (this is the risk).

**The Risk:**
If $|11\rangle$ is actually very low energy (e.g., nearly degenerate with $|10\rangle$ and $|01\rangle$), your simulation will miss it.
*   However, if $|11\rangle$ were that low in energy, your **Monte Carlo** would likely have visited it or found it as a candidate for the Ground State (or close to it).

### 4. Recommendation: The "Star + Loops" Solver

We can compromise. We can stick to the $Hamming=1$ basis (Star Topology) for the *selection*, but use the **General Solver** (Gauss-Jordan) for the *solution*.

**Why?**
Even within the $Hamming=1$ subspace, lateral hopping (if it existed) would create loops.
More importantly, if your MC finds $|10\rangle$ as the GS (e.g., 1 electron on the island), then:
*   $|00\rangle$ is a neighbor (dist=1).
*   $|11\rangle$ is a neighbor (dist=1).
*   Your basis *will* include $|11\rangle$.
*   You *will* capture the path $|00\rangle \leftrightarrow |10\rangle \leftrightarrow |11\rangle$.

**Conclusion on Math/Physics:**
1.  **Basis:** Limiting to $Hamming \le 1$ from GS is consistent with "1st order tunneling dominating transport". If a state is 2 flips away, reaching it requires passing through an intermediate. If that intermediate is the GS, you capture it. If the intermediate is an excited state, reaching the double-flip requires *two* excitation steps, which is exponentially less likely in the low-temperature/perturbative limit.
2.  **Solver:** You **should use the Gauss-Jordan (General) Solver**, not the Analytic Star Solver.
    *   *Reason:* While the basis selection looks like a Star, the physical connectivity might not be. If you have any term that connects two excitations directly (like lateral hopping), or if you define "ground state" locally such that transitions exist between selected neighbors, the Star solver will give $0$ current where there should be flow.
    *   Also, the Gauss-Jordan solver is generic. It doesn't care if the graph is a Star or a Mesh. It is safe. The Analytic solver is risky if your assumptions break.

### Summary Plan

1.  **MC:** Finds Global Minimum (State A).
2.  **Basis Gen:** Finds all states $B$ where $Hamming(A, B) = 1$.
3.  **Filter:** Keeps top 64 states (State A + 63 Bs).
4.  **Matrix Build:**
    *   Calculate rates for $A \leftrightarrow B_i$ (Always exists).
    *   Calculate rates for $B_i \leftrightarrow B_j$ (**Check this!**).
        *   If your Hamiltonian has no lateral hopping ($t_{ij}=0$), these are ALL ZERO.
        *   If so, the matrix is mathematically a Star.
5.  **Solve:** Use **Gauss-Jordan**.
    *   It costs the same $O(N^3)$ for $N=64$ (negligible).
    *   It protects you against future changes (e.g., adding lateral hopping or changing the basis selection logic).

This approach is physically sound for the regime you described and computationally safe.

---

# USER

I think you are doing serious mistake here when you conisder energy and boltzman probability (theram equalibirum distribution) as reliable criterium for selection. In my view it is only rought estimate. But in reality we cannot rely on that because we are out of thermal equalibirum (otherwise we do not need Pauli Master Equation and we are fine with the ground-state solver). In our system we have very low temperature (order of few kelvin or even milikelvin) but we have substantial current, which can creat substantially non-equalibirum distribution of states (especially due to barriers caused by forbiden tuneling between the sites). 

So the doubly excited state |i,j> can be high in energy but that is not saying it will not be populated. 

What is importaint is the conditional probability. And that is what we may use to find relevant states.

We have initial state (grund state) and compute probability to hop to any other singly exited state (this will depend on position of tip which inject electron/holes to the system, states far from tip will have exponentially lower probability to get excited). Now once we have these singly exicted states, we can create doubly excited states from these singly excited states, this again depend on the hopping matrix element to the tip and fermi functions ( we have fermi function of tip and of substrate, and the allowed transiton must fulfil fermi golden rule, conservation of energy we can only hop if change of quantum dot energy equals to chemical poential of adding/removing electron to/from tip/substrate), and all is modulated by t_i=<tip|r_i>~exp( beta|r_i-r_tip| ), where r_i is position of site i. So changes of sites near tip are much more probable than changes of sites far from tip even if the doubly occupied states have higher energy.

We assume that at low current the distribution woul,d be still not so far from equalbitum (we can neglect triply or quadruply occupied states), but distance/hopping rate may be more importaint criterium for non-equalibirum occupancy than energy alone.

---

# Gemini-3-pro

This is a profound insight that changes the architecture of the "Scanner" completely. You are absolutely correct: **In non-equilibrium systems, Kinetics dominates Thermodynamics.**

If the tunneling rate $\Gamma$ varies by orders of magnitude (due to the exponential dependence on the Tip-Site distance $e^{-\kappa r}$), a high-energy state near the tip will be populated orders of magnitude faster than a ground-state-like transition far from the tip.

Using Boltzmann probability ($e^{-E/kT}$) as a sorting criterion is indeed wrong for transport because it assumes the system has time to equilibrate with a bath, whereas here the system is being driven by the tip faster than it can relax.

### The New Strategy: "Kinetic Beam Search"

Instead of searching for "Low Energy," we must search for **"High Flow."**

We need to construct the subspace by estimating the **Stationary Probability** of reachable states.

#### 1. The Physics of the Selection Metric

Let's define the metric for selecting the single-flip state $|i\rangle$ (where site $i$ is flipped relative to Ground State $|0\rangle$).

We approximate the system locally as a two-level system ($|0\rangle \leftrightarrow |i\rangle$) to get a first-order estimate of its population $P_i^{est}$.

$$ P_i^{est} \approx \frac{\Gamma_{0 \to i}}{\Gamma_{0 \to i} + \Gamma_{i \to 0}} $$

Where:
*   $\Gamma_{0 \to i}$: Total rate to enter state $|i\rangle$ (sum of Tip and Substrate contributions).
*   $\Gamma_{i \to 0}$: Total rate to leave state $|i\rangle$ back to $|0\rangle$.

**Why this captures your intuition:**
1.  **Tip Proximity:** If site $i$ is near the tip, the tunneling matrix element $T_{tip}$ is huge. This boosts $\Gamma_{0 \to i}$ (if energy allows).
2.  **Energy Blockade:** The rates $\Gamma$ include the Fermi functions $f(E)$ and $1-f(E)$. If the transition is energetically forbidden (Coulomb blockade), the Fermi factor is zero, $\Gamma \to 0$, and the state is correctly discarded.
3.  **Trapping:** If a state is hard to enter but *very* hard to leave (a "dark state" or trap), $\Gamma_{i \to 0}$ is small, making $P_i^{est}$ large.

**The Criterion:**
We calculate $P_i^{est}$ for all single-flip candidates and select the Top $N$.

#### 2. Handling Double Excitations (The "Active Pair" Heuristic)

You raised the importance of paths like $|0\rangle \to |i\rangle \to |i,j\rangle \to |j\rangle \to |0\rangle$.
Scanning all $N(N-1)/2$ double excitations is too slow ($128^2 \approx 16,000$).

However, we can use the result of the Single-Flip scan to intelligently guess the relevant double excitations.

**The Logic:**
If state $|i\rangle$ has high population, and state $|j\rangle$ has high population, then the double state $|i,j\rangle$ is the most likely 2nd-order state to be populated. It represents the system trying to load electrons onto *both* "active" sites simultaneously.

**Proposed Algorithm:**
1.  **Scan Single Flips:** Calculate $P_k^{est}$ for all 128 sites.
2.  **Select Top Single:** Pick the top $K_{single}$ (e.g., 48) sites.
3.  **Construct Double Candidates:**
    *   Take the very top $K_{top}$ (e.g., top 4) sites from that list.
    *   Form all pairs among them ($4 \times 3 / 2 = 6$ pairs).
    *   Add these $|i,j\rangle$ states to the subspace.
    *   *Physics Check:* This captures the "sequential loading" of the most active sites near the tip.
4.  **Total Subspace:** $\{GS\} \cup \{ \text{Top Singles} \} \cup \{ \text{Top Active Doubles} \}$.

#### 3. The Math for the Kernel (Transition Rates)

To implement this efficient "Kinetic Scanner," we need the explicit rate formulas in the kernel.

**Definitions:**
*   $\Delta E = E_{excited} - E_{GS}$.
*   $T_{tip}$: Tunneling factor (computed in Kernel 1, depends on $r_{tip}$).
*   $T_{sub}$: Substrate coupling (constant parameter).
*   $\mu_{tip} = V_{bias}$, $\mu_{sub} = 0$.

**Rates (Fermi's Golden Rule):**
*   **In-rate ($\Gamma_{in}$):** Rate to *create* the excitation (e.g., add electron).
    $$ \Gamma_{in} = 2\pi [ |T_{sub}|^2 f(\Delta E - \mu_{sub}) + |T_{tip}|^2 f(\Delta E - \mu_{tip}) ] $$
*   **Out-rate ($\Gamma_{out}$):** Rate to *destroy* the excitation (e.g., remove electron).
    $$ \Gamma_{out} = 2\pi [ |T_{sub}|^2 (1 - f(\Delta E - \mu_{sub})) + |T_{tip}|^2 (1 - f(\Delta E - \mu_{tip})) ] $$

*   *Note: If the excitation involves removing an electron (hole), flip the $f$ and $1-f$ logic.*

#### 4. The Revised Solver Architecture

We now have a complete picture. We need a 3-stage process within the Workgroup.

**Kernel: `solve_pme_kinetic`**

**Stage 1: Kinetic Scan (Parallel)**
*   Each thread picks a set of sites (strided loop).
*   Compute $\Delta E$ (Energy cost).
*   Compute $\Gamma_{in}, \Gamma_{out}$.
*   Compute Score $S = \Gamma_{in} / (\Gamma_{in} + \Gamma_{out} + \epsilon)$.
*   Store `(Score, BitIndex)` in LDS.

**Stage 2: Sort & Expand (Parallel)**
*   **Sort:** Parallel bitonic sort of the 128 scores in LDS.
*   **Select:** Keep top $N_{basis}$ indices.
*   **Augment (Doubles):**
    *   Take indices of Top 4.
    *   Compute bitmasks for pairs (e.g., `mask[0] | mask[1]`).
    *   Append these to the basis list.
    *   *Total Basis Size:* $\approx 50-60$. Fits in 64x64 matrix.

**Stage 3: Matrix Solve (Parallel)**
*   Construct the Rate Matrix $K$ for this basis.
    *   *Crucial:* Now that we have Doubles $|i,j\rangle$, we must calculate rates between $|i\rangle \leftrightarrow |i,j\rangle$ and $|j\rangle \leftrightarrow |i,j\rangle$.
    *   This naturally builds the "Diamond/Mesh" topology you asked for.
*   Run Gauss-Jordan Elimination.
*   Compute Current.

### Why this addresses your concern

This approach explicitly solves the **Non-Equilibrium Selection** problem.

1.  **Tip Dominance:** A state with $\Delta E = 1 eV$ (huge) but effectively infinite tunneling barrier will have Score $\approx 0$. A state with $\Delta E = 1 eV$ but right under the tip will have a non-zero Score because the driving force is high.
2.  **Barrier Detection:** If a transition is blocked by the Fermi window (e.g., bias is too small to overcome $\Delta E$), the Fermi function is 0, Rate is 0, Score is 0. The state is ignored.
3.  **Active Cycles:** By including the pairs of the most active single sites, we capture the dominant cyclic current paths $|0\rangle \to |i\rangle \to |i,j\rangle \to |j\rangle \to |0\rangle$ without calculating the full exponentially large space.

This is the most robust approximate solver for large-scale Pauli Master Equations on a GPU. It combines Monte Carlo (for the global vacuum state) with a Kinetics-aware subspace expansion.

---

# USER

I do not like your Kernel 2: Kinetic Subspace Solver (The "Big Gun")
it is too monolithic. We should definitely split the search for relevant sub-space and the solver (linear algebra). This will be both more debugeable - we can check which states we includied what is theyr estimated probability, where they are located (1st pass - sach of substapce) and then separately do the solver (linear algebra) it will also reduce register pressure, and allow us better optimize dimensiton and partitioning. It should be monudlar so we can echvinage different version of the matrix solver and subspace-search independently.
maybe for the search of subspace - we can do two passes,
in 1st pass we will find the singly occupied stats and sort them by estimated probability. (even before sorting, we should filter them by some probability cutoff - if the transition probability is too low, we do not include them in the list for sorting)
in second pass we take top n generated in the sorted list from the first pass and generate new (double exxited states). We must make sure we do not duplicate any states I mean |0>->|i>->|i,j> and  |0>->|j>->|i,j> produce same state |i,j> by two different paths
Perhaps 1 and 2 can be done by the same universal kernel which always takes list of intiial state (can be just one ground state, but also more than one intital states), and generate all relevant exitations, filtered by some trashold probability, and do deduplication
Another possibility is to write single kernel which systematically try to generate all double excitations (i,j), by product of probability pij~pi*pj, that may work as well. Think it through.
please write new modular 2 or 3 pass version
Even before you start writing the kernel write clear physical equation reich documentation bout how to compute transition probability from fermi golder rule considering energy conservation and fermi-dirac distribution fuction of electron in substrate and the tip (tip is shifted wrt substrate by VBias), matrix hoping element from tip depend exponentially on distance, hopping from substrate is constant. When I was solvinf this system importaint insight is that typically voltage is so high that current flow only in one direction. E.g. if voltage shift tip distribution function to higher energies over substrate, then electron always hop from tip to QD and then to sample never opposite. At inverted bias voltage the flow is possoite but again unidirectional. We almost never see hopping electron from tip to QD and back to tip.
Write all this in physically rigroous way with all equations. You can again check out current implementation in OpenCL

---

# Gemini-3-pro

Here is the comprehensive response, starting with the Rigorous Physics Documentation, followed by the Modular OpenCL implementation.

### Part 1: Physics & Equations (Rigorous Documentation)

**1. System Hamiltonian**
The system consists of $N$ sites (quantum dots) occupied by electrons. The many-body state is defined by the occupancy vector $\vec{n} = \{n_1, n_2, \dots, n_N\}$, where $n_i \in \{0, 1\}$.
$$ H_{sys}(\vec{n}) = \sum_{i} \epsilon_i n_i + \sum_{i < j} W_{ij} n_i n_j $$
*   $\epsilon_i$: On-site energy (includes intrinsic level + electrostatic potential from the STM Tip).
*   $W_{ij}$: Coulomb interaction energy between electrons at sites $i$ and $j$.

**2. Tunneling & Rates (Fermi's Golden Rule)**
Transitions occur via single-electron tunneling events between a site $i$ and one of the leads (Tip or Substrate).
*   **Initial State:** $|a\rangle$, Energy $E_a$, Occupancy $\vec{n}_a$.
*   **Final State:** $|b\rangle$, Energy $E_b$, Occupancy $\vec{n}_b$.
*   **Transition Energy:** $\Delta E = E_b - E_a$.

**Tunneling Matrix Elements:**
*   **Substrate:** Assumed constant coupling $\Gamma_{sub}$.
*   **Tip:** Depends exponentially on tip position $\vec{r}_{tip}$.
    $$ \Gamma_{tip, i}(\vec{r}_{tip}) = C_{pre} \cdot e^{-2 \beta |\vec{r}_{tip} - \vec{r}_i|} $$
    *(Note: In code, we often pass `T_factor` $= e^{-\beta r}$, so rate $\propto T_{factor}^2$).*

**Rate Equations:**
The rate $W_{a \to b}$ depends on whether an electron is **added** ($N_b = N_a + 1$) or **removed** ($N_b = N_a - 1$).

*   **Case 1: Adding an Electron (Tunneling IN)**
    The lead must provide an electron with energy $\Delta E$. The probability is proportional to the occupancy of the lead at that energy (Fermi function).
    $$ W_{a \to b}^{L} = \frac{1}{\hbar} \Gamma_{L,i} \cdot f(\Delta E - \mu_L) $$

*   **Case 2: Removing an Electron (Tunneling OUT)**
    The lead must accept an electron with energy $-\Delta E$ (or electron leaves with energy $\Delta E$). The probability is proportional to the *vacancy* of the lead ($1-f$).
    $$ W_{a \to b}^{L} = \frac{1}{\hbar} \Gamma_{L,i} \cdot [1 - f(\Delta E - \mu_L)] $$

*   **Fermi Function:** $f(E) = \frac{1}{1 + e^{E/k_B T}}$.
*   **Chemical Potentials:** $\mu_{sub} = 0$, $\mu_{tip} = e V_{bias}$.

**3. Kinetic Selection Criterion (Non-Equilibrium)**
In a high-bias regime ($eV_{bias} \gg k_B T$), current is unidirectional.
*   **Example ($V_{bias} > 0$):**
    1.  Tip has high $\mu$. Electron hops Tip $\to$ Dot ($i$).
        *   Requires: $\Delta E < \mu_{tip}$ (energetically allowed) AND site $i$ is near Tip.
    2.  Substrate has low $\mu$. Electron hops Dot ($i$) $\to$ Substrate.
        *   Requires: $\Delta E > \mu_{sub}$.
*   **Selection Metric:** We estimate the population $P_i$ of a single excitation in a two-level approximation:
    $$ P_i^{est} \approx \frac{\Gamma_{in}}{\Gamma_{in} + \Gamma_{out}} $$
    *   If $P_i^{est}$ is very low, the state is effectively unreachable (blocked).
    *   If $P_i^{est}$ is high, the state is a "waiting room" for transport.
    *   **Double Excitations:** A double state $|i,j\rangle$ is formed by a transition from $|i\rangle$ (adding $j$) or from $|j\rangle$ (adding $i$). Its relevance is proportional to the product of its parents' activity.

### Part 2: Modular OpenCL Implementation

Here are the modular kernels.
1.  **`kinetic_basis_scanner`**: Finds relevant subspace (Single & Double) based on physics, writes bitmasks to global memory.
2.  **`solve_pme_dense_batch`**: Reads bitmasks, builds matrix, solves PME, calculates current.

#### Kernel 1: The Scanner (Subspace Generator)

#### Kernel 2: The Solver (Linear Algebra)

### Analysis of this approach

1.  **Debuggability:** You can verify the `basis_out` array in Python. If the solver fails, you check: "Did the scanner pick the right states?" "Are the doubles correctly formed?"
2.  **Performance:**
    *   **Scanner:** Uses negligible registers. Just iterating energies and sorting floats.
    *   **Solver:** Matrix fits in LDS ($64 \times 64 = 16\text{KB}$). Threads focus purely on FMA operations.
3.  **Physics Compliance:**
    *   The `calc_total_rate` function correctly incorporates the Fermi blockade logic.
    *   The "Kinetic Score" in the scanner ensures that high-energy states are only included if they are very close to the tip (high $T_{factor}$ overcoming the Energy penalty).

---

# GPT-5.2 Low

## What’s implemented now (end-to-end)

### **Pipeline**
- **MC ground-state**: unchanged, still [solve_mc_2phase()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:721:4-820:19) producing `occupation` (bitmask bytes per tip/pixel).
- **PME postprocess A (fast, single-flip “star”)**:
  - Kernel: `solve_pme_star_analytic`
  - Host: [HubbardSolver.solve_pme_star(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:303:4-311:19)
- **PME postprocess B (kinetic non-eq subspace + dense solver)**:
  - Kernel 1 (basis): `kinetic_basis_scanner` (GS + top kinetic singles + doubles from top parents)
  - Kernel 2 (solver): `solve_pme_dense_batch` (build rate matrix, Gauss-Jordan, compute tip current)
  - Host: [HubbardSolver.solve_pme_dense(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:313:4-328:19)

### **Test/demo**
- [run_Hubbard.py::demo_local_update()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/run_Hubbard.py:535:0-806:24) now does:
  1) precalc `Esite,Tsite`
  2) MC ([solve_mc_2phase](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:721:4-820:19))
  3) PME currents:
     - [pme_curr_star = solver.solve_pme_star(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:303:4-311:19)
     - [pme_curr_dense = solver.solve_pme_dense(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:313:4-328:19)
  4) Saves figure: `./figs/pme_currents_V0.50.png`
- I also changed `__main__` to default `show_plots=False` to avoid hanging in headless / background runs.

### **Ran successfully**
I ran:
```bash
python -u -m pyBall.OCL.run_Hubbard
```
Result:
- Exit code `0`
- Produced:
  - `./figs/pme_currents_V0.50.png`
  - plus existing outputs (`site_maps_imshow_...`, `monte_carlo_results_...`)

## Kernels added/modified in [pyBall/OCL/cl/hubbard.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/hubbard.cl:0:0-0:0)

### **1) `solve_pme_star_analytic`**
- **Purpose**: analytic steady-state for **GS + all single flips**, assuming star topology (no direct excitation↔excitation transitions).
- **Uses**:
  - `deltaE_flip_u4()` (ΔE for flip from sparse W + Esite)
  - `calc_total_rate()` (tip+sub)
  - `rate_tip_only()` (tip-only contribution)

### **2) `kinetic_basis_scanner`**
- **Purpose**: modular *subspace selection* under **non-equilibrium kinetic criterion**.
- **Output**:
  - `basis_out[nTips, MAX_BASIS, 4]` stored as 4×`uint32` words (packed as `int` in global buffer)
  - `n_basis_out[nTips]`
- **Selection metric**:
  - For each site flip `i`, compute ΔE and rates `r_in`, `r_out`
  - score = `r_in/(r_in+r_out)` (0 if blocked)
- **Selection method**:
  - Serial top-K selection on `lid==0` (robust; not fastest, but works and debuggable)
  - Adds doubles from top `N_PARENTS_DOUBLE` parents (`i<j` ⇒ no duplicates)

### **3) `solve_pme_dense_batch`**
- **Purpose**: modular *linear algebra stage*:
  - loads basis masks
  - computes state energies
  - builds generator matrix (only Hamming distance 1 transitions allowed)
  - enforces conservation on diagonal
  - replaces row0 by normalization `sum(P)=1`
  - solves by Gauss-Jordan (no pivoting currently)
  - computes a tip-current-like quantity by summing tip-only rates weighted by probabilities

## pyOpenCL harness in [pyBall/OCL/HubbardSolver.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:0:0-0:0)

Added:
- **Buffers** via [realloc_pme_buffers()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:100:4-111:36):
  - `occ_gs`, `tipTf`, `pme_curr`, `basis`, `n_basis`
- **Kernel setup wrappers**:
  - [setup_solve_pme_star_analytic](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:217:4-238:27)
  - [setup_kinetic_basis_scanner](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:240:4-260:27)
  - [setup_solve_pme_dense_batch](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:262:4-282:27)
- **High-level APIs**:
  - [solve_pme_star(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:303:4-311:19)
  - [solve_pme_dense(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:313:4-328:19)
- **Input prep**: [pme_prepare_inputs(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:284:4-301:43)
  - currently reuses `Tsite` as `tip_T_factors` (this matches your earlier approach; later we can pass a true `exp(-beta r)` buffer if you prefer)

## What’s still missing / needs follow-up (important)

### **1) `solve_pme_dense_batch` correctness / physics conventions**
Right now it compiles and runs, but there are physics details that likely need refinement:
- **Add vs remove logic** in dense solver:
  - I used “adding = bit in final state”.
  - This is a simplification; strictly you want to decide whether the transition is electron-in or electron-out relative to leads, and use `f` vs `(1-f)` per lead consistently.
- **Current expression**:
  - currently sums `P_c * rate_tip_only(c→r)` over directed edges.
  - You probably want **net tip current**: `Σ_{c,r} ( P_c W_tip(c→r) - P_r W_tip(r→c) )` with sign convention.
- **No pivoting** in Gauss-Jordan:
  - can fail / produce `NAN` on near-singular matrices. For robustness you likely want partial pivoting (still doable in LDS).

### **2) Duplicate helper block cleanup**
I disabled the second “prototype helper” section with `#if 0` to stop OpenCL compile errors (duplicate `fermi`, `PI`, etc.). If you want, I can properly remove or refactor it later (but per your rules I avoided larger restructuring beyond what was necessary to compile).

### **3) Performance**
- `kinetic_basis_scanner` currently does top-K selection serially on `lid==0`. Works, but isn’t optimal.
- Once correctness is confirmed, we can switch to:
  - bitonic sort in local memory, or
  - parallel reduction-based top-K (better register pressure control)

## Status summary
- **Implemented**: kernels, Python harness, demo plotting, and **successfully ran** the test producing current maps.
- **Next**: validate/adjust dense PME rate formulas + net-current definition + add pivoting, then optimize scanner sorting.

---

# USER

OK, I see that the output in 
/home/prokop/git/FireCore/figs/pme_currents_V0.50.png
is a bit strange. We should make more rigorous tests. We should compare partity with full solver for small system (2-sites, 4-sites)

look here @pauli_ocl.py @PME.cl this is full solver for small systems, we should reproduce it
1) first we must check if our tip potential and tip hopping is the same (perhaps in @PME.cl is more advanced model of tip potential which we should reimplement in @hubbard.cl )
Analyse the details
2) after that we should check if the site energies are the same for some 1D scan of tip positions (or even one tip position so we can use debung prnt inside the kernels, we need to achive partity there)
3) then we evaluate many body state energies, making sure they are the same btween the monte carlo solver and the full solver
4) then we should check if the groun sate is the same (again for 1D scan)
5) then we compute the excitations, and include all of them in the active sub-space, and solve them, we compare if the probabilities are the same.
6) Only then it make sense to compare the current.

analyse the differences between current codes in @hubbard.cl @PME.cl and discribe them clearly. 

make detailed notes. 

Then make comprehansive plan how to do the testing in multiple steps


---

# GPT-5.1-Codex

Below is a concise comparison and a proposed parity-test plan.

Key model differences: hubbard.cl vs PME.cl vs pauli_ocl.py
1) Tip electrostatics / on-site energy (Esite)
- PME.cl (kernel build_tip_arrays / helpers):
  - Uses multipole expansion with per-site rotation: d, d_mir rotated by a 3×3 R; eval_multipole(d_rot, order, multipole_cs); optional mirror term; ramp term applied to monopole.
  - Parameters: Rtip, zV0, zVd, beta, bMirror, bRamp; uses v_tips per pixel; rotations rots array.
- hubbard.cl (precalc_Esite_Thop):
  - Also uses Emultipole(cs, pTip, nMulti) with per-site coefficients; mirror toggle present but commented out; rotation of pTip/pTipMirror is currently commented out; ramp applied (monopole only) if enabled.
  - Rtip and zV (zMirror + offset) passed; uses tip_data.w as Vbias.
- Likely mismatch: PME.cl rotates into site frame; hubbard.cl currently skips rotation. Mirror handling differs (PME uses bMirror flag; hubbard currently comments out mirror subtraction). Ramp uses zV1=z_tip+zVd in PME; hubbard uses zV array similarly but confirm zVd vs zRampOffset naming. Need to confirm eval_multipole vs Emultipole parity (order vs nMulti, and cs layout).

2) Tip hopping (Thop / T_factors)
- PME.cl: T_factor = exp(-beta * |tip-site|).
- hubbard.cl: Thop = Thop_amp * exp(-Thop_decay * r).
- Potential mismatch: naming (beta vs Thop_decay), amplitude (Thop_amp vs implicit 1.0), and whether squared in rates (PME full solver uses tip factor in transition amplitude? verify pauli_ocl.py usage).
- Also PME.cl may use per-pixel v_tips and params[4]=beta; hubbard uses params dict.

3) Transition rates / Fermi factors
- PME full solver (pauli_ocl.py + PME.cl): likely uses add/remove logic consistent with lead chemical potentials and uses f/deltaE per lead separately.
- hubbard.cl PME kernels:
  - calc_total_rate uses G_sub + G_tip_pre * t_factor^2 with fermi(dE, mu, T) and (1-f) via is_add; rate_tip_only similar.
  - Star solver uses deltaE_flip_u4 with Esite_tip, W, etc.
- Possible sign/convention differences: add vs remove, t_factor squared, and whether current uses net tip current or only forward edges.

4) Many-body state handling
- PME full solver enumerates all states for small systems (pauli_ocl.py), solving full PME matrix.
- hubbard.cl uses Monte Carlo ground state, builds reduced basis (kinetic_basis_scanner) and dense solver limited to MAX_BASIS=64; star analytic assumes only GS-connected transitions.

Comprehensive parity test plan (small systems: 2-site, 4-site)
0) Align parameters and kernels
   - Match Rtip, zMirror/zV0, zRampOffset/zVd, beta vs Thop_decay, Thop_amp (set to 1.0), multipole order/coeffs, bMirror/bRamp, rotations.
   - Ensure tip potentials/hoppings computed identically: enable rotations and mirror in hubbard.cl to mirror PME.cl; set ramp flags same.

1) Single-tip, single-position Esite/Thop parity
   - Use one tip position; run PME.cl build_tip_arrays outputs (E shifts, T_factors) vs hubbard.cl precalc_Esite_Thop (Esite, Tsite).
   - Add debug prints in kernels (if needed) for one tip/site to check d, d_mir, rotated vectors, multipole value, final E, and hopping.
   - Compare arrays numerically in Python harness.

2) 1D tip scan parity (Esite/Thop)
   - Scan along x (or z) for fixed system; compare curves of Esite and Thop per site between PME.cl and hubbard.cl outputs.

3) Many-body energies parity
   - For small N (2, 4 sites), enumerate all occupations.
   - Compute many-body energies using hubbard.cl path (MC ground-state energy function) vs pauli_ocl.py full solver (or PME.cl host path).
   - Confirm GS energy and ordering match.

4) Ground-state parity
   - For scanned tips, confirm GS occupancy matches between MC ground-state search (or brute force in Python) and full solver enumeration.

5) Excitations and basis
   - Enumerate all single/double excitations; compare energies and deltaEs between the two implementations.
   - Temporarily set hubbard dense solver basis to include all states (since small N) to avoid truncation differences.

6) Transition rates parity
   - For each connected pair, compare lead-specific rates (tip/substrate) using the same dE, t_factor, mu, T, Gamma.
   - Align add/remove convention; check whether t_factor is squared in both places.

7) PME solution parity (probabilities)
   - Solve full PME in pauli_ocl.py/PME.cl and hubbard dense solver with full basis; compare stationary probabilities.

8) Current parity
   - Compute net tip current with matched formula: Σ(P_a W_tip(a→b) - P_b W_tip(b→a)).
   - Ensure sign conventions match.

Immediate actionable checks
- Confirm and, if needed, re-enable rotation and mirror handling in hubbard.cl precalc_Esite_Thop to match PME.cl.
- Match hopping decay: set Thop_amp=1 and Thop_decay=beta used in PME.cl for tests.
- Use a Python harness to call both: hubbard solver precalc and pauli_ocl full solver for the same small system and tip position(s); dump Esite/Thop arrays side-by-side.

If you want, I can draft the exact test script comparing steps (1) and (2) and outline the kernel changes (rotate + mirror) to align models.

---

# GPT-5.2

## Implemented tests (2-site embedded in 4-site) + debug gates

### What I added/changed
- **New test scripts**
  - [pyBall/OCL/test_pme_parity_single.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_single.py:0:0-0:0)
    - Single tip position + single bias.
    - Compares **Esite** and **Tsite** between:
      - full solver: [PauliSolverCL](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/pauli_ocl.py:28:0-289:12) ([pyBall/OCL/pauli_ocl.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/pauli_ocl.py:0:0-0:0) + [pyBall/OCL/cl/PME.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/PME.cl:0:0-0:0))
      - hubbard precalc: [HubbardSolver.precalc_esite_thop](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:512:4-570:79) ([pyBall/OCL/cl/hubbard.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/hubbard.cl:0:0-0:0))
    - Also compares **many-body energies** (derived from Esite) vs full solver `StateEs`.
    - **Enables GPU debug prints** via build option `-DDBG_PRECALC=1`.
  - [pyBall/OCL/test_pme_parity_sweep.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_sweep.py:0:0-0:0)
    - 1D **x-sweep** at fixed bias.
    - Plots PME vs Hubbard **Esite(x)** and **T(x)** and saves:
      - [/home/prokop/git/FireCore/figs/parity_tip_xsweep.png](cci:7://file:///home/prokop/git/FireCore/figs/parity_tip_xsweep.png:0:0-0:0)

- **Debug-print gates + build options**
  - [HubbardSolver.__init__](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:63:4-79:83) now accepts `build_options` and forwards it to OpenCL build (so we can do `-DDBG_PRECALC=1`).
  - [hubbard.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/hubbard.cl:0:0-0:0): added `DBG_PRECALC` macro gate (default `0`) to avoid unconditional `printf`.
  - Fixed [pauli_ocl.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/pauli_ocl.py:0:0-0:0) kernel path bug: it now correctly loads [pyBall/OCL/cl/PME.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/PME.cl:0:0-0:0) (it previously looked for `pyBall/cl/PME.cl`).

- **Critical parity fix in [hubbard.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/hubbard.cl:0:0-0:0)**
  To match [PME.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/PME.cl:0:0-0:0) I updated `precalc_Esite_Thop`:
  - **Enabled rotation** (transform `pTip`/`pTipMirror` into site frame using `rotSite`) — previously commented out.
  - **Enabled mirror subtraction** `E -= Emultipole(...pTipMirror...)` when `bMirror != 0`.
  - Fixed `Emultipole()` coefficient-count logic (dipole/quadrupole thresholds).
  - Zero-initialized local `cs[16]` so unused coeffs don’t contain garbage.

## Test results I ran (on your machine)

### 1) Single tip position parity
Command:
```bash
python -u pyBall/OCL/test_pme_parity_single.py
```

Result:
- **Esite parity:** `max|dE| = 0.0`
- **T parity:** `max|dT| = 0.0`
- **Many-body energies derived vs full solver StateEs:** small float noise (`~9e-5 eV` max), OK.

This means step **(1)** and a big chunk of **(2)** in your plan (site shifts + many-body energies) are now matched.

### 2) 1D x-sweep parity
Command:
```bash
python -u pyBall/OCL/test_pme_parity_sweep.py
```

Result:
- Figure saved: [/home/prokop/git/FireCore/figs/parity_tip_xsweep.png](cci:7://file:///home/prokop/git/FireCore/figs/parity_tip_xsweep.png:0:0-0:0)
- **Esite sweep parity:** `max|dE| = 0.0`
- **T sweep parity:** `max|dT| ~ 4.47e-08` (float32-level)

## What’s still missing (next steps for your list)
Your requested steps **3-5** (GS, probabilities, current) need **PME solve parity**, not just precalc parity.

- **Immediate next task** (already added to TODO):
  - Implement a test that solves PME using the Hubbard PME kernels (likely `solve_pme_dense_batch` with a full basis for the small system) and compare:
    - **probabilities** vs [PauliSolverCL](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/pauli_ocl.py:28:0-289:12) `out_probs`
    - **current** vs [PauliSolverCL](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/pauli_ocl.py:28:0-289:12) `out_current`
  - Then do the same over a **bias sweep** (or x sweep).

If you want strict 2-site (N=2) rather than “2-site embedded in 4-site”, we’ll need to either:
- generalize [PME.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/PME.cl:0:0-0:0) away from `#define N_SITES 4`, `N_STATES 16`, or
- keep the embedding approach (which is fine for parity as long as spectator sites are truly inert).

## Status
- **Done:** scripts for single-point + 1D scan; debug print gating; Esite/Thop parity fixed and verified.
- **Next:** implement PME(probabilities/current) parity tests using Hubbard PME solver vs full solver, then extend to bias sweep, then move to 4-site.

---

# USER

WOW, this is serious work you did at once, that is great. GOOD JOB! I want you now do all the requested parity checks, put all the results in single directory both as text file (e.g. csv/dat? or something human readable and as .png plots), make sure the results ar not matching only trivially (e.g. all zero) for that reason we should alaways print spread (min,max) of the variables along the scan. We should do 1D scans along x and along VBias and match the whole scans. 

MAke usre we use exactly the same parameters, and report the exactly parameters we are using, for sweeps report in which range of x_tip or VBias we do the scan.

The properties we want to output/plot and check parity are:
1) site energies (after gating by the tip potential VBias dependent) and hoppings (exponential form tip)
2) ground state energy and energy of all many body states
3) probabilites of many body states after solving PME
4) total current and currents passing through individual sites.

First do this all for 2-sites (can be embeded in 4), ouputing it in one directroy

Then do it for 4-sites (output it in other directory)

clearly refere what are the reslts and directrories where I can find the results.

Now make deteiled plan of this

---

# GPT-5.1-Codex

Here’s a concrete test plan to finish all parity checks and artifacts. I’m in Ask mode; to implement and run, please switch to Code mode.

Plan (2-site first, then 4-site)
1) Shared fixtures (both solvers)
   - Systems: 2-site (embedded in 4 as before with two far-away spectators), then true 4-site.
   - Common params: Rtip=3.0, zV0=-0.5, zVd=0.0, beta=0.5, Gamma=1.0, W_scalar=0 (unless you want a finite W), bMirror=1, bRamp=1, multipole cs[0]=1.0, order=0, rotations=identity.
   - Leads: mu0=0, T0=4 K (kB*4), mu1=Vbias, T1=4 K.
   - Tip Z: zTip=3.0.

2) Scans to run
   - X-scan: tip_x ∈ [-2, 7], y=0, z=3, Vbias fixed (e.g., 0.50 V).
   - V-scan: Vbias ∈ [0.0, 1.0], tip at x=2.5, y=0, z=3.
   - Both for 2-site and 4-site.

3) Outputs per scan (placed in dedicated directories)
   - Directory layout:
     - figs/parity_2site/xscan/
     - figs/parity_2site/Vscan/
     - figs/parity_4site/xscan/
     - figs/parity_4site/Vscan/
   - Files:
     - CSV (or TSV) with columns:
       - x or Vbias
       - Esite_full[i], Esite_hub[i] for all sites
       - T_full[i], T_hub[i]
       - GS_energy_full, GS_energy_hub (if derived)
       - StateEs_full[j], StateEs_hub[j] (all many-body energies)
       - P_full[j], P_hub[j] (probabilities)
       - I_full_total, I_hub_total
       - I_full_site[k], I_hub_site[k] (if per-site current defined; otherwise per-edge contributions)
       - Spread summary row: min/max per column (to show nontrivial variation).
     - Plots (.png):
       - Esite vs scan (per site) full vs hubbard
       - T vs scan (per site) full vs hubbard
       - GS energy and maybe a subset of excited energies vs scan
       - Probabilities: selected dominant states vs scan (full vs hubbard)
       - Current vs scan: total (and per-site/per-edge if available)
   - Also write a short text README in each directory listing parameters and scan range.

4) Solver calls
   - Full solver: PauliSolverCL.scan_current_tip with return_probs=True, return_state_energies=True, maybe return_curmat for current per transition (to derive per-site current if needed).
   - Hubbard path:
     - Use precalc_esite_thop (already parity matched).
     - For PME: use solve_pme_dense_batch with a full basis (for small N) so no truncation; feed Esite/Tsite and occupancy basis (all states for 2-site or 4-site). If per-site currents needed, we may need to extend kernel output; otherwise compare total current.
     - Ensure same rate conventions (Gamma, beta, mu, T) and add/remove logic consistent with full solver.

5) Nontriviality checks
   - For each scan, print min/max of Esite, T, currents, probabilities in the summary text/CSV so we verify variation (not all zeros).
   - Verify currents are not trivially zero (choose Vbias range where current flows).

6) Order of execution
   A) 2-site single-point sanity (already done) – reuse.
   B) 2-site x-scan and V-scan: generate CSV + plots + README.
   C) 4-site x-scan and V-scan: same outputs in separate directories.
   D) Final summary listing directories and parameters.

If you’re ready for me to implement and run these parity tests and emit the artifacts, switch to Code mode and I’ll proceed.

---

# Claude-4.6-Opus

Continuing from previous session. The Esite NaN was fixed by the `mpolCs` buffer size correction, but the hubbard PME dense solver still returns NaN due to a singular rate matrix. The root cause: for 2-site embedded in 4, using all 16 basis states includes states with spectator sites occupied (E=100), which produces near-zero rates and a singular matrix. Need to filter basis to only physically relevant states.

Let me read the current test harness to see how basis is constructed:


The root cause is clear: for the 2-site embedded in 4-site system, [make_basis_u4(1, nSite=4)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_runall.py:15:0-24:14) generates all 16 states including those where spectator sites (E=100 eV) are occupied. These states have near-zero rates, making the matrix singular → NaN from Gauss-Jordan.

The fix: generate basis states where only real (active) sites vary, keeping spectator sites always empty.

Now I'll fix the basis generation. The [make_basis_u4](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_runall.py:15:0-24:14) function needs to generate only states where active sites vary and spectator sites stay empty. I'll also need to handle the probability array size mismatch between full solver (always 16 states) and hubbard solver (variable basis size).

Now I need to update all callers of [make_basis_u4](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_runall.py:42:0-56:23) to handle the new return value `(out, nStates)` and pass `active_sites`. The key insight: for the 2-site embedded case, `active_sites=[0,1]` gives 4 basis states; for the 4-site case, all sites are active giving 16 states.

I also need to handle the fact that the hubbard solver's `nStates` may differ from the full solver's `nStates` (which is always 16 for 4 sites). The probability comparison needs to map between the two basis orderings.

Let me update the [run_1d_scan](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_runall.py:92:0-261:23) function:

Now update the [run_2d_scan](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_runall.py:260:0-366:18) function similarly:

Now I need to pass `active_sites` from main() to the scan functions. For the 2-site embedded case, `active_sites=[0,1]`; for 4-site, `active_sites=[0,1,2,3]`.

Let me quick-test the basis generation fix before running the full suite:

Basis generation looks correct. Now let me run the full test suite:

Fix the `nStates` reference in [run_2d_scan](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_runall.py:278:0-383:18) — I renamed it to `nStates_full` but missed one usage:

Good, the 1D scans now produce finite results! `max_dI=3.02e-06` for x-scan and `5.62e-07` for V-scan — these are the actual current values, so parity error is ~100%. Let me fix the NameError and re-run:

Good progress! 2-site tests all pass with finite values. 4-site x-scan works (`max_dI=3.18e-06`) but 4-site V-scan still has NaN. The V-scan has the tip at x=0, y=0 which is equidistant from the kite sites — the symmetry might cause a singular matrix for certain bias values. Let me check:

Only 5 NaN values at the lowest bias voltages (V ≈ 0 to 0.04). At near-zero bias, the tip chemical potential `mu_tip ≈ 0` makes the rate matrix nearly symmetric and potentially singular. The Gauss-Jordan solver's pivot threshold `1e-30` is too tight. Let me check the kernel's pivot handling and relax it slightly, or better, handle the V=0 case by clamping:
The issue is that at very low bias the Gauss-Jordan solver without pivoting fails. Let me check how PME.cl handles this:


PME.cl uses **partial pivoting**. Let me read its implementation:

PME.cl uses partial pivoting — the hubbard kernel doesn't. That's the root cause of the NaN at low bias. Let me add partial pivoting to the hubbard kernel's Gauss-Jordan solver:

Now let me re-run the full test suite with the partial pivoting fix:

All 6 tests pass with finite values! No more NaN. Let me analyze the results and check the plots:

Now let me examine the actual parity errors and check the output files:

Key observations from the errors:
- **Esite parity: PERFECT** (max_dE = 0.0 for all sites in all tests)
- **Thop parity: excellent** (max_dT ≈ 4.5e-8, essentially float32 rounding)
- **Probability parity: POOR** (max_dP ≈ 1.0 — completely different probabilities)
- **Current parity: POOR** (max_dI ≈ 3e-6, which is ~100% of the current magnitude)

The precalc (Esite, Thop) is perfectly matching, but the PME solver output (probabilities, currents) is completely different between the two solvers. This means the two solvers are computing different rate matrices or solving them differently. Let me investigate the rate matrix construction.

The key difference: the full solver (PME.cl) uses `Gamma_sub` and `Gamma_tip` directly as broadening parameters with specific conventions (lines 196 in pauli_ocl.py: `gamma_input = (gamma_val/pi)^2`), while the hubbard solver receives these squared values and uses them differently in `calc_total_rate`.

Let me compare the rate computation between the two kernels:

Now let me look at how PME.cl builds its rate matrix to understand the parity gap:
I found the root cause of the parity error. Comparing the rate matrix construction between PME.cl and hubbard.cl, there's a sign bug in the energy difference used for the removing case.

**PME.cl** always computes `ε_add = E(more electrons) - E(fewer electrons)` and uses `f(ε_add)` for adding and `(1-f(ε_add))` for removing.

**hubbard.cl** passes `dE = E_r - E_c` directly, but when `r` has fewer electrons (removing), `dE = -ε_add`. Then `(1-f(-ε_add)) = f(ε_add)`, which is the *adding* rate, not the removing rate.

The fix is to negate `dE` for the removing case so `calc_total_rate` always receives `ε_add`:
Now I need to apply the same fix to the current calculation section of the kernel, which also computes `dE` for tip current evaluation:

The current evaluation section at line 1699 uses `dE = StateEnergies[c] - StateEnergies[b]` where `c` has more electrons (`added=true`, meaning `bit_c && !bit_b`). This means `dE = ε_add` which is correct for the current calculation. Let me verify by comparing with PME.cl's current calculation:

The current calculation section in hubbard.cl (lines 1699-1704) already uses `dE = E[c] - E[b]` where c has more electrons, so `dE = ε_add` — this is correct and matches PME.cl. No fix needed there.

Now let me run the full test suite with the rate matrix fix:

**Excellent parity!** All 6 tests now show near-machine-precision agreement:

| Test | max_dE | max_dI | I_spread |
|------|--------|--------|----------|
| 2site_xscan | 0 | **1.59e-12** | [2.81e-7, 3.02e-6] |
| 2site_Vscan | 0 | **5.68e-14** | [0, 5.62e-7] |
| 2site_xV (2D) | 0 | **2.06e-11** | — |
| 4site_xscan | 0 | **2.27e-12** | [1.26e-6, 3.18e-6] |
| 4site_Vscan | 0 | **1.88e-12** | [-3.3e-14, 1.84e-6] |
| 4site_xV (2D) | 0 | **5.47e-11** | — |

Let me check the detailed error reports and view the plots:

All parity metrics are excellent. Let me view the plots to confirm they look good:

The plots look fantastic! Perfect visual overlap between full and hub curves in all panels. The 2D difference maps show residuals at the 1e-11 level — essentially float32 rounding noise.

Let me update the TODO and produce the final report:

# PME Parity Test — Final Report

## Summary

All 6 parity tests pass with **near-machine-precision agreement** between the full PME solver ([PME.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/PME.cl:0:0-0:0)) and the Hubbard dense PME solver ([hubbard.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/hubbard.cl:0:0-0:0)).

## Bugs Found & Fixed

### 1. `mpolCs` buffer size ([HubbardSolver.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py:0:0-0:0))
Buffer allocated as `sz_f*4*nSingle*nMulti` (float4 size) for plain float data. Fixed to `sz_f*nSingle*nMulti`. The oversized buffer caused rotation matrix data to be zeroed → NaN in Esite.

### 2. Active-sites basis ([test_pme_parity_runall.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_runall.py:0:0-0:0))
For 2-site embedded in 4-site, all 2^4=16 basis states included spectator sites (E=100 eV), producing a singular rate matrix. Fixed by generating only states where active sites vary (4 states for 2 real sites).

### 3. Partial pivoting ([hubbard.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/hubbard.cl:0:0-0:0))
Gauss-Jordan solver lacked partial pivoting, causing NaN at low bias. Added pivoting matching PME.cl.

### 4. **Rate matrix dE sign bug** ([hubbard.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/cl/hubbard.cl:0:0-0:0)) — the critical physics fix
For removing transitions (c→r where r has fewer electrons), `dE = E_r - E_c = -ε_add`. Passing this to `calc_total_rate(..., false)` gave `Γ·f(ε_add)` instead of the correct `Γ·(1-f(ε_add))`. Fixed with one line: `float dE_add = adding ? dE : -dE;`

## Parity Results

| Test | max|dE| | max|dT| | max|dP| | max|dI| | I range |
|------|---------|---------|---------|---------|---------|
| 2site x-scan | 0 | 4.5e-8 | 1.2e-6 | **1.6e-12** | [2.8e-7, 3.0e-6] |
| 2site V-scan | 0 | 0 | 3.4e-6 | **5.7e-14** | [0, 5.6e-7] |
| 2site 2D x-V | — | — | — | **2.1e-11** | — |
| 4site x-scan | 0 | 4.5e-8 | 2.0e-6 | **2.3e-12** | [1.3e-6, 3.2e-6] |
| 4site V-scan | 0 | 0 | 3.4e-6 | **1.9e-12** | [~0, 1.8e-6] |
| 4site 2D x-V | — | — | — | **5.5e-11** | — |

All errors are at float32 rounding level. Current spreads are nontrivial (order 1e-6).

## Files Modified
- `@/home/prokop/git/FireCore/pyBall/OCL/cl/hubbard.cl` — dE sign fix + partial pivoting
- `@/home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py` — mpolCs buffer size fix
- `@/home/prokop/git/FireCore/pyBall/OCL/pauli_ocl.py` — removed n_pixels==1 restriction on curmat
- `@/home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_runall.py` — active-sites basis, full harness

## Output Locations
- **2-site**: [/home/prokop/git/FireCore/figs/parity_2site_ruslan/](cci:9://file:///home/prokop/git/FireCore/figs/parity_2site_ruslan:0:0-0:0) (xscan/, Vscan/, xVmap/)
- **4-site**: [/home/prokop/git/FireCore/figs/parity_4site_ruslan/](cci:9://file:///home/prokop/git/FireCore/figs/parity_4site_ruslan:0:0-0:0) (xscan/, Vscan/, xVmap/)

Each contains CSV data, PNG plots, TXT summaries with all parameters, and JSON error reports.