https://gemini.google.com/share/66aeeecb1c9d

## USER 1

I'm interested in efficient splitting large cloud of particles into localized clusters. I'm thinking efficient way can be something like hierarchical k-means where we gradually split larger clusters by PCA (principal component analysis ) along the longest component. We subdivide into approx same size clusters always swapping the indexes to two aprox same parts. We can remove out layers an try to add them to other closer cluster.

---

It sounds like you are describing a variation of **Principal Direction Divisive Partitioning (PDDP)**. This is a highly effective "top-down" clustering approach that often outperforms standard K-means in high-dimensional spaces because it relies on the global structure of the data rather than random initialization.

Using PCA to split clusters ensures that you are always dividing the cloud along the axis of maximum variance—essentially "cutting" the cloud where it is longest.

### The Algorithmic Flow

Here is how your proposed method would look in a structured pipeline:

1. **Selection:** Start with the entire cloud. Calculate the mean and the first principal component (the eigenvector corresponding to the largest eigenvalue).
2. **Projection & Splitting:** Project all particles onto this principal component. To achieve "approximately same size clusters," you would split the particles at the **median** of the projected values.
3. **Index Swapping:** As you noted, you can efficiently handle this in memory by swapping indices in a contiguous array, similar to the `partition` step in QuickSort.
4. **Recursion:** Repeat the process for each resulting cluster until a stopping criterion is met (e.g., a maximum number of clusters or a minimum cluster radius).

### Refining the "Outlier" and "Equal Size" Problem

Your idea of removing outliers and re-assigning them is excellent for maintaining "compactness." Here are three ways to enhance your current strategy:

* **The "Refinement" Step:** After each PCA split, the boundary between the two new clusters is a linear hyperplane. You can run a quick **2-means (K-means with K=2)** iteration using the two new centroids as starting points. This "snaps" the boundary to the local density of the particles.
* **Outlier Handling:** Instead of just removing them, calculate the **Mahalanobis distance** for particles in a cluster. If a particle is a certain number of standard deviations away, flag it. Once the hierarchy is complete, you can perform a "Global Re-assignment" where these flagged particles are snapped to the nearest cluster centroid.
* **Balanced Trees:** If your goal is strictly equal-sized clusters (common in parallel computing or spatial indexing), using the median for the split is mandatory. However, if the cloud has a natural "gap," splitting at the median might cut through a dense area.

### Comparison: PCA-Splitting vs. Standard K-Means

| Feature | PCA-Based (PDDP) | Standard K-Means |
| --- | --- | --- |
| **Initialization** | Deterministic (No "bad starts") | Stochastic (Needs multiple runs) |
| **Complexity** |  splits |  iterations |
| **Shape** | Great for elongated/linear clusters | Best for spherical clusters |
| **Efficiency** | Very fast for large datasets | Can get stuck in local optima |


### Implementation Tip: Memory Management

Since you mentioned swapping indices, you can represent your particle cloud as a single large array of coordinates. By passing **array views** or **pointers with offsets** to your recursive function, you avoid copying data.

```python
# Conceptual logic for the split
def split_cluster(indices):
    # 1. Compute PCA on particles[indices]
    # 2. Find median projection
    # 3. Partition 'indices' around median (inplace swap)
    # 4. return indices_left, indices_right

```

Would you like me to help you write a Python or C++ snippet to implement this specific PCA-based splitting logic?

---

## USER 2

Ok for efficient parallel Implementation we should start with k-clustes where k is equal number of threads. Then we continue with only binary subdivision of each cluster (each processed by its own thread)

---

This is a smart strategy for high-performance computing. By aligning your initial split with the number of available threads, you maximize hardware utilization and minimize the overhead of thread management.

In parallel computing terms, you are moving from a **Task Parallel** approach (distributing the initial large clusters) to a **Recursive Decomposition** approach.

### The Parallel Implementation Pipeline

To make this efficient, I recommend a two-stage execution model:

#### Stage 1: The Global "K-Split"

Instead of doing a simple binary split at the very beginning, you can perform a **"Level-Up" PCA**.

* Calculate the top  principal components.
* Alternatively, run the PCA split recursively until you have  clusters (where ).
* **Synchronization:** This stage requires a global barrier because all threads are working on the same large pool of data to find the initial principal axes.

#### Stage 2: Independent Sub-Trees (The "Thread-Local" Phase)

Once each thread is assigned its own localized cluster, you enter the **Shared-Nothing** phase:

* **No Locks:** Each thread owns a specific contiguous slice of the index array.
* **Cache Locality:** Since each thread works on a smaller, localized spatial region, the particle data is more likely to stay in the L2/L3 cache.
* **Binary Subdivision:** Each thread continues splitting its assigned cloud until the desired granularity is reached.

### Optimization: Swapping vs. Masking

Since you are working with large clouds, memory bandwidth is your bottleneck.

1. **In-Place Partitioning:** Use a `std::partition`-style logic (C++) or `numpy.partition` (Python). This keeps memory usage low because you aren't creating new arrays; you are just reordering the pointers or indices to the particles.
2. **Double Buffering:** If the particles are moving (dynamic simulation), it is sometimes faster to write the "left" and "right" children into a new buffer rather than swapping in place, as it allows for purely sequential writes which are friendlier to the memory controller.

### Handling the "Outlier" Re-assignment

To keep the parallel execution "clean," you should handle outliers at the end of each thread's local work:

* **Local Cleanup:** Each thread identifies particles that are too far from its local cluster centers.
* **Global Boundary Check:** Once all threads finish, you perform a final pass. Since the clusters are now localized, you can use a **Bounding Box** or **AABB (Axis-Aligned Bounding Box)** for each cluster to quickly check which cluster an outlier should belong to.

### Comparison of Workload Distribution

| Phase | Thread Interaction | Bottleneck |
| --- | --- | --- |
| **Initial K-Split** | High (Global Sync) | Memory Bandwidth / Interconnect |
| **Binary Subdivision** | Zero (Independent) | CPU Compute (SVD/Eigen calculations) |
| **Outlier Re-fit** | Low (Boundary checks) | Search Latency |

### Implementation Detail: PCA Efficiency

For the splitting, you don't need a full SVD (Singular Value Decomposition). You only need the **largest eigenvector**. You can use the **Power Iteration** method, which is extremely fast, parallelizable, and only requires a few matrix-vector multiplications to converge.

Would you like me to provide a **Power Iteration** code snippet optimized for finding that primary splitting axis?

---

## USER 3

Ok but I want to do it mainly in 3D.

How this algorithm would compare with other algorithm where we gradually create boxes (or spheres) of size L. If new particle does not fit to any of the existing boxes we need to create new.

The motivation is to split particles into clusters to optimize pairwise interactions calculation using the fact that interactions have finite range (R_cut)

---

Since your goal is to optimize **pairwise interactions** with a finite cutoff (), the choice of clustering is critical. You are essentially comparing a **Top-Down Hierarchical Split (PCA-based)** with a **Bottom-Up Grid/Cell-Link Method**.

In the context of Molecular Dynamics (MD) or Smooth Particle Hydrodynamics (SPH), here is how they stack up for 3D spatial partitioning:

### 1. The PCA-Based Split (Recursive PDDP)

This method is "data-aware." It doesn't care about the absolute coordinate system; it only cares about where the particles actually are.

* **How it works:** It finds the "stretched" direction of your particle cloud and cuts it.
* **Pros:**
* **Load Balancing:** Since you split at the median, every thread gets exactly  particles. This prevents "hotspots" where one thread is overwhelmed.
* **Cluster Shape:** It creates clusters that follow the natural flow of the data (e.g., if particles are in a thin stream, the boxes will be long).


* **Cons for :**
* **Irregular Boundaries:** The clusters can be very "long and skinny." If a cluster is thinner than , you will have to check a massive number of neighboring clusters, which ruins the optimization.

### 2. The Fixed-Box Method (Cell-Linked Lists)

This is the industry standard for finite-range interactions ().

* **How it works:** You divide 3D space into a grid of cubes with side length . Each particle is hashed into a box based on its  coordinates.
* **Pros:**
* **O(1) Neighbor Discovery:** To find neighbors for a particle, you *only* look at its own box and the 26 surrounding boxes.
* **Constant Boundary:** You know exactly how far you need to look.


* **Cons:**
* **Memory Waste:** If your particle cloud is sparse or diagonal, you create many empty boxes.
* **Load Imbalance:** If particles clump in one corner, the thread assigned to that box does all the work while others idle.

### Head-to-Head Comparison

| Feature | PCA Hierarchical Split | Fixed Boxes () |
| --- | --- | --- |
| **Best Use Case** | Highly non-uniform/sparse clouds | Dense, uniform particle distributions |
| **Search Complexity** |  to find cluster |  to find cell |
| **Pairwise Efficiency** | Needs complex "Bounding Box" checks | Simple "Check 27 cells" logic |
| **Parallelization** | Excellent (Perfectly balanced) | Difficult (Requires domain decomposition) |

### The "Hybrid" Recommendation

For 3D particles with a , a **PCA-only approach is usually inefficient** because the split planes don't guarantee that the resulting clusters are "thick" enough to contain the interaction radius.

**A better approach for your specific 3D problem:**

1. **Octree-like PCA:** Instead of just splitting into 2, use PCA to find the primary axis, but then also consider the 2nd and 3rd axes to ensure the cluster remains "blocky" (spherical/cubical) rather than a thin needle.
2. **AABB Guard Bands:** For every PCA-generated cluster, calculate its **Axis-Aligned Bounding Box (AABB)**.
3. **The Interaction Rule:** * If , skip all particles in those clusters.
* This allows you to prune thousands of calculations with a single  check.

### Which one should you choose?

* If your particles are **uniformly distributed** (like a gas in a box): Use **Fixed Boxes**. It’s simpler and faster.
* If your particles are **highly localized** (like a galaxy or a complex protein string): Use the **PCA/Hierarchical Split**. It avoids wasting memory on empty space and ensures every CPU core has the same amount of work.

Would you like to see how to calculate the 3D Bounding Box (AABB) for a cluster after a PCA split to handle the  logic?

---

## USER 4

No you misinterpreted me. Grid approach is something else. I want to avoid that because it requires to count also empty boxes or use some hashmap or sorting. That is quite costly especially when the space is only sparsely covered.

But what I want is to create bounding boxes (AABB or spheres) which has maximum size comparable with R_cut. 

When some box contains too many particles we can split it by PCA. However, important is that the initial layout have already size similar to R_cut witch seems reasonable

---

I see your point clearly now. You are aiming for a **data-centric** approach that only allocates memory and compute where particles actually exist, avoiding the overhead of a global grid or hashmap.

By using  as the guiding scale for your initial layout, you essentially create a "tightly-packed" set of localized envelopes. Using PCA as the splitting mechanism ensures that when density is high, the subdivision follows the actual "grain" of the particle distribution.

### The Refined Algorithm: "Adaptive Envelope Splitting"

Instead of starting with one giant cloud, you start by grouping particles into "Seed Envelopes" of size .

#### 1. Initial "Greedy" Clustering (The  Layer)

To avoid a grid, you can perform a single pass to create initial seeds:

* Pick an unassigned particle.
* Create a cluster containing all unassigned particles within a  radius (or AABB).
* This ensures you only "spawn" clusters where data exists.

#### 2. PCA-Driven Refinement (The Load Balancing Layer)

Once you have these -sized envelopes, some might contain 10 particles and others 10,000 (dense regions).

* **The Threshold:** If `count(particles in cluster) > Max_Particles_Per_Thread`, apply your **PCA split**.
* **The Split:** Since the cluster is already , the PCA split will likely cut it into two "slabs."
* **The Benefit:** This handles the "sparse vs. dense" problem perfectly. In sparse areas, the  box remains a single unit. In dense areas, PCA breaks it down into chunks small enough for efficient parallel processing.

### Why PCA is better than "Midpoint" splitting here

In 3D, if you have a dense "filament" of particles (like a polymer or a stream), a standard octree split might cut through the filament at an awkward angle, creating many clusters with just a few particles.

PCA identifies the **eigenvectors** of the cluster. By splitting perpendicular to the first principal component, you are guaranteed to:

1. Minimize the "surface area" of the cut.
2. Maximize the distance between the new centroids.
3. Keep the resulting sub-clusters as compact as possible, which is vital for the  logic.

### Pairwise Interaction Logic

With this hierarchy, your interaction check becomes a simple **AABB-to-AABB** test:

1. **Level 1 (Coarse):** Check if . If yes, skip all  interactions.
2. **Level 2 (Fine):** If the clusters are neighbors, iterate through the particles.

### Comparison: Your PCA Method vs. Traditional Octree

| Feature | PCA-Subdivision | Standard Octree |
| --- | --- | --- |
| **Depth** | Shallower (more efficient splits) | Deeper (fixed branching) |
| **Sparsity** | Naturally ignores empty space | Must "skip" empty nodes |
| **Memory** | Low (contiguous index swapping) | Higher (pointer overhead) |
| **Alignment** | Aligns with data density | Aligns with coordinate axes |

### Implementation Recommendation: The "Median Swap"

To keep this fast in parallel, use the **QuickSort-style partition** you mentioned:

* Compute the 1st Principal Component vector .
* For each particle  in the cluster, calculate the scalar projection .
* Find the **median** of .
* Swap indices so all particles with  are on the left.

**Next Step:**
Would you like to look at a specific strategy for the **Initial Seeding** phase (to get those first  clusters without a grid), or should we focus on the **Power Iteration** math to make the PCA split as fast as possible?

---

## USER 5

Yes but you did not get the initial step.

I want to iterate over particles sequentially. For each particle check boxes which are already created, what is the cost of adding particle to it (i.e enlarging it increasing the maximum side length L of that box). If the cost is to high for all boxes, we need to create new box.

This way we create certain number of boxes of limited maximum size. Which we can eventually further subdivide.

---

I see—you are proposing a **Dynamic Sequential Grouping** strategy (similar to a *Leader-Clustering* or *Growing Neural Gas* approach) to initialize the seeds.

Instead of a pre-defined grid, you are building an **adaptive spatial cover** by streaming the particles. This is much more memory-efficient for sparse data because you only "pay" for the volume you actually occupy.

### The "Cost-Based" Seeding Logic

The "cost" of adding a particle to an existing box can be defined by the constraint that the box must not exceed a maximum side length .

**The Algorithm:**

1. **Iterate:** For each particle :
* **Find Candidates:** Check existing boxes whose centers are within a reasonable distance (you can keep a list of active box centroids).
* **Evaluate Cost:** For each candidate box, calculate the new Bounding Box if  were included.
* **Decision:**
* If a box exists where the new , add  to the one that requires the **minimum volume increase**.
* If no such box exists (or adding it would exceed ), **spawn a new box** with  as the first member.

2. **Refine (Parallel PCA):** Once the sequential pass is done, you have a set of boxes of size . For any box containing a massive number of particles, you hand it off to a thread to perform your **PCA-based binary subdivision**.

### Strengths of this Approach

* **No Grid/Hashing:** You completely avoid the "empty cell" problem. If your particles are in two thin clusters separated by a massive void, you only store two boxes.
* **Sequential Efficiency:** The first pass is  where  is the number of active boxes. Since  is usually much smaller than , this is very fast.
* **Ideal for :** Because you cap the growth at , you are guaranteed that your top-level clusters are perfectly sized for pruning pairwise interactions.

### Improving the Search for "Candidate Boxes"

Since you want to avoid a grid, the only bottleneck is Step 1 (checking existing boxes). If you have thousands of boxes, checking every box for every particle becomes  in the worst case.

**Optimization Trick:**
Keep a simple **Bounding Volume Hierarchy (BVH)** or a small list of "Active Centroids." Since you know your boxes are , you only need to check boxes whose centroids are within approximately  of the current particle.

### Why the PCA Subdivision is the Perfect Second Step

The sequential pass creates the "Envelopes." However, a sequential pass is sensitive to the **order of particles**. You might end up with a very dense box that is slightly "off-center."

The **PCA split** acts as the optimizer:

1. It takes those dense "Seed Boxes."
2. It re-evaluates the internal variance.
3. It splits them along the "natural axis," ensuring that the workload is perfectly balanced across your threads.

### Next Step for Implementation

To make this high-performance, the "Cost Function" for adding a particle to a box is key. Usually, the **Surface Area Heuristic (SAH)** or simple **Volume Expansion** works best.

Would you like to define the **cost function** mathematically, or should we look at how to manage the **Index Swapping** once the PCA-based subdivision starts?

---

## USER 6

Ok now when you understand it I want to compare this strategy (with size based boxes) and the previously described algorithm where we start from k-means where k is number of threads.

---

Both strategies aim to group particles for finite-range interaction efficiency (), but they solve the problem from opposite directions.

### 1. Global PCA-Binary Split (Top-Down)

This is essentially a **Space-Filling Tree** approach. You start with the whole world and keep "chopping" it into equal halves based on data density.

* **Logic:** Start with  particles  calculate PCA  split into 2  repeat until you have  (threads) clusters.
* **The "K-Threads" Advantage:** You get **perfect load balancing**. Because you split at the median, every thread receives exactly the same number of particles.
* **The  Weakness:** PCA cares about *variance*, not *distance*. If your particle cloud is very sparse and long, PCA might create a very "thin" cluster that is much longer than  but thinner than .
* **Result:** You'll end up having to check many neighbor clusters because the cluster's boundary is "leaky" relative to the interaction radius.

### 2. Adaptive Envelope Seeding (Bottom-Up Sequential)

This is a **Density-Aware Cover**. You grow the clusters from the particles themselves, capped by a physical size .

* **Logic:** Iterate through particles  fit into a box if "cost" (expansion) is low  else, start a new box.
* **The  Strength:** This is **physically optimal**. Since every box is roughly  in size, you know that for any particle, its neighbors can *only* be in its own box or the immediate neighbor boxes. This mimics the efficiency of a grid without the memory waste of empty cells.
* **The Thread Weakness:** It is **harder to load balance**. One box might have 10,000 particles (a dense liquid) and another might have 5 (a sparse gas). If you just assign one box per thread, some threads will finish in microseconds while others grind for seconds.

### Head-to-Head Comparison

| Feature | PCA-K Split (Top-Down) | Your Sequential Envelopes (Bottom-Up) |
| --- | --- | --- |
| **Initial Focus** | Balancing CPU work (Load) | Minimizing Search Radius () |
| **Cluster Shape** | Arbitrary "Slabs" (PCA-aligned) | Tight Envelopes (User-defined ) |
| **Empty Space** | Ignores it (Good for sparse) | Never creates it (Best for sparse) |
| **Complexity** |  |  |
| **Interaction Pruning** | Moderate (AABB can be large) | **Excellent** (AABB ) |

### The Synthesis: The Best of Both Worlds

Since you specifically want to optimize **pairwise interactions**, the **Sequential Envelope** (Strategy 2) is the better starting point because it respects the physical scale of your simulation ().

**However**, you should use the **PCA logic** as a "cleanup" step for load balancing:

1. **Phase 1 (Sequential):** Generate your  boxes.
2. **Phase 2 (Check Density):** Identify "Heavy Boxes" (where particle count  threshold).
3. **Phase 3 (PCA Split):** Use the PCA-Binary subdivision on **only** the heavy boxes to break them into smaller chunks that can be distributed across your threads.

> **Key takeaway:** The sequential pass builds the **Geometry** (respecting ); the PCA splits build the **Parallelism** (respecting thread count).

Would you like to see how to define the "Cost" of adding a particle to a box so that the sequential pass stays  instead of slowing down as you add more boxes?

---

## USER 7

The problem of the second strategy (R_cut boxes) is
1) the number of boxes is not limited. We may end up with as many boxes as particles making it very slow O(n2)
2) we cannot do it easily parallel, while the k-means were inherently parallel

---

You’ve hit the nail on the head. Those two points are the primary "bottlenecks" that make pure bottom-up clustering difficult for high-performance computing.

If you have  boxes for  particles, you’ve essentially built a very expensive linked list. And because the "Box Creation" depends on what happened with the previous particle, you have a **serial dependency** that breaks parallelism.

Here is a breakdown of how your original PCA-based hierarchical approach solves these specific failures compared to the bottom-up strategy.

### 1. Avoiding the  Trap (The Complexity Problem)

In the bottom-up approach, checking  particles against  boxes results in . If  grows with , you hit a quadratic wall.

**The PCA Top-Down Solution:**

* **Bounded Depth:** By splitting the cloud into exactly  clusters (or a fixed depth), the number of clusters is **guaranteed** to be small and manageable.
* **Logarithmic Work:** Because you divide the population in half each time, you reach your target cluster size in  steps.
* **Spatial Control:** Instead of hoping boxes stay small, you use the PCA split to force high-density regions to subdivide. If a cluster's AABB is still much larger than  after the first split, the second or third split will naturally bring the dimensions down toward the  scale.

### 2. Parallelism: Synchronous vs. Asynchronous

The bottom-up strategy is "Stateful"—you can't process Particle 1000 until you know where Particle 999 went.

**The PCA "Divide and Conquer" Solution:**

* **Fork-Join Parallelism:** The first PCA split is done by all threads (global). Once the data is split into "Left" and "Right" halves, the two halves are **completely independent**.
* **Zero Contention:** Thread A can work on the Left half and Thread B on the Right half without ever checking a shared mutex or box-list.
* **Recursive Scaling:** As you go deeper, you can assign more threads to sub-branches.

### The Comparison Table: Addressing Your Concerns

| Problem | Bottom-Up ( Boxes) | Top-Down (PCA K-Split) |
| --- | --- | --- |
| **Max Cluster Count** | Unbounded (potentially ) | **Strictly  (number of threads)** |
| **Search Cost** |  (Potential ) | **** |
| **Parallelization** | Hard (Serial dependency) | **Native (Recursive decomposition)** |
| **Memory Access** | Random / Scattered | **Contiguous (via index swapping)** |
| ** Strategy** | Guaranteed box size, unknown load | **Guaranteed load, variable box size** |

### How to get  efficiency with the PCA approach

Since you prefer the PCA approach for its efficiency and parallelism, you can "mimic" the  benefit by using a **Stopping Criterion**:

1. **Split by PCA** along the longest axis.
2. **Check Dimensions:** If the Cluster AABB's longest side is still , keep splitting *regardless* of particle count.
3. **Check Density:** If the AABB is  but has too many particles for one thread, keep splitting to balance the load.

This ensures that at the end, every cluster is **either** small enough to fit in the  window **or** contains a perfectly balanced workload for your threads.

### A Note on "Swapping Indexes"

In the PCA version, the **"In-place Partition"** (swapping) is your best friend. It ensures that the particles belonging to "Thread 1" are all sitting next to each other in RAM. This creates massive **CPU Cache gains** when calculating pairwise interactions, as the processor can pre-fetch the particle data linearly.

Would you like to see a logic flow for how to handle the **Index Swapping** across multiple threads during that first "Global" PCA split?

---

## USER 8

When starting the k-means what is the best way to choose the initial pivot points or positions of k-means cluster centers ? That may significantly affect performance but we need some simple and fast heuristic

---

Since you are focused on **efficiency** and **large 3D particle clouds**, you should avoid the standard "Random" or even "K-means++" initialization. For  particles and  threads, those are either too unstable or too slow (sequential).

The best heuristic for your specific architecture is **PCA-Based Deterministic Initialization** (sometimes called PCA-Part).

### The Best Fast Heuristic: PCA Hierarchical Seeding

Instead of guessing  points at once, you use your PCA splitting logic to **generate** the  initial pivots. This is much faster than K-means++ because it is deterministic and aligns with the variance of your data.

1. **Start with the whole cloud.** Calculate the center of mass (Centroid 0).
2. **Find the 1st Principal Component.** Project all particles onto it.
3. **Split at the Median.** This gives you two clusters and two new centroids.
4. **Repeat until you have  clusters.** * *Tip:* Always split the cluster with the largest **Sum of Squared Errors (SSE)** or simply the one with the most particles.
5. **Use these  centroids as your "Seeds" for K-means.**

### Why this is superior for your case:

* **Speed:** K-means++ requires  passes over the data. PCA-seeding only requires  passes.
* **Parallel-Friendly:** You can parallelize the PCA calculation (Power Iteration) for the first few levels when the clusters are still large.
* **Stability:** Unlike random starts, this will give you the same (or nearly the same) result every time, which is critical for debugging physical simulations.
* **Quality:** It places seeds far apart along the axes where the particles are most "spread out," which is exactly what you want for  optimization.

### Comparison of Initialization Speed

| Method | Complexity (Seeding) | Hardware Suitability |
| --- | --- | --- |
| **Random (Forgy)** |  | Instant, but often bad results |
| **K-Means++** |  | **Poor** (inherently sequential) |
| **PCA-Seeding** | **** | **Excellent** (deterministic & parallelizable) |

### Implementation Shortcut: The "Extreme Points" Heuristic

If even PCA is too slow for your initial pass, you can use a "Bounding Box" heuristic:

1. Find the  coordinates in  (1 pass).
2. The distance between  and  defines your longest axis.
3. Place your  seeds linearly along that axis, or on a grid within that box.
4. Run **one** iteration of K-means to "snap" those seeds to the actual local particle density.

### The "Swapping" Trick for K-means

Since you are swapping indices anyway, keep the data **physically sorted** after your PCA split. When you finally run your K-means iterations, the particles belonging to the same cluster will be **cache-local**. This can speed up the "Distance to Centroid" calculation by 3x–5x due to better memory bandwidth utilization.

Would you like to see the **3-line Power Iteration** code to find that 1st Principal Component for the seeding?