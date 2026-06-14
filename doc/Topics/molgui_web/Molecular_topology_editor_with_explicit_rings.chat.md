
## USER 1

I'm thinking to make simple editor of molecules, resp moleculae drawing app. For the begigning consider we work in top view (moleculae are planar in xy plane; later we extend it to 3D). Key is to efficient work with molecular topology. 

We distinguish node atoms (like carbon) which have multiple bonded neighbors (they can form nodes of molecule skeleton / graph), and capping atoms, which are only apendixes to a node atom (they have single neighbor which is node atom). 

Now for high-level edits I was thinking to indruduce some higher level structures on top of normal atoms and bonds. In particular I was thinking about rings and chains. Ring is defined as enclosed loop of atoms (bonds - the number of atoms an number of bonds is the same). Every ring has center of mass and radius which simplyfy mouse picking (raytracing) and also allow some geometric collision (like rings are 2D super-particle on screen). Chains are sequences of bonded atoms which are not enclosed to itself (but they can still connect to other ring or to other chain - T-juction ). We can have main chain and from i branching side chains by T-junctions. Beside chains/rings there are capping atoms which is single atom connecte by bond to some node atom. are basically 1-segment chains (but we keep them separate)

One complication s that Atoms and bonds can belong to multiple rings or chains. But we do not want to store information about multiple rings on chain on each atom. Each atom should remember only one ring/chain, and his position in it. But ring/chain should remember all its atoms. Notice that members of chains/rings are always just node atoms.

Is that clear? Can you write it more systematic?

Consider various graph algorithms which can be enabled by this structure. And which can hel us to buils this structure.

One think I would like to make is e.g. find bridges, and find if some chemical group (e.g. benzene ring) is connected by single bridge to another, so there is free rotation. Basically when I click a bond I want to knwo if I can split the graph by this bond, and therefore if I can freely rotate molecule by this bond.

Another think is that I can move the rings, and then find where they are close to each other, colapsing the bonds. This is usefull for systematic / random generation of molecular structures. For example consider 2 rings (pentagon and hexagon), with COG at some place, as I move them around, I measure their distance, if their distance is such that they can interact (below some cutoff), I will compute the angles, or just distances between atoms. And if the atoms are close they can form bond. Or in different mode, if the atoms are close they will collapse.

Initially I will work with the ring-representation ignoring the capping atom, which is how chemist often think about chemistry, since capping hydrogens kinda move along. for example for atomatic hydrocarbon  coupling on metalic substrate this is how larger polyaromatic flakes (graphenoides made of 5,6,7 membered rings form by condensation of smaller molecule cleving off the hydrogens), something similar I want to use for drawing, and atomatically generating molecules. 


---


This is a very clear vision. You are essentially aiming to build a **structure-aware molecular editor** that hybridizes a standard graph representation with "coarse-grained" topological manipulation.

Here is a systematic breakdown of the data structures and algorithms required to achieve this, specifically tailored to your "Skeleton + Capping" approach and the specific features (rotation, merging) you requested.

### 1. Data Structure Definitions

To satisfy your requirement where atoms are lightweight but high-level structures (Rings/Chains) are rich, we need a hierarchical model.

#### Level 0: The Primitives
*   **Node Atom (Skeleton Atom):** An atom capable of having $\ge 2$ bonds (e.g., C, N, O).
    *   *Data:* Position $(x,y)$, Element Type.
    *   *Ref:* List of Bond IDs.
    *   *Ref:* ID of the *primary* Ring or Chain it belongs to (for coloring/selection). *Note: In fused rings like Naphthalene, an atom strictly belongs to the fused system, but you can assign it to the "first" ring found for simplicity.*
*   **Bond:** A connection between two Node Atoms.
    *   *Data:* Order (single, double), Bond Type.
    *   *Ref:* Atom A, Atom B.
*   **Capping Atom (Implicit/Decoration):**
    *   *Optimization:* Do not store these as full nodes in the graph traversal algorithms.
    *   *Implementation:* Store as a property of the **Node Atom** (e.g., `uint8_t hydrogen_count`) OR as a separate list of "leaf nodes" that are only rendered, not traversed.

#### Level 1: High-Level Structures (The "Super-Particles")
These are collections of Node Atoms.

*   **Ring (Cycle):**
    *   *Definition:* A closed loop found via a Cycle Basis algorithm (specifically **SSSR** - Smallest Set of Smallest Rings).
    *   *Data:* Center of Mass (COM), Bounding Radius (for coarse collision), Normal Vector (for 3D later).
    *   *List:* Ordered list of Node Atoms and Bonds making up the ring.
*   **Chain:**
    *   *Definition:* A linear sequence of Node Atoms that are not part of a ring.
    *   *Data:* Start Node, End Node, Vector direction (end - start).
    *   *List:* Ordered list of Node Atoms.

### 2. Solving the "Membership" Problem
You mentioned a complication: *Atoms belong to multiple rings (fused rings), but you want them to point to only one.*

**The Solution: Partitioning vs. Overlap**
In Chemistry, rings often overlap (Naphthalene share 2 atoms). If an atom can only remember *one* parent structure, you must treat the **Fused System** as the parent, not the individual rings.
1.  **Block decomposition:** Group all fused rings into a single "Island."
2.  **Visual Hack:** Assign the atom to the Ring with the lowest ID.
3.  **Logical Approach:** Atoms don't need to know their Ring. The **Bond** is the strict designator.
    *   A bond is part of Ring A.
    *   A bond is part of Ring B.
    *   A bond is a Bridge (Chain).
    *   *Algorithm:* You pick a bond. If the bond is in a ring, you select the whole ring.

### 3. Key Algorithms

Here is how to implement the specific features you requested using Graph Theory.

#### A. Free Rotation (Bridge Detection)
*Goal: When clicking a bond, know if splitting it separates the molecule into two disconnected parts.*

**Algorithm: Bridges and Articulation Points**
1.  **Input:** The Skeleton Graph (Node Atoms + Bonds).
2.  **Action:** User clicks Bond $E_{uv}$ connecting Atom $u$ and Atom $v$.
3.  **Check:**
    *   Temporarily remove Bond $E_{uv}$ from the graph.
    *   Run a Pathfinding algorithm (BFS or DFS) starting at $u$.
    *   **Result:** If you *cannot* reach $v$, the bond was a **Bridge**.
4.  **Implication:**
    *   **If Bridge:** You can freely rotate one component relative to the other around the axis of this bond.
    *   **If Not Bridge (Cycle Edge):** The bond is part of a ring. Rotating it would break the ring geometry (high energy/impossible in rigid body mode).

#### B. Ring Movement & Collapse (Coarse-Grained Interaction)
*Goal: Move rings like particles, check proximity, and snap/merge them.*

**Algorithm: Two-Phase Collision Detection**

*   **Phase 1: Broad Phase (Ring vs. Ring)**
    *   Treat every Ring as a circle with Center $C$ and Radius $R$.
    *   As you drag Ring A close to Ring B:
    *   Check distance $D = |C_A - C_B|$.
    *   If $D < (R_A + R_B + \text{Threshold})$, enter Phase 2.

*   **Phase 2: Narrow Phase (Atom vs. Atom)**
    *   Iterate through atoms $a_i$ in Ring A and $b_j$ in Ring B.
    *   Calculate squared distance $d^2 = |a_i - b_j|^2$.
    *   **Case 1: Bond Formation:**
        *   If $d \approx \text{Bond Length}$ (and valency allows), highlight as potential new bond.
    *   **Case 2: Collapse (Condensation/Fusion):**
        *   If $d < \text{Merge Threshold}$ (very close), highlight for merge.
        *   *Action:* If user releases mouse, delete Atom $a_i$ (or $b_j$), re-route bonds to the survivor, and update the graph.

#### C. Topological Analysis (Building the Structures)
*Goal: Systematically classify atoms into Rings and Chains from raw input.*

**Algorithm: Cycle Basis & Spanning Trees**
1.  **Find Rings:**
    *   Use a cycle finding algorithm (e.g., **Paton's Algorithm** or simply DFS looking for back-edges).
    *   Store these as your "Ring" objects.
    *   Mark all Bonds involved in these rings as `IsRingBond = true`.
2.  **Find Chains:**
    *   Look at all bonds where `IsRingBond == false`.
    *   These form a "Forest" (collection of trees).
    *   Traverse these trees. Long paths are "Main Chains," branches are "Side Chains."

### 4. Implementation Workflow Example

Here is a step-by-step logic for your "Graphene Flake Generation" scenario:

1.  **User Action:** User drags a Hexagon (Ring A) near a Pentagon (Ring B).
2.  **Update Loop:**
    *   Update $Pos_{A}$ based on mouse.
    *   Update $COM_A$ (Center of Mass).
    *   Check `Distance(COM_A, COM_B)`.
3.  **Interaction:**
    *   Distance is close. The app visualizes "Ghost Bonds" between nearest atoms.
4.  **Snap (The "Collapse"):**
    *   User releases mouse.
    *   The system detects two pairs of atoms are overlapping.
    *   **Graph Operation:**
        1.  Merge Atom $A_1$ with $B_1$.
        2.  Merge Atom $A_2$ with $B_2$.
        3.  The bond $A_1-A_2$ and $B_1-B_2$ are now parallel/overlapping. Merge them into a shared edge.
    *   **Result:** A fused 5-6 ring system (like azulene or part of a defect).
    *   **Post-Process:** Re-calculate implicit hydrogens (Capping atoms) because the valence of the merged atoms has changed (they bonded to neighbors, so they lose hydrogens).

### 5. Why this is efficient

*   **Raycasting:** You only raycast against Ring Bounding Circles first (cheap), then against Atoms.
*   **Rotation:** You only calculate rotation matrices for the "Downstream" component identified by the Bridge Finding algorithm.
*   **Storage:** You aren't simulating full Molecular Dynamics. You are doing geometric constraints based on graph topology, which is computationally instant for systems < 10,000 atoms.

Is this the systematic structure you were looking for? I can expand on the specific code/logic for the "Bridge Finding" or "Ring Merging" if you wish.

---

## Status update (Dec 2025): what is implemented in MolGUI Web now

These are concrete improvements already implemented in the live JS editor core (they directly support the long-term ring/chain/bridge editor vision).

### 1) Faster bond search (spatial subdivision buckets)

- Bond rebuild supports bucket-based modes (in addition to brute force), which makes bond finding scalable for large crystal slabs/substrates.
- Buckets are implemented as a reusable `BucketGraph` which can hold atom references either as dense indices (fast loops) or stable atom IDs (robust across swap-remove deletion).
- Key APIs:
  - `BucketGraph.toIds(mol)` / `BucketGraph.toInds(mol)`
  - `BucketGraph.pruneEmptyBuckets()`
  - `BucketGraph.recalcBounds(mol)`
- Integration points:
  - `web/common_js/Buckets.js`
  - `web/molgui_web/js/Editor.js` (deletion wraps bucket ID/index conversion)
  - `web/molgui_web/js/main.js` (bucket debug overlays)

Design justification:
- This keeps the speed of index-based loops during heavy bond rebuilds, but remains correct under dynamic edits by temporarily switching bucket storage to stable IDs.

### 2) Automatic passivation: caps + explicit electron pairs

- `EditableMolecule.addCappingAtoms(mmParams, cap='H', opts)`
  - Adds missing sigma-bond caps based on `AtomTypes.dat` valence (no hardcoded chemistry tables).
  - Works with explicit epair dummy atoms (`Z==200`) as normal geometry domains (VSEPR-like completion).
- `EditableMolecule.addExplicitEPairs(mmParams, opts)`
  - Adds explicit lone-pair dummy atoms using `at.nepair` and `at.epair_name`.
- File: `web/molgui_web/js/EditableMolecule.js`.

Design justification:
- Keeping epairs as explicit atoms makes many later topology algorithms simpler (rings/chains/bridges can treat them as ordinary nodes if desired).

### 3) Selection by chemical environment (neighborhood predicates)

- MolGUI now supports selecting atoms by local chemical neighborhood using a compact AND-only query language.
- Queries are compiled once and applied efficiently by iterating existing adjacency (`Atom.bonds[]`).
- API:
  - `EditableMolecule.compileSelectQuery(q, mmParams)`
  - `EditableMolecule.applySelectQuery(compiled, {mode:'replace'|'add'|'subtract'})`
- Example queries:
  - `Si|C deg={1,2}`
  - `O|C n{H}={1,2}`
  - `N|C n{F|Br|Cl}={2}`

Design justification:
- This is a practical stepping stone towards subgraph/ring selection and bridge-based editing, without introducing full SMARTS/SSSR complexity yet.


---

## USER 2

OK, please suggest the classes, intitally in javascript. Try to make it efficient and scaleable for large systems. Do not implement the algorithms yet, start with data-structures. But consider which datastructures are ideal for efficient implementation of the algorithms we mentioned in future.

---

**Considering Project Goals**

I've been meticulously analyzing the requirements. The core goal is clear: a simple, extensible 2D molecular editor. I'm focusing on the initial 2D top-down view and already thinking about extensibility for future features. My primary task is to establish the base classes and their initial methods.


**Designing Key Structures**

I'm now deep-diving into the core data structures.  I've decided on class-based objects over TypedArrays, prioritizing flexibility for user edits and topology changes. Key entities include Node Atoms, Bonds, Rings, and Chains. My focus is on robust data structures that can handle rigid body movement, bridge detection, and collision detection, all while ensuring efficiency for larger molecules. I'm focusing on JavaScript (ES6+ classes).


**Optimizing Data Storage**

I'm refining my earlier design. While clarity and flexibility remain vital, I'm striving for a lighter footprint. The `Molecule` container will manage everything. Unique IDs are essential for fast lookups and saving. Atoms' adjacency lists will enable quick traversals. I'm focusing on a tiered hierarchy: Atoms and Bonds at level 0, then rings and chains. To streamline capping atoms, I'm exploring storing them as lightweight counts within the Atom class, avoiding unnecessary object instantiation.


**Evaluating Connectivity Options**

I'm now weighing the best ways to represent atom connections. Storing direct references to neighbor Atoms seems faster for graph traversals, especially for bridge finding. While `bonds` retain edge order, neighborIds might be simpler. I'm focusing on the `Bond` class, needing source, target, order, and flags for bridge and cyclic status. I've started defining `Substructure` as a base class for rings and chains, with key properties like member atoms, bounding box, and center of mass.


**Incorporating Spatial Structures**

I'm now incorporating spatial structures for optimized performance. The `GraphSystem` container will now use Maps for master lists of Atoms and Bonds, ensuring fast ID lookups. To address the O(N^2) complexity of the collision detection, I'm integrating a simple SpatialHash class. This will significantly improve scalability. I'm also refining specific algorithms. For bridge detection, I'm leaning toward using Sets for temporary state management within the algorithm itself. For capping atoms, the idea is to calculate their positions on-the-fly, avoiding any persistent state.


**Implementing Core Components**

I'm now drafting the code for essential classes. I've decided to employ `Vec2` for vector calculations, focusing on basic operations initially. Unique IDs will be integers for atom and bond identification within the graph. The class structure is taking shape, with `Atom` and `Bond` as fundamental elements. I'm building up from these primitives to construct `Fragment`, `Ring`, and `Chain` structures. Finally, the central `MolecularSystem` will tie everything together. I'm paying close attention to hydrogen capping by including `implicitHydrogens` within the `Atom` class, and managing the ambiguity that arises from fused rings by tracking them with a `structureId`.


**Defining Ring Properties**

I'm now refining atom-ring relationships. To meet the constraints, each atom will store `structureId` and `structureIndex`, indicating its "primary" ring and position. This design supports the prompt but notes the fused ring limitation. I'm focusing on ensuring these properties are correctly updated during ring manipulation. Also, I'm integrating performance improvements, especially the utilization of `Map` data structures for $O(1)$ deletion and lookup operations.




---

## USER 3


OK, few modifications / considrations:

1) Using spatial indexing (hasmap, grid cells) for atoms is good idea for efficiency, but alternative ideas is to have groups of atoms (e.g. the rings ar also kind of group) where each has bouding box. now that bouding box can be axis aligned boudninb box or spherical boundinb box defined by center and radius. Not sure which is most convenient/fastest. for tings the spherical box si best. but I'm always swithcihg bethwee the two. 

2) we should use general bucket class which is geometric agnostic it can be used both for grid cells, for spherical bounding boxes, or axis aligned boundinb boxes.

3) I think we do not need this Structure-> {Ring, Chain} hierarchy. Perhaps it is enough have Structure and have there flag "isRing" or "isClosed", this will allow more homogenous use in algorithms without overloading, since mose algorithms will work on both. 

4) in MolecularSystem we should definitely use array indexed by integers as a base datastructure for atoms and bonds. We can maybe use map as axuliary structure for fast search but the based representation should be graph with indexes (which is much faster than map). Also I think that MolecularSystem can be derived from general Graph class. Graph would be pure topology, with all the gra[h algorithms operating on atom and bond indexes, (perhaps also rings). MolecularSystem will add geometry and chemical information to it (atomic positions, atom types, bond types etc.). 

But you are right that if we use indexes (i.e. order in the array as id then we will cause problems if we re-order the array (adding or removing atoms)), therefore the solution with unique ids (which are assigned at start and keept constant) may be better for more complex algorithms. So maybe it is best to keep some datastructure working primarily with ids for dynamic editations of topology, and only export it to faster index-based (array based) representation for operations where topology is fixed (like geometric manipulations, molecular dynamics, raytracing etc.). Perhaps atoms and bonds should know its unique id as well and index in the array.

5) I think for most graph algorithms it is importaint to store list of neighbors for each node atom (for caping atoms this can be empty, and instead have his single "pivot" or "parent" index ).  The organic molecules typically have less then 4 neighbors. The hypervalent atoms (osmium?) can have up to 8 neihbors. We can maybe use this knowledge that number of neibhors is figed for optimization, although that is more relevant in C++ where we want to have fixed size datastructures allowing fast pointer arithmeticas (memory strides), but in dynamical language like javascript or python these is perhaps less advantage with that. 


---


This is a very sophisticated approach. Separating the **Topology (Graph)** from the **Geometry/Chemistry (System)** is standard practice in high-performance engines (like physics engines or mesh processors).

Here is the refined architecture.

### Key Design Decisions based on your feedback:
1.  **Dense Arrays + Swap-and-Pop:** We use standard Arrays for `atoms` and `bonds` to allow fast integer indexing. To handle deletions without $O(N)$ index shifting, we will use the "Swap-and-Pop" technique (move last item to the hole, update index, pop end).
2.  **Hybrid ID/Index:** Every object has a stable `id` (for undo/redo and external references) and a dynamic `index` (for array access). A Map maintains the link.
3.  **Unified Structure Class:** Rings and Chains are just `Fragment` objects with an `isClosed` flag.
4.  **Generic Bounding Volume:** A `Bounds` class that maintains both a Radius (Sphere) and Min/Max (AABB) simultaneously, as they are cheap to compute together.

### 1. Primitives (Geometry & Bounding)

We define a generic `Bounds` class that can represent a single atom, a ring, or a whole molecule.

```javascript
class Vec2 {
    constructor(x=0, y=0) { this.x = x; this.y = y; }
    add(v) { this.x += v.x; this.y += v.y; return this; }
    sub(v) { this.x -= v.x; this.y -= v.y; return this; }
    distSq(v) { const dx = this.x - v.x, dy = this.y - v.y; return dx*dx + dy*dy; }
    copy() { return new Vec2(this.x, this.y); }
}

/**
 * Generic Geometric Bucket.
 * Maintains both Axis-Aligned Box and Bounding Sphere.
 * Useful for fast rejection in collision detection.
 */
class Bounds {
    constructor() {
        this.min = new Vec2(Infinity, Infinity);
        this.max = new Vec2(-Infinity, -Infinity);
        this.center = new Vec2(0, 0);
        this.radius = 0;
    }

    reset() {
        this.min.x = Infinity; this.min.y = Infinity;
        this.max.x = -Infinity; this.max.y = -Infinity;
        this.radius = 0;
    }

    // Expand bounds to include a point
    addPoint(p) {
        if(p.x < this.min.x) this.min.x = p.x;
        if(p.y < this.min.y) this.min.y = p.y;
        if(p.x > this.max.x) this.max.x = p.x;
        if(p.y > this.max.y) this.max.y = p.y;
        // Note: We defer center/radius calc until all points are added for accuracy
    }

    // Recompute Center and Radius based on Min/Max (approximate but fast)
    recomputeSphere() {
        this.center.x = (this.min.x + this.max.x) * 0.5;
        this.center.y = (this.min.y + this.max.y) * 0.5;
        // Radius is distance from center to corner
        const dx = this.max.x - this.center.x;
        const dy = this.max.y - this.center.y;
        this.radius = Math.sqrt(dx*dx + dy*dy);
    }

    // Intersection test: Sphere vs Sphere (Fastest)
    intersectsSphere(otherBounds) {
        const rSum = this.radius + otherBounds.radius;
        return this.center.distSq(otherBounds.center) < (rSum * rSum);
    }

    // Intersection test: AABB vs AABB (More precise for long chains)
    intersectsBox(other) {
        if (this.max.x < other.min.x || this.min.x > other.max.x) return false;
        if (this.max.y < other.min.y || this.min.y > other.max.y) return false;
        return true;
    }
}
```

### 2. The Abstract Graph (Topology Layer)

This layer knows nothing about Chemistry. It deals with IDs, Indices, and Connectivity.

```javascript
class GraphNode {
    constructor(id, index) {
        this.id = id;          // Unique immutable ID
        this.index = index;    // Current index in the master array
        
        // Adjacency list: Stores indices of connected Edges.
        // We store Edges (not just Nodes) because we need to know *which* bond connects them.
        this.adj = []; 
    }
}

class GraphEdge {
    constructor(id, index, nodeIndexA, nodeIndexB) {
        this.id = id;
        this.index = index;
        this.a = nodeIndexA; // Index of Node A
        this.b = nodeIndexB; // Index of Node B
    }
    
    getOther(myNodeIndex) {
        return (this.a === myNodeIndex) ? this.b : this.a;
    }
}

class Graph {
    constructor() {
        this.nodes = []; // Dense Array
        this.edges = []; // Dense Array
        
        // Maps Unique ID -> Array Index for fast lookups
        this.nodeMap = new Map(); 
        this.edgeMap = new Map();
        
        this._idCounter = 1;
    }

    // --- Core Operations ---

    createNode() {
        const id = this._idCounter++;
        const index = this.nodes.length;
        const node = new GraphNode(id, index); // In subclass this will be Atom
        
        this.nodes.push(node);
        this.nodeMap.set(id, index);
        return node;
    }

    createEdge(nodeIndexA, nodeIndexB) {
        const id = this._idCounter++;
        const index = this.edges.length;
        const edge = new GraphEdge(id, index, nodeIndexA, nodeIndexB); // In subclass this will be Bond
        
        this.edges.push(edge);
        this.edgeMap.set(id, index);
        
        // Update topology
        this.nodes[nodeIndexA].adj.push(index);
        this.nodes[nodeIndexB].adj.push(index);
        
        return edge;
    }

    /**
     * Remove Node using Swap-and-Pop to keep array dense.
     * O(1) removal, but changes the index of the last element.
     */
    removeNode(indexToRemove) {
        // 1. Remove all attached edges first (omitted for brevity, but necessary)
        
        const lastIndex = this.nodes.length - 1;
        const nodeToRemove = this.nodes[indexToRemove];
        
        if (indexToRemove === lastIndex) {
            this.nodes.pop();
        } else {
            // Swap last element into the hole
            const nodeToMove = this.nodes[lastIndex];
            this.nodes[indexToRemove] = nodeToMove;
            nodeToMove.index = indexToRemove; // Update its internal index
            this.nodes.pop();
            
            // Update the Map
            this.nodeMap.set(nodeToMove.id, indexToRemove);
        }
        
        this.nodeMap.delete(nodeToRemove.id);
        
        // Note: Edges connected to 'nodeToMove' still point to the old index (lastIndex).
        // We must update those edges to point to 'indexToRemove'.
        // This is the cost of Swap-and-Pop in a graph.
        this._fixEdgeReferences(lastIndex, indexToRemove);
    }
    
    _fixEdgeReferences(oldIdx, newIdx) {
        const node = this.nodes[newIdx];
        for(let edgeIdx of node.adj) {
            const edge = this.edges[edgeIdx];
            if(edge.a === oldIdx) edge.a = newIdx;
            else if(edge.b === oldIdx) edge.b = newIdx;
        }
    }
}
```

### 3. The Molecular System (Chemical Layer)

This extends the Graph. `Atom` extends `GraphNode` and adds position/element. `Bond` extends `GraphEdge`.

```javascript
class Atom extends GraphNode {
    constructor(id, index, element, x, y) {
        super(id, index);
        this.element = element; // "C", "N", etc.
        this.pos = new Vec2(x, y);
        
        // Capping atoms are just a property (number of implicit Hydrogens)
        // Or a list of lightweight descriptors: [{element: "H"}, {element: "Cl"}]
        this.caps = []; 
        
        // Quick reference to high-level structure
        this.fragmentId = -1; // ID of the Ring/Chain this belongs to
    }
    
    addCap(element="H") {
        this.caps.push({ element });
    }
}

class Bond extends GraphEdge {
    constructor(id, index, u, v, order) {
        super(id, index, u, v);
        this.order = order; // 1, 2, 3
        this.isBridge = false; // Computed by algorithm
    }
}

/**
 * Represents a Ring or a Chain.
 * It is a "Group" of atoms with a Bounding Volume.
 */
class Fragment {
    constructor(id) {
        this.id = id;
        this.atomIndices = []; // List of integer indices
        this.isClosed = false; // True = Ring, False = Chain
        this.bounds = new Bounds();
    }

    // Recomputes bounds based on current atom positions
    updateBounds(system) {
        this.bounds.reset();
        for(let i of this.atomIndices) {
            const atom = system.nodes[i]; // Access via dense array
            this.bounds.addPoint(atom.pos);
        }
        this.bounds.recomputeSphere();
    }
}

class MolecularSystem extends Graph {
    constructor() {
        super();
        this.fragments = []; // Array of Fragment objects
    }

    // Override factory method
    createNode(element, x, y) {
        const id = this._idCounter++;
        const index = this.nodes.length;
        const atom = new Atom(id, index, element, x, y);
        
        this.nodes.push(atom);
        this.nodeMap.set(id, index);
        return atom;
    }

    // Override factory method
    createEdge(idxA, idxB, order=1) {
        const id = this._idCounter++;
        const index = this.edges.length;
        const bond = new Bond(id, index, idxA, idxB, order);
        
        this.edges.push(bond);
        this.edgeMap.set(id, index);
        
        this.nodes[idxA].adj.push(index);
        this.nodes[idxB].adj.push(index);
        return bond;
    }

    /**
     * Example Interaction: Find fragment under mouse.
     * Uses the optimized Bounding Sphere check first.
     */
    pickFragment(mousePos) {
        // 1. Broad Phase: Check Bounds
        const candidates = [];
        const mouseBounds = new Bounds(); 
        mouseBounds.center = mousePos; 
        mouseBounds.radius = 5.0; // Click tolerance

        for(let frag of this.fragments) {
            // Fast reject
            if(frag.bounds.intersectsSphere(mouseBounds)) {
                candidates.push(frag);
            }
        }

        // 2. Narrow Phase: Check Atoms inside candidates
        // ...
    }
}
```

### Why this fits your requirements:

1.  **Efficiency:** The core data is `this.nodes` (an Array). Loops for physics/drawing are `for(let i=0; i<n; i++) nodes[i]`. This is the fastest possible iteration in JS.
2.  **Scalability:** The `Bounds` class allows us to have 100 rings, and checking collisions only involves 100 simple math checks, not 100 * (atoms per ring) checks.
3.  **Topology:** By storing `adj` (adjacency) as indices on the `GraphNode`, implementing "Bridge Finding" or "Path Finding" is standard Graph Theory.
4.  **Flexibility:** `Fragment` handles both Rings and Chains uniformly. You just check `isClosed` if you specifically need ring logic (like aromaticity).
5.  **Capping:** Handled as `atom.caps`, keeping the main graph skeleton pure.

In the next steps (algorithms), you can write functions that take `MolecularSystem` as input, iterate the integer arrays, and populate `isBridge` flags or generate `fragments`.