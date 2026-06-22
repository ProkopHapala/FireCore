/// Build simple neighbor-index adjacency (no l0/K) from mol.bonds. Each entry is array of neighbor indices.
export function buildBondAdjacencySimple(nAtoms, bonds, mol) {
    const adj = new Array(nAtoms);
    const counts = new Int32Array(nAtoms);
    for (let i = 0; i < nAtoms; i++) adj[i] = [];
    for (let k = 0; k < bonds.length; k++) {
        const b = bonds[k];
        if (typeof b.ensureIndices === 'function') b.ensureIndices(mol);
        const a = b.a | 0, c = b.b | 0;
        if (a < 0 || a >= nAtoms || c < 0 || c >= nAtoms) throw new Error(`buildBondAdjacencySimple: bond index out of range (${a},${c}) nAtoms=${nAtoms}`);
        adj[a].push(c);
        adj[c].push(a);
        counts[a]++; counts[c]++;
    }
    return { adj, counts };
}

/// Canonical key for a ring: rotate so smallest atom is first, pick lexicographically smaller of forward/reversed.
function canonicalRingKey(ring) {
    const n = ring.length;
    let minIdx = 0, minVal = ring[0];
    for (let i = 1; i < n; i++) { if (ring[i] < minVal) { minVal = ring[i]; minIdx = i; } }
    const fwd = [], rev = [];
    for (let i = 0; i < n; i++) { fwd.push(ring[(minIdx + i) % n]); rev.push(ring[(minIdx - i + n) % n]); }
    let key = '';
    for (let i = 0; i < n; i++) {
        if (fwd[i] < rev[i]) { key = fwd.join(','); break; }
        if (fwd[i] > rev[i]) { key = rev.join(','); break; }
    }
    if (!key) key = fwd.join(',');
    return key;
}

/// Bounded BFS finding all shortest paths from src to dst, avoiding forbiddenVertex.
/// Returns array of paths (each is array of atom indices from src to dst inclusive).
function boundedShortestPaths(adj, src, dst, forbiddenVertex, maxDepth) {
    if (src === forbiddenVertex || dst === forbiddenVertex) return [];
    if (src === dst) return [[src]];
    const n = adj.length;
    const parent = new Int32Array(n).fill(-1);
    const depthArr = new Int32Array(n).fill(-1);
    const visited = new Uint8Array(n);
    visited[forbiddenVertex] = 1;
    visited[src] = 1;
    depthArr[src] = 0;
    const q = [src];
    let qi = 0;
    let foundDepth = -1;
    const hits = [];
    while (qi < q.length) {
        const u = q[qi++];
        const du = depthArr[u];
        if (foundDepth >= 0 && du >= foundDepth) continue;
        if (du >= maxDepth) continue;
        const neighs = adj[u];
        for (let k = 0; k < neighs.length; k++) {
            const v = neighs[k] | 0;
            if (visited[v]) continue;
            visited[v] = 1;
            parent[v] = u;
            depthArr[v] = du + 1;
            if (v === dst) {
                foundDepth = du + 1;
                const path = [];
                let p = v;
                while (p >= 0) { path.push(p); p = parent[p]; }
                path.reverse();
                hits.push(path);
            } else {
                q.push(v);
            }
        }
    }
    return hits;
}

/// Check if a ring is primitive: no pair of non-adjacent ring vertices has a graph shortcut.
function isPrimitiveRing(adj, ring) {
    const n = ring.length;
    if (n < 4) return true; // 3-rings are always primitive
    const ringSet = new Set(ring);
    for (let i = 0; i < n; i++) {
        for (let j = i + 2; j < n; j++) {
            if (i === 0 && j === n - 1) continue; // adjacent on cycle
            const dcycle = Math.min(j - i, n - (j - i));
            if (dcycle <= 1) continue;
            // BFS from ring[i] to ring[j], bounded by dcycle-1, avoiding ring vertices (except endpoints)
            const forbidden = new Set();
            for (let k = 0; k < n; k++) { if (k !== i && k !== j) forbidden.add(ring[k]); }
            const dist = bfsDistanceBounded(adj, ring[i], ring[j], forbidden, dcycle - 1);
            if (dist >= 0 && dist < dcycle) return false;
        }
    }
    return true;
}

function bfsDistanceBounded(adj, src, dst, forbiddenSet, maxDepth) {
    if (src === dst) return 0;
    const n = adj.length;
    const visited = new Uint8Array(n);
    if (forbiddenSet) for (const v of forbiddenSet) visited[v] = 1;
    visited[src] = 1;
    const q = [src];
    let qi = 0;
    let depth = 0;
    while (qi < q.length) {
        const levelSize = q.length - qi;
        if (depth >= maxDepth) return -1;
        for (let l = 0; l < levelSize; l++) {
            const u = q[qi++];
            const neighs = adj[u];
            for (let k = 0; k < neighs.length; k++) {
                const v = neighs[k] | 0;
                if (visited[v]) continue;
                if (v === dst) return depth + 1;
                visited[v] = 1;
                q.push(v);
            }
        }
        depth++;
    }
    return -1;
}

/// Find primitive rings using King/Franzblau method: for each atom c, for each pair of bonded neighbors (a,b),
/// find shortest path from a to b avoiding c. Ring = [c, a, ...path..., b, c].
/// @param adj — neighbor-index adjacency array (from buildBondAdjacencySimple or bondsAdj1 mapped)
/// @param opts — { kMin=3, kMax=8, lMax=2.5, pos=null, checkPrimitive=true }
/// @returns { rings, hist, ringOfAtom, ringOfBond }
export function findPrimitiveRings(adj, opts = {}) {
    const n = adj.length;
    const kMin = (opts.kMin !== undefined) ? (opts.kMin | 0) : 3;
    const kMax = (opts.kMax !== undefined) ? (opts.kMax | 0) : 8;
    const lMax = (opts.lMax !== undefined) ? +opts.lMax : 2.5;
    const checkPrim = (opts.checkPrimitive !== undefined) ? !!opts.checkPrimitive : true;
    const pos = opts.pos || null; // Float64Array or array of [x,y,z] for distance pruning
    const lMax2 = lMax * lMax;

    const seen = new Set();
    const rings = [];

    for (let c = 0; c < n; c++) {
        const neighs = adj[c];
        if (neighs.length < 2) continue;
        for (let ia = 0; ia < neighs.length; ia++) {
            const a = neighs[ia] | 0;
            for (let ib = ia + 1; ib < neighs.length; ib++) {
                const b = neighs[ib] | 0;
                // Ring size = pathLen(a→b) + 2. We need pathLen <= kMax - 2.
                const maxPathLen = kMax - 2;
                const paths = boundedShortestPaths(adj, a, b, c, maxPathLen);
                for (let pi = 0; pi < paths.length; pi++) {
                    const path = paths[pi];
                    const ringSize = path.length + 1; // +1 for c
                    if (ringSize < kMin || ringSize > kMax) continue;
                    // Geometric distance prune: |pos[a] - pos[b]| must be <= (ringSize-1)*lMax
                    if (pos) {
                        const pa = pos[a], pb = pos[b];
                        const dx = pa[0] - pb[0], dy = pa[1] - pb[1], dz = pa[2] - pb[2];
                        const d2 = dx*dx + dy*dy + dz*dz;
                        const maxReach = (ringSize - 1) * lMax;
                        if (d2 > maxReach * maxReach) continue;
                    }
                    const ring = [c, ...path];
                    if (ring.length !== ringSize) throw new Error(`findPrimitiveRings: ring length ${ring.length} != expected ${ringSize}`);
                    const key = canonicalRingKey(ring);
                    if (seen.has(key)) continue;
                    if (checkPrim && !isPrimitiveRing(adj, ring)) continue;
                    seen.add(key);
                    rings.push(ring);
                }
            }
        }
    }

    // Build histograms and per-atom / per-bond ring membership
    const hist = new Array(kMax + 1).fill(0);
    const ringOfAtom = new Array(n).fill(null).map(() => []);
    const ringOfBond = new Map();
    for (let r = 0; r < rings.length; r++) {
        const ring = rings[r];
        const sz = ring.length;
        hist[sz]++;
        for (let i = 0; i < sz; i++) {
            ringOfAtom[ring[i]].push(r);
            const a = ring[i], b = ring[(i + 1) % sz];
            const bk = a < b ? `${a},${b}` : `${b},${a}`;
            if (!ringOfBond.has(bk)) ringOfBond.set(bk, []);
            ringOfBond.get(bk).push(r);
        }
    }
    return { rings, hist, ringOfAtom, ringOfBond, natoms: n };
}

/// Classify rings by size. Returns { n5, n6, n7, defectRings, summary }
export function classifyRings(rings) {
    const hist = {};
    const defectRings = [];
    for (let r = 0; r < rings.length; r++) {
        const sz = rings[r].length;
        hist[sz] = (hist[sz] || 0) + 1;
        if (sz !== 6) defectRings.push({ index: r, size: sz, atoms: rings[r] });
    }
    const summary = Object.entries(hist).sort((a, b) => +a[0] - +b[0])
        .map(([k, v]) => `${k}-ring: ${v}`).join(', ');
    return { hist, defectRings, summary, nRings: rings.length };
}

/// Build adjacency from mol.bonds (convenience wrapper around buildBondAdjacencySimple).
export function buildAdjFromMol(mol) {
    const n = mol.atoms.length;
    return buildBondAdjacencySimple(n, mol.bonds, mol);
}

/// Filter adjacency to only heavy atoms (Z > 1), remapping indices.
/// Returns { adj, pos, oldIdx, nHeavy } where oldIdx maps new→old indices.
export function filterAdjToHeavy(mol, adjAll) {
    const heavyMap = new Map();
    const heavyList = [];
    for (let i = 0; i < mol.atoms.length; i++) {
        if (mol.atoms[i].Z > 1) {
            heavyMap.set(i, heavyList.length);
            heavyList.push(i);
        }
    }
    const nHeavy = heavyList.length;
    const adj = new Array(nHeavy);
    for (let i = 0; i < nHeavy; i++) adj[i] = [];
    for (let i = 0; i < mol.atoms.length; i++) {
        if (!heavyMap.has(i)) continue;
        const ni = heavyMap.get(i);
        for (const j of adjAll[i]) {
            if (!heavyMap.has(j)) continue;
            adj[ni].push(heavyMap.get(j));
        }
    }
    const pos = heavyList.map(i => [mol.atoms[i].pos.x, mol.atoms[i].pos.y, mol.atoms[i].pos.z]);
    return { adj, pos, oldIdx: heavyList, nHeavy };
}

/// Remap rings from heavy-only indices to mol indices.
export function remapRings(rings, oldIdx) {
    return rings.map(ring => ring.map(i => oldIdx[i]));
}

/// Map ringOfBond keys from heavy-only indices back to mol indices.
export function remapRingOfBond(ringOfBond, oldIdx) {
    const remapped = new Map();
    for (const [key, ringIds] of ringOfBond.entries()) {
        const [a, b] = key.split(',').map(Number);
        const na = oldIdx[a], nb = oldIdx[b];
        const nkey = na < nb ? `${na},${nb}` : `${nb},${na}`;
        remapped.set(nkey, ringIds);
    }
    return remapped;
}

/// Complete ring detection pipeline on a molecule: build adj, filter to heavy, find primitive rings, classify, remap.
/// Returns { ringsMol, ringOfBondMol, cls, result, oldIdx }.
export function runRingsOnMol(mol, opts = {}) {
    const { adj: adjAll } = buildAdjFromMol(mol);
    const { adj, pos, oldIdx, nHeavy } = filterAdjToHeavy(mol, adjAll);
    const kMin = (opts.kMin !== undefined) ? (opts.kMin | 0) : 3;
    const kMax = (opts.kMax !== undefined) ? (opts.kMax | 0) : 8;
    const lMax = (opts.lMax !== undefined) ? +opts.lMax : 2.0;
    const checkPrim = opts.checkPrimitive !== false;
    const result = findPrimitiveRings(adj, { kMin, kMax, pos, lMax, checkPrimitive: checkPrim });
    const cls = classifyRings(result.rings);
    const ringsMol = remapRings(result.rings, oldIdx);
    const ringOfBondMol = remapRingOfBond(result.ringOfBond, oldIdx);
    return { ringsMol, ringOfBondMol, cls, result, oldIdx, nHeavy };
}