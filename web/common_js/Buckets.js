/// Generic spatial partitioning data structure — buckets with AABB bounds and neighbor lists.
/// Not tied to any specific geometry or partitioning strategy.
import { Vec3 } from './Vec3.js';

/// Standalone AABB overlap test for flat array / typed array usage.
/// aMin/aMax/bMin/bMax can be [x,y,z] arrays or Vec3-like objects with numeric indices.
export function aabbOverlap3D(aMin, aMax, bMin, bMax, margin = 0.0) {
    return (aMax[0] + margin >= bMin[0] && aMin[0] <= bMax[0] + margin && aMax[1] + margin >= bMin[1] && aMin[1] <= bMax[1] + margin && aMax[2] + margin >= bMin[2] && aMin[2] <= bMax[2] + margin);
}

/// Squared distance from point p to AABB [bMin, bMax] (0 if inside).
export function dist2ToAabb(p, bMin, bMax) {
    const dx = Math.max(0.0, Math.max(bMin[0] - p[0], p[0] - bMax[0]));
    const dy = Math.max(0.0, Math.max(bMin[1] - p[1], p[1] - bMax[1]));
    const dz = Math.max(0.0, Math.max(bMin[2] - p[2], p[2] - bMax[2]));
    return dx * dx + dy * dy + dz * dz;
}

export class Bucket {
    constructor() {
        this.atoms = [];
        this.neigh = [];
        this.pmin = new Vec3(+Infinity, +Infinity, +Infinity);
        this.pmax = new Vec3(-Infinity, -Infinity, -Infinity);
        this.meta = null;
    }

    resetBounds() {
        this.pmin.x = +Infinity; this.pmin.y = +Infinity; this.pmin.z = +Infinity;
        this.pmax.x = -Infinity; this.pmax.y = -Infinity; this.pmax.z = -Infinity;
    }

    /// Add atom index i to bucket. If pos provided, updates AABB bounds.
    /// Optional radius expands bounds by r (for VdW-aware collision AABBs).
    addAtomIndex(i, pos = null, radius = 0.0) {
        this.atoms.push(i | 0);
        if (pos) {
            const r = +radius;
            if (pos.x - r < this.pmin.x) this.pmin.x = pos.x - r;
            if (pos.y - r < this.pmin.y) this.pmin.y = pos.y - r;
            if (pos.z - r < this.pmin.z) this.pmin.z = pos.z - r;
            if (pos.x + r > this.pmax.x) this.pmax.x = pos.x + r;
            if (pos.y + r > this.pmax.y) this.pmax.y = pos.y + r;
            if (pos.z + r > this.pmax.z) this.pmax.z = pos.z + r;
        }
    }

    overlapAABB(other, margin = 0.0) {
        return aabbOverlap3D(this.pmin, this.pmax, other.pmin, other.pmax, margin);
    }
}

export class BucketGraph {
    constructor() {
        this.buckets = [];
        this.meta = null;
        this.bStoreId = false;
    }

    clear() {
        this.buckets.length = 0;
        this.meta = null;
    }

    get n() { return this.buckets.length; }

    addBucket() {
        const b = new Bucket();
        this.buckets.push(b);
        return b;
    }

    toIds(mol) {
        if (this.bStoreId) return;
        if (!mol || !mol.atoms) throw new Error('BucketGraph.toIds: mol required');
        for (let ib = 0; ib < this.buckets.length; ib++) {
            const bs = this.buckets[ib].atoms;
            for (let k = 0; k < bs.length; k++) {
                const ia = bs[k] | 0;
                const a = mol.atoms[ia];
                if (!a) throw new Error(`BucketGraph.toIds: atom missing for index ia=${ia} in bucket ib=${ib}`);
                bs[k] = a.id;
            }
        }
        this.bStoreId = true;
    }

    toInds(mol) {
        if (!this.bStoreId) return;
        if (!mol || !mol.atoms || !mol.getAtomIndex) throw new Error('BucketGraph.toInds: mol required');
        for (let ib = 0; ib < this.buckets.length; ib++) {
            const bs = this.buckets[ib].atoms;
            let n = 0;
            for (let k = 0; k < bs.length; k++) {
                const id = bs[k];
                const ia = mol.getAtomIndex(id);
                if (ia < 0) continue;
                bs[n++] = ia;
            }
            bs.length = n;
        }
        this.bStoreId = false;
    }

    recalcBounds(mol, radius = null) {
        if (this.bStoreId) throw new Error('BucketGraph.recalcBounds: buckets are in ID mode; call toInds(mol) first');
        if (!mol || !mol.atoms) throw new Error('BucketGraph.recalcBounds: mol required');
        for (let ib = 0; ib < this.buckets.length; ib++) {
            const b = this.buckets[ib];
            b.resetBounds();
            const bs = b.atoms;
            for (let k = 0; k < bs.length; k++) {
                const ia = bs[k] | 0;
                const a = mol.atoms[ia];
                if (!a) throw new Error(`BucketGraph.recalcBounds: atom missing for index ia=${ia} in bucket ib=${ib}`);
                const r = radius ? +radius[ia] : 0.0;
                const p = a.pos;
                if (p.x - r < b.pmin.x) b.pmin.x = p.x - r;
                if (p.y - r < b.pmin.y) b.pmin.y = p.y - r;
                if (p.z - r < b.pmin.z) b.pmin.z = p.z - r;
                if (p.x + r > b.pmax.x) b.pmax.x = p.x + r;
                if (p.y + r > b.pmax.y) b.pmax.y = p.y + r;
                if (p.z + r > b.pmax.z) b.pmax.z = p.z + r;
            }
        }
    }

    pruneEmptyBuckets() {
        if (this.bStoreId) throw new Error('BucketGraph.pruneEmptyBuckets: buckets are in ID mode; call toInds(mol) first');
        const n0 = this.buckets.length | 0;
        const map = new Int32Array(n0);
        let n1 = 0;
        for (let ib = 0; ib < n0; ib++) {
            const b = this.buckets[ib];
            if (b && b.atoms && b.atoms.length) map[ib] = n1++;
            else map[ib] = -1;
        }
        if (n1 === n0) return map;

        const bsNew = new Array(n1);
        let j = 0;
        for (let ib = 0; ib < n0; ib++) {
            const jb = map[ib];
            if (jb < 0) continue;
            const b = this.buckets[ib];
            const neighNew = [];
            for (let k = 0; k < b.neigh.length; k++) {
                const kn = map[b.neigh[k] | 0] | 0;
                if (kn < 0) continue;
                neighNew.push(kn);
            }
            b.neigh = neighNew;
            bsNew[j++] = b;
        }
        this.buckets = bsNew;
        if (this.meta && this.meta.dense2bucket) this.meta.dense2bucket = null;
        return map;
    }

    getBucketCenterFromBounds(ib, out = null) {
        const b = this.buckets[ib | 0];
        if (!b) throw new Error(`BucketGraph.getBucketCenterFromBounds: bucket missing ib=${ib}`);
        if (!out) out = new Vec3();
        out.x = 0.5 * (b.pmin.x + b.pmax.x);
        out.y = 0.5 * (b.pmin.y + b.pmax.y);
        out.z = 0.5 * (b.pmin.z + b.pmax.z);
        return out;
    }

    /// Find overlapping buckets via AABB intersection and populate each bucket's neigh[] list.
    /// Optional margin expands overlap test. Symmetric: if A overlaps B, both get each other.
    findOverlapNeighbors(margin = 0.0) {
        const n = this.buckets.length;
        for (let i = 0; i < n; i++) this.buckets[i].neigh.length = 0;
        for (let i = 0; i < n; i++) {
            const bi = this.buckets[i];
            for (let j = i + 1; j < n; j++) {
                const bj = this.buckets[j];
                if (bi.overlapAABB(bj, margin)) {
                    bi.neigh.push(j);
                    bj.neigh.push(i);
                }
            }
        }
    }

    /// Export to flat typed arrays for GPU / NPZ packing.
    /// Returns { group_atoms, group_nAtoms, bbox_min, bbox_max, icol, icolGroup, nGroups, groupCap }
    /// where groupCap is max atoms per bucket. icol/icolGroup are per-atom (nAtoms total, -1 if not in any bucket).
    exportFlat(nAtoms, groupCap = 0) {
        const nGroups = this.buckets.length;
        if (groupCap <= 0) {
            for (let g = 0; g < nGroups; g++) {
                const c = this.buckets[g].atoms.length;
                if (c > groupCap) groupCap = c;
            }
        }
        const group_atoms = new Int32Array(nGroups * groupCap);
        group_atoms.fill(-1);
        const group_nAtoms = new Int32Array(nGroups);
        const bbox_min = new Float64Array(nGroups * 3);
        const bbox_max = new Float64Array(nGroups * 3);
        const icol = new Int32Array(nAtoms);
        const icolGroup = new Int32Array(nAtoms);
        icol.fill(-1);
        icolGroup.fill(-1);

        for (let g = 0; g < nGroups; g++) {
            const b = this.buckets[g];
            const atoms = b.atoms;
            const n = atoms.length;
            if (n > groupCap) throw new Error(`BucketGraph.exportFlat: group ${g} has ${n} atoms > groupCap ${groupCap}`);
            group_nAtoms[g] = n;
            for (let k = 0; k < n; k++) {
                const ia = atoms[k] | 0;
                group_atoms[g * groupCap + k] = ia;
                icol[ia] = (g * groupCap + k) | 0;
                icolGroup[ia] = g;
            }
            bbox_min[g * 3 + 0] = b.pmin.x; bbox_min[g * 3 + 1] = b.pmin.y; bbox_min[g * 3 + 2] = b.pmin.z;
            bbox_max[g * 3 + 0] = b.pmax.x; bbox_max[g * 3 + 1] = b.pmax.y; bbox_max[g * 3 + 2] = b.pmax.z;
        }

        return { group_atoms, group_nAtoms, bbox_min, bbox_max, icol, icolGroup, nGroups, groupCap: groupCap | 0 };
    }
}
