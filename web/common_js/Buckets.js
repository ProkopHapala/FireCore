import { Vec3 } from './Vec3.js';

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

    addAtomIndex(i, pos = null) {
        this.atoms.push(i | 0);
        if (pos) {
            if (pos.x < this.pmin.x) this.pmin.x = pos.x;
            if (pos.y < this.pmin.y) this.pmin.y = pos.y;
            if (pos.z < this.pmin.z) this.pmin.z = pos.z;
            if (pos.x > this.pmax.x) this.pmax.x = pos.x;
            if (pos.y > this.pmax.y) this.pmax.y = pos.y;
            if (pos.z > this.pmax.z) this.pmax.z = pos.z;
        }
    }

    overlapAABB(other, margin = 0.0) {
        const m = (margin !== undefined) ? +margin : 0.0;
        const ax0 = this.pmin.x - m, ay0 = this.pmin.y - m, az0 = this.pmin.z - m;
        const ax1 = this.pmax.x + m, ay1 = this.pmax.y + m, az1 = this.pmax.z + m;
        const bx0 = other.pmin.x - m, by0 = other.pmin.y - m, bz0 = other.pmin.z - m;
        const bx1 = other.pmax.x + m, by1 = other.pmax.y + m, bz1 = other.pmax.z + m;
        return (ax0 <= bx1) && (ax1 >= bx0) && (ay0 <= by1) && (ay1 >= by0) && (az0 <= bz1) && (az1 >= bz0);
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

    recalcBounds(mol) {
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
                const p = a.pos;
                if (p.x < b.pmin.x) b.pmin.x = p.x;
                if (p.y < b.pmin.y) b.pmin.y = p.y;
                if (p.z < b.pmin.z) b.pmin.z = p.z;
                if (p.x > b.pmax.x) b.pmax.x = p.x;
                if (p.y > b.pmax.y) b.pmax.y = p.y;
                if (p.z > b.pmax.z) b.pmax.z = p.z;
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
}

export function crystalCellKey(ix, iy, iz, na, nb) { return ((iz * (nb | 0) + (iy | 0)) * (na | 0) + (ix | 0)) | 0; }

export function buildCrystalCellBucketsFromMol(mol, na, nb, nc, lvec, origin = new Vec3(0, 0, 0)) {
    if (!mol || !mol.atoms) throw new Error('buildCrystalCellBucketsFromMol: mol required');
    const A = lvec[0], B = lvec[1], C = lvec[2];
    const g = new BucketGraph();
    const ncell = (na | 0) * (nb | 0) * (nc | 0);
    const dense = new Int32Array(ncell);
    dense.fill(-1);

    for (let ia = 0; ia < mol.atoms.length; ia++) {
        const a = mol.atoms[ia];
        const m = a.cellIndex;
        if (m === undefined || m === null) continue;
        const ic = m | 0;
        if (ic < 0 || ic >= ncell) continue;
        let ib = dense[ic];
        if (ib < 0) {
            ib = g.buckets.length;
            dense[ic] = ib;
            const b = g.addBucket();
            const ix = ic % (na | 0);
            const iy = ((ic / (na | 0)) | 0) % (nb | 0);
            const iz = (ic / ((na | 0) * (nb | 0))) | 0;
            b.meta = { ix, iy, iz, icellDense: ic };
        }
        g.buckets[ib].addAtomIndex(ia, a.pos);
    }

    const bucketOfDense = dense;
    for (let ib = 0; ib < g.buckets.length; ib++) {
        const meta = g.buckets[ib].meta;
        const ix = meta.ix | 0, iy = meta.iy | 0, iz = meta.iz | 0;
        for (let dz = -1; dz <= 1; dz++) {
            for (let dy = -1; dy <= 1; dy++) {
                for (let dx = -1; dx <= 1; dx++) {
                    const jx = ix + dx, jy = iy + dy, jz = iz + dz;
                    if (jx < 0 || jy < 0 || jz < 0 || jx >= (na | 0) || jy >= (nb | 0) || jz >= (nc | 0)) continue;
                    const jc = crystalCellKey(jx, jy, jz, na, nb);
                    const jb = bucketOfDense[jc] | 0;
                    if (jb < 0) continue;
                    if (jb < ib) continue;
                    g.buckets[ib].neigh.push(jb);
                }
            }
        }
    }

    g.meta = { kind: 'crystal_cells', na: na | 0, nb: nb | 0, nc: nc | 0, origin: origin.clone(), lvec: [A.clone(), B.clone(), C.clone()], dense2bucket: bucketOfDense };
    return g;
}

export function buildWireframeCellVerts(A, B, C, O, out = null, i0 = 0) {
    const o = O ? O : new Vec3(0, 0, 0);
    const p100 = new Vec3().setV(o).add(A);
    const p010 = new Vec3().setV(o).add(B);
    const p001 = new Vec3().setV(o).add(C);
    const p110 = new Vec3().setV(o).add(A).add(B);
    const p101 = new Vec3().setV(o).add(A).add(C);
    const p011 = new Vec3().setV(o).add(B).add(C);
    const p111 = new Vec3().setV(o).add(A).add(B).add(C);
    const edges = [
        o, p100, o, p010, o, p001,
        p100, p110, p100, p101,
        p010, p110, p010, p011,
        p001, p101, p001, p011,
        p110, p111, p101, p111, p011, p111,
    ];
    const n = edges.length * 3;
    if (!out) out = new Float32Array(n + (i0 | 0));
    let k = i0 | 0;
    for (let i = 0; i < edges.length; i++) {
        const v = edges[i];
        out[k++] = v.x; out[k++] = v.y; out[k++] = v.z;
    }
    return out;
}

export function buildWireframeAABBVerts(pmin, pmax, out = null, i0 = 0) {
    const x0 = pmin.x, y0 = pmin.y, z0 = pmin.z;
    const x1 = pmax.x, y1 = pmax.y, z1 = pmax.z;
    const ps = [
        [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
        [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1],
    ];
    const e = [
        [0, 1], [1, 2], [2, 3], [3, 0],
        [4, 5], [5, 6], [6, 7], [7, 4],
        [0, 4], [1, 5], [2, 6], [3, 7],
    ];
    const n = e.length * 2 * 3;
    if (!out) out = new Float32Array(n + (i0 | 0));
    let k = i0 | 0;
    for (let i = 0; i < e.length; i++) {
        const a = ps[e[i][0]];
        const b = ps[e[i][1]];
        out[k++] = a[0]; out[k++] = a[1]; out[k++] = a[2];
        out[k++] = b[0]; out[k++] = b[1]; out[k++] = b[2];
    }
    return out;
}

export function _testBucketGraphBasic() {
    const g = new BucketGraph();
    const b0 = g.addBucket();
    b0.addAtomIndex(0, new Vec3(0, 0, 0));
    b0.addAtomIndex(1, new Vec3(1, 1, 1));
    const b1 = g.addBucket();
    b1.addAtomIndex(2, new Vec3(2, 2, 2));
    b0.neigh.push(1);
    if (!(b0.atoms.length === 2 && b1.atoms.length === 1)) throw new Error('_testBucketGraphBasic failed');
    return true;
}
