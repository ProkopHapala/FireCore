/// Uniform 3D grid partitioning for crystal unit cells.
/// Uses BucketGraph from Buckets.js as the generic data structure.
import { Vec3 } from './Vec3.js';
import { BucketGraph } from './Buckets.js';

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
