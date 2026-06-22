import { Vec3 } from './Vec3.js';

export class Buckets {
    constructor(ncell = 0, nobj = 0, { withObj2cell = true } = {}) {
        this.maxInBucket = 0;
        this.nobjSize = 0;
        this.ncell    = 0;
        this.nobj     = 0;
        this.cellNs   = new Int32Array(0);
        this.cellI0s  = new Int32Array(0);
        this.cell2obj = new Int32Array(0);
        this.obj2cell = withObj2cell ? new Int32Array(0) : null;
        if (ncell > 0) this.resizeCells(ncell);
        if (nobj > 0) this.resizeObjs(nobj, { withObj2cell });
    }

    resizeCells(ncell) {
        ncell = ncell | 0;
        if (!(ncell >= 0)) throw new Error('Buckets.resizeCells: ncell must be >=0');
        if (ncell === this.ncell) return false;
        const cellNs = new Int32Array(ncell);
        const cellI0s = new Int32Array(ncell);
        const n0 = Math.min(this.ncell | 0, ncell);
        if (n0 > 0) {
            cellNs.set(this.cellNs.subarray(0, n0));
            cellI0s.set(this.cellI0s.subarray(0, n0));
        }
        this.cellNs = cellNs;
        this.cellI0s = cellI0s;
        this.ncell = ncell;
        return true;
    }

    resizeObjs(nobj, { withObj2cell = false } = {}) {
        nobj = nobj | 0;
        if (!(nobj >= 0)) throw new Error('Buckets.resizeObjs: nobj must be >=0');
        this.nobj = nobj;
        const need = nobj > (this.nobjSize | 0);
        if (need) {
            this.nobjSize = nobj;
            this.cell2obj = new Int32Array(this.nobjSize);
            this.cell2obj.fill(-1);
            if (withObj2cell) {
                this.obj2cell = new Int32Array(nobj);
                this.obj2cell.fill(-1);
            }
        } else {
            if (withObj2cell && (!this.obj2cell || (this.obj2cell.length | 0) < nobj)) {
                this.obj2cell = new Int32Array(nobj);
                this.obj2cell.fill(-1);
            }
        }
        return need;
    }

    bindObjs(nobj, obj2cell) {
        if (!obj2cell) throw new Error('Buckets.bindObjs: obj2cell required');
        this.resizeObjs(nobj, { withObj2cell: false });
        this.obj2cell = obj2cell;
    }

    clean() {
        this.cellNs.fill(0);
    }

    cleanO2C(fill = -1) {
        if (!this.obj2cell) throw new Error('Buckets.cleanO2C: obj2cell not allocated');
        this.obj2cell.fill(fill | 0);
    }

    count(nobj = -1, obj2cell = null) {
        if (nobj < 0) nobj = this.nobj | 0;
        if (!obj2cell) obj2cell = this.obj2cell;
        if (!obj2cell) throw new Error('Buckets.count: obj2cell required');
        if ((nobj | 0) > (obj2cell.length | 0)) throw new Error('Buckets.count: nobj exceeds obj2cell length');
        const ncell = this.ncell | 0;
        for (let i = 0; i < (nobj | 0); i++) {
            const ic = obj2cell[i] | 0;
            if (ic >= 0) {
                if (ic >= ncell) throw new Error(`Buckets.count: obj2cell[${i}]=${ic} out of range ncell=${ncell}`);
                this.cellNs[ic]++;
            }
        }
    }

    updateOffsets() {
        const ncell = this.ncell | 0;
        let ntot = 0;
        let maxInBucket = 0;
        for (let k = 0; k < ncell; k++) {
            this.cellI0s[k] = ntot;
            const ni = this.cellNs[k] | 0;
            ntot += ni;
            if (ni > maxInBucket) maxInBucket = ni;
            this.cellNs[k] = 0;
        }
        this.maxInBucket = maxInBucket;
        if (ntot > (this.nobjSize | 0)) throw new Error(`Buckets.updateOffsets: total items ntot=${ntot} exceeds nobjSize=${this.nobjSize}`);
        return ntot;
    }

    objectsToCells(nobj = -1, obj2cell = null) {
        if (nobj < 0) nobj = this.nobj | 0;
        if (!obj2cell) obj2cell = this.obj2cell;
        if (!obj2cell) throw new Error('Buckets.objectsToCells: obj2cell required');
        const ncell = this.ncell | 0;
        const nobjSize = this.nobjSize | 0;
        for (let i = 0; i < (nobj | 0); i++) {
            const ic = obj2cell[i] | 0;
            if (ic < 0) continue;
            if (ic >= ncell) throw new Error(`Buckets.objectsToCells: obj2cell[${i}]=${ic} out of range ncell=${ncell}`);
            const ni = this.cellNs[ic] | 0;
            const j = (this.cellI0s[ic] | 0) + ni;
            if (j >= nobjSize) throw new Error(`Buckets.objectsToCells: j=${j} out of range nobjSize=${nobjSize}`);
            this.cell2obj[j] = i;
            this.cellNs[ic] = ni + 1;
        }
    }

    updateCells(nobj = -1, obj2cell = null) {
        if (obj2cell === null) obj2cell = this.obj2cell;
        if (!obj2cell) throw new Error('Buckets.updateCells: obj2cell required');
        if (nobj < 0) nobj = this.nobj | 0;
        for (let i = 0; i < (this.nobjSize | 0); i++) this.cell2obj[i] = -1;
        this.clean();
        this.count(nobj, obj2cell);
        this.updateOffsets();
        this.objectsToCells(nobj, obj2cell);
    }

    inCellRange(icell) {
        icell = icell | 0;
        if (icell < 0 || icell >= (this.ncell | 0)) return { i0: 0, n: 0 };
        const i0 = this.cellI0s[icell] | 0;
        const n = this.cellNs[icell] | 0;
        return { i0, n };
    }

    inCell(icell) {
        const r = this.inCellRange(icell);
        return this.cell2obj.subarray(r.i0, r.i0 + r.n);
    }

    forEachInCell(icell, fn) {
        const r = this.inCellRange(icell);
        const i0 = r.i0 | 0;
        const i1 = i0 + (r.n | 0);
        for (let i = i0; i < i1; i++) fn(this.cell2obj[i] | 0);
    }
}

export function overlapAABB6(a6, i, b6, j, margin = 0.0) {
    const mi = (margin !== undefined) ? +margin : 0.0;
    const ia = (i * 6) | 0;
    const ib = (j * 6) | 0;
    const ax0 = a6[ia  ] - mi, ay0 = a6[ia+1] - mi, az0 = a6[ia+2] - mi;
    const ax1 = a6[ia+3] + mi, ay1 = a6[ia+4] + mi, az1 = a6[ia+5] + mi;
    const bx0 = b6[ib  ] - mi, by0 = b6[ib+1] - mi, bz0 = b6[ib+2] - mi;
    const bx1 = b6[ib+3] + mi, by1 = b6[ib+4] + mi, bz1 = b6[ib+5] + mi;
    return (ax0 <= bx1) && (ax1 >= bx0) && (ay0 <= by1) && (ay1 >= by0) && (az0 <= bz1) && (az1 >= bz0);
}

export function fitAABB6FromPoints(pos, indices, out6, ibox = 0) {
    const i6 = (ibox * 6) | 0;
    let x0 = +Infinity, y0 = +Infinity, z0 = +Infinity;
    let x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
    for (let ii = 0; ii < (indices.length | 0); ii++) {
        const i = indices[ii] | 0;
        const p = pos[i];
        const x = p.x, y = p.y, z = p.z;
        if (x < x0) x0 = x;
        if (y < y0) y0 = y;
        if (z < z0) z0 = z;
        if (x > x1) x1 = x;
        if (y > y1) y1 = y;
        if (z > z1) z1 = z;
    }
    out6[i6] = x0; out6[i6 + 1] = y0; out6[i6 + 2] = z0;
    out6[i6 + 3] = x1; out6[i6 + 4] = y1; out6[i6 + 5] = z1;
    return out6;
}

export function _testBucketsBasic() {
    const b = new Buckets(4, 10, { withObj2cell: true });
    for (let i = 0; i < 10; i++) b.obj2cell[i] = (i % 4);
    b.updateCells();
    let s = 0;
    for (let ic = 0; ic < 4; ic++) {
        const a = b.inCell(ic);
        for (let j = 0; j < a.length; j++) s += a[j] + 1;
    }
    if (s <= 0) throw new Error('_testBucketsBasic failed');
    return true;
}

export { Vec3 };
