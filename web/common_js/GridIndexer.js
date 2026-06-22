import { Vec3 } from './Vec3.js';

export class DenseGridIndexer {
    constructor(p0 = new Vec3(0, 0, 0), h = 1.0, nx = 1, ny = 1, nz = 1) {
        this.p0 = p0.clone ? p0.clone() : new Vec3(p0.x, p0.y, p0.z);
        this.h = +h;
        this.nx = nx | 0;
        this.ny = ny | 0;
        this.nz = nz | 0;
        if (!(this.h > 0)) throw new Error('DenseGridIndexer: h must be >0');
        if (this.nx <= 0 || this.ny <= 0 || this.nz <= 0) throw new Error('DenseGridIndexer: nx,ny,nz must be >0');
        this.ncell = (this.nx * this.ny * this.nz) | 0;
        this._strideY = this.nx | 0;
        this._strideZ = (this.nx * this.ny) | 0;
    }

    cellOfIJK(ix, iy, iz) {
        ix |= 0; iy |= 0; iz |= 0;
        if (ix < 0 || ix >= this.nx) return -1;
        if (iy < 0 || iy >= this.ny) return -1;
        if (iz < 0 || iz >= this.nz) return -1;
        return (ix + this._strideY * (iy + this.ny * iz)) | 0;
    }

    ijkOfCell(icell, out = null) {
        icell |= 0;
        if (!out) out = [0, 0, 0];
        if (icell < 0 || icell >= (this.ncell | 0)) { out[0] = 0; out[1] = 0; out[2] = 0; return out; }
        const iz = (icell / this._strideZ) | 0;
        const r0 = (icell - iz * this._strideZ) | 0;
        const iy = (r0 / this._strideY) | 0;
        const ix = (r0 - iy * this._strideY) | 0;
        out[0] = ix; out[1] = iy; out[2] = iz;
        return out;
    }

    cellOfPos(p) {
        const ix = Math.floor((p.x - this.p0.x) / this.h) | 0;
        const iy = Math.floor((p.y - this.p0.y) / this.h) | 0;
        const iz = Math.floor((p.z - this.p0.z) / this.h) | 0;
        return this.cellOfIJK(ix, iy, iz);
    }

    forEachNeighborCell(icell, fn) {
        const ijk = this.ijkOfCell(icell, _tmpIJK);
        const ix0 = ijk[0] | 0;
        const iy0 = ijk[1] | 0;
        const iz0 = ijk[2] | 0;
        for (let dz = -1; dz <= 1; dz++) {
            const iz = iz0 + dz;
            if (iz < 0 || iz >= this.nz) continue;
            for (let dy = -1; dy <= 1; dy++) {
                const iy = iy0 + dy;
                if (iy < 0 || iy >= this.ny) continue;
                for (let dx = -1; dx <= 1; dx++) {
                    const ix = ix0 + dx;
                    if (ix < 0 || ix >= this.nx) continue;
                    fn((ix + this._strideY * (iy + this.ny * iz)) | 0);
                }
            }
        }
    }

    fillObj2Cell(pos, obj2cell) {
        const n = pos.length | 0;
        if ((obj2cell.length | 0) < n) throw new Error('DenseGridIndexer.fillObj2Cell: obj2cell too small');
        for (let i = 0; i < n; i++) obj2cell[i] = this.cellOfPos(pos[i]);
        return obj2cell;
    }
}

function _packKey(ix, iy, iz) { return `${ix},${iy},${iz}`; }

export class SparseGridIndexer {
    constructor(p0 = new Vec3(0, 0, 0), h = 1.0) {
        this.p0 = p0.clone ? p0.clone() : new Vec3(p0.x, p0.y, p0.z);
        this.h = +h;
        if (!(this.h > 0)) throw new Error('SparseGridIndexer: h must be >0');
        this.map = new Map();
        this.cellX = new Int32Array(0);
        this.cellY = new Int32Array(0);
        this.cellZ = new Int32Array(0);
        this.ncell = 0;
    }

    clear() {
        this.map.clear();
        this.ncell = 0;
    }

    _ensureCells(n) {
        const n0 = this.cellX.length | 0;
        if (n <= n0) return;
        const n1 = Math.max(n, Math.max(16, (n0 * 2) | 0)) | 0;
        const x = new Int32Array(n1);
        const y = new Int32Array(n1);
        const z = new Int32Array(n1);
        if (n0 > 0) { x.set(this.cellX); y.set(this.cellY); z.set(this.cellZ); }
        this.cellX = x; this.cellY = y; this.cellZ = z;
    }

    ijkOfPos(p, out = null) {
        if (!out) out = [0, 0, 0];
        out[0] = Math.floor((p.x - this.p0.x) / this.h) | 0;
        out[1] = Math.floor((p.y - this.p0.y) / this.h) | 0;
        out[2] = Math.floor((p.z - this.p0.z) / this.h) | 0;
        return out;
    }

    cellOfIJK(ix, iy, iz) {
        const key = _packKey(ix | 0, iy | 0, iz | 0);
        const v = this.map.get(key);
        return (v === undefined) ? -1 : (v | 0);
    }

    cellOfPos(p) {
        const ijk = this.ijkOfPos(p, _tmpIJK);
        return this.cellOfIJK(ijk[0], ijk[1], ijk[2]);
    }

    buildFromPositions(pos, obj2cell) {
        this.clear();
        const n = pos.length | 0;
        if ((obj2cell.length | 0) < n) throw new Error('SparseGridIndexer.buildFromPositions: obj2cell too small');
        for (let i = 0; i < n; i++) {
            const ijk = this.ijkOfPos(pos[i], _tmpIJK);
            const ix = ijk[0] | 0, iy = ijk[1] | 0, iz = ijk[2] | 0;
            const key = _packKey(ix, iy, iz);
            let ic = this.map.get(key);
            if (ic === undefined) {
                ic = this.ncell | 0;
                this._ensureCells((ic + 1) | 0);
                this.cellX[ic] = ix;
                this.cellY[ic] = iy;
                this.cellZ[ic] = iz;
                this.map.set(key, ic);
                this.ncell = (ic + 1) | 0;
            }
            obj2cell[i] = ic | 0;
        }
        return obj2cell;
    }

    forEachNeighborCell(icell, fn) {
        icell |= 0;
        if (icell < 0 || icell >= (this.ncell | 0)) return;
        const ix0 = this.cellX[icell] | 0;
        const iy0 = this.cellY[icell] | 0;
        const iz0 = this.cellZ[icell] | 0;
        for (let dz = -1; dz <= 1; dz++) {
            const iz = (iz0 + dz) | 0;
            for (let dy = -1; dy <= 1; dy++) {
                const iy = (iy0 + dy) | 0;
                for (let dx = -1; dx <= 1; dx++) {
                    const ix = (ix0 + dx) | 0;
                    const key = _packKey(ix, iy, iz);
                    const jc = this.map.get(key);
                    if (jc !== undefined) fn(jc | 0);
                }
            }
        }
    }
}

export function chooseGridIndexer(p0, h, nx, ny, nz, nobj, { occMin = 0.1, maxDenseCells = 200000 } = {}) {
    const Nd = (nx * ny * nz) | 0;
    const occ = (Nd > 0) ? (nobj / Nd) : 0.0;
    if (Nd <= (maxDenseCells | 0) && occ >= +occMin) return new DenseGridIndexer(p0, h, nx, ny, nz);
    return new SparseGridIndexer(p0, h);
}

export function gridDimsFromAABB(pmin, pmax, h) {
    const hx = +h;
    if (!(hx > 0)) throw new Error('gridDimsFromAABB: h must be >0');
    const nx = Math.max(1, Math.ceil((pmax.x - pmin.x) / hx) | 0);
    const ny = Math.max(1, Math.ceil((pmax.y - pmin.y) / hx) | 0);
    const nz = Math.max(1, Math.ceil((pmax.z - pmin.z) / hx) | 0);
    return [nx, ny, nz];
}

export function fitAABBFromVec3(pos) {
    const n = pos.length | 0;
    const pmin = new Vec3(+Infinity, +Infinity, +Infinity);
    const pmax = new Vec3(-Infinity, -Infinity, -Infinity);
    for (let i = 0; i < n; i++) {
        const p = pos[i];
        if (p.x < pmin.x) pmin.x = p.x;
        if (p.y < pmin.y) pmin.y = p.y;
        if (p.z < pmin.z) pmin.z = p.z;
        if (p.x > pmax.x) pmax.x = p.x;
        if (p.y > pmax.y) pmax.y = p.y;
        if (p.z > pmax.z) pmax.z = p.z;
    }
    return { pmin, pmax };
}

export function _testGridIndexerBasic() {
    const pos = [new Vec3(0.1, 0.1, 0.1), new Vec3(1.1, 0.1, 0.1), new Vec3(0.1, 1.1, 0.1), new Vec3(0.1, 0.1, 1.1)];
    const obj2cell = new Int32Array(pos.length);
    const dense = new DenseGridIndexer(new Vec3(0, 0, 0), 1.0, 4, 4, 4);
    dense.fillObj2Cell(pos, obj2cell);
    let c = 0;
    dense.forEachNeighborCell(obj2cell[0], () => { c++; });
    if (!(c > 0)) throw new Error('_testGridIndexerBasic dense failed');

    const sparse = new SparseGridIndexer(new Vec3(0, 0, 0), 1.0);
    sparse.buildFromPositions(pos, obj2cell);
    c = 0;
    sparse.forEachNeighborCell(obj2cell[0], () => { c++; });
    if (!(c > 0)) throw new Error('_testGridIndexerBasic sparse failed');
    return true;
}

const _tmpIJK = [0, 0, 0];

export { Vec3 };
