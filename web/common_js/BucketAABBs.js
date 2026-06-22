import { Vec3 } from './Vec3.js';

export class AABB {
    constructor(pmin = null, pmax = null) {
        this.pmin = pmin ? pmin : new Vec3(+Infinity, +Infinity, +Infinity);
        this.pmax = pmax ? pmax : new Vec3(-Infinity, -Infinity, -Infinity);
    }

    setEmpty() {
        this.pmin.x = +Infinity; this.pmin.y = +Infinity; this.pmin.z = +Infinity;
        this.pmax.x = -Infinity; this.pmax.y = -Infinity; this.pmax.z = -Infinity;
        return this;
    }

    setAABB(pmin, pmax) {
        this.pmin.x = +pmin.x; this.pmin.y = +pmin.y; this.pmin.z = +pmin.z;
        this.pmax.x = +pmax.x; this.pmax.y = +pmax.y; this.pmax.z = +pmax.z;
        return this;
    }

    overlap(other, margin = 0.0) {
        const m = (margin !== undefined) ? +margin : 0.0;
        const ax0 = this.pmin.x - m, ay0 = this.pmin.y - m, az0 = this.pmin.z - m;
        const ax1 = this.pmax.x + m, ay1 = this.pmax.y + m, az1 = this.pmax.z + m;
        const bx0 = other.pmin.x - m, by0 = other.pmin.y - m, bz0 = other.pmin.z - m;
        const bx1 = other.pmax.x + m, by1 = other.pmax.y + m, bz1 = other.pmax.z + m;
        return (ax0 <= bx1) && (ax1 >= bx0) && (ay0 <= by1) && (ay1 >= by0) && (az0 <= bz1) && (az1 >= bz0);
    }
}

export class AABBs {
    constructor(ncell = 0) {
        this.ncell = 0;
        this.pmins = [];
        this.pmaxs = [];
        if (ncell > 0) this.resize(ncell);
    }

    resize(ncell) {
        ncell = ncell | 0;
        if (!(ncell >= 0)) throw new Error('AABBs.resize: ncell must be >=0');
        if (ncell === (this.ncell | 0)) return false;
        const n0 = this.ncell | 0;
        if (ncell > n0) {
            for (let i = n0; i < ncell; i++) {
                this.pmins[i] = new Vec3(+Infinity, +Infinity, +Infinity);
                this.pmaxs[i] = new Vec3(-Infinity, -Infinity, -Infinity);
            }
        } else {
            this.pmins.length = ncell;
            this.pmaxs.length = ncell;
        }
        this.ncell = ncell;
        return true;
    }

    setEmpty(icell) {
        const pmin = this.pmins[icell | 0];
        const pmax = this.pmaxs[icell | 0];
        pmin.x = +Infinity; pmin.y = +Infinity; pmin.z = +Infinity;
        pmax.x = -Infinity; pmax.y = -Infinity; pmax.z = -Infinity;
    }

    setAABB(icell, pmin, pmax) {
        const i = icell | 0;
        const a = this.pmins[i];
        const b = this.pmaxs[i];
        a.x = +pmin.x; a.y = +pmin.y; a.z = +pmin.z;
        b.x = +pmax.x; b.y = +pmax.y; b.z = +pmax.z;
    }

    getPMin(icell, out = null) {
        const p = this.pmins[icell | 0];
        if (!out) out = new Vec3();
        out.x = p.x; out.y = p.y; out.z = p.z;
        return out;
    }

    getPMax(icell, out = null) {
        const p = this.pmaxs[icell | 0];
        if (!out) out = new Vec3();
        out.x = p.x; out.y = p.y; out.z = p.z;
        return out;
    }

    overlap(i, j, margin = 0.0) {
        const m = (margin !== undefined) ? +margin : 0.0;
        const a0 = this.pmins[i | 0], a1 = this.pmaxs[i | 0];
        const b0 = this.pmins[j | 0], b1 = this.pmaxs[j | 0];
        const ax0 = a0.x - m, ay0 = a0.y - m, az0 = a0.z - m;
        const ax1 = a1.x + m, ay1 = a1.y + m, az1 = a1.z + m;
        const bx0 = b0.x - m, by0 = b0.y - m, bz0 = b0.z - m;
        const bx1 = b1.x + m, by1 = b1.y + m, bz1 = b1.z + m;
        return (ax0 <= bx1) && (ax1 >= bx0) && (ay0 <= by1) && (ay1 >= by0) && (az0 <= bz1) && (az1 >= bz0);
    }
}

export class BucketAABBs extends AABBs {}

export function forEachAABBOverlappingPair(aabbs, fn, { margin = 0.0 } = {}) {
    const n = aabbs.ncell | 0;
    const m = (margin !== undefined) ? +margin : 0.0;
    const pmins = aabbs.pmins;
    const pmaxs = aabbs.pmaxs;
    for (let i = 0; i < n; i++) {
        const a0 = pmins[i], a1 = pmaxs[i];
        const ax0 = a0.x - m, ay0 = a0.y - m, az0 = a0.z - m;
        const ax1 = a1.x + m, ay1 = a1.y + m, az1 = a1.z + m;
        for (let j = i + 1; j < n; j++) {
            const b0 = pmins[j], b1 = pmaxs[j];
            const bx0 = b0.x - m, by0 = b0.y - m, bz0 = b0.z - m;
            const bx1 = b1.x + m, by1 = b1.y + m, bz1 = b1.z + m;
            const ok = (ax0 <= bx1) && (ax1 >= bx0) && (ay0 <= by1) && (ay1 >= by0) && (az0 <= bz1) && (az1 >= bz0);
            if (ok) fn(i, j);
        }
    }
}

export function _testBucketAABBsBasic() {
    const a = new BucketAABBs(2);
    a.setAABB(0, new Vec3(0, 0, 0), new Vec3(1, 1, 1));
    a.setAABB(1, new Vec3(0.5, 0.5, 0.5), new Vec3(2, 2, 2));
    let n = 0;
    forEachAABBOverlappingPair(a, () => { n++; });
    if (n !== 1) throw new Error('_testBucketAABBsBasic failed');
    return true;
}

export { Vec3 };
