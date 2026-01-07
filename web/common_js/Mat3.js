import { Vec3 } from "./Vec3.js";

export class Mat3 {
    constructor(a = null, b = null, c = null) {
        this.a = a ? a.clone() : new Vec3(1, 0, 0);
        this.b = b ? b.clone() : new Vec3(0, 1, 0);
        this.c = c ? c.clone() : new Vec3(0, 0, 1);
    }

    setABC(a, b, c) { this.a.setV(a); this.b.setV(b); this.c.setV(c); return this; }

    clone() { return new Mat3(this.a, this.b, this.c); }

    mulVec(v, out = new Vec3()) {
        out.set(
            this.a.x * v.x + this.b.x * v.y + this.c.x * v.z,
            this.a.y * v.x + this.b.y * v.y + this.c.y * v.z,
            this.a.z * v.x + this.b.z * v.y + this.c.z * v.z
        );
        return out;
    }

    transpose() {
        const ax = this.a.x, ay = this.a.y, az = this.a.z;
        const bx = this.b.x, by = this.b.y, bz = this.b.z;
        const cx = this.c.x, cy = this.c.y, cz = this.c.z;
        this.a.set(ax, bx, cx);
        this.b.set(ay, by, cy);
        this.c.set(az, bz, cz);
        return this;
    }

    static mul(A, B) {
        const a = new Vec3(), b = new Vec3(), c = new Vec3();
        A.mulVec(B.a, a);
        A.mulVec(B.b, b);
        A.mulVec(B.c, c);
        return new Mat3(a, b, c);
    }

    static fromAxisAngle(axis, angle) {
        const u = axis.clone();
        const l = u.normalize();
        if (l < 1e-12) return new Mat3();
        const ca = Math.cos(angle);
        const sa = Math.sin(angle);

        const ex = new Vec3(1, 0, 0).rotateCSA(ca, sa, u);
        const ey = new Vec3(0, 1, 0).rotateCSA(ca, sa, u);
        const ez = new Vec3(0, 0, 1).rotateCSA(ca, sa, u);
        return new Mat3(ex, ey, ez);
    }

    static alignVectorToZ(v) {
        const n = v.clone();
        const l = n.normalize();
        if (l < 1e-12) throw new Error('Mat3.alignVectorToZ: zero vector');
        const z = new Vec3(0, 0, 1);
        const c = n.dot(z);
        if (c > 1.0 - 1e-12) return new Mat3();
        if (c < -1.0 + 1e-12) return new Mat3(new Vec3(-1, 0, 0), new Vec3(0, -1, 0), new Vec3(0, 0, 1));
        const axis = new Vec3().setCross(n, z);
        const s = axis.norm();
        if (s < 1e-12) return new Mat3();
        axis.mulScalar(1.0 / s);
        const ang = Math.atan2(s, c);
        return Mat3.fromAxisAngle(axis, ang);
    }

    /// Build an orthonormal frame from forward and up vectors.
    static fromForwardUp(forward, up) {
        const f = forward.clone().normalize();
        const u = up.clone().normalize();
        const s = u.clone().setCross(f, u).normalize();
        const v = f.clone().setCross(s, f).normalize();
        return new Mat3(s, v, f);
    }

    /// Apply rotation about the forward axis of an existing frame.
    static rotateAroundForward(M, angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        const R = new Mat3(
            new Vec3(c, -s, 0),
            new Vec3(s,  c, 0),
            new Vec3(0,  0, 1)
        );
        return Mat3.mul(M, R);
    }
}

if (typeof window !== 'undefined') {
    window.Mat3 = Mat3;
}
