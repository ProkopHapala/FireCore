import { Vec3 } from "../../common_js/Vec3.js";
import { Mat3 } from "../../common_js/Mat3.js";

export class Atom {
    constructor(id, i, Z = 6, x = 0, y = 0, z = 0) {
        this.id = id;
        this.i = i;
        this.Z = Z;
        this.charge = 0;
        this.flags = 0;
        this.pos = new Vec3(x, y, z);
        this.bonds = [];
        this.frag = -1;
        this.fragSlot = -1;
    }
}

export class Bond {
    constructor(id, i, aId, bId, order = 1, type = 1) {
        this.id = id;
        this.i = i;
        this.aId = aId;
        this.bId = bId;
        this.a = -1;
        this.b = -1;
        this.order = order;
        this.type = type;
        this.isBridge = false;
        this.isRingEdge = false;
        this.topoVersionCached = -1;
    }

    other(ai) { return (ai === this.a) ? this.b : this.a; }

    ensureIndices(mol) {
        if (this.topoVersionCached === mol.topoVersion) return;
        const a = mol.getAtomIndex(this.aId);
        const b = mol.getAtomIndex(this.bId);
        if (a < 0 || b < 0) throw new Error(`Bond.ensureIndices(): atom missing for bond id=${this.id} aId=${this.aId} bId=${this.bId}`);
        this.a = a;
        this.b = b;
        this.topoVersionCached = mol.topoVersion;
    }
}

export class Bounds {
    constructor() {
        this.min = new Vec3(+Infinity, +Infinity, +Infinity);
        this.max = new Vec3(-Infinity, -Infinity, -Infinity);
        this.center = new Vec3(0, 0, 0);
        this.radius = 0;
    }

    reset() {
        this.min.set(+Infinity, +Infinity, +Infinity);
        this.max.set(-Infinity, -Infinity, -Infinity);
        this.center.set(0, 0, 0);
        this.radius = 0;
    }

    addPoint(p) {
        if (p.x < this.min.x) this.min.x = p.x;
        if (p.y < this.min.y) this.min.y = p.y;
        if (p.z < this.min.z) this.min.z = p.z;
        if (p.x > this.max.x) this.max.x = p.x;
        if (p.y > this.max.y) this.max.y = p.y;
        if (p.z > this.max.z) this.max.z = p.z;
    }

    finalize() {
        this.center.set(
            (this.min.x + this.max.x) * 0.5,
            (this.min.y + this.max.y) * 0.5,
            (this.min.z + this.max.z) * 0.5,
        );
        const dx = this.max.x - this.center.x;
        const dy = this.max.y - this.center.y;
        const dz = this.max.z - this.center.z;
        this.radius = Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    intersectsSphere(other) {
        const r = this.radius + other.radius;
        const dx = this.center.x - other.center.x;
        const dy = this.center.y - other.center.y;
        const dz = this.center.z - other.center.z;
        return (dx * dx + dy * dy + dz * dz) <= (r * r);
    }

    intersectsAABB(other) {
        if (this.max.x < other.min.x || this.min.x > other.max.x) return false;
        if (this.max.y < other.min.y || this.min.y > other.max.y) return false;
        if (this.max.z < other.min.z || this.min.z > other.max.z) return false;
        return true;
    }
}

export class Fragment {
    constructor(id, isClosed = false) {
        this.id = id;
        this.isClosed = isClosed;
        this.atomIds = [];
        this.bondIds = [];
        this.bounds = new Bounds();
    }

    updateBounds(mol) {
        this.bounds.reset();
        for (const id of this.atomIds) {
            const ia = mol.getAtomIndex(id);
            if (ia < 0) continue;
            this.bounds.addPoint(mol.atoms[ia].pos);
        }
        this.bounds.finalize();
    }
}

export class EditableMolecule {
    constructor() {
        this.atoms = [];
        this.bonds = [];
        this.fragments = [];

        this.id2atom = new Map();
        this.id2bond = new Map();

        this.lastAtomId = 0;
        this.lastBondId = 0;

        this.selection = new Set();

        this.dirtyTopo = true;
        this.dirtyGeom = true;
        this.dirtyFrags = true;
        this.dirtyExport = true;

        this.topoVersion = 0;
        this.lockDepth = 0;
        this._topoVersionLocked = -1;
    }

    get nAtoms() { return this.atoms.length; }

    _touchTopo() {
        this.topoVersion++;
        this.dirtyTopo = true;
        this.dirtyFrags = true;
        this.dirtyExport = true;
    }

    _touchGeom() {
        this.dirtyGeom = true;
        this.dirtyFrags = true;
        this.dirtyExport = true;
    }

    _assertUnlocked(op = "edit") {
        if (this.lockDepth > 0) throw new Error(`EditableMolecule: topology is locked, cannot ${op}`);
    }

    lockTopology() {
        this.lockDepth++;
        if (this.lockDepth === 1) this._topoVersionLocked = this.topoVersion;
        return this._topoVersionLocked;
    }

    unlockTopology() {
        if (this.lockDepth <= 0) throw new Error('EditableMolecule.unlockTopology(): not locked');
        this.lockDepth--;
        if (this.lockDepth === 0) this._topoVersionLocked = -1;
    }

    assertLocked() {
        if (this.lockDepth <= 0) throw new Error('EditableMolecule: expected topology lock');
        if (this._topoVersionLocked !== this.topoVersion) throw new Error('EditableMolecule: topoVersion changed while locked');
    }

    getAtomIndex(id) {
        const i = this.id2atom.get(id);
        return (i === undefined) ? -1 : i;
    }

    getBondIndex(id) {
        const i = this.id2bond.get(id);
        return (i === undefined) ? -1 : i;
    }

    addAtom(x = 0, y = 0, z = 0, Z = 6) { // legacy signature compatible with MoleculeSystem.addAtom(x,y,z,type)
        this._assertUnlocked('addAtom');
        const id = ++this.lastAtomId;
        const i = this.atoms.length;
        const a = new Atom(id, i, Z, x, y, z);
        this.atoms.push(a);
        this.id2atom.set(id, i);
        this._touchTopo();
        return id;
    }

    addAtomZ(Z = 6, x = 0, y = 0, z = 0) { return this.addAtom(x, y, z, Z); }

    setAtomPosById(id, x, y, z) {
        const i = this.getAtomIndex(id);
        if (i < 0) throw new Error(`setAtomPosById(): atom not found id=${id}`);
        const p = this.atoms[i].pos;
        p.x = x; p.y = y; p.z = z;
        this._touchGeom();
    }

    addBond(aId, bId, order = 1, type = 1) {
        this._assertUnlocked('addBond');
        const a = this.getAtomIndex(aId);
        const b = this.getAtomIndex(bId);
        if (a < 0 || b < 0) throw new Error(`addBond(): atom missing aId=${aId} bId=${bId}`);
        const id = ++this.lastBondId;
        const i = this.bonds.length;
        const bond = new Bond(id, i, aId, bId, order, type);
        bond.a = a;
        bond.b = b;
        this.bonds.push(bond);
        this.id2bond.set(id, i);
        this.atoms[a].bonds.push(i);
        this.atoms[b].bonds.push(i);
        this._touchTopo();
        bond.topoVersionCached = this.topoVersion;
        return id;
    }

    selectAtom(id, mode = 'replace') {
        if (mode === 'replace') {
            this.selection.clear();
            this.selection.add(id);
        } else if (mode === 'add') {
            this.selection.add(id);
        } else if (mode === 'subtract') {
            this.selection.delete(id);
        }
        this.dirtyExport = true;
    }

    select(idOrIndex, mode = 'replace') { // compatibility with MoleculeSystem.select()
        const id = this.id2atom.has(idOrIndex) ? idOrIndex : (this.atoms[idOrIndex] ? this.atoms[idOrIndex].id : -1);
        if (id < 0) throw new Error(`select(): invalid atom id/index=${idOrIndex}`);
        this.selectAtom(id, mode);
    }

    clearSelection() { this.selection.clear(); this.dirtyExport = true; }

    deleteSelectedAtoms() {
        this._assertUnlocked('deleteSelectedAtoms');
        if (this.selection.size === 0) return;
        const ids = Array.from(this.selection);
        this.selection.clear();
        for (const id of ids) this.removeAtomById(id);
        this.dirtyExport = true;
    }

    _removeBondIndexFromAtom(ai, ib) {
        const bs = this.atoms[ai].bonds;
        for (let k = 0; k < bs.length; k++) {
            if (bs[k] === ib) { bs[k] = bs[bs.length - 1]; bs.pop(); return; }
        }
    }

    removeBondByIndex(ib) {
        this._assertUnlocked('removeBond');
        const nb = this.bonds.length;
        if (ib < 0 || ib >= nb) throw new Error(`removeBondByIndex(): out of range ib=${ib}`);

        const bond = this.bonds[ib];
        bond.ensureIndices(this);

        this._removeBondIndexFromAtom(bond.a, ib);
        this._removeBondIndexFromAtom(bond.b, ib);

        const last = nb - 1;
        if (ib !== last) {
            const moved = this.bonds[last];
            this.bonds[ib] = moved;
            moved.i = ib;
            this.id2bond.set(moved.id, ib);

            moved.ensureIndices(this);
            this._replaceBondIndexInAtom(moved.a, last, ib);
            this._replaceBondIndexInAtom(moved.b, last, ib);
        }
        this.bonds.pop();
        this.id2bond.delete(bond.id);
        this._touchTopo();
    }

    _replaceBondIndexInAtom(ai, fromIb, toIb) {
        const bs = this.atoms[ai].bonds;
        for (let k = 0; k < bs.length; k++) {
            if (bs[k] === fromIb) { bs[k] = toIb; return; }
        }
    }

    removeBondById(id) {
        const ib = this.getBondIndex(id);
        if (ib < 0) return;
        this.removeBondByIndex(ib);
    }

    removeAtomByIndex(i) {
        this._assertUnlocked('removeAtom');
        const na = this.atoms.length;
        if (i < 0 || i >= na) throw new Error(`removeAtomByIndex(): out of range i=${i}`);

        const a = this.atoms[i];
        const bondIdsToRemove = a.bonds.map(ib => this.bonds[ib] ? this.bonds[ib].id : -1).filter(id => id >= 0);
        for (const id of bondIdsToRemove) this.removeBondById(id);

        const last = na - 1;
        if (i !== last) {
            const moved = this.atoms[last];
            this.atoms[i] = moved;
            moved.i = i;
            this.id2atom.set(moved.id, i);
        }
        this.atoms.pop();
        this.id2atom.delete(a.id);
        this.selection.delete(a.id);
        this._touchTopo();
    }

    removeAtomById(id) {
        const i = this.getAtomIndex(id);
        if (i < 0) return;
        this.removeAtomByIndex(i);
    }

    clear() {
        this._assertUnlocked('clear');
        this.atoms.length = 0;
        this.bonds.length = 0;
        this.fragments.length = 0;
        this.id2atom.clear();
        this.id2bond.clear();
        this.selection.clear();
        this._touchTopo();
    }

    addAtomsFromArrays(pos3, types1) {
        if (!pos3 || !types1) throw new Error('addAtomsFromArrays: pos3/types1 must be provided');
        const n = (types1.length | 0);
        if ((pos3.length | 0) !== n * 3) throw new Error(`addAtomsFromArrays: pos3.length=${pos3.length} != 3*types.length=${n * 3}`);
        const ids = new Array(n);
        for (let i = 0; i < n; i++) {
            const i3 = i * 3;
            ids[i] = this.addAtom(pos3[i3], pos3[i3 + 1], pos3[i3 + 2], types1[i]);
        }
        return ids;
    }

    updateNeighborList() { /* adjacency is maintained incrementally in atoms[].bonds */ }

    attachGroupByMarker() { throw new Error('EditableMolecule.attachGroupByMarker(): not implemented yet (needs port from legacy MoleculeSystem)'); }
    attachParsedByDirection() { throw new Error('EditableMolecule.attachParsedByDirection(): not implemented yet (needs port from legacy MoleculeSystem)'); }

    appendParsedSystem(other, opts = {}) {
        if (!other || !other.pos || !other.types) throw new Error('appendParsedSystem: other must have pos/types');
        const n = other.types.length | 0;
        const pos = (opts.pos !== undefined) ? opts.pos : new Vec3(0, 0, 0);
        const rot = (opts.rot !== undefined) ? opts.rot : null;
        const bPosVec = (pos instanceof Vec3);
        if (!bPosVec && !(Array.isArray(pos) && pos.length >= 3)) throw new Error('appendParsedSystem: opts.pos must be Vec3 or [x,y,z]');
        const px = bPosVec ? pos.x : pos[0];
        const py = bPosVec ? pos.y : pos[1];
        const pz = bPosVec ? pos.z : pos[2];

        const bRotMat3 = (rot instanceof Mat3);
        const bRotFlat = (!!rot && !bRotMat3 && (rot.length === 9));
        if (rot && !bRotMat3 && !bRotFlat) throw new Error('appendParsedSystem: opts.rot must be Mat3 or flat[9]');
        const tmp = new Vec3();
        const tmp2 = new Vec3();
        const ids = new Array(n);
        for (let i = 0; i < n; i++) {
            const i3 = i * 3;
            let x = other.pos[i3];
            let y = other.pos[i3 + 1];
            let z = other.pos[i3 + 2];
            if (bRotMat3) {
                tmp.set(x, y, z);
                rot.mulVec(tmp, tmp2);
                x = tmp2.x; y = tmp2.y; z = tmp2.z;
            } else if (bRotFlat) {
                const rx = rot[0] * x + rot[1] * y + rot[2] * z;
                const ry = rot[3] * x + rot[4] * y + rot[5] * z;
                const rz = rot[6] * x + rot[7] * y + rot[8] * z;
                x = rx; y = ry; z = rz;
            }
            ids[i] = this.addAtom(x + px, y + py, z + pz, other.types[i]);
        }
        if (other.bonds && other.bonds.length) {
            for (const [a0, b0] of other.bonds) {
                const aId = ids[a0 | 0];
                const bId = ids[b0 | 0];
                if (aId === undefined || bId === undefined) throw new Error(`appendParsedSystem: bond refers to out-of-range atom a0=${a0} b0=${b0}`);
                this.addBond(aId, bId);
            }
        }
        return ids;
    }

    recalculateBonds(mmParams = null) {
        this._assertUnlocked('recalculateBonds');
        // remove all bonds
        while (this.bonds.length > 0) this.removeBondByIndex(this.bonds.length - 1);

        const defaultRcut = 1.6;
        const defaultRcut2 = defaultRcut * defaultRcut;
        const bondFactor = 1.3;

        const n = this.atoms.length;
        for (let i = 0; i < n; i++) {
            const ai = this.atoms[i];
            const zi = ai.Z;
            let r1 = 0.7;
            if (mmParams && mmParams.byAtomicNumber && mmParams.byAtomicNumber[zi]) {
                r1 = mmParams.byAtomicNumber[zi].Rcov;
            }
            const xi = ai.pos.x, yi = ai.pos.y, zi0 = ai.pos.z;
            for (let j = i + 1; j < n; j++) {
                const aj = this.atoms[j];
                const zj = aj.Z;
                let rCut2 = defaultRcut2;
                if (mmParams && mmParams.byAtomicNumber && mmParams.byAtomicNumber[zj]) {
                    const r2 = mmParams.byAtomicNumber[zj].Rcov;
                    const rSum = (r1 + r2) * bondFactor;
                    rCut2 = rSum * rSum;
                }
                const dx = xi - aj.pos.x;
                const dy = yi - aj.pos.y;
                const dz = zi0 - aj.pos.z;
                const dist2 = dx * dx + dy * dy + dz * dz;
                if (dist2 < rCut2) {
                    this.addBond(ai.id, aj.id);
                }
            }
        }
    }

    exportToMoleculeSystem(ms) {
        const n = this.atoms.length;
        if (ms.capacity < n) ms.resize(Math.max(n, ms.capacity * 2));

        if (!ms.atomIds || (ms.atomIds.length | 0) < (ms.capacity | 0)) ms.atomIds = new Int32Array(ms.capacity);
        ms._nextAtomId = (this.lastAtomId + 1) | 0;

        ms.nAtoms = n;
        for (let i = 0; i < n; i++) {
            const a = this.atoms[i];
            const i3 = i * 3;
            ms.pos[i3] = a.pos.x;
            ms.pos[i3 + 1] = a.pos.y;
            ms.pos[i3 + 2] = a.pos.z;
            ms.types[i] = a.Z;
            ms.atomIds[i] = a.id;
        }

        const bonds = new Array(this.bonds.length);
        for (let i = 0; i < this.bonds.length; i++) {
            const b = this.bonds[i];
            b.ensureIndices(this);
            bonds[i] = [b.a, b.b];
        }
        ms.setBonds(bonds);

        const sel = new Set();
        for (const id of this.selection) {
            const i = this.getAtomIndex(id);
            if (i >= 0) sel.add(i);
        }
        ms.selection = sel;

        ms.isDirty = true;
        this.dirtyTopo = false;
        this.dirtyGeom = false;
        this.dirtyExport = false;
    }

    toXYZString(opts = {}) {
        const qs = opts.qs || null;
        const bQ = !!qs;
        const lvec = (opts.lvec !== undefined) ? opts.lvec : (this.lvec || null);
        const out = [];
        out.push(String(this.atoms.length));
        if (lvec && lvec.length === 3) {
            const a = lvec[0], b = lvec[1], c = lvec[2];
            out.push(`lvs ${a.x} ${a.y} ${a.z}   ${b.x} ${b.y} ${b.z}   ${c.x} ${c.y} ${c.z}`);
        } else {
            out.push('Generated by EditableMolecule');
        }
        for (let i = 0; i < this.atoms.length; i++) {
            const at = this.atoms[i];
            const sym = EditableMolecule.Z_TO_SYMBOL[at.Z] || 'X';
            const x = at.pos.x.toFixed(6);
            const y = at.pos.y.toFixed(6);
            const z = at.pos.z.toFixed(6);
            if (bQ) out.push(`${sym} ${x} ${y} ${z} ${(qs[i] !== undefined) ? qs[i] : 0.0}`);
            else out.push(`${sym} ${x} ${y} ${z}`);
        }
        return out.join('\n') + '\n';
    }

    static normalizeSymbol(s) {
        if (!s) return s;
        const a = s.trim();
        if (a.length === 1) return a.toUpperCase();
        return a[0].toUpperCase() + a.slice(1).toLowerCase();
    }

    static symbolToZ(sym) {
        const s = EditableMolecule.normalizeSymbol(sym);
        const z = EditableMolecule.SYMBOL_TO_Z[s];
        if (!z) throw new Error(`symbolToZ: unknown symbol '${sym}'`);
        return z;
    }

    static asZ(x) {
        if (typeof x === 'number') return x | 0;
        if (typeof x !== 'string') throw new Error(`asZ: unsupported type ${typeof x}`);
        return EditableMolecule.symbolToZ(x);
    }

    static parseMol2(text) {
        const lines = text.split(/\r?\n/);
        let mode = '';
        let lvec = null;
        const apos = [];
        const types = [];
        const bonds = [];
        for (let il = 0; il < lines.length; il++) {
            const line = lines[il].trim();
            if (!line) continue;
            if (line.startsWith('@lvs')) {
                const p = line.split(/\s+/).slice(1).map(parseFloat);
                if (p.length >= 9) {
                    lvec = [
                        new Vec3(p[0], p[1], p[2]),
                        new Vec3(p[3], p[4], p[5]),
                        new Vec3(p[6], p[7], p[8])
                    ];
                }
                continue;
            }
            if (line.startsWith('@<TRIPOS>ATOM')) { mode = 'ATOM'; continue; }
            if (line.startsWith('@<TRIPOS>BOND')) { mode = 'BOND'; continue; }
            if (line.startsWith('@<TRIPOS>')) { mode = ''; continue; }
            if (mode === 'ATOM') {
                const parts = line.split(/\s+/);
                if (parts.length < 6) continue;
                const sym = EditableMolecule.normalizeSymbol(parts[1].replace(/[^A-Za-z].*$/, ''));
                const x = parseFloat(parts[2]);
                const y = parseFloat(parts[3]);
                const z = parseFloat(parts[4]);
                const Z = EditableMolecule.symbolToZ(sym);
                apos.push(x, y, z);
                types.push(Z);
            } else if (mode === 'BOND') {
                const parts = line.split(/\s+/);
                if (parts.length < 4) continue;
                const a = (parseInt(parts[1]) | 0) - 1;
                const b = (parseInt(parts[2]) | 0) - 1;
                if (a >= 0 && b >= 0) bonds.push([a, b]);
            }
        }
        return {
            pos: new Float32Array(apos),
            types: new Uint8Array(types),
            bonds,
            lvec
        };
    }

    static parseXYZ(text) {
        const lines = text.split(/\r?\n/);
        let i0 = 0;
        if (lines.length >= 2) i0 = 2;
        const pos = [];
        const types = [];
        for (let i = i0; i < lines.length; i++) {
            const line = lines[i].trim();
            if (!line) continue;
            const parts = line.split(/\s+/);
            if (parts.length < 4) continue;
            const sym = parts[0];
            const x = parseFloat(parts[1]);
            const y = parseFloat(parts[2]);
            const z = parseFloat(parts[3]);
            if (!isFinite(x) || !isFinite(y) || !isFinite(z)) continue;
            const Z = EditableMolecule.asZ(sym);
            pos.push(x, y, z);
            types.push(Z);
        }
        return { pos: new Float32Array(pos), types: new Uint8Array(types), bonds: [], lvec: null };
    }

    // These are intentionally small; expand only when needed.
    static SYMBOL_TO_Z = {
        H: 1, He: 2, C: 6, N: 7, O: 8, F: 9, Na: 11, Mg: 12, Al: 13, Si: 14, P: 15, S: 16, Cl: 17, K: 19, Ca: 20,
        Fe: 26, Cu: 29, Zn: 30, Se: 34, Br: 35, I: 53, Xe: 54
    };

    static Z_TO_SYMBOL = {
        1: 'H', 2: 'He', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl',
        19: 'K', 20: 'Ca', 26: 'Fe', 29: 'Cu', 30: 'Zn', 34: 'Se', 35: 'Br', 53: 'I', 54: 'Xe'
    };
}
