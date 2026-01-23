import { Vec3 } from "../common_js/Vec3.js";
import { Mat3 } from "../common_js/Mat3.js";
import { MMParams } from "./MMParams.js";

/// Build orthonormal basis perpendicular to direction dir (outputs u,v).
function orthonormalBasisFromDir(dir, outU, outV) {
    const a = Math.abs(dir.x), b = Math.abs(dir.y), c = Math.abs(dir.z);
    const tmp = (a < 0.9) ? new Vec3(1, 0, 0) : ((b < 0.9) ? new Vec3(0, 1, 0) : new Vec3(0, 0, 1));
    outU.setV(tmp);
    outU.subMul(dir, outU.dot(dir));
    if (!(outU.normalize() > 0)) throw new Error('orthonormalBasisFromDir: failed to build basis');
    outV.setCross(dir, outU);
    if (!(outV.normalize() > 0)) throw new Error('orthonormalBasisFromDir: failed to build basis (v)');
}

/// Compute missing VSEPR directions from existing neighbor vectors.
function missingDirsVSEPR(vs, nMissing, totalDomains, outDirs) {
    outDirs.length = 0;
    const nb = vs.length | 0;
    if (nMissing <= 0) return outDirs;
    if (nb === 0) {
        if (totalDomains === 4) {
            const a = 0.57735026919;
            outDirs.push(new Vec3( a, a, a));
            outDirs.push(new Vec3( a,-a,-a));
            outDirs.push(new Vec3(-a, a,-a));
            outDirs.push(new Vec3(-a,-a, a));
        } else if (totalDomains === 3) {
            outDirs.push(new Vec3(1, 0, 0));
            outDirs.push(new Vec3(-0.5, 0.86602540378, 0));
            outDirs.push(new Vec3(-0.5, -0.86602540378, 0));
        } else if (totalDomains === 2) {
            outDirs.push(new Vec3(1, 0, 0));
            outDirs.push(new Vec3(-1, 0, 0));
        } else {
            throw new Error(`missingDirsVSEPR: unsupported totalDomains=${totalDomains} for nb=0`);
        }
        while (outDirs.length > nMissing) outDirs.pop();
        return outDirs;
    }

    if (totalDomains === 2) {
        if (nb !== 1 || nMissing !== 1) throw new Error(`missingDirsVSEPR: linear expects nb=1,nMissing=1 got nb=${nb} nMissing=${nMissing}`);
        outDirs.push(vs[0].clone().mulScalar(-1));
        return outDirs;
    }

    if (totalDomains === 3) {
        if (nb === 2 && nMissing === 1) {
            const m = vs[0].clone().add(vs[1]);
            if (!(m.normalize() > 0)) throw new Error('missingDirsVSEPR: nb=2 planar but v1+v2 is zero');
            outDirs.push(m.mulScalar(-1));
            return outDirs;
        }
        if (nb === 1 && nMissing === 2) {
            const axis = vs[0].clone().mulScalar(-1);
            if (!(axis.normalize() > 0)) throw new Error('missingDirsVSEPR: nb=1 planar axis zero');
            const u = new Vec3();
            const v = new Vec3();
            orthonormalBasisFromDir(axis, u, v);
            const ca = -0.5;
            const sa = 0.86602540378;
            outDirs.push(axis.clone().mulScalar(ca).addMul(u, sa));
            outDirs.push(axis.clone().mulScalar(ca).addMul(u, -sa));
            return outDirs;
        }
        throw new Error(`missingDirsVSEPR: unsupported planar nb=${nb} nMissing=${nMissing}`);
    }

    if (totalDomains === 4) {
        if (nb === 3 && nMissing === 1) {
            const m = vs[0].clone().add(vs[1]).add(vs[2]);
            if (!(m.normalize() > 0)) throw new Error('missingDirsVSEPR: nb=3 tetra but sum is zero');
            outDirs.push(m.mulScalar(-1));
            return outDirs;
        }
        if (nb === 2 && nMissing === 2) {
            const m_c = vs[0].clone().add(vs[1]);
            if (!(m_c.normalize() > 0)) throw new Error('missingDirsVSEPR: nb=2 tetra but v1+v2 is zero');
            const m_b = new Vec3().setCross(vs[0], vs[1]);
            if (!(m_b.normalize() > 0)) {
                const u = new Vec3();
                const v = new Vec3();
                orthonormalBasisFromDir(m_c, u, v);
                m_b.setV(u);
            }
            const cc = 0.57735026919;
            const cb = 0.81649658092;
            outDirs.push(m_c.clone().mulScalar(-cc).addMul(m_b, cb).normalize() ? m_c.clone().mulScalar(-cc).addMul(m_b, cb).normalize() : null);
            outDirs.pop();
            const d1 = m_c.clone().mulScalar(-cc).addMul(m_b, cb);
            if (!(d1.normalize() > 0)) throw new Error('missingDirsVSEPR: failed normalize tetra dir1');
            const d2 = m_c.clone().mulScalar(-cc).addMul(m_b, -cb);
            if (!(d2.normalize() > 0)) throw new Error('missingDirsVSEPR: failed normalize tetra dir2');
            outDirs.push(d1);
            outDirs.push(d2);
            return outDirs;
        }
        if (nb === 1 && nMissing === 3) {
            const v1 = vs[0];
            const u = new Vec3();
            const v = new Vec3();
            orthonormalBasisFromDir(v1, u, v);
            const a = -1.0 / 3.0;
            const b = Math.sqrt(8.0 / 9.0);
            const c120 = -0.5;
            const s120 = 0.86602540378;
            const u2 = u.clone().mulScalar(c120).addMul(v, s120);
            const u3 = u.clone().mulScalar(c120).addMul(v, -s120);
            const d1 = v1.clone().mulScalar(a).addMul(u, b);
            const d2 = v1.clone().mulScalar(a).addMul(u2, b);
            const d3 = v1.clone().mulScalar(a).addMul(u3, b);
            if (!(d1.normalize() > 0 && d2.normalize() > 0 && d3.normalize() > 0)) throw new Error('missingDirsVSEPR: failed normalize tetra nb=1');
            outDirs.push(d1);
            outDirs.push(d2);
            outDirs.push(d3);
            return outDirs;
        }
        throw new Error(`missingDirsVSEPR: unsupported tetra nb=${nb} nMissing=${nMissing}`);
    }

    throw new Error(`missingDirsVSEPR: unsupported totalDomains=${totalDomains}`);
}

export class Atom {
    constructor(id, i, Z = 6, x = 0, y = 0, z = 0) {
        this.id = id;
        this.i = i;
        this.Z = Z;
        this.atype = -1;
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
        this.mmParams = null;

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

        this.lvec = [
            new Vec3(1, 0, 0),
            new Vec3(0, 1, 0),
            new Vec3(0, 0, 1)
        ];
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

    setAtomTypeByName(id, typeName, mmParams) {
        const ia = this.getAtomIndex(id);
        if (ia < 0) throw new Error(`setAtomTypeByName: atom not found id=${id}`);
        const t = mmParams.resolveTypeOrElementToAtomType(typeName);
        const a = this.atoms[ia];
        a.atype = t.atype;
        if (t.iZ > 0) a.Z = t.iZ;
        this._touchTopo();
        return t;
    }

    addCappingAtoms(mmParams, cap = 'H', opts = {}) {
        const bBond = (opts.bBond !== undefined) ? !!opts.bBond : true;
        const onlySelection = (opts.onlySelection !== undefined) ? !!opts.onlySelection : true;
        const bondFactor = (opts.bondFactor !== undefined) ? +opts.bondFactor : 1.1;
        const bPrint = !!opts.bPrint;
        const sel = Array.from(this.selection);
        if (onlySelection && sel.length === 0) throw new Error('addCappingAtoms: selection is empty');
        const ids = onlySelection ? sel : this.atoms.map(a => a.id);
        const capT = mmParams.resolveTypeOrElementToAtomType(cap);
        const out = { nAdded: 0, capIds: [] };

        const vs = [];
        const dirs = [];
        const tmp = new Vec3();
        for (let ii = 0; ii < ids.length; ii++) {
            const id = ids[ii];
            const ia = this.getAtomIndex(id);
            if (ia < 0) continue;
            const a = this.atoms[ia];
            const at = mmParams.getAtomTypeForAtom(a);
            const sigmaMax = (at.valence | 0);
            if (sigmaMax < 0) throw new Error(`addCappingAtoms: valence<0 for type='${at.name || '?'}' Z=${a.Z}`);
            let nbReal = 0;
            let nEp = 0;
            for (let k = 0; k < a.bonds.length; k++) {
                const ib = a.bonds[k] | 0;
                const bnd = this.bonds[ib];
                if (!bnd) continue;
                bnd.ensureIndices(this);
                const jb = bnd.other(ia);
                if (jb < 0) continue;
                const z = this.atoms[jb].Z | 0;
                if (z === 200) nEp++; else nbReal++;
            }
            const nDang = sigmaMax - nbReal;
            const nbAll = (nbReal + nEp) | 0;
            const totalDomains = (nbAll + nDang) | 0;
            if (bPrint) {
                const tname = (at && at.name) ? at.name : '?';
                console.log(`addCappingAtoms id=${a.id} i=${ia} Z=${a.Z} type=${tname} valence=${at.valence|0} nepair=${at.nepair|0} npi=${at.npi|0} nbReal=${nbReal} nbAll=${nbAll} nEp=${nEp} nDang=${nDang} totalDomains=${totalDomains} onlySelection=${onlySelection}`);
            }
            if (nDang <= 0) continue;

            if (totalDomains !== 4 && totalDomains !== 3 && totalDomains !== 2) {
                throw new Error(`addCappingAtoms: unsupported totalDomains=${totalDomains} (id=${a.id} Z=${a.Z} type='${at.name || '?'}' nbReal=${nbReal} nbAll=${nbAll} nDang=${nDang})`);
            }

            vs.length = 0;
            for (let k = 0; k < a.bonds.length; k++) {
                const ib = a.bonds[k] | 0;
                const bnd = this.bonds[ib];
                if (!bnd) continue;
                bnd.ensureIndices(this);
                const jb = bnd.other(ia);
                if (jb < 0) continue;
                tmp.setV(this.atoms[jb].pos);
                tmp.sub(a.pos);
                if (!(tmp.normalize() > 0)) continue;
                vs.push(tmp.clone());
            }

            missingDirsVSEPR(vs, nDang, totalDomains, dirs);
            const r = mmParams.bondLengthEstimate(a.Z | 0, capT.iZ | 0, bondFactor, 1.0);
            for (let j = 0; j < dirs.length; j++) {
                const d = dirs[j];
                const p = a.pos;
                const capId = this.addAtom(p.x + d.x * r, p.y + d.y * r, p.z + d.z * r, capT.iZ | 0);
                const iCap = this.getAtomIndex(capId);
                if (iCap < 0) throw new Error('addCappingAtoms: internal error (cap missing)');
                this.atoms[iCap].atype = capT.atype;
                out.capIds.push(capId);
                out.nAdded++;
                if (bBond) this.addBond(a.id, capId, 1, 1);
            }
        }
        return out;
    }

    addExplicitEPairs(mmParams, opts = {}) {
        const bBond = (opts.bBond !== undefined) ? !!opts.bBond : true;
        const onlySelection = (opts.onlySelection !== undefined) ? !!opts.onlySelection : true;
        const bondFactor = (opts.bondFactor !== undefined) ? +opts.bondFactor : 1.0;
        const sel = Array.from(this.selection);
        if (onlySelection && sel.length === 0) throw new Error('addExplicitEPairs: selection is empty');
        const ids = onlySelection ? sel : this.atoms.map(a => a.id);
        const out = { nAdded: 0, epairIds: [] };

        const vs = [];
        const dirs = [];
        const tmp = new Vec3();
        for (let ii = 0; ii < ids.length; ii++) {
            const id = ids[ii];
            const ia = this.getAtomIndex(id);
            if (ia < 0) continue;
            const a = this.atoms[ia];
            const at = getAtomTypeForAtom(mmParams, a);
            const ne = (at.nepair | 0);
            if (ne <= 0) continue;

            let nHave = 0;
            for (let k = 0; k < a.bonds.length; k++) {
                const ib = a.bonds[k] | 0;
                const bnd = this.bonds[ib];
                if (!bnd) continue;
                bnd.ensureIndices(this);
                const jb = bnd.other(ia);
                if (jb < 0) continue;
                if ((this.atoms[jb].Z | 0) === 200) nHave++;
            }
            const nAdd = ne - nHave;
            if (nAdd <= 0) continue;

            vs.length = 0;
            for (let k = 0; k < a.bonds.length; k++) {
                const ib = a.bonds[k] | 0;
                const bnd = this.bonds[ib];
                if (!bnd) continue;
                bnd.ensureIndices(this);
                const jb = bnd.other(ia);
                if (jb < 0) continue;
                tmp.setV(this.atoms[jb].pos);
                tmp.sub(a.pos);
                if (!(tmp.normalize() > 0)) continue;
                vs.push(tmp.clone());
            }

            const totalDomains = ((at.valence | 0) + (at.nepair | 0) + (at.npi | 0)) | 0;
            missingDirsVSEPR(vs, nAdd, totalDomains, dirs);

            const epName = (at.epair_name && at.epair_name !== '*') ? at.epair_name : 'E';
            const epT = mmParams.resolveTypeOrElementToAtomType(epName);
            const r = mmParams.bondLengthEstimate(a.Z | 0, epT.iZ | 0, bondFactor, 0.5);
            for (let j = 0; j < dirs.length; j++) {
                const d = dirs[j];
                const p = a.pos;
                const eId = this.addAtom(p.x + d.x * r, p.y + d.y * r, p.z + d.z * r, epT.iZ | 0);
                const iE = this.getAtomIndex(eId);
                if (iE < 0) throw new Error('addExplicitEPairs: internal error (epair missing)');
                this.atoms[iE].atype = epT.atype;
                out.epairIds.push(eId);
                out.nAdded++;
                if (bBond) this.addBond(a.id, eId, 1, 1);
            }
        }
        return out;
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

    translateAtoms(ids, vec) {
        if (!ids || ids.length === 0) return;
        if (!(vec instanceof Vec3)) throw new Error('translateAtoms: vec must be Vec3');
        for (const id of ids) {
            const ia = this.id2atom.get(id);
            const a = (ia !== undefined) ? this.atoms[ia] : null;
            if (a && a.pos) a.pos.add(vec);
            //console.log(`translateAtoms: translated atom ${id} to ${a.pos}`);
        }
        this._touchGeom();
    }

    rotateAtoms(ids, axis, deg, center) {
        if (!ids || ids.length === 0) return;
        if (!(axis instanceof Vec3)) throw new Error('rotateAtoms: axis must be Vec3');
        if (!(center instanceof Vec3)) throw new Error('rotateAtoms: center must be Vec3');
        const vAxis = axis.clone();
        const ln = vAxis.normalize();
        if (!(ln > 0)) throw new Error('rotateAtoms: axis length is zero');
        const rad = deg * Math.PI / 180.0;
        const R = Mat3.fromAxisAngle(vAxis, rad);
        for (const id of ids) {
            const ia = this.id2atom.get(id);
            const a = (ia !== undefined) ? this.atoms[ia] : null;
            if (a && a.pos) {
                a.pos.sub(center);
                R.mulVec(a.pos, a.pos);
                a.pos.add(center);
            }
        }
        this._touchGeom();
    }

    clearSelection() {
        this.selection.clear();
        this.isDirty = true;
        this.dirtyExport = true;
    }

    selectAll() {
        this.selection.clear();
        for (let i = 0; i < this.atoms.length; i++) {
            const a = this.atoms[i];
            if (a) this.selection.add(a.id);
        }
        this.isDirty = true;
        this.dirtyExport = true;
    }

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
        this.lastAtomId = 0;
        this.lastBondId = 0;
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

    /// Rebuild atom->bond adjacency from current bonds (safety net / bulk rebuild).
    updateNeighborList() {
        const n = this.atoms.length;
        for (let i = 0; i < n; i++) {
            const a = this.atoms[i];
            if (!a) throw new Error(`updateNeighborList: missing atom at index ${i}`);
            a.bonds.length = 0;
        }
        for (let ib = 0; ib < this.bonds.length; ib++) {
            const b = this.bonds[ib];
            if (!b) throw new Error(`updateNeighborList: missing bond at index ${ib}`);
            b.ensureIndices(this);
            const ia = b.a | 0;
            const ja = b.b | 0;
            if (ia < 0 || ia >= n || ja < 0 || ja >= n) throw new Error(`updateNeighborList: bond ${ib} has invalid atom indices a=${ia} b=${ja}`);
            this.atoms[ia].bonds.push(ib);
            this.atoms[ja].bonds.push(ib);
        }
    }

    printSizes(label = '') {
        const prefix = label ? `[EditableMolecule:${label}]` : '[EditableMolecule]';
        console.log(`${prefix} atoms=${this.atoms.length} bonds=${this.bonds.length} fragments=${this.fragments.length} selection=${this.selection.size}`);
        if (this.atoms.length > 0) {
            const sample = Math.min(5, this.atoms.length);
            for (let i = 0; i < sample; i++) {
                const a = this.atoms[i];
                if (!a) continue;
                const bondCount = a.bonds ? a.bonds.length : 0;
                console.log(`${prefix} atom[${i}] id=${a.id} Z=${a.Z} degree=${bondCount}`);
            }
        }
        if (this.bonds.length > 0) {
            const sampleB = Math.min(5, this.bonds.length);
            for (let ib = 0; ib < sampleB; ib++) {
                const b = this.bonds[ib];
                if (!b) continue;
                b.ensureIndices(this);
                console.log(`${prefix} bond[${ib}] id=${b.id} a=${b.a} b=${b.b}`);
            }
        }
    }

    _findBondId( ia, ib) {
        const a = this.atoms[ia];
        for (const ibond of a.bonds) {
            const b = this.bonds[ibond];
            if (!b) continue;
            b.ensureIndices(this);
            const other = b.other(ia);
            if (other === ib) return b.id;
        }
        return null;
    }

    _atomHeavyHydCounts( ia) {
        const a = this.atoms[ia];
        let heavy = 0, hyd = 0;
        for (const ib of a.bonds) {
            const b = this.bonds[ib];
            if (!b) continue;
            b.ensureIndices(this);
            const jb = b.other(ia);
            if (jb < 0 || jb >= this.atoms.length) continue;
            const nb = this.atoms[jb];
            if (!nb) continue;
            if (nb.Z === 1) hyd++;
            else heavy++;
        }
        return { heavy, hyd };
    }

    _sumHydrogenDirs( ia) {
        const a = this.atoms[ia];
        const acc = new Vec3();
        for (const ib of a.bonds) {
            const b = this.bonds[ib];
            if (!b) continue;
            b.ensureIndices(this);
            const jb = b.other(ia);
            if (jb < 0 || jb >= this.atoms.length) continue;
            const nb = this.atoms[jb];
            if (!nb || nb.Z !== 1) continue;
            const dir = new Vec3().setSub(nb.pos, a.pos);
            if (dir.normalize() > 1e-12) acc.add(dir);
        }
        return acc;
    }



    static _buildFrame(forward, up) {
        const f = forward.clone();
        if (f.normalize() < 1e-12) throw new Error('_buildFrame: zero forward');
        const u = up.clone();
        u.makeOrthoU(f);
        if (u.normalize() < 1e-12) throw new Error('_buildFrame: up parallel to forward');
        const l = new Vec3().setCross(u, f);
        if (l.normalize() < 1e-12) throw new Error('_buildFrame: failed to build left');
        return new Mat3(f, u, l);
    }

    static _rotateFrameAroundForward(M, phi) {
        const f = M.a.clone();
        const u = M.b.clone();
        const l = M.c.clone();
        const c = Math.cos(phi);
        const s = Math.sin(phi);
        const u2 = u.clone().mulScalar(c).addMul(l, s);
        const l2 = l.clone().mulScalar(c).addMul(u, -s);
        return new Mat3(f, u2, l2);
    }

    replicate(nrep, lvec = null) {
        this._assertUnlocked('replicate');
        const nx = (nrep[0] | 0);
        const ny = (nrep[1] | 0);
        const nz = (nrep[2] | 0);
        const nxyz = nx * ny * nz;
        if (nxyz <= 1) return;
        const lv = lvec || this.lvec;
        if (!lv || lv.length < 3) throw new Error('replicate: lattice vectors (lvec) required');

        const oldAtoms = this.atoms.slice();
        const oldBonds = this.bonds.slice();
        const na = oldAtoms.length;

        for (let iz = 0; iz < nz; iz++) {
            for (let iy = 0; iy < ny; iy++) {
                for (let ix = 0; ix < nx; ix++) {
                    if (ix === 0 && iy === 0 && iz === 0) continue;
                    const shift = new Vec3();
                    shift.addMul(lv[0], ix);
                    shift.addMul(lv[1], iy);
                    shift.addMul(lv[2], iz);
                    
                    const idMap = new Map();
                    for (let i = 0; i < na; i++) {
                        const a = oldAtoms[i];
                        const newId = this.addAtom(a.pos.x + shift.x, a.pos.y + shift.y, a.pos.z + shift.z, a.Z);
                        const iNew = this.getAtomIndex(newId);
                        this.atoms[iNew].atype = a.atype;
                        this.atoms[iNew].charge = a.charge;
                        if (a.cellIndex !== undefined) {
                            // Update cellIndex for the replicated atom to stay unique if it's a crystal replication
                            // This depends on how the caller expects cellIndex to behave.
                            // For simple replication, we might just copy it or offset it.
                            this.atoms[iNew].cellIndex = a.cellIndex; 
                        }
                        idMap.set(a.id, newId);
                    }
                    for (let i = 0; i < oldBonds.length; i++) {
                        const b = oldBonds[i];
                        this.addBond(idMap.get(b.aId), idMap.get(b.bId), b.order, b.type);
                    }
                }
            }
        }
        this.lvec = [
            lv[0].clone().mulScalar(nx),
            lv[1].clone().mulScalar(ny),
            lv[2].clone().mulScalar(nz)
        ];
        this._touchTopo();
    }

    recalculateBonds(mmParams = null, opts = {}) {
        this._assertUnlocked('recalculateBonds');
        // remove all bonds
        while (this.bonds.length > 0) this.removeBondByIndex(this.bonds.length - 1);

        const defaultRcut = (opts.defaultRcut !== undefined) ? +opts.defaultRcut : 1.6;
        const defaultRcut2 = defaultRcut * defaultRcut;
        const bondFactor = (opts.bondFactor !== undefined) ? +opts.bondFactor : 1.3;
        const stats = opts.stats ? opts.stats : null;

        const t0 = (opts.time !== false && typeof performance !== 'undefined') ? performance.now() : 0;
        const n = this.atoms.length;
        for (let i = 0; i < n; i++) {
            const ai = this.atoms[i];
            const xi = ai.pos.x, yi = ai.pos.y, zi0 = ai.pos.z;
            for (let j = i + 1; j < n; j++) {
                const aj = this.atoms[j];
                const rCut2 = mmParams.bondCutoff2(ai.Z, aj.Z, defaultRcut2, bondFactor, stats);
                if (stats) stats.nAtomPairs = (stats.nAtomPairs | 0) + 1;
                const dx = xi - aj.pos.x;
                const dy = yi - aj.pos.y;
                const dz = zi0 - aj.pos.z;
                const dist2 = dx * dx + dy * dy + dz * dz;
                if (stats) stats.nDist2 = (stats.nDist2 | 0) + 1;
                if (dist2 < rCut2) {
                    this.addBond(ai.id, aj.id);
                    if (stats) stats.nBondsAdded = (stats.nBondsAdded | 0) + 1;
                }
            }
        }
        if (opts.time !== false && typeof performance !== 'undefined' && stats) stats.timeMs = (performance.now() - t0);
        return stats;
    }

    recalculateBondsBucketNeighbors(mmParams, buckets, opts = {}) {
        this._assertUnlocked('recalculateBondsBucketNeighbors');
        if (!buckets || !buckets.buckets) throw new Error('recalculateBondsBucketNeighbors: buckets (BucketGraph) required');

        while (this.bonds.length > 0) this.removeBondByIndex(this.bonds.length - 1);

        const defaultRcut = (opts.defaultRcut !== undefined) ? +opts.defaultRcut : 1.6;
        const defaultRcut2 = defaultRcut * defaultRcut;
        const bondFactor = (opts.bondFactor !== undefined) ? +opts.bondFactor : 1.3;
        const stats = opts.stats ? opts.stats : null;
        const t0 = (opts.time !== false && typeof performance !== 'undefined') ? performance.now() : 0;

        const bs = buckets.buckets;
        for (let ib = 0; ib < bs.length; ib++) {
            const bi = bs[ib];
            const as = bi.atoms;
            const neigh = bi.neigh;
            for (let kk = 0; kk < neigh.length; kk++) {
                const jb = neigh[kk] | 0;
                const bj = bs[jb];
                const bsj = bj.atoms;
                if (stats) stats.nBucketPairs = (stats.nBucketPairs | 0) + 1;
                for (let ia0 = 0; ia0 < as.length; ia0++) {
                    const ia = as[ia0] | 0;
                    const ai = this.atoms[ia];
                    const xi = ai.pos.x, yi = ai.pos.y, zi0 = ai.pos.z;
                    const j0 = (jb === ib) ? (ia0 + 1) : 0;
                    for (let jb0 = j0; jb0 < bsj.length; jb0++) {
                        const ja = bsj[jb0] | 0;
                        const aj = this.atoms[ja];
                        const rCut2 = mmParams.bondCutoff2(ai.Z, aj.Z, defaultRcut2, bondFactor, stats);
                        if (stats) stats.nAtomPairs = (stats.nAtomPairs | 0) + 1;
                        const dx = xi - aj.pos.x;
                        const dy = yi - aj.pos.y;
                        const dz = zi0 - aj.pos.z;
                        const dist2 = dx * dx + dy * dy + dz * dz;
                        if (stats) stats.nDist2 = (stats.nDist2 | 0) + 1;
                        if (dist2 < rCut2) {
                            this.addBond(ai.id, aj.id);
                            if (stats) stats.nBondsAdded = (stats.nBondsAdded | 0) + 1;
                        }
                    }
                }
            }
        }
        if (opts.time !== false && typeof performance !== 'undefined' && stats) stats.timeMs = (performance.now() - t0);
        return stats;
    }

    recalculateBondsBucketAllPairsAABB(mmParams, buckets, opts = {}) {
        this._assertUnlocked('recalculateBondsBucketAllPairsAABB');
        if (!buckets || !buckets.buckets) throw new Error('recalculateBondsBucketAllPairsAABB: buckets (BucketGraph) required');

        while (this.bonds.length > 0) this.removeBondByIndex(this.bonds.length - 1);

        const defaultRcut = (opts.defaultRcut !== undefined) ? +opts.defaultRcut : 1.6;
        const defaultRcut2 = defaultRcut * defaultRcut;
        const bondFactor = (opts.bondFactor !== undefined) ? +opts.bondFactor : 1.3;
        const stats = opts.stats ? opts.stats : null;
        const t0 = (opts.time !== false && typeof performance !== 'undefined') ? performance.now() : 0;

        const maxR = mmParams.maxRcovFromAtoms(this.atoms, 0.7);
        const margin = (opts.margin !== undefined) ? +opts.margin : ((2.0 * maxR) * bondFactor);

        const bs = buckets.buckets;
        for (let ib = 0; ib < bs.length; ib++) {
            const bi = bs[ib];
            const as = bi.atoms;
            for (let jb = 0; jb <= ib; jb++) {
                const bj = bs[jb];
                if (stats) stats.nAABBTests = (stats.nAABBTests | 0) + 1;
                if (!bj.overlapAABB(bi, margin)) continue;
                if (stats) stats.nBucketPairs = (stats.nBucketPairs | 0) + 1;
                const bsj = bj.atoms;
                for (let ia0 = 0; ia0 < as.length; ia0++) {
                    const ia = as[ia0] | 0;
                    const ai = this.atoms[ia];
                    const xi = ai.pos.x, yi = ai.pos.y, zi0 = ai.pos.z;
                    const j0 = (jb === ib) ? (ia0 + 1) : 0;
                    for (let jb0 = j0; jb0 < bsj.length; jb0++) {
                        const ja = bsj[jb0] | 0;
                        const aj = this.atoms[ja];
                        const rCut2 = mmParams.bondCutoff2(ai.Z, aj.Z, defaultRcut2, bondFactor, stats);
                        if (stats) stats.nAtomPairs = (stats.nAtomPairs | 0) + 1;
                        const dx = xi - aj.pos.x;
                        const dy = yi - aj.pos.y;
                        const dz = zi0 - aj.pos.z;
                        const dist2 = dx * dx + dy * dy + dz * dz;
                        if (stats) stats.nDist2 = (stats.nDist2 | 0) + 1;
                        if (dist2 < rCut2) {
                            this.addBond(ai.id, aj.id);
                            if (stats) stats.nBondsAdded = (stats.nBondsAdded | 0) + 1;
                        }
                    }
                }
            }
        }
        if (opts.time !== false && typeof performance !== 'undefined' && stats) stats.timeMs = (performance.now() - t0);
        return stats;
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

    /// --- Moved to MoleculeSelection ---
    static compileSelectQuery(_q, _mmParams, _asZFn) { throw new Error('compileSelectQuery: install MoleculeSelection (installMoleculeSelectionMethods)'); }
    applySelectQuery(_compiled, _opts = {}) { throw new Error('applySelectQuery: install MoleculeSelection (installMoleculeSelectionMethods)'); }

    //  -- Moved to MoleculeUtils.js
    attachGroupByMarker          (_groupParsed, _markerX = 'Xe', _markerY = 'He', _opts = {} ) { throw new Error('attachGroupByMarker: install MoleculeUtils (installMoleculeUtilsMethods)');  }
    attachParsedByDirection      (_capAtom, _groupParsed, _params = {} )                       { throw new Error('attachParsedByDirection: install MoleculeUtils (installMoleculeUtilsMethods)'); }
    static _getParsedPos         (_pos3, _i)                                                   { throw new Error('_getParsedPos: install MoleculeUtils (installMoleculeUtilsMethods)'); }
    static _findMarkerPairsMol   (_mol, _zX, _zY)                                              { throw new Error('_findMarkerPairsMol: install MoleculeUtils (installMoleculeUtilsMethods)');}
    static _findMarkerPairsParsed(_parsed, _zX, _zY)                                           { throw new Error('_findMarkerPairsParsed: install MoleculeUtils (installMoleculeUtilsMethods)');}
    static _computeMarkerAttachRotation(_backboneMol, _groupParsed, _bind, _gind, _zBX, _zBY, _zGX, _zGY) { throw new Error('_computeMarkerAttachRotation: install MoleculeUtils (installMoleculeUtilsMethods)'); }
    static _transformParsed(_parsed, _R, _Xb, _A2)                                                        { throw new Error('_transformParsed: install MoleculeUtils (installMoleculeUtilsMethods)'); }
    static _removeAtomsFromParsed(_parsed, _toRemoveSet)                                                  {  throw new Error('_removeAtomsFromParsed: install MoleculeUtils (installMoleculeUtilsMethods)');  }

    // -- Moved to MoleculeIO.js
    static normalizeSymbol(_s)     { throw new Error('normalizeSymbol: install MoleculeIO (installMoleculeIOMethods)'); }
    static symbolToZ(_sym)         { throw new Error('symbolToZ: install MoleculeIO (installMoleculeIOMethods)'); }
    static asZ(_x)                 { throw new Error('asZ: install MoleculeIO (installMoleculeIOMethods)'); }
    static parseMol2(text)         { throw new Error('parseMol2: install MoleculeIO (installMoleculeIOMethods)'); }
    static parseXYZ(text)          { throw new Error('parseXYZ: install MoleculeIO (installMoleculeIOMethods)'); }

    // Symbol helpers provided by MoleculeIO installer.
    static SYMBOL_TO_Z = {};
    static Z_TO_SYMBOL = {};
}
