import { Vec3 } from "../../common_js/Vec3.js";
import { Mat3 } from "../../common_js/Mat3.js";

function _bondCut2(z1, z2, mmParams, defaultRcut2, bondFactor, stats) {
    if (stats) stats.nBondCut2 = (stats.nBondCut2 | 0) + 1;
    let r1 = 0.7;
    let r2 = 0.7;
    if (mmParams && mmParams.byAtomicNumber) {
        const e1 = mmParams.byAtomicNumber[z1 | 0];
        const e2 = mmParams.byAtomicNumber[z2 | 0];
        if (e1 && (e1.Rcov > 0)) r1 = e1.Rcov;
        if (e2 && (e2.Rcov > 0)) r2 = e2.Rcov;
        if (e1 && e2) {
            const rSum = (r1 + r2) * bondFactor;
            if (stats) stats.nRcovCutoff = (stats.nRcovCutoff | 0) + 1;
            return rSum * rSum;
        }
    }
    if (stats) stats.nDefaultCutoff = (stats.nDefaultCutoff | 0) + 1;
    return defaultRcut2;
}

function _maxRcovFromMolAtoms(atoms, mmParams, fallback = 0.7) {
    let r = +fallback;
    if (mmParams && mmParams.byAtomicNumber) {
        for (let i = 0; i < atoms.length; i++) {
            const z = atoms[i].Z | 0;
            const et = mmParams.byAtomicNumber[z];
            if (et && (et.Rcov > r)) r = et.Rcov;
        }
    }
    return r;
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

function _ensureMMTypeIndex(mmParams) {
    if (!mmParams || !mmParams.atomTypes) throw new Error('_ensureMMTypeIndex: mmParams with atomTypes required');
    if (mmParams._atomTypeNameToIndex && mmParams._atomTypeIndexToName) return;
    const names = Object.keys(mmParams.atomTypes);
    names.sort();
    const n2i = {};
    const i2n = names.slice();
    for (let i = 0; i < i2n.length; i++) n2i[i2n[i]] = i;
    mmParams._atomTypeNameToIndex = n2i;
    mmParams._atomTypeIndexToName = i2n;
}

function _atomTypeNameToIndex(mmParams, name) {
    _ensureMMTypeIndex(mmParams);
    const i = mmParams._atomTypeNameToIndex[name];
    return (i === undefined) ? -1 : (i | 0);
}

function _atomTypeIndexToName(mmParams, i) {
    _ensureMMTypeIndex(mmParams);
    const n = mmParams._atomTypeIndexToName[i | 0];
    return n ? n : null;
}

function _resolveTypeOrElementToAtomType(mmParams, nameOrEl) {
    if (!mmParams) throw new Error('_resolveTypeOrElementToAtomType: mmParams required');
    const s = String(nameOrEl || '').trim();
    if (!s) throw new Error('_resolveTypeOrElementToAtomType: empty type');
    if (mmParams.atomTypes && mmParams.atomTypes[s]) {
        const at = mmParams.atomTypes[s];
        return { name: s, atype: _atomTypeNameToIndex(mmParams, s), iZ: at.iZ | 0 };
    }
    if (mmParams.elementTypes && mmParams.elementTypes[s]) {
        const el = mmParams.elementTypes[s];
        const iZ = el.iZ | 0;
        if (mmParams.atomTypes && mmParams.atomTypes[s]) {
            const at = mmParams.atomTypes[s];
            return { name: s, atype: _atomTypeNameToIndex(mmParams, s), iZ: at.iZ | 0 };
        }
        if (mmParams.atomTypes) {
            for (const k in mmParams.atomTypes) {
                const at = mmParams.atomTypes[k];
                if (at && at.element_name === s) return { name: k, atype: _atomTypeNameToIndex(mmParams, k), iZ };
            }
        }
        return { name: s, atype: -1, iZ };
    }
    throw new Error(`Unknown cap type/element '${s}'`);
}

function _parseCountSet(s) {
    const t = String(s || '').trim();
    if (!t.startsWith('{') || !t.endsWith('}')) throw new Error(`selectQuery: expected count set {..}, got '${t}'`);
    const inner = t.slice(1, -1).trim();
    if (!inner) return new Set();
    const parts = inner.split(',');
    const out = new Set();
    for (let i = 0; i < parts.length; i++) {
        const x = parseInt(parts[i].trim(), 10);
        if (!isFinite(x)) throw new Error(`selectQuery: invalid count '${parts[i]}'`);
        out.add(x | 0);
    }
    return out;
}

function _compileTokenSetToMatcher(mmParams, tokSetStr) {
    const s = String(tokSetStr || '').trim();
    if (!s) throw new Error('selectQuery: empty token set');
    if (s === '*') return { bAny: true, zSet: null, tSet: null, fnMatch: (_a) => true };
    const toks = s.split('|').map(x => x.trim()).filter(x => x.length > 0);
    if (toks.length === 0) throw new Error('selectQuery: empty token set');
    const zSet = new Set();
    const tSet = new Set();
    for (let i = 0; i < toks.length; i++) {
        const tok = toks[i];
        if (tok === '*') return { bAny: true, zSet: null, tSet: null, fnMatch: (_a) => true };
        if (tok.indexOf('_') >= 0) {
            const it = _atomTypeNameToIndex(mmParams, tok);
            if (it < 0) throw new Error(`selectQuery: unknown atom type '${tok}'`);
            tSet.add(it | 0);
        } else {
            const z = EditableMolecule.asZ(tok);
            if (!(z > 0)) throw new Error(`selectQuery: invalid element '${tok}'`);
            zSet.add(z | 0);
        }
    }
    const fnMatch = (a) => {
        if (!a) return false;
        const iz = a.Z | 0;
        if (zSet.size && zSet.has(iz)) return true;
        const it = (a.atype !== undefined) ? (a.atype | 0) : -1;
        if (tSet.size && it >= 0 && tSet.has(it)) return true;
        return false;
    };
    return { bAny: false, zSet, tSet, fnMatch };
}

function _compileSelectQuery(mmParams, q) {
    const s = String(q || '').trim();
    if (!s) throw new Error('selectQuery: empty query');
    const parts = s.split(/\s+/).filter(x => x.length > 0);
    if (parts.length === 0) throw new Error('selectQuery: empty query');
    const atomTok = parts[0];
    const atomMatcher = _compileTokenSetToMatcher(mmParams, atomTok);
    const constraints = [];
    for (let i = 1; i < parts.length; i++) {
        const p = parts[i];
        if (!p) continue;
        if (p.startsWith('n{')) {
            const j = p.indexOf('}');
            if (j < 0) throw new Error(`selectQuery: missing '}' in '${p}'`);
            const neiSetStr = p.slice(2, j + 1);
            const rest = p.slice(j + 1);
            if (!rest.startsWith('={')) throw new Error(`selectQuery: expected '={' after n{..}, got '${p}'`);
            const cntSetStr = rest.slice(1);
            const cntSet = _parseCountSet(cntSetStr);
            const neiInner = neiSetStr.slice(2, -1).trim();
            const neiMatcher = _compileTokenSetToMatcher(mmParams, neiInner || '*');
            constraints.push({ neiMatcher, cntSet, src: p });
            continue;
        }
        if (p.startsWith('deg')) {
            const rest = p.slice(3);
            if (!rest.startsWith('={')) throw new Error(`selectQuery: expected deg={..}, got '${p}'`);
            const cntSet = _parseCountSet(rest.slice(1));
            const neiMatcher = _compileTokenSetToMatcher(mmParams, '*');
            constraints.push({ neiMatcher, cntSet, src: p });
            continue;
        }
        throw new Error(`selectQuery: cannot parse token '${p}'`);
    }
    return { q: s, atomMatcher, constraints };
}

function _getAtomTypeForAtom(mmParams, atom) {
    if (!mmParams) throw new Error('_getAtomTypeForAtom: mmParams required');
    if (!atom) throw new Error('_getAtomTypeForAtom: atom required');
    if (atom.atype !== undefined && (atom.atype | 0) >= 0) {
        const name = _atomTypeIndexToName(mmParams, atom.atype | 0);
        if (name && mmParams.atomTypes && mmParams.atomTypes[name]) return mmParams.atomTypes[name];
    }
    if (mmParams.byAtomicNumber) {
        const et = mmParams.byAtomicNumber[atom.Z | 0];
        if (et && mmParams.atomTypes && mmParams.atomTypes[et.name]) return mmParams.atomTypes[et.name];
    }
    if (mmParams.atomTypes && mmParams.atomTypes['*']) return mmParams.atomTypes['*'];
    throw new Error(`_getAtomTypeForAtom: cannot resolve atom type for Z=${atom.Z}`);
}

function _bondLengthEst(zA, zB, mmParams, bondFactor = 1.1, fallback = 1.0) {
    if (mmParams && mmParams.byAtomicNumber) {
        const e1 = mmParams.byAtomicNumber[zA | 0];
        const e2 = mmParams.byAtomicNumber[zB | 0];
        const r1 = (e1 && (e1.Rcov > 0)) ? e1.Rcov : 0.0;
        const r2 = (e2 && (e2.Rcov > 0)) ? e2.Rcov : 0.0;
        if (r1 > 0 && r2 > 0) return (r1 + r2) * bondFactor;
    }
    return +fallback;
}

function _orthonormalBasisFromDir(dir, outU, outV) {
    const a = Math.abs(dir.x), b = Math.abs(dir.y), c = Math.abs(dir.z);
    const tmp = (a < 0.9) ? new Vec3(1, 0, 0) : ((b < 0.9) ? new Vec3(0, 1, 0) : new Vec3(0, 0, 1));
    outU.setV(tmp);
    outU.subMul(dir, outU.dot(dir));
    if (!(outU.normalize() > 0)) throw new Error('_orthonormalBasisFromDir: failed to build basis');
    outV.setCross(dir, outU);
    if (!(outV.normalize() > 0)) throw new Error('_orthonormalBasisFromDir: failed to build basis (v)');
}

function _missingDirsVSEPR(vs, nMissing, totalDomains, outDirs) {
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
            throw new Error(`_missingDirsVSEPR: unsupported totalDomains=${totalDomains} for nb=0`);
        }
        while (outDirs.length > nMissing) outDirs.pop();
        return outDirs;
    }

    if (totalDomains === 2) {
        if (nb !== 1 || nMissing !== 1) throw new Error(`_missingDirsVSEPR: linear expects nb=1,nMissing=1 got nb=${nb} nMissing=${nMissing}`);
        outDirs.push(vs[0].clone().mulScalar(-1));
        return outDirs;
    }

    if (totalDomains === 3) {
        if (nb === 2 && nMissing === 1) {
            const m = vs[0].clone().add(vs[1]);
            if (!(m.normalize() > 0)) throw new Error('_missingDirsVSEPR: nb=2 planar but v1+v2 is zero');
            outDirs.push(m.mulScalar(-1));
            return outDirs;
        }
        if (nb === 1 && nMissing === 2) {
            const axis = vs[0].clone().mulScalar(-1);
            if (!(axis.normalize() > 0)) throw new Error('_missingDirsVSEPR: nb=1 planar axis zero');
            const u = new Vec3();
            const v = new Vec3();
            _orthonormalBasisFromDir(axis, u, v);
            const ca = -0.5;
            const sa = 0.86602540378;
            outDirs.push(axis.clone().mulScalar(ca).addMul(u, sa));
            outDirs.push(axis.clone().mulScalar(ca).addMul(u, -sa));
            return outDirs;
        }
        throw new Error(`_missingDirsVSEPR: unsupported planar nb=${nb} nMissing=${nMissing}`);
    }

    if (totalDomains === 4) {
        if (nb === 3 && nMissing === 1) {
            const m = vs[0].clone().add(vs[1]).add(vs[2]);
            if (!(m.normalize() > 0)) throw new Error('_missingDirsVSEPR: nb=3 tetra but sum is zero');
            outDirs.push(m.mulScalar(-1));
            return outDirs;
        }
        if (nb === 2 && nMissing === 2) {
            const m_c = vs[0].clone().add(vs[1]);
            if (!(m_c.normalize() > 0)) throw new Error('_missingDirsVSEPR: nb=2 tetra but v1+v2 is zero');
            const m_b = new Vec3().setCross(vs[0], vs[1]);
            if (!(m_b.normalize() > 0)) {
                const u = new Vec3();
                const v = new Vec3();
                _orthonormalBasisFromDir(m_c, u, v);
                m_b.setV(u);
            }
            const cc = 0.57735026919;
            const cb = 0.81649658092;
            outDirs.push(m_c.clone().mulScalar(-cc).addMul(m_b, cb).normalize() ? m_c.clone().mulScalar(-cc).addMul(m_b, cb).normalize() : null);
            outDirs.pop();
            const d1 = m_c.clone().mulScalar(-cc).addMul(m_b, cb);
            if (!(d1.normalize() > 0)) throw new Error('_missingDirsVSEPR: failed normalize tetra dir1');
            const d2 = m_c.clone().mulScalar(-cc).addMul(m_b, -cb);
            if (!(d2.normalize() > 0)) throw new Error('_missingDirsVSEPR: failed normalize tetra dir2');
            outDirs.push(d1);
            outDirs.push(d2);
            return outDirs;
        }
        if (nb === 1 && nMissing === 3) {
            const v1 = vs[0];
            const u = new Vec3();
            const v = new Vec3();
            _orthonormalBasisFromDir(v1, u, v);
            const a = -1.0 / 3.0;
            const b = Math.sqrt(8.0 / 9.0);
            const c120 = -0.5;
            const s120 = 0.86602540378;
            const u2 = u.clone().mulScalar(c120).addMul(v, s120);
            const u3 = u.clone().mulScalar(c120).addMul(v, -s120);
            const d1 = v1.clone().mulScalar(a).addMul(u, b);
            const d2 = v1.clone().mulScalar(a).addMul(u2, b);
            const d3 = v1.clone().mulScalar(a).addMul(u3, b);
            if (!(d1.normalize() > 0 && d2.normalize() > 0 && d3.normalize() > 0)) throw new Error('_missingDirsVSEPR: failed normalize tetra nb=1');
            outDirs.push(d1);
            outDirs.push(d2);
            outDirs.push(d3);
            return outDirs;
        }
        throw new Error(`_missingDirsVSEPR: unsupported tetra nb=${nb} nMissing=${nMissing}`);
    }

    throw new Error(`_missingDirsVSEPR: unsupported totalDomains=${totalDomains}`);
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

    setAtomTypeByName(id, typeName, mmParams) {
        const ia = this.getAtomIndex(id);
        if (ia < 0) throw new Error(`setAtomTypeByName: atom not found id=${id}`);
        const t = _resolveTypeOrElementToAtomType(mmParams, typeName);
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
        const capT = _resolveTypeOrElementToAtomType(mmParams, cap);
        const out = { nAdded: 0, capIds: [] };

        const vs = [];
        const dirs = [];
        const tmp = new Vec3();
        for (let ii = 0; ii < ids.length; ii++) {
            const id = ids[ii];
            const ia = this.getAtomIndex(id);
            if (ia < 0) continue;
            const a = this.atoms[ia];
            const at = _getAtomTypeForAtom(mmParams, a);
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

            _missingDirsVSEPR(vs, nDang, totalDomains, dirs);
            const r = _bondLengthEst(a.Z | 0, capT.iZ | 0, mmParams, bondFactor, 1.0);
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
            const at = _getAtomTypeForAtom(mmParams, a);
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
            _missingDirsVSEPR(vs, nAdd, totalDomains, dirs);

            const epName = (at.epair_name && at.epair_name !== '*') ? at.epair_name : 'E';
            const epT = _resolveTypeOrElementToAtomType(mmParams, epName);
            const r = _bondLengthEst(a.Z | 0, epT.iZ | 0, mmParams, bondFactor, 0.5);
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
        const tmp = vec;
        for (const id of ids) {
            const a = this.id2atom.get(id);
            if (a && a.pos) a.pos.add(tmp);
        }
        this._touchGeom();
    }

    rotateAtoms(ids, axis, deg, center = null) {
        if (!ids || ids.length === 0) return;
        const vAxis = (axis instanceof Vec3) ? axis.clone() : new Vec3(axis[0], axis[1], axis[2]);
        const ln = vAxis.normalize();
        if (!(ln > 0)) throw new Error('rotateAtoms: axis length is zero');
        const rad = deg * Math.PI / 180.0;
        const R = new Mat3();
        R.setRotate(rad, vAxis);
        let ctr = null;
        if (center instanceof Vec3) {
            ctr = center;
        } else if (center && center.length >= 3) {
            ctr = new Vec3(center[0], center[1], center[2]);
        } else {
            ctr = new Vec3();
            for (const id of ids) {
                const a = this.id2atom.get(id);
                if (a) ctr.add(a.pos);
            }
            ctr.mul(1.0 / ids.length);
        }
        for (const id of ids) {
            const a = this.id2atom.get(id);
            if (!a) continue;
            a.pos.sub(ctr);
            R.dotVec(a.pos, a.pos);
            a.pos.add(ctr);
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

    static compileSelectQuery(q, mmParams) {
        if (!mmParams) throw new Error('compileSelectQuery: mmParams required');
        return _compileSelectQuery(mmParams, q);
    }

    applySelectQuery(compiled, opts = {}) {
        const mode = (opts.mode !== undefined) ? String(opts.mode) : 'replace';
        const bPrint = !!opts.bPrint;
        if (!compiled || !compiled.atomMatcher || !compiled.constraints) throw new Error('applySelectQuery: invalid compiled query');
        const atomMatch = compiled.atomMatcher.fnMatch;
        const constraints = compiled.constraints;
        const sel = this.selection;
        if (mode === 'replace') sel.clear();
        const nAtoms = this.atoms.length | 0;
        const counts = new Int32Array(constraints.length | 0);
        let nHit = 0;
        for (let ia = 0; ia < nAtoms; ia++) {
            const a = this.atoms[ia];
            if (!a) continue;
            if (!atomMatch(a)) continue;
            for (let k = 0; k < counts.length; k++) counts[k] = 0;
            const bs = a.bonds;
            for (let ib0 = 0; ib0 < bs.length; ib0++) {
                const ib = bs[ib0] | 0;
                const bnd = this.bonds[ib];
                if (!bnd) continue;
                bnd.ensureIndices(this);
                const jb = bnd.other(ia);
                if (jb < 0) continue;
                const nb = this.atoms[jb];
                if (!nb) continue;
                for (let ic = 0; ic < constraints.length; ic++) {
                    if (constraints[ic].neiMatcher.fnMatch(nb)) counts[ic]++;
                }
            }
            let ok = true;
            for (let ic = 0; ic < constraints.length; ic++) {
                const cs = constraints[ic].cntSet;
                if (cs.size && !cs.has(counts[ic] | 0)) { ok = false; break; }
            }
            if (!ok) continue;
            const id = a.id;
            if (mode === 'subtract') sel.delete(id);
            else sel.add(id);
            nHit++;
            if (bPrint) console.log(`selectQuery hit id=${id} i=${ia} Z=${a.Z}`);
        }
        this.dirtyExport = true;
        return { nHit };
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

    attachGroupByMarker(groupParsed, markerX = 'Xe', markerY = 'He', opts = {}) {
        const maxIter = (opts.maxIter !== undefined) ? (opts.maxIter | 0) : 10000;
        const zX = EditableMolecule.asZ(markerX);
        const zY = EditableMolecule.asZ(markerY);
        const groupMarkerX = (opts.groupMarkerX !== undefined) ? opts.groupMarkerX : markerX;
        const groupMarkerY = (opts.groupMarkerY !== undefined) ? opts.groupMarkerY : markerY;
        const zGX = EditableMolecule.asZ(groupMarkerX);
        const zGY = EditableMolecule.asZ(groupMarkerY);
        if (!groupParsed || !groupParsed.pos || !groupParsed.types || !groupParsed.bonds) throw new Error('attachGroupByMarker: groupParsed must have pos/types/bonds');
        const gPairs = EditableMolecule._findMarkerPairsParsed(groupParsed, zGX, zGY);
        if (gPairs.length !== 1) throw new Error(`attachGroupByMarker: group must have exactly one marker pair, got ${gPairs.length}`);
        const gind = gPairs[0]; // [iX, iY, iA] in groupParsed local indices

        const pairs0 = EditableMolecule._findMarkerPairsMol(this, zX, zY);
        if (pairs0.length === 0) throw new Error(`attachGroupByMarker: backbone has no marker pairs X='${markerX}' Y='${markerY}'`);

        let it = 0;
        while (it < maxIter) {
            it++;
            const pairs = EditableMolecule._findMarkerPairsMol(this, zX, zY);
            if (pairs.length === 0) break;

            const bind = pairs[0]; // {xId, yId, aId}
            const R = EditableMolecule._computeMarkerAttachRotation(this, groupParsed, bind, gind, zX, zY, zGX, zGY);

            const Xb = this.atoms[this.getAtomIndex(bind.xId)].pos;
            const A2 = EditableMolecule._getParsedPos(groupParsed.pos, gind[2]);
            const tr = EditableMolecule._transformParsed(groupParsed, R, Xb, A2);
            const removed = EditableMolecule._removeAtomsFromParsed(tr, new Set([gind[0], gind[1]]));
            const idxAnchorGroup = removed.oldToNew[gind[2]];
            if (idxAnchorGroup < 0) throw new Error('attachGroupByMarker: group anchor unexpectedly removed');

            const ids = this.appendParsedSystem(removed, { pos: new Vec3(0, 0, 0) });
            const groupAnchorId = ids[idxAnchorGroup];
            if (groupAnchorId === undefined) throw new Error('attachGroupByMarker: groupAnchorId missing after append');
            this.addBond(bind.aId, groupAnchorId);

            // delete backbone marker atoms (order does not matter when deleting by ID)
            this.removeAtomById(bind.xId);
            this.removeAtomById(bind.yId);
        }
        if (it >= maxIter) throw new Error(`attachGroupByMarker: exceeded maxIter=${maxIter}`);
    }

    attachParsedByDirection(capAtom, groupParsed, params = {}) {
        const capId = capAtom | 0;
        const iCap = this.getAtomIndex(capId);
        if (iCap < 0) throw new Error(`attachParsedByDirection: capAtom not found id=${capAtom}`);
        if (!groupParsed || !groupParsed.pos || !groupParsed.types || !groupParsed.bonds) throw new Error('attachParsedByDirection: groupParsed must have pos/types/bonds');

        const cap = this.atoms[iCap];
        if (!cap.bonds || cap.bonds.length === 0) throw new Error('attachParsedByDirection: capAtom has no neighbors');

        const backId = (params.backAtom !== undefined) ? (params.backAtom | 0) : (() => {
            const ib = cap.bonds[0];
            const b = this.bonds[ib];
            b.ensureIndices(this);
            return this.atoms[b.other(iCap)].id;
        })();
        const iBack = this.getAtomIndex(backId);
        if (iBack < 0) throw new Error(`attachParsedByDirection: backAtom not found id=${backId}`);

        const bondLen = (params.bondLen !== undefined) ? +params.bondLen : 1.5;
        const upArr = (params.up !== undefined) ? params.up : [0, 0, 1];
        const up = new Vec3(+upArr[0], +upArr[1], +upArr[2]);
        const twistDeg = (params.twistDeg !== undefined) ? +params.twistDeg : 0.0;

        const iGa = (params.groupAnchor !== undefined) ? ((params.groupAnchor | 0) - 1) : 0;
        const iGf = (params.groupForwardRef !== undefined) ? ((params.groupForwardRef | 0) - 1) : 1;
        const iGu = (params.groupUpRef !== undefined) ? ((params.groupUpRef | 0) - 1) : -1;
        if (iGa < 0 || iGa >= groupParsed.types.length) throw new Error(`attachParsedByDirection: groupAnchor out of range ${iGa}`);
        if (iGf < 0 || iGf >= groupParsed.types.length) throw new Error(`attachParsedByDirection: groupForwardRef out of range ${iGf}`);
        if (iGu >= groupParsed.types.length) throw new Error(`attachParsedByDirection: groupUpRef out of range ${iGu}`);

        const Xb = this.atoms[iBack].pos;
        const Xcap = cap.pos;
        const f = new Vec3().setSub(Xcap, Xb);
        if (f.normalize() < 1e-12) throw new Error('attachParsedByDirection: cap-back vector is zero');

        const Mb0 = EditableMolecule._buildFrame(f, up);
        const Mb = (Math.abs(twistDeg) > 1e-9) ? EditableMolecule._rotateFrameAroundForward(Mb0, twistDeg * Math.PI / 180.0) : Mb0;

        const Xg = EditableMolecule._getParsedPos(groupParsed.pos, iGa);
        const PgF = EditableMolecule._getParsedPos(groupParsed.pos, iGf);
        const fg = new Vec3().setSub(PgF, Xg);

        let ug = null;
        if (iGu >= 0) {
            const PgU = EditableMolecule._getParsedPos(groupParsed.pos, iGu);
            ug = new Vec3().setSub(PgU, Xg);
        } else {
            ug = new Vec3(0, 0, 1);
        }

        const Mg = EditableMolecule._buildFrame(fg, ug);
        const R = Mat3.mul(Mb, Mg.clone().transpose());

        const Xattach = Xb.clone().addMul(f, bondLen);
        const tr = EditableMolecule._transformParsed(groupParsed, R, Xattach, Xg);
        const ids = this.appendParsedSystem(tr);
        const groupAnchorId = ids[iGa];
        if (groupAnchorId === undefined) throw new Error('attachParsedByDirection: groupAnchorId missing after append');
        this.addBond(backId, groupAnchorId);
        this.removeAtomById(capId);
    }

    static _getParsedPos(pos3, i) {
        const i3 = (i | 0) * 3;
        return new Vec3(pos3[i3], pos3[i3 + 1], pos3[i3 + 2]);
    }

    static _findMarkerPairsMol(mol, zX, zY) {
        const out = [];
        for (let i = 0; i < mol.atoms.length; i++) {
            const aX = mol.atoms[i];
            if (aX.Z !== zX) continue;
            let yId = -1;
            let aId = -1;
            for (const ib of aX.bonds) {
                const b = mol.bonds[ib];
                b.ensureIndices(mol);
                const j = b.other(i);
                const aj = mol.atoms[j];
                if (aj.Z === zY) yId = aj.id;
                else if (aj.Z !== zX) aId = aj.id;
            }
            if (yId >= 0 && aId >= 0) out.push({ xId: aX.id, yId, aId });
        }
        return out;
    }

    static _findMarkerPairsParsed(parsed, zX, zY) {
        const n = parsed.types.length | 0;
        const ngs = new Array(n).fill(null).map(() => []);
        for (const [a, b] of parsed.bonds) { ngs[a | 0].push(b | 0); ngs[b | 0].push(a | 0); }
        const out = [];
        for (let i = 0; i < n; i++) {
            if (parsed.types[i] !== zX) continue;
            let iY = -1;
            let iA = -1;
            for (const j of ngs[i]) {
                const tj = parsed.types[j];
                if (tj === zY) iY = j;
                else if (tj !== zX) iA = j;
            }
            if (iY >= 0 && iA >= 0) out.push([i, iY, iA]);
        }
        return out;
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

    static _computeMarkerAttachRotation(backboneMol, groupParsed, bind, gind, zBX, zBY, zGX, zGY) {
        const iXb = backboneMol.getAtomIndex(bind.xId);
        const iYb = backboneMol.getAtomIndex(bind.yId);
        const iAb = backboneMol.getAtomIndex(bind.aId);
        if (iXb < 0 || iYb < 0 || iAb < 0) throw new Error('_computeMarkerAttachRotation: backbone indices missing');
        const Xb = backboneMol.atoms[iXb].pos;
        const Yb = backboneMol.atoms[iYb].pos;
        const Ab = backboneMol.atoms[iAb].pos;

        const fb = new Vec3().setSub(Ab, Xb);
        const ub = new Vec3().setSub(Yb, Xb);
        const Mb = EditableMolecule._buildFrame(fb, ub);

        const Xg = EditableMolecule._getParsedPos(groupParsed.pos, gind[0]);
        const Yg = EditableMolecule._getParsedPos(groupParsed.pos, gind[1]);
        const Ag = EditableMolecule._getParsedPos(groupParsed.pos, gind[2]);

        const fg0 = new Vec3().setSub(Ag, Xg);
        if (fg0.normalize() < 1e-12) throw new Error('_computeMarkerAttachRotation: group forward is zero');
        const fg = fg0.mulScalar(-1.0);
        const ug = new Vec3().setSub(Yg, Xg);
        const Mg = EditableMolecule._buildFrame(fg, ug);

        return Mat3.mul(Mb, Mg.clone().transpose());
    }

    static _transformParsed(parsed, R, Xb, A2) {
        const n = parsed.types.length | 0;
        const pos = new Float32Array(n * 3);
        const tmp = new Vec3();
        const tmp2 = new Vec3();
        for (let i = 0; i < n; i++) {
            const i3 = i * 3;
            tmp.set(parsed.pos[i3] - A2.x, parsed.pos[i3 + 1] - A2.y, parsed.pos[i3 + 2] - A2.z);
            R.mulVec(tmp, tmp2);
            tmp2.add(Xb);
            pos[i3] = tmp2.x;
            pos[i3 + 1] = tmp2.y;
            pos[i3 + 2] = tmp2.z;
        }
        return { ...parsed, pos };
    }

    static _removeAtomsFromParsed(parsed, toRemoveSet) {
        const n = parsed.types.length | 0;
        const oldToNew = new Int32Array(n).fill(-1);
        let nNew = 0;
        for (let i = 0; i < n; i++) {
            if (!toRemoveSet.has(i)) { oldToNew[i] = nNew; nNew++; }
        }
        const pos = new Float32Array(nNew * 3);
        const types = new Uint8Array(nNew);
        for (let i = 0; i < n; i++) {
            const j = oldToNew[i];
            if (j < 0) continue;
            const i3 = i * 3;
            const j3 = j * 3;
            pos[j3] = parsed.pos[i3];
            pos[j3 + 1] = parsed.pos[i3 + 1];
            pos[j3 + 2] = parsed.pos[i3 + 2];
            types[j] = parsed.types[i];
        }
        const bonds = [];
        for (const [a, b] of parsed.bonds) {
            const na = oldToNew[a | 0];
            const nb = oldToNew[b | 0];
            if (na >= 0 && nb >= 0) bonds.push([na, nb]);
        }
        return { pos, types, bonds, lvec: parsed.lvec, oldToNew };
    }

    exportAsParsed(mol = this) {
        const n = mol.atoms.length;
        const pos = new Float32Array(n * 3);
        const types = new Uint8Array(n);
        const bonds = [];
        for (let i = 0; i < n; i++) {
            const a = mol.atoms[i];
            const i3 = i * 3;
            pos[i3] = a.pos.x; pos[i3 + 1] = a.pos.y; pos[i3 + 2] = a.pos.z;
            types[i] = a.Z;
        }
        for (let i = 0; i < mol.bonds.length; i++) {
            const b = mol.bonds[i];
            b.ensureIndices(mol);
            bonds.push([b.a, b.b]);
        }
        return { pos, types, bonds, lvec: mol.lvec ? [mol.lvec[0].clone(), mol.lvec[1].clone(), mol.lvec[2].clone()] : null };
    }

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
                const rCut2 = _bondCut2(ai.Z, aj.Z, mmParams, defaultRcut2, bondFactor, stats);
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
                        const rCut2 = _bondCut2(ai.Z, aj.Z, mmParams, defaultRcut2, bondFactor, stats);
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

        const maxR = _maxRcovFromMolAtoms(this.atoms, mmParams, 0.7);
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
                        const rCut2 = _bondCut2(ai.Z, aj.Z, mmParams, defaultRcut2, bondFactor, stats);
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

    toMol2String(opts = {}) {
        const lvec = (opts.lvec !== undefined) ? opts.lvec : (this.lvec || null);
        const name = (opts.name !== undefined) ? String(opts.name) : 'EditableMolecule';
        const out = [];
        out.push('@<TRIPOS>MOLECULE');
        out.push(name);
        out.push(` ${this.atoms.length} ${this.bonds.length} 0 0 0`);
        out.push('SMALL');
        out.push('GASTEIGER');
        if (lvec && lvec.length === 3) {
            const a = lvec[0], b = lvec[1], c = lvec[2];
            out.push(`@lvs ${a.x} ${a.y} ${a.z}    ${b.x} ${b.y} ${b.z}   ${c.x} ${c.y} ${c.z}`);
        }
        out.push('');
        out.push('@<TRIPOS>ATOM');
        for (let i = 0; i < this.atoms.length; i++) {
            const at = this.atoms[i];
            const sym = EditableMolecule.Z_TO_SYMBOL[at.Z] || 'X';
            const aname = `${sym}${i + 1}`;
            const x = at.pos.x.toFixed(4);
            const y = at.pos.y.toFixed(4);
            const z = at.pos.z.toFixed(4);
            const q = (at.charge !== undefined && at.charge !== null) ? (+at.charge).toFixed(4) : '0.0000';
            out.push(`${String(i + 1).padStart(7)} ${aname.padEnd(6)} ${x.padStart(10)} ${y.padStart(10)} ${z.padStart(10)} ${sym.padEnd(5)} 1  UNL1  ${q.padStart(10)}`);
        }
        out.push('@<TRIPOS>BOND');
        for (let i = 0; i < this.bonds.length; i++) {
            const b = this.bonds[i];
            b.ensureIndices(this);
            const a = b.a + 1;
            const c = b.b + 1;
            const ord = (b.order !== undefined && b.order !== null) ? String(b.order) : '1';
            out.push(`${String(i + 1).padStart(6)} ${String(a).padStart(5)} ${String(c).padStart(5)} ${ord}`);
        }
        out.push('');
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
