export class MoleculeSystem {
    constructor(capacity = 1000) {
        this.capacity = capacity;
        this.nAtoms = 0;

        // Structure of Arrays (SoA)
        this.pos = new Float32Array(this.capacity * 3); // x, y, z
        this.types = new Uint8Array(this.capacity);     // element type (atomic number)

        // Bonds (Simple list for now, can be optimized later)
        this.bonds = []; // Array of [id1, id2]

        // Dirty flags
        this.isDirty = true;

        // Selection
        this.selection = new Set();

        // Neighbor List (Adjacency List)
        // Array of Set<int> or Array<int>
        this.neighborList = [];
    }

    select(id, mode = 'replace') {
        if (mode === 'replace') {
            this.selection.clear();
            this.selection.add(id);
        } else if (mode === 'add') {
            this.selection.add(id);
        } else if (mode === 'subtract') {
            this.selection.delete(id);
        }
        this.isDirty = true; // Trigger redraw to update colors
    }

    clearSelection() {
        this.selection.clear();
        this.isDirty = true;
    }

    setBonds(bonds) {
        this.bonds = bonds ? bonds.map(b => [b[0], b[1]]) : [];
        this.updateNeighborList();
        this.isDirty = true;
    }

    addAtom(x, y, z, type) {
        if (this.nAtoms >= this.capacity) {
            this.resize(this.capacity * 2);
        }

        const i = this.nAtoms;
        this.pos[i * 3] = x;
        this.pos[i * 3 + 1] = y;
        this.pos[i * 3 + 2] = z;
        this.types[i] = type;

        // Initialize neighbor list for new atom
        if (this.neighborList.length <= i) {
            this.neighborList.push([]);
        } else {
            this.neighborList[i] = [];
        }

        this.nAtoms++;
        this.isDirty = true;
        return i;
    }

    addBond(id1, id2) {
        this.bonds.push([id1, id2]);

        // Update neighbor list
        if (this.neighborList[id1]) this.neighborList[id1].push(id2);
        if (this.neighborList[id2]) this.neighborList[id2].push(id1);

        this.isDirty = true;
    }

    updateNeighborList() {
        this.neighborList = new Array(this.nAtoms).fill(null).map(() => []);
        for (const [id1, id2] of this.bonds) {
            if (id1 < this.nAtoms && id2 < this.nAtoms) {
                this.neighborList[id1].push(id2);
                this.neighborList[id2].push(id1);
            }
        }
    }

    recalculateBonds(mmParams = null) {
        this.bonds = [];

        // Optimization: Use Spatial Hash if N > 1000 (TODO)

        const defaultRcut = 1.6;
        const defaultRcut2 = defaultRcut * defaultRcut;
        const bondFactor = 1.3; // Tolerance factor for covalent bonds

        for (let i = 0; i < this.nAtoms; i++) {
            const type1 = this.types[i];
            let r1 = 0.7; // Default ~Carbon
            if (mmParams && mmParams.byAtomicNumber[type1]) {
                r1 = mmParams.byAtomicNumber[type1].Rcov;
            }

            for (let j = i + 1; j < this.nAtoms; j++) {
                const type2 = this.types[j];
                let r2 = 0.7;
                let rCut2 = defaultRcut2;

                if (mmParams) {
                    if (mmParams.byAtomicNumber[type2]) {
                        r2 = mmParams.byAtomicNumber[type2].Rcov;
                    }
                    const rSum = (r1 + r2) * bondFactor;
                    rCut2 = rSum * rSum;
                }

                const dx = this.pos[i * 3] - this.pos[j * 3];
                const dy = this.pos[i * 3 + 1] - this.pos[j * 3 + 1];
                const dz = this.pos[i * 3 + 2] - this.pos[j * 3 + 2];

                const dist2 = dx * dx + dy * dy + dz * dz;
                if (dist2 < rCut2) {
                    this.bonds.push([i, j]);
                }
            }
        }
        this.updateNeighborList();
        this.isDirty = true;
        MoleculeSystem._log('info', `Recalculated bonds. Found ${this.bonds.length} bonds.`);
    }

    deleteSelectedAtoms() {
        if (this.selection.size === 0) return;

        const toDelete = Array.from(this.selection).sort((a, b) => b - a); // Delete from end to avoid index shift issues? 
        // Actually, swap-remove changes indices, so we need to be careful.
        // Better strategy: 
        // 1. Mark atoms to keep.
        // 2. Rebuild arrays.
        // 3. Re-map selection (clear it).

        const oldToNew = new Int32Array(this.nAtoms).fill(-1);
        let newCount = 0;

        // 1. Calculate new indices
        for (let i = 0; i < this.nAtoms; i++) {
            if (!this.selection.has(i)) {
                oldToNew[i] = newCount;
                newCount++;
            }
        }

        // 2. Compact Arrays
        const newPos = new Float32Array(this.capacity * 3);
        const newTypes = new Uint8Array(this.capacity);

        for (let i = 0; i < this.nAtoms; i++) {
            if (oldToNew[i] !== -1) {
                const newIdx = oldToNew[i];
                newPos[newIdx * 3] = this.pos[i * 3];
                newPos[newIdx * 3 + 1] = this.pos[i * 3 + 1];
                newPos[newIdx * 3 + 2] = this.pos[i * 3 + 2];
                newTypes[newIdx] = this.types[i];
            }
        }

        this.pos = newPos;
        this.types = newTypes;

        // 3. Rebuild Bonds
        const newBonds = [];
        for (const [id1, id2] of this.bonds) {
            if (oldToNew[id1] !== -1 && oldToNew[id2] !== -1) {
                newBonds.push([oldToNew[id1], oldToNew[id2]]);
            }
        }
        this.bonds = newBonds;

        this.nAtoms = newCount;
        this.selection.clear();
        this.updateNeighborList();
        this.isDirty = true;

        MoleculeSystem._log('info', `Deleted atoms. New count: ${this.nAtoms}`);
    }

    deleteAtoms(indices) {
        if (!indices) return;
        const rem = Array.isArray(indices) ? indices : Array.from(indices);
        if (rem.length === 0) return;
        const toRemove = new Set(rem.map(i => i | 0));
        const oldToNew = new Int32Array(this.nAtoms).fill(-1);
        let newCount = 0;
        for (let i = 0; i < this.nAtoms; i++) {
            if (!toRemove.has(i)) {
                oldToNew[i] = newCount;
                newCount++;
            }
        }
        const newPos = new Float32Array(this.capacity * 3);
        const newTypes = new Uint8Array(this.capacity);
        for (let i = 0; i < this.nAtoms; i++) {
            const j = oldToNew[i];
            if (j < 0) continue;
            const i3 = i * 3;
            const j3 = j * 3;
            newPos[j3] = this.pos[i3];
            newPos[j3 + 1] = this.pos[i3 + 1];
            newPos[j3 + 2] = this.pos[i3 + 2];
            newTypes[j] = this.types[i];
        }
        const newBonds = [];
        for (const [a, b] of this.bonds) {
            const na = oldToNew[a];
            const nb = oldToNew[b];
            if (na >= 0 && nb >= 0) newBonds.push([na, nb]);
        }
        this.pos = newPos;
        this.types = newTypes;
        this.bonds = newBonds;
        this.nAtoms = newCount;
        const newSel = new Set();
        for (const i of this.selection) {
            const j = oldToNew[i];
            if (j >= 0) newSel.add(j);
        }
        this.selection = newSel;
        this.updateNeighborList();
        this.isDirty = true;
        MoleculeSystem._log('info', `Deleted atoms. New count: ${this.nAtoms}`);
    }

    resize(newCapacity) {
        MoleculeSystem._log('info', `Resizing MoleculeSystem to ${newCapacity}`);
        const newPos = new Float32Array(newCapacity * 3);
        const newTypes = new Uint8Array(newCapacity);

        newPos.set(this.pos);
        newTypes.set(this.types);

        this.pos = newPos;
        this.types = newTypes;
        this.capacity = newCapacity;
    }

    clear() {
        this.nAtoms = 0;
        this.bonds = [];
        this.isDirty = true;
    }

    static _log(level, msg) {
        const g = (typeof globalThis !== 'undefined') ? globalThis : null;
        const l = g && g.logger ? g.logger : null;
        if (l) {
            if (level === 'error' && l.error) { l.error(msg); return; }
            if (level === 'warn' && l.warn) { l.warn(msg); return; }
            if (l.info) { l.info(msg); return; }
        }
        if (level === 'error') { console.error(msg); return; }
        if (level === 'warn') { console.warn(msg); return; }
        console.log(msg);
    }

    ensureAdditionalCapacity(nAdd) {
        const need = this.nAtoms + nAdd;
        if (need <= this.capacity) return;
        let cap = this.capacity;
        while (cap < need) cap = Math.max(1024, cap * 2);
        this.resize(cap);
    }

    addAtomsFromArrays(pos3, types1) {
        if (!pos3 || !types1) throw new Error('addAtomsFromArrays: pos3/types1 must be provided');
        const n = (types1.length | 0);
        if ((pos3.length | 0) !== n * 3) throw new Error(`addAtomsFromArrays: pos3.length=${pos3.length} != 3*types.length=${n * 3}`);
        this.ensureAdditionalCapacity(n);
        for (let i = 0; i < n; i++) {
            const i3 = i * 3;
            this.addAtom(pos3[i3], pos3[i3 + 1], pos3[i3 + 2], types1[i]);
        }
        this.isDirty = true;
    }

    appendParsedSystem(other, opts = {}) {
        if (!other || !other.pos || !other.types) throw new Error('appendParsedSystem: other must have pos/types');
        const n = other.types.length | 0;
        const offset = this.nAtoms;
        const pos = (opts.pos !== undefined) ? opts.pos : [0, 0, 0];
        const rot = (opts.rot !== undefined) ? opts.rot : null; // 3x3 row-major
        this.ensureAdditionalCapacity(n);
        for (let i = 0; i < n; i++) {
            const i3 = i * 3;
            let x = other.pos[i3];
            let y = other.pos[i3 + 1];
            let z = other.pos[i3 + 2];
            if (rot) {
                const rx = rot[0] * x + rot[1] * y + rot[2] * z;
                const ry = rot[3] * x + rot[4] * y + rot[5] * z;
                const rz = rot[6] * x + rot[7] * y + rot[8] * z;
                x = rx; y = ry; z = rz;
            }
            this.addAtom(x + pos[0], y + pos[1], z + pos[2], other.types[i]);
        }
        if (other.bonds && other.bonds.length) {
            for (const [a0, b0] of other.bonds) {
                this.bonds.push([(a0 | 0) + offset, (b0 | 0) + offset]);
            }
        }
        this.updateNeighborList();
        this.isDirty = true;
        return offset;
    }

    findMarkerPairs(markerX, markerY) {
        const zX = MoleculeSystem._asZ(markerX);
        const zY = MoleculeSystem._asZ(markerY);
        const out = [];
        for (let i = 0; i < this.nAtoms; i++) {
            if (this.types[i] !== zX) continue;
            const ngs = this.neighborList[i] || [];
            let iY = -1;
            let iA = -1;
            for (const j of ngs) {
                const tj = this.types[j];
                if (tj === zY) iY = j;
                else if (tj !== zX) iA = j;
            }
            if (iY >= 0 && iA >= 0) out.push([i, iY, iA]);
        }
        return out;
    }

    attachGroupByMarker(groupParsed, markerX = 'Xe', markerY = 'He', opts = {}) {
        const maxIter = (opts.maxIter !== undefined) ? (opts.maxIter | 0) : 10000;
        const zX = MoleculeSystem._asZ(markerX);
        const zY = MoleculeSystem._asZ(markerY);
        if (!groupParsed || !groupParsed.pos || !groupParsed.types || !groupParsed.bonds) throw new Error('attachGroupByMarker: groupParsed must have pos/types/bonds');
        const groupPairs = MoleculeSystem.findMarkerPairsParsed(groupParsed, zX, zY);
        if (groupPairs.length !== 1) throw new Error(`attachGroupByMarker: group must have exactly one marker pair, got ${groupPairs.length}`);
        const gind = groupPairs[0];
        let it = 0;
        while (it < maxIter) {
            it++;
            const pairs = this.findMarkerPairs(zX, zY);
            if (pairs.length === 0) break;
            const bind = pairs[0];
            const R = MoleculeSystem.computeGroupRotation(this, groupParsed, bind, gind);
            const Xb = MoleculeSystem._getPos3(this.pos, bind[0]);
            const A2 = MoleculeSystem._getPos3(groupParsed.pos, gind[2]);
            const tr = MoleculeSystem._transformParsed(groupParsed, R, Xb, A2);
            const removed = MoleculeSystem._removeAtomsFromParsed(tr, new Set([gind[0], gind[1]]));
            const idxAnchorGroup = removed.oldToNew[gind[2]];
            if (idxAnchorGroup < 0) throw new Error('attachGroupByMarker: group anchor unexpectedly removed');
            const offset = this.appendParsedSystem(removed);
            this.addBond(bind[2], offset + idxAnchorGroup);
            this.deleteAtoms([bind[0], bind[1]]);
        }
        if (it >= maxIter) throw new Error(`attachGroupByMarker: exceeded maxIter=${maxIter}`);
        this.updateNeighborList();
        this.isDirty = true;
    }

    static genReplicatedCell(params) {
        const { lvec, basisPos, basisTypes, nRep = [1, 1, 1], origin = [0, 0, 0] } = params || {};
        if (!lvec || !basisPos || !basisTypes) throw new Error('genReplicatedCell: lvec, basisPos, basisTypes required');
        const na = nRep[0] | 0;
        const nb = nRep[1] | 0;
        const nc = nRep[2] | 0;
        if (na <= 0 || nb <= 0 || nc <= 0) throw new Error(`genReplicatedCell: invalid nRep=${nRep}`);
        const nBasis = basisTypes.length | 0;
        if ((basisPos.length | 0) !== nBasis * 3) throw new Error('genReplicatedCell: basisPos must be flat array length 3*nBasis');
        const nTot = nBasis * na * nb * nc;
        const pos = new Float32Array(nTot * 3);
        const types = new Uint8Array(nTot);
        const ax = lvec[0][0], ay = lvec[0][1], az = lvec[0][2];
        const bx = lvec[1][0], by = lvec[1][1], bz = lvec[1][2];
        const cx = lvec[2][0], cy = lvec[2][1], cz = lvec[2][2];
        const ox = origin[0], oy = origin[1], oz = origin[2];
        let ia = 0;
        for (let iz = 0; iz < nc; iz++) {
            for (let iy = 0; iy < nb; iy++) {
                for (let ix = 0; ix < na; ix++) {
                    const sx = ox + ax * ix + bx * iy + cx * iz;
                    const sy = oy + ay * ix + by * iy + cy * iz;
                    const sz = oz + az * ix + bz * iy + cz * iz;
                    for (let ib = 0; ib < nBasis; ib++) {
                        const o3 = (ia * 3);
                        const b3 = (ib * 3);
                        pos[o3] = sx + basisPos[b3];
                        pos[o3 + 1] = sy + basisPos[b3 + 1];
                        pos[o3 + 2] = sz + basisPos[b3 + 2];
                        types[ia] = basisTypes[ib];
                        ia++;
                    }
                }
            }
        }
        return { pos, types, nTot, lvec: [
            [ax * na, ay * na, az * na],
            [bx * nb, by * nb, bz * nb],
            [cx * nc, cy * nc, cz * nc]
        ] };
    }

    static genNaClStep(params = {}) {
        const a = (params.a !== undefined) ? +params.a : (5.6413 / 2);
        const nx = (params.nx !== undefined) ? (params.nx | 0) : 13;
        const ny = (params.ny !== undefined) ? (params.ny | 0) : 12;
        const nz = (params.nz !== undefined) ? (params.nz | 0) : 3;
        const Q0 = (params.Q0 !== undefined) ? +params.Q0 : 0.7;
        if (nx <= 0 || ny <= 0 || nz <= 0) throw new Error(`genNaClStep: invalid nx,ny,nz = ${nx},${ny},${nz}`);
        const Lx = a * nx;
        const LxStep = Lx * 0.5;
        const ax = -1.0 / nx;
        const nTot = nx * ny * nz;
        const pos = new Float32Array(nTot * 3);
        const types = new Uint8Array(nTot);
        const qs = new Float32Array(nTot);
        let ia = 0;
        for (let iz = 0; iz < nz; iz++) {
            for (let ix = 0; ix < nx; ix++) {
                for (let iy = 0; iy < ny; iy++) {
                    let x = ix * a;
                    const y = iy * a;
                    let z = iz * a;
                    let i = ix + iy + iz;
                    const z_ = z + ax * x;
                    x -= ax * z;
                    z = z_;
                    if (x > LxStep) {
                        z += a;
                        i += 1;
                    }
                    const o3 = ia * 3;
                    pos[o3] = x;
                    pos[o3 + 1] = y;
                    pos[o3 + 2] = z;
                    if ((i & 1) === 0) {
                        types[ia] = 11;
                        qs[ia] = +Q0;
                    } else {
                        types[ia] = 17;
                        qs[ia] = -Q0;
                    }
                    ia++;
                }
            }
        }
        const lvec = [
            [nx * a, 0, 0],
            [0, ny * a, 0],
            [0, 0, nz * a]
        ];
        return { pos, types, qs, nTot, lvec, a, nx, ny, nz, Q0 };
    }

    static assemblePolymerSequence(sequence, monomers, opts = {}) {
        const _0 = (opts._0 !== undefined) ? (opts._0 | 0) : 1;
        const defaultCap = (opts.capacity !== undefined) ? (opts.capacity | 0) : 4096;
        if (!sequence) throw new Error('assemblePolymerSequence: sequence is empty');
        const sys = new MoleculeSystem(defaultCap);
        let pos = [0, 0, 0];
        let prev = null;
        let prevOffset = 0;
        for (let i = 0; i < sequence.length; i++) {
            const key = sequence[i];
            if (key === '_') continue;
            const rec = monomers[key];
            if (!rec || !rec.parsed) throw new Error(`assemblePolymerSequence: missing monomer '${key}' or its parsed data`);
            if (!rec.anchors || rec.anchors.length < 2) throw new Error(`assemblePolymerSequence: monomer '${key}' missing anchors [head,tail] (1-based)`);
            if (!prev) {
                sys.appendParsedSystem(rec.parsed);
                prev = rec;
                prevOffset = 0;
            } else {
                const lv = (rec.parsed.lvec && rec.parsed.lvec[1]) ? rec.parsed.lvec[1] : [0, 0, 0];
                pos = [pos[0] + lv[0], pos[1] + lv[1], pos[2] + lv[2]];
                const off = sys.appendParsedSystem(rec.parsed, { pos });
                const iPrevTail = prevOffset + ((prev.anchors[1] | 0) - _0);
                const iNewHead = off + ((rec.anchors[0] | 0) - _0);
                sys.addBond(iPrevTail, iNewHead);
                prev = rec;
                prevOffset = off;
            }
        }
        sys.updateNeighborList();
        sys.isDirty = true;
        return sys;
    }

    static assemblePolymerFromTokens(tokens, monomers, opts = {}) {
        const _0 = (opts._0 !== undefined) ? (opts._0 | 0) : 1;
        const defaultCap = (opts.capacity !== undefined) ? (opts.capacity | 0) : 4096;
        if (!tokens || tokens.length === 0) throw new Error('assemblePolymerFromTokens: tokens is empty');
        const sys = new MoleculeSystem(defaultCap);
        let pos = [0, 0, 0];
        let prev = null;
        let prevOffset = 0;
        for (let i = 0; i < tokens.length; i++) {
            const key = tokens[i];
            if (!key) continue;
            const rec = monomers[key];
            if (!rec || !rec.parsed) throw new Error(`assemblePolymerFromTokens: missing monomer '${key}' or its parsed data`);
            if (!rec.anchors || rec.anchors.length < 2) throw new Error(`assemblePolymerFromTokens: monomer '${key}' missing anchors [head,tail] (1-based)`);
            if (!prev) {
                sys.appendParsedSystem(rec.parsed);
                prev = rec;
                prevOffset = 0;
            } else {
                const lv = (rec.parsed.lvec && rec.parsed.lvec[1]) ? rec.parsed.lvec[1] : [0, 0, 0];
                pos = [pos[0] + lv[0], pos[1] + lv[1], pos[2] + lv[2]];
                const off = sys.appendParsedSystem(rec.parsed, { pos });
                const iPrevTail = prevOffset + ((prev.anchors[1] | 0) - _0);
                const iNewHead = off + ((rec.anchors[0] | 0) - _0);
                sys.addBond(iPrevTail, iNewHead);
                prev = rec;
                prevOffset = off;
            }
        }
        sys.updateNeighborList();
        sys.isDirty = true;
        return sys;
    }

    attachParsedByDirection(capAtom, groupParsed, params = {}) {
        const iCap = capAtom | 0;
        if (iCap < 0 || iCap >= this.nAtoms) throw new Error(`attachParsedByDirection: capAtom out of range ${capAtom}`);
        if (!groupParsed || !groupParsed.pos || !groupParsed.types || !groupParsed.bonds) throw new Error('attachParsedByDirection: groupParsed must have pos/types/bonds');
        this.updateNeighborList();
        const ngs = this.neighborList[iCap];
        if (!ngs || ngs.length === 0) throw new Error('attachParsedByDirection: capAtom has no neighbors');
        const iBack = (params.backAtom !== undefined) ? (params.backAtom | 0) : (ngs[0] | 0);
        if (iBack < 0 || iBack >= this.nAtoms) throw new Error(`attachParsedByDirection: backAtom out of range ${iBack}`);
        const bondLen = (params.bondLen !== undefined) ? +params.bondLen : 1.5;
        const up = (params.up !== undefined) ? params.up : [0, 0, 1];
        const twistDeg = (params.twistDeg !== undefined) ? +params.twistDeg : 0.0;
        const iGa = (params.groupAnchor !== undefined) ? ((params.groupAnchor | 0) - 1) : 0;
        const iGf = (params.groupForwardRef !== undefined) ? ((params.groupForwardRef | 0) - 1) : 1;
        const iGu = (params.groupUpRef !== undefined) ? ((params.groupUpRef | 0) - 1) : -1;
        if (iGa < 0 || iGa >= groupParsed.types.length) throw new Error(`attachParsedByDirection: groupAnchor out of range ${iGa}`);
        if (iGf < 0 || iGf >= groupParsed.types.length) throw new Error(`attachParsedByDirection: groupForwardRef out of range ${iGf}`);
        if (iGu >= groupParsed.types.length) throw new Error(`attachParsedByDirection: groupUpRef out of range ${iGu}`);
        const Xb = MoleculeSystem._getPos3(this.pos, iBack);
        const Xcap = MoleculeSystem._getPos3(this.pos, iCap);
        const f = MoleculeSystem._unit3(MoleculeSystem._sub3(Xcap, Xb));
        const Mb0 = MoleculeSystem.buildFrame(f, up);
        const Mb = (Math.abs(twistDeg) > 1e-9) ? MoleculeSystem._rotateFrameAroundForward(Mb0, twistDeg * Math.PI / 180.0) : Mb0;
        const Xg = MoleculeSystem._getPos3(groupParsed.pos, iGa);
        const PgF = MoleculeSystem._getPos3(groupParsed.pos, iGf);
        let fg = MoleculeSystem._sub3(PgF, Xg);
        let ug = null;
        if (iGu >= 0) {
            const PgU = MoleculeSystem._getPos3(groupParsed.pos, iGu);
            ug = MoleculeSystem._sub3(PgU, Xg);
        } else {
            ug = [0, 0, 1];
        }
        const Mg = MoleculeSystem.buildFrame(fg, ug);
        const R = MoleculeSystem._matMul3(Mb, MoleculeSystem._matT3(Mg));
        const Xattach = [Xb[0] + f[0] * bondLen, Xb[1] + f[1] * bondLen, Xb[2] + f[2] * bondLen];
        const tr = MoleculeSystem._transformParsed(groupParsed, R, Xattach, Xg);
        const offset = this.appendParsedSystem(tr);
        this.addBond(iBack, offset + iGa);
        this.deleteAtoms([iCap]);
        this.updateNeighborList();
        this.isDirty = true;
    }

    static _asZ(x) {
        if (typeof x === 'number') return x | 0;
        if (typeof x !== 'string') throw new Error(`_asZ: unsupported type ${typeof x}`);
        const s = x.trim();
        const z = MoleculeSystem.SYMBOL_TO_Z[s] || MoleculeSystem.SYMBOL_TO_Z[MoleculeSystem._normalizeSymbol(s)];
        if (!z) throw new Error(`Unknown element symbol '${x}'`);
        return z;
    }

    static _normalizeSymbol(s) {
        if (!s) return s;
        const a = s.trim();
        if (a.length === 1) return a.toUpperCase();
        return a[0].toUpperCase() + a.slice(1).toLowerCase();
    }

    static _getPos3(pos, i) {
        const i3 = (i | 0) * 3;
        return [pos[i3], pos[i3 + 1], pos[i3 + 2]];
    }

    static _sub3(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
    static _add3(a, b) { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
    static _dot3(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
    static _cross3(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
    static _norm3(a) { return Math.hypot(a[0], a[1], a[2]); }
    static _mul3s(a, s) { return [a[0] * s, a[1] * s, a[2] * s]; }
    static _unit3(a) { const r = MoleculeSystem._norm3(a); if (r < 1e-12) throw new Error('zero-length vector'); return [a[0] / r, a[1] / r, a[2] / r]; }

    static buildFrame(forward, up) {
        const f = MoleculeSystem._unit3(forward);
        let u = up;
        const fu = MoleculeSystem._dot3(u, f);
        u = MoleculeSystem._sub3(u, MoleculeSystem._mul3s(f, fu));
        u = MoleculeSystem._unit3(u);
        const l = MoleculeSystem._unit3(MoleculeSystem._cross3(u, f));
        return [f[0], u[0], l[0], f[1], u[1], l[1], f[2], u[2], l[2]];
    }

    static _rotateFrameAroundForward(M, phi) {
        const f = [M[0], M[3], M[6]];
        const u = [M[1], M[4], M[7]];
        const l = [M[2], M[5], M[8]];
        const c = Math.cos(phi);
        const s = Math.sin(phi);
        const u2 = [u[0] * c + l[0] * s, u[1] * c + l[1] * s, u[2] * c + l[2] * s];
        const l2 = [l[0] * c - u[0] * s, l[1] * c - u[1] * s, l[2] * c - u[2] * s];
        return [f[0], u2[0], l2[0], f[1], u2[1], l2[1], f[2], u2[2], l2[2]];
    }

    static _matMul3(A, B) {
        return [
            A[0] * B[0] + A[1] * B[3] + A[2] * B[6],
            A[0] * B[1] + A[1] * B[4] + A[2] * B[7],
            A[0] * B[2] + A[1] * B[5] + A[2] * B[8],
            A[3] * B[0] + A[4] * B[3] + A[5] * B[6],
            A[3] * B[1] + A[4] * B[4] + A[5] * B[7],
            A[3] * B[2] + A[4] * B[5] + A[5] * B[8],
            A[6] * B[0] + A[7] * B[3] + A[8] * B[6],
            A[6] * B[1] + A[7] * B[4] + A[8] * B[7],
            A[6] * B[2] + A[7] * B[5] + A[8] * B[8]
        ];
    }

    static _matT3(A) { return [A[0], A[3], A[6], A[1], A[4], A[7], A[2], A[5], A[8]]; }

    static computeGroupRotation(backboneSys, groupParsed, bind, gind) {
        const Xb = MoleculeSystem._getPos3(backboneSys.pos, bind[0]);
        const Ab = MoleculeSystem._getPos3(backboneSys.pos, bind[2]);
        const Yb = MoleculeSystem._getPos3(backboneSys.pos, bind[1]);
        const fb = MoleculeSystem._unit3(MoleculeSystem._sub3(Ab, Xb));
        const ub = MoleculeSystem._unit3(MoleculeSystem._sub3(Yb, Xb));
        const Mb = MoleculeSystem.buildFrame(fb, ub);
        const Xg = MoleculeSystem._getPos3(groupParsed.pos, gind[0]);
        const Ag = MoleculeSystem._getPos3(groupParsed.pos, gind[2]);
        const Yg = MoleculeSystem._getPos3(groupParsed.pos, gind[1]);
        const fg0 = MoleculeSystem._unit3(MoleculeSystem._sub3(Ag, Xg));
        const fg = MoleculeSystem._mul3s(fg0, -1.0);
        const ug = MoleculeSystem._unit3(MoleculeSystem._sub3(Yg, Xg));
        const Mg = MoleculeSystem.buildFrame(fg, ug);
        return MoleculeSystem._matMul3(Mb, MoleculeSystem._matT3(Mg));
    }

    static findMarkerPairsParsed(parsed, zX, zY) {
        const n = parsed.types.length | 0;
        const ngs = new Array(n).fill(null).map(() => []);
        for (const [a, b] of parsed.bonds) { ngs[a].push(b); ngs[b].push(a); }
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

    static _transformParsed(parsed, R, Xb, A2) {
        const n = parsed.types.length | 0;
        const pos = new Float32Array(n * 3);
        for (let i = 0; i < n; i++) {
            const i3 = i * 3;
            const x0 = parsed.pos[i3] - A2[0];
            const y0 = parsed.pos[i3 + 1] - A2[1];
            const z0 = parsed.pos[i3 + 2] - A2[2];
            pos[i3] = R[0] * x0 + R[1] * y0 + R[2] * z0 + Xb[0];
            pos[i3 + 1] = R[3] * x0 + R[4] * y0 + R[5] * z0 + Xb[1];
            pos[i3 + 2] = R[6] * x0 + R[7] * y0 + R[8] * z0 + Xb[2];
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
            const na = oldToNew[a];
            const nb = oldToNew[b];
            if (na >= 0 && nb >= 0) bonds.push([na, nb]);
        }
        return { pos, types, bonds, lvec: parsed.lvec, oldToNew };
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
                        [p[0], p[1], p[2]],
                        [p[3], p[4], p[5]],
                        [p[6], p[7], p[8]]
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
                const sym = MoleculeSystem._normalizeSymbol(parts[1].replace(/[^A-Za-z].*$/, ''));
                const x = parseFloat(parts[2]);
                const y = parseFloat(parts[3]);
                const z = parseFloat(parts[4]);
                const Z = MoleculeSystem._asZ(sym);
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

    toXYZString(opts = {}) {
        const bQ = !!opts.qs;
        const qs = opts.qs || null;
        const n = this.nAtoms;
        const out = [];
        out.push(String(n));
        if (opts.lvec && opts.lvec.length === 3) {
            const lv = opts.lvec;
            out.push(`lvs ${lv[0][0]} ${lv[0][1]} ${lv[0][2]}   ${lv[1][0]} ${lv[1][1]} ${lv[1][2]}   ${lv[2][0]} ${lv[2][1]} ${lv[2][2]}`);
        } else {
            out.push('Generated by MoleculeSystem');
        }
        for (let i = 0; i < n; i++) {
            const i3 = i * 3;
            const sym = MoleculeSystem.Z_TO_SYMBOL[this.types[i]] || 'X';
            const x = this.pos[i3].toFixed(6);
            const y = this.pos[i3 + 1].toFixed(6);
            const z = this.pos[i3 + 2].toFixed(6);
            if (bQ) out.push(`${sym} ${x} ${y} ${z} ${qs ? qs[i] : 0.0}`);
            else out.push(`${sym} ${x} ${y} ${z}`);
        }
        return out.join('\n') + '\n';
    }

    static SYMBOL_TO_Z = {
        H: 1, C: 6, N: 7, O: 8, F: 9, Na: 11, Mg: 12, Al: 13, Si: 14, P: 15, S: 16, Cl: 17, K: 19, Ca: 20,
        Fe: 26, Cu: 29, Zn: 30, Br: 35, I: 53, Se: 34, Xe: 54, He: 2
    };

    static Z_TO_SYMBOL = {
        1: 'H', 2: 'He', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl',
        19: 'K', 20: 'Ca', 26: 'Fe', 29: 'Cu', 30: 'Zn', 34: 'Se', 35: 'Br', 53: 'I', 54: 'Xe'
    };
}
