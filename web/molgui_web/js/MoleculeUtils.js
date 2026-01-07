import { EditableMolecule } from "./EditableMolecule.js";
import { Vec3 } from "../../common_js/Vec3.js";
import { Mat3 } from "../../common_js/Mat3.js";


/// Assemble a polymer molecule from monomer tokens and anchor definitions.
export function assemblePolymerFromTokens(tokens, monomers, opts = {}) {
    const _0 = (opts._0 !== undefined) ? (opts._0 | 0) : 1;
    if (!tokens || tokens.length === 0) throw new Error('assemblePolymerFromTokens: tokens is empty');

    const mol = new EditableMolecule();

    const pos = new Vec3(0, 0, 0);
    let prev = null;
    let prevTailId = -1;

    const _getParsedPos = (parsed, i) => {
        const i3 = (i | 0) * 3;
        return new Vec3(parsed.pos[i3], parsed.pos[i3 + 1], parsed.pos[i3 + 2]);
    };

    for (let i = 0; i < tokens.length; i++) {
        const key = tokens[i];
        if (!key) continue;
        const rec = monomers[key];
        if (!rec || !rec.parsed) throw new Error(`assemblePolymerFromTokens: missing monomer '${key}' or its parsed data`);
        if (!rec.anchors || rec.anchors.length < 2) throw new Error(`assemblePolymerFromTokens: monomer '${key}' missing anchors [head,tail] (1-based)`);

        if (!prev) {
            const ids = mol.appendParsedSystem(rec.parsed, { pos });
            const iTail = ((rec.anchors[1] | 0) - _0);
            if (iTail < 0 || iTail >= ids.length) throw new Error(`assemblePolymerFromTokens: tail anchor out of range for '${key}' iTail=${iTail}`);
            prevTailId = ids[iTail];
            prev = rec;
        } else {
            const lv = (rec.parsed.lvec && rec.parsed.lvec[1]) ? rec.parsed.lvec[1] : null;
            if (lv) {
                if (!(lv instanceof Vec3)) throw new Error('assemblePolymerFromTokens: parsed.lvec must be Vec3[3]');
                const iHead = ((rec.anchors[0] | 0) - _0);
                if (iHead < 0 || iHead >= rec.parsed.types.length) throw new Error(`assemblePolymerFromTokens: head anchor out of range for '${key}' iHead=${iHead}`);
                const iPrevTail = mol.getAtomIndex(prevTailId);
                if (iPrevTail < 0) throw new Error('assemblePolymerFromTokens: prevTailId not found');
                const pTail = mol.atoms[iPrevTail].pos;
                const pHeadLocal = _getParsedPos(rec.parsed, iHead);
                const pHeadPlus = pHeadLocal.clone().add(pos).add(lv);
                const pHeadMinus = pHeadLocal.clone().add(pos).sub(lv);
                const dPlus = pTail.dist2(pHeadPlus);
                const dMinus = pTail.dist2(pHeadMinus);
                if (dMinus < dPlus) pos.sub(lv);
                else pos.add(lv);
            }
            const ids = mol.appendParsedSystem(rec.parsed, { pos });
            const iHead = ((rec.anchors[0] | 0) - _0);
            const iTail = ((rec.anchors[1] | 0) - _0);
            if (iHead < 0 || iHead >= ids.length) throw new Error(`assemblePolymerFromTokens: head anchor out of range for '${key}' iHead=${iHead}`);
            if (iTail < 0 || iTail >= ids.length) throw new Error(`assemblePolymerFromTokens: tail anchor out of range for '${key}' iTail=${iTail}`);
            const headId = ids[iHead];
            mol.addBond(prevTailId, headId);
            prevTailId = ids[iTail];
            prev = rec;
        }
    }

    return mol;
}

/// Attach a parsed group to backbone markers X/Y (installs bonds and removes markers).
export function attachGroupByMarker(mol, groupParsed, markerX = 'Xe', markerY = 'He', opts = {}) {
    const maxIter = (opts.maxIter !== undefined) ? (opts.maxIter | 0) : 10000;
    const zX = EditableMolecule.asZ(markerX);
    const zY = EditableMolecule.asZ(markerY);
    const groupMarkerX = (opts.groupMarkerX !== undefined) ? opts.groupMarkerX : markerX;
    const groupMarkerY = (opts.groupMarkerY !== undefined) ? opts.groupMarkerY : markerY;
    const zGX = EditableMolecule.asZ(groupMarkerX);
    const zGY = EditableMolecule.asZ(groupMarkerY);
    if (!groupParsed || !groupParsed.pos || !groupParsed.types || !groupParsed.bonds) throw new Error('attachGroupByMarker: groupParsed must have pos/types/bonds');
    const gPairs = _findMarkerPairsParsed(groupParsed, zGX, zGY);
    if (gPairs.length !== 1) throw new Error(`attachGroupByMarker: group must have exactly one marker pair, got ${gPairs.length}`);
    const gind = gPairs[0];

    const pairs0 = _findMarkerPairsMol(mol, zX, zY);
    if (pairs0.length === 0) throw new Error(`attachGroupByMarker: backbone has no marker pairs X='${markerX}' Y='${markerY}'`);

    let it = 0;
    while (it < maxIter) {
        it++;
        const pairs = _findMarkerPairsMol(mol, zX, zY);
        if (pairs.length === 0) break;

        const bind = pairs[0];
        const R = _computeMarkerAttachRotation(mol, groupParsed, bind, gind, zX, zY, zGX, zGY);

        const Xb = mol.atoms[mol.getAtomIndex(bind.xId)].pos;
        const A2 = _getParsedPos(groupParsed.pos, gind[2]);
        const tr = _transformParsed(groupParsed, R, Xb, A2);
        const removed = _removeAtomsFromParsed(tr, new Set([gind[0], gind[1]]));
        const idxAnchorGroup = removed.oldToNew[gind[2]];
        if (idxAnchorGroup < 0) throw new Error('attachGroupByMarker: group anchor unexpectedly removed');

        const ids = mol.appendParsedSystem(removed, { pos: new Vec3(0, 0, 0) });
        const groupAnchorId = ids[idxAnchorGroup];
        if (groupAnchorId === undefined) throw new Error('attachGroupByMarker: groupAnchorId missing after append');
        mol.addBond(bind.aId, groupAnchorId);

        mol.removeAtomById(bind.xId);
        mol.removeAtomById(bind.yId);
    }
    if (it >= maxIter) throw new Error(`attachGroupByMarker: exceeded maxIter=${maxIter}`);
}

/// Append another EditableMolecule by exporting its parsed data and reusing appendParsedSystem.
export function appendMolecule(mol, otherMol, opts = {}) {
    if (!otherMol || !otherMol.exportAsParsed) throw new Error('appendMolecule: otherMol must be EditableMolecule with exportAsParsed');
    const parsed = otherMol.exportAsParsed();
    return mol.appendParsedSystem(parsed, opts);
}

/// Attach another EditableMolecule via marker matching (exports parsed first).
export function attachMoleculeByMarker(mol, otherMol, markerX = 'Xe', markerY = 'He', opts = {}) {
    if (!otherMol || !otherMol.exportAsParsed) throw new Error('attachMoleculeByMarker: otherMol must be EditableMolecule with exportAsParsed');
    const parsed = otherMol.exportAsParsed();
    return attachGroupByMarker(mol, parsed, markerX, markerY, opts);
}

/// Attach another EditableMolecule via directional frame (exports parsed first).
export function attachMoleculeByDirection(mol, otherMol, capAtom, params = {}) {
    if (!otherMol || !otherMol.exportAsParsed) throw new Error('attachMoleculeByDirection: otherMol must be EditableMolecule with exportAsParsed');
    const parsed = otherMol.exportAsParsed();
    return attachParsedByDirection(mol, capAtom, parsed, params);
}

/// Attach parsed fragment to a cap atom following forward/up references.
export function attachParsedByDirection(mol, capAtom, groupParsed, params = {}) {
    const capId = capAtom | 0;
    const iCap = mol.getAtomIndex(capId);
    if (iCap < 0) throw new Error(`attachParsedByDirection: capAtom not found id=${capAtom}`);
    if (!groupParsed || !groupParsed.pos || !groupParsed.types || !groupParsed.bonds) throw new Error('attachParsedByDirection: groupParsed must have pos/types/bonds');

    const cap = mol.atoms[iCap];
    if (!cap.bonds || cap.bonds.length === 0) throw new Error('attachParsedByDirection: capAtom has no neighbors');

    const backId = (params.backAtom !== undefined) ? (params.backAtom | 0) : (() => {
        const ib = cap.bonds[0];
        const b = mol.bonds[ib];
        b.ensureIndices(mol);
        return mol.atoms[b.other(iCap)].id;
    })();
    const iBack = mol.getAtomIndex(backId);
    if (iBack < 0) throw new Error(`attachParsedByDirection: backAtom not found id=${backId}`);

    const upId = (params.upAtom !== undefined) ? (params.upAtom | 0) : backId;
    const iUp = mol.getAtomIndex(upId);
    if (iUp < 0) throw new Error(`attachParsedByDirection: upAtom not found id=${upId}`);

    const forwardRef = (params.forwardRef !== undefined) ? params.forwardRef : 0;
    const upRef = (params.upRef !== undefined) ? params.upRef : 1;
    const anchorRef = (params.anchorRef !== undefined) ? params.anchorRef : 2;
    const twist = (params.twist !== undefined) ? +params.twist : 0.0;
    const zBack = (params.zBack !== undefined) ? params.zBack : 6;

    const iForward = forwardRef | 0;
    const iUpRef = upRef | 0;
    const iAnchor = anchorRef | 0;
    const pAnchor = _getParsedPos(groupParsed.pos, iAnchor);
    const pForward = _getParsedPos(groupParsed.pos, iForward);
    const pUp = _getParsedPos(groupParsed.pos, iUpRef);
    const fg = pForward.clone().sub(pAnchor).normalize();
    const ug = pUp.clone().sub(pAnchor).normalize();

    const Xb = mol.atoms[iCap].pos;
    const Ab = mol.atoms[iBack].pos;
    const Ub = mol.atoms[iUp].pos;
    const fb = Ab.clone().sub(Xb).normalize();
    const ub = Ub.clone().sub(Xb).normalize();

    let R = Mat3.fromForwardUp(fb, ub);
    const Mg = Mat3.fromForwardUp(fg, ug);
    R = Mat3.mul(R, Mg.clone().transpose());
    if (Math.abs(twist) > 1e-12) R = Mat3.rotateAroundForward(R, twist);

    const tr = _transformParsed(groupParsed, R, Xb, pAnchor);
    const ids = mol.appendParsedSystem(tr, { pos: new Vec3(0, 0, 0) });
    mol.addBond(capId, ids[iForward]);
    if (zBack > 0) mol.addBond(backId, ids[iAnchor]);
    return ids;
}

/// Replicate molecule over lattice repeats (nrep[3]) updating lvec.
export function replicateMolecule(mol, nrep, lvec = null) {
    mol._assertUnlocked('replicate');
    const nx = (nrep[0] | 0);
    const ny = (nrep[1] | 0);
    const nz = (nrep[2] | 0);
    const nxyz = nx * ny * nz;
    if (nxyz <= 1) return;
    const lv = lvec || mol.lvec;
    if (!lv || lv.length < 3) throw new Error('replicate: lattice vectors (lvec) required');

    const oldAtoms = mol.atoms.slice();
    const oldBonds = mol.bonds.slice();
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
                    const newId = mol.addAtom(a.pos.x + shift.x, a.pos.y + shift.y, a.pos.z + shift.z, a.Z);
                    const iNew = mol.getAtomIndex(newId);
                    mol.atoms[iNew].atype = a.atype;
                    mol.atoms[iNew].charge = a.charge;
                    idMap.set(a.id, newId);
                }
                for (let i = 0; i < oldBonds.length; i++) {
                    const b = oldBonds[i];
                    mol.addBond(idMap.get(b.aId), idMap.get(b.bId), b.order, b.type);
                }
            }
        }
    }
    mol.lvec = [
        lv[0].clone().mulScalar(nx),
        lv[1].clone().mulScalar(ny),
        lv[2].clone().mulScalar(nz)
    ];
    mol._touchTopo();
}

function _getParsedPos(parsedPos, i) {
    const i3 = (i | 0) * 3;
    return new Vec3(parsedPos[i3], parsedPos[i3 + 1], parsedPos[i3 + 2]);
}

function _findMarkerPairsMol(mol, zX, zY) {
    const out = [];
    for (let ia = 0; ia < mol.atoms.length; ia++) {
        const a = mol.atoms[ia];
        if (!a || a.Z !== zX) continue;
        const nb = a.bonds;
        for (let ib0 = 0; ib0 < nb.length; ib0++) {
            const ib = nb[ib0] | 0;
            const bnd = mol.bonds[ib];
            if (!bnd) continue;
            bnd.ensureIndices(mol);
            const jb = bnd.other(ia);
            if (jb < 0) continue;
            const nbAtom = mol.atoms[jb];
            if (!nbAtom || nbAtom.Z !== zY) continue;
            const nnb = nbAtom.bonds;
            for (let ic0 = 0; ic0 < nnb.length; ic0++) {
                const ic = nnb[ic0] | 0;
                const b2 = mol.bonds[ic];
                if (!b2) continue;
                b2.ensureIndices(mol);
                const jc = b2.other(jb);
                if (jc < 0 || jc === ia) continue;
                const aId = mol.atoms[jc].id;
                out.push({ xId: a.id, yId: nbAtom.id, aId });
            }
        }
    }
    return out;
}

function _findMarkerPairsParsed(parsed, zX, zY) {
    const out = [];
    const n = parsed.types.length | 0;
    for (let i = 0; i < n; i++) {
        if (parsed.types[i] !== zX) continue;
        for (const [a, b] of parsed.bonds) {
            const ia = a | 0;
            const ib = b | 0;
            if (ia !== i && ib !== i) continue;
            const j = (ia === i) ? ib : ia;
            if ((parsed.types[j] | 0) !== zY) continue;
            for (const [c, d] of parsed.bonds) {
                const ic = c | 0;
                const id = d | 0;
                if (ic === j && id !== i) { out.push([i, j, id]); break; }
                if (id === j && ic !== i) { out.push([i, j, ic]); break; }
            }
        }
    }
    return out;
}

function _computeMarkerAttachRotation(mol, groupParsed, bind, gind, zX, zY, zGX, zGY) {
    const Xb = mol.atoms[mol.getAtomIndex(bind.xId)].pos;
    const Yb = mol.atoms[mol.getAtomIndex(bind.yId)].pos;
    const Ab = mol.atoms[mol.getAtomIndex(bind.aId)].pos;
    const fg = _getParsedPos(groupParsed.pos, gind[0]).clone().sub(_getParsedPos(groupParsed.pos, gind[2])).normalize();
    const ug = _getParsedPos(groupParsed.pos, gind[1]).clone().sub(_getParsedPos(groupParsed.pos, gind[2])).normalize();
    const fb = Xb.clone().sub(Ab).normalize();
    const ub = Yb.clone().sub(Ab).normalize();
    const Mb = Mat3.fromForwardUp(fb, ub);
    const Mg = Mat3.fromForwardUp(fg, ug);
    return Mat3.mul(Mb, Mg.clone().transpose());
}

function _transformParsed(parsed, R, Xb, A2) {
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

function _removeAtomsFromParsed(parsed, toRemoveSet) {
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

/// Install utility helpers onto EditableMolecule (polymer/attach/replicate).
export function installMoleculeUtilsMethods(cls = EditableMolecule) {
    console.log('[installMoleculeUtilsMethods] installing on', cls && cls.name);
    if (!cls || !cls.prototype) throw new Error('installMoleculeUtilsMethods: invalid class');
    cls.assemblePolymerFromTokens = assemblePolymerFromTokens;
    cls.prototype.assemblePolymerFromTokens = function (tokens, monomers, opts = {}) {
        return assemblePolymerFromTokens(tokens, monomers, opts);
    };
    cls.prototype.appendMolecule = function (otherMol, opts = {}) {
        return appendMolecule(this, otherMol, opts);
    };
    cls.prototype.attachGroupByMarker = function (groupParsed, markerX = 'Xe', markerY = 'He', opts = {}) {
        return attachGroupByMarker(this, groupParsed, markerX, markerY, opts);
    };
    cls.prototype.attachMoleculeByMarker = function (otherMol, markerX = 'Xe', markerY = 'He', opts = {}) {
        return attachMoleculeByMarker(this, otherMol, markerX, markerY, opts);
    };
    cls.prototype.attachParsedByDirection = function (capAtom, groupParsed, params = {}) {
        return attachParsedByDirection(this, capAtom, groupParsed, params);
    };
    cls.prototype.attachMoleculeByDirection = function (capAtom, otherMol, params = {}) {
        return attachMoleculeByDirection(this, otherMol, capAtom, params);
    };
    cls.prototype.replicate = function (nrep, lvec = null) {
        return replicateMolecule(this, nrep, lvec);
    };
    // expose helpers (optional, low-level)
    cls._getParsedPos = _getParsedPos;
    cls._findMarkerPairsMol = _findMarkerPairsMol;
    cls._findMarkerPairsParsed = _findMarkerPairsParsed;
    cls._computeMarkerAttachRotation = _computeMarkerAttachRotation;
    cls._transformParsed = _transformParsed;
    cls._removeAtomsFromParsed = _removeAtomsFromParsed;
}
