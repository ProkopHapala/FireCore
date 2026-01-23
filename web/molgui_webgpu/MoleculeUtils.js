import { EditableMolecule } from "./EditableMolecule.js";
import { Vec3 } from "../common_js/Vec3.js";
import { Mat3 } from "../common_js/Mat3.js";


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

// --- Bridge utilities (moved from ScriptRunner) ---

export function collapseBridgeAt(mol, idBridge) {
    if (!mol || !mol.atoms) throw new Error('collapse_bridge_at: molecule missing');
    const iBridge = mol.getAtomIndex(idBridge);
    if (iBridge < 0) throw new Error(`collapse_bridge_at: atom id=${idBridge} not found`);
    const a = mol.atoms[iBridge];
    if (!a) throw new Error(`collapse_bridge_at: atom missing at index ${iBridge}`);
    const neighIds = [];
    const hydIds = [];
    for (const ib of a.bonds) {
        const b = mol.bonds[ib];
        if (!b) continue;
        b.ensureIndices(mol);
        const jb = b.other(iBridge);
        if (jb < 0 || jb >= mol.atoms.length) continue;
        const nb = mol.atoms[jb];
        if (!nb) continue;
        if (nb.Z === 1) hydIds.push(nb.id);
        else neighIds.push(nb.id);
    }
    if (neighIds.length !== 2) throw new Error(`collapse_bridge_at: expected 2 non-H neighbors, got ${neighIds.length}`);
    if (neighIds[0] === neighIds[1]) throw new Error('collapse_bridge_at: neighbors identical, cannot collapse');
    for (const hid of hydIds) mol.removeAtomById(hid);
    mol.removeAtomById(idBridge);
    const i0 = mol.getAtomIndex(neighIds[0]);
    const i1 = mol.getAtomIndex(neighIds[1]);
    if (i0 < 0 || i1 < 0) throw new Error('collapse_bridge_at: neighbor atom missing after removal');
    let alreadyBonded = false;
    for (const ib of mol.atoms[i0].bonds) {
        const b = mol.bonds[ib];
        if (!b) continue;
        b.ensureIndices(mol);
        const ja = b.other(i0);
        if (ja < 0 || ja >= mol.atoms.length) continue;
        if (mol.atoms[ja].id === neighIds[1]) { alreadyBonded = true; break; }
    }
    if (!alreadyBonded) mol.addBond(neighIds[0], neighIds[1]);
    return neighIds;
}

export function collapseBridgeRandom(mol) {
    if (!mol || !mol.selection) throw new Error('collapse_bridge: selection missing');
    if (mol.selection.size === 0) throw new Error('collapse_bridge: selection empty');
    const ids = Array.from(mol.selection);
    const pick = Math.floor(Math.random() * ids.length);
    return collapseBridgeAt(mol, ids[pick]);
}

export function collapseAllBridges(mol, opts = {}) {
    const collapseRound = (requireH2) => {
        const nsel = mol.selectBridgeCandidates ? mol.selectBridgeCandidates({ requireH2 }) : 0;
        let nCollapsed = 0;
        for (let i = 0; i < nsel; i++) {
            const ids = Array.from(mol.selection);
            if (ids.length === 0) break;
            const pick = Math.floor(Math.random() * ids.length);
            collapseBridgeAt(mol, ids[pick]);
            nCollapsed++;
            if (mol.selectBridgeCandidates) mol.selectBridgeCandidates({ requireH2 });
            if (mol.selection.size === 0) break;
        }
        return nCollapsed;
    };
    let total = 0;
    total += collapseRound(true);
    total += collapseRound(false);
    return total;
}

export function insertBridge(mol, aId, bId, params = {}) {
    if (!mol || !mol.atoms) throw new Error('insert_bridge: molecule missing');
    const dbg = params.dbg || (()=>{});
    const ia = mol.getAtomIndex(aId);
    const ib = mol.getAtomIndex(bId);
    if (ia < 0 || ib < 0) throw new Error('insert_bridge: atom not found');
    const pa = mol.atoms[ia].pos;
    const pb = mol.atoms[ib].pos;
    const bondId = mol._findBondId(ia, ib);
    if (bondId === null) throw new Error('insert_bridge: no bond between given atoms');
    mol.removeBondById(bondId);
    const mid = new Vec3().setAdd(pa, pb).mulScalar(0.5);
    const dir = new Vec3().setSub(pb, pa);
    const L = dir.normalize();
    if (L < 1e-12) throw new Error('insert_bridge: zero-length bond');
    const upA = mol._sumHydrogenDirs(ia);
    const upB = mol._sumHydrogenDirs(ib);
    const up = new Vec3().setAdd(upA, upB);
    let upLen = up.normalize();
    if (upLen < 1e-6) {
        up.set(Math.abs(dir.x) < 0.9 ? 1 : 0, Math.abs(dir.x) < 0.9 ? 0 : 1, 0);
        up.setCross(dir, up).normalize();
    }
    const upOffset = (params.upOffsetFactor !== undefined ? +params.upOffsetFactor : 0.5) * L;
    const cPos = new Vec3().setAdd(mid, new Vec3().setV(up).mulScalar(upOffset));
    const side = new Vec3().setCross(dir, up);
    if (side.normalize() < 1e-12) {
        side.set(Math.abs(up.x) < 0.9 ? 1 : 0, Math.abs(up.x) < 0.9 ? 0 : 1, 0);
        side.setCross(dir, side).normalize();
    }
    const hDist = (params.hDist !== undefined) ? +params.hDist : 1.09;
    const hUpOffset = (params.hUpOffsetFactor !== undefined) ? +params.hUpOffsetFactor : 0.5 * L;
    const cId = mol.addAtom(cPos.x, cPos.y, cPos.z, 6);
    if (!(cId >= 0)) throw new Error('insert_bridge: failed to add carbon');
    mol.addBond(aId, cId);
    mol.addBond(cId, bId);
    dbg(`[insert_bridge] c=${cId} dir=(${dir.x.toFixed(3)},${dir.y.toFixed(3)},${dir.z.toFixed(3)}) up=(${up.x.toFixed(3)},${up.y.toFixed(3)},${up.z.toFixed(3)}) side=(${side.x.toFixed(3)},${side.y.toFixed(3)},${side.z.toFixed(3)})`);
    if (params.addHydrogens !== false) {
        const base = new Vec3().setAdd(cPos, new Vec3().setV(up).mulScalar(hUpOffset));
        const h1p = new Vec3().setAdd(base, new Vec3().setV(side).mulScalar(hDist));
        const h2p = new Vec3().setAdd(base, new Vec3().setV(side).mulScalar(-hDist));
        const h1 = mol.addAtom(h1p.x, h1p.y, h1p.z, 1);
        const h2 = mol.addAtom(h2p.x, h2p.y, h2p.z, 1);
        if (!(h1 >= 0 && h2 >= 0)) throw new Error('insert_bridge: failed to add hydrogens');
        mol.addBond(cId, h1);
        mol.addBond(cId, h2);
        dbg(`[insert_bridge] h1=${h1} h2=${h2}`);
    }
    return cId;
}

export function insertBridgeRandom(mol, params = {}) {
    if (!mol || !mol.atoms) throw new Error('insert_bridge_random: molecule missing');
    const minHeavy = (params.minHeavy !== undefined) ? params.minHeavy | 0 : 3;
    const minHyd = (params.minHyd !== undefined) ? params.minHyd | 0 : 1;
    const candidates = [];
    for (let ib = 0; ib < mol.bonds.length; ib++) {
        const b = mol.bonds[ib];
        if (!b) continue;
        b.ensureIndices(mol);
        const ia = b.a;
        const ja = b.b;
        const a = mol.atoms[ia];
        const c = mol.atoms[ja];
        if (!a || !c) continue;
        if (a.Z !== 6 || c.Z !== 6) continue;
        const stA = mol._atomHeavyHydCounts(ia);
        const stB = mol._atomHeavyHydCounts(ja);
        const okA = (stA.heavy >= minHeavy) && (stA.hyd >= minHyd);
        const okB = (stB.heavy >= minHeavy) && (stB.hyd >= minHyd);
        if (okA && okB) candidates.push({ aId: a.id, bId: c.id });
    }
    if (candidates.length === 0) throw new Error('insert_bridge_random: no candidate bonds found');
    const pick = (params.pickIndex !== undefined) ? (params.pickIndex | 0) : Math.floor(Math.random() * candidates.length);
    if (pick < 0 || pick >= candidates.length) throw new Error(`insert_bridge_random: pickIndex out of range (0..${candidates.length - 1})`);
    const chosen = candidates[pick];
    const dbg = (msg)=>{};
    return insertBridge(mol, chosen.aId, chosen.bId, { ...params, dbg });
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
