export function surfaceAtomIndices(mol) {
    if (!mol || !mol.atoms) throw new Error('surfaceAtomIndices: mol required');
    const out = [];
    for (let ia = 0; ia < mol.atoms.length; ia++) {
        const a = mol.atoms[ia];
        const z = a.Z | 0;
        if (z === 1) { out.push(ia); continue; }
        if (z <= 1) continue;
        let nHeavy = 0;
        for (let k = 0; k < (a.bonds ? a.bonds.length : 0); k++) {
            const ib = a.bonds[k] | 0;
            const b = mol.bonds[ib];
            if (!b) continue;
            b.ensureIndices(mol);
            const ja = (b.a === ia) ? (b.b | 0) : (b.a | 0);
            if ((mol.atoms[ja].Z | 0) > 1) nHeavy++;
        }
        if (nHeavy < 4) out.push(ia);
    }
    return out;
}

function dist2(pos, ia, ja) {
    const dx = pos[ia * 3 + 0] - pos[ja * 3 + 0];
    const dy = pos[ia * 3 + 1] - pos[ja * 3 + 1];
    const dz = pos[ia * 3 + 2] - pos[ja * 3 + 2];
    return dx * dx + dy * dy + dz * dz;
}

function chooseFarthestPivot(pos, candidates, pivots) {
    let best = -1;
    let bestMinD2 = -1;
    for (let ii = 0; ii < candidates.length; ii++) {
        const ia = candidates[ii] | 0;
        let minD2 = Infinity;
        for (let k = 0; k < pivots.length; k++) {
            const d2 = dist2(pos, ia, pivots[k] | 0);
            if (d2 < minD2) minD2 = d2;
        }
        if (minD2 > bestMinD2) { bestMinD2 = minD2; best = ia; }
    }
    return best;
}

function assignToPivots(pos, surf, pivots) {
    const nG = pivots.length | 0;
    const gidOf = new Int32Array(surf.length);
    const counts = new Int32Array(nG);
    for (let ii = 0; ii < surf.length; ii++) {
        const ia = surf[ii] | 0;
        let bestG = 0;
        let bestD2 = dist2(pos, ia, pivots[0] | 0);
        for (let g = 1; g < nG; g++) {
            const d2 = dist2(pos, ia, pivots[g] | 0);
            if (d2 < bestD2) { bestD2 = d2; bestG = g; }
        }
        gidOf[ii] = bestG;
        counts[bestG]++;
    }
    return { gidOf, counts };
}

export function buildCollisionWorkgroups({ pos, mol, groupCap = 32, fillFactor = 0.8, maxSplitIter = 64 }) {
    if (!pos || !(pos.length > 0)) throw new Error('buildCollisionWorkgroups: pos required');
    if (!mol || !mol.atoms) throw new Error('buildCollisionWorkgroups: mol required');
    const nAtoms = mol.atoms.length | 0;
    if (pos.length !== nAtoms * 3) throw new Error(`buildCollisionWorkgroups: pos length mismatch pos=${pos.length} nAtoms=${nAtoms}`);
    const surf = surfaceAtomIndices(mol);
    const nSurf = surf.length | 0;
    if (nSurf === 0) throw new Error('buildCollisionWorkgroups: no surface atoms');
    const tgt = Math.max(1, Math.floor((groupCap | 0) * (+fillFactor)));
    let nGroups = Math.ceil(nSurf / tgt);
    if (nGroups < 1) nGroups = 1;

    const pivots = [surf[0] | 0];
    while (pivots.length < nGroups) {
        const p = chooseFarthestPivot(pos, surf, pivots);
        if (p < 0) throw new Error('buildCollisionWorkgroups: failed to choose pivot');
        pivots.push(p | 0);
    }

    let assign = assignToPivots(pos, surf, pivots);
    let splitIter = 0;
    while (true) {
        let gMax = -1;
        let cMax = 0;
        for (let g = 0; g < assign.counts.length; g++) {
            const c = assign.counts[g] | 0;
            if (c > cMax) { cMax = c; gMax = g; }
        }
        if (cMax <= groupCap) break;
        if (splitIter >= maxSplitIter) throw new Error(`buildCollisionWorkgroups: capacity split did not converge (cMax=${cMax} groupCap=${groupCap})`);
        const members = [];
        for (let ii = 0; ii < surf.length; ii++) if ((assign.gidOf[ii] | 0) === gMax) members.push(surf[ii] | 0);
        const pNew = chooseFarthestPivot(pos, members, pivots);
        if (pNew < 0) throw new Error('buildCollisionWorkgroups: failed to split pivot');
        pivots.push(pNew | 0);
        assign = assignToPivots(pos, surf, pivots);
        splitIter++;
    }

    nGroups = pivots.length | 0;
    const group_nAtoms = new Int32Array(nGroups);
    for (let g = 0; g < nGroups; g++) group_nAtoms[g] = assign.counts[g] | 0;

    const group_atoms = new Int32Array(nGroups * groupCap);
    group_atoms.fill(-1);
    const writePtr = new Int32Array(nGroups);
    for (let ii = 0; ii < surf.length; ii++) {
        const ia = surf[ii] | 0;
        const g = assign.gidOf[ii] | 0;
        const k = writePtr[g] | 0;
        if (k >= groupCap) throw new Error(`buildCollisionWorkgroups: group overflow g=${g} k=${k} groupCap=${groupCap}`);
        group_atoms[g * groupCap + k] = ia;
        writePtr[g] = (k + 1) | 0;
    }

    const icol = new Int32Array(nAtoms);
    const icolGroup = new Int32Array(nAtoms);
    icol.fill(-1);
    icolGroup.fill(-1);
    for (let g = 0; g < nGroups; g++) {
        for (let il = 0; il < groupCap; il++) {
            const ia = group_atoms[g * groupCap + il] | 0;
            if (ia < 0) continue;
            icol[ia] = (g * groupCap + il) | 0;
            icolGroup[ia] = g | 0;
        }
    }

    return { nAtoms, surface: surf, pivots, groupCap: groupCap | 0, nGroups, group_atoms, group_nAtoms, icol, icolGroup };
}

export function computeGroupAABBs({ pos, radius, group_atoms, group_nAtoms, groupCap }) {
    if (!pos || !radius) throw new Error('computeGroupAABBs: pos and radius required');
    const nGroups = (group_nAtoms.length | 0);
    const bbox_min = new Float64Array(nGroups * 3);
    const bbox_max = new Float64Array(nGroups * 3);
    for (let g = 0; g < nGroups; g++) {
        let xmin = Infinity, ymin = Infinity, zmin = Infinity;
        let xmax = -Infinity, ymax = -Infinity, zmax = -Infinity;
        const base = g * (groupCap | 0);
        const n = group_nAtoms[g] | 0;
        for (let il = 0; il < n; il++) {
            const ia = group_atoms[base + il] | 0;
            if (ia < 0) continue;
            const r = +radius[ia];
            const x = +pos[ia * 3 + 0], y = +pos[ia * 3 + 1], z = +pos[ia * 3 + 2];
            xmin = Math.min(xmin, x - r); ymin = Math.min(ymin, y - r); zmin = Math.min(zmin, z - r);
            xmax = Math.max(xmax, x + r); ymax = Math.max(ymax, y + r); zmax = Math.max(zmax, z + r);
        }
        if (!(xmin < Infinity)) throw new Error(`computeGroupAABBs: empty group g=${g}`);
        bbox_min[g * 3 + 0] = xmin; bbox_min[g * 3 + 1] = ymin; bbox_min[g * 3 + 2] = zmin;
        bbox_max[g * 3 + 0] = xmax; bbox_max[g * 3 + 1] = ymax; bbox_max[g * 3 + 2] = zmax;
    }
    return { bbox_min, bbox_max };
}

export function buildExclIcol_1_2_3({ mol, icol, EXCL_MAX = 16, ipbc = 0 }) {
    if (!mol || !mol.atoms) throw new Error('buildExclIcol_1_2_3: mol required');
    if (!icol) throw new Error('buildExclIcol_1_2_3: icol required');
    const nAtoms = mol.atoms.length | 0;
    const excl = new Int32Array(nAtoms * (EXCL_MAX | 0));
    const excl_count = new Int32Array(nAtoms);
    excl.fill(-1);

    for (let ia = 0; ia < nAtoms; ia++) {
        if ((icol[ia] | 0) < 0) continue;
        const m = new Map();
        const add = (ja) => {
            const j = ja | 0;
            if (j === ia) return;
            const icj = icol[j] | 0;
            if (icj < 0) return;
            if (!m.has(icj)) m.set(icj, true);
        };
        const a = mol.atoms[ia];
        for (let k = 0; k < (a.bonds ? a.bonds.length : 0); k++) {
            const ib = a.bonds[k] | 0;
            const b = mol.bonds[ib];
            if (!b) continue;
            b.ensureIndices(mol);
            const ja = (b.a === ia) ? (b.b | 0) : (b.a | 0);
            add(ja);
        }
        for (let k = 0; k < (a.bonds ? a.bonds.length : 0); k++) {
            const ib = a.bonds[k] | 0;
            const b = mol.bonds[ib];
            if (!b) continue;
            b.ensureIndices(mol);
            const ja = (b.a === ia) ? (b.b | 0) : (b.a | 0);
            const aj = mol.atoms[ja];
            for (let kk = 0; kk < (aj.bonds ? aj.bonds.length : 0); kk++) {
                const ib2 = aj.bonds[kk] | 0;
                const b2 = mol.bonds[ib2];
                if (!b2) continue;
                b2.ensureIndices(mol);
                const jb = (b2.a === ja) ? (b2.b | 0) : (b2.a | 0);
                add(jb);
            }
        }
        const keys = Array.from(m.keys());
        keys.sort((a, b) => (a | 0) - (b | 0));
        if (keys.length > (EXCL_MAX | 0)) throw new Error(`buildExclIcol_1_2_3: overflow ia=${ia} n=${keys.length} EXCL_MAX=${EXCL_MAX}`);
        const base = ia * (EXCL_MAX | 0);
        for (let i = 0; i < keys.length; i++) excl[base + i] = ((ipbc << 24) | (keys[i] | 0)) | 0;
        excl_count[ia] = keys.length | 0;
    }

    return { excl, excl_count, EXCL_MAX: EXCL_MAX | 0 };
}
