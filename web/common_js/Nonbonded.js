import { surfaceAtomIndices } from './CollisionWorkgroups.js';
import { aabbOverlap3D, dist2ToAabb } from './Buckets.js';

export function combineREQ(mmParams, atomI, atomJ) {
    const ti = mmParams.getAtomTypeForAtom(atomI);
    const tj = mmParams.getAtomTypeForAtom(atomJ);
    if (!ti || !tj) throw new Error('combineREQ: missing atom type');
    const R0 = +ti.RvdW + +tj.RvdW;
    const E0 = Math.sqrt(Math.max(0.0, +ti.EvdW)) * Math.sqrt(Math.max(0.0, +tj.EvdW));
    const Q = (+ti.Qbase) * (+tj.Qbase);
    const H = (+ti.Hb) * (+tj.Hb);
    return { R0, E0, Q, H };
}

export function collisionPair(dp, R_min, E_min, R_cut, R_cut2) {
    const dx = dp[0], dy = dp[1], dz = dp[2];
    const r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > R_cut2 * R_cut2) return { fx: 0, fy: 0, fz: 0, E: 0, inRange: false };
    const r = Math.sqrt(r2);
    const d1 = R_cut - R_min;
    const d2 = R_cut2 - R_cut;
    const k1 = -2.0 * E_min / (d1 * (d1 + d2));
    const k2 = -k1 * (d1 / d2);
    let E, F;
    if (r < R_cut) {
        const x = r - R_min;
        E = 0.5 * k1 * x * x + E_min;
        F = -k1 * x;
    } else {
        const x = r - R_cut2;
        E = 0.5 * k2 * x * x;
        F = -k2 * x;
    }
    const fr = F / r;
    return { fx: dx * fr, fy: dy * fr, fz: dz * fr, E, inRange: true };
}

export function ljqh_pair(dp, REQ, R2damp = 1e-8) {
    const dx = dp[0], dy = dp[1], dz = dp[2];
    const r2 = dx * dx + dy * dy + dz * dz;
    const ir2_ = 1.0 / (r2 + R2damp);
    const Ec = 14.3996448915 * REQ.Q * Math.sqrt(ir2_);
    const ir2 = 1.0 / r2;
    const u2 = (REQ.R0 * REQ.R0) * ir2;
    const u6 = u2 * u2 * u2;
    const vdW = u6 * REQ.E0;
    const E = (u6 - 2.0) * vdW + Ec;
    const fr = -12.0 * (u6 - 1.0) * vdW * ir2 - Ec * ir2_;
    return { fx: dx * fr, fy: dy * fr, fz: dz * fr, E };
}

function bonded12and13Sets(mol, atomSet) {
    const nAtoms = mol.atoms.length | 0;
    const out = new Array(nAtoms);
    for (let i = 0; i < nAtoms; i++) out[i] = null;
    for (let ia = 0; ia < nAtoms; ia++) {
        if (atomSet && !atomSet.has(ia)) continue;
        const a = mol.atoms[ia];
        const s = new Set();
        const n1 = [];
        for (let k = 0; k < (a.bonds ? a.bonds.length : 0); k++) {
            const b = mol.bonds[a.bonds[k] | 0];
            if (!b) continue;
            b.ensureIndices(mol);
            const ja = (b.a === ia) ? (b.b | 0) : (b.a | 0);
            if (ja === ia) continue;
            n1.push(ja);
            if (atomSet && !atomSet.has(ja)) continue;
            s.add(ja);
        }
        for (let i1 = 0; i1 < n1.length; i1++) {
            const ja = n1[i1] | 0;
            const aj = mol.atoms[ja];
            for (let k = 0; k < (aj.bonds ? aj.bonds.length : 0); k++) {
                const b = mol.bonds[aj.bonds[k] | 0];
                if (!b) continue;
                b.ensureIndices(mol);
                const jb = (b.a === ja) ? (b.b | 0) : (b.a | 0);
                if (jb === ia) continue;
                if (atomSet && !atomSet.has(jb)) continue;
                s.add(jb);
            }
        }
        out[ia] = s;
    }
    return out;
}

export function computeNonbondBruteForceKernelStyle({ pos, mol, mmParams, surfaceOnly = true, exclude12and13 = true, R2damp = 1e-8, collectPairs = false }) {
    if (!pos || !mol || !mmParams) throw new Error('computeNonbondBruteForceKernelStyle: pos/mol/mmParams required');
    const nAtoms = mol.atoms.length | 0;
    if (pos.length !== nAtoms * 3) throw new Error(`computeNonbondBruteForceKernelStyle: pos length mismatch pos=${pos.length} nAtoms=${nAtoms}`);
    const surfSet = surfaceOnly ? new Set(surfaceAtomIndices(mol)) : null;
    const exSets = exclude12and13 ? bonded12and13Sets(mol, surfSet) : null;
    const f = new Float64Array(nAtoms * 3);
    let Etot = 0.0;
    const pairs = collectPairs ? [] : null;
    for (let ia = 0; ia < nAtoms; ia++) {
        if (surfSet && !surfSet.has(ia)) continue;
        const ex = exSets ? exSets[ia] : null;
        const i0 = ia * 3;
        const xi = pos[i0 + 0], yi = pos[i0 + 1], zi = pos[i0 + 2];
        for (let ja = 0; ja < nAtoms; ja++) {
            if (ja === ia) continue;
            if (surfSet && !surfSet.has(ja)) continue;
            if (ex && ex.has(ja)) continue;
            const j0 = ja * 3;
            const dx = pos[j0 + 0] - xi;
            const dy = pos[j0 + 1] - yi;
            const dz = pos[j0 + 2] - zi;
            const REQ = combineREQ(mmParams, mol.atoms[ia], mol.atoms[ja]);
            const { fx, fy, fz, E } = ljqh_pair([dx, dy, dz], REQ, R2damp);
            f[i0 + 0] += fx; f[i0 + 1] += fy; f[i0 + 2] += fz;
            Etot += E;
            if (pairs) pairs.push([ia, ja]);
        }
    }
    return { f, Etot, pairs };
}

export function computeNonbondByGroups({ pos, mol, mmParams, group_atoms, group_nAtoms, group_bbox_min, group_bbox_max, radius, excl_icol, EXCL_MAX = 16, groupCap = 32, icol = null, surfaceOnly = true, margin = 0.0, R2damp = 1e-8, collectPairs = false, collisionMode = false, R_cut = 0.0, R_cut2 = 0.0 }) {
    if (!pos || !mol || !mmParams) throw new Error('computeNonbondByGroups: pos/mol/mmParams required');
    const nAtoms = mol.atoms.length | 0;
    if (pos.length !== nAtoms * 3) throw new Error(`computeNonbondByGroups: pos length mismatch pos=${pos.length} nAtoms=${nAtoms}`);
    if (!group_atoms || !group_nAtoms || !group_bbox_min || !group_bbox_max) throw new Error('computeNonbondByGroups: group arrays required');
    if (!radius) throw new Error('computeNonbondByGroups: radius required');
    if (!excl_icol) throw new Error('computeNonbondByGroups: excl_icol required');
    const nGroups = group_nAtoms.length | 0;
    const f = new Float64Array(nAtoms * 3);
    let Etot = 0.0;
    const pairs = collectPairs ? [] : null;
    const icolByAtom = icol;
    if (!icolByAtom) throw new Error('computeNonbondByGroups: icol required');

    const margin2 = (+margin) * (+margin);
    for (let g = 0; g < nGroups; g++) {
        const aMin = [group_bbox_min[g * 3 + 0], group_bbox_min[g * 3 + 1], group_bbox_min[g * 3 + 2]];
        const aMax = [group_bbox_max[g * 3 + 0], group_bbox_max[g * 3 + 1], group_bbox_max[g * 3 + 2]];

        const overlapGroups = [];
        for (let h = 0; h < nGroups; h++) {
            if (!aabbOverlap3D(aMin, aMax, [group_bbox_min[h * 3 + 0], group_bbox_min[h * 3 + 1], group_bbox_min[h * 3 + 2]], [group_bbox_max[h * 3 + 0], group_bbox_max[h * 3 + 1], group_bbox_max[h * 3 + 2]], margin)) continue;
            overlapGroups.push(h | 0);
        }
        overlapGroups.sort((x, y) => x - y);

        const ghostIcols = [];
        for (let idx = 0; idx < overlapGroups.length; idx++) {
            const h = overlapGroups[idx] | 0;
            const nH = group_nAtoms[h] | 0;
            const baseH = h * (groupCap | 0);
            for (let il = 0; il < nH; il++) {
                const ja = group_atoms[baseH + il] | 0;
                if (ja < 0) continue;
                const pj = [pos[ja * 3 + 0], pos[ja * 3 + 1], pos[ja * 3 + 2]];
                if ((margin > 0) && (dist2ToAabb(pj, aMin, aMax) >= margin2)) continue;
                ghostIcols.push((h * (groupCap | 0) + il) | 0);
            }
        }

        const nG = group_nAtoms[g] | 0;
        const baseG = g * (groupCap | 0);
        for (let il = 0; il < nG; il++) {
            const ia = group_atoms[baseG + il] | 0;
            if (ia < 0) continue;
            if (surfaceOnly && (icolByAtom[ia] | 0) < 0) continue;
            const ic_i = icolByAtom[ia] | 0;
            if (ic_i < 0) continue;

            const i0 = ia * 3;
            const xi = pos[i0 + 0], yi = pos[i0 + 1], zi = pos[i0 + 2];

            let iex = ia * (EXCL_MAX | 0);
            const iex_end = iex + (EXCL_MAX | 0) - 1;
            let jex = excl_icol[iex] | 0;

            for (let jj = 0; jj < ghostIcols.length; jj++) {
                const ic_j = ghostIcols[jj] | 0;
                if (ic_j === ic_i) continue;
                const ja = group_atoms[ic_j] | 0;
                if (ja < 0) continue;

                if (jex !== -1) {
                    if ((iex < iex_end) && ((jex & 0xFFFFFF) < ic_j)) iex++;
                    jex = excl_icol[iex] | 0;
                }
                if (jex === ((0 << 24) | ic_j)) continue;

                const j0 = ja * 3;
                const dx = pos[j0 + 0] - xi;
                const dy = pos[j0 + 1] - yi;
                const dz = pos[j0 + 2] - zi;
                if (collisionMode) {
                    const REQ = combineREQ(mmParams, mol.atoms[ia], mol.atoms[ja]);
                    const R_min = REQ.R0;
                    const E_min = -REQ.E0;
                    const res = collisionPair([dx, dy, dz], R_min, E_min, R_cut > 0 ? R_cut : R_min * 1.2, R_cut2 > 0 ? R_cut2 : R_min * 1.5);
                    if (!res.inRange) continue;
                    f[i0 + 0] += res.fx; f[i0 + 1] += res.fy; f[i0 + 2] += res.fz;
                    Etot += res.E;
                    if (pairs) pairs.push([ia, ja]);
                } else {
                    const REQ = combineREQ(mmParams, mol.atoms[ia], mol.atoms[ja]);
                    const { fx, fy, fz, E } = ljqh_pair([dx, dy, dz], REQ, R2damp);
                    f[i0 + 0] += fx; f[i0 + 1] += fy; f[i0 + 2] += fz;
                    Etot += E;
                    if (pairs) pairs.push([ia, ja]);
                }
            }
        }
    }

    return { f, Etot, pairs };
}

export function maxAbsDiff(a, b) {
    if (a.length !== b.length) throw new Error(`maxAbsDiff: length mismatch ${a.length} vs ${b.length}`);
    let m = 0.0;
    for (let i = 0; i < a.length; i++) {
        const d = Math.abs(a[i] - b[i]);
        if (d > m) m = d;
    }
    return m;
}
