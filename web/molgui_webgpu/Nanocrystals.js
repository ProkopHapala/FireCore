/// Reusable nanocrystal generation module: CIF loading, sphere/plane cutting, pruning, capping, defect operators.
/// Used by CLI scripts (gen_nanocrystals.mjs) and test scripts (test_ring_detection.mjs).

import fs from 'node:fs';
import path from 'node:path';
import { spawnSync } from 'node:child_process';

import { Vec3 } from '../common_js/Vec3.js';
import { MMParams } from './MMParams.js';
import { EditableMolecule } from './EditableMolecule.js';
import * as CrystalUtils from './CrystalUtils.js';
import { installMoleculeIOMethods } from './MoleculeIO.js';
import { selectBridgeCandidates } from './MoleculeSelection.js';
import { collapseBridgeAt, collapseAllBridges, insertBridge, fuseSiH2ClashPairs, atomDist, atomAngleDeg, neighborIndices, bondedPairSet } from './MoleculeUtils.js';

installMoleculeIOMethods(EditableMolecule);

/// Default args object for programmatic use (no CLI parsing needed).
/// Pass overrides to customize. resourceDir should point to cpp/common_resources.
export function defaultArgs(resourceDir, overrides = {}) {
    const o = {
        cif: path.resolve(resourceDir, 'crystals/Si-sym.cif'),
        applySymmetry: true,
        dedupSymmetry: true,
        dedupTol: 0.1,
        samples: 0,
        nx: [3, 3], ny: [3, 3], nz: [3, 3],
        centered: true,
        planeTemplates: ['a111'],
        planeSymC: 6.0,
        planeMode: 'ang',
        planeC0: 0.0,
        planeCScale: 0.55,
        planeCJitter: 0.10,
        bonds: true,
        elementTypes: path.join(resourceDir, 'ElementTypes.dat'),
        atomTypes: path.join(resourceDir, 'AtomTypes.dat'),
        bondTypes: path.join(resourceDir, 'BondTypes.dat'),
        angleTypes: path.join(resourceDir, 'AngleTypes.dat'),
        bondFactor: 1.30,
        defaultRcut: 1.60,
        caps: 'H',
        minSiDegree: 2,
        minHeavyDegree: 2,
        pruneMaxIter: 10,
        collapse: 0,
        collapseAll: false,
        requireH2: true,
        collapseProb: 0.0,
        insertProb: 0.0,
        fuseProb: 0.0,
        fuseHClashMax: 2.0,
        fuseSiMax: 4.5,
        outwardBias: 0.35,
        resolveClashes: 1,
        capHHBonds: 0,
        capHHBondDist: 1.8,
        cutMode: 'planes',
        sphereR: 6.0,
        sphereNrep: 5,
        rcutHeavy: 0.0,
        E_SiH2: 1.0, E_SiH3: 4.0, E_bare: 10.0, E_bridge: 2.0, muH: 0.0,
        prefix: 'nc',
        seed: 0,
    };
    return Object.assign(o, overrides);
}

export function loadMMParams(args) {
    const mm = new MMParams();
    const readText = (p) => fs.readFileSync(p, 'utf8');
    mm.parseElementTypes(readText(args.elementTypes));
    mm.parseAtomTypes(readText(args.atomTypes));
    mm.parseBondTypes(readText(args.bondTypes));
    mm.parseAngleTypes(readText(args.angleTypes));
    return mm;
}

export function buildCrystalFromCIFText(cifText, opts) {
    const crystal = CrystalUtils.cifToCrystalData(cifText);
    const lvecCIF = CrystalUtils.latticeVectorsFromParams(crystal.lattice);
    let sites = crystal.sites ? crystal.sites.slice() : [];
    if (opts.applySymmetry) {
        const symOps = crystal.symOps || [];
        if (!symOps.length) throw new Error('CIF requested symmetry application but file lacks symOps');
        sites = CrystalUtils.applySymmetryOpsFracSites(sites, symOps, { tol: 1e-6 });
    }
    if (opts.dedupSymmetry) {
        sites = CrystalUtils.dedupFracSitesByTolA(sites, lvecCIF, opts.dedupTol);
    }
    return CrystalUtils.cellDataFromFracSites(lvecCIF, sites);
}

export function buildPlanesFromTemplates(lvec, templates, planeSymC, planeMode, cmin, cmax) {
    const expanded = CrystalUtils.expandPlaneTemplates(templates, planeSymC);
    const b = CrystalUtils.reciprocalLattice(lvec);
    const plsFinal = [];
    for (const p of expanded) {
        const n = new Vec3().setLincomb3(p.h || 0, b[0], p.k || 0, b[1], p.l || 0, b[2]);
        plsFinal.push({ n, cmin, cmax });
    }
    return { plsFinal, planeMode };
}

export function collapseProbabilisticBridges(mol, prob, requireH2) {
    if (!(prob > 0)) return 0;
    const nsel = selectBridgeCandidates(mol, { requireH2 });
    if (nsel <= 0) return 0;
    const ids = Array.from(mol.selection);
    let nCollapsed = 0;
    for (let i = 0; i < ids.length; i++) {
        if (Math.random() >= prob) continue;
        const id = ids[i];
        if (mol.getAtomIndex(id) < 0) continue;
        collapseBridgeAt(mol, id);
        nCollapsed++;
    }
    return nCollapsed;
}

export function randInt(a, b) {
    const i0 = a | 0;
    const i1 = b | 0;
    if (i1 < i0) throw new Error(`randInt: invalid range ${a}..${b}`);
    const u = Math.random();
    const x = i0 + Math.floor(u * (i1 - i0 + 1));
    return (x > i1) ? i1 : x;
}

export function atomNeighborsCounts(mol, ia, heavyZ = 14) {
    const a = mol.atoms[ia];
    if (!a) return { nH: 0, nHeavy: 0, nHeavySame: 0 };
    let nH = 0, nHeavy = 0, nHeavySame = 0;
    const hz = heavyZ | 0;
    for (let k = 0; k < a.bonds.length; k++) {
        const ib = a.bonds[k] | 0;
        const b = mol.bonds[ib];
        if (!b) continue;
        b.ensureIndices(mol);
        const ja = b.other(ia);
        if (ja < 0 || ja >= mol.atoms.length) continue;
        const nb = mol.atoms[ja];
        if (!nb) continue;
        const z = nb.Z | 0;
        if (z === 1) nH++;
        else {
            nHeavy++;
            if (z === hz) nHeavySame++;
        }
    }
    return { nH, nHeavy, nHeavySame, nSi: nHeavySame };
}

export function pruneUndercoordinatedHeavyIter(mol, heavyZ, minHeavyDegree, maxIter) {
    const hz = heavyZ | 0;
    const minD = minHeavyDegree | 0;
    const toRemove = [];
    const toRemoveH = [];
    let totalRemoved = 0;
    let iter = 0;
    while (iter < maxIter) {
        iter++;
        toRemove.length = 0;
        toRemoveH.length = 0;
        for (let ia = 0; ia < mol.atoms.length; ia++) {
            const a = mol.atoms[ia];
            if (!a || (a.Z | 0) !== hz) continue;
            const c = atomNeighborsCounts(mol, ia, hz);
            if (c.nHeavySame < minD) {
                toRemove.push(a.id);
                for (let k = 0; k < a.bonds.length; k++) {
                    const ib = a.bonds[k] | 0;
                    const b = mol.bonds[ib];
                    if (!b) continue;
                    b.ensureIndices(mol);
                    const ja = b.other(ia);
                    if (ja < 0 || ja >= mol.atoms.length) continue;
                    const nb = mol.atoms[ja];
                    if (nb && (nb.Z | 0) === 1) toRemoveH.push(nb.id);
                }
            }
        }
        if (toRemove.length === 0) break;
        for (let i = 0; i < toRemoveH.length; i++) mol.removeAtomById(toRemoveH[i]);
        for (let i = 0; i < toRemove.length; i++) mol.removeAtomById(toRemove[i]);
        totalRemoved += toRemove.length;
    }
    return { totalRemoved, iter };
}

export function classifySurfaceCounts(mol, heavyZ = 14) {
    const hz = heavyZ | 0;
    let nSi = 0;
    let nH = 0;
    let nSiH = 0;
    let nSiH2 = 0;
    let nSiH3 = 0;
    let nBare = 0;
    let nBridge = 0;
    for (let ia = 0; ia < mol.atoms.length; ia++) {
        const a = mol.atoms[ia];
        if (!a) continue;
        const z = a.Z | 0;
        if (z === 1) { nH++; continue; }
        if (z !== hz) continue;
        nSi++;
        const c = atomNeighborsCounts(mol, ia, hz);
        const isSurface = (c.nHeavySame < 4);
        if (!isSurface) continue;
        if (c.nH === 0) {
            if (c.nHeavySame === 2) nBridge++;
            else nBare++;
        } else if (c.nH === 1) nSiH++;
        else if (c.nH === 2) nSiH2++;
        else if (c.nH >= 3) nSiH3++;
    }
    return { nSi, nH, nSiH, nSiH2, nSiH3, nBare, nBridge };
}

export function computeEnergy(cnt, args, extra = {}) {
    const nCollapsed = extra.nCollapsed || 0;
    const nInserted = extra.nInserted || 0;
    const E = args.E_SiH2 * cnt.nSiH2 + args.E_SiH3 * cnt.nSiH3 + args.E_bare * cnt.nBare + args.E_bridge * cnt.nBridge + args.muH * cnt.nH;
    return { E, nCollapsed, nInserted };
}

export function defaultRcutHeavy(heavyZ) {
    return ((heavyZ | 0) === 14) ? 2.55 : 1.75;
}

export function getHeavyZ(cell) {
    if (!cell.basisTypes || cell.basisTypes.length === 0) throw new Error('getHeavyZ: empty basisTypes');
    return cell.basisTypes[0] | 0;
}

export function buildSphereCutMol(cell, nrep, R, centered, dedupTol) {
    const nr = centered ? (2 * (nrep | 0) + 1) : (nrep | 0);
    const origin = centered
        ? new Vec3().setLincomb3(-(nrep | 0), cell.lvec[0], -(nrep | 0), cell.lvec[1], -(nrep | 0), cell.lvec[2])
        : new Vec3(0, 0, 0);
    const mol = CrystalUtils.genReplicatedCell({
        lvec: cell.lvec,
        basisPos: cell.basisPos,
        basisTypes: cell.basisTypes,
        basisCharges: cell.basisCharges,
        nRep: [nr, nr, nr],
        origin,
        dedup: true,
        dedupTol,
    });
    const cog = new Vec3(0, 0, 0);
    for (let ia = 0; ia < mol.atoms.length; ia++) cog.add(mol.atoms[ia].pos);
    if (mol.atoms.length === 0) throw new Error('buildSphereCutMol: no atoms after replication');
    cog.mulScalar(1.0 / mol.atoms.length);
    for (let ia = 0; ia < mol.atoms.length; ia++) mol.atoms[ia].pos.sub(cog);
    const R2 = R * R;
    const removeIds = [];
    for (let ia = 0; ia < mol.atoms.length; ia++) {
        const p = mol.atoms[ia].pos;
        if (p.norm2() > R2) removeIds.push(mol.atoms[ia].id);
    }
    for (let i = 0; i < removeIds.length; i++) mol.removeAtomById(removeIds[i]);
    if (mol.atoms.length < 10) throw new Error(`buildSphereCutMol: too few atoms (${mol.atoms.length}) for R=${R} nrep=${nrep}`);
    return mol;
}

export function assignCrystalAtomTypes(mol, mm) {
    for (let ia = 0; ia < mol.atoms.length; ia++) {
        const a = mol.atoms[ia];
        const z = a.Z | 0;
        if (z === 14) mol.setAtomTypeByName(a.id, 'Si', mm);
        else if (z === 6) mol.setAtomTypeByName(a.id, 'C', mm);
        else if (z === 1) mol.setAtomTypeByName(a.id, 'H', mm);
    }
}

export function mulberry32(seed) {
    let a = seed | 0;
    return function () {
        a |= 0;
        a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

/// Generate a single nanocrystal from args + cell. Returns { mol, nCollapsed, nInserted, nCaps, nPruned, nHHBonds, cnt, enr, name }.
export function generateNanocrystal(args, cell, mm, heavyZ, iout) {
    const nx = (args.cutMode === 'sphere') ? 0 : randInt(args.nx[0], args.nx[1]);
    const ny = (args.cutMode === 'sphere') ? 0 : randInt(args.ny[0], args.ny[1]);
    const nz = (args.cutMode === 'sphere') ? 0 : randInt(args.nz[0], args.nz[1]);

    const naEff = args.centered ? (2 * nx + 1) : nx;
    const nbEff = args.centered ? (2 * ny + 1) : ny;
    const ncEff = args.centered ? (2 * nz + 1) : nz;
    const La = cell.lvec[0].norm() * (naEff || 1);
    const Lb = cell.lvec[1].norm() * (nbEff || 1);
    const Lc = cell.lvec[2].norm() * (ncEff || 1);
    const Lavg = (La + Lb + Lc) / 3.0;

    const cBase = args.planeC0 + args.planeCScale * Lavg;
    const dj = args.planeCJitter * cBase;
    const jitter = (args.cutMode === 'sphere') ? 0 : dj * (2.0 * Math.random() - 1.0);
    const cmin = -cBase + jitter;
    const cmax = +cBase + jitter;

    const { plsFinal, planeMode } = (args.cutMode === 'planes')
        ? buildPlanesFromTemplates(cell.lvec, args.planeTemplates, args.planeSymC, args.planeMode, cmin, cmax)
        : { plsFinal: null, planeMode: args.planeMode };

    const mol = (args.cutMode === 'sphere')
        ? buildSphereCutMol(cell, args.sphereNrep, args.sphereR, true, args.dedupTol)
        : CrystalUtils.genReplicatedCellCutPlanes({
            lvec: cell.lvec,
            basisPos: cell.basisPos,
            basisTypes: cell.basisTypes,
            basisCharges: cell.basisCharges,
            nRep: [nx, ny, nz],
            origin: new Vec3(0, 0, 0),
            planes: plsFinal,
            planeMode,
            centered: args.centered,
            dedup: true,
            dedupTol: args.dedupTol,
        });

    if (!mm) throw new Error('generateNanocrystal: mmParams missing');

    mol.recalculateBonds(mm, { defaultRcut: args.defaultRcut, bondFactor: args.bondFactor });

    const pr = pruneUndercoordinatedHeavyIter(mol, heavyZ, args.minHeavyDegree, args.pruneMaxIter);
    if (pr.totalRemoved > 0) mol.recalculateBonds(mm, { defaultRcut: args.defaultRcut, bondFactor: args.bondFactor });

    assignCrystalAtomTypes(mol, mm);

    let nCaps = 0;
    if (args.caps && args.caps !== '0' && args.caps !== 'none') {
        const r = mol.addCappingAtoms(mm, args.caps, { onlySelection: false, bBond: true, bondFactor: args.bondFactor, outwardBias: args.outwardBias, resolveClashes: !!args.resolveClashes });
        nCaps = r.nAdded | 0;
        if (nCaps > 0) mol.recalculateBonds(mm, { defaultRcut: args.defaultRcut, bondFactor: args.bondFactor });
    }

    const pr2 = pruneUndercoordinatedHeavyIter(mol, heavyZ, args.minHeavyDegree, args.pruneMaxIter);
    if (pr2.totalRemoved > 0) mol.recalculateBonds(mm, { defaultRcut: args.defaultRcut, bondFactor: args.bondFactor });

    let nCollapsed = 0;
    let nInserted = 0;
    let nFused = 0;
    const surfaceFilter = (m, ia) => {
        const c = atomNeighborsCounts(m, ia, heavyZ);
        return c.nHeavySame < 4;
    };
    if (args.collapseAll) {
        nCollapsed = collapseAllBridges(mol);
    } else {
        if (args.insertProb > 0) {
            const candidates = [];
            for (const b of mol.bonds) {
                if (!b) continue; b.ensureIndices(mol);
                const ia = b.a, ja = b.b;
                const a = mol.atoms[ia], c = mol.atoms[ja];
                if (!a || !c) continue;
                if ((a.Z | 0) !== heavyZ || (c.Z | 0) !== heavyZ) continue;
                if (!surfaceFilter(mol, ia) || !surfaceFilter(mol, ja)) continue;
                candidates.push([a.id, c.id]);
            }

            for (const [aId, bId] of candidates) {
                if (Math.random() >= args.insertProb) continue;
                const cid = insertBridge(mol, aId, bId, { addHydrogens: true, hDist: 1.3, upOffsetFactor: 0.5 });
                const ic = mol.getAtomIndex(cid);
                if (ic >= 0) mol.atoms[ic].Z = heavyZ;
                nInserted++;
            }
        }

        if (args.collapseProb > 0) {
            mol.selection.clear();
            selectBridgeCandidates(mol, { z: heavyZ, minHeavy: 2, minHyd: 2, requireH2: true, surfaceFilter });
            for (const id of Array.from(mol.selection)) {
                if (Math.random() < args.collapseProb) {
                    collapseBridgeAt(mol, id);
                    nCollapsed++;
                }
            }
        }
    }

    if (args.fuseProb > 0) {
        nFused = fuseSiH2ClashPairs(mol, {
            heavyZ, surfaceFilter, prob: args.fuseProb,
            hClashMax: args.fuseHClashMax, siMax: args.fuseSiMax,
        });
    }

    let nHHBonds = 0;
    if (args.capHHBonds) {
        const hr = mol.addCapHHBonds(args.capHHBondDist);
        nHHBonds = hr.nAdded | 0;
    }

    const cnt = classifySurfaceCounts(mol, heavyZ);
    const enr = computeEnergy(cnt, args, { nCollapsed, nInserted });
    const nPruned = (pr.totalRemoved + pr2.totalRemoved) | 0;

    const cutTag = (args.cutMode === 'sphere')
        ? `sphereR${args.sphereR.toFixed(1)}_nrep${args.sphereNrep}`
        : `nx${nx}_ny${ny}_nz${nz}_c${cBase.toFixed(3)}_tpl${args.planeTemplates.join('-')}`;
    const name = `${args.prefix}_i${String(iout).padStart(4, '0')}_${cutTag}_pr${nPruned}_cap${nCaps}_col${nCollapsed}_ins${nInserted}_fus${nFused}`;

    return { mol, nCollapsed, nInserted, nFused, nCaps, nPruned, nHHBonds, cnt, enr, name, nx, ny, nz };
}

// ============================================================================
// NC geometry audit (moved from test_nanocrystal_geometry.mjs)
// ============================================================================

export function isBulkSi(mol, ia) {
    const a = mol.atoms[ia];
    if ((a.Z | 0) !== 14) return false;
    let nSi = 0;
    for (const jb of neighborIndices(mol, ia)) if ((mol.atoms[jb].Z | 0) === 14) nSi++;
    return nSi === 4;
}

export function auditGeometry(mol, mm, opts = {}) {
    const siHLo = opts.siHLo ?? 1.40, siHHi = opts.siHHi ?? 1.55;
    const angLo = opts.angLo ?? 100, angHi = opts.angHi ?? 120;
    const clashCut = opts.clashCut ?? 1.8;
    const bonded = bondedPairSet(mol);
    const report = { natoms: mol.atoms.length, nbonds: mol.bonds.length, siH: { n: 0, min: null, max: null, bad: [] }, angles: { n: 0, bad: [], surfaceBad: [] }, clashes: [], valence: { bad: [] }, summary: {} };
    for (let ia = 0; ia < mol.atoms.length; ia++) {
        const a = mol.atoms[ia];
        const nbs = neighborIndices(mol, ia);
        const z = a.Z | 0;
        if (z === 14) {
            const nSi = nbs.filter(jb => (mol.atoms[jb].Z | 0) === 14).length;
            const nH = nbs.filter(jb => (mol.atoms[jb].Z | 0) === 1).length;
            const expectedH = Math.max(0, 4 - nSi);
            if (nH !== expectedH) report.valence.bad.push({ id: a.id, ia, nSi, nH, expectedH });
            for (const jb of nbs) {
                if ((mol.atoms[jb].Z | 0) !== 1) continue;
                const d = atomDist(a, mol.atoms[jb]);
                report.siH.n++;
                if (report.siH.min === null || d < report.siH.min) report.siH.min = d;
                if (report.siH.max === null || d > report.siH.max) report.siH.max = d;
                if (d < siHLo || d > siHHi) report.siH.bad.push({ siId: a.id, hId: mol.atoms[jb].id, d });
            }
            if (nbs.length >= 3) {
                const bulk = isBulkSi(mol, ia);
                for (let i = 0; i < nbs.length; i++) for (let j = i + 1; j < nbs.length; j++) {
                    const ang = atomAngleDeg(a, mol.atoms[nbs[i]], mol.atoms[nbs[j]]);
                    report.angles.n++;
                    if (!isFinite(ang) || ang < angLo || ang > angHi) {
                        const rec = { siId: a.id, ia, ang, neighbors: [mol.atoms[nbs[i]].Z, mol.atoms[nbs[j]].Z], bulk };
                        if (bulk) report.angles.bad.push(rec); else report.angles.surfaceBad.push(rec);
                    }
                }
            }
        }
    }
    for (let i = 0; i < mol.atoms.length; i++) {
        for (let j = i + 1; j < mol.atoms.length; j++) {
            const lo = Math.min(i, j), hi = Math.max(i, j);
            if (bonded.has(`${lo}-${hi}`)) continue;
            const d = atomDist(mol.atoms[i], mol.atoms[j]);
            if (d < clashCut) {
                const zi = mol.atoms[i].Z | 0, zj = mol.atoms[j].Z | 0;
                if (zi === 1 || zj === 1) report.clashes.push({ i, j, zi, zj, d });
            }
        }
    }
    report.summary = {
        siH_ok: report.siH.bad.length === 0,
        bulk_angles_ok: report.angles.bad.length === 0,
        clashes_ok: report.clashes.length === 0,
        valence_ok: report.valence.bad.length === 0,
        pass: report.siH.bad.length === 0 && report.angles.bad.length === 0 && report.clashes.length === 0 && report.valence.bad.length === 0,
        n_surface_angle_exceptions: report.angles.surfaceBad.length,
    };
    return report;
}

export function oneLineSummary(name, r) {
    const s = r.summary;
    return `[${name}] atoms=${r.natoms} Si-H=[${r.siH.min?.toFixed(3) ?? '?'},${r.siH.max?.toFixed(3) ?? '?'}] clashes=${r.clashes.length} bulkAngBad=${r.angles.bad.length} surfAngExc=${r.angles.surfaceBad.length} valenceBad=${r.valence.bad.length} PASS=${s.pass}`;
}

export function buildSiH4() {
    const mol = new EditableMolecule();
    const si = mol.addAtom(0, 0, 0, 14);
    const r = 1.48;
    const s = 1 / Math.sqrt(3);
    const dirs = [[s, s, s], [s, -s, -s], [-s, s, -s], [-s, -s, s]];
    for (const [x, y, z] of dirs) {
        const h = mol.addAtom(x * r, y * r, z * r, 1);
        mol.addBond(si, h, 1, 1);
    }
    return mol;
}

// ============================================================================
// Ensemble pipeline utilities (moved from run_nanocrystal_ensemble.mjs)
// ============================================================================

export function sampleValue(spec, rnd) {
    if (spec == null) return null;
    if (typeof spec === 'number' || typeof spec === 'string' || typeof spec === 'boolean') return spec;
    if (spec.uniform) return spec.uniform[0] + rnd() * (spec.uniform[1] - spec.uniform[0]);
    if (spec.uniform_int) return (spec.uniform_int[0] + Math.floor(rnd() * (spec.uniform_int[1] - spec.uniform_int[0] + 1))) | 0;
    if (spec.choice) return spec.choice[Math.floor(rnd() * spec.choice.length)];
    throw new Error(`sampleValue: unknown spec ${JSON.stringify(spec)}`);
}

export function sampleCutSpec(cutMixture, rnd) {
    const weights = cutMixture.map(c => +c.weight);
    const sum = weights.reduce((a, b) => a + b, 0);
    let u = rnd() * sum;
    let idx = 0;
    for (let i = 0; i < weights.length; i++) {
        u -= weights[i];
        if (u <= 0) { idx = i; break; }
        idx = i;
    }
    const comp = { ...cutMixture[idx] };
    delete comp.weight;
    const out = { cutMode: comp.cutMode };
    for (const [k, v] of Object.entries(comp)) {
        if (k === 'cutMode') continue;
        out[k] = sampleValue(v, rnd);
    }
    return out;
}

export function crystalId(seed, index) {
    return `${String(seed >>> 0).padStart(8, '0')}_${String(index).padStart(6, '0')}`;
}

export function appendManifest(manifestPath, rec) {
    fs.appendFileSync(manifestPath, JSON.stringify(rec) + '\n');
}

export function pyEnv(repoRoot) {
    return { ...process.env, PYTHONPATH: [repoRoot, process.env.PYTHONPATH].filter(Boolean).join(path.delimiter) };
}

export function runCmd(cmd, args, { cwd, label = '', statusPath = null, python = 'python3', repoRoot } = {}) {
    const t0 = performance.now();
    const r = spawnSync(cmd, args, { cwd, encoding: 'utf8', stdio: ['ignore', 'ignore', 'pipe'], env: pyEnv(repoRoot || cwd) });
    const timing_ms = performance.now() - t0;
    if (r.status !== 0) {
        const err = new Error(`[${label}] failed exit=${r.status}\n${r.stderr}`);
        err.timing_ms = timing_ms;
        throw err;
    }
    let info = { status: 'ok', timing_ms };
    if (statusPath && fs.existsSync(statusPath)) {
        info = { ...JSON.parse(fs.readFileSync(statusPath, 'utf8')), timing_ms };
    }
    return { stderr: r.stderr, timing_ms, info };
}

export function runPy(subcmd, pyArgs, opts) {
    const py = opts.python || 'python3';
    return runCmd(py, ['-m', 'pyBall.nanocrystal_pipeline', subcmd, ...pyArgs], opts);
}

export function writeTimingReport(outputDir, stageTimings) {
    const stages = {};
    for (const [name, arr] of Object.entries(stageTimings)) {
        if (!arr.length) continue;
        arr.sort((a, b) => a - b);
        const sum = arr.reduce((a, b) => a + b, 0);
        stages[name] = { min_ms: arr[0], median_ms: arr[Math.floor(arr.length / 2)], max_ms: arr[arr.length - 1], total_ms: sum, n: arr.length };
    }
    fs.mkdirSync(outputDir, { recursive: true });
    fs.writeFileSync(path.join(outputDir, 'timing_report.json'), JSON.stringify({ stages, notes: 'Measured wall times; hardware-specific.' }, null, 2));
    const lines = ['# Nanocrystal ensemble timing report', ''];
    for (const [k, v] of Object.entries(stages)) {
        lines.push(`## ${k}`);
        lines.push(`- n=${v.n} min=${v.min_ms.toFixed(1)} ms median=${v.median_ms.toFixed(1)} ms max=${v.max_ms.toFixed(1)} ms total=${v.total_ms.toFixed(1)} ms`);
        lines.push('');
    }
    fs.writeFileSync(path.join(outputDir, 'timing_report.md'), lines.join('\n'));
}

/// Convert ensemble config format into defaultArgs for generateNanocrystal().
export function genArgsFromConfig(config, resourceDir, overrides = {}) {
    const cut = config.cut || {};
    const rep = config.replication || {};
    const pass = config.passivation || {};
    const def = config.defects || {};
    const o = {
        cif: config.cif || path.join(resourceDir, 'crystals/Si-sym.cif'),
        applySymmetry: config.applySymmetry ?? 0,
        cutMode: cut.cutMode || 'sphere',
        caps: pass.caps || 'H',
        resolveClashes: pass.resolveClashes ?? false,
        outwardBias: pass.outwardBias ?? 0,
        insertProb: def.insertProb ?? 0,
        collapseProb: def.collapseProb ?? 0,
        bondFactor: 1.30,
        defaultRcut: 1.60,
        minHeavyDegree: 2,
        pruneMaxIter: 3,
        samples: 1,
        seed: 42,
        maxFiles: 1,
    };
    if (pass.capHHBonds) { o.capHHBonds = true; o.capHHBondDist = pass.capHHBondDist ?? 1.8; }
    if (cut.cutMode === 'sphere') {
        o.sphereR = cut.sphereR ?? 5.0;
        o.sphereNrep = cut.sphereNrep ?? 5;
        o.rcutHeavy = cut.rcutHeavy ?? 0;
    } else if (cut.cutMode === 'planes') {
        o.nx = [rep.nx ?? 2, rep.nx ?? 2];
        o.ny = [rep.ny ?? 2, rep.ny ?? 2];
        o.nz = [rep.nz ?? 2, rep.nz ?? 2];
        o.centered = 1;
        o.planeTemplates = cut.planeTemplates || ['a111'];
        o.planeSymC = cut.planeSymC ?? 6;
        o.planeMode = cut.planeMode || 'ang';
        o.planeC0 = cut.planeC0 ?? 0;
        o.planeCScale = cut.planeCScale ?? 0.45;
        o.planeCJitter = cut.planeCJitter ?? 0;
    }
    return defaultArgs(resourceDir, { ...o, ...overrides });
}

export function genParamsFromAtlasEntry(atlasBase, entry) {
    return {
        cif: atlasBase.cif,
        applySymmetry: atlasBase.applySymmetry ?? 0,
        cut: entry.cut,
        replication: entry.replication || { nx: 2, ny: 2, nz: 2 },
        passivation: atlasBase.passivation,
        defects: atlasBase.defects || { insertProb: 0, collapseProb: 0 },
    };
}
