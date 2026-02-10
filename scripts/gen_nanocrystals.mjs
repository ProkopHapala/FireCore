import fs from 'node:fs';
import path from 'node:path';

import { Vec3 } from '../web/common_js/Vec3.js';

import { MMParams } from '../web/molgui_webgpu/MMParams.js';
import { EditableMolecule } from '../web/molgui_webgpu/EditableMolecule.js';
import * as CrystalUtils from '../web/molgui_webgpu/CrystalUtils.js';
import { toMol2String, toXYZString, installMoleculeIOMethods } from '../web/molgui_webgpu/MoleculeIO.js';
import { selectBridgeCandidates } from '../web/molgui_webgpu/MoleculeSelection.js';
import { collapseBridgeAt, collapseAllBridges, insertBridge } from '../web/molgui_webgpu/MoleculeUtils.js';

installMoleculeIOMethods(EditableMolecule);

function parseProb(s, name) {
    const v = +String(s).trim();
    if (!isFinite(v) || v < 0 || v > 1) throw new Error(`${name}: expected probability in [0,1], got '${s}'`);
    return v;
}

function parseCSVRange(s) {
    const t = String(s || '').trim();
    if (!t) throw new Error('parseCSVRange: empty');
    const parts = t.split(',').map(x => x.trim()).filter(x => x.length > 0);
    if (parts.length === 1) {
        const v = parseInt(parts[0], 10);
        if (!isFinite(v)) throw new Error(`parseCSVRange: bad int '${s}'`);
        return [v | 0, v | 0];
    }
    if (parts.length === 2) {
        const a = parseInt(parts[0], 10);
        const b = parseInt(parts[1], 10);
        if (!isFinite(a) || !isFinite(b)) throw new Error(`parseCSVRange: bad range '${s}'`);
        return [a | 0, b | 0];
    }
    throw new Error(`parseCSVRange: expected 'n' or 'a,b', got '${s}'`);
}

function rangeInclusive(a, b) {
    const out = [];
    const i0 = a | 0;
    const i1 = b | 0;
    if (i1 < i0) throw new Error(`rangeInclusive: invalid range ${a}..${b}`);
    for (let i = i0; i <= i1; i++) out.push(i);
    return out;
}

function parsePlaneTemplates(s) {
    const t = String(s || '').trim();
    if (!t) return [];
    return t.split(',').map(x => x.trim()).filter(x => x.length > 0);
}

function parsePlaneSymC(s) {
    const v = +String(s || '').trim();
    if (!isFinite(v) || !(v > 0)) throw new Error(`bad --planeSymC '${s}'`);
    return v;
}

function parseArgs(argv) {
    const repoRoot = path.resolve(path.dirname(new URL(import.meta.url).pathname), '..');
    const defaultRes = path.resolve(repoRoot, 'cpp/common_resources');

    const out = {
        cif: path.resolve(defaultRes, 'crystals/Si-sym.cif'),
        applySymmetry: true,
        dedupSymmetry: true,
        dedupTol: 0.1,

        samples: 0,
        nx: [3, 3],
        ny: [3, 3],
        nz: [3, 3],
        centered: true,

        planeTemplates: ['a111'],
        planeSymC: 6.0,
        planeMode: 'ang',

        planeC0: 0.0,
        planeCScale: 0.55,
        planeCJitter: 0.10,

        bonds: true,
        elementTypes: path.join(defaultRes, 'ElementTypes.dat'),
        atomTypes: path.join(defaultRes, 'AtomTypes.dat'),
        bondTypes: path.join(defaultRes, 'BondTypes.dat'),
        angleTypes: path.join(defaultRes, 'AngleTypes.dat'),
        bondFactor: 1.30,
        defaultRcut: 1.60,

        caps: 'H',
        minSiDegree: 2,
        pruneMaxIter: 10,

        collapse: 0,
        collapseAll: false,
        requireH2: true,

        collapseProb: 0.10,
        insertProb: 0.0,

        E_SiH2: 1.0,
        E_SiH3: 4.0,
        E_bare: 10.0,
        E_bridge: 2.0,
        muH: 0.0,

        outDir: path.resolve(repoRoot, 'OUT_nanocrystals'),
        prefix: 'si_nc',
        maxFiles: 0,
        seed: 0,
        statsCsv: null,
        stackedXyzOut: null,
    };

    for (let i = 2; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => {
            if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`);
            i++;
            return String(argv[i]);
        };

        if (a === '--cif') out.cif = path.resolve(nxt());
        else if (a === '--applySymmetry') out.applySymmetry = (nxt() !== '0');
        else if (a === '--dedupSymmetry') out.dedupSymmetry = (nxt() !== '0');
        else if (a === '--dedupTol') out.dedupTol = +nxt();

        else if (a === '--samples') out.samples = parseInt(nxt(), 10) | 0;
        else if (a === '--nx' || a === '--nx-range') out.nx = parseCSVRange(nxt());
        else if (a === '--ny' || a === '--ny-range') out.ny = parseCSVRange(nxt());
        else if (a === '--nz' || a === '--nz-range') out.nz = parseCSVRange(nxt());
        else if (a === '--centered') out.centered = (nxt() !== '0');

        else if (a === '--planeTemplates') out.planeTemplates = parsePlaneTemplates(nxt());
        else if (a === '--planeSymC') out.planeSymC = parsePlaneSymC(nxt());
        else if (a === '--planeMode') out.planeMode = nxt();

        else if (a === '--planeC0') out.planeC0 = +nxt();
        else if (a === '--planeCScale') out.planeCScale = +nxt();
        else if (a === '--planeCJitter') out.planeCJitter = +nxt();

        else if (a === '--bonds') out.bonds = (nxt() !== '0');
        else if (a === '--elementTypes') out.elementTypes = path.resolve(nxt());
        else if (a === '--atomTypes') out.atomTypes = path.resolve(nxt());
        else if (a === '--bondTypes') out.bondTypes = path.resolve(nxt());
        else if (a === '--angleTypes') out.angleTypes = path.resolve(nxt());
        else if (a === '--bondFactor') out.bondFactor = +nxt();
        else if (a === '--defaultRcut') out.defaultRcut = +nxt();

        else if (a === '--caps') out.caps = nxt();
        else if (a === '--minSiDegree') out.minSiDegree = parseInt(nxt(), 10) | 0;
        else if (a === '--pruneMaxIter') out.pruneMaxIter = parseInt(nxt(), 10) | 0;

        else if (a === '--collapse') out.collapse = parseInt(nxt(), 10) | 0;
        else if (a === '--collapseAll') out.collapseAll = (nxt() !== '0');
        else if (a === '--requireH2') out.requireH2 = (nxt() !== '0');

        else if (a === '--collapseProb') out.collapseProb = parseProb(nxt(), '--collapseProb');
        else if (a === '--insertProb') out.insertProb = parseProb(nxt(), '--insertProb');

        else if (a === '--E_SiH2') out.E_SiH2 = +nxt();
        else if (a === '--E_SiH3') out.E_SiH3 = +nxt();
        else if (a === '--E_bare') out.E_bare = +nxt();
        else if (a === '--E_bridge') out.E_bridge = +nxt();
        else if (a === '--muH') out.muH = +nxt();

        else if (a === '--outDir') out.outDir = path.resolve(nxt());
        else if (a === '--prefix') out.prefix = nxt();
        else if (a === '--maxFiles') out.maxFiles = parseInt(nxt(), 10) | 0;
        else if (a === '--seed') out.seed = parseInt(nxt(), 10) | 0;

        else if (a === '--statsCsv') out.statsCsv = path.resolve(nxt());
        else if (a === '--stackedXyzOut') out.stackedXyzOut = path.resolve(nxt());

        else throw new Error(`Unknown arg ${a}`);
    }

    if (out.planeMode !== 'ang' && out.planeMode !== 'frac') throw new Error("--planeMode must be 'ang' or 'frac'");
    if (!Array.isArray(out.planeTemplates) || out.planeTemplates.length === 0) throw new Error('--planeTemplates must be non-empty (e.g. a111 or a100,a110,a111)');
    if (!(out.dedupTol > 0)) throw new Error('--dedupTol must be >0');
    if (out.collapse < 0) throw new Error('--collapse must be >=0');
    if (out.samples < 0) throw new Error('--samples must be >=0');
    if (out.minSiDegree < 0) throw new Error('--minSiDegree must be >=0');
    if (out.pruneMaxIter < 0) throw new Error('--pruneMaxIter must be >=0');
    if (!(isFinite(out.planeCScale) && out.planeCScale >= 0)) throw new Error('--planeCScale must be >=0');
    if (!(isFinite(out.planeCJitter) && out.planeCJitter >= 0)) throw new Error('--planeCJitter must be >=0');
    if (!isFinite(out.E_SiH2) || !isFinite(out.E_SiH3) || !isFinite(out.E_bare) || !isFinite(out.E_bridge) || !isFinite(out.muH)) throw new Error('Energy coefficients must be finite');

    return out;
}

function loadMMParams(args) {
    const mm = new MMParams();
    const readText = (p) => fs.readFileSync(p, 'utf8');
    mm.parseElementTypes(readText(args.elementTypes));
    mm.parseAtomTypes(readText(args.atomTypes));
    mm.parseBondTypes(readText(args.bondTypes));
    mm.parseAngleTypes(readText(args.angleTypes));
    return mm;
}

function buildCrystalFromCIFText(cifText, opts) {
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

function buildPlanesFromTemplates(lvec, templates, planeSymC, planeMode, cmin, cmax) {
    const expanded = CrystalUtils.expandPlaneTemplates(templates, planeSymC);
    const b = CrystalUtils.reciprocalLattice(lvec);
    const plsFinal = [];
    for (const p of expanded) {
        const n = new Vec3().setLincomb3(p.h || 0, b[0], p.k || 0, b[1], p.l || 0, b[2]);
        plsFinal.push({ n, cmin, cmax });
    }
    return { plsFinal, planeMode };
}

function collapseProbabilisticBridges(mol, prob, requireH2) {
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

function randInt(a, b) {
    const i0 = a | 0;
    const i1 = b | 0;
    if (i1 < i0) throw new Error(`randInt: invalid range ${a}..${b}`);
    const u = Math.random();
    const x = i0 + Math.floor(u * (i1 - i0 + 1));
    return (x > i1) ? i1 : x;
}

function _atomNeighborsCounts(mol, ia) {
    const a = mol.atoms[ia];
    if (!a) return { nH: 0, nSi: 0, nHeavy: 0 };
    let nH = 0, nSi = 0, nHeavy = 0;
    for (let k = 0; k < a.bonds.length; k++) {
        const ib = a.bonds[k] | 0;
        const b = mol.bonds[ib];
        if (!b) continue;
        b.ensureIndices(mol);
        const ja = b.other(ia);
        if (ja < 0 || ja >= mol.atoms.length) continue;
        const nb = mol.atoms[ja];
        if (!nb) continue;
        if ((nb.Z | 0) === 1) nH++;
        else {
            nHeavy++;
            if ((nb.Z | 0) === 14) nSi++;
        }
    }
    return { nH, nSi, nHeavy };
}

function pruneUndercoordinatedSiIter(mol, minSiDegree, maxIter) {
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
            if (!a || (a.Z | 0) !== 14) continue;
            const c = _atomNeighborsCounts(mol, ia);
            if (c.nSi < (minSiDegree | 0)) {
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

function classifySurfaceCounts(mol) {
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
        if (z !== 14) continue;
        nSi++;
        const c = _atomNeighborsCounts(mol, ia);
        const isSurface = (c.nSi < 4); // diamond bulk has 4 Si neighbors
        if (!isSurface) continue; // bulk atoms ignored in surface classification
        if (c.nH === 0) {
            if (c.nSi === 2) nBridge++;
            else nBare++;
        } else if (c.nH === 1) nSiH++;
        else if (c.nH === 2) nSiH2++;
        else if (c.nH >= 3) nSiH3++;
    }
    return { nSi, nH, nSiH, nSiH2, nSiH3, nBare, nBridge };
}

function computeEnergy(cnt, args, extra = {}) {
    const nCollapsed = extra.nCollapsed || 0;
    const nInserted = extra.nInserted || 0;
    const E = args.E_SiH2 * cnt.nSiH2 + args.E_SiH3 * cnt.nSiH3 + args.E_bare * cnt.nBare + args.E_bridge * cnt.nBridge + args.muH * cnt.nH;
    return { E, nCollapsed, nInserted };
}

function mulberry32(seed) {
    let a = seed | 0;
    return function () {
        a |= 0;
        a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

async function main() {
    const args = parseArgs(process.argv);

    if (args.seed !== 0) {
        const rnd = mulberry32(args.seed);
        Math.random = rnd;
    }

    fs.mkdirSync(args.outDir, { recursive: true });

    const cifText = fs.readFileSync(args.cif, 'utf8');
    const cell = buildCrystalFromCIFText(cifText, args);

    const mm = (args.bonds || (args.caps && args.caps !== '0' && args.caps !== 'none')) ? loadMMParams(args) : null;

    const nTarget = (args.samples > 0) ? (args.samples | 0) : 0;
    const nLimit = (args.maxFiles > 0) ? (args.maxFiles | 0) : 0;
    const nOut = (nLimit > 0) ? (nTarget > 0 ? Math.min(nTarget, nLimit) : nLimit) : (nTarget > 0 ? nTarget : 0);

    if (args.samples <= 0 && args.maxFiles <= 0) {
        console.log('[gen_nanocrystals] NOTE: neither --samples nor --maxFiles specified; defaulting to 1 output');
    }
    const nGen = (nOut > 0) ? nOut : 1;

    if (args.statsCsv) {
        if (!fs.existsSync(args.statsCsv)) {
            fs.writeFileSync(args.statsCsv, 'i,nx,ny,nz,nAtoms,nBonds,nPruned,nCaps,nCollapsed,nInserted,nSi,nH,nSiH,nSiH2,nSiH3,nBare,nBridge,E\n');
        }
    }

    const stackedXYZ = args.stackedXyzOut ? [] : null;

    for (let iout = 1; iout <= nGen; iout++) {
        const nx = randInt(args.nx[0], args.nx[1]);
        const ny = randInt(args.ny[0], args.ny[1]);
        const nz = randInt(args.nz[0], args.nz[1]);

        const naEff = args.centered ? (2 * nx + 1) : nx;
        const nbEff = args.centered ? (2 * ny + 1) : ny;
        const ncEff = args.centered ? (2 * nz + 1) : nz;
        const La = cell.lvec[0].norm() * naEff;
        const Lb = cell.lvec[1].norm() * nbEff;
        const Lc = cell.lvec[2].norm() * ncEff;
        const Lavg = (La + Lb + Lc) / 3.0;

        const cBase = args.planeC0 + args.planeCScale * Lavg;
        const dj = args.planeCJitter * cBase;
        const jitter = dj * (2.0 * Math.random() - 1.0);
        const cmin = -cBase + jitter;
        const cmax = +cBase + jitter;

        const { plsFinal, planeMode } = buildPlanesFromTemplates(cell.lvec, args.planeTemplates, args.planeSymC, args.planeMode, cmin, cmax);

        const mol = CrystalUtils.genReplicatedCellCutPlanes({
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

        if (!mm) throw new Error('Internal: mmParams missing');

        mol.recalculateBonds(mm, { defaultRcut: args.defaultRcut, bondFactor: args.bondFactor });

        const pr = pruneUndercoordinatedSiIter(mol, args.minSiDegree, args.pruneMaxIter);
        if (pr.totalRemoved > 0) mol.recalculateBonds(mm, { defaultRcut: args.defaultRcut, bondFactor: args.bondFactor });

        let nCaps = 0;
        if (args.caps && args.caps !== '0' && args.caps !== 'none') {
            const r = mol.addCappingAtoms(mm, args.caps, { onlySelection: false, bBond: true });
            nCaps = r.nAdded | 0;
            if (nCaps > 0) mol.recalculateBonds(mm, { defaultRcut: args.defaultRcut, bondFactor: args.bondFactor });
        }

        const pr2 = pruneUndercoordinatedSiIter(mol, args.minSiDegree, args.pruneMaxIter);
        if (pr2.totalRemoved > 0) mol.recalculateBonds(mm, { defaultRcut: args.defaultRcut, bondFactor: args.bondFactor });

        // Surface SiH2 insert then collapse using generalized helpers
        let nCollapsed = 0;
        let nInserted = 0;
        if (args.collapseAll) {
            nCollapsed = collapseAllBridges(mol); // legacy carbon path
        } else {
            const surfaceFilter = (m, ia) => {
                const c = _atomNeighborsCounts(m, ia);
                return c.nSi < 4;
            };

            // Insert first to create bridgeable sites (surface-only)
            if (args.insertProb > 0) {
                const candidates = [];
                for (const b of mol.bonds) {
                    if (!b) continue; b.ensureIndices(mol);
                    const ia = b.a, ja = b.b;
                    const a = mol.atoms[ia], c = mol.atoms[ja];
                    if (!a || !c) continue;
                    if ((a.Z | 0) !== 14 || (c.Z | 0) !== 14) continue;
                    if (!surfaceFilter(mol, ia, a) || !surfaceFilter(mol, ja, c)) continue;
                    candidates.push([a.id, c.id]);
                }
                console.log(`[insert debug] surface Si-Si candidate bonds: ${candidates.length}`);

                for (const [aId, bId] of candidates) {
                    if (Math.random() >= args.insertProb) continue;
                    const cid = insertBridge(mol, aId, bId, { addHydrogens: true, hDist: 1.3, upOffsetFactor: 0.5 });
                    const ic = mol.getAtomIndex(cid);
                    if (ic >= 0) mol.atoms[ic].Z = 14; // ensure Si bridge
                    nInserted++;
                }
            }

            if (args.collapseProb > 0) {
                // select SiH2 bridge-like: need two heavy neighbors and >=2 H on surface
                mol.selection.clear();
                selectBridgeCandidates(mol, { z: 14, minHeavy: 2, minHyd: 2, requireH2: true, surfaceFilter });
                for (const id of Array.from(mol.selection)) {
                    if (Math.random() < args.collapseProb) {
                        collapseBridgeAt(mol, id);
                        nCollapsed++;
                    }
                }
            }
        }

        const cnt = classifySurfaceCounts(mol);
        const enr = computeEnergy(cnt, args, { nCollapsed, nInserted });

        const name = `${args.prefix}_i${String(iout).padStart(4, '0')}_nx${nx}_ny${ny}_nz${nz}_c${cBase.toFixed(3)}_tpl${args.planeTemplates.join('-')}_pr${(pr.totalRemoved + pr2.totalRemoved)}_cap${nCaps}_col${nCollapsed}_ins${nInserted}`;
        const mol2 = toMol2String(mol, { name });
        const outPath = path.join(args.outDir, name + '.mol2');
        fs.writeFileSync(outPath, mol2);

        if (stackedXYZ) {
            const xyz = toXYZString(mol, { lvec: mol.lvec });
            stackedXYZ.push(xyz.trimEnd());
        }

        const nPruned = (pr.totalRemoved + pr2.totalRemoved) | 0;
        console.log(`[gen_nanocrystals] wrote ${outPath} atoms=${mol.atoms.length} bonds=${mol.bonds.length} pruned=${nPruned} caps=${nCaps} collapsed=${nCollapsed} inserted=${nInserted} E=${enr.E.toFixed(6)} SiH=${cnt.nSiH} SiH2=${cnt.nSiH2} SiH3=${cnt.nSiH3} bare=${cnt.nBare}`);

        if (args.statsCsv) {
            const line = [
                iout, nx, ny, nz,
                mol.atoms.length, mol.bonds.length,
                nPruned, nCaps, nCollapsed, nInserted,
                cnt.nSi, cnt.nH, cnt.nSiH, cnt.nSiH2, cnt.nSiH3, cnt.nBare, cnt.nBridge,
                enr.E
            ].join(',') + '\n';
            fs.appendFileSync(args.statsCsv, line);
        }
    }

    if (stackedXYZ) {
        const xyzPath = args.stackedXyzOut;
        fs.writeFileSync(xyzPath, stackedXYZ.join('\n'));
        console.log(`[gen_nanocrystals] wrote stacked XYZ ${xyzPath} with ${stackedXYZ.length} frames`);
    }
}

main().catch((e) => {
    console.error(e);
    process.exit(1);
});
