/// @deprecated Use `node tests/tSiNCs/nanocrystals.mjs generate` instead. This script is kept for backward compatibility.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import { toMol2String, toXYZString } from '../../web/molgui_webgpu/MoleculeIO.js';
import {
    loadMMParams, buildCrystalFromCIFText, getHeavyZ, generateNanocrystal, mulberry32,
} from '../../web/molgui_webgpu/Nanocrystals.js';

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
    const __dirname = path.dirname(fileURLToPath(import.meta.url));
const TEST_DIR = __dirname;
const repoRoot = path.resolve(__dirname, '../..');
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

        collapseProb: 0.0,
        insertProb: 0.0,
        fuseProb: 0.0,
        fuseHClashMax: 2.0,
        outwardBias: 0.35,
        resolveClashes: 1,
        capHHBonds: 0,
        capHHBondDist: 1.8,

        cutMode: 'planes',
        sphereR: 6.0,
        sphereNrep: 5,
        rcutHeavy: 0.0,
        minHeavyDegree: 2,

        E_SiH2: 1.0,
        E_SiH3: 4.0,
        E_bare: 10.0,
        E_bridge: 2.0,
        muH: 0.0,

        outDir: path.resolve(TEST_DIR, 'OUT_nanocrystals'),
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
        else if (a === '--minSiDegree') out.minHeavyDegree = parseInt(nxt(), 10) | 0;
        else if (a === '--pruneMaxIter') out.pruneMaxIter = parseInt(nxt(), 10) | 0;

        else if (a === '--collapse') out.collapse = parseInt(nxt(), 10) | 0;
        else if (a === '--collapseAll') out.collapseAll = (nxt() !== '0');
        else if (a === '--requireH2') out.requireH2 = (nxt() !== '0');

        else if (a === '--collapseProb') out.collapseProb = parseProb(nxt(), '--collapseProb');
        else if (a === '--insertProb') out.insertProb = parseProb(nxt(), '--insertProb');
        else if (a === '--fuseProb') out.fuseProb = parseProb(nxt(), '--fuseProb');
        else if (a === '--fuseHClashMax') out.fuseHClashMax = +nxt();
        else if (a === '--outwardBias') out.outwardBias = +nxt();
        else if (a === '--resolveClashes') out.resolveClashes = (nxt() !== '0');
        else if (a === '--capHHBonds') out.capHHBonds = (nxt() !== '0');
        else if (a === '--capHHBondDist') out.capHHBondDist = +nxt();

        else if (a === '--cutMode') out.cutMode = nxt();
        else if (a === '--sphereR' || a === '--sphere-r') out.sphereR = +nxt();
        else if (a === '--sphereNrep' || a === '--sphere-nrep') out.sphereNrep = parseInt(nxt(), 10) | 0;
        else if (a === '--rcutHeavy' || a === '--rcut-heavy') out.rcutHeavy = +nxt();
        else if (a === '--minHeavyDegree') out.minHeavyDegree = parseInt(nxt(), 10) | 0;

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
    if (out.cutMode === 'planes' && (!Array.isArray(out.planeTemplates) || out.planeTemplates.length === 0)) throw new Error('--planeTemplates must be non-empty for cutMode=planes');
    if (!(out.dedupTol > 0)) throw new Error('--dedupTol must be >0');
    if (out.collapse < 0) throw new Error('--collapse must be >=0');
    if (out.samples < 0) throw new Error('--samples must be >=0');
    if (out.minSiDegree < 0) throw new Error('--minSiDegree must be >=0');
    if (out.pruneMaxIter < 0) throw new Error('--pruneMaxIter must be >=0');
    if (!(isFinite(out.planeCScale) && out.planeCScale >= 0)) throw new Error('--planeCScale must be >=0');
    if (!(isFinite(out.planeCJitter) && out.planeCJitter >= 0)) throw new Error('--planeCJitter must be >=0');
    if (!(isFinite(out.capHHBondDist) && out.capHHBondDist > 0)) throw new Error('--capHHBondDist must be >0');
    if (out.cutMode !== 'planes' && out.cutMode !== 'sphere') throw new Error("--cutMode must be 'planes' or 'sphere'");
    out.minSiDegree = out.minHeavyDegree;
    if (!isFinite(out.outwardBias) || out.outwardBias < 0 || out.outwardBias > 1) throw new Error('--outwardBias must be in [0,1]');
    if (out.rcutHeavy < 0) throw new Error('--rcutHeavy must be >=0 (0=auto)');
    if (!(out.sphereR > 0)) throw new Error('--sphereR must be >0');
    if (out.sphereNrep < 0) throw new Error('--sphereNrep must be >=0');
    if (!isFinite(out.E_SiH2) || !isFinite(out.E_SiH3) || !isFinite(out.E_bare) || !isFinite(out.E_bridge) || !isFinite(out.muH)) throw new Error('Energy coefficients must be finite');

    return out;
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
    const heavyZ = getHeavyZ(cell);

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
            fs.writeFileSync(args.statsCsv, 'i,nx,ny,nz,nAtoms,nBonds,nPruned,nCaps,nHHBonds,nCollapsed,nInserted,nSi,nH,nSiH,nSiH2,nSiH3,nBare,nBridge,E\n');
        }
    }

    const stackedXYZ = args.stackedXyzOut ? [] : null;

    for (let iout = 1; iout <= nGen; iout++) {
        const result = generateNanocrystal(args, cell, mm, heavyZ, iout);
        const { mol, nCollapsed, nInserted, nFused, nCaps, nPruned, nHHBonds, cnt, enr, name, nx, ny, nz } = result;

        const mol2 = toMol2String(mol, { name });
        const outPath = path.join(args.outDir, name + '.mol2');
        fs.writeFileSync(outPath, mol2);
        fs.writeFileSync(path.join(args.outDir, name + '.xyz'), toXYZString(mol, { lvec: mol.lvec }));

        if (stackedXYZ) {
            const xyz = toXYZString(mol, { lvec: mol.lvec });
            stackedXYZ.push(xyz.trimEnd());
        }

        console.log(`[gen_nanocrystals] wrote ${outPath} atoms=${mol.atoms.length} bonds=${mol.bonds.length} pruned=${nPruned} caps=${nCaps} hhBonds=${nHHBonds} collapsed=${nCollapsed} inserted=${nInserted} fused=${nFused} E=${enr.E.toFixed(6)} SiH=${cnt.nSiH} SiH2=${cnt.nSiH2} SiH3=${cnt.nSiH3} bare=${cnt.nBare}`);

        if (args.statsCsv) {
            const line = [
                iout, nx, ny, nz,
                mol.atoms.length, mol.bonds.length,
                nPruned, nCaps, nHHBonds, nCollapsed, nInserted,
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
