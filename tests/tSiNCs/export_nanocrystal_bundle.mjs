#!/usr/bin/env node
/// Generate nanocrystal + export crystal/topology NPZ bundle for viewer fixtures and pipeline QA.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import { exportNanocrystalBundle } from '../../web/molgui_webgpu/NanocrystalExport.js';
import {
    loadMMParams, buildCrystalFromCIFText, getHeavyZ, generateNanocrystal, mulberry32, defaultArgs,
} from '../../web/molgui_webgpu/Nanocrystals.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const TEST_DIR = __dirname;
const REPO = path.resolve(__dirname, '../..');
const RES = path.join(REPO, 'cpp/common_resources');

function parseProb(s, name) {
    const v = +String(s).trim();
    if (!isFinite(v) || v < 0 || v > 1) throw new Error(`${name}: expected probability in [0,1], got '${s}'`);
    return v;
}

function parseArgs(argv) {
    const out = defaultArgs(RES, {
        cif: path.join(RES, 'crystals/diamond_primitive.cif'),
        outDir: path.join(TEST_DIR, 'OUT_nanocrystal_bundle'),
        id: 'nanocrystal',
        seed: 42,
        insertProb: 0.0,
        collapseProb: 0.0,
        fuseProb: 0.0,
        fuseHClashMax: 2.0,
        fuseSiMax: 4.5,
        exportTopology: true,
        exportSurfaceGroups: true,
        addAngle: true,
        addDihedral: true,
        K12: 0, K13: 0, K14: 0,
        groupCap: 32,
        maxNeighbors: 48,
        crystalOnly: false,
    });
    for (let i = 2; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => { if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`); i++; return String(argv[i]); };
        if (a === '--cif') out.cif = path.resolve(nxt());
        else if (a === '--applySymmetry') out.applySymmetry = (nxt() !== '0');
        else if (a === '--out' || a === '--outDir') out.outDir = path.resolve(nxt());
        else if (a === '--id') out.id = nxt();
        else if (a === '--seed') out.seed = parseInt(nxt(), 10) | 0;
        else if (a === '--heavyZ') { /* CIF defines heavy element; kept for CLI parity */ nxt(); }
        else if (a === '--sphere' || a === '--sphereR') { out.cutMode = 'sphere'; out.sphereR = +nxt(); }
        else if (a === '--sphereNrep') out.sphereNrep = parseInt(nxt(), 10) | 0;
        else if (a === '--planes') { out.cutMode = 'planes'; out.planeTemplates = nxt().split(',').map(s => s.trim()).filter(Boolean); }
        else if (a === '--nx') { const v = parseInt(nxt(), 10); out.nx = [v, v]; out.ny = [v, v]; out.nz = [v, v]; }
        else if (a === '--planeCScale') out.planeCScale = +nxt();
        else if (a === '--planeCJitter') out.planeCJitter = +nxt();
        else if (a === '--insertProb') out.insertProb = parseProb(nxt(), '--insertProb');
        else if (a === '--collapseProb') out.collapseProb = parseProb(nxt(), '--collapseProb');
        else if (a === '--fuseProb') out.fuseProb = parseProb(nxt(), '--fuseProb');
        else if (a === '--fuseHClashMax') out.fuseHClashMax = +nxt();
        else if (a === '--fuseSiMax') out.fuseSiMax = +nxt();
        else if (a === '--collapseAll') out.collapseAll = (nxt() !== '0');
        else if (a === '--caps') out.caps = nxt();
        else if (a === '--resolveClashes') out.resolveClashes = (nxt() !== '0');
        else if (a === '--exportTopology') out.exportTopology = (nxt() !== '0');
        else if (a === '--exportSurfaceGroups') out.exportSurfaceGroups = (nxt() !== '0');
        else if (a === '--addAngle') out.addAngle = (nxt() !== '0');
        else if (a === '--addDihedral') out.addDihedral = (nxt() !== '0');
        else if (a === '--K12') out.K12 = +nxt();
        else if (a === '--K13') out.K13 = +nxt();
        else if (a === '--K14') out.K14 = +nxt();
        else if (a === '--groupCap') out.groupCap = parseInt(nxt(), 10) | 0;
        else if (a === '--maxNeighbors') out.maxNeighbors = parseInt(nxt(), 10) | 0;
        else if (a === '--crystalOnly') out.crystalOnly = (nxt() !== '0');
        else if (a === '--help' || a === '-h') { console.log(usage()); process.exit(0); }
        else throw new Error(`Unknown arg ${a}`);
    }
    if (!(out.sphereR > 0) && out.cutMode === 'sphere') throw new Error('--sphereR must be >0 for sphere cut');
    if (path.basename(out.cif).includes('diamond_primitive')) out.applySymmetry = false;
    return out;
}

function usage() {
    return `Usage: node tests/tSiNCs/export_nanocrystal_bundle.mjs [options]

Generate nanocrystal + write 01_crystal.npz, 03_topology.npz, meta.json under --out.

Options:
  --cif PATH              CIF file (default: diamond_primitive.cif)
  --out PATH              output directory
  --id STR                fixture id / meta label
  --seed N                RNG seed (default 42)
  --sphere R              sphere cut, radius R Å (sets cutMode=sphere)
  --sphereNrep N          sphere replication (default 5)
  --planes a111,...       plane cut templates (sets cutMode=planes)
  --nx N                  cubic replication N×N×N for plane cuts
  --insertProb P          bridge insert probability [0,1]
  --collapseProb P        bridge collapse probability [0,1]
  --caps H|none           capping element
  --exportTopology 0|1    write 03_topology.npz (default 1)
  --crystalOnly 1         alias: --exportTopology 0

Example:
  node scripts/export_nanocrystal_bundle.mjs \\
    --cif cpp/common_resources/crystals/diamond_primitive.cif \\
    --sphere 6 --insertProb 0.15 --collapseProb 0.10 \\
    --out tests/tSiNCs/fixtures/npz_viewer/diamond_sphere_R6_defects`;
}

function main() {
    const args = parseArgs(process.argv);
    if (args.crystalOnly) args.exportTopology = false;
    if (args.seed !== 0) Math.random = mulberry32(args.seed);

    const cifText = fs.readFileSync(args.cif, 'utf8');
    const cell = buildCrystalFromCIFText(cifText, args);
    const heavyZ = getHeavyZ(cell);
    const mm = loadMMParams(args);

    const result = generateNanocrystal(args, cell, mm, heavyZ, 1);
    const { mol, nCollapsed, nInserted, nFused, cnt, name } = result;

    const genParams = {
        cif: args.cif,
        cutMode: args.cutMode,
        sphereR: args.sphereR,
        sphereNrep: args.sphereNrep,
        planeTemplates: args.planeTemplates,
        nx: args.nx, ny: args.ny, nz: args.nz,
        insertProb: args.insertProb,
        collapseProb: args.collapseProb,
        fuseProb: args.fuseProb,
        fuseHClashMax: args.fuseHClashMax,
        caps: args.caps,
        seed: args.seed,
        heavyZ,
        name,
    };

    const bundle = exportNanocrystalBundle(mol, mm, {
        outDir: args.outDir,
        id: args.id,
        genParams,
        exportTopology: args.exportTopology,
        exportSurfaceGroups: args.exportSurfaceGroups,
        addAngle: args.addAngle,
        addDihedral: args.addDihedral,
        K12: args.K12, K13: args.K13, K14: args.K14,
        maxNeighbors: args.maxNeighbors,
        groupCap: args.groupCap,
        defectsMeta: { nCollapsed, nInserted, nFused, insertProb: args.insertProb, collapseProb: args.collapseProb, fuseProb: args.fuseProb },
    });

    console.log(`[export_nanocrystal_bundle] id=${args.id} atoms=${mol.atoms.length} bonds=${mol.bonds.length}`);
    console.log(`  collapsed=${nCollapsed} inserted=${nInserted} fused=${nFused} bridge=${cnt.nBridge} SiH2=${cnt.nSiH2}`);
    console.log(`  crystal: ${bundle.crystalNpz}`);
    if (bundle.topologyNpz) console.log(`  topology: ${bundle.topologyNpz}`);
    console.log(`  meta: ${bundle.metaPath}`);
}

main();
