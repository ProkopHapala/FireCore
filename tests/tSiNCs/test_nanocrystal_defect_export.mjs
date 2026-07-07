#!/usr/bin/env node
/// Defect + topology NPZ export regression harness (C diamond + Si planes).
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import { readNpzFile } from '../../web/common_js/npzIO.js';
import { exportNanocrystalBundle } from '../../web/molgui_webgpu/NanocrystalExport.js';
import {
    loadMMParams, buildCrystalFromCIFText, getHeavyZ, generateNanocrystal, mulberry32, defaultArgs,
} from '../../web/molgui_webgpu/Nanocrystals.js';

function atomDistIdx(mol, i, j) {
    const pi = mol.atoms[i].pos, pj = mol.atoms[j].pos;
    const dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(__dirname, '../..');
const RES = path.join(REPO, 'cpp/common_resources');
const OUT = path.join(REPO, 'tests/tSiNCs/out_defect_export');

function bondedPairSet(mol) {
    const s = new Set();
    for (const b of mol.bonds || []) {
        if (!b) continue;
        b.ensureIndices(mol);
        const ia = b.a | 0, ib = b.b | 0;
        if (ia < 0 || ib < 0) continue;
        const lo = Math.min(ia, ib), hi = Math.max(ia, ib);
        s.add(`${lo}-${hi}`);
    }
    return s;
}

function minHHClashDistance(mol) {
    const n = mol.atoms.length;
    const bonded = bondedPairSet(mol);
    let minD = Infinity;
    for (let i = 0; i < n; i++) {
        if ((mol.atoms[i].Z | 0) !== 1) continue;
        for (let j = i + 1; j < n; j++) {
            if ((mol.atoms[j].Z | 0) !== 1) continue;
            const lo = Math.min(i, j), hi = Math.max(i, j);
            if (bonded.has(`${lo}-${hi}`)) continue;
            const d = atomDistIdx(mol, i, j);
            if (d < minD) minD = d;
        }
    }
    return minD;
}

function countNeighbors(mol, ia) {
    let n = 0;
    const a = mol.atoms[ia];
    if (!a || !a.bonds) return 0;
    for (let k = 0; k < a.bonds.length; k++) {
        const ib = a.bonds[k] | 0;
        const b = mol.bonds[ib];
        if (b) n++;
    }
    return n;
}

function validateTopologyNpz(npzPath, mol) {
    const { arrays, entries } = readNpzFile(fs, npzPath);
    const n = arrays.natoms ? arrays.natoms[0] : arrays.Z.length;
    const nGroups = arrays.n_groups ? arrays.n_groups[0] : 0;
    if (!arrays.group_bbox_min) throw new Error(`${npzPath}: missing group_bbox_min`);
    if (!arrays.neigh_count) throw new Error(`${npzPath}: missing neigh_count`);
    const gRows = entries.get('group_bbox_min')?.shape?.[0] ?? (arrays.group_bbox_min.length / 3);
    if (gRows !== nGroups) throw new Error(`${npzPath}: group_bbox_min rows ${gRows} !== n_groups ${nGroups}`);
    for (let i = 0; i < n; i++) {
        const nc = arrays.neigh_count[i] | 0;
        const molNb = countNeighbors(mol, i);
        if (molNb > 0 && nc < 1) throw new Error(`${npzPath}: atom ${i} mol neighbors=${molNb} but neigh_count=${nc}`);
    }
    // non-isolated: every atom with Z>1 should have neigh_count >= 1
    for (let i = 0; i < n; i++) {
        if ((arrays.Z[i] | 0) <= 1) continue;
        if ((arrays.neigh_count[i] | 0) < 1) throw new Error(`${npzPath}: heavy atom ${i} neigh_count=${arrays.neigh_count[i]}`);
    }
    if (!arrays.schema_version || arrays.schema_version[0] !== 1) throw new Error(`${npzPath}: schema_version must be 1`);
    return { n, nGroups, n_angle: arrays.n_angle?.[0] ?? 0, n_bond: arrays.n_bond?.[0] ?? 0 };
}

function runCase(caseDef) {
    const { name, cif, heavyZ, seed, genOverrides, expect } = caseDef;
    const args = defaultArgs(RES, {
        cif: path.join(RES, cif),
        seed,
        caps: 'H',
        minHeavyDegree: 2,
        ...genOverrides,
    });
    Math.random = mulberry32(seed);
    const cell = buildCrystalFromCIFText(fs.readFileSync(args.cif, 'utf8'), args);
    const hz = getHeavyZ(cell);
    if ((hz | 0) !== (heavyZ | 0)) throw new Error(`${name}: CIF heavyZ=${hz} expected ${heavyZ}`);
    const mm = loadMMParams(args);
    const { mol, nCollapsed, nInserted, cnt } = generateNanocrystal(args, cell, mm, hz, 1);
    if (!(mol.atoms.length > 0)) throw new Error(`${name}: natoms must be > 0`);

    const minHH = minHHClashDistance(mol);
    if (minHH < 1.0) throw new Error(`${name}: H-H clash min dist=${minHH.toFixed(3)} Å < 1.0`);

    const outDir = path.join(OUT, name);
    fs.mkdirSync(outDir, { recursive: true });
    const genParams = { ...genOverrides, insertProb: args.insertProb, collapseProb: args.collapseProb, seed, heavyZ: hz };
    exportNanocrystalBundle(mol, mm, {
        outDir, id: name, genParams,
        defectsMeta: { nCollapsed, nInserted, insertProb: args.insertProb, collapseProb: args.collapseProb },
    });
    const topoInfo = validateTopologyNpz(path.join(outDir, '03_topology.npz'), mol);

    if (expect.nInsertedMin != null && nInserted < expect.nInsertedMin) throw new Error(`${name}: nInserted=${nInserted} < ${expect.nInsertedMin}`);
    if (expect.nCollapsedMin != null && nCollapsed < expect.nCollapsedMin) throw new Error(`${name}: nCollapsed=${nCollapsed} < ${expect.nCollapsedMin}`);
    if (expect.bridgeLtCaps != null && cnt.nBridge >= expect.bridgeLtCaps) throw new Error(`${name}: nBridge=${cnt.nBridge} expected < ${expect.bridgeLtCaps} (caps baseline)`);
    if (expect.nAngleMin != null && topoInfo.n_angle < expect.nAngleMin) throw new Error(`${name}: n_angle=${topoInfo.n_angle} < ${expect.nAngleMin}`);

    console.log(`[PASS] ${name} atoms=${mol.atoms.length} collapsed=${nCollapsed} inserted=${nInserted} bridge=${cnt.nBridge} minHH=${minHH.toFixed(3)} groups=${topoInfo.nGroups}`);
    return { name, nCollapsed, nInserted, cnt, topoInfo, minHH };
}

const CASES = [
    { name: 'C_caps_only', cif: 'crystals/diamond_primitive.cif', heavyZ: 6, seed: 1001,
      genOverrides: { applySymmetry: false, cutMode: 'sphere', sphereR: 6.0, sphereNrep: 5, insertProb: 0, collapseProb: 0 },
      expect: { nAngleMin: 1 } },
    { name: 'C_insert_only', cif: 'crystals/diamond_primitive.cif', heavyZ: 6, seed: 1,
      genOverrides: { applySymmetry: false, cutMode: 'sphere', sphereR: 6.0, sphereNrep: 5, insertProb: 0.2, collapseProb: 0 },
      expect: { nInsertedMin: 1, nAngleMin: 1 } },
    { name: 'C_collapse_only', cif: 'crystals/diamond_primitive.cif', heavyZ: 6, seed: 1,
      genOverrides: { applySymmetry: false, cutMode: 'sphere', sphereR: 6.0, sphereNrep: 5, insertProb: 0, collapseProb: 0.2 },
      expect: { nCollapsedMin: 1, nAngleMin: 1 } },
    { name: 'C_mixed_planes', cif: 'crystals/diamond_primitive.cif', heavyZ: 6, seed: 4,
      genOverrides: { applySymmetry: false, cutMode: 'planes', nx: [2, 2], ny: [2, 2], nz: [2, 2], planeTemplates: ['a111'], planeCScale: 0.40, planeCJitter: 0, insertProb: 0.15, collapseProb: 0.10 },
      expect: { nAngleMin: 1 } },
    { name: 'Si_caps_only', cif: 'crystals/Si-sym.cif', heavyZ: 14, seed: 1,
      genOverrides: { cutMode: 'sphere', sphereR: 6.0, sphereNrep: 5, insertProb: 0, collapseProb: 0 },
      expect: { nAngleMin: 1 } },
    { name: 'Si_insert_only', cif: 'crystals/Si-sym.cif', heavyZ: 14, seed: 1,
      genOverrides: { cutMode: 'sphere', sphereR: 6.0, sphereNrep: 5, insertProb: 0.2, collapseProb: 0 },
      expect: { nInsertedMin: 1, nAngleMin: 1 } },
    { name: 'Si_collapse_only', cif: 'crystals/Si-sym.cif', heavyZ: 14, seed: 1,
      genOverrides: { cutMode: 'sphere', sphereR: 6.0, sphereNrep: 5, insertProb: 0, collapseProb: 0.2 },
      expect: { nCollapsedMin: 1, nAngleMin: 1 } },
    { name: 'Si_mixed_planes', cif: 'crystals/Si-sym.cif', heavyZ: 14, seed: 3,
      genOverrides: { cutMode: 'planes', nx: [2, 2], ny: [2, 2], nz: [2, 2], planeTemplates: ['a111'], planeCScale: 0.40, planeCJitter: 0, insertProb: 0.15, collapseProb: 0.10 },
      expect: { nAngleMin: 1 } },
];

function main() {
    fs.mkdirSync(OUT, { recursive: true });
    const results = new Map();
    for (const c of CASES) {
        if (c.name.endsWith('_collapse_only')) {
            const prefix = c.heavyZ === 6 ? 'C_caps_only' : 'Si_caps_only';
            const base = results.get(prefix);
            if (base) c.expect = { ...c.expect, bridgeLtCaps: base.cnt.nBridge + 1 };
        }
        results.set(c.name, runCase(c));
    }
    console.log(`\n[test_nanocrystal_defect_export] all ${CASES.length} cases passed`);
}

main();
