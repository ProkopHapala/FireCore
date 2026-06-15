#!/usr/bin/env node
/** Geometry audit for Si/C nanocrystal passivation (headless Node). Writes geometry_report.json. */
import fs from 'node:fs';
import path from 'node:path';

import { Vec3 } from '../web/common_js/Vec3.js';
import { MMParams } from '../web/molgui_webgpu/MMParams.js';
import { EditableMolecule } from '../web/molgui_webgpu/EditableMolecule.js';
import { installMoleculeIOMethods } from '../web/molgui_webgpu/MoleculeIO.js';

installMoleculeIOMethods(EditableMolecule);

const REPO = path.resolve(path.dirname(new URL(import.meta.url).pathname), '..');
const FIX = path.join(REPO, 'tests/tSiNCs/fixtures/vibration_parallel/structures');
const RES = path.join(REPO, 'cpp/common_resources');

const PRESETS = {
    G0: { kind: 'build', build: 'sih4' },
    G1: { kind: 'gen', genArgs: ['--nx-range', '1,1', '--ny-range', '1,1', '--nz-range', '1,1', '--centered', '1', '--planeTemplates', 'a111', '--planeSymC', '6', '--planeCScale', '0.50', '--planeCJitter', '0', '--minSiDegree', '2'] },
    G2: { kind: 'gen', genArgs: ['--nx-range', '2,2', '--ny-range', '2,2', '--nz-range', '2,2', '--centered', '1', '--planeTemplates', 'a111', '--planeSymC', '6', '--planeCScale', '0.40', '--planeCJitter', '0', '--minSiDegree', '2'] },
    G3: { kind: 'gen', genArgs: ['--nx-range', '2,2', '--ny-range', '2,2', '--nz-range', '2,2', '--centered', '1', '--planeTemplates', 'a111', '--planeSymC', '6', '--planeCScale', '0.40', '--planeCJitter', '0', '--insertProb', '0.2', '--collapseProb', '0.3'] },
    G4: { kind: 'file', mol2: path.join(FIX, 'adamantane.mol2') },
    G5: { kind: 'file', xyz: path.join(FIX, 'diamond_primitive.xyz') },
};

function loadMM() {
    const mm = new MMParams();
    mm.parseElementTypes(fs.readFileSync(path.join(RES, 'ElementTypes.dat'), 'utf8'));
    mm.parseAtomTypes(fs.readFileSync(path.join(RES, 'AtomTypes.dat'), 'utf8'));
    mm.parseBondTypes(fs.readFileSync(path.join(RES, 'BondTypes.dat'), 'utf8'));
    return mm;
}

function dist(a, b) { return Math.sqrt(a.pos.dist2(b.pos)); }

function angleDeg(oa, ob, oc) {
    const u = new Vec3().setV(ob.pos).sub(oa.pos);
    const v = new Vec3().setV(oc.pos).sub(oa.pos);
    const nu = u.norm(), nv = v.norm();
    if (!(nu > 1e-12 && nv > 1e-12)) return NaN;
    const c = u.dot(v) / (nu * nv);
    return Math.acos(Math.max(-1, Math.min(1, c))) * (180 / Math.PI);
}

function loadMolFromFile(p) {
    const ext = path.extname(p).toLowerCase();
    const text = fs.readFileSync(p, 'utf8');
    const parsed = (ext === '.xyz') ? EditableMolecule.parseXYZ(text) : EditableMolecule.parseMol2(text);
    const mol = new EditableMolecule();
    mol.appendParsedSystem(parsed);
    return mol;
}

function buildSiH4() {
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

function bondedPairSet(mol) {
    const s = new Set();
    for (const b of mol.bonds) {
        if (!b) continue;
        b.ensureIndices(mol);
        const ia = b.a | 0, ib = b.b | 0;
        if (ia < 0 || ib < 0) continue;
        const lo = Math.min(ia, ib), hi = Math.max(ia, ib);
        s.add(`${lo}-${hi}`);
    }
    return s;
}

function neighborIndices(mol, ia) {
    const out = [];
    const a = mol.atoms[ia];
    for (const ib of a.bonds) {
        const b = mol.bonds[ib | 0];
        if (!b) continue;
        b.ensureIndices(mol);
        const jb = b.other(ia);
        if (jb >= 0) out.push(jb);
    }
    return out;
}

function isBulkSi(mol, ia) {
    const a = mol.atoms[ia];
    if ((a.Z | 0) !== 14) return false;
    let nSi = 0;
    for (const jb of neighborIndices(mol, ia)) if ((mol.atoms[jb].Z | 0) === 14) nSi++;
    return nSi === 4;
}

function auditGeometry(mol, mm, opts = {}) {
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
                const d = dist(a, mol.atoms[jb]);
                report.siH.n++;
                if (report.siH.min === null || d < report.siH.min) report.siH.min = d;
                if (report.siH.max === null || d > report.siH.max) report.siH.max = d;
                if (d < siHLo || d > siHHi) report.siH.bad.push({ siId: a.id, hId: mol.atoms[jb].id, d });
            }
            if (nbs.length >= 3) {
                const bulk = isBulkSi(mol, ia);
                for (let i = 0; i < nbs.length; i++) for (let j = i + 1; j < nbs.length; j++) {
                    const ang = angleDeg(a, mol.atoms[nbs[i]], mol.atoms[nbs[j]]);
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
            const d = dist(mol.atoms[i], mol.atoms[j]);
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

function oneLineSummary(name, r) {
    const s = r.summary;
    return `[${name}] atoms=${r.natoms} Si-H=[${r.siH.min?.toFixed(3) ?? '?'},${r.siH.max?.toFixed(3) ?? '?'}] clashes=${r.clashes.length} bulkAngBad=${r.angles.bad.length} surfAngExc=${r.angles.surfaceBad.length} valenceBad=${r.valence.bad.length} PASS=${s.pass}`;
}

async function runGenPreset(genArgs, outDir) {
    const { spawnSync } = await import('node:child_process');
    fs.mkdirSync(outDir, { recursive: true });
    const cmd = ['node', 'scripts/gen_nanocrystals.mjs', '--samples', '1', '--seed', '42', '--maxFiles', '1', '--caps', 'H', '--insertProb', '0', '--collapseProb', '0', '--outDir', outDir, '--prefix', 'audit', ...genArgs];
    const res = spawnSync(cmd[0], cmd.slice(1), { cwd: REPO, encoding: 'utf8' });
    if (res.status !== 0) throw new Error(`gen_nanocrystals failed:\n${res.stdout}\n${res.stderr}`);
    const mol2s = fs.readdirSync(outDir).filter(f => f.endsWith('.mol2')).sort();
    if (!mol2s.length) throw new Error(`no mol2 in ${outDir}`);
    return path.join(outDir, mol2s[mol2s.length - 1]);
}

function parseArgs(argv) {
    const out = { preset: null, input: null, outJson: null, outDir: path.join(REPO, 'tests/tSiNCs/geometry'), bondFactor: 1.30, defaultRcut: 1.60, failOnError: true };
    for (let i = 2; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => { if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`); return String(argv[++i]); };
        if (a === '--preset') out.preset = nxt().toUpperCase();
        else if (a === '--input') out.input = path.resolve(nxt());
        else if (a === '--outJson') out.outJson = path.resolve(nxt());
        else if (a === '--outDir') out.outDir = path.resolve(nxt());
        else if (a === '--bondFactor') out.bondFactor = +nxt();
        else if (a === '--defaultRcut') out.defaultRcut = +nxt();
        else if (a === '--no-fail') out.failOnError = false;
        else if (a === '--all-presets') out.allPresets = true;
        else throw new Error(`Unknown arg ${a}`);
    }
    return out;
}

async function runOne(name, mol, mm, bondOpts, outDir) {
    if (mol.bonds.length === 0) mol.recalculateBonds(mm, bondOpts);
    const report = auditGeometry(mol, mm);
    report.preset = name;
    const line = oneLineSummary(name, report);
    console.log(line);
    fs.mkdirSync(outDir, { recursive: true });
    const jsonPath = path.join(outDir, `${name.toLowerCase()}_geometry_report.json`);
    fs.writeFileSync(jsonPath, JSON.stringify(report, null, 2));
    const xyzPath = path.join(outDir, `${name.toLowerCase()}.xyz`);
    const { toXYZString } = await import('../web/molgui_webgpu/MoleculeIO.js');
    fs.writeFileSync(xyzPath, toXYZString(mol, {}));
    return { report, jsonPath, line };
}

async function main() {
    const args = parseArgs(process.argv);
    const mm = loadMM();
    const bondOpts = { bondFactor: args.bondFactor, defaultRcut: args.defaultRcut };
    const results = [];

    if (args.allPresets) {
        for (const key of Object.keys(PRESETS)) {
            const p = PRESETS[key];
            let mol;
            if (p.kind === 'build') mol = buildSiH4();
            else if (p.kind === 'file') {
                const fp = p.mol2 || p.xyz;
                if (!fs.existsSync(fp)) { console.log(`[${key}] SKIP missing ${fp}`); continue; }
                mol = loadMolFromFile(fp);
            } else if (p.kind === 'gen') {
                const genOut = path.join(args.outDir, `gen_${key.toLowerCase()}`);
                const mol2 = await runGenPreset(p.genArgs, genOut);
                mol = loadMolFromFile(mol2);
            }
            results.push(await runOne(key, mol, mm, bondOpts, args.outDir));
        }
    } else if (args.preset) {
        const p = PRESETS[args.preset];
        if (!p) throw new Error(`Unknown preset ${args.preset}; known: ${Object.keys(PRESETS).join(', ')}`);
        let mol;
        if (p.kind === 'build') mol = buildSiH4();
        else if (p.kind === 'file') mol = loadMolFromFile(p.mol2 || p.xyz);
        else if (p.kind === 'gen') {
            const genOut = path.join(args.outDir, `gen_${args.preset.toLowerCase()}`);
            mol = loadMolFromFile(await runGenPreset(p.genArgs, genOut));
        }
        results.push(await runOne(args.preset, mol, mm, bondOpts, args.outDir));
    } else if (args.input) {
        results.push(await runOne(path.basename(args.input), loadMolFromFile(args.input), mm, bondOpts, args.outDir));
    } else {
        console.error('usage: node scripts/test_nanocrystal_geometry.mjs --preset G2 | --input file.mol2 | --all-presets [--outDir DIR]');
        process.exit(1);
    }

    if (args.outJson && results.length === 1) fs.copyFileSync(results[0].jsonPath, args.outJson);

    const failed = results.filter(r => !r.report.summary.pass);
    const passList = ['G0', 'G1', 'G2'];
    const requiredFailed = failed.filter(r => passList.includes(r.report.preset));
    if (args.failOnError && requiredFailed.length > 0) {
        console.error(`[test_nanocrystal_geometry] FAILED required ${requiredFailed.map(r => r.report.preset).join(', ')}`);
        process.exit(1);
    }
    if (args.failOnError && args.allPresets && failed.length > requiredFailed.length) {
        console.log(`[test_nanocrystal_geometry] optional presets failed: ${failed.filter(r => !passList.includes(r.report.preset)).map(r => r.report.preset).join(', ')}`);
    }
    console.log(`[test_nanocrystal_geometry] done ${results.length} case(s)`);
}

main().catch(e => { console.error(e); process.exit(1); });
