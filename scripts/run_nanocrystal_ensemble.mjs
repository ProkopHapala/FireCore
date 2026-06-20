#!/usr/bin/env node
/// Nanocrystal ensemble pipeline — single entrypoint (JS generate/topology + Python MMFF/eigh).
import fs from 'node:fs';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

import { molToCrystalArrays, writeCrystalNpz, crystalToXYZ, readXyzPositions, crystalToJson, writeCrystalJson, readNpzFile } from '../web/common_js/npzIO.js';
import { exportCrystalCompareSvgViews, atlasIndexHtml } from '../web/common_js/nanocrystalSvg.js';
import { buildTopologyNpz } from '../web/common_js/exportFF.js';
import { loadMMParamsFromDir, loadMolFromMol2, bondsForVisualization, getCrystalBondsFromFiles } from '../web/common_js/MolIO.js';
const __dirname = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(__dirname, '..');

function usage() {
    return `Usage: node scripts/run_nanocrystal_ensemble.mjs [options]

Options:
  --config PATH          ensemble JSON (default: scripts/ensemble.example.json)
  --data-dir PATH        per-crystal NPZ cache (default: OUT_nc_ensemble/data)
  --output-dir PATH      accumulated outputs (default: OUT_nc_ensemble/out)
  --work-dir PATH        scratch (default: OUT_nc_ensemble/work)
  --n-crystals N         override config n_crystals
  --python PATH          Python executable (default: python3 or $PYTHON)
  --force                re-run all stages
  --no-debug             skip per-crystal SVG debug plots
  --atlas PATH           shape showcase atlas (JSON entries, no full ensemble)
  --help                 show this message

Viewer: pipeline writes self-contained out/viewer.html (works from file://) + crystals/*/crystal_*.json`;
}

function parseArgs(argv) {
    const out = {
        config: path.join(REPO, 'scripts/ensemble.example.json'),
        dataDir: path.join(REPO, 'OUT_nc_ensemble/data'),
        outputDir: path.join(REPO, 'OUT_nc_ensemble/out'),
        workDir: path.join(REPO, 'OUT_nc_ensemble/work'),
        nCrystals: null,
        python: process.env.PYTHON || 'python3',
        force: false,
        debug: null,
        atlas: null,
    };
    for (let i = 2; i < argv.length; i++) {
        const a = String(argv[i]);
        if (a === '--help' || a === '-h') { console.log(usage()); process.exit(0); }
        const nxt = () => { if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`); i++; return String(argv[i]); };
        if (a === '--config') out.config = nxt();
        else if (a === '--data-dir') out.dataDir = nxt();
        else if (a === '--output-dir') out.outputDir = nxt();
        else if (a === '--work-dir') out.workDir = nxt();
        else if (a === '--n-crystals') out.nCrystals = parseInt(nxt(), 10);
        else if (a === '--python') out.python = nxt();
        else if (a === '--force') out.force = true;
        else if (a === '--no-debug') out.debug = false;
        else if (a === '--atlas') out.atlas = nxt();
        else throw new Error(`Unknown arg: ${a}\n\n${usage()}`);
    }
    return out;
}

function mulberry32(seed) {
    let a = seed | 0;
    return () => {
        a |= 0;
        a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

function sampleValue(spec, rnd) {
    if (spec == null) return null;
    if (typeof spec === 'number' || typeof spec === 'string' || typeof spec === 'boolean') return spec;
    if (spec.uniform) return spec.uniform[0] + rnd() * (spec.uniform[1] - spec.uniform[0]);
    if (spec.uniform_int) return (spec.uniform_int[0] + Math.floor(rnd() * (spec.uniform_int[1] - spec.uniform_int[0] + 1))) | 0;
    if (spec.choice) return spec.choice[Math.floor(rnd() * spec.choice.length)];
    throw new Error(`sampleValue: unknown spec ${JSON.stringify(spec)}`);
}

function sampleCutSpec(cutMixture, rnd) {
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

function crystalId(seed, index) {
    return `${String(seed >>> 0).padStart(8, '0')}_${String(index).padStart(6, '0')}`;
}

function appendManifest(manifestPath, rec) {
    fs.appendFileSync(manifestPath, JSON.stringify(rec) + '\n');
}

function pyEnv() {
    return { ...process.env, PYTHONPATH: [REPO, process.env.PYTHONPATH].filter(Boolean).join(path.delimiter) };
}

function runCmd(cmd, args, { cwd = REPO, label = '', statusPath = null } = {}) {
    const t0 = performance.now();
    const r = spawnSync(cmd, args, { cwd, encoding: 'utf8', stdio: ['ignore', 'ignore', 'pipe'], env: pyEnv() });
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

function runPy(subcmd, pyArgs, opts) {
    return runCmd(opts.python || 'python3', ['-m', 'pyBall.nanocrystal_pipeline', subcmd, ...pyArgs], opts);
}

function stageExists(p) { return fs.existsSync(p); }

function buildGenCli(genParams, outDir, seed) {
    const g = genParams;
    const args = [
        'scripts/gen_nanocrystals.mjs',
        '--cif', path.join(REPO, g.cif),
        '--applySymmetry', String(g.applySymmetry ?? 0),
        '--samples', '1',
        '--seed', String(seed),
        '--maxFiles', '1',
        '--outDir', outDir,
        '--prefix', 'nc',
        '--caps', g.passivation?.caps ?? 'H',
        '--insertProb', String(g.defects?.insertProb ?? 0),
        '--collapseProb', String(g.defects?.collapseProb ?? 0),
        '--cutMode', g.cut.cutMode,
    ];
    const p = g.passivation || {};
    if (p.capHHBonds) args.push('--capHHBonds', '1', '--capHHBondDist', String(p.capHHBondDist ?? 1.8));
    args.push('--resolveClashes', (p.resolveClashes ? '1' : '0'));
    if (p.outwardBias != null) args.push('--outwardBias', String(p.outwardBias));
    const c = g.cut;
    if (c.cutMode === 'sphere') {
        args.push('--sphereR', String(c.sphereR), '--sphereNrep', String(c.sphereNrep | 0));
    } else if (c.cutMode === 'planes') {
        const nx = g.replication?.nx ?? 2;
        const ny = g.replication?.ny ?? 2;
        const nz = g.replication?.nz ?? 2;
        args.push('--nx-range', `${nx},${nx}`, '--ny-range', `${ny},${ny}`, '--nz-range', `${nz},${nz}`);
        const tpl = c.planeTemplates || ['a111'];
        args.push('--planeTemplates', Array.isArray(tpl) ? tpl.join(',') : String(tpl));
        args.push('--planeSymC', String(c.planeSymC ?? 6));
        args.push('--planeMode', c.planeMode || 'ang');
        args.push('--planeC0', String(c.planeC0 ?? 0));
        args.push('--planeCScale', String(c.planeCScale ?? 0.45));
        args.push('--planeCJitter', String(c.planeCJitter ?? 0));
    } else {
        throw new Error(`cutMode ${c.cutMode} not implemented in pass 1`);
    }
    return args;
}

function writeTimingReport(outputDir, stageTimings) {
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

function writeCrystalViewerJson(plotDir, { id, label, posInit, Zinit, bondsInit, posRel, bondsRel }) {
    fs.mkdirSync(plotDir, { recursive: true });
    writeCrystalJson(fs, path.join(plotDir, 'crystal_init.json'), crystalToJson({ id, label, stage: 'init', pos: posInit, Z: Zinit, bonds_ij: bondsInit }));
    writeCrystalJson(fs, path.join(plotDir, 'crystal_relaxed.json'), crystalToJson({ id, label, stage: 'relaxed', pos: posRel, Z: Zinit, bonds_ij: bondsRel }));
}

function enrichViewerJsonWithTopology(plotDir, topoNpzPath) {
    if (!fs.existsSync(topoNpzPath)) return;
    const { arrays } = readNpzFile(fs, topoNpzPath);
    const icolGroup = arrays.icolGroup;
    const bboxMin = arrays.group_bbox_min;
    const bboxMax = arrays.group_bbox_max;
    if (!icolGroup) throw new Error(`enrichViewerJsonWithTopology: missing icolGroup in ${topoNpzPath}`);
    if (!bboxMin || !bboxMax) throw new Error(`enrichViewerJsonWithTopology: missing group_bbox_min/max in ${topoNpzPath}`);
    for (const name of ['crystal_init.json', 'crystal_relaxed.json']) {
        const p = path.join(plotDir, name);
        if (!fs.existsSync(p)) continue;
        const obj = JSON.parse(fs.readFileSync(p, 'utf8'));
        const pos = Float64Array.from(obj.pos);
        const Z = Int32Array.from(obj.Z);
        const bonds = obj.bonds_ij ? Int32Array.from(obj.bonds_ij) : null;
        const extra = {
            icolGroup: Array.from(icolGroup),
            group_bbox_min: Array.from(bboxMin),
            group_bbox_max: Array.from(bboxMax),
        };
        const upd = crystalToJson({ id: obj.id, label: obj.label, stage: obj.stage, pos, Z, bonds_ij: bonds, extra });
        writeCrystalJson(fs, p, upd);
    }
}

function installViewer(outDir, crystals, title = 'Nanocrystal viewer') {
    const manifest = { title, crystals };
    fs.writeFileSync(path.join(outDir, 'viewer_manifest.json'), JSON.stringify(manifest, null, 2));
    const embedded = {
        title,
        crystals: crystals.map(c => ({
            id: c.id,
            label: c.label,
            init: JSON.parse(fs.readFileSync(path.join(outDir, c.init), 'utf8')),
            relaxed: c.relaxed ? JSON.parse(fs.readFileSync(path.join(outDir, c.relaxed), 'utf8')) : null,
        })),
    };
    const tpl = fs.readFileSync(path.join(REPO, 'web/common_js/nanocrystalViewer.html'), 'utf8');
    const html = tpl.replace('/*__NC_VIEWER_DATA__*/null', JSON.stringify(embedded));
    fs.writeFileSync(path.join(outDir, 'viewer.html'), html);
}

function writeCrystalCompareSvgs(plotDir, { posInit, Zinit, bondsInit, posRel, bondsRel, title, views }) {
    fs.mkdirSync(plotDir, { recursive: true });
    const svgs = exportCrystalCompareSvgViews({ posA: posInit, ZA: Zinit, bondsA: bondsInit, posB: posRel, ZB: Zinit, bondsB: bondsRel, views, title });
    for (const [view, svg] of Object.entries(svgs)) {
        const suffix = views.length === 1 ? '' : `_${view}`;
        fs.writeFileSync(path.join(plotDir, `compare${suffix}.svg`), svg);
    }
}

function genParamsFromAtlasEntry(atlasBase, entry) {
    return {
        cif: atlasBase.cif,
        applySymmetry: atlasBase.applySymmetry ?? 0,
        cut: entry.cut,
        replication: entry.replication || { nx: 2, ny: 2, nz: 2 },
        passivation: atlasBase.passivation,
        defects: atlasBase.defects || { insertProb: 0, collapseProb: 0 },
    };
}

async function runAtlas(atlasPath, outputDir, paths) {
    const atlas = JSON.parse(fs.readFileSync(path.resolve(atlasPath), 'utf8'));
    const atlasDir = path.join(outputDir, 'atlas');
    fs.mkdirSync(atlasDir, { recursive: true });
    const views = atlas.views || ['111'];
    const rows = [];
    const viewerCrystals = [];
    console.log(`[atlas] ${atlas.entries.length} shape entries → ${atlasDir}`);
    for (const entry of atlas.entries) {
        const id = entry.id;
        const edir = path.join(atlasDir, id);
        fs.mkdirSync(edir, { recursive: true });
        const genParams = genParamsFromAtlasEntry(atlas, entry);
        const seed = entry.seed | 0;
        const genWork = path.join(edir, 'work_gen');
        fs.mkdirSync(genWork, { recursive: true });
        console.log(`  [atlas] ${id}: ${entry.label}`);
        runCmd('node', buildGenCli({ ...genParams, cut: genParams.cut }, genWork, seed), { label: `atlas-${id}` });
        const mol2s = fs.readdirSync(genWork).filter(f => f.endsWith('.mol2')).sort();
        if (!mol2s.length) throw new Error(`atlas ${id}: no mol2`);
        const pInitMol2 = path.join(edir, 'init.mol2');
        fs.copyFileSync(path.join(genWork, mol2s[0]), pInitMol2);
        const mol = loadMolFromMol2(pInitMol2, loadMMParamsFromDir());
        const arrays = molToCrystalArrays(mol);
        const pInitXyz = path.join(edir, 'init.xyz');
        fs.writeFileSync(pInitXyz, crystalToXYZ(arrays.pos, arrays.Z));
        fs.writeFileSync(path.join(edir, 'gen_params.json'), JSON.stringify(genParams, null, 2));
        const pRelaxedXyz = path.join(edir, 'relaxed.xyz');
        runPy('relax', ['--init-xyz', pInitXyz, '--out-npz', path.join(edir, '02_relaxed.npz'), '--out-xyz', pRelaxedXyz, '--allow-unconverged'], { label: `atlas-relax-${id}`, python: paths.python });
        const bondsInit = bondsForVisualization(mol);
        const bondsRel = getCrystalBondsFromFiles(pInitMol2, pRelaxedXyz);
        const { pos: posRel } = readXyzPositions(fs, pRelaxedXyz);
        writeCrystalCompareSvgs(edir, { posInit: arrays.pos, Zinit: arrays.Z, bondsInit, posRel, bondsRel, title: `${id}: ${entry.label}`, views });
        writeCrystalViewerJson(edir, { id, label: entry.label, posInit: arrays.pos, Zinit: arrays.Z, bondsInit, posRel, bondsRel });
        viewerCrystals.push({ id, label: entry.label, init: `${id}/crystal_init.json`, relaxed: `${id}/crystal_relaxed.json` });
        const compareSvg = fs.existsSync(path.join(edir, 'compare_111.svg')) ? `${id}/compare_111.svg` : `${id}/compare.svg`;
        rows.push({ id, label: entry.label, genParams, svgRel: compareSvg, natoms: arrays.natoms });
    }
    const html = atlasIndexHtml(rows, 'C diamond nanocrystal shape atlas');
    fs.writeFileSync(path.join(atlasDir, 'index.html'), html);
    installViewer(atlasDir, viewerCrystals, 'C diamond shape atlas');
    const md = ['# C diamond nanocrystal shape atlas', '', '[Interactive 3D viewer](viewer.html)', '', '| Shape | Atoms | Preview | gen_params |', '|---|---:|---|---|'];
    for (const r of rows) {
        md.push(`| **${r.label}** (\`${r.id}\`) | ${r.natoms} | [compare](${r.svgRel}) | \`${JSON.stringify(r.genParams)}\` |`);
    }
    fs.writeFileSync(path.join(atlasDir, 'index.md'), md.join('\n') + '\n');
    fs.writeFileSync(path.join(atlasDir, 'atlas_meta.json'), JSON.stringify({ atlas_path: atlasPath, entries: rows }, null, 2));
    console.log(`[atlas] index: ${path.join(atlasDir, 'index.html')}`);
    return atlasDir;
}

async function processCrystal(index, cfg, paths, stageTimings, manifestPath, doDebug, debugViews, viewerCrystals) {
    const seed = (cfg.seed_base + index) | 0;
    const rnd = mulberry32(seed);
    const id = crystalId(seed, index);
    const cdir = path.join(paths.dataDir, 'crystals', id);
    fs.mkdirSync(cdir, { recursive: true });

    const cut = sampleCutSpec(cfg.cut_mixture, rnd);
    const rep = cfg.replication || {};
    const genParams = {
        cif: cfg.cif,
        applySymmetry: cfg.applySymmetry ?? 0,
        cut,
        replication: { nx: sampleValue(rep.nx, rnd), ny: sampleValue(rep.ny, rnd), nz: sampleValue(rep.nz, rnd) },
        passivation: cfg.passivation,
        defects: {
            insertProb: sampleValue(cfg.defects?.insertProb, rnd) ?? 0,
            collapseProb: sampleValue(cfg.defects?.collapseProb, rnd) ?? 0,
        },
    };
    const genJson = JSON.stringify(genParams);
    const metaObj = { crystal_id: id, seed, index, gen_params: genParams, bonds_ij: null };

    const p01 = path.join(cdir, '01_init.npz');
    const pInitMol2 = path.join(cdir, 'init.mol2');
    const pInitXyz = path.join(cdir, 'init.xyz');
    const p02 = path.join(cdir, '02_relaxed.npz');
    const pRelaxedXyz = path.join(cdir, 'relaxed.xyz');
    const p03 = path.join(cdir, '03_topology.npz');
    const p04 = path.join(cdir, '04_hessian.npz');
    const p05 = path.join(cdir, '05_spectrum.npz');
    const plotDir = path.join(paths.outputDir, 'crystals', id);

    console.log(`\n[ensemble] crystal ${index + 1} id=${id} cutMode=${cut.cutMode}`);

    if (!stageExists(p01) || paths.force) {
        const genWork = path.join(cdir, 'work_gen');
        fs.mkdirSync(genWork, { recursive: true });
        const t0 = performance.now();
        runCmd('node', buildGenCli({ ...genParams, cut }, genWork, seed), { label: 'generate' });
        stageTimings.generate.push(performance.now() - t0);
        const mol2s = fs.readdirSync(genWork).filter(f => f.endsWith('.mol2')).sort();
        if (!mol2s.length) throw new Error(`no mol2 in ${genWork}`);
        fs.copyFileSync(path.join(genWork, mol2s[0]), pInitMol2);
        const mm = loadMMParamsFromDir();
        const mol = loadMolFromMol2(pInitMol2, mm);
        const arrays = molToCrystalArrays(mol);
        const bondsVis = bondsForVisualization(mol);
        metaObj.bonds_ij = bondsVis ? Array.from(bondsVis) : null;
        fs.writeFileSync(path.join(cdir, 'meta.json'), JSON.stringify(metaObj, null, 2));
        writeCrystalNpz(fs, p01, { ...arrays, gen_params: genJson, timing_ms: performance.now() - t0 });
        fs.writeFileSync(pInitXyz, crystalToXYZ(arrays.pos, arrays.Z));
        appendManifest(manifestPath, { crystal_id: id, stage: 'generate', natoms: arrays.natoms, status: 'ok' });
    }

    if (!fs.existsSync(path.join(cdir, 'meta.json'))) {
        const mol = loadMolFromMol2(pInitMol2, loadMMParamsFromDir());
        const bondsVis = bondsForVisualization(mol);
        metaObj.bonds_ij = bondsVis ? Array.from(bondsVis) : null;
        fs.writeFileSync(path.join(cdir, 'meta.json'), JSON.stringify(metaObj, null, 2));
    }

    const { Z: Zinit, pos: posInit } = readXyzPositions(fs, pInitXyz);

    if (!stageExists(p02) || paths.force) {
        const r = runPy('relax', ['--init-xyz', pInitXyz, '--out-npz', p02, '--out-xyz', pRelaxedXyz, '--allow-unconverged'], { label: 'relax', statusPath: p02.replace(/\.npz$/, '.status.json'), python: paths.python });
        stageTimings.relax.push(r.timing_ms);
        appendManifest(manifestPath, { crystal_id: id, stage: 'relax', ...r.info });
    }

    if (doDebug && stageExists(pRelaxedXyz)) {
        const bondsInit = getCrystalBondsFromFiles(pInitMol2, pInitXyz);
        const bondsRel = getCrystalBondsFromFiles(pInitMol2, pRelaxedXyz);
        const { pos: posRel } = readXyzPositions(fs, pRelaxedXyz);
        writeCrystalCompareSvgs(plotDir, { posInit, Zinit, bondsInit, posRel, bondsRel, title: id, views: debugViews });
        writeCrystalViewerJson(plotDir, { id, label: id, posInit, Zinit, bondsInit, posRel, bondsRel });
        viewerCrystals.push({ id, label: id, init: `crystals/${id}/crystal_init.json`, relaxed: `crystals/${id}/crystal_relaxed.json` });
    }

    if (!stageExists(p03) || paths.force) {
        const t0 = performance.now();
        buildTopologyNpz({ mol2Path: pInitMol2, relaxedXyzPath: pRelaxedXyz, outNpzPath: p03, groupCap: cfg.group_cap ?? 32 });
        const timing_ms = performance.now() - t0;
        stageTimings.topology.push(timing_ms);
        appendManifest(manifestPath, { crystal_id: id, stage: 'topology', timing_ms, status: 'ok' });
    }

    if (doDebug && stageExists(p03)) {
        enrichViewerJsonWithTopology(plotDir, p03);
    }

    if (!stageExists(p04) || paths.force) {
        const r = runPy('hessian', ['--relaxed-xyz', pRelaxedXyz, '--out-npz', p04], { label: 'hessian', statusPath: p04.replace(/\.npz$/, '.status.json'), python: paths.python });
        stageTimings.hessian.push(r.timing_ms);
        appendManifest(manifestPath, { crystal_id: id, stage: 'hessian', ...r.info });
    }

    if (!stageExists(p05) || paths.force) {
        const spec = cfg.spectrum || {};
        const pyArgs = ['--hessian-npz', p04, '--out-npz', p05, '--margin-frac', String(spec.margin_frac ?? 0.02), '--min-bins', String(spec.min_bins ?? 64), '--sigma-bins', String(spec.sigma_bins ?? 1.5)];
        if (doDebug) pyArgs.push('--out-plot', path.join(plotDir, 'spectrum.png'));
        const r = runPy('spectrum', pyArgs, { label: 'spectrum', statusPath: p05.replace(/\.npz$/, '.status.json'), python: paths.python });
        stageTimings.spectrum.push(r.timing_ms);
        appendManifest(manifestPath, { crystal_id: id, stage: 'spectrum', ...r.info });
    }

    return p05;
}

async function main() {
    const args = parseArgs(process.argv);
    const paths = { python: args.python, force: args.force };

    if (args.atlas) {
        const outputDir = path.resolve(args.outputDir);
        fs.mkdirSync(outputDir, { recursive: true });
        await runAtlas(args.atlas, outputDir, paths);
        return;
    }

    const cfgPath = path.resolve(args.config);
    const cfg = JSON.parse(fs.readFileSync(cfgPath, 'utf8'));
    const nCrystals = args.nCrystals ?? cfg.n_crystals ?? 3;
    const dataDir = path.resolve(args.dataDir);
    const outputDir = path.resolve(args.outputDir);
    const workDir = path.resolve(args.workDir);
    fs.mkdirSync(dataDir, { recursive: true });
    fs.mkdirSync(outputDir, { recursive: true });
    fs.mkdirSync(workDir, { recursive: true });

    const doDebug = args.debug !== null ? args.debug : !!(cfg.debug?.enabled);
    const maxDebug = cfg.debug?.max_crystals ?? 999;
    const debugViews = cfg.debug?.views?.length ? cfg.debug.views : ['111', '100', '001'];

    fs.writeFileSync(path.join(dataDir, 'ensemble_meta.json'), JSON.stringify({ config: cfg, config_path: cfgPath, n_crystals: nCrystals }, null, 2));
    const manifestPath = path.join(dataDir, 'manifest.jsonl');
    if (!fs.existsSync(manifestPath)) fs.writeFileSync(manifestPath, '');

    const pathsFull = { dataDir, outputDir, workDir, python: args.python, force: args.force };
    const stageTimings = { generate: [], relax: [], topology: [], hessian: [], spectrum: [] };
    const spectrumPaths = [];
    const viewerCrystals = [];

    for (let i = 0; i < nCrystals; i++) {
        const crystalDebug = doDebug && i < maxDebug;
        spectrumPaths.push(await processCrystal(i, cfg, pathsFull, stageTimings, manifestPath, crystalDebug, debugViews, viewerCrystals));
    }

    if (viewerCrystals.length) installViewer(outputDir, viewerCrystals, `Ensemble viewer (n=${viewerCrystals.length})`);

    const accDir = path.join(outputDir, 'accumulated');
    fs.mkdirSync(accDir, { recursive: true });
    runPy('accumulate', ['--inputs', ...spectrumPaths, '--out-dir', accDir], { label: 'accumulate', python: args.python });
    writeTimingReport(outputDir, stageTimings);

    console.log(`\n[ensemble] done: ${nCrystals} crystals`);
    console.log(`  data:   ${dataDir}`);
    console.log(`  output: ${outputDir}`);
    console.log(`  plots:  ${path.join(outputDir, 'plots')}/spectrum_histogram.png (+ pol_x/y/z)`);
    if (viewerCrystals.length) console.log(`  viewer: ${path.join(outputDir, 'viewer.html')}`);
}

main().catch((e) => { console.error(e); process.exit(1); });
