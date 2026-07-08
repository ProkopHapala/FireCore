#!/usr/bin/env node
/// Unified nanocrystal CLI — consolidates gen_nanocrystals, run_nanocrystal_ensemble,
/// build_linearized_topology, test_nanocrystal_geometry, debug_nanocrystal_nonbond_groups,
/// and test_ring_detection into a single entrypoint with subcommands.
/// Old scripts are deprecated but not deleted until this passes equivalent tests.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
 
const __dirname = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(__dirname, '..');
const RES = path.join(REPO, 'cpp/common_resources');
const VIEWER_TPL = path.join(REPO, 'web/common_js/nanocrystalViewer.html');
 
import { toMol2String, toXYZString } from '../web/molgui_webgpu/MoleculeIO.js';
import {
    loadMMParams, buildCrystalFromCIFText, getHeavyZ, generateNanocrystal, mulberry32, defaultArgs,
    auditGeometry, oneLineSummary, buildSiH4, isBulkSi,
    sampleValue, sampleCutSpec, crystalId, appendManifest, runCmd, runPy, writeTimingReport,
    genArgsFromConfig, genParamsFromAtlasEntry,
} from '../web/molgui_webgpu/Nanocrystals.js';
import { loadMMParamsFromDir, loadMolFromMol2, loadMol, applyPositions, bondsForVisualization, getCrystalBondsFromFiles } from '../web/common_js/MolIO.js';
import { molToCrystalArrays, writeCrystalNpz, crystalToXYZ, readXyzPositions, crystalToJson, writeCrystalJson, readNpzFile } from '../web/common_js/npzIO.js';
import { exportCrystalSvg, exportCrystalCompareSvgViews, bondsForSvgMol, atlasIndexHtml, writeCrystalCompareSvgs, writeCrystalViewerJson, enrichViewerJsonWithTopology, installViewer } from '../web/common_js/nanocrystalSvg.js';
import { buildTopologyNpz, buildTopologyFull } from '../web/common_js/exportFF.js';
import { runRingsOnMol } from '../web/common_js/Graph.js';
import { buildCollisionWorkgroups, buildExclIcol_1_2_3 } from '../web/common_js/CollisionWorkgroups.js';
import { computeNonbondBruteForceKernelStyle, computeNonbondByGroups, maxAbsDiff } from '../web/common_js/Nonbonded.js';
 
// ============================================================================
// Shared utilities
// ============================================================================
 
function loadJSON(p) { return JSON.parse(fs.readFileSync(p, 'utf8')); }
 
function mergeConfig(jsonCfg, cliOverrides) {
    if (!jsonCfg) return cliOverrides;
    if (!cliOverrides) return jsonCfg;
    const out = { ...jsonCfg };
    for (const [k, v] of Object.entries(cliOverrides)) if (v !== undefined && v !== null) out[k] = v;
    return out;
}
 
function parseProb(s, name) {
    const v = +String(s).trim();
    if (!isFinite(v) || v < 0 || v > 1) throw new Error(`${name}: expected probability in [0,1], got '${s}'`);
    return v;
}
 
function parseCSVRange(s) {
    const t = String(s || '').trim();
    if (!t) throw new Error('parseCSVRange: empty');
    const parts = t.split(',').map(x => x.trim()).filter(x => x.length > 0);
    if (parts.length === 1) { const v = parseInt(parts[0], 10); return [v | 0, v | 0]; }
    if (parts.length === 2) { return [parseInt(parts[0], 10) | 0, parseInt(parts[1], 10) | 0]; }
    throw new Error(`parseCSVRange: expected 'n' or 'a,b', got '${s}'`);
}
 
function usage() {
    return `Usage: node scripts/nanocrystal.mjs <subcommand> [options]
 
Subcommands:
  generate   — single crystal generation (mol2+xyz output)
  ensemble   — batch ensemble pipeline (gen→relax→topo→hessian→spectrum→accumulate)
  topology   — linearized topology build (mol2/xyz/npz → topology.npz + viewer)
  rings      — ring detection on a crystal (SVG output)
  audit      — geometry audit (bond lengths, angles, clashes, valence)
  nonbond    — nonbond group debug + collision pair export
 
Global options:
  --config PATH          JSON config file (CLI flags override JSON values)
  --help, -h             show this message
 
Run 'node scripts/nanocrystal.mjs <subcommand> --help' for subcommand-specific options.`;
}
 
// ============================================================================
// generate subcommand
// ============================================================================
 
function generateUsage() {
    return `Usage: node scripts/nanocrystal.mjs generate [options]
 
Options:
  --cif PATH             CIF file (default: cpp/common_resources/crystals/Si-sym.cif)
  --cutMode MODE         'planes' or 'sphere' (default: planes)
  --sphereR R            sphere radius for cutMode=sphere
  --sphereNrep N         sphere replication count
  --nx-range A,B         cell replication range for cutMode=planes
  --ny-range A,B         cell replication range
  --nz-range A,B         cell replication range
  --planeTemplates T     comma-separated plane templates (default: a111)
  --planeSymC N          plane symmetry count (default: 6)
  --planeCScale F        plane cutoff scale (default: 0.55)
  --planeCJitter F       plane cutoff jitter (default: 0.10)
  --caps S               cap atom type (default: H)
  --insertProb P         probability of inserting -CH2- bridges
  --collapseProb P       probability of collapsing bridges
  --outwardBias F        outward bias for cap placement (default: 0.35)
  --resolveClashes 0|1   resolve steric clashes (default: 1)
  --seed N               random seed (0 = Math.random)
  --samples N            number of crystals to generate
  --outDir DIR           output directory
  --prefix S             filename prefix
  --config PATH          JSON config file`;
}
 
function cmdGenerate(argv) {
    const cfg = { cif: null, cutMode: null, sphereR: null, sphereNrep: null, nx: null, ny: null, nz: null, planeTemplates: null, planeSymC: null, planeCScale: null, planeCJitter: null, caps: null, insertProb: null, collapseProb: null, outwardBias: null, resolveClashes: null, seed: null, samples: null, outDir: null, prefix: null, config: null };
    for (let i = 0; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => { if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`); return String(argv[++i]); };
        if (a === '--help' || a === '-h') { console.log(generateUsage()); process.exit(0); }
        else if (a === '--config') cfg.config = nxt();
        else if (a === '--cif') cfg.cif = path.resolve(nxt());
        else if (a === '--cutMode') cfg.cutMode = nxt();
        else if (a === '--sphereR') cfg.sphereR = +nxt();
        else if (a === '--sphereNrep') cfg.sphereNrep = parseInt(nxt(), 10) | 0;
        else if (a === '--nx-range' || a === '--nx') cfg.nx = parseCSVRange(nxt());
        else if (a === '--ny-range' || a === '--ny') cfg.ny = parseCSVRange(nxt());
        else if (a === '--nz-range' || a === '--nz') cfg.nz = parseCSVRange(nxt());
        else if (a === '--planeTemplates') cfg.planeTemplates = nxt().split(',').map(s => s.trim()).filter(s => s.length);
        else if (a === '--planeSymC') cfg.planeSymC = +nxt();
        else if (a === '--planeCScale') cfg.planeCScale = +nxt();
        else if (a === '--planeCJitter') cfg.planeCJitter = +nxt();
        else if (a === '--caps') cfg.caps = nxt();
        else if (a === '--insertProb') cfg.insertProb = parseProb(nxt(), '--insertProb');
        else if (a === '--collapseProb') cfg.collapseProb = parseProb(nxt(), '--collapseProb');
        else if (a === '--outwardBias') cfg.outwardBias = +nxt();
        else if (a === '--resolveClashes') cfg.resolveClashes = (nxt() !== '0');
        else if (a === '--seed') cfg.seed = parseInt(nxt(), 10) | 0;
        else if (a === '--samples') cfg.samples = parseInt(nxt(), 10) | 0;
        else if (a === '--outDir') cfg.outDir = path.resolve(nxt());
        else if (a === '--prefix') cfg.prefix = nxt();
        else throw new Error(`Unknown arg: ${a}\n\n${generateUsage()}`);
    }
    let jsonCfg = null;
    if (cfg.config) jsonCfg = loadJSON(cfg.config);
    const merged = mergeConfig(jsonCfg, cfg);
    // Build nested config structure from flat CLI args for genArgsFromConfig
    const nested = { cif: merged.cif, applySymmetry: merged.applySymmetry ?? true };
    nested.cut = {};
    if (merged.cutMode) nested.cut.cutMode = merged.cutMode;
    if (merged.sphereR != null) nested.cut.sphereR = merged.sphereR;
    if (merged.sphereNrep != null) nested.cut.sphereNrep = merged.sphereNrep;
    if (merged.planeTemplates) nested.cut.planeTemplates = merged.planeTemplates;
    if (merged.planeSymC != null) nested.cut.planeSymC = merged.planeSymC;
    if (merged.planeCScale != null) nested.cut.planeCScale = merged.planeCScale;
    if (merged.planeCJitter != null) nested.cut.planeCJitter = merged.planeCJitter;
    if (merged.nx || merged.ny || merged.nz) nested.replication = { nx: merged.nx?.[0], ny: merged.ny?.[0], nz: merged.nz?.[0] };
    nested.passivation = {};
    if (merged.caps) nested.passivation.caps = merged.caps;
    if (merged.outwardBias != null) nested.passivation.outwardBias = merged.outwardBias;
    if (merged.resolveClashes != null) nested.passivation.resolveClashes = merged.resolveClashes;
    nested.defects = {};
    if (merged.insertProb != null) nested.defects.insertProb = merged.insertProb;
    if (merged.collapseProb != null) nested.defects.collapseProb = merged.collapseProb;
    // Merge with jsonCfg if it has nested structure
    if (jsonCfg) {
        if (jsonCfg.cut) nested.cut = { ...jsonCfg.cut, ...nested.cut };
        if (jsonCfg.replication && !nested.replication) nested.replication = jsonCfg.replication;
        if (jsonCfg.passivation) nested.passivation = { ...jsonCfg.passivation, ...nested.passivation };
        if (jsonCfg.defects) nested.defects = { ...jsonCfg.defects, ...nested.defects };
    }
    const args = genArgsFromConfig(nested, RES, { outDir: merged.outDir || path.join(REPO, 'OUT_nanocrystals'), prefix: merged.prefix || 'nc', seed: merged.seed ?? 42, samples: merged.samples ?? 1 });
    if (args.seed !== 0) { const rnd = mulberry32(args.seed); Math.random = rnd; }
    fs.mkdirSync(args.outDir, { recursive: true });
    const cifText = fs.readFileSync(args.cif, 'utf8');
    const cell = buildCrystalFromCIFText(cifText, args);
    const heavyZ = getHeavyZ(cell);
    const mm = loadMMParams(args);
    const nGen = args.samples > 0 ? args.samples : 1;
    for (let iout = 1; iout <= nGen; iout++) {
        const result = generateNanocrystal(args, cell, mm, heavyZ, iout);
        const { mol, name, nCollapsed, nInserted, nCaps, nPruned, nHHBonds, cnt, enr } = result;
        const mol2 = toMol2String(mol, { name });
        const outPath = path.join(args.outDir, name + '.mol2');
        fs.writeFileSync(outPath, mol2);
        fs.writeFileSync(path.join(args.outDir, name + '.xyz'), toXYZString(mol, { lvec: mol.lvec }));
        console.log(`[generate] wrote ${outPath} atoms=${mol.atoms.length} bonds=${mol.bonds.length} pruned=${nPruned} caps=${nCaps} hhBonds=${nHHBonds} collapsed=${nCollapsed} inserted=${nInserted} E=${enr.E.toFixed(6)}`);
    }
}
 
// ============================================================================
// ensemble subcommand
// ============================================================================
 
function ensembleUsage() {
    return `Usage: node scripts/nanocrystal.mjs ensemble [options]
 
Options:
  --config PATH          ensemble JSON (default: scripts/ensemble.example.json)
  --data-dir PATH        per-crystal NPZ cache (default: OUT_nc_ensemble/data)
  --output-dir PATH      accumulated outputs (default: OUT_nc_ensemble/out)
  --work-dir PATH        scratch (default: OUT_nc_ensemble/work)
  --n-crystals N         override config n_crystals
  --python PATH          Python executable (default: python3 or $PYTHON)
  --force                re-run all stages
  --no-debug             skip per-crystal SVG debug plots
  --atlas PATH           shape showcase atlas (JSON entries, no full ensemble)`;
}
 
async function cmdEnsemble(argv) {
    const cfg = { config: path.join(REPO, 'scripts/ensemble.example.json'), dataDir: path.join(REPO, 'OUT_nc_ensemble/data'), outputDir: path.join(REPO, 'OUT_nc_ensemble/out'), workDir: path.join(REPO, 'OUT_nc_ensemble/work'), nCrystals: null, python: process.env.PYTHON || 'python3', force: false, debug: null, atlas: null };
    for (let i = 0; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => { if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`); return String(argv[++i]); };
        if (a === '--help' || a === '-h') { console.log(ensembleUsage()); process.exit(0); }
        else if (a === '--config') cfg.config = nxt();
        else if (a === '--data-dir') cfg.dataDir = nxt();
        else if (a === '--output-dir') cfg.outputDir = nxt();
        else if (a === '--work-dir') cfg.workDir = nxt();
        else if (a === '--n-crystals') cfg.nCrystals = parseInt(nxt(), 10);
        else if (a === '--python') cfg.python = nxt();
        else if (a === '--force') cfg.force = true;
        else if (a === '--no-debug') cfg.debug = false;
        else if (a === '--atlas') cfg.atlas = nxt();
        else throw new Error(`Unknown arg: ${a}\n\n${ensembleUsage()}`);
    }
    const paths = { python: cfg.python, force: cfg.force, repoRoot: REPO };
 
    if (cfg.atlas) {
        const outputDir = path.resolve(cfg.outputDir);
        fs.mkdirSync(outputDir, { recursive: true });
        await runAtlas(cfg.atlas, outputDir, paths);
        return;
    }
 
    const cfgObj = loadJSON(path.resolve(cfg.config));
    const nCrystals = cfg.nCrystals ?? cfgObj.n_crystals ?? 3;
    const dataDir = path.resolve(cfg.dataDir);
    const outputDir = path.resolve(cfg.outputDir);
    fs.mkdirSync(dataDir, { recursive: true });
    fs.mkdirSync(outputDir, { recursive: true });
    const doDebug = cfg.debug !== null ? cfg.debug : !!(cfgObj.debug?.enabled);
    const maxDebug = cfgObj.debug?.max_crystals ?? 999;
    const debugViews = cfgObj.debug?.views?.length ? cfgObj.debug.views : ['111', '100', '001'];
    fs.writeFileSync(path.join(dataDir, 'ensemble_meta.json'), JSON.stringify({ config: cfgObj, config_path: cfg.config, n_crystals: nCrystals }, null, 2));
    const manifestPath = path.join(dataDir, 'manifest.jsonl');
    if (!fs.existsSync(manifestPath)) fs.writeFileSync(manifestPath, '');
    const stageTimings = { generate: [], relax: [], topology: [], hessian: [], spectrum: [] };
    const spectrumPaths = [];
    const viewerCrystals = [];
    for (let i = 0; i < nCrystals; i++) {
        const crystalDebug = doDebug && i < maxDebug;
        spectrumPaths.push(await processCrystal(i, cfgObj, { dataDir, outputDir, ...paths }, stageTimings, manifestPath, crystalDebug, debugViews, viewerCrystals));
    }
    if (viewerCrystals.length) installViewer(outputDir, viewerCrystals, VIEWER_TPL, `Ensemble viewer (n=${viewerCrystals.length})`);
    const accDir = path.join(outputDir, 'accumulated');
    fs.mkdirSync(accDir, { recursive: true });
    runPy('accumulate', ['--inputs', ...spectrumPaths, '--out-dir', accDir], { label: 'accumulate', python: cfg.python, cwd: REPO, repoRoot: REPO });
    writeTimingReport(outputDir, stageTimings);
    console.log(`\n[ensemble] done: ${nCrystals} crystals`);
    console.log(`  data:   ${dataDir}`);
    console.log(`  output: ${outputDir}`);
    if (viewerCrystals.length) console.log(`  viewer: ${path.join(outputDir, 'viewer.html')}`);
}
 
async function processCrystal(index, cfg, paths, stageTimings, manifestPath, doDebug, debugViews, viewerCrystals) {
    const seed = (cfg.seed_base + index) | 0;
    const rnd = mulberry32(seed);
    const id = crystalId(seed, index);
    const cdir = path.join(paths.dataDir, 'crystals', id);
    fs.mkdirSync(cdir, { recursive: true });
    const cut = sampleCutSpec(cfg.cut_mixture, rnd);
    const rep = cfg.replication || {};
    const genParams = { cif: cfg.cif, applySymmetry: cfg.applySymmetry ?? 0, cut, replication: { nx: sampleValue(rep.nx, rnd), ny: sampleValue(rep.ny, rnd), nz: sampleValue(rep.nz, rnd) }, passivation: cfg.passivation, defects: { insertProb: sampleValue(cfg.defects?.insertProb, rnd) ?? 0, collapseProb: sampleValue(cfg.defects?.collapseProb, rnd) ?? 0 } };
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
 
    if (!fs.existsSync(p01) || paths.force) {
        const args = genArgsFromConfig({ ...genParams, cut }, RES, { seed, samples: 1, maxFiles: 1 });
        if (args.seed !== 0) { const r = mulberry32(args.seed); Math.random = r; }
        const cifText = fs.readFileSync(args.cif, 'utf8');
        const cell = buildCrystalFromCIFText(cifText, args);
        const heavyZ = getHeavyZ(cell);
        const mm = loadMMParams(args);
        const t0 = performance.now();
        const result = generateNanocrystal(args, cell, mm, heavyZ, 1);
        stageTimings.generate.push(performance.now() - t0);
        const { mol } = result;
        const mol2 = toMol2String(mol, { name: 'nc' });
        fs.writeFileSync(pInitMol2, mol2);
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
 
    if (!fs.existsSync(p02) || paths.force) {
        const r = runPy('relax', ['--init-xyz', pInitXyz, '--out-npz', p02, '--out-xyz', pRelaxedXyz, '--allow-unconverged'], { label: 'relax', statusPath: p02.replace(/\.npz$/, '.status.json'), python: paths.python, cwd: REPO, repoRoot: REPO });
        stageTimings.relax.push(r.timing_ms);
        appendManifest(manifestPath, { crystal_id: id, stage: 'relax', ...r.info });
    }
    if (doDebug && fs.existsSync(pRelaxedXyz)) {
        const bondsInit = getCrystalBondsFromFiles(pInitMol2, pInitXyz);
        const bondsRel = getCrystalBondsFromFiles(pInitMol2, pRelaxedXyz);
        const { pos: posRel } = readXyzPositions(fs, pRelaxedXyz);
        writeCrystalCompareSvgs(plotDir, { posInit, Zinit, bondsInit, posRel, bondsRel, title: id, views: debugViews });
        writeCrystalViewerJson(plotDir, { id, label: id, posInit, Zinit, bondsInit, posRel, bondsRel });
        viewerCrystals.push({ id, label: id, init: `crystals/${id}/crystal_init.json`, relaxed: `crystals/${id}/crystal_relaxed.json` });
    }
    if (!fs.existsSync(p03) || paths.force) {
        const t0 = performance.now();
        buildTopologyNpz({ mol2Path: pInitMol2, relaxedXyzPath: pRelaxedXyz, outNpzPath: p03, groupCap: cfg.group_cap ?? 32 });
        stageTimings.topology.push(performance.now() - t0);
        appendManifest(manifestPath, { crystal_id: id, stage: 'topology', timing_ms: performance.now() - t0, status: 'ok' });
    }
    if (doDebug && fs.existsSync(p03)) enrichViewerJsonWithTopology(plotDir, p03);
    if (!fs.existsSync(p04) || paths.force) {
        const r = runPy('hessian', ['--relaxed-xyz', pRelaxedXyz, '--out-npz', p04], { label: 'hessian', statusPath: p04.replace(/\.npz$/, '.status.json'), python: paths.python, cwd: REPO, repoRoot: REPO });
        stageTimings.hessian.push(r.timing_ms);
        appendManifest(manifestPath, { crystal_id: id, stage: 'hessian', ...r.info });
    }
    if (!fs.existsSync(p05) || paths.force) {
        const spec = cfg.spectrum || {};
        const pyArgs = ['--hessian-npz', p04, '--out-npz', p05, '--margin-frac', String(spec.margin_frac ?? 0.02), '--min-bins', String(spec.min_bins ?? 64), '--sigma-bins', String(spec.sigma_bins ?? 1.5)];
        if (doDebug) pyArgs.push('--out-plot', path.join(plotDir, 'spectrum.png'));
        const r = runPy('spectrum', pyArgs, { label: 'spectrum', statusPath: p05.replace(/\.npz$/, '.status.json'), python: paths.python, cwd: REPO, repoRoot: REPO });
        stageTimings.spectrum.push(r.timing_ms);
        appendManifest(manifestPath, { crystal_id: id, stage: 'spectrum', ...r.info });
    }
    return p05;
}

async function runAtlas(atlasPath, outputDir, paths) {
    const atlas = loadJSON(path.resolve(atlasPath));
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
        const args = genArgsFromConfig(genParams, RES, { seed, samples: 1, maxFiles: 1 });
        if (args.seed !== 0) { const r = mulberry32(args.seed); Math.random = r; }
        const cifText = fs.readFileSync(args.cif, 'utf8');
        const cell = buildCrystalFromCIFText(cifText, args);
        const heavyZ = getHeavyZ(cell);
        const mm = loadMMParams(args);
        console.log(`  [atlas] ${id}: ${entry.label}`);
        const result = generateNanocrystal(args, cell, mm, heavyZ, 1);
        const { mol } = result;
        const pInitMol2 = path.join(edir, 'init.mol2');
        fs.writeFileSync(pInitMol2, toMol2String(mol, { name: id }));
        const arrays = molToCrystalArrays(mol);
        const pInitXyz = path.join(edir, 'init.xyz');
        fs.writeFileSync(pInitXyz, crystalToXYZ(arrays.pos, arrays.Z));
        fs.writeFileSync(path.join(edir, 'gen_params.json'), JSON.stringify(genParams, null, 2));
        const pRelaxedXyz = path.join(edir, 'relaxed.xyz');
        runPy('relax', ['--init-xyz', pInitXyz, '--out-npz', path.join(edir, '02_relaxed.npz'), '--out-xyz', pRelaxedXyz, '--allow-unconverged'], { label: `atlas-relax-${id}`, python: paths.python, cwd: REPO, repoRoot: REPO });
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
    installViewer(atlasDir, viewerCrystals, VIEWER_TPL, 'C diamond shape atlas');
    console.log(`[atlas] index: ${path.join(atlasDir, 'index.html')}`);
}

// ============================================================================
// topology subcommand
// ============================================================================

function topologyUsage() {
    return `Usage: node scripts/nanocrystal.mjs topology [options]

Options:
  --input PATH           .mol2 or .xyz file
  --out-dir DIR          output directory
  --add-angle 0|1        include 1-3 angle springs (default: 1)
  --add-dihedral 0|1     include 1-4 dihedral springs (default: 1)
  --K12 F                override bond spring constant
  --K13 F                override angle spring constant
  --K14 F                override dihedral spring constant
  --max-neighbors N      max neighbors per atom (default: 48)
  --export-npz 0|1       export NPZ (default: 1)
  --export-json 0|1      export JSON spring list (default: 0)
  --html-viewer 0|1      export HTML viewer (default: 1)
  --npz-viewer 0|1       export NPZ viewer HTML (default: 1)
  --graph-dist 0|1       validate graph distances (default: 1)`;
}

async function cmdTopology(argv) {
    const cfg = { input: null, outDir: path.join(REPO, 'tests/tSiNCs/linearized'), addAngle: true, addDihedral: true, K12: 0, K13: 0, K14: 0, maxNeighbors: 48, exportNpz: true, exportJson: false, htmlViewer: true, npzViewer: true, graphDist: true, config: null };
    for (let i = 0; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => { if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`); return String(argv[++i]); };
        if (a === '--help' || a === '-h') { console.log(topologyUsage()); process.exit(0); }
        else if (a === '--config') cfg.config = nxt();
        else if (a === '--input') cfg.input = path.resolve(nxt());
        else if (a === '--out-dir') cfg.outDir = path.resolve(nxt());
        else if (a === '--add-angle') cfg.addAngle = (nxt() !== '0');
        else if (a === '--add-dihedral') cfg.addDihedral = (nxt() !== '0');
        else if (a === '--K12') cfg.K12 = +nxt();
        else if (a === '--K13') cfg.K13 = +nxt();
        else if (a === '--K14') cfg.K14 = +nxt();
        else if (a === '--max-neighbors') cfg.maxNeighbors = parseInt(nxt(), 10) | 0;
        else if (a === '--export-npz') cfg.exportNpz = (nxt() !== '0');
        else if (a === '--export-json') cfg.exportJson = (nxt() !== '0');
        else if (a === '--html-viewer') cfg.htmlViewer = (nxt() !== '0');
        else if (a === '--npz-viewer') cfg.npzViewer = (nxt() !== '0');
        else if (a === '--graph-dist') cfg.graphDist = (nxt() !== '0');
        else throw new Error(`Unknown arg: ${a}\n\n${topologyUsage()}`);
    }
    let jsonCfg = null;
    if (cfg.config) jsonCfg = loadJSON(cfg.config);
    const merged = mergeConfig(jsonCfg, cfg);
    if (!merged.input) throw new Error('--input is required\n\n' + topologyUsage());
    fs.mkdirSync(merged.outDir, { recursive: true });
    const results = await buildTopologyFull({ ...merged, inputPath: merged.input });
    console.log(`[topology] ${results.base}: atoms=${results.topo.n_real} bonds=${results.nBond} angle_springs=${results.nAng} dihedral_springs=${results.nDih}`);
    if (results.npzPath) console.log(`  NPZ: ${results.npzPath}`);
    if (results.jsonPath) console.log(`  JSON: ${results.jsonPath}`);
    if (results.htmlPath) console.log(`  HTML: ${results.htmlPath}`);
    if (results.npzHtmlPath) console.log(`  NPZ viewer: ${results.npzHtmlPath}`);
}

// ============================================================================
// rings subcommand
// ============================================================================

function ringsUsage() {
    return `Usage: node scripts/nanocrystal.mjs rings [options]

Options:
  --input PATH           .mol2 file
  --out-dir DIR          output directory for SVG
  --label S              label for SVG title`;
}

function cmdRings(argv) {
    const cfg = { input: null, outDir: path.join(REPO, 'tests/tSiNCs/rings'), label: 'rings', config: null };
    for (let i = 0; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => { if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`); return String(argv[++i]); };
        if (a === '--help' || a === '-h') { console.log(ringsUsage()); process.exit(0); }
        else if (a === '--config') cfg.config = nxt();
        else if (a === '--input') cfg.input = path.resolve(nxt());
        else if (a === '--out-dir') cfg.outDir = path.resolve(nxt());
        else if (a === '--label') cfg.label = nxt();
        else throw new Error(`Unknown arg: ${a}\n\n${ringsUsage()}`);
    }
    let jsonCfg = null;
    if (cfg.config) jsonCfg = loadJSON(cfg.config);
    const merged = mergeConfig(jsonCfg, cfg);
    if (!merged.input) throw new Error('--input is required\n\n' + ringsUsage());
    fs.mkdirSync(merged.outDir, { recursive: true });
    const mm = loadMMParamsFromDir();
    const mol = loadMol(merged.input, mm);
    const ringData = runRingsOnMol(mol, {});
    const { cls, result, oldIdx } = ringData;
    console.log(`[rings] heavy atoms: ${ringData.nHeavy} rings: ${cls.nRings} (${cls.summary})`);
    for (const dr of cls.defectRings || []) {
        const molAtoms = dr.atoms.map(i => oldIdx[i]);
        console.log(`  defect ring size=${dr.size} atoms=[${molAtoms.join(',')}]`);
    }
    const svgData = bondsForSvgMol(mol);
    const ringOfBond = new Map();
    const rings = result.rings || [];
    for (let ri = 0; ri < rings.length; ri++) {
        const ring = rings[ri];
        for (let k = 0; k < ring.length; k++) {
            const a = ring[k], b = ring[(k + 1) % ring.length];
            const key = a < b ? `${a},${b}` : `${b},${a}`;
            if (!ringOfBond.has(key)) ringOfBond.set(key, []);
            ringOfBond.get(key).push(ri);
        }
    }
    const svg = exportCrystalSvg({ pos: svgData.pos, Z: svgData.Z, bonds_ij: svgData.bonds_ij, view: '111', title: merged.label, ringOfBond, rings });
    const svgPath = path.join(merged.outDir, `${merged.label}_rings.svg`);
    fs.writeFileSync(svgPath, svg);
    console.log(`[rings] wrote ${svgPath}`);
}

// ============================================================================
// audit subcommand
// ============================================================================

function auditUsage() {
    return `Usage: node scripts/nanocrystal.mjs audit [options]

Options:
  --input PATH           .mol2 or .xyz file to audit
  --preset NAME          built-in preset (G0-G5)
  --all-presets          run all presets
  --out-dir DIR          output directory for reports
  --out-json PATH        write single report to PATH
  --bondFactor F         bond recalculation factor (default: 1.30)
  --defaultRcut F        default bond cutoff (default: 1.60)
  --no-fail              don't exit non-zero on required preset failures`;
}

const AUDIT_PRESETS = {
    G0: { kind: 'build', build: 'sih4' },
    G1: { kind: 'gen', genArgs: ['--nx-range', '1,1', '--ny-range', '1,1', '--nz-range', '1,1', '--centered', '1', '--planeTemplates', 'a111', '--planeSymC', '6', '--planeCScale', '0.50', '--planeCJitter', '0', '--minSiDegree', '2'] },
    G2: { kind: 'gen', genArgs: ['--nx-range', '2,2', '--ny-range', '2,2', '--nz-range', '2,2', '--centered', '1', '--planeTemplates', 'a111', '--planeSymC', '6', '--planeCScale', '0.40', '--planeCJitter', '0', '--minSiDegree', '2'] },
    G3: { kind: 'gen', genArgs: ['--nx-range', '2,2', '--ny-range', '2,2', '--nz-range', '2,2', '--centered', '1', '--planeTemplates', 'a111', '--planeSymC', '6', '--planeCScale', '0.40', '--planeCJitter', '0', '--insertProb', '0.2', '--collapseProb', '0.3'] },
    G4: { kind: 'file', mol2: path.join(REPO, 'tests/tSiNCs/fixtures/vibration_parallel/structures/adamantane.mol2') },
    G5: { kind: 'file', xyz: path.join(REPO, 'tests/tSiNCs/fixtures/vibration_parallel/structures/diamond_primitive.xyz') },
};

async function cmdAudit(argv) {
    const cfg = { preset: null, input: null, outJson: null, outDir: path.join(REPO, 'tests/tSiNCs/geometry'), bondFactor: 1.30, defaultRcut: 1.60, failOnError: true, allPresets: false, config: null };
    for (let i = 0; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => { if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`); return String(argv[++i]); };
        if (a === '--help' || a === '-h') { console.log(auditUsage()); process.exit(0); }
        else if (a === '--config') cfg.config = nxt();
        else if (a === '--preset') cfg.preset = nxt().toUpperCase();
        else if (a === '--input') cfg.input = path.resolve(nxt());
        else if (a === '--out-json') cfg.outJson = path.resolve(nxt());
        else if (a === '--out-dir') cfg.outDir = path.resolve(nxt());
        else if (a === '--bondFactor') cfg.bondFactor = +nxt();
        else if (a === '--defaultRcut') cfg.defaultRcut = +nxt();
        else if (a === '--no-fail') cfg.failOnError = false;
        else if (a === '--all-presets') cfg.allPresets = true;
        else throw new Error(`Unknown arg: ${a}\n\n${auditUsage()}`);
    }
    let jsonCfg = null;
    if (cfg.config) jsonCfg = loadJSON(cfg.config);
    const merged = mergeConfig(jsonCfg, cfg);
    const mm = loadMMParamsFromDir();
    const bondOpts = { bondFactor: merged.bondFactor, defaultRcut: merged.defaultRcut };
    const results = [];
    const runOne = (name, mol) => {
        if (mol.bonds.length === 0) mol.recalculateBonds(mm, bondOpts);
        const report = auditGeometry(mol, mm);
        report.preset = name;
        const line = oneLineSummary(name, report);
        console.log(line);
        fs.mkdirSync(merged.outDir, { recursive: true });
        const jsonPath = path.join(merged.outDir, `${name.toLowerCase()}_geometry_report.json`);
        fs.writeFileSync(jsonPath, JSON.stringify(report, null, 2));
        const xyzPath = path.join(merged.outDir, `${name.toLowerCase()}.xyz`);
        fs.writeFileSync(xyzPath, toXYZString(mol, {}));
        return { report, jsonPath, line };
    };
    const loadMolForAudit = (p) => loadMol(p, mm);
    const runGenPreset = (genArgs, outDir) => {
        fs.mkdirSync(outDir, { recursive: true });
        // Parse preset genArgs (old CLI format) into nested config
        const presetCfg = { applySymmetry: true, cut: { cutMode: 'planes', planeTemplates: ['a111'], planeSymC: 6, planeCScale: 0.50, planeCJitter: 0 }, replication: { nx: 1, ny: 1, nz: 1 }, passivation: { caps: 'H' }, defects: { insertProb: 0, collapseProb: 0 } };
        for (let i = 0; i < genArgs.length; i++) {
            const a = genArgs[i];
            if (a === '--nx-range' || a === '--nx') { const v = genArgs[++i].split(',').map(x => parseInt(x, 10)); presetCfg.replication.nx = v[0]; }
            else if (a === '--ny-range' || a === '--ny') { const v = genArgs[++i].split(',').map(x => parseInt(x, 10)); presetCfg.replication.ny = v[0]; }
            else if (a === '--nz-range' || a === '--nz') { const v = genArgs[++i].split(',').map(x => parseInt(x, 10)); presetCfg.replication.nz = v[0]; }
            else if (a === '--planeTemplates') presetCfg.cut.planeTemplates = genArgs[++i].split(',');
            else if (a === '--planeSymC') presetCfg.cut.planeSymC = +genArgs[++i];
            else if (a === '--planeCScale') presetCfg.cut.planeCScale = +genArgs[++i];
            else if (a === '--planeCJitter') presetCfg.cut.planeCJitter = +genArgs[++i];
            else if (a === '--insertProb') presetCfg.defects.insertProb = +genArgs[++i];
            else if (a === '--collapseProb') presetCfg.defects.collapseProb = +genArgs[++i];
            else if (a === '--centered') { i++; /* skip value */ }
            else if (a === '--minSiDegree') { i++; /* skip, handled by genArgsFromConfig */ }
        }
        const args = genArgsFromConfig(presetCfg, RES, { seed: 42, samples: 1, maxFiles: 1 });
        const rnd = mulberry32(42); Math.random = rnd;
        const cifText = fs.readFileSync(args.cif, 'utf8');
        const cell = buildCrystalFromCIFText(cifText, args);
        const heavyZ = getHeavyZ(cell);
        const mm2 = loadMMParams(args);
        const result = generateNanocrystal(args, cell, mm2, heavyZ, 1);
        const mol2Path = path.join(outDir, 'audit_0001.mol2');
        fs.writeFileSync(mol2Path, toMol2String(result.mol, { name: 'audit' }));
        return mol2Path;
    };
    if (merged.allPresets) {
        for (const key of Object.keys(AUDIT_PRESETS)) {
            const p = AUDIT_PRESETS[key];
            let mol;
            if (p.kind === 'build') mol = buildSiH4();
            else if (p.kind === 'file') {
                const fp = p.mol2 || p.xyz;
                if (!fs.existsSync(fp)) { console.log(`[${key}] SKIP missing ${fp}`); continue; }
                mol = loadMolForAudit(fp);
            } else if (p.kind === 'gen') {
                const genOut = path.join(merged.outDir, `gen_${key.toLowerCase()}`);
                const mol2 = runGenPreset(p.genArgs, genOut);
                mol = loadMolForAudit(mol2);
            }
            results.push(runOne(key, mol));
        }
    } else if (merged.preset) {
        const p = AUDIT_PRESETS[merged.preset];
        if (!p) throw new Error(`Unknown preset ${merged.preset}; known: ${Object.keys(AUDIT_PRESETS).join(', ')}`);
        let mol;
        if (p.kind === 'build') mol = buildSiH4();
        else if (p.kind === 'file') mol = loadMolForAudit(p.mol2 || p.xyz);
        else if (p.kind === 'gen') {
            const genOut = path.join(merged.outDir, `gen_${merged.preset.toLowerCase()}`);
            mol = loadMolForAudit(runGenPreset(p.genArgs, genOut));
        }
        results.push(runOne(merged.preset, mol));
    } else if (merged.input) {
        results.push(runOne(path.basename(merged.input), loadMolForAudit(merged.input)));
    } else {
        console.error('usage: node scripts/nanocrystal.mjs audit --preset G2 | --input file.mol2 | --all-presets [--out-dir DIR]');
        process.exit(1);
    }
    if (merged.outJson && results.length === 1) fs.copyFileSync(results[0].jsonPath, merged.outJson);
    const failed = results.filter(r => !r.report.summary.pass);
    const passList = ['G0', 'G1', 'G2'];
    const requiredFailed = failed.filter(r => passList.includes(r.report.preset));
    if (merged.failOnError && requiredFailed.length > 0) {
        console.error(`[audit] FAILED required ${requiredFailed.map(r => r.report.preset).join(', ')}`);
        process.exit(1);
    }
    console.log(`[audit] done ${results.length} case(s)`);
}

// ============================================================================
// nonbond subcommand
// ============================================================================

function nonbondUsage() {
    return `Usage: node scripts/nanocrystal.mjs nonbond [options]

Options:
  --mol2 PATH            mol2 file
  --xyz PATH             relaxed xyz file
  --topo PATH            optional topology NPZ (for pre-built collision groups)
  --margin F             collision margin (default: 0.0)
  --collision            use collision mode (getSR_x2_smooth)
  --rcut R               cutoff radius 1
  --rcut2 R              cutoff radius 2
  --out-pairs PATH       write collision pairs JSON`;
}

function cmdNonbond(argv) {
    const cfg = { mol2: null, xyz: null, topo: null, margin: 0.0, outPairs: null, collision: false, rcut: 0, rcut2: 0, config: null };
    for (let i = 0; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => { if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`); return String(argv[++i]); };
        if (a === '--help' || a === '-h') { console.log(nonbondUsage()); process.exit(0); }
        else if (a === '--config') cfg.config = nxt();
        else if (a === '--mol2') cfg.mol2 = nxt();
        else if (a === '--xyz') cfg.xyz = nxt();
        else if (a === '--topo') cfg.topo = nxt();
        else if (a === '--margin') cfg.margin = parseFloat(nxt());
        else if (a === '--out-pairs') cfg.outPairs = nxt();
        else if (a === '--collision') cfg.collision = true;
        else if (a === '--rcut') cfg.rcut = parseFloat(nxt());
        else if (a === '--rcut2') cfg.rcut2 = parseFloat(nxt());
        else throw new Error(`Unknown arg: ${a}\n\n${nonbondUsage()}`);
    }
    let jsonCfg = null;
    if (cfg.config) jsonCfg = loadJSON(cfg.config);
    const merged = mergeConfig(jsonCfg, cfg);
    if (!merged.mol2 || !merged.xyz) throw new Error('Missing --mol2 or --xyz\n\n' + nonbondUsage());
    const mm = loadMMParamsFromDir();
    const mol = loadMolFromMol2(merged.mol2, mm);
    const { pos } = readXyzPositions(fs, merged.xyz);
    applyPositions(mol, pos);
    const nAtoms = mol.atoms.length | 0;
    const posFlat = new Float64Array(nAtoms * 3);
    for (let i = 0; i < nAtoms; i++) { const a = mol.atoms[i]; posFlat[i*3]=+a.pos.x; posFlat[i*3+1]=+a.pos.y; posFlat[i*3+2]=+a.pos.z; }
    let radius = new Float64Array(nAtoms);
    for (let i = 0; i < nAtoms; i++) {
        const at = mm.getAtomTypeForAtom(mol.atoms[i]);
        const r = (at && at.RvdW > 0) ? at.RvdW : 1.5;
        if (!(r > 0)) throw new Error(`invalid RvdW ${r} at atom ${i}`);
        radius[i] = +r;
    }
    let icol, group_atoms, group_nAtoms, bbox_min, bbox_max, excl_icol, groupCap = 32;
    if (merged.topo) {
        const { arrays } = readNpzFile(fs, merged.topo);
        if (!arrays.icol || !arrays.group_atoms) throw new Error(`topology npz missing required arrays: ${merged.topo}`);
        icol = arrays.icol; group_atoms = arrays.group_atoms; group_nAtoms = arrays.group_nAtoms;
        bbox_min = arrays.group_bbox_min; bbox_max = arrays.group_bbox_max; excl_icol = arrays.excl_icol;
        groupCap = arrays.group_cap ? (arrays.group_cap[0] | 0) : 32;
    } else {
        const wg = buildCollisionWorkgroups({ pos: posFlat, mol, radius, groupCap: 32, fillFactor: 0.8 });
        icol = wg.icol; group_atoms = wg.group_atoms; group_nAtoms = wg.group_nAtoms; groupCap = wg.groupCap;
        bbox_min = wg.bbox_min; bbox_max = wg.bbox_max;
        const ex2 = buildExclIcol_1_2_3({ mol, icol, EXCL_MAX: 16, ipbc: 0 });
        excl_icol = ex2.excl;
    }
    const brute = computeNonbondBruteForceKernelStyle({ pos: posFlat, mol, mmParams: mm, surfaceOnly: true, exclude12and13: true, collectPairs: true });
    const grouped = computeNonbondByGroups({ pos: posFlat, mol, mmParams: mm, group_atoms, group_nAtoms, group_bbox_min: bbox_min, group_bbox_max: bbox_max, radius, excl_icol, EXCL_MAX: 16, groupCap, icol, surfaceOnly: true, margin: merged.margin, collectPairs: true, collisionMode: merged.collision, R_cut: merged.rcut, R_cut2: merged.rcut2 });
    const df = maxAbsDiff(brute.f, grouped.f);
    const dE = Math.abs(brute.Etot - grouped.Etot);
    const bruteSet = new Set(brute.pairs.map(([a, b]) => a < b ? `${a},${b}` : `${b},${a}`));
    const groupedSet = new Set(grouped.pairs.map(([a, b]) => a < b ? `${a},${b}` : `${b},${a}`));
    const onlyInBrute = [...bruteSet].filter(k => !groupedSet.has(k));
    const onlyInGrouped = [...groupedSet].filter(k => !bruteSet.has(k));
    console.log(`# nonbond parity`);
    console.log(`mol2: ${path.resolve(merged.mol2)}`);
    console.log(`xyz:  ${path.resolve(merged.xyz)}`);
    if (merged.topo) console.log(`topo: ${path.resolve(merged.topo)}`);
    console.log(`mode: ${merged.collision ? 'collision (getSR_x2_smooth)' : 'full LJ+Q'}`);
    console.log(`margin: ${merged.margin}`);
    if (merged.collision) console.log(`R_cut: ${merged.rcut}, R_cut2: ${merged.rcut2}`);
    console.log(`maxAbsDiff(force): ${df}`);
    console.log(`absDiff(E): ${dE}`);
    console.log(`brute pairs: ${brute.pairs.length}, grouped pairs: ${grouped.pairs.length}`);
    console.log(`only in brute: ${onlyInBrute.length}, only in grouped: ${onlyInGrouped.length}`);
    if (onlyInBrute.length > 0) console.log(`  onlyInBrute (first 10): ${onlyInBrute.slice(0, 10).join(' ')}`);
    if (onlyInGrouped.length > 0) console.log(`  onlyInGrouped (first 10): ${onlyInGrouped.slice(0, 10).join(' ')}`);
    if (merged.outPairs) {
        const icolGroupArr = merged.topo ? readNpzFile(fs, merged.topo).arrays.icolGroup : icol;
        const out = { nAtoms, groupCap, icolGroup: Array.from(icolGroupArr), pairs: grouped.pairs.map(([a, b]) => [a, b]), pairGroups: grouped.pairs.map(([a, b]) => [icolGroupArr[a] | 0, icolGroupArr[b] | 0]) };
        fs.writeFileSync(merged.outPairs, JSON.stringify(out));
        console.log(`wrote collision pairs to ${merged.outPairs}`);
    }
}

// ============================================================================
// main dispatch
// ============================================================================

async function main() {
    const argv = process.argv.slice(2);
    if (argv.length === 0 || argv[0] === '--help' || argv[0] === '-h') {
        console.log(usage());
        process.exit(0);
    }
    const subcmd = argv[0];
    const rest = argv.slice(1);
    switch (subcmd) {
        case 'generate': cmdGenerate(rest); break;
        case 'ensemble': await cmdEnsemble(rest); break;
        case 'topology': await cmdTopology(rest); break;
        case 'rings': cmdRings(rest); break;
        case 'audit': await cmdAudit(rest); break;
        case 'nonbond': cmdNonbond(rest); break;
        default: console.error(`Unknown subcommand '${subcmd}'\n\n${usage()}`); process.exit(1);
    }
}

main().catch((err) => { console.error(err); process.exit(1); });
