#!/usr/bin/env node
/// Generate small symmetric H-passivated Si and C diamond nanocrystals (sphere cuts, <100 atoms).
/// Config: tests/tSiNCs/small_symmetric_nc.json → OUT_small_nc/<material>/<id>/{init.mol2,init.xyz,gen_params.json}
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { toMol2String, toXYZString } from '../../web/molgui_webgpu/MoleculeIO.js';
import { loadMMParams, buildCrystalFromCIFText, getHeavyZ, generateNanocrystal, genArgsFromConfig, mulberry32 } from '../../web/molgui_webgpu/Nanocrystals.js';
import { molToCrystalArrays, crystalToXYZ } from '../../web/common_js/npzIO.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const TEST_DIR = __dirname;
const REPO = path.resolve(__dirname, '../..');
const RES = path.join(REPO, 'cpp/common_resources');

function loadCfg(p) {
    const cfgPath = p ? path.resolve(p) : path.join(TEST_DIR, 'small_symmetric_nc.json');
    return { cfg: JSON.parse(fs.readFileSync(cfgPath, 'utf8')), cfgPath };
}

function main() {
    const cfgArg = process.argv[2] && !process.argv[2].startsWith('-') ? process.argv[2] : null;
    const { cfg, cfgPath } = loadCfg(cfgArg);
    const outRoot = path.resolve(TEST_DIR, cfg.outDir || 'OUT_small_nc');
    const seed = cfg.seed | 0;
    const rows = [];
    console.log(`[small_nc] config ${cfgPath} → ${outRoot}`);
    for (const mat of cfg.materials) {
        for (const cut of mat.cuts) {
            const id = `${mat.id}_${cut.id}`;
            const edir = path.join(outRoot, mat.id, cut.id);
            fs.mkdirSync(edir, { recursive: true });
            const genParams = {
                cif: mat.cif,
                applySymmetry: mat.applySymmetry ?? 0,
                cut: { cutMode: 'sphere', sphereR: cut.sphereR, sphereNrep: cut.sphereNrep ?? 4 },
                replication: cfg.replication || { nx: 2, ny: 2, nz: 2 },
                passivation: cfg.passivation,
                defects: cfg.defects || { insertProb: 0, collapseProb: 0 },
            };
            const args = genArgsFromConfig(genParams, RES, { seed, samples: 1, maxFiles: 1 });
            if (seed !== 0) { const rnd = mulberry32(seed); Math.random = rnd; }
            const cell = buildCrystalFromCIFText(fs.readFileSync(args.cif, 'utf8'), args);
            const mm = loadMMParams(args);
            const { mol } = generateNanocrystal(args, cell, mm, getHeavyZ(cell), 1);
            const arrays = molToCrystalArrays(mol);
            fs.writeFileSync(path.join(edir, 'init.mol2'), toMol2String(mol, { name: id }));
            fs.writeFileSync(path.join(edir, 'init.xyz'), crystalToXYZ(arrays.pos, arrays.Z));
            fs.writeFileSync(path.join(edir, 'gen_params.json'), JSON.stringify({ id, label: `${mat.label} R=${cut.sphereR}`, ...genParams, natoms: arrays.natoms }, null, 2));
            rows.push({ id, mat: mat.id, sphereR: cut.sphereR, natoms: arrays.natoms, path: edir });
            console.log(`  ${id}: ${arrays.natoms} atoms (R=${cut.sphereR} Å)`);
        }
    }
    const md = ['# Small symmetric H-passivated nanocrystals (<100 atoms)', '', `Config: \`${path.relative(REPO, cfgPath)}\``, '', '| id | material | R (Å) | atoms |', '|---|---|---:|---:|'];
    for (const r of rows) md.push(`| \`${r.id}\` | ${r.mat} | ${r.sphereR} | ${r.natoms} |`);
    fs.writeFileSync(path.join(outRoot, 'index.md'), md.join('\n') + '\n');
    const over = rows.filter(r => r.natoms >= 100);
    if (over.length) console.warn(`[small_nc] WARNING: ${over.length} structure(s) have natoms >= 100: ${over.map(r => r.id).join(', ')}`);
    console.log(`[small_nc] done: ${rows.length} structures, index ${path.join(outRoot, 'index.md')}`);
}

main();
