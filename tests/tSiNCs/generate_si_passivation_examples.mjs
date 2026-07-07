#!/usr/bin/env node
/// Generate ~1 nm Si nanocrystal NPZ bundles with varied surface passivation (viewer / FTIR examples).
import fs from 'node:fs';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(__dirname, '../..');
const OUT_ROOT = path.join(__dirname, 'fixtures/si_1nm_passivation');
const CIF = 'cpp/common_resources/crystals/Si-sym.cif';
const CLI = path.join(REPO, 'scripts/export_nanocrystal_bundle.mjs');

/// ~1–1.3 nm Si NCs: sphere diameter 11–13 Å, or small {111}/{100} plane cuts.
const CASES = [
    { dir: '01_sphere_11A_H_standard', id: 'si_sphere_11A_H', desc: 'Spherical ~1.1 nm (R=5.5 Å), H passivation, no defects', args: ['--sphere', '5.5', '--insertProb', '0', '--collapseProb', '0', '--seed', '42'] },
    { dir: '02_sphere_13A_H_standard', id: 'si_sphere_13A_H', desc: 'Spherical ~1.3 nm (R=6.5 Å), H passivation, no defects', args: ['--sphere', '6.5', '--insertProb', '0', '--collapseProb', '0', '--seed', '42'] },
    { dir: '03_sphere_13A_H_bridge_insert', id: 'si_sphere_13A_insert', desc: 'R=6.5 Å, H caps + stochastic -SiH2- bridge insert (insertProb=0.18)', args: ['--sphere', '6.5', '--insertProb', '0.18', '--collapseProb', '0', '--seed', '1'] },
    { dir: '04_sphere_13A_H_bridge_collapse', id: 'si_sphere_13A_collapse', desc: 'R=6.5 Å, H caps + bridge collapse → fused surface rings (collapseProb=0.18)', args: ['--sphere', '6.5', '--insertProb', '0', '--collapseProb', '0.18', '--seed', '1'] },
    { dir: '05_sphere_13A_H_mixed_defects', id: 'si_sphere_13A_mixed', desc: 'R=6.5 Å, H caps + insert + collapse (mixed surface defects)', args: ['--sphere', '6.5', '--insertProb', '0.12', '--collapseProb', '0.10', '--seed', '1'] },
    { dir: '06_sphere_13A_bare_surface', id: 'si_sphere_13A_bare', desc: 'R=6.5 Å, no H caps — bare / under-coordinated surface Si', args: ['--sphere', '6.5', '--caps', 'none', '--insertProb', '0', '--collapseProb', '0', '--seed', '42'] },
    { dir: '07_faceted_111_H', id: 'si_faceted_111_H', desc: '{111}-faceted ~1 nm chunk (G1 plane cut), H passivation', args: ['--planes', 'a111', '--nx', '1', '--planeCScale', '0.42', '--planeCJitter', '0', '--insertProb', '0', '--collapseProb', '0', '--seed', '42'] },
    { dir: '09_sphere_13A_H_fuse', id: 'si_sphere_13A_fuse', desc: 'R=6.5 Å, H caps + SiH2 clash fusion (cleave 2 H, bond Si–Si)', args: ['--sphere', '6.5', '--resolveClashes', '0', '--fuseProb', '1', '--insertProb', '0', '--collapseProb', '0', '--seed', '42'] },
];

function runCase(c) {
    const outDir = path.join(OUT_ROOT, c.dir);
    fs.mkdirSync(outDir, { recursive: true });
    const argv = [CLI, '--cif', CIF, '--out', outDir, '--id', c.id, ...c.args];
    const r = spawnSync('node', argv, { cwd: REPO, encoding: 'utf8' });
    if (r.status !== 0) {
        console.error(r.stdout);
        console.error(r.stderr);
        throw new Error(`generate_si_passivation_examples: failed ${c.dir} (exit ${r.status})`);
    }
    const metaPath = path.join(outDir, 'meta.json');
    const meta = JSON.parse(fs.readFileSync(metaPath, 'utf8'));
    return { ...c, outDir, meta, natoms: meta.natoms, defects: meta.defects };
}

function main() {
    fs.mkdirSync(OUT_ROOT, { recursive: true });
    const index = [];
    for (const c of CASES) {
        const row = runCase(c);
        console.log(`[ok] ${c.dir}  natoms=${row.natoms}  defects=${JSON.stringify(row.defects)}`);
        index.push({
            dir: c.dir,
            id: c.id,
            description: c.desc,
            natoms: row.natoms,
            nbonds: row.meta.nbonds,
            defects: row.defects,
            gen_params: row.meta.gen_params,
            files: ['01_crystal.npz', '03_topology.npz', 'meta.json'],
        });
    }
  const indexPath = path.join(OUT_ROOT, 'index.json');
    fs.writeFileSync(indexPath, JSON.stringify({ schema: 'si_1nm_passivation_v1', cases: index }, null, 2));
    console.log(`\n[generate_si_passivation_examples] wrote ${CASES.length} cases under ${OUT_ROOT}`);
    console.log(`  index: ${indexPath}`);
}

main();
