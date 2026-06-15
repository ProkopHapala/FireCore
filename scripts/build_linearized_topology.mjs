#!/usr/bin/env node
/// Headless MMFFL linearized topology builder (K12/K13/K14). Job 2: Linearized_topology.progress.md
import fs from 'node:fs';
import path from 'node:path';

import { MMParams } from '../web/molgui_webgpu/MMParams.js';
import { EditableMolecule } from '../web/molgui_webgpu/EditableMolecule.js';
import { installMoleculeIOMethods, Z_TO_SYMBOL } from '../web/molgui_webgpu/MoleculeIO.js';
import { installMoleculeUtilsMethods } from '../web/molgui_webgpu/MoleculeUtils.js';
import {
    buildMMFFLTopology,
    enumerateGraphDistances,
    exportTopologySpringList,
    packLinearTopologyForGPU,
} from '../web/molgui_webgpu/MMFFLTopology.js';
import { writeTopologyNpzFile } from '../web/molgui_webgpu/LinearizedTopologyNpz.js';
import { exportStickViewerHTML, exportStickViewerNpzLoaderHTML } from '../web/molgui_webgpu/LinearizedTopologyViewer.js';

installMoleculeIOMethods(EditableMolecule);
installMoleculeUtilsMethods(EditableMolecule);

const __dirname = path.dirname(new URL(import.meta.url).pathname);
const repoRoot = path.resolve(__dirname, '..');
const defaultRes = path.join(repoRoot, 'cpp/common_resources');

function parseArgs(argv) {
    const out = {
        input: null,
        outDir: path.join(repoRoot, 'tests/tSiNCs/linearized'),
        addAngle: true,
        addDihedral: true,
        htmlViewer: true,
        exportNpz: true,
        exportJson: false,
        npzViewer: true,
        K12: 0,
        K13: 0,
        K14: 0,
        maxNeighbors: 48,
        graphDist: true,
        elementTypes: path.join(defaultRes, 'ElementTypes.dat'),
        atomTypes: path.join(defaultRes, 'AtomTypes.dat'),
        bondTypes: path.join(defaultRes, 'BondTypes.dat'),
        angleTypes: path.join(defaultRes, 'AngleTypes.dat'),
        dihedralTypes: path.join(defaultRes, 'DihedralTypes.dat'),
    };
    for (let i = 2; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => { if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`); i++; return String(argv[i]); };
        if (a === '--input') out.input = nxt();
        else if (a === '--out-dir') out.outDir = nxt();
        else if (a === '--add-angle') out.addAngle = (nxt() !== '0');
        else if (a === '--add-dihedral') out.addDihedral = (nxt() !== '0');
        else if (a === '--html-viewer') out.htmlViewer = (nxt() !== '0');
        else if (a === '--export-npz') out.exportNpz = (nxt() !== '0');
        else if (a === '--export-json') out.exportJson = (nxt() !== '0');
        else if (a === '--npz-viewer') out.npzViewer = (nxt() !== '0');
        else if (a === '--K12') out.K12 = +nxt();
        else if (a === '--K13') out.K13 = +nxt();
        else if (a === '--K14') out.K14 = +nxt();
        else if (a === '--max-neighbors') out.maxNeighbors = parseInt(nxt(), 10) | 0;
        else if (a === '--graph-dist') out.graphDist = (nxt() !== '0');
        else if (a === '--elementTypes') out.elementTypes = nxt();
        else if (a === '--atomTypes') out.atomTypes = nxt();
        else if (a === '--bondTypes') out.bondTypes = nxt();
        else if (a === '--angleTypes') out.angleTypes = nxt();
        else if (a === '--dihedralTypes') out.dihedralTypes = nxt();
        else throw new Error(`Unknown arg ${a}`);
    }
    if (!out.input) throw new Error('Usage: node scripts/build_linearized_topology.mjs --input path.mol2|xyz [--out-dir dir] [--export-json 0|1] [--export-npz 0|1]');
    return out;
}

function resolveInputPath(pIn) {
    const p = String(pIn).trim();
    const tried = [];
    const tryPath = (q) => { const qq = path.resolve(q); tried.push(qq); return fs.existsSync(qq) ? qq : null; };
    let ok = tryPath(p);
    if (ok) return ok;
    ok = tryPath(path.join(repoRoot, p));
    if (ok) return ok;
    throw new Error(`--input not found: '${p}'. Tried:\n${tried.map(s => '  ' + s).join('\n')}`);
}

function loadMolecule(inputPath, mm) {
    const ext = path.extname(inputPath).toLowerCase();
    const text = fs.readFileSync(inputPath, 'utf8');
    const mol = new EditableMolecule();
    mol.clear();
    if (ext === '.mol2') {
        mol.appendParsedSystem(EditableMolecule.parseMol2(text));
    } else if (ext === '.xyz') {
        mol.appendParsedSystem(EditableMolecule.parseXYZ(text));
        mol.recalculateBonds(mm);
    } else {
        throw new Error(`Unsupported input extension '${ext}' (use .mol2 or .xyz)`);
    }
    if (!mol.bonds || mol.bonds.length === 0) mol.recalculateBonds(mm);
    return mol;
}

function buildPrimaryAdj(mol, mm, typeNames) {
    const n = mol.atoms.length;
    const adj = new Array(n);
    for (let i = 0; i < n; i++) adj[i] = [];
    for (let i = 0; i < mol.bonds.length; i++) {
        const b = mol.bonds[i];
        b.ensureIndices(mol);
        const a = b.a | 0, c = b.b | 0;
        const bp = mm.getBondParams(typeNames[a], typeNames[c]);
        const l0 = +bp.l0, K = +bp.k;
        adj[a].push([c, l0, K]);
        adj[c].push([a, l0, K]);
    }
    return adj;
}

function validateSpringGraphDistances(topo, bondsAdj1) {
    const gdist = enumerateGraphDistances(bondsAdj1, 6);
    const pairDist = (a, b) => {
        const n = bondsAdj1.length;
        const dist = new Int32Array(n);
        dist.fill(-1);
        const q = [a | 0];
        dist[a] = 0;
        let qi = 0;
        while (qi < q.length) {
            const u = q[qi++];
            const du = dist[u] | 0;
            for (const nb of bondsAdj1[u] || []) {
                const v = nb[0] | 0;
                if (dist[v] >= 0) continue;
                dist[v] = du + 1;
                q.push(v);
            }
        }
        return dist[b | 0] | 0;
    };
    const errors = [];
    for (let k = 0; k < topo.bonds_primary.length; k++) {
        const e = topo.bonds_primary[k];
        const d = pairDist(e[0], e[1]);
        if (d !== 1) errors.push(`bond (${e[0]},${e[1]}) graph_dist=${d} expected 1`);
    }
    for (const e of topo.bonds_angle) {
        const d = pairDist(e[0], e[1]);
        if (d !== 2) errors.push(`angle spring (${e[0]},${e[1]}) graph_dist=${d} expected 2`);
    }
    for (const e of topo.bonds_dihedral) {
        const d = pairDist(e[0], e[1]);
        if (d !== 3) errors.push(`dihedral spring (${e[0]},${e[1]}) graph_dist=${d} expected 3 (MMFF 1-4)`);
    }
    return { gdist, errors };
}

function stemName(inputPath) {
    return path.basename(inputPath, path.extname(inputPath));
}

async function main() {
    const args = parseArgs(process.argv);
    const inputPath = resolveInputPath(args.input);
    const outDir = path.resolve(args.outDir);
    fs.mkdirSync(outDir, { recursive: true });

    const readText = (p) => fs.readFileSync(p, 'utf8');
    const mm = new MMParams();
    mm.parseElementTypes(readText(args.elementTypes));
    mm.parseAtomTypes(readText(args.atomTypes));
    mm.parseBondTypes(readText(args.bondTypes));
    mm.parseAngleTypes(readText(args.angleTypes));
    mm.parseDihedralTypes(readText(args.dihedralTypes));

    const mol = loadMolecule(inputPath, mm);
    const topoOpts = {
        type_source: 'table',
        add_angle: args.addAngle,
        add_dihedral: args.addDihedral,
        add_pi: false,
        add_epair: false,
        K12: args.K12,
        K13: args.K13,
        K14: args.K14,
    };
    const topo = buildMMFFLTopology(mol, mm, topoOpts);

    const nBond = topo.bonds_primary.length | 0;
    const nAng = topo.bonds_angle.length | 0;
    const nDih = (topo.bonds_dihedral ? topo.bonds_dihedral.length : 0) | 0;
    console.log(`[build_linearized_topology] ${path.basename(inputPath)}: atoms=${topo.n_real} bonds=${nBond} angle_springs=${nAng} dihedral_springs=${nDih}`);

    const typeNames = [];
    for (let i = 0; i < topo.n_real; i++) typeNames[i] = mm.resolveTypeNameTable(Z_TO_SYMBOL[mol.atoms[i].Z] || 'X');
    const bondsAdj1 = buildPrimaryAdj(mol, mm, typeNames);

    let graphReport = null;
    if (args.graphDist) {
        const { gdist, errors } = validateSpringGraphDistances(topo, bondsAdj1);
        graphReport = { hist_undirected_pairs: gdist.hist, per_atom_max_dist: Array.from(gdist.perAtomMax), validation_errors: errors };
        console.log('[graph_dist] undirected pair histogram (dist 1..4):', gdist.hist.slice(1, 5).join(', '));
        if (errors.length > 0) {
            console.error('[graph_dist] VALIDATION ERRORS:');
            for (const e of errors.slice(0, 20)) console.error('  ', e);
            if (errors.length > 20) console.error(`  ... and ${errors.length - 20} more`);
            throw new Error(`graph distance validation failed (${errors.length} errors)`);
        }
    }

    const packing = packLinearTopologyForGPU(topo, { maxNeighbors: args.maxNeighbors });
    console.log(`[packing] maxDegree=${packing.maxDegree} maxNeighbors=${packing.maxNeighbors}`);

    const base = stemName(inputPath);
    const meta = {
        source: inputPath,
        n_bond: nBond,
        n_angle: nAng,
        n_dihedral: nDih,
    };

    if (args.exportJson) {
        const springs = exportTopologySpringList(topo);
        const topoJson = {
            meta: {
                ...meta,
                natoms: topo.n_real | 0,
                K12_override: args.K12,
                K13_override: args.K13,
                K14_override: args.K14,
                max_neighbors_required: packing.maxDegree,
                max_neighbors_budget: args.maxNeighbors,
                graph_dist: graphReport,
            },
            springs,
        };
        const jsonPath = path.join(outDir, `${base}_topology.json`);
        fs.writeFileSync(jsonPath, JSON.stringify(topoJson, null, 2));
        console.log(`Wrote ${jsonPath}`);
    }

    if (args.htmlViewer) {
        const html = exportStickViewerHTML(topo, `${base} MMFFL`);
        const htmlPath = path.join(outDir, `${base}_debug_viewer.html`);
        fs.writeFileSync(htmlPath, html);
        console.log(`Wrote ${htmlPath}`);
    }

    if (args.exportNpz) {
        const npzPath = path.join(outDir, `${base}_topology.npz`);
        writeTopologyNpzFile(fs, npzPath, topo, packing, mol, meta);
        console.log(`Wrote ${npzPath}`);
    }

    if (args.npzViewer) {
        const npzHtml = exportStickViewerNpzLoaderHTML(`${base} MMFFL`, base);
        const npzHtmlPath = path.join(outDir, `${base}_npz_viewer.html`);
        fs.writeFileSync(npzHtmlPath, npzHtml);
        console.log(`Wrote ${npzHtmlPath}`);
    }
}

main().catch((err) => { console.error(err); process.exit(1); });
