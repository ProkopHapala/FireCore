import fs from 'node:fs';
import path from 'node:path';

import { MMParams } from './MMParams.js';
import { EditableMolecule } from './EditableMolecule.js';
import { installMoleculeIOMethods } from './MoleculeIO.js';
import { installMoleculeUtilsMethods } from './MoleculeUtils.js';
import { installMoleculeSelectionMethods } from './MoleculeSelection.js';
import { buildXPDBTopology } from './XPDBTopology.js';
import { buildMMFFLTopology } from './MMFFLTopology.js';

installMoleculeIOMethods(EditableMolecule);
installMoleculeUtilsMethods(EditableMolecule);
installMoleculeSelectionMethods(EditableMolecule);

const __dirname = path.dirname(new URL(import.meta.url).pathname);
const defaultRes = path.resolve(__dirname, '../../cpp/common_resources');
const repoRoot = path.resolve(__dirname, '../..');


/*

Run like this:

node dump_xpdb_topology.mjs --xyz cpp/common_resources/xyz/OCH2.xyz --useMMFFL 1 --angles 1 --addPi 0 --addEpair 0 --reportTypes 0 --mol2Out OUT_js_OCH2.mol2


guanine.xyz
pyridine.xyz
pyrrol.xyz


node web/molgui_webgpu/dump_xpdb_topology.mjs --xyz          ./cpp/common_resources/xyz/guanine.xyz --mol2Out  OUT_js_guanine.mol2
python -u test_TiledJacobi_molecules.py       --molecule ../../cpp/common_resources/xyz/guanine.xyz --mol2_out OUT_py_guanine.mol2 --print_buffers 1 --topology_only 1 



node web/molgui_webgpu/dump_xpdb_topology.mjs --xyz          ./cpp/common_resources/xyz/pyridine.xyz --mol2Out  OUT_js_pyridine.mol2
python -u test_TiledJacobi_molecules.py       --molecule ../../cpp/common_resources/xyz/pyridine.xyz --mol2_out OUT_py_pyridine.mol2 --print_buffers 1 --topology_only 1 


node web/molgui_webgpu/dump_xpdb_topology.mjs --xyz          ./cpp/common_resources/xyz/pyrrol.xyz --mol2Out  OUT_js_pyrrol.mol2
python -u test_TiledJacobi_molecules.py       --molecule ../../cpp/common_resources/xyz/pyrrol.xyz --mol2_out OUT_py_pyrrol.mol2 --print_buffers 1 --topology_only 1 

*/

function parseArgs(argv) {
    const out = {
        xyz: null,
        enableAngles: true,
        useMMFFL: true,
        addPi: true,
        twoPi: true,
        addPiAlign: false,
        addEpair: true,
        addEpairPairs: false,
        reportTypes: true,
        L_pi: 1.0,
        L_epair: 0.5,
        // Match Python defaults from test_TiledJacobi_molecules.py
        k_angle: 100.0,
        k_pi: 50.0,
        k_pi_orth: 30.0,
        k_pi_align: 15.0,
        k_ep: 40.0,
        k_ep_orth: 25.0,
        k_ep_pair: 10.0,
        maxBonds: 16,
        defaultL: 1.3,
        defaultK: 200.0,
        mol2Out: null,
        elementTypes: path.join(defaultRes, 'ElementTypes.dat'),
        atomTypes: path.join(defaultRes, 'AtomTypes.dat'),
        bondTypes: path.join(defaultRes, 'BondTypes.dat'),
        angleTypes: path.join(defaultRes, 'AngleTypes.dat'),
    };
    for (let i = 2; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => {
            if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`);
            i++;
            return String(argv[i]);
        };
        if (a === '--xyz') out.xyz = nxt();
        else if (a === '--angles') out.enableAngles = (nxt() !== '0');
        else if (a === '--useMMFFL') out.useMMFFL = (nxt() !== '0');
        else if (a === '--addPi') out.addPi = (nxt() !== '0');
        else if (a === '--twoPi') out.twoPi = (nxt() !== '0');
        else if (a === '--addPiAlign') out.addPiAlign = (nxt() !== '0');
        else if (a === '--addEpair') out.addEpair = (nxt() !== '0');
        else if (a === '--addEpairPairs') out.addEpairPairs = (nxt() !== '0');
        else if (a === '--reportTypes') out.reportTypes = (nxt() !== '0');
        else if (a === '--L_pi') out.L_pi = +nxt();
        else if (a === '--L_epair') out.L_epair = +nxt();
        else if (a === '--k_angle') out.k_angle = +nxt();
        else if (a === '--k_pi') out.k_pi = +nxt();
        else if (a === '--k_pi_orth') out.k_pi_orth = +nxt();
        else if (a === '--k_pi_align') out.k_pi_align = +nxt();
        else if (a === '--k_ep') out.k_ep = +nxt();
        else if (a === '--k_ep_orth') out.k_ep_orth = +nxt();
        else if (a === '--k_ep_pair') out.k_ep_pair = +nxt();
        else if (a === '--mol2Out') out.mol2Out = nxt();
        else if (a === '--maxBonds') out.maxBonds = parseInt(nxt(), 10) | 0;
        else if (a === '--defaultL') out.defaultL = +nxt();
        else if (a === '--defaultK') out.defaultK = +nxt();
        else if (a === '--elementTypes') out.elementTypes = nxt();
        else if (a === '--atomTypes') out.atomTypes = nxt();
        else if (a === '--bondTypes') out.bondTypes = nxt();
        else if (a === '--angleTypes') out.angleTypes = nxt();
        else throw new Error(`Unknown arg ${a}`);
    }
    if (!out.xyz) throw new Error('Usage: node dump_xpdb_topology.mjs --xyz path/to/file.xyz [--useMMFFL 0|1] [--addPi 0|1] [--addEpair 0|1] [--mol2Out out.mol2]');
    return out;
}

function packBondArrays(bondsAdj, nAtoms, nMaxBonded) {
    const bondIndices = new Int32Array(nAtoms * nMaxBonded);
    bondIndices.fill(-1);
    const bondLenStiff = new Float32Array(nAtoms * nMaxBonded * 2);
    bondLenStiff.fill(0);

    const addBond = (i, j, L, K) => {
        let slot = -1;
        const base = i * nMaxBonded;
        for (let k = 0; k < nMaxBonded; k++) {
            const jj = bondIndices[base + k];
            if (jj === -1 || jj === j) { slot = k; break; }
        }
        if (slot === -1) throw new Error(`packBondArrays: atom ${i} max bonds exceeded (nMaxBonded=${nMaxBonded})`);
        bondIndices[base + slot] = j;
        bondLenStiff[(base + slot) * 2 + 0] = L;
        bondLenStiff[(base + slot) * 2 + 1] = K;
    };

    for (let i = 0; i < nAtoms; i++) {
        const neighs = bondsAdj[i] || [];
        for (let b of neighs) {
            addBond(i, b[0], +b[1], +b[2]);
        }
    }

    return { bondIndices, bondLenStiff };
}

function formatArray2D_Int(arr, nRows, nCols) {
    const lines = [];
    for (let i = 0; i < nRows; i++) {
        const row = [];
        const base = i * nCols;
        for (let j = 0; j < nCols; j++) row.push(String(arr[base + j] | 0));
        lines.push('[' + row.join(' ') + ']');
    }
    return lines.join('\n');
}

function formatArray2D_FloatPairs(arr, nRows, nCols) {
    const lines = [];
    for (let i = 0; i < nRows; i++) {
        const row = [];
        const base = i * nCols;
        for (let j = 0; j < nCols; j++) {
            const L = arr[(base + j) * 2 + 0];
            const K = arr[(base + j) * 2 + 1];
            row.push(`${L.toFixed(6)},${K.toFixed(6)}`);
        }
        lines.push('[' + row.join(' ') + ']');
    }
    return lines.join('\n');
}

function mol2FromTopology(topo, name = 'Molecule') {
    const fmt9_4 = (x) => (+x).toFixed(4).padStart(9);
    const fmt10_4 = (x) => (+x).toFixed(4).padStart(10);

    const nAtoms = topo.n_all | 0;
    const bondsSet = new Set();
    const bondsAll = [];
    const addPair = (a, b) => {
        let i = a | 0;
        let j = b | 0;
        if (j < i) { const t = i; i = j; j = t; }
        const key = `${i}|${j}`;
        if (bondsSet.has(key)) return;
        bondsSet.add(key);
        bondsAll.push([i, j]);
    };

    // Python parity: bonds_all = sorted(set(bonds_primary + bonds_derived))
    for (const [a, b] of topo.bonds_primary) addPair(a, b);
    if (topo.bonds_angle) for (const [a, b] of topo.bonds_angle) addPair(a, b);
    if (topo.bonds_pi) for (const [a, b] of topo.bonds_pi) addPair(a, b);
    if (topo.bonds_pi_align) for (const [a, b] of topo.bonds_pi_align) addPair(a, b);
    if (topo.bonds_epair) for (const [a, b] of topo.bonds_epair) addPair(a, b);
    if (topo.bonds_epair_pair) for (const [a, b] of topo.bonds_epair_pair) addPair(a, b);

    bondsAll.sort((p, q) => (p[0] - q[0]) || (p[1] - q[1]));

    const out = [];
    out.push('@<TRIPOS>MOLECULE');
    out.push(name);
    out.push(`${String(nAtoms).padStart(3)} ${String(bondsAll.length).padStart(3)} 0 0 0`);
    out.push('SMALL');
    out.push('GASTEIGER');
    out.push('');
    out.push('@<TRIPOS>ATOM');

    for (let i = 0; i < nAtoms; i++) {
        const tname = String(topo.type_names[i]);
        let ename = tname.replace(/\./g, '_');
        ename = ename.split('_')[0];
        const p = topo.apos[i];
        const x = p[0], y = p[1], z = p[2];
        const q = 0.0;
        out.push(`${String(i + 1).padStart(7)} ${ename.padEnd(8)} ${fmt9_4(x)} ${fmt9_4(y)} ${fmt9_4(z)} ${tname.padEnd(5)} ${String(1).padStart(3)}  ${'UNL1'.padEnd(7)} ${fmt10_4(q)}`);
    }

    out.push('@<TRIPOS>BOND');
    for (let i = 0; i < bondsAll.length; i++) {
        const a1 = (bondsAll[i][0] | 0) + 1;
        const a2 = (bondsAll[i][1] | 0) + 1;
        out.push(`${String(i + 1).padStart(6)} ${String(a1).padStart(5)} ${String(a2).padStart(5)} ${String(1).padStart(4)}`);
    }
    return out.join('\n') + '\n';
}

async function main() {
    const args = parseArgs(process.argv);

    const mm = new MMParams();
    const readText = (p) => fs.readFileSync(p, 'utf8');

    mm.parseElementTypes(readText(args.elementTypes));
    mm.parseAtomTypes(readText(args.atomTypes));
    mm.parseBondTypes(readText(args.bondTypes));
    mm.parseAngleTypes(readText(args.angleTypes));

    const resolveXYZPath = (pIn) => {
        const p = String(pIn || '').trim();
        if (!p) throw new Error('Empty --xyz path');
        const tried = [];
        const tryPath = (q) => {
            const qq = path.resolve(q);
            tried.push(qq);
            if (fs.existsSync(qq)) return qq;
            return null;
        };

        // 1) as provided (relative to cwd or absolute)
        let ok = tryPath(p);
        if (ok) return ok;

        // 2) treat '/cpp/..' as repo-root relative
        if (p.startsWith('/cpp/')) {
            ok = tryPath(path.join(repoRoot, p.slice(1))); // remove leading '/'
            if (ok) return ok;
        }
        // 3) common case: 'cpp/..' already relative, but passed from different cwd
        if (p.startsWith('cpp/')) {
            ok = tryPath(path.join(repoRoot, p));
            if (ok) return ok;
        }

        throw new Error(`--xyz file not found. Provided='${p}'. Tried:\n${tried.map(s => '  ' + s).join('\n')}`);
    };

    const xyzPath = resolveXYZPath(args.xyz);
    const xyzText = fs.readFileSync(xyzPath, 'utf8');
    const parsed = EditableMolecule.parseXYZ(xyzText);

    const mol = new EditableMolecule();
    mol.clear();
    mol.appendParsedSystem(parsed);
    mol.recalculateBonds(mm);

    let bondsAdj = null;
    let stats = null;
    let topo = null;

    if (args.useMMFFL) {
        topo = buildMMFFLTopology(mol, mm, {
            report_types: !!args.reportTypes,
            add_angle: !!args.enableAngles,
            add_pi: !!args.addPi,
            two_pi: !!args.twoPi,
            add_pi_align: !!args.addPiAlign,
            add_epair: !!args.addEpair,
            add_epair_pairs: !!args.addEpairPairs,
            L_pi: args.L_pi,
            L_epair: args.L_epair,
            k_angle: args.k_angle,
            k_pi: args.k_pi,
            k_pi_orth: args.k_pi_orth,
            k_pi_align: args.k_pi_align,
            k_ep: args.k_ep,
            k_ep_orth: args.k_ep_orth,
            k_ep_pair: args.k_ep_pair,
            defaultL: args.defaultL,
            defaultK: args.defaultK,
        });
        // attach typeName to atoms for MOL2 export parity
        for (let i = 0; i < mol.atoms.length; i++) mol.atoms[i].typeName = topo.type_names[i];

        // Build XPDB adjacency from linear bonds only (primary+linear). Keep simple for parity.
        const nAtoms2 = mol.atoms.length;
        bondsAdj = new Array(nAtoms2);
        for (let i = 0; i < nAtoms2; i++) bondsAdj[i] = [];
        const maxB = args.maxBonds | 0;
        const addEdge = (i, j, L, K) => {
            const arr = bondsAdj[i];
            if (arr.length >= maxB) return;
            arr.push([j, +L, +K]);
        };
        // primary
        for (const [a, b] of topo.bonds_primary) {
            const ta = topo.type_names[a];
            const tb = topo.type_names[b];
            const bp = mm.getBondParams(ta, tb);
            const l0 = (bp && bp.l0 > 0) ? bp.l0 : args.defaultL;
            const K = (bp && bp.k > 0) ? bp.k : args.defaultK;
            addEdge(a, b, l0, K);
            addEdge(b, a, l0, K);
        }
        // linear
        for (const lb of topo.bonds_linear) {
            const a = lb[0] | 0;
            const b = lb[1] | 0;
            const l0 = +lb[2];
            const K = +lb[3];
            addEdge(a, b, l0, K);
            addEdge(b, a, l0, K);
        }

        stats = {
            atoms: topo.n_all,
            n_real: topo.n_real,
            n_all: topo.n_all,
            bonds_primary: topo.bonds_primary.length,
            bonds_linear: topo.bonds_linear.length,
        };
    } else {
        const r = buildXPDBTopology(mol, mm, {
            includeAngleConstraints: !!args.enableAngles,
            maxBonds: args.maxBonds,
            defaultL: args.defaultL,
            defaultK: args.defaultK,
        });
        bondsAdj = r.bondsAdj;
        stats = r.stats;
    }

    const nAtoms = mol.nAtoms || mol.atoms.length;
    const nMaxBonded = 16;
    const { bondIndices, bondLenStiff } = packBondArrays(bondsAdj, nAtoms, nMaxBonded);

    const outLines = [];
    outLines.push(`[DEBUG] dump_xpdb_topology.mjs xyz=${args.xyz}`);
    outLines.push(`[DEBUG] stats=${JSON.stringify(stats)}`);
    outLines.push(`[DEBUG] nAtoms=${nAtoms} nMaxBonded=${nMaxBonded} enableAngles=${args.enableAngles}`);
    outLines.push(`[DEBUG] useMMFFL=${args.useMMFFL} addPi=${args.addPi} addEpair=${args.addEpair} twoPi=${args.twoPi} addPiAlign=${args.addPiAlign} addEpairPairs=${args.addEpairPairs}`);
    if (topo) {
        outLines.push(`[TOPO] n_real=${topo.n_real} n_all=${topo.n_all}`);
        outLines.push(`[TOPO] bonds_primary=${topo.bonds_primary.length} bonds_angle=${topo.bonds_angle.length} bonds_pi=${topo.bonds_pi.length} bonds_pi_align=${topo.bonds_pi_align.length} bonds_epair=${topo.bonds_epair.length} bonds_epair_pair=${topo.bonds_epair_pair.length}`);
    }
    outLines.push('');
    outLines.push('[DEBUG] pos0:');
    for (let i = 0; i < nAtoms; i++) {
        const a = mol.atoms[i];
        outLines.push(`${a.pos.x.toFixed(6)} ${a.pos.y.toFixed(6)} ${a.pos.z.toFixed(6)}`);
    }
    outLines.push('');
    outLines.push('[DEBUG] bond_indices:');
    outLines.push(formatArray2D_Int(bondIndices, nAtoms, nMaxBonded));
    outLines.push('');
    outLines.push('[DEBUG] bondLenStiff (L,K pairs):');
    outLines.push(formatArray2D_FloatPairs(bondLenStiff, nAtoms, nMaxBonded));

    const outPath = path.resolve(process.cwd(), 'OUT-dump_xpdb_topology.txt');
    fs.writeFileSync(outPath, outLines.join('\n') + '\n');
    console.log(`Wrote ${outPath}`);

    if (args.mol2Out) {
        const mol2 = topo ? mol2FromTopology(topo, 'Molecule') : mol.toMol2String({ name: path.basename(args.mol2Out).replace(/\.mol2$/i, '') });
        fs.writeFileSync(args.mol2Out, mol2);
        console.log(`Wrote ${args.mol2Out}`);
    }
}

main();
