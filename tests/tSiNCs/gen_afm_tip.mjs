#!/usr/bin/env node
/// Generate H-passivated Si or Diamond nanocrystal with a spike (AFM-tip-like)
/// apex pointing along the [111] direction, then rotate so [111] → z-axis.
/// The shape is a truncated tetrahedron: 3 side facets from {111} planes +
/// 1 bottom truncation plane. The apex naturally terminates at a single atom.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(__dirname, '../..');
const RES = path.join(REPO, 'cpp/common_resources');

import { Vec3 } from '../../web/common_js/Vec3.js';
import { Mat3 } from '../../web/common_js/Mat3.js';
import { toMol2String, toXYZString } from '../../web/molgui_webgpu/MoleculeIO.js';
import {
    loadMMParams, buildCrystalFromCIFText, getHeavyZ,
    pruneUndercoordinatedHeavyIter, assignCrystalAtomTypes,
} from '../../web/molgui_webgpu/Nanocrystals.js';
import * as CrystalUtils from '../../web/molgui_webgpu/CrystalUtils.js';

function parseArgs() {
    const argv = process.argv.slice(2);
    const opts = {
        material: 'Si',           // 'Si' or 'C' (diamond)
        sideHalf: 5.0,            // Å — controls apex height (apex at 3*sideHalf along [111])
        bottomDepth: 6.0,         // Å — controls base position below origin
        nRep: 5,                  // cell replication (centered)
        outDir: path.join(__dirname, 'OUT_afm_tip'),
        prefix: null,             // auto from material
    };
    for (let i = 0; i < argv.length; i++) {
        const a = argv[i];
        if (a === '--material' || a === '-m') opts.material = argv[++i];
        else if (a === '--sideHalf') opts.sideHalf = parseFloat(argv[++i]);
        else if (a === '--bottomDepth') opts.bottomDepth = parseFloat(argv[++i]);
        else if (a === '--nRep') opts.nRep = parseInt(argv[++i], 10) | 0;
        else if (a === '--outDir') opts.outDir = argv[++i];
        else if (a === '--help' || a === '-h') {
            console.log('Usage: gen_afm_tip.mjs [--material Si|C] [--sideHalf 5.0] [--bottomDepth 6.0] [--nRep 5] [--outDir DIR]');
            process.exit(0);
        }
    }
    return opts;
}

/// Build custom planes for a tetrahedral tip along [111].
/// Side facets: 3 planes with normals along [1,1,-1], [1,-1,1], [-1,1,1]
///   (the 3 {111} directions that are NOT [111]).
/// Each side plane keeps atoms where n·r <= sideHalf (interior of tetrahedron).
/// Bottom plane: normal [-1,-1,-1], keeps atoms where n·r <= bottomDepth.
/// Apex forms naturally at distance 3*sideHalf along [111].
function buildTipPlanes(lvec, sideHalf, bottomDepth) {
    const b = CrystalUtils.reciprocalLattice(lvec);
    // Side facet normals: [1,1,-1], [1,-1,1], [-1,1,1] in reciprocal space
    const sideHKLs = [[1,1,-1], [1,-1,1], [-1,1,1]];
    // Bottom normal: [-1,-1,-1]
    const bottomHKL = [-1,-1,-1];
    const planes = [];
    const BIG = 1e4;  // effectively unbounded on the "open" side
    for (const [h,k,l] of sideHKLs) {
        const n = new Vec3().setLincomb3(h, b[0], k, b[1], l, b[2]);
        const ln = n.normalize();
        if (!(ln > 0)) throw new Error(`buildTipPlanes: zero normal for [${h},${k},${l}]`);
        planes.push({ n, cmin: -BIG, cmax: sideHalf });
    }
    {
        const [h,k,l] = bottomHKL;
        const n = new Vec3().setLincomb3(h, b[0], k, b[1], l, b[2]);
        const ln = n.normalize();
        if (!(ln > 0)) throw new Error(`buildTipPlanes: zero normal for bottom`);
        planes.push({ n, cmin: -BIG, cmax: bottomDepth });
    }
    return planes;
}

function main() {
    const opts = parseArgs();
    const material = opts.material;
    const cifFile = (material === 'C' || material === 'diamond' || material === 'Diamond')
        ? path.join(RES, 'crystals/C_diamond_sym.cif')
        : path.join(RES, 'crystals/Si-sym.cif');
    const prefix = opts.prefix || `afm_tip_${material}`;
    const heavyZ = (material === 'C' || material === 'diamond' || material === 'Diamond') ? 6 : 14;

    console.log(`=== AFM Tip Nanocrystal Generator ===`);
    console.log(`Material: ${material} (Z=${heavyZ})`);
    console.log(`CIF: ${cifFile}`);
    console.log(`sideHalf=${opts.sideHalf} Å, bottomDepth=${opts.bottomDepth} Å, nRep=${opts.nRep}`);

    // Load CIF and build crystal with symmetry
    const cifText = fs.readFileSync(cifFile, 'utf8');
    const cell = buildCrystalFromCIFText(cifText, {
        applySymmetry: true,
        dedupSymmetry: true,
        dedupTol: 0.1,
    });
    console.log(`Cell: ${cell.basisTypes.length} basis atoms, lvec[0]=${cell.lvec[0]}`);

    // Load MM params
    const args = {
        elementTypes: path.join(RES, 'ElementTypes.dat'),
        atomTypes: path.join(RES, 'AtomTypes.dat'),
        bondTypes: path.join(RES, 'BondTypes.dat'),
        angleTypes: path.join(RES, 'AngleTypes.dat'),
    };
    const mm = loadMMParams(args);

    // Build custom tip planes
    const planes = buildTipPlanes(cell.lvec, opts.sideHalf, opts.bottomDepth);
    console.log(`Planes: ${planes.length} (3 side + 1 bottom)`);

    // Cut crystal with custom planes
    const n = opts.nRep;
    const mol = CrystalUtils.genReplicatedCellCutPlanes({
        lvec: cell.lvec,
        basisPos: cell.basisPos,
        basisTypes: cell.basisTypes,
        basisCharges: cell.basisCharges,
        nRep: [n, n, n],
        origin: new Vec3(0, 0, 0),
        planes,
        planeMode: 'ang',
        centered: true,
        dedup: true,
        dedupTol: 0.1,
    });
    console.log(`After cut: ${mol.atoms.length} atoms`);

    // Recalculate bonds
    mol.recalculateBonds(mm, { defaultRcut: 1.60, bondFactor: 1.30 });

    // Prune undercoordinated heavy atoms
    const pr = pruneUndercoordinatedHeavyIter(mol, heavyZ, 2, 10);
    if (pr.totalRemoved > 0) {
        console.log(`Pruned ${pr.totalRemoved} undercoordinated atoms`);
        mol.recalculateBonds(mm, { defaultRcut: 1.60, bondFactor: 1.30 });
    }

    // Assign atom types
    assignCrystalAtomTypes(mol, mm);

    // Cap with H
    const capResult = mol.addCappingAtoms(mm, 'H', {
        onlySelection: false, bBond: true, bondFactor: 1.30,
        outwardBias: 0.35, resolveClashes: true,
    });
    console.log(`Added ${capResult.nAdded} H caps`);
    if (capResult.nAdded > 0) mol.recalculateBonds(mm, { defaultRcut: 1.60, bondFactor: 1.30 });

    // Second prune pass (in case capping revealed new undercoordinated atoms)
    const pr2 = pruneUndercoordinatedHeavyIter(mol, heavyZ, 2, 10);
    if (pr2.totalRemoved > 0) {
        console.log(`Pruned ${pr2.totalRemoved} more undercoordinated atoms`);
        mol.recalculateBonds(mm, { defaultRcut: 1.60, bondFactor: 1.30 });
    }

    // Rotate so [111] → z-axis
    const dir111 = new Vec3(1, 1, 1);
    const R = Mat3.alignVectorToZ(dir111);
    CrystalUtils.rotateMoleculeInPlace(mol, R);
    console.log(`Rotated [111] → z-axis`);

    // Count atoms
    let nHeavy = 0, nH = 0;
    for (const a of mol.atoms) {
        if ((a.Z | 0) === heavyZ) nHeavy++;
        else if ((a.Z | 0) === 1) nH++;
    }
    console.log(`Final: ${mol.atoms.length} atoms (${nHeavy} ${material}, ${nH} H)`);

    // Find apex atom (highest z among heavy atoms)
    let apexZ = -Infinity, apexIdx = -1;
    for (let i = 0; i < mol.atoms.length; i++) {
        const a = mol.atoms[i];
        if ((a.Z | 0) === heavyZ && a.pos.z > apexZ) { apexZ = a.pos.z; apexIdx = i; }
    }
    if (apexIdx >= 0) {
        const apex = mol.atoms[apexIdx];
        console.log(`Apex atom at z=${apexZ.toFixed(3)} Å, pos=(${apex.pos.x.toFixed(3)}, ${apex.pos.y.toFixed(3)}, ${apex.pos.z.toFixed(3)})`);
    }

    // Write output
    fs.mkdirSync(opts.outDir, { recursive: true });
    const baseName = `${prefix}_s${opts.sideHalf.toFixed(1)}_b${opts.bottomDepth.toFixed(1)}_n${opts.nRep}`;
    const mol2Path = path.join(opts.outDir, `${baseName}.mol2`);
    const xyzPath = path.join(opts.outDir, `${baseName}.xyz`);
    fs.writeFileSync(mol2Path, toMol2String(mol));
    fs.writeFileSync(xyzPath, toXYZString(mol));
    console.log(`Written: ${mol2Path}`);
    console.log(`Written: ${xyzPath}`);
    console.log(`Done.`);
}

main();
