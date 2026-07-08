#!/usr/bin/env node
/// @deprecated Use `node tests/tSiNCs/nanocrystals.mjs nonbond` instead. This script is kept for backward compatibility.
import fs from 'node:fs';
import path from 'node:path';

import { loadMMParamsFromDir, loadMolFromMol2, applyPositions } from '../../web/common_js/MolIO.js';
import { readXyzPositions, readNpzFile } from '../../web/common_js/npzIO.js';
import { buildCollisionWorkgroups, buildExclIcol_1_2_3 } from '../../web/common_js/CollisionWorkgroups.js';
import { computeNonbondBruteForceKernelStyle, computeNonbondByGroups, maxAbsDiff } from '../../web/common_js/Nonbonded.js';

function usage() {
    console.log(`Usage: debug_nanocrystal_nonbond_groups.mjs --mol2 MOL2 --xyz relaxed.xyz [--topo topology.npz] [--margin 0.0] [--collision] [--rcut R] [--rcut2 R] [--out-pairs FILE]`);
}

function parseArgs(argv) {
    const out = { mol2: null, xyz: null, topo: null, margin: 0.0, outPairs: null, collision: false, rcut: 0, rcut2: 0 };
    for (let i = 2; i < argv.length; i++) {
        const a = argv[i];
        if (a === '--mol2') out.mol2 = argv[++i];
        else if (a === '--xyz') out.xyz = argv[++i];
        else if (a === '--topo') out.topo = argv[++i];
        else if (a === '--margin') out.margin = parseFloat(argv[++i]);
        else if (a === '--out-pairs') out.outPairs = argv[++i];
        else if (a === '--collision') out.collision = true;
        else if (a === '--rcut') out.rcut = parseFloat(argv[++i]);
        else if (a === '--rcut2') out.rcut2 = parseFloat(argv[++i]);
        else if (a === '--help') { usage(); process.exit(0); }
        else throw new Error(`Unknown arg '${a}'`);
    }
    if (!out.mol2 || !out.xyz) throw new Error('Missing --mol2 or --xyz');
    return out;
}

function main() {
    const args = parseArgs(process.argv);
    const mm = loadMMParamsFromDir();
    const mol = loadMolFromMol2(args.mol2, mm);
    const { pos } = readXyzPositions(fs, args.xyz);
    applyPositions(mol, pos);

    const nAtoms = mol.atoms.length | 0;
    const posFlat = new Float64Array(nAtoms * 3);
    for (let i = 0; i < nAtoms; i++) {
        const a = mol.atoms[i];
        posFlat[i * 3 + 0] = +a.pos.x;
        posFlat[i * 3 + 1] = +a.pos.y;
        posFlat[i * 3 + 2] = +a.pos.z;
    }

    let radius = new Float64Array(nAtoms);
    for (let i = 0; i < nAtoms; i++) {
        const at = mm.getAtomTypeForAtom(mol.atoms[i]);
        const r = (at && at.RvdW > 0) ? at.RvdW : 1.5;
        if (!(r > 0)) throw new Error(`invalid RvdW ${r} at atom ${i}`);
        radius[i] = +r;
    }

    let icol, group_atoms, group_nAtoms, bbox_min, bbox_max, excl_icol;
    let groupCap = 32;
    if (args.topo) {
        const { arrays } = readNpzFile(fs, args.topo);
        if (!arrays.icol || !arrays.group_atoms || !arrays.group_nAtoms || !arrays.group_bbox_min || !arrays.group_bbox_max || !arrays.excl_icol) {
            throw new Error(`topology npz missing required arrays: ${args.topo}`);
        }
        icol = arrays.icol;
        group_atoms = arrays.group_atoms;
        group_nAtoms = arrays.group_nAtoms;
        bbox_min = arrays.group_bbox_min;
        bbox_max = arrays.group_bbox_max;
        excl_icol = arrays.excl_icol;
        groupCap = arrays.group_cap ? (arrays.group_cap[0] | 0) : 32;
    } else {
        const wg = buildCollisionWorkgroups({ pos: posFlat, mol, radius, groupCap: 32, fillFactor: 0.8 });
        icol = wg.icol;
        group_atoms = wg.group_atoms;
        group_nAtoms = wg.group_nAtoms;
        groupCap = wg.groupCap;
        bbox_min = wg.bbox_min;
        bbox_max = wg.bbox_max;
        const ex2 = buildExclIcol_1_2_3({ mol, icol, EXCL_MAX: 16, ipbc: 0 });
        excl_icol = ex2.excl;
    }

    const brute = computeNonbondBruteForceKernelStyle({ pos: posFlat, mol, mmParams: mm, surfaceOnly: true, exclude12and13: true, collectPairs: true });
    const grouped = computeNonbondByGroups({ pos: posFlat, mol, mmParams: mm, group_atoms, group_nAtoms, group_bbox_min: bbox_min, group_bbox_max: bbox_max, radius, excl_icol, EXCL_MAX: 16, groupCap, icol, surfaceOnly: true, margin: args.margin, collectPairs: true, collisionMode: args.collision, R_cut: args.rcut, R_cut2: args.rcut2 });

    const df = maxAbsDiff(brute.f, grouped.f);
    const dE = Math.abs(brute.Etot - grouped.Etot);

    const bruteSet = new Set(brute.pairs.map(([a, b]) => a < b ? `${a},${b}` : `${b},${a}`));
    const groupedSet = new Set(grouped.pairs.map(([a, b]) => a < b ? `${a},${b}` : `${b},${a}`));
    const onlyInBrute = [...bruteSet].filter(k => !groupedSet.has(k));
    const onlyInGrouped = [...groupedSet].filter(k => !bruteSet.has(k));

    console.log(`# nonbond parity`);
    console.log(`mol2: ${path.resolve(args.mol2)}`);
    console.log(`xyz:  ${path.resolve(args.xyz)}`);
    if (args.topo) console.log(`topo: ${path.resolve(args.topo)}`);
    console.log(`mode: ${args.collision ? 'collision (getSR_x2_smooth)' : 'full LJ+Q'}`);
    console.log(`margin: ${args.margin}`);
    if (args.collision) console.log(`R_cut: ${args.rcut}, R_cut2: ${args.rcut2}`);
    console.log(`maxAbsDiff(force): ${df}`);
    console.log(`absDiff(E): ${dE}`);
    console.log(`brute pairs: ${brute.pairs.length}, grouped pairs: ${grouped.pairs.length}`);
    console.log(`only in brute: ${onlyInBrute.length}, only in grouped: ${onlyInGrouped.length}`);
    if (onlyInBrute.length > 0) console.log(`  onlyInBrute (first 10): ${onlyInBrute.slice(0, 10).join(' ')}`);
    if (onlyInGrouped.length > 0) console.log(`  onlyInGrouped (first 10): ${onlyInGrouped.slice(0, 10).join(' ')}`);

    if (args.outPairs) {
        const icolGroupArr = args.topo ? readNpzFile(fs, args.topo).arrays.icolGroup : icol;
        const out = {
            nAtoms,
            groupCap,
            icolGroup: Array.from(icolGroupArr),
            pairs: grouped.pairs.map(([a, b]) => [a, b]),
            pairGroups: grouped.pairs.map(([a, b]) => [icolGroupArr[a] | 0, icolGroupArr[b] | 0]),
        };
        fs.writeFileSync(args.outPairs, JSON.stringify(out));
        console.log(`wrote collision pairs to ${args.outPairs}`);
    }
}

main();
