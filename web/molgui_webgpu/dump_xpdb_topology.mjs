import fs from 'node:fs';
import path from 'node:path';

import { MMParams } from './MMParams.js';
import { EditableMolecule } from './EditableMolecule.js';
import { installMoleculeIOMethods } from './MoleculeIO.js';
import { installMoleculeUtilsMethods } from './MoleculeUtils.js';
import { installMoleculeSelectionMethods } from './MoleculeSelection.js';
import { buildXPDBTopology } from './XPDBTopology.js';

installMoleculeIOMethods(EditableMolecule);
installMoleculeUtilsMethods(EditableMolecule);
installMoleculeSelectionMethods(EditableMolecule);

function parseArgs(argv) {
    const out = {
        xyz: null,
        enableAngles: true,
        maxBonds: 16,
        defaultL: 1.5,
        defaultK: 100.0,
        elementTypes: '../../cpp/common_resources/ElementTypes.dat',
        atomTypes: '../../cpp/common_resources/AtomTypes.dat',
        bondTypes: '../../cpp/common_resources/BondTypes.dat',
        angleTypes: '../../cpp/common_resources/AngleTypes.dat',
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
        else if (a === '--maxBonds') out.maxBonds = parseInt(nxt(), 10) | 0;
        else if (a === '--defaultL') out.defaultL = +nxt();
        else if (a === '--defaultK') out.defaultK = +nxt();
        else if (a === '--elementTypes') out.elementTypes = nxt();
        else if (a === '--atomTypes') out.atomTypes = nxt();
        else if (a === '--bondTypes') out.bondTypes = nxt();
        else if (a === '--angleTypes') out.angleTypes = nxt();
        else throw new Error(`Unknown arg ${a}`);
    }
    if (!out.xyz) throw new Error('Usage: node dump_xpdb_topology.mjs --xyz path/to/file.xyz [--angles 0|1]');
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

async function main() {
    const args = parseArgs(process.argv);

    const mm = new MMParams();
    const readText = (p) => fs.readFileSync(p, 'utf8');

    mm.parseElementTypes(readText(args.elementTypes));
    mm.parseAtomTypes(readText(args.atomTypes));
    mm.parseBondTypes(readText(args.bondTypes));
    mm.parseAngleTypes(readText(args.angleTypes));

    const xyzText = fs.readFileSync(args.xyz, 'utf8');
    const parsed = EditableMolecule.parseXYZ(xyzText);

    const mol = new EditableMolecule();
    mol.clear();
    mol.appendParsedSystem(parsed);
    mol.recalculateBonds(mm);

    const { bondsAdj, stats } = buildXPDBTopology(mol, mm, {
        includeAngleConstraints: !!args.enableAngles,
        maxBonds: args.maxBonds,
        defaultL: args.defaultL,
        defaultK: args.defaultK,
    });

    const nAtoms = mol.nAtoms || mol.atoms.length;
    const nMaxBonded = 16;
    const { bondIndices, bondLenStiff } = packBondArrays(bondsAdj, nAtoms, nMaxBonded);

    const outLines = [];
    outLines.push(`[DEBUG] dump_xpdb_topology.mjs xyz=${args.xyz}`);
    outLines.push(`[DEBUG] stats=${JSON.stringify(stats)}`);
    outLines.push(`[DEBUG] nAtoms=${nAtoms} nMaxBonded=${nMaxBonded} enableAngles=${args.enableAngles}`);
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
}

main();
