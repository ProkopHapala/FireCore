/// Build linearized topology NPZ from mol2 + relaxed xyz (nanocrystal pipeline).
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import { MMParams } from '../molgui_webgpu/MMParams.js';
import { EditableMolecule } from '../molgui_webgpu/EditableMolecule.js';
import { installMoleculeIOMethods } from '../molgui_webgpu/MoleculeIO.js';
import { installMoleculeUtilsMethods } from '../molgui_webgpu/MoleculeUtils.js';
import { buildMMFFLTopology, packLinearTopologyForGPU } from '../molgui_webgpu/MMFFLTopology.js';
import { writeTopologyNpzFile } from '../molgui_webgpu/LinearizedTopologyNpz.js';
import { readXyzPositions, molToCrystalArrays } from './npzIO.js';

installMoleculeIOMethods(EditableMolecule);
installMoleculeUtilsMethods(EditableMolecule);

const defaultRes = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../../cpp/common_resources');

export function loadMMParamsFromDir(resDir = defaultRes) {
    const readText = (p) => fs.readFileSync(p, 'utf8');
    const mm = new MMParams();
    mm.parseElementTypes(readText(path.join(resDir, 'ElementTypes.dat')));
    mm.parseAtomTypes(readText(path.join(resDir, 'AtomTypes.dat')));
    mm.parseBondTypes(readText(path.join(resDir, 'BondTypes.dat')));
    mm.parseAngleTypes(readText(path.join(resDir, 'AngleTypes.dat')));
    mm.parseDihedralTypes(readText(path.join(resDir, 'DihedralTypes.dat')));
    return mm;
}

export function loadMolFromMol2(mol2Path, mm) {
    const text = fs.readFileSync(mol2Path, 'utf8');
    const mol = new EditableMolecule();
    mol.clear();
    mol.appendParsedSystem(EditableMolecule.parseMol2(text));
    if (!mol.bonds || mol.bonds.length === 0) mol.recalculateBonds(mm);
    return mol;
}

export function applyPositions(mol, pos) {
    const n = mol.atoms.length;
    if (pos.length < n * 3) throw new Error(`applyPositions: pos length ${pos.length} < ${n * 3}`);
    for (let i = 0; i < n; i++) {
        const a = mol.atoms[i];
        a.pos.x = pos[i * 3];
        a.pos.y = pos[i * 3 + 1];
        a.pos.z = pos[i * 3 + 2];
    }
}

/// Bonds for SVG: all mol2 bonds except H–H (capHH topology bonds confuse the view).
export function bondsForVisualization(mol) {
    const { Z, bonds_ij } = molToCrystalArrays(mol);
    if (!bonds_ij || !bonds_ij.length) return null;
    const out = [];
    const nb = bonds_ij.length / 2;
    for (let k = 0; k < nb; k++) {
        const a = bonds_ij[k * 2] | 0, b = bonds_ij[k * 2 + 1] | 0;
        if ((Z[a] | 0) === 1 && (Z[b] | 0) === 1) continue;
        out.push(a, b);
    }
    return out.length ? Int32Array.from(out) : null;
}

export function getCrystalBondsFromFiles(mol2Path, xyzPath, mm = null) {
    const mmParams = mm || loadMMParamsFromDir();
    const mol = loadMolFromMol2(mol2Path, mmParams);
    const { pos } = readXyzPositions(fs, xyzPath);
    applyPositions(mol, pos);
    return bondsForVisualization(mol);
}

export function buildTopologyNpz({ mol2Path, relaxedXyzPath, outNpzPath, mm = null, maxNeighbors = 48, addAngle = false, addDihedral = false }) {
    const mmParams = mm || loadMMParamsFromDir();
    const mol = loadMolFromMol2(mol2Path, mmParams);
    if (relaxedXyzPath) {
        const { pos } = readXyzPositions(fs, relaxedXyzPath);
        applyPositions(mol, pos);
    }
    const topo = buildMMFFLTopology(mol, mmParams, {
        type_source: 'table',
        add_angle: addAngle,
        add_dihedral: addDihedral,
        add_pi: false,
        add_epair: false,
        K12: 0, K13: 0, K14: 0,
    });
    const packing = packLinearTopologyForGPU(topo, { maxNeighbors });
    const meta = {
        source: mol2Path,
        n_bond: topo.bonds_primary.length | 0,
        n_angle: topo.bonds_angle.length | 0,
        n_dihedral: (topo.bonds_dihedral ? topo.bonds_dihedral.length : 0) | 0,
    };
    writeTopologyNpzFile(fs, outNpzPath, topo, packing, mol, meta);
    return { natoms: topo.n_real | 0, meta, packing };
}
