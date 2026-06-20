/// Build topology NPZ from mol2 + relaxed xyz: MMFFL topology + collision workgroups + exclusions.
import fs from 'node:fs';
import { buildMMFFLTopology, packLinearTopologyForGPU } from '../molgui_webgpu/MMFFLTopology.js';
import { writeTopologyNpzFile } from '../molgui_webgpu/LinearizedTopologyNpz.js';
import { buildCollisionWorkgroups, buildExclIcol_1_2_3 } from './CollisionWorkgroups.js';
import { readXyzPositions } from './npzIO.js';
import { loadMMParamsFromDir, loadMolFromMol2, applyPositions } from './MolIO.js';

export function buildTopologyNpz({ mol2Path, relaxedXyzPath, outNpzPath, mm = null, maxNeighbors = 48, addAngle = false, addDihedral = false, groupCap = 32 }) {
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

    const nAtoms = topo.n_real | 0;
    const posFlat = new Float64Array(nAtoms * 3);
    for (let i = 0; i < nAtoms; i++) {
        const q = topo.apos[i];
        posFlat[i * 3 + 0] = +q[0];
        posFlat[i * 3 + 1] = +q[1];
        posFlat[i * 3 + 2] = +q[2];
    }
    const radius = new Float64Array(nAtoms);
    for (let i = 0; i < nAtoms; i++) {
        const a = mol.atoms[i];
        const at = mmParams.getAtomTypeForAtom(a);
        const r = (at && at.RvdW > 0) ? at.RvdW : ((at && at.element && at.element.RvdW > 0) ? at.element.RvdW : 1.5);
        if (!(r > 0)) throw new Error(`buildTopologyNpz: invalid RvdW=${r} at atom ${i}`);
        radius[i] = +r;
    }
    radius._shape = [nAtoms];

    const wg = buildCollisionWorkgroups({ pos: posFlat, mol, radius, groupCap, fillFactor: 0.8 });
    const ex2 = buildExclIcol_1_2_3({ mol, icol: wg.icol, EXCL_MAX: 16, ipbc: 0 });
    wg.group_atoms._shape = [wg.nGroups, wg.groupCap];
    wg.group_nAtoms._shape = [wg.nGroups];
    wg.icol._shape = [nAtoms];
    wg.icolGroup._shape = [nAtoms];
    wg.bbox_min._shape = [wg.nGroups, 3];
    wg.bbox_max._shape = [wg.nGroups, 3];
    ex2.excl._shape = [nAtoms, ex2.EXCL_MAX];
    ex2.excl_count._shape = [nAtoms];
    const meta = {
        source: mol2Path,
        n_bond: topo.bonds_primary.length | 0,
        n_angle: topo.bonds_angle.length | 0,
        n_dihedral: (topo.bonds_dihedral ? topo.bonds_dihedral.length : 0) | 0,
    };
    const extra = {
        radius,
        icol: wg.icol,
        icolGroup: wg.icolGroup,
        group_atoms: wg.group_atoms,
        group_nAtoms: wg.group_nAtoms,
        group_bbox_min: wg.bbox_min,
        group_bbox_max: wg.bbox_max,
        excl_icol: ex2.excl,
        excl_count: ex2.excl_count,
        excl_max: new Int32Array([ex2.EXCL_MAX | 0]),
        group_cap: new Int32Array([wg.groupCap | 0]),
        n_groups: new Int32Array([wg.nGroups | 0]),
    };
    writeTopologyNpzFile(fs, outNpzPath, topo, packing, mol, meta, extra);
    return { natoms: topo.n_real | 0, meta, packing };
}
