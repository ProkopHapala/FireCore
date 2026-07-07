/// Bootstrap NPZ viewer fixtures (Agent A) until Agent C delivers production exports.
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { writeCrystalNpz, writeNpzCompressed } from '../../../../web/common_js/npzIO.js';

const __dir = path.dirname(fileURLToPath(import.meta.url));

// 2-atom Si dimer crystal
const pos = new Float64Array([0, 0, 0, 0, 0, 2.35]);
const Z = new Int32Array([14, 14]);
const bonds_ij = new Int32Array([0, 1]);
writeCrystalNpz(fs, path.join(__dir, '01_init.npz'), { pos, Z, bonds_ij, gen_params: '{"fixture":"agent_a_bootstrap"}' });

// topology NPZ with AABB overlay (2 surface groups)
const n = 4;
const pos4 = new Float64Array([
    0, 0, 0,  1.5, 0, 0,
    0, 1.5, 0,  1.5, 1.5, 0,
]);
const Z4 = new Int32Array([14, 14, 14, 14]);
pos4._shape = [n, 3];
Z4._shape = [n];
const icolGroup = new Int32Array([0, 0, 1, 1]);
icolGroup._shape = [n];
const group_bbox_min = new Float64Array([
    -0.5, -0.5, -0.5,
    1.0, 1.0, -0.5,
]);
const group_bbox_max = new Float64Array([
    2.0, 0.5, 0.5,
    2.5, 2.5, 0.5,
]);
group_bbox_min._shape = [2, 3];
group_bbox_max._shape = [2, 3];
const n_groups = new Int32Array([2]);
const natoms = new Int32Array([n]);
writeNpzCompressed(path.join(__dir, '03_topology.npz'), {
    pos: pos4, Z: Z4, natoms, icolGroup, group_bbox_min, group_bbox_max, n_groups,
}, fs);

// loose .npy fast path
import { encodeNpy } from '../../../../web/common_js/npzIO.js';
fs.writeFileSync(path.join(__dir, 'pos.npy'), encodeNpy(pos4, [n, 3]));
fs.writeFileSync(path.join(__dir, 'Z.npy'), encodeNpy(Z4, [n]));

console.log('bootstrap_fixtures: wrote 01_init.npz, 03_topology.npz, pos.npy, Z.npy');
