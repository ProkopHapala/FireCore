/// Pack linearized MMFFL topology into numpy .npz (Node, binary — no Python).
import { writeNpzCompressed } from '../common_js/npzWrite.js';

/// Typed arrays ready for encodeNpzCompressed. Keeps packing buffers by reference where possible.
export function buildTopologyNpzArrays(topo, packing, mol, meta = {}) {
    const n = packing.nAtoms | 0;
    const m = packing.maxNeighbors | 0;
    if (!(n > 0) || !(m > 0)) throw new Error(`buildTopologyNpzArrays: invalid n=${n} m=${m}`);
    const pos = new Float64Array(n * 3);
    const Z = new Int32Array(n);
    for (let i = 0; i < n; i++) {
        const p = topo.apos[i];
        pos[i * 3 + 0] = +p[0];
        pos[i * 3 + 1] = +p[1];
        pos[i * 3 + 2] = +p[2];
        Z[i] = mol.atoms[i].Z | 0;
    }
    const neigh_idx = packing.neighs instanceof Int32Array ? packing.neighs : Int32Array.from(packing.neighs);
    const stick_class = packing.stick_class instanceof Int32Array ? packing.stick_class : Int32Array.from(packing.stick_class);
    const KLsFlat = packing.KLs instanceof Float32Array ? packing.KLs : Float32Array.from(packing.KLs);
    const bond_l0 = new Float32Array(n * m);
    const bond_k = new Float32Array(n * m);
    const neigh_count = new Int32Array(n);
    for (let i = 0; i < n; i++) {
        let c = 0;
        const base = i * m;
        for (let k = 0; k < m; k++) {
            const j = neigh_idx[base + k] | 0;
            bond_l0[base + k] = KLsFlat[(base + k) * 2 + 0];
            bond_k[base + k] = KLsFlat[(base + k) * 2 + 1];
            if (j >= 0) c++;
        }
        neigh_count[i] = c;
    }
    const KLs = KLsFlat;
    KLs._shape = [n, m, 2];
    pos._shape = [n, 3];
    Z._shape = [n];
    neigh_idx._shape = [n, m];
    stick_class._shape = [n, m];
    bond_l0._shape = [n, m];
    bond_k._shape = [n, m];
    neigh_count._shape = [n];
    const natoms = new Int32Array([n]);
    const max_neighbors = new Int32Array([m]);
    const n_bond = new Int32Array([meta.n_bond | 0]);
    const n_angle = new Int32Array([meta.n_angle | 0]);
    const n_dihedral = new Int32Array([meta.n_dihedral | 0]);
    return {
        pos, Z, neigh_idx, neigh_count, bond_l0, bond_k, KLs, stick_class,
        natoms, max_neighbors, n_bond, n_angle, n_dihedral,
    };
}

export function writeTopologyNpzFile(fs, outPath, topo, packing, mol, meta = {}) {
    const arrays = buildTopologyNpzArrays(topo, packing, mol, meta);
    writeNpzCompressed(outPath, arrays, fs);
}
