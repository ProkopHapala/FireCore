import fs from 'node:fs';
import path from 'node:path';

import './headless_init.mjs';
import { buildXPDBInputsFromXYZArgs } from './dump_xpdb_topology.mjs';

import { XPDB_WebGPU } from './XPDB_WebGPU.js';

function parseArgs(argv) {
    const out = {
        xyz: null,
        out: null,
        // keep consistent defaults with dump_xpdb_topology.mjs
        enableAngles: true,
        useMMFFL: true,
        addPi: true,
        twoPi: true,
        addPiAlign: true,
        addEpair: true,
        addEpairPairs: true,
        reportTypes: false,
        L_pi: 1.0,
        L_epair: 0.5,
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
        atom_rad: 0.2,
        atom_mass: 1.0,
        dump_fixed: 6,
        alignPiVectors: false,
        elementTypes: path.resolve('cpp/common_resources/ElementTypes.dat'),
        atomTypes: path.resolve('cpp/common_resources/AtomTypes.dat'),
        bondTypes: path.resolve('cpp/common_resources/BondTypes.dat'),
        angleTypes: path.resolve('cpp/common_resources/AngleTypes.dat'),
        dt: 0.01,
        maxRadius: 0.2,
        k_coll: 0.0,
        omega: 0.0,
        momentum_beta: 0.0,
        inner_iters: 0,
    };

    for (let i = 2; i < argv.length; i++) {
        const a = String(argv[i]);
        const nxt = () => {
            if (i + 1 >= argv.length) throw new Error(`Missing value after ${a}`);
            i++;
            return String(argv[i]);
        };
        if (a === '--xyz') out.xyz = nxt();
        else if (a === '--out') out.out = nxt();
        else if (a === '--angles') out.enableAngles = (nxt() !== '0');
        else if (a === '--useMMFFL') out.useMMFFL = (nxt() !== '0');
        else if (a === '--addPi') out.addPi = (nxt() !== '0');
        else if (a === '--twoPi') out.twoPi = (nxt() !== '0');
        else if (a === '--addPiAlign') out.addPiAlign = (nxt() !== '0');
        else if (a === '--addEpair') out.addEpair = (nxt() !== '0');
        else if (a === '--addEpairPairs') out.addEpairPairs = (nxt() !== '0');
        else if (a === '--align_pi_vectors') out.alignPiVectors = (nxt() !== '0');
        else if (a === '--atom_rad') out.atom_rad = +nxt();
        else if (a === '--atom_mass') out.atom_mass = +nxt();
        else if (a === '--defaultK') out.defaultK = +nxt();
        else if (a === '--dt') out.dt = +nxt();
        else if (a === '--maxRadius') out.maxRadius = +nxt();
        else if (a === '--k_coll') out.k_coll = +nxt();
        else if (a === '--omega') out.omega = +nxt();
        else if (a === '--momentum_beta') out.momentum_beta = +nxt();
        else if (a === '--inner_iters') out.inner_iters = parseInt(nxt(), 10) | 0;
        else throw new Error(`Unknown arg ${a}`);
    }
    if (!out.xyz) throw new Error('Usage: node run_xpdb_webgpu_headless.mjs --xyz <molecule.xyz> --out <dump.txt>');
    if (!out.out) throw new Error('Usage: node run_xpdb_webgpu_headless.mjs --xyz <molecule.xyz> --out <dump.txt>');
    return out;
}

async function main() {
    const args = parseArgs(process.argv);

    const inputs = await buildXPDBInputsFromXYZArgs({
        xyz: args.xyz,
        enableAngles: args.enableAngles,
        useMMFFL: args.useMMFFL,
        addPi: args.addPi,
        twoPi: args.twoPi,
        addPiAlign: args.addPiAlign,
        addEpair: args.addEpair,
        addEpairPairs: args.addEpairPairs,
        reportTypes: args.reportTypes,
        L_pi: args.L_pi,
        L_epair: args.L_epair,
        k_angle: args.k_angle,
        k_pi: args.k_pi,
        k_pi_orth: args.k_pi_orth,
        k_pi_align: args.k_pi_align,
        k_ep: args.k_ep,
        k_ep_orth: args.k_ep_orth,
        k_ep_pair: args.k_ep_pair,
        maxBonds: args.maxBonds,
        defaultL: args.defaultL,
        defaultK: args.defaultK,
        atom_rad: args.atom_rad,
        atom_mass: args.atom_mass,
        dump_fixed: args.dump_fixed,
        alignPiVectors: args.alignPiVectors,
        elementTypes: args.elementTypes,
        atomTypes: args.atomTypes,
        bondTypes: args.bondTypes,
        angleTypes: args.angleTypes,
        dump_topo_debug: null,
    });

    const sim = new XPDB_WebGPU(inputs.nAtoms, 64);
    await sim.init();

    // Upload buffers exactly as built by buildXPDBInputsFromXYZArgs()
    sim.device.queue.writeBuffer(sim.buffers.pos, 0, inputs.pos4);
    sim.device.queue.writeBuffer(sim.buffers.predPos, 0, inputs.pos4);
    sim.device.queue.writeBuffer(sim.buffers.atomParams, 0, inputs.params4);
    sim.device.queue.writeBuffer(sim.buffers.bondIdxGlobal, 0, inputs.bondIndices);
    sim.device.queue.writeBuffer(sim.buffers.bondLenStiff, 0, inputs.bondLenStiff);

    // Run: update_bboxes + build_local_topology (+ solve with 0 iters, no effect)
    sim.step(args.dt, args.inner_iters, args.k_coll, args.omega, args.momentum_beta, null, -1, args.maxRadius);

    await XPDB_WebGPU.dumpTypedState(sim, {
        label: `XPDB_WebGPU headless xyz=${inputs.xyzPath}`,
        fixed: args.dump_fixed,
        writeText: async (text) => {
            fs.writeFileSync(args.out, text, 'utf8');
        }
    });
    console.log(`Wrote ${args.out}`);
}

main();
