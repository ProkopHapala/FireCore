import fs from 'node:fs';
import path from 'node:path';

import './headless_init.mjs';
import { buildXPDBInputsFromXYZArgs } from './dump_xpdb_topology.mjs';

import { XPDB_WebGPU } from './XPDB_WebGPU.js';


/*

run like this:

node web/molgui_webgpu/run_xpdb_webgpu_headless.mjs --xyz tests/tDFT_pentacene/pentacene.xyz --out web/molgui_webgpu/out_webgpu_pentacene_distort.txt --inner_iters 10 --n_steps 10 --traj_xyz web/molgui_webgpu/traj_webgpu_pentacene.xyz --dt 0.01 --omega 0.8 --momentum_beta 0 --k_coll 300 --coll_scale 1.2 --bbox_scale 1.2 --maxRadius 0.2 --init_scale 1.05 --init_noise 0.02 --init_seed 12345 --twoPi 0 --addPiAlign 0

PYTHONPATH=/home/prokop/git/FireCore:$PYTHONPATH python3 pyBall/XPDB_AVBD/test_XPDB_new_dump_headless.py --molecule tests/tDFT_pentacene/pentacene.xyz --dump_out pyBall/XPDB_AVBD/out_ocl_pentacene_distort.txt --inner_iters 10 --n_steps 10 --traj_xyz pyBall/XPDB_AVBD/traj_ocl_pentacene.xyz --traj_real_only 1 --dt 0.01 --omega 0.8 --momentum_beta 0 --k_coll 300 --coll_scale 1.2 --bbox_scale 1.2 --init_scale 1.05 --init_noise 0.02 --init_seed 12345 --two_pi 0 --add_pi_align 0

python3 pyBall/XPDB_AVBD/compare_xyz_trajectories.py pyBall/XPDB_AVBD/traj_ocl_pentacene.xyz web/molgui_webgpu/traj_webgpu_pentacene.xyz --atol 1e-5

*/

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
        coll_scale: 2.0,
        bbox_scale: 2.0,
        k_coll: 200.0,
        omega: 0.0,
        momentum_beta: 0.0,
        inner_iters: 0,

        n_steps: 0,
        traj_xyz: null,
        traj_real_only: true,
        init_scale: 1.0,
        init_noise: 0.0,
        init_seed: 1337,
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
        else if (a === '--coll_scale') out.coll_scale = +nxt();
        else if (a === '--bbox_scale') out.bbox_scale = +nxt();
        else if (a === '--k_coll') out.k_coll = +nxt();
        else if (a === '--omega') out.omega = +nxt();
        else if (a === '--momentum_beta') out.momentum_beta = +nxt();
        else if (a === '--inner_iters') out.inner_iters = parseInt(nxt(), 10) | 0;
        else if (a === '--n_steps') out.n_steps = parseInt(nxt(), 10) | 0;
        else if (a === '--traj_xyz') out.traj_xyz = nxt();
        else if (a === '--traj_real_only') out.traj_real_only = (nxt() !== '0');
        else if (a === '--init_scale') out.init_scale = +nxt();
        else if (a === '--init_noise') out.init_noise = +nxt();
        else if (a === '--init_seed') out.init_seed = parseInt(nxt(), 10) | 0;
        else throw new Error(`Unknown arg ${a}`);
    }
    if (!out.xyz) throw new Error('Usage: node run_xpdb_webgpu_headless.mjs --xyz <molecule.xyz> --out <dump.txt>');
    if (!out.out) throw new Error('Usage: node run_xpdb_webgpu_headless.mjs --xyz <molecule.xyz> --out <dump.txt>');
    return out;
}

function writeXYZFrame(stream, symbols, pos4, nAtoms, step, dt, title) {
    stream.write(`${nAtoms}\n`);
    stream.write(`${title} step=${step} time=${(step * dt).toFixed(6)}\n`);
    for (let i = 0; i < nAtoms; i++) {
        const s = symbols[i];
        const x = pos4[i * 4 + 0];
        const y = pos4[i * 4 + 1];
        const z = pos4[i * 4 + 2];
        stream.write(`${s} ${x.toFixed(8)} ${y.toFixed(8)} ${z.toFixed(8)}\n`);
    }
}

function deterministicNoise(idx, axis, seed) {
    const x = Math.sin(seed * 12.9898 + idx * 78.233 + axis * 37.719) * 43758.5453;
    const frac = x - Math.floor(x);
    return frac * 2.0 - 1.0;
}

function applyInitialDistortion(pos4, nAtoms, { scale = 1.0, noise = 0.0, seed = 1337 } = {}) {
    if (scale !== 1.0) {
        for (let i = 0; i < nAtoms; i++) {
            const base = i * 4;
            pos4[base + 0] *= scale;
            pos4[base + 1] *= scale;
            pos4[base + 2] *= scale;
        }
    }
    if (noise !== 0.0) {
        for (let i = 0; i < nAtoms; i++) {
            const base = i * 4;
            pos4[base + 0] += noise * deterministicNoise(i, 0, seed);
            pos4[base + 1] += noise * deterministicNoise(i, 1, seed);
            pos4[base + 2] += noise * deterministicNoise(i, 2, seed);
        }
    }
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

    // Apply deterministic distortion to initial coordinates while keeping bond rest lengths from original inputs
    const pos4 = new Float32Array(inputs.pos4);
    applyInitialDistortion(pos4, inputs.nAtoms, { scale: args.init_scale, noise: args.init_noise, seed: args.init_seed });

    // Upload buffers with distorted positions
    sim.device.queue.writeBuffer(sim.buffers.pos, 0, pos4);
    sim.device.queue.writeBuffer(sim.buffers.predPos, 0, pos4);
    sim.device.queue.writeBuffer(sim.buffers.atomParams, 0, inputs.params4);
    sim.device.queue.writeBuffer(sim.buffers.bondIdxGlobal, 0, inputs.bondIndices);
    sim.device.queue.writeBuffer(sim.buffers.bondLenStiff, 0, inputs.bondLenStiff);

    // Optional trajectory: run multiple steps and dump pos each frame.
    if ((args.n_steps | 0) > 0) {
        if (!args.traj_xyz) throw new Error('--n_steps > 0 requires --traj_xyz PATH');
        const Z_TO_SYMBOL = inputs.mol.constructor.Z_TO_SYMBOL;
        if (!Z_TO_SYMBOL) throw new Error('Z_TO_SYMBOL missing; MoleculeIO not installed?');
        const topo = inputs.topo;
        const nReal = topo && (topo.n_real !== undefined) ? (topo.n_real | 0) : inputs.nAtoms;
        const nWrite = args.traj_real_only ? nReal : inputs.nAtoms;
        const symbols = new Array(nWrite);
        for (let i = 0; i < nWrite; i++) {
            const a = inputs.mol.atoms[i];
            const Z = a ? (a.Z | 0) : 0;
            symbols[i] = Z_TO_SYMBOL[Z] || 'X';
        }
        const stream = fs.createWriteStream(args.traj_xyz, { encoding: 'utf8' });
        for (let istep = 0; istep < (args.n_steps | 0); istep++) {
            sim.step({
                dt: args.dt,
                iterations: args.inner_iters,
                k_coll: args.k_coll,
                omega: args.omega,
                momentum_beta: args.momentum_beta,
                maxRadius: args.maxRadius,
                coll_scale: args.coll_scale,
                bbox_scale: args.bbox_scale,
            });
            const pos = await sim.readBuffer(sim.buffers.pos).then(buf => new Float32Array(buf));
            writeXYZFrame(stream, symbols, pos, nWrite, istep, args.dt, 'XPDB_WebGPU');
        }
        await new Promise((resolve, reject) => {
            stream.end(() => resolve());
            stream.on('error', (e) => reject(e));
        });
    } else {
        // Run: update_bboxes + build_local_topology (+ solve with 0 iters, no effect)
        sim.step({
            dt: args.dt,
            iterations: args.inner_iters,
            k_coll: args.k_coll,
            omega: args.omega,
            momentum_beta: args.momentum_beta,
            maxRadius: args.maxRadius,
            coll_scale: args.coll_scale,
            bbox_scale: args.bbox_scale,
        });
    }

    await XPDB_WebGPU.dumpTypedState(sim, {
        label: `XPDB_WebGPU headless xyz=${inputs.xyzPath}`,
        fixed: args.dump_fixed,
        writeText: async (text) => {
            fs.writeFileSync(args.out, text, 'utf8');
        }
    });
    console.log(`Wrote ${args.out}`);
    if ((args.n_steps | 0) > 0) console.log(`Wrote ${args.traj_xyz}`);
}

main();
