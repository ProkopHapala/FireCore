#!/usr/bin/env python3
"""Consolidated headless test script for RRsp3.

Replaces: test_RRsp3_momentum.py, test_RRsp3_smoke.py, test_RRsp3_debug.py,
          test_RRsp3_debug_runner.py, test_RRsp3_convergence.py

Subcommands:
    momentum    -- test linear and angular momentum conservation
    topology    -- test collision/exclusion mapping via kernel debug prints
    smoke       -- basic smoke test (step completes, no NaN, ranges OK)
    convergence -- relax a stretched XYZ sheet and compare kernel convergence
    suite       -- run full test suite across systems and kernels

Examples:
    python test_RRsp3_headless.py momentum --system L2 --kernel all --steps 50
    python test_RRsp3_headless.py topology --system L1 --shift 1.5
    python test_RRsp3_headless.py smoke --system L0
    python test_RRsp3_headless.py convergence --xyz pentacene.xyz --kernel all --stretch 1.05
    python test_RRsp3_headless.py suite --level L3 --kernels current,shapematch
"""

import sys
import os
import argparse
import subprocess
import json
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_SHARED_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..', 'shared'))
for _p in (_THIS_DIR, _SHARED_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from RRsp3 import RRsp3, build_neighs_bk_from_bonds, make_bk_slots_clustered, make_exclusions_1st_2nd
from RRsp3_utils import (
    build_packed_system,
    build_packed_system_from_xyz,
    make_ports_from_neighs,
    gather_dx_port,
    compute_momentum_deltas,
    compute_mean_constraint_error,
    parse_kernel_debug_logs,
    validate_interactions,
    build_interaction_truth,
    check_local_ranges,
    write_xyz_frame,
    write_test_report,
    print_report,
    make_planar_fixmask,
    perturb_state_planar,
)

# ------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------

ALL_KERNELS = ['current', 'orig', 'substep_optimized', 'shapematch', 'eigen']

SYSTEM_SPECS = {
    'L0': {'names': ['h2o'], 'shift': None},
    'L1': {'names': ['h2o', 'h2o'], 'shift': None},
    'L2': {'names': ['ch3oh'], 'shift': None},
    'L3': {'names': ['ch3oh', 'ch3oh'], 'shift': None},
    'L4': {'names': ['ch2nh', 'h2o'], 'shift': None},
    'L5': {'names': ['tri3'], 'shift': None},
}

# ------------------------------------------------------------------
# Momentum Test
# ------------------------------------------------------------------

def cmd_momentum(args):
    kernels = _parse_kernels(args.kernel)
    results_all = []
    any_fail = False

    for kname in kernels:
        print(f"\n{'='*50}")
        print(f"MOMENTUM TEST  system={args.system}  kernel={kname}")
        print(f"{'='*50}")

        spec = SYSTEM_SPECS[args.system]
        sysdata = build_packed_system(spec['names'], shift=spec['shift'], group_size=64)
        natoms = sysdata['natoms']

        sim = RRsp3(natoms, group_size=64, prefer_gpu=True)

        # Perturb state
        rng = np.random.default_rng(int(args.seed))
        if args.plane:
            pos_p, quat_p = perturb_state_planar(sysdata['pos'], sysdata['quat'], pos_scale=float(args.pos_scale), rot_scale=float(args.rot_scale), rng=rng, plane=args.plane)
            fixmask = make_planar_fixmask(natoms, plane=args.plane)
            fixmask[sysdata['is_pad']] |= (1 | 2 | 4)
        else:
            from XPTB_utils import perturb_state
            pos_p, quat_p = perturb_state(sysdata['pos'], sysdata['quat'], pos_scale=float(args.pos_scale), rot_scale=float(args.rot_scale), rng=rng)
            fixmask = sysdata['fixmask'].copy()

        sim.upload_state(pos_p, sysdata['invm'], quat=quat_p)
        sim.upload_fixmask(fixmask)
        sim.upload_radius(sysdata['rad'])
        sim.upload_neighs_and_exclusions(sysdata['neighs'], sysdata['excl1'], sysdata['excl2'])
        sim.upload_cluster_ports(sysdata['port_local_atoms'], sysdata['K_atoms'], nnode_per_group=sysdata['nnode_per_group'])
        sim.upload_bkSlots(sysdata['bkSlots'])
        sim.reset_momentum()

        dP_list = []
        dL_list = []
        traj_path = None
        fxyz = None
        if args.traj:
            traj_path = args.traj.replace('.xyz', f'_{args.system}_{kname}.xyz')
            fxyz = open(traj_path, 'w')

        for istep in range(int(args.steps)):
            pos4, quat4 = sim.download_pos_quat()
            pos0 = pos4[:, :3].copy()

            sim.step_cluster(
                nnode_per_group=sysdata['nnode_per_group'],
                dt=float(args.dt),
                k_coll=float(args.k_coll),
                relaxation=float(args.relax),
                bbox_margin=float(args.bbox_margin),
                momentum_beta=float(args.momentum_beta),
                port_kernel=kname,
                rot_mass_scale=float(args.rot_mass_scale),
                n_rot_substeps=int(args.n_rot_substeps),
                rot_eps=float(args.rot_eps),
                theta_max=float(args.theta_max),
            )

            # Download deltas for momentum check
            dpos_coll = sim.download_dpos_coll()
            dpos_node = sim.download_dpos_node(nnode_per_group=sysdata['nnode_per_group'])
            drot_node = sim.download_drot_node(nnode_per_group=sysdata['nnode_per_group'])
            dpos_neigh = sim.download_dpos_neigh(nnode_per_group=sysdata['nnode_per_group'])

            dx_port = gather_dx_port(natoms, sysdata['nnode_per_group'], 64, dpos_node, dpos_neigh, sysdata['bkSlots'])
            dx_coll = dpos_coll[:, :3]

            dP, dL = compute_momentum_deltas(
                pos0, sysdata['m'], dx_coll, dx_port, dpos_node, drot_node,
                sysdata['nnode_per_group'], 64, sysdata['bkSlots'],
                relax=float(args.relax), massless_rot=_is_massless_kernel(kname)
            )
            nP = float(np.linalg.norm(dP))
            nL = float(np.linalg.norm(dL))
            dP_list.append(nP)
            dL_list.append(nL)

            if args.verbose:
                print(f"  step {istep:4d} |dP|={nP:.3e} |dL|={nL:.3e}")

            if nP > float(args.epsP) or nL > float(args.epsL):
                print(f"FAIL: momentum conservation at step={istep}: |dP|={nP:.3e} |dL|={nL:.3e}")
                any_fail = True
                break

            if fxyz is not None and (istep % int(args.traj_every)) == 0:
                real = sysdata['real']
                elems_real = [e for e in sysdata['elems_all'] if e != 'X']
                comment = f"step={istep} |dP|={nP:.3e} |dL|={nL:.3e}"
                write_xyz_frame(fxyz, elems_real, pos0[real], comment=comment)

        if fxyz is not None:
            fxyz.close()

        if len(dP_list) == int(args.steps):
            print(f"[PASS] All {args.steps} steps within tolerance  |dP|<{args.epsP}  |dL|<{args.epsL}")

        results_all.append({
            'test_type': 'momentum',
            'system': args.system,
            'kernel': kname,
            'steps': len(dP_list),
            'dP_max': float(max(dP_list)) if dP_list else None,
            'dL_max': float(max(dL_list)) if dL_list else None,
            'dP_final': float(dP_list[-1]) if dP_list else None,
            'dL_final': float(dL_list[-1]) if dL_list else None,
            'pass': len(dP_list) == int(args.steps),
            'trajectory': traj_path,
        })

    if args.report:
        write_test_report(args.report, {'momentum_tests': results_all})

    return 1 if any_fail else 0


# ------------------------------------------------------------------
# Topology Test
# ------------------------------------------------------------------

def cmd_topology(args):
    spec = SYSTEM_SPECS[args.system]
    shift_val = args.shift if args.shift is not None else 1.5

    # Build shift list from single value (offset per molecule)
    nmol = len(spec['names'])
    shift = [np.array([float(i) * shift_val, 0.0, 0.0], dtype=np.float32) for i in range(nmol)]

    sysdata = build_packed_system(spec['names'], shift=shift, group_size=64)
    natoms = sysdata['natoms']

    # Build truth table
    truth = build_interaction_truth(spec['names'], sysdata['packed'], shift, collision_radius=1.0)

    # Determine debug window from system size
    nreal = int(np.sum(sysdata['real']))
    dbg_end = min(nreal, 8)
    dbg_gids = set(range(dbg_end))

    # Compile with debug prints
    build_opts = [
        "-DENABLE_DEBUG_PRINTS",
        f"-DDEBUG_VERBOSITY=3",
        f"-DDEBUG_COMPONENTS=1",      # collision only
        f"-DDEBUG_TARGET_WG=-1",
        f"-DDEBUG_GID_START=0",
        f"-DDEBUG_GID_END={dbg_end}",
        "-DENABLE_COLL=1",
        "-DENABLE_PORT=0",
    ]

    sim = RRsp3(natoms, group_size=64, prefer_gpu=True, build_options=build_opts)
    sim.upload_state(sysdata['pos'], sysdata['invm'], quat=sysdata['quat'])
    sim.upload_fixmask(sysdata['fixmask'])
    sim.upload_radius(sysdata['rad'])
    sim.upload_neighs_and_exclusions(sysdata['neighs'], sysdata['excl1'], sysdata['excl2'])

    # Run one step via debug runner subprocess to capture stdout
    tmpdir = os.path.join(_THIS_DIR, "_tmp_rrsp3_headless")
    os.makedirs(tmpdir, exist_ok=True)

    np.save(os.path.join(tmpdir, "pos.npy"), sysdata['pos'])
    np.save(os.path.join(tmpdir, "invm.npy"), sysdata['invm'])
    np.save(os.path.join(tmpdir, "quat.npy"), sysdata['quat'])
    np.save(os.path.join(tmpdir, "rad.npy"), sysdata['rad'])
    np.save(os.path.join(tmpdir, "neighs.npy"), sysdata['neighs'])
    np.save(os.path.join(tmpdir, "excl1.npy"), sysdata['excl1'])
    np.save(os.path.join(tmpdir, "excl2.npy"), sysdata['excl2'])

    env = os.environ.copy()
    env['RRSP3_TMPDIR'] = tmpdir
    env['RRSP3_NNODE_PER_GROUP'] = str(sysdata['nnode_per_group'])
    env['RRSP3_BUILD_OPTS'] = " ".join(build_opts)

    runner = os.path.join(_THIS_DIR, "test_RRsp3_debug_runner.py")
    if not os.path.exists(runner):
        print(f"WARNING: debug runner not found at {runner}; using inline execution (stdout may not capture kernel prints)")
        # Inline fallback
        sim.step_cluster(
            nnode_per_group=sysdata['nnode_per_group'],
            dt=0.1, k_coll=50.0, relaxation=0.5, bbox_margin=0.5,
            port_kernel='current'
        )
        print("[SKIP] Cannot validate topology without debug runner subprocess for stdout capture")
        return 0

    cmd = [sys.executable, runner]
    out = subprocess.run(cmd, cwd=os.path.abspath(os.path.join(_THIS_DIR, "..", "..")), capture_output=True, text=True, env=env)

    parsed = parse_kernel_debug_logs(out.stdout)
    passed, report = validate_interactions(parsed, truth, group_size=64, dbg_gids=dbg_gids)

    print("\n--- Topology Validation Report ---")
    for line in report:
        print(line)
    print(f"\nOverall: {'PASS' if passed else 'FAIL'}")

    if args.report:
        write_test_report(args.report, {
            'test_type': 'topology',
            'system': args.system,
            'pass': passed,
            'report': report,
            'n_coll_logs': len(parsed.get('coll', [])),
            'n_topo_logs': len(parsed.get('topology', [])),
        })

    return 0 if passed else 1


# ------------------------------------------------------------------
# Smoke Test
# ------------------------------------------------------------------

def cmd_smoke(args):
    spec = SYSTEM_SPECS[args.system]
    sysdata = build_packed_system(spec['names'], shift=spec['shift'], group_size=64)
    natoms = sysdata['natoms']

    sim = RRsp3(natoms, group_size=64, prefer_gpu=True)
    sim.upload_state(sysdata['pos'], sysdata['invm'], quat=sysdata['quat'])
    sim.upload_fixmask(sysdata['fixmask'])
    sim.upload_radius(sysdata['rad'])
    sim.upload_neighs_and_exclusions(sysdata['neighs'], sysdata['excl1'], sysdata['excl2'])
    sim.upload_cluster_ports(sysdata['port_local_atoms'], sysdata['K_atoms'], nnode_per_group=sysdata['nnode_per_group'])
    sim.upload_bkSlots(sysdata['bkSlots'])

    # Run one step
    sim.step_cluster(
        nnode_per_group=sysdata['nnode_per_group'],
        dt=0.1, k_coll=50.0, relaxation=0.5, bbox_margin=0.5,
        port_kernel='current'
    )

    # Check outputs (only real atoms, padding has NaN by design)
    pos4, quat4 = sim.download_pos_quat()
    if np.any(np.isnan(pos4[sysdata['real']])):
        raise RuntimeError("SMOKE FAIL: NaN in positions after step")
    if np.any(np.isnan(quat4[sysdata['real']])):
        raise RuntimeError("SMOKE FAIL: NaN in quaternions after step")

    ghost_counts = sim.download(sim.cl_ghost_counts, (natoms // 64,), np.int32)
    ghost_indices = sim.download(sim.cl_ghost_indices, ((natoms // 64) * sim.max_ghosts,), np.int32)
    neighs_local = sim.download(sim.cl_neighs_local, (natoms, 4), np.int32)
    excl1_local = sim.download(sim.cl_excl1_local, (natoms, 4), np.int32)
    excl2_local = sim.download(sim.cl_excl2_local, (natoms, 4), np.int32)

    check_local_ranges(neighs_local, excl1_local, excl2_local, ghost_counts, group_size=64)

    print("ghost_counts:", ghost_counts.tolist())
    for ig in range(natoms // 64):
        g = int(ghost_counts[ig])
        print(f"  group {ig}: {g} ghosts (first few: {ghost_indices[ig * sim.max_ghosts:ig * sim.max_ghosts + min(g, 8)].tolist()})")

    print("[PASS] Smoke test completed without NaN or range errors")
    return 0


# ------------------------------------------------------------------
# Collision / Reindexing Validation
# ------------------------------------------------------------------

def cmd_collision(args):
    """Run GPU collision pipeline and compare against vectorized brute-force reference.

    Builds a multi-molecule system with forced overlaps, runs broad+narrow phase,
    downloads dpos_coll, ghosts, and local indices, then:
      1. Compares GPU dpos_coll to brute_force_collision_dpos()
      2. Validates ghost list (no duplicates, within AABB margin)
      3. Validates neighs/excl round-trip (global -> local -> global)
      4. Optionally plots atoms colored by workgroup
    """
    from RRsp3_utils import (
        brute_force_collision_dpos,
        validate_ghost_list,
        validate_neighs_local,
        plot_workgroups,
    )

    # Build overlapping system: 4 H2O molecules close together
    # With group_size=8 -> 2 workgroups (8+8), forcing ghost-atom exchange.
    group_size = int(args.group_size)
    sysdata = build_packed_system(
        ['h2o'] * 4,
        shift=[[i * 1.2, 0.0, 0.0] for i in range(4)],
        group_size=group_size,
    )
    natoms = sysdata['natoms']
    sim = RRsp3(natoms, group_size=group_size, prefer_gpu=True)

    # Perturb positions to create overlaps (radius=1.0 for all real atoms)
    pos0 = sysdata['pos'].copy()
    rng = np.random.default_rng(int(args.seed))
    pos0 += rng.normal(scale=0.3, size=pos0.shape).astype(np.float32)

    sim.upload_state(pos0, sysdata['invm'], quat=sysdata['quat'])
    sim.upload_fixmask(sysdata['fixmask'])
    sim.upload_radius(sysdata['rad'])
    sim.upload_neighs_and_exclusions(sysdata['neighs'], sysdata['excl1'], sysdata['excl2'])
    sim.reset_momentum()

    # Run broad phase + topology to populate ghosts and local indices
    sim.run_bboxes_and_topology(bbox_margin=args.bbox_margin)

    # Download topology data for validation
    ghost_indices, ghost_counts = sim.download_ghosts()
    neighs_local = sim.download_neighs_local()
    excl1_local, excl2_local = sim.download_excl_local()
    bmin, bmax = sim.download_bboxes()

    # Validate ghosts
    rmax = float(np.max(sysdata['rad'][sysdata['real']])) if np.any(sysdata['real']) else 0.0
    margin_sq = float((2.0 * rmax + float(args.bbox_margin)) ** 2)
    ok_ghost, report_ghost = validate_ghost_list(
        pos0, sysdata['rad'], ghost_indices.ravel(), ghost_counts, group_size, margin_sq
    )
    for line in report_ghost:
        print(line)

    # Validate reindexing
    ok_reidx, report_reidx = validate_neighs_local(
        sysdata['neighs'], neighs_local, ghost_indices.ravel(), ghost_counts, group_size
    )
    for line in report_reidx:
        print(line)

    # Also validate exclusions (using same round-trip checker)
    ok_excl, report_excl = validate_neighs_local(
        sysdata['excl1'], excl1_local, ghost_indices.ravel(), ghost_counts, group_size
    )
    for line in report_excl:
        print(line)
    ok_excl2, report_excl2 = validate_neighs_local(
        sysdata['excl2'], excl2_local, ghost_indices.ravel(), ghost_counts, group_size
    )
    for line in report_excl2:
        print(line)

    # Run collision kernel (skip_ports avoids needing port buffers)
    dpos_coll_gpu4, _, _, _ = sim.compute_cluster_deltas(
        nnode_per_group=sysdata['nnode_per_group'],
        dt=args.dt,
        k_coll=args.k_coll,
        bbox_margin=args.bbox_margin,
        skip_ports=True,
    )
    dpos_coll_gpu = dpos_coll_gpu4[:, :3]

    # Brute-force reference from ORIGINAL positions (GPU collision reads same pos)
    dpos_coll_ref = brute_force_collision_dpos(
        pos0, sysdata['rad'], sysdata['excl1'], sysdata['excl2'], invm=sysdata['invm']
    )

    # Compare only real atoms
    real = sysdata['real']
    diff = np.abs(dpos_coll_gpu[real] - dpos_coll_ref[real])
    max_err = float(np.max(diff))
    mean_err = float(np.mean(diff))
    tol = float(args.tol)

    print(f"Collision comparison: max_err={max_err:.6e}  mean_err={mean_err:.6e}  tol={tol}")
    if max_err > tol:
        print(f"[FAIL] GPU collision deviates from brute-force reference")
        any_fail = True
    else:
        print(f"[PASS] GPU collision matches brute-force reference")
        any_fail = False

    if not (ok_ghost and ok_reidx and ok_excl and ok_excl2):
        any_fail = True
        print("[FAIL] Reindexing validation failed")
    else:
        print("[PASS] Reindexing validation passed")

    if args.plot:
        plot_workgroups(pos0, sysdata['packed']['group_id'], sysdata['is_pad'], outpath=args.plot)

    return 1 if any_fail else 0


# ------------------------------------------------------------------
# Convergence Test
# ------------------------------------------------------------------

def cmd_convergence(args):
    kernels = _parse_kernels(args.kernel)
    betas = [float(v) for v in args.beta_list.split(',')] if args.beta_list else [float(args.momentum_beta)]
    dts   = [float(v) for v in args.dt_list.split(',')] if args.dt_list else [float(args.dt)]
    all_results = []
    any_fail = False

    for kname in kernels:
        for beta in betas:
            for dt in dts:
                print(f"\n{'='*50}")
                print(f"CONVERGENCE TEST  xyz={args.xyz}  kernel={kname}  beta={beta}  dt={dt}")
                print(f"{'='*50}")

                sysdata = build_packed_system_from_xyz(args.xyz, group_size=64)
                natoms = sysdata['natoms']

                sim = RRsp3(natoms, group_size=64, prefer_gpu=True)

                # Compute ports from ORIGINAL (unstretched) geometry so the solver
                # tries to relax back to the original bond lengths.
                port_local_orig, K_atoms = make_ports_from_neighs(sysdata['pos'], sysdata['neighs'], K=200.0)

                # Stretch / perturb initial state
                pos0 = sysdata['pos'].copy()
                quat0 = sysdata['quat'].copy()
                if args.stretch != 1.0:
                    # Scale x,y in-plane (planar sheet)
                    pos0[:, 0] *= args.stretch
                    pos0[:, 1] *= args.stretch
                if args.pos_scale > 0.0:
                    rng = np.random.default_rng(int(args.seed))
                    noise = rng.normal(scale=args.pos_scale, size=pos0.shape).astype(np.float32)
                    pos0 += noise

                sim.upload_state(pos0, sysdata['invm'], quat=quat0)
                sim.upload_fixmask(sysdata['fixmask'])
                sim.upload_radius(sysdata['rad'])
                sim.upload_neighs_and_exclusions(sysdata['neighs'], sysdata['excl1'], sysdata['excl2'])
                sim.upload_cluster_ports(port_local_orig, K_atoms, nnode_per_group=sysdata['nnode_per_group'])
                sim.upload_bkSlots(sysdata['bkSlots'])
                sim.reset_momentum()

                errors = []
                geo_errs = []
                traj_path = None
                fxyz = None
                if args.traj:
                    tag = f"{kname}_b{beta}_d{dt}"
                    traj_path = args.traj.replace('.xyz', f'_{tag}.xyz')
                    fxyz = open(traj_path, 'w')

                for istep in range(int(args.max_iter)):
                    sim.step_cluster(
                        nnode_per_group=sysdata['nnode_per_group'],
                        dt=dt,
                        k_coll=float(args.k_coll),
                        relaxation=float(args.relax),
                        bbox_margin=float(args.bbox_margin),
                        momentum_beta=beta,
                        port_kernel=kname,
                        rot_mass_scale=float(args.rot_mass_scale),
                        n_rot_substeps=int(args.n_rot_substeps),
                        rot_eps=float(args.rot_eps),
                        theta_max=float(args.theta_max),
                    )

                    dpos_node = sim.download_dpos_node(nnode_per_group=sysdata['nnode_per_group'])
                    dpos_neigh = sim.download_dpos_neigh(nnode_per_group=sysdata['nnode_per_group'])
                    step_err = float(np.mean(np.abs(dpos_node[:, :3])) + np.mean(np.abs(dpos_neigh[:, :3])))
                    errors.append(step_err)

                    pos4, quat4 = sim.download_pos_quat()
                    gerr = compute_mean_constraint_error(pos4, quat4, sysdata['neighs'], port_local_orig)
                    geo_errs.append(gerr)

                    if args.verbose:
                        print(f"  step {istep:4d}  mean_step={step_err:.6e}  geo_err={gerr:.6e}")

                    if fxyz is not None and (istep % int(args.traj_every)) == 0:
                        real = sysdata['real']
                        elems_real = [e for e in sysdata['elems_all'] if e != 'X']
                        comment = f"step={istep} mean_step={step_err:.6e} geo_err={gerr:.6e}"
                        write_xyz_frame(fxyz, elems_real, pos4[real, :3], comment=comment)

                    if step_err < float(args.tol):
                        print(f"[PASS] Converged at step {istep}: mean_step={step_err:.3e} geo_err={gerr:.3e} < tol={args.tol}")
                        break
                else:
                    print(f"[FAIL] Did not converge within {args.max_iter} iterations: mean_step={errors[-1]:.3e} geo_err={geo_errs[-1]:.3e}")
                    any_fail = True

                if fxyz is not None:
                    fxyz.close()

                label = kname
                if len(betas) > 1:
                    label += f' β={beta}'
                if len(dts) > 1:
                    label += f' dt={dt}'

                all_results.append({
                    'test_type': 'convergence',
                    'xyz': args.xyz,
                    'kernel': kname,
                    'label': label,
                    'stretch': args.stretch,
                    'beta': beta,
                    'dt': dt,
                    'iterations': len(errors),
                    'final_error': float(errors[-1]),
                    'converged': errors[-1] < float(args.tol),
                    'trajectory': traj_path,
                    'step_errs': errors[:],
                    'geo_errs': geo_errs[:],
                })

    if args.save_npz:
        data = {}
        for r in all_results:
            k = r['kernel']
            data[f'{k}_step'] = np.array(r.get('step_errs', []), dtype=np.float32)
            data[f'{k}_geo'] = np.array(r.get('geo_errs', []), dtype=np.float32)
        np.savez(args.save_npz, **data)
        print(f"Saved convergence data to {args.save_npz}")

    if args.plot:
        _plot_convergence(all_results, args)

    if args.report:
        write_test_report(args.report, {'convergence_tests': all_results})

    return 1 if any_fail else 0


# ------------------------------------------------------------------
# Plot Convergence
# ------------------------------------------------------------------

def _plot_convergence(all_results, args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    betas = [float(v) for v in args.beta_list.split(',')] if args.beta_list else [float(args.momentum_beta)]
    dts   = [float(v) for v in args.dt_list.split(',')] if args.dt_list else [float(args.dt)]
    btxt = f"beta=[{','.join(str(b) for b in betas)}]" if len(betas) > 1 else f"beta={betas[0]}"
    dtxt = f"dt=[{','.join(str(d) for d in dts)}]" if len(dts) > 1 else f"dt={dts[0]}"
    title = f"kernel={args.kernel}  stretch={args.stretch}  {dtxt}  {btxt}  relax={args.relax}"
    fig.suptitle(title, fontsize=10)

    for r in all_results:
        label = r.get('label', r['kernel'])
        steps = np.arange(r['iterations'])
        ax1.plot(steps, r['step_errs'], label=label, linewidth=1.5)
        ax2.plot(steps, r['geo_errs'], label=label, linewidth=1.5)

    ax1.set_yscale('log')
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Mean step size')
    ax1.set_title('Solver step size')
    ax1.legend()
    ax1.grid(True, which='both', ls='--', alpha=0.4)

    ax2.set_yscale('log')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Mean geometric error')
    ax2.set_title('Tip-to-atom constraint violation')
    ax2.legend()
    ax2.grid(True, which='both', ls='--', alpha=0.4)

    fig.tight_layout()
    out = args.plot if isinstance(args.plot, str) and args.plot else 'convergence.png'
    fig.savefig(out, dpi=150)
    print(f"Saved plot to {out}")


# ------------------------------------------------------------------
# Suite Test
# ------------------------------------------------------------------

def cmd_suite(args):
    levels = args.level.split(',') if args.level else ['L0', 'L1', 'L2']
    kernels = _parse_kernels(args.kernel)
    all_results = []
    any_fail = False

    for level in levels:
        for kname in kernels:
            print(f"\n{'='*50}")
            print(f"SUITE: {level} / {kname}")
            print(f"{'='*50}")

            # Reconstruct args-like dict for momentum test
            sub_args = argparse.Namespace(
                system=level,
                kernel=kname,
                steps=args.steps,
                dt=args.dt,
                k_coll=args.k_coll,
                relax=args.relax,
                bbox_margin=args.bbox_margin,
                momentum_beta=args.momentum_beta,
                rot_mass_scale=args.rot_mass_scale,
                n_rot_substeps=args.n_rot_substeps,
                rot_eps=args.rot_eps,
                theta_max=args.theta_max,
                pos_scale=args.pos_scale,
                rot_scale=args.rot_scale,
                seed=args.seed,
                epsP=args.epsP,
                epsL=args.epsL,
                plane=args.plane,
                traj='',
                traj_every=1,
                verbose=args.verbose,
                report='',
            )
            try:
                rc = cmd_momentum(sub_args)
            except Exception as e:
                print(f"EXCEPTION: {e}")
                rc = 1

            if rc != 0:
                any_fail = True

    print(f"\n{'='*50}")
    print(f"SUITE FINAL: {'FAIL' if any_fail else 'PASS'}")
    print(f"{'='*50}")
    return 1 if any_fail else 0


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def _parse_kernels(s):
    if s == 'all':
        return ALL_KERNELS[:]
    return [k.strip() for k in str(s).split(',') if k.strip()]


def _is_massless_kernel(kname):
    """Return True for massless geometric-alignment kernels (no physical rotational inertia)."""
    return kname in ('substep_optimized', 'shapematch', 'eigen')


# ------------------------------------------------------------------
# Main / CLI
# ------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="RRsp3 consolidated headless test")
    sub = ap.add_subparsers(dest='command', required=True)

    # Common args
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--system', type=str, default='L0', help='Test system: L0,L1,L2,L3,L4')
    common.add_argument('--report', type=str, default='', help='Write JSON report to path')
    common.add_argument('--verbose', action='store_true', help='Verbose per-step output')

    # momentum
    p_mom = sub.add_parser('momentum', parents=[common])
    p_mom.add_argument('--kernel', type=str, default='current', help='Kernel name or "all"')
    p_mom.add_argument('--steps', type=int, default=50)
    p_mom.add_argument('--dt', type=float, default=0.1)
    p_mom.add_argument('--k_coll', type=float, default=50.0)
    p_mom.add_argument('--relax', type=float, default=0.5)
    p_mom.add_argument('--bbox_margin', type=float, default=0.5)
    p_mom.add_argument('--momentum_beta', type=float, default=0.0)
    p_mom.add_argument('--rot_mass_scale', type=float, default=1.0)
    p_mom.add_argument('--n_rot_substeps', type=int, default=5)
    p_mom.add_argument('--rot_eps', type=float, default=0.0)
    p_mom.add_argument('--theta_max', type=float, default=0.0)
    p_mom.add_argument('--pos_scale', type=float, default=0.05)
    p_mom.add_argument('--rot_scale', type=float, default=0.0)
    p_mom.add_argument('--seed', type=int, default=0)
    p_mom.add_argument('--epsP', type=float, default=1e-4)
    p_mom.add_argument('--epsL', type=float, default=1e-3)
    p_mom.add_argument('--plane', type=str, default='', help="Constrain to plane: xy, xz, yz")
    p_mom.add_argument('--traj', type=str, default='', help='Write trajectory to .xyz')
    p_mom.add_argument('--traj_every', type=int, default=1)
    p_mom.set_defaults(func=cmd_momentum)

    # topology
    p_topo = sub.add_parser('topology', parents=[common])
    p_topo.add_argument('--shift', type=float, default=None, help='Inter-molecular spacing')
    p_topo.set_defaults(func=cmd_topology)

    # smoke
    p_smoke = sub.add_parser('smoke', parents=[common])
    p_smoke.set_defaults(func=cmd_smoke)

    # collision / reindexing validation
    p_coll = sub.add_parser('collision', parents=[common])
    p_coll.add_argument('--group_size', type=int, default=8, help='Workgroup size (must divide natoms)')
    p_coll.add_argument('--seed', type=int, default=0)
    p_coll.add_argument('--dt', type=float, default=0.1)
    p_coll.add_argument('--k_coll', type=float, default=50.0)
    p_coll.add_argument('--bbox_margin', type=float, default=0.5)
    p_coll.add_argument('--tol', type=float, default=1e-4, help='Max error tolerance vs brute-force')
    p_coll.add_argument('--plot', type=str, default='', help='Path to save workgroup plot PNG')
    p_coll.set_defaults(func=cmd_collision)

    # convergence
    p_conv = sub.add_parser('convergence', parents=[common])
    p_conv.add_argument('--xyz', type=str, default='/home/prokophapala/git/FireCore/cpp/common_resources/xyz/pentacene.xyz', help='Input XYZ file')
    p_conv.add_argument('--kernel', type=str, default='current', help='Kernel name or "all"')
    p_conv.add_argument('--max_iter', type=int, default=500)
    p_conv.add_argument('--dt', type=float, default=0.1)
    p_conv.add_argument('--dt_list', type=str, default='', help='Comma-separated dt values for param sweep')
    p_conv.add_argument('--k_coll', type=float, default=50.0)
    p_conv.add_argument('--relax', type=float, default=0.5)
    p_conv.add_argument('--bbox_margin', type=float, default=0.5)
    p_conv.add_argument('--momentum_beta', type=float, default=0.0)
    p_conv.add_argument('--beta_list', type=str, default='', help='Comma-separated beta values for param sweep')
    p_conv.add_argument('--rot_mass_scale', type=float, default=1.0)
    p_conv.add_argument('--n_rot_substeps', type=int, default=5)
    p_conv.add_argument('--rot_eps', type=float, default=0.0)
    p_conv.add_argument('--theta_max', type=float, default=0.0)
    p_conv.add_argument('--pos_scale', type=float, default=0.0)
    p_conv.add_argument('--stretch', type=float, default=2.0, help='In-plane stretch factor (1.0 = no stretch)')
    p_conv.add_argument('--seed', type=int, default=0)
    p_conv.add_argument('--tol', type=float, default=1e-5, help='Mean step size tolerance')
    p_conv.add_argument('--traj', type=str, default='', help='Write trajectory to .xyz')
    p_conv.add_argument('--traj_every', type=int, default=1)
    p_conv.add_argument('--save_npz', type=str, default='', help='Save step/geo error arrays to .npz')
    p_conv.add_argument('--plot', type=str, nargs='?', const='convergence.png', default='', help='Plot convergence curves (optional path)')
    p_conv.set_defaults(func=cmd_convergence)

    # suite
    p_suite = sub.add_parser('suite', parents=[common])
    p_suite.add_argument('--level', type=str, default='L0,L1,L2', help='Comma-separated levels')
    p_suite.add_argument('--kernel', type=str, default='current', help='Kernel name or "all"')
    p_suite.add_argument('--steps', type=int, default=50)
    p_suite.add_argument('--dt', type=float, default=0.1)
    p_suite.add_argument('--k_coll', type=float, default=50.0)
    p_suite.add_argument('--relax', type=float, default=0.5)
    p_suite.add_argument('--bbox_margin', type=float, default=0.5)
    p_suite.add_argument('--momentum_beta', type=float, default=0.0)
    p_suite.add_argument('--rot_mass_scale', type=float, default=1.0)
    p_suite.add_argument('--n_rot_substeps', type=int, default=5)
    p_suite.add_argument('--rot_eps', type=float, default=0.0)
    p_suite.add_argument('--theta_max', type=float, default=0.0)
    p_suite.add_argument('--pos_scale', type=float, default=0.05)
    p_suite.add_argument('--rot_scale', type=float, default=0.0)
    p_suite.add_argument('--seed', type=int, default=0)
    p_suite.add_argument('--epsP', type=float, default=1e-4)
    p_suite.add_argument('--epsL', type=float, default=1e-3)
    p_suite.add_argument('--plane', type=str, default='')
    p_suite.set_defaults(func=cmd_suite)

    args = ap.parse_args()
    rc = args.func(args)
    sys.exit(rc)


if __name__ == '__main__':
    main()
