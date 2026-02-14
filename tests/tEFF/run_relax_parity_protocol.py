#!/usr/bin/env python3
import os, re, argparse, numpy as np
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

from pyBall import eFF as eff
from pyBall.OCL import eFF_ocl


def parse_core_mode(xyz_path):
    with open(xyz_path, 'r') as f:
        f.readline()
        c = f.readline().strip()
    m = re.search(r"na,ne,core\s+\d+\s+\d+\s+(\w)", c)
    if not m:
        raise RuntimeError(f"Cannot parse core mode from {xyz_path}")
    return m.group(1)


def cpu_relax_single(xyz_in, xyz_out, trj_out, nsteps, dt, fconv, opt_alg):
    eff.setVerbosity(0, 0)
    eff.setTrjName(trj_out, savePerNsteps=1, bDel=True)
    Es, _, _ = eff.processXYZ_e(xyz_in, nstepMax=nsteps, dt=dt, Fconv=fconv, optAlg=opt_alg, xyz_out=xyz_out, fgo_out=None, bOutputs=(1, 0, 0))
    eff.getBuffs()
    na, ne = eff.na, eff.ne
    pos = np.zeros((na + ne, 4), dtype=np.float64)
    pos[:na, :3] = eff.apos
    pos[na:, :3] = eff.epos
    pos[na:, 3] = eff.esize
    frc = np.zeros((na + ne, 4), dtype=np.float64)
    frc[:na, :3] = eff.aforce
    frc[na:, :3] = eff.eforce
    frc[na:, 3] = eff.fsize
    E5 = np.array(Es, dtype=np.float64)
    if E5.ndim == 2:
        E5 = E5[-1]
    return dict(na=na, ne=ne, pos=pos, frc=frc, E5=E5[:5])


def _write_frame(fout, sysdef, pblock):
    na = sysdef['na']
    ne = sysdef['ne']
    core_mode = sysdef.get('core_mode', 'f')
    fout.write(f"{na+ne}\n")
    fout.write(f"na,ne,core {na} {ne} {core_mode} | GPU localMD relaxed\n")
    for ia in range(na):
        name = sysdef['ions'][ia][0]
        x, y, z = pblock[ia, :3]
        fout.write(f"{name} {x:.8f} {y:.8f} {z:.8f}\n")
    for ie in range(ne):
        spin = sysdef['electrons'][ie][1]
        x, y, z, s = pblock[na + ie]
        if spin == 0:
            nm, q = 'e2', 2.0
        elif spin > 0:
            nm, q = 'e+', 1.0
        else:
            nm, q = 'e-', -1.0
        fout.write(f"{nm} {x:.8f} {y:.8f} {z:.8f} {q:.1f} {s:.8f}\n")


def gpu_relax(xyz_in, xyz_out, trj_out, nsteps, dt, damping, nloc, device, fix_atoms=False):
    ocl = eFF_ocl.EFF_OCL(nloc=nloc, device_index=device, bEnergyKernel=True)
    ocl.load_xyzs(xyz_in, bPrint=False)
    if len(ocl.systems) == 0:
        raise RuntimeError(f"No systems parsed from {xyz_in}")

    # Optional: fix all nuclei (ions) using fixmask bits {x|y|z}=7; relax electrons only
    # (CPU uses constraints; GPU uses fixmask buffer)
    if fix_atoms:
        for sysdef in ocl.systems:
            na = int(sysdef['na']); ne = int(sysdef['ne'])
            fm = np.zeros((na + ne,), dtype=np.uint8)
            fm[:na] = 7
            sysdef['fixmask'] = fm
    ocl.realloc_buffers()
    ocl.upload_data()

    b_frozen = all(s.get('core_mode', 'f') == 'f' for s in ocl.systems)

    pos = None
    if trj_out and nsteps > 0:
        open(trj_out, 'w').close()
        for _ in range(nsteps):
            pos = ocl.relax_systems(n_steps=1, dt=dt, damping=damping, bFrozenCore=b_frozen)
            with open(trj_out, 'a') as fout:
                i0 = 0
                for sysdef in ocl.systems:
                    nt = sysdef['na'] + sysdef['ne']
                    _write_frame(fout, sysdef, pos[i0:i0+nt])
                    i0 += nt
    else:
        pos = ocl.relax_systems(n_steps=max(1, nsteps), dt=dt, damping=damping, bFrozenCore=b_frozen)

    forces = ocl.get_forces().copy()
    E5 = ocl.eval_energies_localmd(bFrozenCore=b_frozen)

    with open(xyz_out, 'w') as fout:
        i0 = 0
        for sysdef in ocl.systems:
            nt = sysdef['na'] + sysdef['ne']
            _write_frame(fout, sysdef, pos[i0:i0+nt])
            i0 += nt

    return dict(pos=pos.astype(np.float64), frc=forces.astype(np.float64), E5=E5.astype(np.float64), systems=ocl.systems)


def _write_xyz_like_input(xyz_template, xyz_out, pos):
    ocl = eFF_ocl.EFF_OCL(nloc=max(32, pos.shape[0]), device_index=0, bEnergyKernel=False)
    ocl.load_xyzs(xyz_template, bPrint=False)
    if len(ocl.systems) != 1:
        raise RuntimeError(f"_write_xyz_like_input expects 1 system in template {xyz_template}, got {len(ocl.systems)}")
    sysdef = ocl.systems[0]
    nt = int(sysdef['na']) + int(sysdef['ne'])
    if pos.shape[0] < nt:
        raise RuntimeError(f"_write_xyz_like_input: pos has {pos.shape[0]} rows but nt={nt}")
    with open(xyz_out, 'w') as fout:
        _write_frame(fout, sysdef, pos[:nt])


def _systems_to_fixed_arrays(systems):
    if len(systems) == 0:
        raise RuntimeError("No systems")
    nconf = len(systems)
    na0 = int(systems[0]['na']); ne0 = int(systems[0]['ne'])
    for i, s in enumerate(systems):
        if int(s['na']) != na0 or int(s['ne']) != ne0:
            raise RuntimeError(f"[scan] inconsistent na/ne across confs: conf0 na={na0} ne={ne0} vs conf{i} na={s['na']} ne={s['ne']}")
    fixed_poss = np.zeros((nconf, na0, 4), dtype=np.float64)
    fixed_inds = np.zeros((nconf, na0, 2), dtype=np.int32)
    fixed_inds[:, :, 0] = np.arange(na0, dtype=np.int32)[None, :]
    fixed_inds[:, :, 1] = 7
    for ic, sysdef in enumerate(systems):
        for ia in range(na0):
            fixed_poss[ic, ia, :3] = sysdef['ions'][ia][1]
    return fixed_poss, fixed_inds, na0, ne0


def _systems_to_fixed_arrays(systems):
    if len(systems) == 0:
        raise RuntimeError("No systems")
    nconf = len(systems)
    na0 = int(systems[0]['na']); ne0 = int(systems[0]['ne'])
    for i, s in enumerate(systems):
        if int(s['na']) != na0 or int(s['ne']) != ne0:
            raise RuntimeError(f"[scan] inconsistent na/ne across confs: conf0 na={na0} ne={ne0} vs conf{i} na={s['na']} ne={s['ne']}")
    fixed_poss = np.zeros((nconf, na0, 4), dtype=np.float64)
    fixed_inds = np.zeros((nconf, na0, 2), dtype=np.int32)
    fixed_inds[:, :, 0] = np.arange(na0, dtype=np.int32)[None, :]
    fixed_inds[:, :, 1] = 7
    for ic, sysdef in enumerate(systems):
        for ia in range(na0):
            fixed_poss[ic, ia, :3] = sysdef['ions'][ia][1]
    return fixed_poss, fixed_inds, na0, ne0


def _maxabs_electron_force(frc, na, ne):
    if ne <= 0:
        return 0.0
    f = frc[na:na+ne, :3]
    return float(np.max(np.abs(f)))


def _assert_ions_fixed(pos0, pos1, na, tol=1e-8, label=''):
    if na <= 0:
        return
    d = pos1[:na, :3] - pos0[:na, :3]
    m = float(np.max(np.abs(d)))
    if m > tol:
        raise RuntimeError(f"[fixmask] ions moved! {label} max|dR|={m} tol={tol}")


def cpu_localmd_traj_cpp(xyz_in, trj_out, nsteps, dt, damping, fix_atoms=False, size_min=0.001):
    # Use C++ integrator (DynamicOpt::move_MD_dbg via ialg=3) to mirror GPU localMD, and let C++ write trajectory.
    # Note: run() currently resets velocities for ialg>0, so this is meant to match the current GPU path used in tests.
    eff.setVerbosity(0, 0)
    eff.setEFFDbgPair(0, 0, 0, 0, 1)

    # Initialize buffers and write initial frame into trajectory
    # We write the initial frame explicitly so that CPU/GPU trajectories are aligned.
    eff.processXYZ_e(xyz_in, nstepMax=0, dt=dt, Fconv=1e-3, optAlg=3, xyz_out=trj_out, fgo_out=None, bOutputs=(0, 0, 0))
    eff.getBuffs()
    na = int(eff.na); ne = int(eff.ne)

    if fix_atoms and na > 0:
        fixed_poss = np.zeros((na, 4), dtype=np.float64)
        fixed_inds = np.zeros((na, 2), dtype=np.int32)
        fixed_poss[:, :3] = eff.apos[:, :3]
        fixed_inds[:, 0] = np.arange(na, dtype=np.int32)
        fixed_inds[:, 1] = 7
        eff.set_constrains(na, fixed_poss, fixed_inds, True)

    # C++ trajectory (do not delete; initial frame already written)
    eff.setTrjName(trj_out, savePerNsteps=1, bDel=False)
    eff.initOpt(dt=dt, damping=damping, f_limit=1000.0, bMass=False)
    eff.run(nstepMax=int(nsteps), dt=dt, Fconv=0.0, ialg=3, bOutE=False, bOutF=False)

    eff.eval(); eff.getBuffs()
    E7 = eff.getEnergyTerms(sh=(7,)).astype(np.float64)
    Ek      = E7[0]
    Eee     = E7[1] + E7[2]
    Eae     = E7[4] + E7[5]
    Eaa     = E7[6]
    Etot    = Ek + Eee + Eae + Eaa
    E5 = np.array([Etot, Ek, Eee, Eae, Eaa], dtype=np.float64)

    pos = np.zeros((na + ne, 4), dtype=np.float64)
    if na > 0:
        pos[:na, :3] = eff.apos[:, :3]
    if ne > 0:
        pos[na:, :3] = eff.epos[:, :3]
        pos[na:, 3] = eff.esize[:]
        pos[na:, 3] = np.maximum(pos[na:, 3], float(size_min))
    return dict(na=na, ne=ne, pos=pos, E5=E5)


def cpu_md_mirror_traj(xyz_in, trj_out, nsteps, dt, damping, fix_atoms=False, size_min=0.001):
    """CPU-side integrator intended to mirror GPU localMD update:
    v = v*damping + f*dt
    p = p + v*dt
    with ion fix (if requested) and electron size clamp.

    Uses eff.eval() for forces at current positions.
    """
    eff.setVerbosity(0, 0)
    eff.setEFFDbgPair(0, 0, 0, 0, 1)
    eff.processXYZ_e(xyz_in, nstepMax=0, dt=dt, Fconv=1e-3, optAlg=1, xyz_out=None, fgo_out=None, bOutputs=(0, 0, 0))
    eff.getBuffs()
    na = int(eff.na); ne = int(eff.ne)
    nt = na + ne

    # local velocities (CPU-only); C++ buffers don't need them for force evaluation
    v = np.zeros((nt, 4), dtype=np.float64)

    # template for writing frames with correct ion names/spins/core_mode
    ocl_tmp = eFF_ocl.EFF_OCL(nloc=max(32, nt), device_index=0, bEnergyKernel=False)
    ocl_tmp.load_xyzs(xyz_in, bPrint=False)
    if len(ocl_tmp.systems) != 1:
        raise RuntimeError(f"cpu_md_mirror_traj: cannot parse template {xyz_in}")
    sysdef = ocl_tmp.systems[0]

    open(trj_out, 'w').close()

    apos_fix = None
    if fix_atoms and na > 0:
        apos_fix = np.array(eff.apos[:, :3], copy=True)

    # write initial frame
    pos0 = np.zeros((nt, 4), dtype=np.float64)
    pos0[:na, :3] = eff.apos
    pos0[na:, :3] = eff.epos
    pos0[na:, 3] = eff.esize
    with open(trj_out, 'a') as fout:
        _write_frame(fout, sysdef, pos0)

    for istep in range(int(nsteps)):
        eff.eval()
        eff.getBuffs()

        # assemble force in {x,y,z,w}
        f = np.zeros((nt, 4), dtype=np.float64)
        f[:na, :3] = eff.aforce
        f[na:, :3] = eff.eforce
        f[na:, 3] = eff.fsize

        # integrate
        v = v * damping + f * dt
        pos0 = pos0 + v * dt

        if fix_atoms and na > 0:
            pos0[:na, :3] = apos_fix  # keep ions fixed at initial positions
            v[:na, :] = 0.0

        if ne > 0:
            pos0[na:, 3] = np.maximum(pos0[na:, 3], float(size_min))

        # write back to C++ buffers for next force evaluation
        if na > 0:
            if fix_atoms:
                eff.apos[:, :3] = apos_fix
            else:
                eff.apos[:, :3] = pos0[:na, :3]
        if ne > 0:
            eff.epos[:, :3] = pos0[na:, :3]
            eff.esize[:] = pos0[na:, 3]

        with open(trj_out, 'a') as fout:
            _write_frame(fout, sysdef, pos0)

    # final snapshot
    eff.eval(); eff.getBuffs()
    E7 = eff.getEnergyTerms(sh=(7,)).astype(np.float64)
    Ek      = E7[0]
    Eee     = E7[1] + E7[2]
    Eae     = E7[4] + E7[5]
    Eaa     = E7[6]
    Etot    = Ek + Eee + Eae + Eaa
    E5 = np.array([Etot, Ek, Eee, Eae, Eaa], dtype=np.float64)
    return dict(na=na, ne=ne, pos=pos0.copy(), E5=E5)


def gpu_localmd_traj(xyz_in, trj_out, nsteps, dt, damping, nloc, device, fix_atoms=False):
    ocl = eFF_ocl.EFF_OCL(nloc=nloc, device_index=device, bEnergyKernel=True)
    ocl.load_xyzs(xyz_in, bPrint=False)
    if len(ocl.systems) != 1:
        raise RuntimeError(f"gpu_localmd_traj expects 1 system, got {len(ocl.systems)} from {xyz_in}")
    if fix_atoms:
        sysdef = ocl.systems[0]
        na = int(sysdef['na']); ne = int(sysdef['ne'])
        fm = np.zeros((na + ne,), dtype=np.uint8); fm[:na] = 7
        sysdef['fixmask'] = fm
    ocl.realloc_buffers(); ocl.upload_data()
    b_frozen = all(s.get('core_mode', 'f') == 'f' for s in ocl.systems)

    open(trj_out, 'w').close()
    na = int(ocl.systems[0]['na']); ne = int(ocl.systems[0]['ne']); nt = na + ne
    with open(trj_out, 'a') as fout:
        _write_frame(fout, ocl.systems[0], ocl.pos_h[:nt].astype(np.float64))
    for _ in range(int(nsteps)):
        pos = ocl.relax_systems(n_steps=1, dt=dt, damping=damping, bFrozenCore=b_frozen)
        ocl.pos_h[:, :] = pos
        with open(trj_out, 'a') as fout:
            _write_frame(fout, ocl.systems[0], ocl.pos_h[:nt].astype(np.float64))
    E5 = ocl.eval_energies_localmd(bFrozenCore=b_frozen).astype(np.float64)
    return dict(pos=ocl.pos_h[:nt].astype(np.float64), E5=E5[0, :5].copy(), na=na, ne=ne)


def _traj_first_divergence(cpu_xyz, gpu_xyz, na, ne, tol_pos=1e-3, tol_size=1e-3):
    # lightweight parser: assumes frames are well-formed and have fixed particle ordering
    def _read_frames(path):
        frames = []
        with open(path, 'r') as f:
            while True:
                l = f.readline()
                if not l:
                    break
                ls = l.strip()
                if not ls:
                    continue
                nt = int(ls.split()[0])
                _ = f.readline()
                data = []
                nread = 0
                while nread < nt:
                    l2 = f.readline()
                    if not l2:
                        break
                    s = l2.split()
                    if len(s) == 0 or s[0].startswith('#'):
                        continue
                    if s[0].startswith('e'):
                        x,y,z = float(s[1]), float(s[2]), float(s[3])
                        sz = float(s[5]) if len(s) > 5 else 0.0
                        data.append([x,y,z,sz])
                    else:
                        x,y,z = float(s[1]), float(s[2]), float(s[3])
                        data.append([x,y,z,0.0])
                    nread += 1
                frames.append(np.array(data, dtype=np.float64))
        return frames

    c = _read_frames(cpu_xyz)
    g = _read_frames(gpu_xyz)
    n = min(len(c), len(g))
    if n == 0:
        return None
    for i in range(n):
        dc = c[i][na:na+ne, :]
        dg = g[i][na:na+ne, :]
        dxyz = np.max(np.abs(dc[:, :3] - dg[:, :3]))
        ds   = np.max(np.abs(dc[:, 3]  - dg[:, 3])) if ne > 0 else 0.0
        if (dxyz > tol_pos) or (ds > tol_size):
            return dict(step=i, pos_xyz_max=float(dxyz), pos_size_max=float(ds))
    return None


def gpu_relax_converge(xyz_in, xyz_out, trj_out, max_iter, dt, damping, nloc, device, fix_atoms=False, ftol=1e-2, check_fixmask=True):
    ocl = eFF_ocl.EFF_OCL(nloc=nloc, device_index=device, bEnergyKernel=True)
    ocl.load_xyzs(xyz_in, bPrint=False)
    if len(ocl.systems) != 1:
        raise RuntimeError(f"gpu_relax_converge expects 1 system, got {len(ocl.systems)} from {xyz_in}")

    if fix_atoms:
        sysdef = ocl.systems[0]
        na = int(sysdef['na']); ne = int(sysdef['ne'])
        fm = np.zeros((na + ne,), dtype=np.uint8)
        fm[:na] = 7
        sysdef['fixmask'] = fm

    ocl.realloc_buffers(); ocl.upload_data()
    b_frozen = all(s.get('core_mode', 'f') == 'f' for s in ocl.systems)
    na = int(ocl.systems[0]['na']); ne = int(ocl.systems[0]['ne'])

    pos0 = ocl.pos_h.copy().astype(np.float64)
    if trj_out:
        open(trj_out, 'w').close()
        with open(trj_out, 'a') as fout:
            _write_frame(fout, ocl.systems[0], pos0[:na+ne])

    nstep = 0
    fe = np.inf
    frc = None
    while nstep < int(max_iter):
        pos = ocl.relax_systems(n_steps=1, dt=dt, damping=damping, bFrozenCore=b_frozen)
        ocl.pos_h[:, :] = pos
        frc = ocl.get_forces().copy().astype(np.float64)
        fe = _maxabs_electron_force(frc, na, ne)
        nstep += 1
        if trj_out:
            with open(trj_out, 'a') as fout:
                _write_frame(fout, ocl.systems[0], ocl.pos_h[:na+ne])
        if fe < ftol:
            break

    posf = ocl.pos_h.copy().astype(np.float64)
    if check_fixmask and fix_atoms:
        _assert_ions_fixed(pos0, posf, na, tol=1e-8, label=os.path.basename(xyz_in))

    E5 = ocl.eval_energies_localmd(bFrozenCore=b_frozen).astype(np.float64)
    with open(xyz_out, 'w') as fout:
        _write_frame(fout, ocl.systems[0], posf[:na+ne])
    return dict(pos=posf, frc=frc, E5=E5, na=na, ne=ne, nstep=nstep, fe=fe)


def cpu_relax_converge_single(xyz_in, xyz_out, trj_out, max_iter, dt, damping, fconv, opt_alg, fix_atoms=False, ftol=1e-2):
    eff.setVerbosity(0, 0)
    eff.setEFFDbgPair(0, 0, 0, 0, 1)
    eff.processXYZ_e(xyz_in, nstepMax=0, dt=dt, Fconv=fconv, optAlg=opt_alg, xyz_out=None, fgo_out=None, bOutputs=(0, 0, 0))
    eff.getBuffs()
    na = int(eff.na); ne = int(eff.ne)

    pos0 = np.zeros((na + ne, 4), dtype=np.float64)
    pos0[:na, :3] = eff.apos
    pos0[na:, :3] = eff.epos
    pos0[na:, 3] = eff.esize

    if fix_atoms:
        fixed_poss = np.zeros((na, 4), dtype=np.float64)
        fixed_inds = np.zeros((na, 2), dtype=np.int32)
        fixed_poss[:, :3] = eff.apos[:, :3]
        fixed_inds[:, 0] = np.arange(na, dtype=np.int32)
        fixed_inds[:, 1] = 7
        eff.set_constrains(na, fixed_poss, fixed_inds, True)

    if trj_out:
        open(trj_out, 'w').close()
    eff.initOpt(dt=dt, damping=damping, f_limit=1000.0, bMass=False)

    nstep = 0
    fe = np.inf
    while nstep < int(max_iter):
        eff.run(nstepMax=1, dt=dt, Fconv=fconv, ialg=opt_alg, bOutE=False, bOutF=False)
        eff.eval()
        eff.getBuffs()
        fe = float(np.max(np.abs(eff.eforce[:, :3]))) if ne > 0 else 0.0
        nstep += 1
        if trj_out:
            pos_step = np.zeros((na + ne, 4), dtype=np.float64)
            pos_step[:na, :3] = eff.apos
            pos_step[na:, :3] = eff.epos
            pos_step[na:, 3] = eff.esize
            # append frame in GPU-readable dialect (na,ne,core + electron records)
            ocl_tmp = eFF_ocl.EFF_OCL(nloc=max(32, na + ne), device_index=0, bEnergyKernel=False)
            ocl_tmp.load_xyzs(xyz_in, bPrint=False)
            if len(ocl_tmp.systems) != 1:
                raise RuntimeError(f"cpu_relax_converge_single: cannot parse template {xyz_in} for trajectory")
            with open(trj_out, 'a') as fout:
                _write_frame(fout, ocl_tmp.systems[0], pos_step[:na+ne])
        if fe < ftol:
            break

    eff.getBuffs()
    posf = np.zeros((na + ne, 4), dtype=np.float64)
    posf[:na, :3] = eff.apos
    posf[na:, :3] = eff.epos
    posf[na:, 3] = eff.esize
    if fix_atoms:
        _assert_ions_fixed(pos0, posf, na, tol=1e-8, label=os.path.basename(xyz_in))

    # Write output in the same XYZ dialect that the GPU parser understands (na,ne,core + e+/e- records)
    _write_xyz_like_input(xyz_in, xyz_out, posf)
    E7 = eff.getEnergyTerms(sh=(7,)).astype(np.float64)
    Ek      = E7[0]
    Eee     = E7[1] + E7[2]
    Eae     = E7[4] + E7[5]
    Eaa     = E7[6]
    Etot    = Ek + Eee + Eae + Eaa
    E5 = np.array([Etot, Ek, Eee, Eae, Eaa], dtype=np.float64)
    return dict(pos=posf, na=na, ne=ne, nstep=nstep, fe=fe, E5=E5)


def _perturb_electrons_in_xyz(xyz_in, xyz_out, dmax=0.05, seed=0):
    rng = np.random.default_rng(int(seed))
    with open(xyz_in, 'r') as f:
        lines = f.readlines()
    if len(lines) < 2:
        raise RuntimeError(f"Bad xyz {xyz_in}")
    nt = int(lines[0].strip().split()[0])

    out = []
    out.append(lines[0])
    out.append(lines[1].rstrip() + f" | e-perturb dmax={dmax} seed={seed}" + "\n")
    ip = 2
    nread = 0
    while nread < nt:
        s = lines[ip].split()
        if len(s) == 0 or s[0].startswith('#'):
            out.append(lines[ip]); ip += 1; continue
        if s[0].startswith('e'):
            r = np.array([float(s[1]), float(s[2]), float(s[3])], dtype=float)
            r += rng.uniform(-dmax, +dmax, size=3)
            s[1] = f"{r[0]:.8f}"; s[2] = f"{r[1]:.8f}"; s[3] = f"{r[2]:.8f}"
            out.append(" ".join(s) + "\n")
        else:
            out.append(lines[ip])
        ip += 1
        nread += 1
    with open(xyz_out, 'w') as f:
        f.writelines(out)


def compare_single(cpu, gpu):
    na, ne = cpu['na'], cpu['ne']
    nt = na + ne
    gpos = gpu['pos'][:nt]
    gfrc = gpu['frc'][:nt]

    dpos_xyz = cpu['pos'][:, :3] - gpos[:, :3]
    dpos_s = cpu['pos'][na:, 3] - gpos[na:, 3]
    dfrc_xyz = cpu['frc'][:, :3] - gfrc[:, :3]
    dfrc_s = cpu['frc'][na:, 3] - gfrc[na:, 3]
    dE = cpu['E5'] - gpu['E5'][0, :5]

    return {
        'pos_xyz_max': float(np.max(np.abs(dpos_xyz))),
        'pos_size_max': float(np.max(np.abs(dpos_s))) if ne > 0 else 0.0,
        'frc_xyz_max': float(np.max(np.abs(dfrc_xyz))),
        'frc_size_max': float(np.max(np.abs(dfrc_s))) if ne > 0 else 0.0,
        'E5_max': float(np.max(np.abs(dE))),
        'dE': dE,
    }


def _perturb_pos(pos, na, ne, sigma=0.05, seed=0, bPerturbAtoms=True, bPerturbElectrons=True, bPerturbSizes=False):
    rng = np.random.default_rng(int(seed))
    out = np.array(pos, copy=True)
    if bPerturbAtoms and na > 0:
        out[:na, :3] += rng.normal(0.0, sigma, size=(na, 3))
    if bPerturbElectrons and ne > 0:
        out[na:na+ne, :3] += rng.normal(0.0, sigma, size=(ne, 3))
        if bPerturbSizes:
            out[na:na+ne, 3] += rng.normal(0.0, sigma, size=(ne,))
    return out


def _write_xyz_from_ocl_systems(fname, systems, pos, comment=""):
    with open(fname, 'w') as fout:
        i0 = 0
        for sysdef in systems:
            nt = int(sysdef['na']) + int(sysdef['ne'])
            if comment:
                fout.write(f"{nt}\n")
                fout.write(comment + "\n")
                for ia in range(int(sysdef['na'])):
                    name = sysdef['ions'][ia][0]
                    x, y, z = pos[i0 + ia, :3]
                    fout.write(f"{name} {x:.8f} {y:.8f} {z:.8f}\n")
                na = int(sysdef['na']); ne = int(sysdef['ne'])
                for ie in range(ne):
                    spin = sysdef['electrons'][ie][1]
                    x, y, z, s = pos[i0 + na + ie]
                    if spin == 0:
                        nm, q = 'e2', 2.0
                    elif spin > 0:
                        nm, q = 'e+', 1.0
                    else:
                        nm, q = 'e-', -1.0
                    fout.write(f"{nm} {x:.8f} {y:.8f} {z:.8f} {q:.1f} {s:.8f}\n")
            else:
                _write_frame(fout, sysdef, pos[i0:i0+nt])
            i0 += nt


def _maxabs(a):
    return float(np.max(np.abs(a))) if a.size else 0.0


def _snapshot_cpu():
    eff.getBuffs()
    na, ne = eff.na, eff.ne
    pos = np.zeros((na + ne, 4), dtype=np.float64)
    pos[:na, :3] = eff.apos
    pos[na:, :3] = eff.epos
    pos[na:, 3] = eff.esize
    frc = np.zeros((na + ne, 4), dtype=np.float64)
    frc[:na, :3] = eff.aforce
    frc[na:, :3] = eff.eforce
    frc[na:, 3] = eff.fsize
    # getEnergyTerms(): [Ek, Eee, EeePaul, EeeExch, Eae, EaePaul, Eaa]
    E7 = eff.getEnergyTerms(sh=(7,)).astype(np.float64)
    Ek      = E7[0]
    Eee     = E7[1] + E7[2]
    Eae     = E7[4] + E7[5]
    Eaa     = E7[6]
    Etot    = Ek + Eee + Eae + Eaa
    E5 = np.array([Etot, Ek, Eee, Eae, Eaa], dtype=np.float64)
    return dict(na=na, ne=ne, pos=pos, frc=frc, E5=E5)


def _count_xyz_confs(xyz_path):
    nconf = 0
    with open(xyz_path, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            ls = line.strip()
            if not ls:
                continue
            try:
                ntot = int(ls.split()[0])
            except Exception:
                continue
            comment = f.readline()
            if not comment:
                break
            # consume exactly ntot non-empty, non-comment particle lines
            nread = 0
            while nread < ntot:
                l2 = f.readline()
                if not l2:
                    break
                s2 = l2.strip().split()
                if (len(s2) == 0) or s2[0].startswith('#'):
                    continue
                nread += 1
            nconf += 1
    return nconf


def _extract_xyz_conf(xyz_in, iconf, xyz_out):
    iconf = int(iconf)
    if iconf < 0:
        raise ValueError(f"iconf must be >=0, got {iconf}")
    with open(xyz_in, 'r') as fin:
        ic = -1
        while True:
            l = fin.readline()
            if not l:
                break
            ls = l.strip()
            if not ls:
                continue
            try:
                nt = int(ls.split()[0])
            except Exception:
                raise RuntimeError(f"Failed to parse XYZ natom/e count line: {l!r} in {xyz_in}")
            comment = fin.readline()
            if not comment:
                raise RuntimeError(f"Unexpected EOF after comment line in {xyz_in}")
            lines = [l, comment]
            for _ in range(nt):
                li = fin.readline()
                if not li:
                    raise RuntimeError(f"Unexpected EOF while reading frame nt={nt} in {xyz_in}")
                lines.append(li)
            ic += 1
            if ic == iconf:
                with open(xyz_out, 'w') as fout:
                    fout.writelines(lines)
                return xyz_out
    raise RuntimeError(f"iconf={iconf} out of range for {xyz_in}")


def _snapshot_gpu(ocl):
    forces = ocl.get_forces().copy().astype(np.float64)
    pos = ocl.pos_h.copy().astype(np.float64)
    E5 = ocl.eval_energies_localmd(bFrozenCore=any(s.get('core_mode','f') == 'f' for s in ocl.systems)).astype(np.float64)
    return dict(pos=pos, frc=forces, E5=E5)


def cpu_relax_scan_loop_cpp(xyz_template, systems, pos_all, nsteps, dt, damping, fix_atoms=True, size_min=0.001):
    # CPU scan relaxation using C++ integrator ialg=3 (move_MD_dbg) with per-frame initialization from pos_all.
    # systems: list of parsed system defs (from EFF_OCL)
    # pos_all: float64 [sum_nt,4] concatenated blocks matching systems ordering
    assert len(systems) > 0
    na = int(systems[0]['na']); ne = int(systems[0]['ne']); nt = na + ne
    for s in systems:
        if int(s['na']) != na or int(s['ne']) != ne:
            raise RuntimeError("cpu_relax_scan_loop_cpp expects constant na/ne across scan")

    # init CPU once from a single-config template to allocate buffers and set coreMode
    eff.setVerbosity(0, 0)
    eff.setEFFDbgPair(0, 0, 0, 0, 1)
    eff.processXYZ_e(xyz_template, nstepMax=0, dt=dt, Fconv=1e-3, optAlg=3, xyz_out=None, fgo_out=None, bOutputs=(0, 0, 0))
    eff.getBuffs()
    if int(eff.na) != na or int(eff.ne) != ne:
        raise RuntimeError(f"CPU init mismatch na/ne: cpu=({eff.na},{eff.ne}) scan=({na},{ne})")

    fixed_inds = None
    if fix_atoms and na > 0:
        fixed_inds = np.zeros((na, 2), dtype=np.int32)
        fixed_inds[:, 0] = np.arange(na, dtype=np.int32)
        fixed_inds[:, 1] = 7

    nconf = len(systems)
    E5s = np.zeros((nconf, 5), dtype=np.float64)
    pos_fin = np.zeros((nconf, nt, 4), dtype=np.float64)
    i0 = 0
    eff.initOpt(dt=dt, damping=damping, f_limit=1000.0, bMass=False)
    for ic in range(nconf):
        blk = pos_all[i0:i0+nt].astype(np.float64, copy=False)
        i0 += nt
        if na > 0:
            eff.apos[:, :3] = blk[:na, :3]
        if ne > 0:
            eff.epos[:, :3] = blk[na:, :3]
            eff.esize[:] = np.maximum(blk[na:, 3], float(size_min))

        if fix_atoms and na > 0:
            fixed_poss = np.zeros((na, 4), dtype=np.float64)
            fixed_poss[:, :3] = eff.apos[:, :3]
            eff.set_constrains(na, fixed_poss, fixed_inds, True)

        eff.run(nstepMax=int(nsteps), dt=dt, Fconv=0.0, ialg=3, bOutE=False, bOutF=False)
        eff.eval(); eff.getBuffs()
        E7 = eff.getEnergyTerms(sh=(7,)).astype(np.float64)
        Ek      = E7[0]
        Eee     = E7[1] + E7[2]
        Eae     = E7[4] + E7[5]
        Eaa     = E7[6]
        Etot    = Ek + Eee + Eae + Eaa
        E5s[ic, :] = (Etot, Ek, Eee, Eae, Eaa)
        if na > 0:
            pos_fin[ic, :na, :3] = eff.apos[:, :3]
        if ne > 0:
            pos_fin[ic, na:, :3] = eff.epos[:, :3]
            pos_fin[ic, na:, 3] = eff.esize[:]
    return E5s, pos_fin


def _write_xyz_frame(f, pos, na, ne, comment=""):
    nt = na + ne
    assert pos.shape[0] >= nt
    f.write(f"{nt}\n")
    f.write(comment + "\n")
    for ia in range(na):
        x,y,z = pos[ia,0], pos[ia,1], pos[ia,2]
        # element symbols are not carried here; write generic 'A'
        f.write(f"A {x:.16f} {y:.16f} {z:.16f}\n")
    for ie in range(ne):
        ip = na + ie
        x,y,z,s = pos[ip,0], pos[ip,1], pos[ip,2], pos[ip,3]
        f.write(f"e {x:.16f} {y:.16f} {z:.16f} 0 {s:.16f}\n")


def short_stepwise_debug(xyz_in, outdir, nsteps, dt, damping, fconv, opt_alg, nloc, device, tol_f=1e-4, tol_e=1e-4, tol_p=1e-3, tol_s=1e-3, dbg_kind=1, dbg_i=0, dbg_j=1, dbg_allpairs=False, fix_atoms=False, xyz_traj_cpu=None, xyz_traj_gpu=None):
    os.makedirs(outdir, exist_ok=True)

    nconf = _count_xyz_confs(xyz_in)
    if nconf != 1:
        raise RuntimeError(f"--stepwise requires single-config xyz, got nconf={nconf} in {xyz_in}")

    # --- init CPU with electrons
    eff.setVerbosity(0, 0)
    eff.setEFFDbgPair(0, 0, 0, 0, 1)
    eff.processXYZ_e(xyz_in, nstepMax=0, dt=dt, Fconv=fconv, optAlg=opt_alg, xyz_out=None, fgo_out=None, bOutputs=(0, 0, 0))
    if fix_atoms:
        eff.getBuffs()
        na_ = int(eff.na); ne_ = int(eff.ne)
        fixed_poss = np.zeros((na_, 4), dtype=np.float64)
        fixed_inds = np.zeros((na_, 2), dtype=np.int32)
        fixed_poss[:, :3] = eff.apos[:, :3]
        fixed_inds[:, 0] = np.arange(na_, dtype=np.int32)
        fixed_inds[:, 1] = 7
        eff.set_constrains(na_, fixed_poss, fixed_inds, True)
    eff.initOpt(dt=dt, damping=damping, f_limit=1000.0, bMass=False)
    cpu0 = _snapshot_cpu()
    na, ne = cpu0['na'], cpu0['ne']
    nt = na + ne

    # --- init GPU (no debug build first)
    ocl = eFF_ocl.EFF_OCL(nloc=nloc, device_index=device, bEnergyKernel=True)
    ocl.load_xyzs(xyz_in, bPrint=False)
    if fix_atoms:
        for sys in ocl.systems:
            na = int(sys['na']); ne = int(sys['ne'])
            fm = np.zeros((na+ne,), dtype=np.uint8)
            fm[:na] = np.uint8(7)  # fix xyz for all atoms
            sys['fixmask'] = fm
    if len(ocl.systems) != 1:
        raise RuntimeError(f"short_stepwise_debug expects 1 system, got {len(ocl.systems)}")
    ocl.realloc_buffers()
    ocl.upload_data()
    b_frozen = all(s.get('core_mode', 'f') == 'f' for s in ocl.systems)

    diverge_step = None
    diverge_metric = None

    fcpu = open(xyz_traj_cpu, 'w') if xyz_traj_cpu else None
    fgpu = open(xyz_traj_gpu, 'w') if xyz_traj_gpu else None

    # write initial frame
    if fcpu: _write_xyz_frame(fcpu, cpu0['pos'], na, ne, comment=f"cpu step=-1")
    if fgpu: _write_xyz_frame(fgpu, ocl.pos_h.astype(np.float64), na, ne, comment=f"gpu step=-1")

    for istep in range(nsteps):
        # --- Pre-move compare: evaluate energies/forces at the same geometry
        eff.eval()
        cpu_pre = _snapshot_cpu()

        # Force path parity: use localMD to compute forces with dt=0 (no movement) and damping=1 (no vel change)
        pos_tmp = ocl.relax_systems(n_steps=1, dt=0.0, damping=1.0, bFrozenCore=b_frozen)
        ocl.pos_h[:, :] = pos_tmp
        E5g_pre = ocl.eval_energies_localmd(bFrozenCore=b_frozen)
        gpu_pre = dict(pos=ocl.pos_h.copy().astype(np.float64), frc=ocl.get_forces().copy().astype(np.float64), E5=E5g_pre.astype(np.float64))

        if fix_atoms:
            dfrc_pre = cpu_pre['frc'][na:nt] - gpu_pre['frc'][na:nt]
        else:
            dfrc_pre = cpu_pre['frc'][:nt] - gpu_pre['frc'][:nt]
        dE_pre = cpu_pre['E5'][:5] - gpu_pre['E5'][0, :5]

        # --- Do one integration step on both sides
        eff.run(nstepMax=1, dt=dt, Fconv=fconv, ialg=opt_alg, bOutE=False, bOutF=False)
        cpu_post = _snapshot_cpu()

        pos_out = ocl.relax_systems(n_steps=1, dt=dt, damping=damping, bFrozenCore=b_frozen)
        ocl.pos_h[:, :] = pos_out
        gpu_post = dict(pos=ocl.pos_h.copy().astype(np.float64))

        if fcpu: _write_xyz_frame(fcpu, cpu_post['pos'], na, ne, comment=f"cpu step={istep}")
        if fgpu: _write_xyz_frame(fgpu, gpu_post['pos'], na, ne, comment=f"gpu step={istep}")

        dpos_post = cpu_post['pos'][:nt] - gpu_post['pos'][:nt]

        metrics = {
            'istep': istep,
            'pos_xyz_max': _maxabs(dpos_post[:, :3]),
            'pos_s_max': _maxabs(dpos_post[na:, 3]) if ne > 0 else 0.0,
            'frc_xyz_max': _maxabs(dfrc_pre[:, :3]),
            'frc_s_max': _maxabs(dfrc_pre[:, 3]) if ne > 0 else 0.0,
            'E5_max': _maxabs(dE_pre),
            'dE': dE_pre,
        }

        with open(os.path.join(outdir, 'stepwise_metrics.txt'), 'a') as f:
            f.write(repr(metrics) + "\n")

        # Protocol: for short trajectories, treat energies/forces as *diagnostics*, not hard acceptance criteria.
        # Divergence trigger is positions and electron sizes only (chaos/noise grows quickly).
        diverged = (metrics['pos_xyz_max'] > tol_p) or (metrics['pos_s_max'] > tol_s)
        if diverged:
            diverge_step = istep
            diverge_metric = metrics
            np.savez(
                os.path.join(outdir, f'diverge_step_{istep:04d}.npz'),
                cpu_pre_pos=cpu_pre['pos'], cpu_pre_frc=cpu_pre['frc'], cpu_pre_E5=cpu_pre['E5'],
                cpu_post_pos=cpu_post['pos'], cpu_post_frc=cpu_post['frc'], cpu_post_E5=cpu_post['E5'],
                gpu_pre_pos=gpu_pre['pos'], gpu_pre_frc=gpu_pre['frc'], gpu_pre_E5=gpu_pre['E5'],
                gpu_post_pos=gpu_post['pos'],
            )
            break

    if diverge_step is None:
        if fcpu: fcpu.close()
        if fgpu: fgpu.close()
        return dict(status='ok', diverge_step=None)

    # --- rerun with prints at the diverging step only
    # Choose dbg_i/dbg_j in the shared indexing: ions [0..na-1], electrons [na..na+ne-1]
    eff.setVerbosity(0, 0)
    if dbg_allpairs:
        eff.setEFFDbgAllPairs(1)
    eff.setEFFDbgPair(1, diverge_step, dbg_kind, dbg_i, dbg_j)
    eff.processXYZ_e(xyz_in, nstepMax=diverge_step + 1, dt=dt, Fconv=fconv, optAlg=opt_alg, xyz_out=None, fgo_out=None, bOutputs=(0, 0, 0))
    eff.setEFFDbgPair(0, 0, 0, 0, 1)
    if dbg_allpairs:
        eff.setEFFDbgAllPairs(0)

    ocl_dbg = eFF_ocl.EFF_OCL(nloc=nloc, device_index=device, bEnergyKernel=True, dbg_pair=True, idbg_sys=0, idbg_step=diverge_step, idbg_i=dbg_i, idbg_j=dbg_j, dbg_allpairs=dbg_allpairs)
    ocl_dbg.load_xyzs(xyz_in, bPrint=False)
    if fix_atoms:
        for sys in ocl_dbg.systems:
            na = int(sys['na']); ne = int(sys['ne'])
            fm = np.zeros((na+ne,), dtype=np.uint8)
            fm[:na] = np.uint8(7)
            sys['fixmask'] = fm
    ocl_dbg.realloc_buffers(); ocl_dbg.upload_data()
    b_frozen = all(s.get('core_mode', 'f') == 'f' for s in ocl_dbg.systems)
    pos_out = ocl_dbg.relax_systems(n_steps=diverge_step + 1, dt=dt, damping=damping, bFrozenCore=b_frozen)
    ocl_dbg.pos_h[:, :] = pos_out

    return dict(status='diverged', diverge_step=diverge_step, metrics=diverge_metric)


def main():
    ap = argparse.ArgumentParser(description='Run eFF CPU/GPU relaxation parity protocol: short, long, scan')
    ap.add_argument('--single-xyz', default='../../cpp/sketches_SDL/Molecular/data/H2O_fixcore.xyz')
    ap.add_argument('--scan-xyz', default='export/scan_data/distscan_H2O__spins_fc.xyz')
    ap.add_argument('--outdir', default='export/relax_parity')
    ap.add_argument('--short-steps', type=int, default=10)
    ap.add_argument('--long-steps', type=int, default=1000)
    ap.add_argument('--scan-steps', type=int, default=100)
    ap.add_argument('--dt', type=float, default=0.001)
    ap.add_argument('--damping', type=float, default=0.1)
    ap.add_argument('--fconv', type=float, default=1e-3)
    ap.add_argument('--opt-alg', type=int, default=1, help='CPU optimizer: 1=MDdamp, 2=FIRE. For trajectory parity with GPU localMD use 1.')
    ap.add_argument('--nloc', type=int, default=32)
    ap.add_argument('--device', type=int, default=0)
    ap.add_argument('--skip-scan', action='store_true')
    ap.add_argument('--stepwise', action='store_true')
    ap.add_argument('--tol-f', type=float, default=1e-4)
    ap.add_argument('--tol-e', type=float, default=1e-4)
    ap.add_argument('--tol-p', type=float, default=1e-3, help='Position tolerance [A] for stepwise divergence trigger (use ~1e-3 for 10 steps).')
    ap.add_argument('--tol-s', type=float, default=1e-3, help='Electron size tolerance for stepwise divergence trigger (use ~1e-3).')
    ap.add_argument('--dbg-kind', type=int, default=1, help='1=EE 2=AE 3=AA')
    ap.add_argument('--dbg-i', type=int, default=0)
    ap.add_argument('--dbg-j', type=int, default=1)
    ap.add_argument('--dbg-allpairs', action='store_true', help='On divergence rerun, dump all pairs for selected kind on both CPU/GPU (CH4-size manageable).')
    ap.add_argument('--fix-atoms', action='store_true', help='Fix all nuclei (atoms) on GPU (xyz mask=7); relax only electrons. CPU already supports fixed DOFs separately.')
    ap.add_argument('--xyz-traj-cpu', default=None, help='Write CPU stepwise trajectory xyz here (only in --stepwise).')
    ap.add_argument('--xyz-traj-gpu', default=None, help='Write GPU stepwise trajectory xyz here (only in --stepwise).')
    ap.add_argument('--perturb-sigma', type=float, default=0.05, help='Perturbation sigma [A] for long-minimum test (Gaussian).')
    ap.add_argument('--perturb-seed', type=int, default=0, help='RNG seed for perturbation.')
    ap.add_argument('--perturb-electrons-only', action='store_true', help='Only perturb electrons (not atoms) in long-minimum test.')
    ap.add_argument('--scan-converge', action='store_true', help='Run 1D-scan convergence protocol: CPU relax -> GPU relax -> GPU relax from CPU-min -> GPU relax from perturbed CPU-min. Requires --fix-atoms and single-config scan.')
    ap.add_argument('--max-iter', type=int, default=10000)
    ap.add_argument('--ftol-elec', type=float, default=1e-2, help='Convergence threshold for electron xyz force component (max abs).')
    ap.add_argument('--perturb-dmax', type=float, default=0.05, help='Uniform perturbation box half-width [A] for electron positions.')
    ap.add_argument('--long-parity', action='store_true', help='Run long trajectory parity using CPU mirror-integrator vs GPU localMD from the same initial xyz. Writes cpu/gpu trajectories and a divergence report.')
    ap.add_argument('--parity-steps', type=int, default=2000)
    ap.add_argument('--tol-parity-pos', type=float, default=1e-3)
    ap.add_argument('--tol-parity-size', type=float, default=1e-3)
    ap.add_argument('--iconf', type=int, default=0, help='Configuration index to use from multi-frame XYZ (only used for --long-parity and --stepwise).')
    ap.add_argument('--cpu-parity-backend', choices=['cpp', 'py'], default='cpp', help='CPU backend for --long-parity: cpp uses ialg=3 (move_MD_dbg) and C++ trajectory; py uses python mirror loop (debug/reference).')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    single_xyz = os.path.abspath(os.path.join(os.path.dirname(__file__), args.single_xyz)) if not os.path.isabs(args.single_xyz) else args.single_xyz
    scan_xyz = os.path.abspath(os.path.join(os.path.dirname(__file__), args.scan_xyz)) if not os.path.isabs(args.scan_xyz) else args.scan_xyz

    if args.stepwise:
        if args.iconf != 0:
            tmp_xyz = os.path.join(args.outdir, f'single_iconf_{args.iconf}.xyz')
            _extract_xyz_conf(single_xyz, args.iconf, tmp_xyz)
            single_xyz = tmp_xyz
        # Short, rigorous mode: detect first divergence and print one selected pair/kind
        res = short_stepwise_debug(single_xyz, args.outdir, args.short_steps, args.dt, args.damping, args.fconv, args.opt_alg, args.nloc, args.device, tol_f=args.tol_f, tol_e=args.tol_e, tol_p=args.tol_p, tol_s=args.tol_s, dbg_kind=args.dbg_kind, dbg_i=args.dbg_i, dbg_j=args.dbg_j, dbg_allpairs=args.dbg_allpairs, fix_atoms=args.fix_atoms, xyz_traj_cpu=args.xyz_traj_cpu, xyz_traj_gpu=args.xyz_traj_gpu)
        print('STEPWISE', res)
        return

    if args.scan_converge:
        if not args.fix_atoms:
            raise RuntimeError("--scan-converge requires --fix-atoms")
        nconf = _count_xyz_confs(scan_xyz)
        if nconf != 1:
            raise RuntimeError(f"--scan-converge requires single-config xyz (1 frame), got nconf={nconf} in {scan_xyz}")

        cpu_trj = os.path.join(args.outdir, 'scanconv_cpu_trj.xyz')
        cpu_fin = os.path.join(args.outdir, 'scanconv_cpu_final.xyz')
        gpu_trj = os.path.join(args.outdir, 'scanconv_gpu_trj.xyz')
        gpu_fin = os.path.join(args.outdir, 'scanconv_gpu_final.xyz')
        gpu_from_cpu_trj = os.path.join(args.outdir, 'scanconv_gpu_from_cpu_trj.xyz')
        gpu_from_cpu_fin = os.path.join(args.outdir, 'scanconv_gpu_from_cpu_final.xyz')
        pert_xyz = os.path.join(args.outdir, 'scanconv_cpu_final_pert.xyz')
        gpu_from_pert_trj = os.path.join(args.outdir, 'scanconv_gpu_from_pert_trj.xyz')
        gpu_from_pert_fin = os.path.join(args.outdir, 'scanconv_gpu_from_pert_final.xyz')

        # 1) CPU converge
        cpu = cpu_relax_converge_single(scan_xyz, cpu_fin, cpu_trj, args.max_iter, args.dt, args.damping, args.fconv, args.opt_alg, fix_atoms=True, ftol=args.ftol_elec)

        # 2) GPU converge from original
        gpu = gpu_relax_converge(scan_xyz, gpu_fin, gpu_trj, args.max_iter, args.dt, args.damping, args.nloc, args.device, fix_atoms=True, ftol=args.ftol_elec, check_fixmask=True)

        # 3) GPU converge starting from CPU relaxed geometry
        gpu_from_cpu = gpu_relax_converge(cpu_fin, gpu_from_cpu_fin, gpu_from_cpu_trj, args.max_iter, args.dt, args.damping, args.nloc, args.device, fix_atoms=True, ftol=args.ftol_elec, check_fixmask=True)

        # 4) Perturb CPU relaxed electrons and relax on GPU
        _perturb_electrons_in_xyz(cpu_fin, pert_xyz, dmax=args.perturb_dmax, seed=args.perturb_seed)
        gpu_from_pert = gpu_relax_converge(pert_xyz, gpu_from_pert_fin, gpu_from_pert_trj, args.max_iter, args.dt, args.damping, args.nloc, args.device, fix_atoms=True, ftol=args.ftol_elec, check_fixmask=True)

        rep = os.path.join(args.outdir, 'scanconv_report.txt')
        with open(rep, 'w') as f:
            f.write('=== SCAN-CONVERGE (single frame) ===\n')
            f.write(f"CPU: nstep={cpu['nstep']} fe={cpu['fe']} E5={cpu['E5']}\n")
            f.write(f"GPU: nstep={gpu['nstep']} fe={gpu['fe']} E5={gpu['E5'][0,:5]}\n")
            f.write(f"GPU_from_CPU: nstep={gpu_from_cpu['nstep']} fe={gpu_from_cpu['fe']} E5={gpu_from_cpu['E5'][0,:5]}\n")
            f.write(f"GPU_from_pert: nstep={gpu_from_pert['nstep']} fe={gpu_from_pert['fe']} E5={gpu_from_pert['E5'][0,:5]}\n")
        print('WROTE', rep)
        return

    if args.long_parity:
        if not args.fix_atoms:
            raise RuntimeError("--long-parity requires --fix-atoms (electron-only relax parity)")

        nconf = _count_xyz_confs(single_xyz)
        if nconf != 1:
            tmp_xyz = os.path.join(args.outdir, f'longpar_iconf_{args.iconf}.xyz')
            _extract_xyz_conf(single_xyz, args.iconf, tmp_xyz)
            single_xyz = tmp_xyz

        cpu_trj = os.path.join(args.outdir, 'longpar_cpu_trj.xyz')
        gpu_trj = os.path.join(args.outdir, 'longpar_gpu_trj.xyz')
        if args.cpu_parity_backend == 'cpp':
            cpu = cpu_localmd_traj_cpp(single_xyz, cpu_trj, args.parity_steps, args.dt, args.damping, fix_atoms=True)
        else:
            cpu = cpu_md_mirror_traj(single_xyz, cpu_trj, args.parity_steps, args.dt, args.damping, fix_atoms=True)
        gpu = gpu_localmd_traj(single_xyz, gpu_trj, args.parity_steps, args.dt, args.damping, args.nloc, args.device, fix_atoms=True)
        na = int(cpu['na']); ne = int(cpu['ne'])
        div = _traj_first_divergence(cpu_trj, gpu_trj, na, ne, tol_pos=args.tol_parity_pos, tol_size=args.tol_parity_size)
        rep = os.path.join(args.outdir, 'longpar_report.txt')
        with open(rep, 'w') as f:
            f.write('=== LONG-PARITY (CPU mirror vs GPU localMD) ===\n')
            f.write(f"steps={args.parity_steps} dt={args.dt} damping={args.damping}\n")
            f.write(f"CPU_E5={cpu['E5']}\n")
            f.write(f"GPU_E5={gpu['E5']}\n")
            f.write(f"tol_pos={args.tol_parity_pos} tol_size={args.tol_parity_size}\n")
            f.write(f"divergence={div}\n")
        print('WROTE', rep)
        return

    cmp_short = None
    if args.short_steps > 0:
        # ---- Step 1: short trajectory parity
        cpu_short = cpu_relax_single(single_xyz, os.path.join(args.outdir, 'cpu_short_final.xyz'), os.path.join(args.outdir, 'cpu_short_trj.xyz'), args.short_steps, args.dt, args.fconv, args.opt_alg)
        gpu_short = gpu_relax(single_xyz, os.path.join(args.outdir, 'gpu_short_final.xyz'), os.path.join(args.outdir, 'gpu_short_trj.xyz'), args.short_steps, args.dt, args.damping, args.nloc, args.device, fix_atoms=args.fix_atoms)
        cmp_short = compare_single(cpu_short, gpu_short)

    cmp_long = None
    cmp_from_pert = None
    if args.long_steps > 0:
        # ---- Step 2: long run parity
        cpu_long = cpu_relax_single(single_xyz, os.path.join(args.outdir, 'cpu_long_final.xyz'), os.path.join(args.outdir, 'cpu_long_trj.xyz'), args.long_steps, args.dt, args.fconv, args.opt_alg)
        gpu_long = gpu_relax(single_xyz, os.path.join(args.outdir, 'gpu_long_final.xyz'), os.path.join(args.outdir, 'gpu_long_trj.xyz'), args.long_steps, args.dt, args.damping, args.nloc, args.device, fix_atoms=args.fix_atoms)
        cmp_long = compare_single(cpu_long, gpu_long)

        # ---- Step 2b: long-minimum robustness test (CPU minimum -> perturb -> GPU relax)
        # Use CPU-relaxed geometry as a starting point, perturb it, then relax on GPU.
        ocl_tmp = eFF_ocl.EFF_OCL(nloc=args.nloc, device_index=args.device, bEnergyKernel=True)
        cpu_min_xyz = os.path.join(args.outdir, 'cpu_long_final.xyz')
        ocl_tmp.load_xyzs(cpu_min_xyz, bPrint=False)
        if len(ocl_tmp.systems) != 1:
            raise RuntimeError(f"Expected exactly 1 system in {cpu_min_xyz}, got {len(ocl_tmp.systems)}")
        ocl_tmp.realloc_buffers(); ocl_tmp.upload_data()
        na = int(ocl_tmp.systems[0]['na']); ne = int(ocl_tmp.systems[0]['ne'])
        pos0 = ocl_tmp.pos_h.astype(np.float64)
        posp = _perturb_pos(pos0, na, ne, sigma=args.perturb_sigma, seed=args.perturb_seed, bPerturbAtoms=(not args.perturb_electrons_only), bPerturbElectrons=True, bPerturbSizes=False)
        xyz_pert = os.path.join(args.outdir, 'cpu_long_final_perturbed.xyz')
        _write_xyz_from_ocl_systems(xyz_pert, ocl_tmp.systems, posp, comment=f"na,ne,core {na} {ne} {ocl_tmp.systems[0].get('core_mode','f')} | perturbed sigma={args.perturb_sigma} seed={args.perturb_seed}")
        gpu_from_pert = gpu_relax(xyz_pert, os.path.join(args.outdir, 'gpu_from_pert_final.xyz'), os.path.join(args.outdir, 'gpu_from_pert_trj.xyz'), args.long_steps, args.dt, args.damping, args.nloc, args.device, fix_atoms=args.fix_atoms)
        cmp_from_pert = compare_single(cpu_long, gpu_from_pert)

    scan_stats = {'nconf': 0, 'Etot_max': np.nan, 'Ek_max': np.nan, 'Eee_max': np.nan, 'Eae_max': np.nan, 'Eaa_max': np.nan, 'pos_xyz_max': np.nan, 'pos_size_max': np.nan}
    if not args.skip_scan:
        # ---- Step 3: scan-level parallel GPU and CPU reference
        eff.setVerbosity(0, 0)
        eff.setTrjName(os.path.join(args.outdir, 'cpu_scan_trj.xyz'), savePerNsteps=10, bDel=True)

        if args.fix_atoms:
            ocl_scan_ref = eFF_ocl.EFF_OCL(nloc=args.nloc, device_index=args.device, bEnergyKernel=True)
            ocl_scan_ref.load_xyzs(scan_xyz, bPrint=False)
            if len(ocl_scan_ref.systems) == 0:
                raise RuntimeError(f"No systems in scan {scan_xyz}")

            ocl_scan_ref.realloc_buffers()
            ocl_scan_ref.upload_data()

            # Initialize CPU from extracted frame0 to avoid very verbose multi-frame processXYZ_e
            tmp0 = os.path.join(args.outdir, 'scan_init_iconf0.xyz')
            _extract_xyz_conf(scan_xyz, 0, tmp0)

            # CPU relax each frame with C++ integrator ialg=3
            Es5_cpu, pos_cpu = cpu_relax_scan_loop_cpp(tmp0, ocl_scan_ref.systems, ocl_scan_ref.pos_h.astype(np.float64), args.scan_steps, args.dt, args.damping, fix_atoms=True)

            # GPU relax all frames in parallel
            gpu_scan = gpu_relax(scan_xyz, os.path.join(args.outdir, 'gpu_scan_relaxed.xyz'), None, args.scan_steps, args.dt, args.damping, args.nloc, args.device, fix_atoms=True)
            Es5_gpu = gpu_scan['E5'][:, :5]

            n = min(len(Es5_cpu), len(Es5_gpu))
            dscan = Es5_cpu[:n, :5] - Es5_gpu[:n, :5]
            scan_stats.update({
                'nconf': int(n),
                'Etot_max': float(np.max(np.abs(dscan[:, 0]))) if n > 0 else np.nan,
                'Ek_max': float(np.max(np.abs(dscan[:, 1]))) if n > 0 else np.nan,
                'Eee_max': float(np.max(np.abs(dscan[:, 2]))) if n > 0 else np.nan,
                'Eae_max': float(np.max(np.abs(dscan[:, 3]))) if n > 0 else np.nan,
                'Eaa_max': float(np.max(np.abs(dscan[:, 4]))) if n > 0 else np.nan,
            })

            na_scan = int(ocl_scan_ref.systems[0]['na']); ne_scan = int(ocl_scan_ref.systems[0]['ne'])
            pos_gpu = gpu_scan['pos'][:n*(na_scan+ne_scan), :].reshape((n, na_scan+ne_scan, 4))
            dpos_xyz = pos_cpu[:n, :, :3] - pos_gpu[:, :, :3]
            dpos_s = pos_cpu[:n, na_scan:, 3] - pos_gpu[:, na_scan:, 3]
            scan_stats['pos_xyz_max'] = float(np.max(np.abs(dpos_xyz))) if n > 0 else np.nan
            scan_stats['pos_size_max'] = float(np.max(np.abs(dpos_s))) if n > 0 else np.nan

            np.save(os.path.join(args.outdir, 'Es5_cpu.npy'), Es5_cpu)
            np.save(os.path.join(args.outdir, 'Es5_gpu.npy'), Es5_gpu)
            np.save(os.path.join(args.outdir, 'Es5_diff.npy'), dscan)
        else:
            Es_cpu_scan, _, _ = eff.processXYZ_e(scan_xyz, nstepMax=args.scan_steps, dt=args.dt, Fconv=args.fconv, optAlg=args.opt_alg, xyz_out=os.path.join(args.outdir, 'cpu_scan_relaxed.xyz'), fgo_out=None, bOutputs=(1, 0, 0))
            Es_cpu_scan = np.array(Es_cpu_scan, dtype=np.float64)

            gpu_scan = gpu_relax(scan_xyz, os.path.join(args.outdir, 'gpu_scan_relaxed.xyz'), None, args.scan_steps, args.dt, args.damping, args.nloc, args.device, fix_atoms=False)
            Es_gpu_scan = gpu_scan['E5']

            n = min(len(Es_cpu_scan), len(Es_gpu_scan))
            dscan = Es_cpu_scan[:n, :5] - Es_gpu_scan[:n, :5]
            scan_stats.update({
                'nconf': int(n),
                'Etot_max': float(np.max(np.abs(dscan[:, 0]))) if n > 0 else np.nan,
                'Ek_max': float(np.max(np.abs(dscan[:, 1]))) if n > 0 else np.nan,
                'Eee_max': float(np.max(np.abs(dscan[:, 2]))) if n > 0 else np.nan,
                'Eae_max': float(np.max(np.abs(dscan[:, 3]))) if n > 0 else np.nan,
                'Eaa_max': float(np.max(np.abs(dscan[:, 4]))) if n > 0 else np.nan,
            })

    with open(os.path.join(args.outdir, 'relax_parity_report.txt'), 'w') as f:
        f.write('=== SHORT ===\n')
        if cmp_short is not None:
            for k, v in cmp_short.items():
                f.write(f"{k}: {v}\n")
        else:
            f.write('SKIPPED\n')

        f.write('\n=== LONG ===\n')
        if cmp_long is not None:
            for k, v in cmp_long.items():
                f.write(f"{k}: {v}\n")
        else:
            f.write('SKIPPED\n')

        f.write('\n=== LONG_MINIMUM_ROBUSTNESS ===\n')
        if cmp_from_pert is not None:
            for k, v in cmp_from_pert.items():
                f.write(f"{k}: {v}\n")
        else:
            f.write('SKIPPED\n')

        f.write('\n=== SCAN_PARITY ===\n')
        for k, v in scan_stats.items():
            f.write(f"{k}: {v}\n")

    print('WROTE', os.path.join(args.outdir, 'relax_parity_report.txt'))
    print('SHORT', cmp_short)
    if cmp_long is not None:
        print('LONG ', cmp_long)
    if cmp_from_pert is not None:
        print('PERT ', cmp_from_pert)
    print('SCAN ', scan_stats)


if __name__ == '__main__':
    main()
