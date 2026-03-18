import ctypes
import ctypes.util
import json
import math
import re
import subprocess
import sys
from copy import deepcopy
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest


def _opencl_available() -> bool:
    libname = ctypes.util.find_library("OpenCL")
    if not libname:
        return False
    try:
        cl = ctypes.CDLL(libname)
    except OSError:
        return False
    cl.clGetPlatformIDs.argtypes = [
        ctypes.c_uint,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.POINTER(ctypes.c_uint),
    ]
    cl.clGetPlatformIDs.restype = ctypes.c_int
    num = ctypes.c_uint(0)
    err = cl.clGetPlatformIDs(0, None, ctypes.byref(num))
    if err == -1001:  # CL_PLATFORM_NOT_FOUND_KHR
        return False
    if err != 0:
        return False
    return num.value > 0


OPENCL_AVAILABLE = _opencl_available()
pytestmark = pytest.mark.skipif(not OPENCL_AVAILABLE, reason="OpenCL platform not available")


REPO_ROOT = Path(__file__).resolve().parents[2]
CPP_DIR = REPO_ROOT / "cpp"


CASE_SCRIPT = r'''
import json, math, os, pathlib, sys, tempfile
import numpy as np
repo_root = pathlib.Path(sys.argv[1])
mode = sys.argv[2]
nsys = int(sys.argv[3])
temperature = float(sys.argv[4])
sys.path.insert(0, str(repo_root))
from pyBall import MMFF_multi as mm

base = repo_root / "cpp" / "common_resources"
xyz = base / "entropic_spring_20.xyz"
DT = np.float32(0.05)
FLIMIT = np.float32(10.0)

full_initial = [0.5, 0.0, 0.0, -0.5, 0.0, 0.0]
full_final = [1.5, 0.0, 0.0, -1.5, 0.0, 0.0]
tail_initial = [0.5 + (1.5 - 0.5) * (2.0 / 3.0), 0.0, 0.0, -0.5 + (-1.5 + 0.5) * (2.0 / 3.0), 0.0, 0.0]
equilibrated_mid = [1.0, 0.0, 0.0, -1.0, 0.0, 0.0]
estimator_fd_half_width = 1.0 / 16.0

def interp_lambda(lam):
    return [full_initial[i] + lam * (full_final[i] - full_initial[i]) for i in range(len(full_initial))]

fd_quarter_initial = interp_lambda(0.25 - estimator_fd_half_width)
fd_quarter_final = interp_lambda(0.25 + estimator_fd_half_width)
fd_threequarter_initial = interp_lambda(0.75 - estimator_fd_half_width)
fd_threequarter_final = interp_lambda(0.75 + estimator_fd_half_width)

def run_case(initial, final, nlambda=4, nmdsteps=400, neqsteps=40):
    tmp = pathlib.Path(tempfile.mkdtemp(prefix=f"ti_{mode}_{nsys}_"))
    os.chdir(tmp)
    for name in ("common_resources", "data"):
        (tmp / name).symlink_to(base, target_is_directory=True)
    mm.setVerbosity(-1, 0)
    mm.init(
        nSys_=nsys,
        xyz_name=str(xyz),
        sElementTypes=str(base / "ElementTypes.dat"),
        sAtomTypes=str(base / "AtomTypes.dat"),
        sBondTypes=str(base / "BondTypes.dat"),
        sAngleTypes=str(base / "AngleTypes.dat"),
        bMMFF=True,
        bEpairs=False,
        T=temperature,
        gamma=1.0 / (200.0 * 0.05),
    )
    result = mm.computeFreeEnergy(
        2,
        initial,
        final,
        nLambda=nlambda,
        nMDsteps=nmdsteps,
        nEQsteps=neqsteps,
        Fconv=1e-6,
        mode=0,
        K=5.0,
    )
    rows = []
    for line in (tmp / "entropic_spring_20_free_energy.dat").read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            rows.append([float(x) for x in line.split()])
    print("JSON_RESULT=" + json.dumps({"free_energy": result, "rows": rows}))


def _hardcoded_bond_forces(positions, neighs, l0=1.198, k=80.0):
    ref = np.zeros_like(positions, dtype=np.float32)
    natoms = positions.shape[0]
    for i in range(natoms):
        for j in neighs[i]:
            j = int(j)
            if j < 0 or j >= natoms:
                continue
            if i < j:
                d = positions[j] - positions[i]
                l = np.linalg.norm(d)
                f = d / l * ((l - l0) * k)
                ref[i] += f
                ref[j] -= f
    return ref


def _hardcoded_bond_energy(positions, neighs, l0=1.198, k=80.0):
    energy = 0.0
    natoms = positions.shape[0]
    for i in range(natoms):
        for j in neighs[i]:
            j = int(j)
            if j < 0 or j >= natoms:
                continue
            if i < j:
                d = np.asarray(positions[j], dtype=np.float64) - np.asarray(positions[i], dtype=np.float64)
                l = math.sqrt(float(np.dot(d, d)))
                dl = l - l0
                energy += 0.5 * k * dl * dl
    return float(energy)


def _limit_force_vectors(forces, flimit=FLIMIT):
    out = np.array(forces, copy=True)
    norms = np.linalg.norm(out, axis=1)
    mask = norms > flimit
    out[mask] *= (flimit / norms[mask])[:, None]
    return out


def run_static_force_audit():
    tmp = pathlib.Path(tempfile.mkdtemp(prefix=f"force_{mode}_{nsys}_"))
    os.chdir(tmp)
    for name in ("common_resources", "data"):
        (tmp / name).symlink_to(base, target_is_directory=True)
    mm.setVerbosity(-1, 0)
    mm.init(
        nSys_=1,
        xyz_name=str(xyz),
        sElementTypes=str(base / "ElementTypes.dat"),
        sAtomTypes=str(base / "AtomTypes.dat"),
        sBondTypes=str(base / "BondTypes.dat"),
        sAngleTypes=str(base / "AngleTypes.dat"),
        bMMFF=True,
        bEpairs=False,
        nPBC=(0, 0, 0),
        gridnPBC=(0, 0, 0),
        T=-1.0,
        gamma=-1.0,
        nExplore=0,
        nRelax=0,
    )
    mm.getBuffs()

    natoms = int(mm.natoms)
    nnode = int(mm.nnode)
    positions = np.array(mm.gpu_atoms[0, :natoms, :3], copy=True)
    neighs = np.array(mm.gpu_neighs[0, :natoms, :], copy=True)

    center = positions.mean(axis=0)
    positions[:, 0] = center[0] + 1.12 * (positions[:, 0] - center[0])
    positions[:, 1] = center[1] + 0.05 * np.sin(np.linspace(0.0, np.pi, natoms, dtype=np.float32))

    mm.gpu_atoms[0, :natoms, :3] = positions
    mm.gpu_avel[0, :, :] = 0.0
    mm.gpu_aforces[0, :, :] = 0.0
    mm.upload(bParams=False, bForces=True, bVel=True)

    ref = _hardcoded_bond_forces(positions, neighs)

    mm.eval_getMMFFf4_ocl()
    gpu1 = np.array(mm.gpu_aforces[0, :natoms, :3], copy=True)
    pi1 = np.array(mm.gpu_aforces[0, natoms:natoms + nnode, :3], copy=True)

    mm.gpu_aforces[0, :natoms, 0] = np.linspace(1.0, 2.0, natoms, dtype=np.float32)
    mm.gpu_aforces[0, :natoms, 1] = np.linspace(-0.5, 0.5, natoms, dtype=np.float32)
    mm.gpu_aforces[0, :natoms, 2] = 3.0
    mm.upload(bParams=False, bForces=True, bVel=False)

    mm.eval_getMMFFf4_ocl()
    gpu2 = np.array(mm.gpu_aforces[0, :natoms, :3], copy=True)

    print(
        "JSON_RESULT="
        + json.dumps(
            {
                "max_abs_ref_gpu1": float(np.max(np.abs(ref - gpu1))),
                "max_abs_gpu1_gpu2": float(np.max(np.abs(gpu1 - gpu2))),
                "max_pi_force": float(np.max(np.abs(pi1))) if pi1.size else 0.0,
                "max_force_norm": float(np.max(np.linalg.norm(gpu1, axis=1))),
                "sum_force": [float(x) for x in gpu1.sum(axis=0)],
                "first_atom_ref": [float(x) for x in ref[0]],
                "first_atom_gpu": [float(x) for x in gpu1[0]],
            }
        )
    )


def run_dynamic_force_audit():
    tmp = pathlib.Path(tempfile.mkdtemp(prefix=f"force_dyn_{mode}_{nsys}_"))
    os.chdir(tmp)
    for name in ("common_resources", "data"):
        (tmp / name).symlink_to(base, target_is_directory=True)
    mm.setVerbosity(-1, 0)
    mm.init(
        nSys_=1,
        xyz_name=str(xyz),
        sElementTypes=str(base / "ElementTypes.dat"),
        sAtomTypes=str(base / "AtomTypes.dat"),
        sBondTypes=str(base / "BondTypes.dat"),
        sAngleTypes=str(base / "AngleTypes.dat"),
        bMMFF=True,
        bEpairs=False,
        nPBC=(0, 0, 0),
        gridnPBC=(0, 0, 0),
        T=-1.0,
        gamma=-1.0,
        nExplore=0,
        nRelax=0,
    )
    mm.getBuffs()

    mm.MDloop(perframe=1, Ftol=0.0, iParalel=3, perVF=1)
    mm.download(bForces=True, bVel=True)

    natoms = int(mm.natoms)
    positions0 = np.array(mm.gpu_atoms[0, :natoms, :3], copy=True)
    neighs = np.array(mm.gpu_neighs[0, :natoms, :], copy=True)
    center = positions0.mean(axis=0)

    positions = np.array(positions0, copy=True)
    positions[:, 0] = center[0] + 1.12 * (positions[:, 0] - center[0])
    positions[:, 1] = center[1] + 0.05 * np.sin(np.linspace(0.0, np.pi, natoms, dtype=np.float32))
    velocities0 = np.zeros_like(positions, dtype=np.float32)

    ref_force = _limit_force_vectors(_hardcoded_bond_forces(positions, neighs))
    ref_vel1 = velocities0 + ref_force * DT
    ref_pos1 = positions + ref_vel1 * DT

    mm.gpu_atoms[0, :natoms, :3] = positions
    mm.gpu_avel[0, :, :] = 0.0
    mm.upload(bParams=False, bForces=False, bVel=True)

    mm.MDloop(perframe=1, Ftol=0.0, iParalel=3, perVF=1)
    mm.download(bForces=True, bVel=True)

    pos_gpu = np.array(mm.gpu_atoms[0, :natoms, :3], copy=True)
    vel_gpu = np.array(mm.gpu_avel[0, :natoms, :3], copy=True)
    eff_force = vel_gpu / DT

    print(
        "JSON_RESULT="
        + json.dumps(
            {
                "max_force_err_from_velocity": float(np.max(np.abs(eff_force - ref_force))),
                "max_pos_err": float(np.max(np.abs(pos_gpu - ref_pos1))),
                "max_vel_err": float(np.max(np.abs(vel_gpu - ref_vel1))),
                "max_force_norm": float(np.max(np.linalg.norm(ref_force, axis=1))),
                "first_atom_ref_force": [float(x) for x in ref_force[0]],
                "first_atom_gpu_velocity": [float(x) for x in vel_gpu[0]],
                "first_atom_gpu_position": [float(x) for x in pos_gpu[0]],
            }
        )
    )


def _init_mm_for_ti(run_nsys, run_temperature):
    tmp = pathlib.Path(tempfile.mkdtemp(prefix=f"ti_whitebox_{mode}_{run_nsys}_"))
    os.chdir(tmp)
    for name in ("common_resources", "data"):
        (tmp / name).symlink_to(base, target_is_directory=True)
    mm.setVerbosity(-1, 0)
    mm.init(
        nSys_=run_nsys,
        xyz_name=str(xyz),
        sElementTypes=str(base / "ElementTypes.dat"),
        sAtomTypes=str(base / "AtomTypes.dat"),
        sBondTypes=str(base / "BondTypes.dat"),
        sAngleTypes=str(base / "AngleTypes.dat"),
        bMMFF=True,
        bEpairs=False,
        T=run_temperature,
        gamma=1.0 / (200.0 * 0.05),
    )
    mm.getBuffs()
    return tmp


def _run_single_lambda_estimator(batch, nlambda=5, nmdsteps=20000, neqsteps=2000, prod_blocks=None):
    _init_mm_for_ti(1, temperature)
    nprod = nmdsteps // nlambda
    blocks = prod_blocks if prod_blocks is not None else [nprod]
    mm.setupTIConstraintsBatch(2, full_initial, full_final, nlambda, batch)
    mm.MDloop(perframe=neqsteps, Ftol=1e-6, iParalel=3, perVF=100)
    mm.zeroTIAverageForces()
    for block in blocks:
        mm.MDloop(perframe=block, Ftol=1e-6, iParalel=3, perVF=100)
    mm.downloadTIAverageForces()
    est = mm.getTIHostEstimator(0, nlambda, nmdsteps)
    avg = np.array(mm.gpu_averageForces[0], copy=True)
    print(
        "JSON_RESULT="
        + json.dumps(
            {
                "average_forces": [float(x) for x in avg],
                "force_sum": float(est[0]),
                "dE": float(est[1]),
                "mean_sq": float(est[2]),
                "sem": float(est[3]),
                "scale": float(est[4]),
                "nProdSteps": int(round(est[5])),
            }
        )
    )


def _run_batch_boundary_state_probe(target_batch=1, nlambda=5, nmdsteps=20000, neqsteps=2000):
    _init_mm_for_ti(1, temperature)
    natoms = int(mm.natoms)
    nprod = nmdsteps // nlambda

    mm.resetTIBatchState()
    mm.setupTIConstraintsBatch(2, full_initial, full_final, nlambda, target_batch)
    mm.download_sys(0, bForces=False, bVel=True)
    fresh_pos = np.array(mm.gpu_atoms[0, :natoms, :3], copy=True)
    fresh_vel = np.array(mm.gpu_avel[0, :natoms, :3], copy=True)

    mm.setupTIConstraintsBatch(2, full_initial, full_final, nlambda, 0)
    mm.MDloop(perframe=neqsteps, Ftol=1e-6, iParalel=3, perVF=100)
    mm.zeroTIAverageForces()
    mm.MDloop(perframe=nprod, Ftol=1e-6, iParalel=3, perVF=100)

    mm.resetTIBatchState()
    mm.setupTIConstraintsBatch(2, full_initial, full_final, nlambda, target_batch)
    mm.download_sys(0, bForces=False, bVel=True)
    cont_pos = np.array(mm.gpu_atoms[0, :natoms, :3], copy=True)
    cont_vel = np.array(mm.gpu_avel[0, :natoms, :3], copy=True)

    print(
        "JSON_RESULT="
        + json.dumps(
            {
                "max_pos_diff": float(np.max(np.abs(cont_pos - fresh_pos))),
                "max_vel_diff": float(np.max(np.abs(cont_vel - fresh_vel))),
                "fresh_first_atom": [float(x) for x in fresh_pos[0]],
                "continued_first_atom": [float(x) for x in cont_pos[0]],
                "fresh_first_vel": [float(x) for x in fresh_vel[0]],
                "continued_first_vel": [float(x) for x in cont_vel[0]],
            }
        )
    )


def _run_nondivisible_nmdsteps_guard_probe():
    tmp = _init_mm_for_ti(1, temperature)
    result = mm.computeFreeEnergy(
        2,
        full_initial,
        full_final,
        nLambda=5,
        nMDsteps=20003,
        nEQsteps=2000,
        Fconv=1e-6,
        mode=0,
        K=5.0,
    )
    print(
        "JSON_RESULT="
        + json.dumps(
            {
                "free_energy": None if math.isnan(result) else float(result),
                "is_nan": bool(math.isnan(result)),
                "file_exists": bool((tmp / "entropic_spring_20_free_energy.dat").exists()),
            }
        )
    )


def _run_snapshot_fd_integrand_probe(batch=2, nlambda=5, neqsteps=2000, nsamples=5, eps=1.0e-4):
    _init_mm_for_ti(1, temperature)
    natoms = int(mm.natoms)
    path_initial = np.array(full_initial, dtype=np.float64).reshape(2, 3)
    path_final = np.array(full_final, dtype=np.float64).reshape(2, 3)
    mm.resetTIBatchState()
    mm.setupTIConstraintsBatch(2, full_initial, full_final, nlambda, batch)
    mm.MDloop(perframe=neqsteps, Ftol=1e-6, iParalel=3, perVF=100)

    lam = float(batch) / float(nlambda - 1)
    target_plus = path_initial + (path_final - path_initial) * (lam + eps)
    target_minus = path_initial + (path_final - path_initial) * (lam - eps)
    neighs = np.array(mm.gpu_neighs[0, :natoms, :], copy=True)

    constrained = np.where(np.array(mm.gpu_constr[0, :natoms, 3], copy=True) > 1.0e5)[0]
    assert constrained.shape == (2,), constrained
    cks = np.array(mm.gpu_constrK[0, constrained, :3], copy=True)

    samples = []
    for isample in range(nsamples):
        mm.download_sys(0, bForces=False, bVel=True)
        snap_pos = np.array(mm.gpu_atoms[0, :natoms, :3], dtype=np.float64, copy=True)

        mm.eval_getMMFFf4_ocl()
        forces = np.array(mm.gpu_aforces[0, :natoms, :3], copy=True)
        proj0 = float(np.dot(forces[constrained[0]], cks[0]))
        proj1 = float(np.dot(forces[constrained[1]], cks[1]))
        inst_sum = proj0 + proj1
        inst_dE = -inst_sum

        plus_pos = np.array(snap_pos, copy=True)
        minus_pos = np.array(snap_pos, copy=True)
        plus_pos[constrained[0]] = target_plus[0]
        plus_pos[constrained[1]] = target_plus[1]
        minus_pos[constrained[0]] = target_minus[0]
        minus_pos[constrained[1]] = target_minus[1]

        e_plus = _hardcoded_bond_energy(plus_pos, neighs)
        e_minus = _hardcoded_bond_energy(minus_pos, neighs)

        fd_dudl = (e_plus - e_minus) / (2.0 * eps)
        print(
            f"TI_DEBUG_SNAPSHOT_FD sample={isample} lambda={lam:.17g} eps={eps:.17g} model=hardcoded_bond "
            f"Eplus={e_plus:.17g} Eminus={e_minus:.17g} fd_dU_dlambda={fd_dudl:.17g} "
            f"proj0={proj0:.17g} proj1={proj1:.17g} inst_sum={inst_sum:.17g} inst_dE={inst_dE:.17g}"
        )
        samples.append(
            {
                "sample": isample,
                "lambda": lam,
                "eps": eps,
                "model": "hardcoded_bond",
                "e_plus": e_plus,
                "e_minus": e_minus,
                "fd_dudl": fd_dudl,
                "proj0": proj0,
                "proj1": proj1,
                "inst_sum": inst_sum,
                "inst_dE": inst_dE,
                "abs_diff": abs(fd_dudl - inst_dE),
            }
        )

        if isample + 1 < nsamples:
            mm.MDloop(perframe=1, Ftol=0.0, iParalel=3, perVF=1)

    print(
        "JSON_RESULT="
        + json.dumps(
            {
                "lambda": lam,
                "eps": eps,
                "samples": samples,
                "max_abs_diff": float(max(s["abs_diff"] for s in samples)),
            }
        )
    )

try:
    if mode == "zero_path":
        run_case(full_initial, full_initial)
    elif mode == "stretch":
        run_case(full_initial, full_final)
    elif mode == "stretch_long":
        run_case(full_initial, full_final, nlambda=4, nmdsteps=12000, neqsteps=400)
    elif mode == "stretch_second_half":
        run_case(tail_initial, full_final, nlambda=2, nmdsteps=200)
    elif mode == "fd_lambda_1_3_long":
        run_case(interp_lambda(1.0 / 3.0 - 1.0 / 18.0), interp_lambda(1.0 / 3.0 + 1.0 / 18.0), nlambda=2, nmdsteps=6000, neqsteps=400)
    elif mode == "fd_lambda_2_3_long":
        run_case(interp_lambda(2.0 / 3.0 - 1.0 / 18.0), interp_lambda(2.0 / 3.0 + 1.0 / 18.0), nlambda=2, nmdsteps=6000, neqsteps=400)
    elif mode == "gpu_accounting":
        run_case(full_initial, full_final, nlambda=2, nmdsteps=2, neqsteps=0)
    elif mode == "gpu_accounting_equilibrated":
        run_case(full_initial, full_final, nlambda=2, nmdsteps=40, neqsteps=2000)
    elif mode == "stretch_equilibrated_estimator":
        run_case(full_initial, full_final, nlambda=5, nmdsteps=20000, neqsteps=2000)
    elif mode == "stretch_reverse_equilibrated_estimator":
        run_case(full_final, full_initial, nlambda=5, nmdsteps=20000, neqsteps=2000)
    elif mode == "stretch_first_half_equilibrated_estimator":
        run_case(full_initial, equilibrated_mid, nlambda=5, nmdsteps=20000, neqsteps=2000)
    elif mode == "stretch_second_half_equilibrated_estimator":
        run_case(equilibrated_mid, full_final, nlambda=5, nmdsteps=20000, neqsteps=2000)
    elif mode == "fd_lambda_1_4_equilibrated":
        run_case(fd_quarter_initial, fd_quarter_final, nlambda=2, nmdsteps=20000, neqsteps=2000)
    elif mode == "fd_lambda_3_4_equilibrated":
        run_case(fd_threequarter_initial, fd_threequarter_final, nlambda=2, nmdsteps=20000, neqsteps=2000)
    elif mode == "stretch_equilibrated_analytic":
        run_case(full_initial, full_final, nlambda=5, nmdsteps=20000, neqsteps=2000)
    elif mode == "bond_force_static":
        run_static_force_audit()
    elif mode == "bond_force_dynamic":
        run_dynamic_force_audit()
    elif mode == "ti_single_lambda_normalization":
        _run_single_lambda_estimator(batch=2)
    elif mode == "ti_single_lambda_block_single":
        _run_single_lambda_estimator(batch=2)
    elif mode == "ti_single_lambda_block_split":
        _run_single_lambda_estimator(batch=2, prod_blocks=[2000, 2000])
    elif mode == "ti_batch_boundary_state":
        _run_batch_boundary_state_probe(target_batch=1)
    elif mode == "ti_nondivisible_nmdsteps_guard":
        _run_nondivisible_nmdsteps_guard_probe()
    elif mode == "ti_snapshot_fd_integrand":
        _run_snapshot_fd_integrand_probe(batch=2)
    else:
        raise ValueError(mode)
finally:
    try:
        mm.clear()
    except Exception:
        pass
'''


def _run_case_impl(mode: str, nsys: int, temperature: float = 300.0, stdout_path: Path | None = None) -> dict:
    proc = subprocess.run(
        [sys.executable, "-c", CASE_SCRIPT, str(REPO_ROOT), mode, str(nsys), str(temperature)],
        cwd=CPP_DIR,
        capture_output=True,
        text=True,
    )
    if stdout_path is not None:
        stdout_path.write_text(proc.stdout)
    assert proc.returncode == 0, (
        f"computeFreeEnergy subprocess failed for {mode} nSys={nsys}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )
    for line in reversed(proc.stdout.splitlines()):
        if line.startswith("JSON_RESULT="):
            return json.loads(line.split("=", 1)[1])
    raise AssertionError(f"JSON result missing for {mode} nSys={nsys}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")


@lru_cache(maxsize=None)
def _run_case_cached(mode: str, nsys: int, temperature: float = 300.0) -> dict:
    return _run_case_impl(mode, nsys, temperature, stdout_path=None)


def _run_case(mode: str, nsys: int, temperature: float = 300.0, stdout_path: Path | None = None) -> dict:
    if stdout_path is not None:
        return _run_case_impl(mode, nsys, temperature, stdout_path=stdout_path)
    return deepcopy(_run_case_cached(mode, nsys, temperature))


def _ti_step_series(samples):
    cumulative_fe = 0.0
    cumulative_var = 0.0
    prev_dE = 0.0
    prev_sigma = 0.0
    out = []
    for lam, dE, sigma in samples:
        if lam > 1e-9:
            d_lambda = lam - out[-1][0]
            cumulative_fe += 0.5 * (dE + prev_dE) * d_lambda
            cumulative_var += 0.25 * d_lambda * d_lambda * (sigma * sigma + prev_sigma * prev_sigma)
        prev_dE = dE
        prev_sigma = sigma
        out.append((lam, cumulative_fe, math.sqrt(cumulative_var)))
    return out


@lru_cache(maxsize=None)
def _bond_only_gaussian_chain_mean_l2(
    temperature: float = 300.0,
    bond_length: float = 1.198,
    bond_k: float = 80.0,
) -> float:
    k_b = 8.617333262145e-5
    if temperature <= 0.0:
        return bond_length * bond_length

    kbt = k_b * temperature
    sigma = math.sqrt(kbt / bond_k)
    l_max = max(4.0 * bond_length, bond_length + 40.0 * sigma)
    grid = np.linspace(0.0, l_max, 20001, dtype=np.float64)
    boltz = np.exp(-0.5 * bond_k * (grid - bond_length) ** 2 / kbt)
    denom = np.trapz((grid**2) * boltz, grid)
    numer = np.trapz((grid**4) * boltz, grid)
    return float(numer / denom)


def _bond_only_gaussian_chain_reference_force_constant(
    temperature: float = 300.0,
    n_segments: int = 20 - 1,
    bond_length: float = 1.198,
    bond_k: float = 80.0,
) -> float:
    k_b = 8.617333262145e-5
    mean_l2 = _bond_only_gaussian_chain_mean_l2(temperature, bond_length, bond_k)
    return 3.0 * k_b * temperature / (n_segments * mean_l2)


def _bond_only_gaussian_chain_reference_free_energy(distance: float, temperature: float = 300.0) -> float:
    k_eff = _bond_only_gaussian_chain_reference_force_constant(temperature=temperature)
    return 0.5 * k_eff * distance * distance


def _bond_only_gaussian_chain_reference_dedlambda(distance: float, dr_dlambda: float, temperature: float = 300.0) -> float:
    k_eff = _bond_only_gaussian_chain_reference_force_constant(temperature=temperature)
    return k_eff * distance * dr_dlambda


def _rows_by_lambda(rows: list[list[float]]) -> dict[float, list[float]]:
    return {row[0]: row for row in rows}


def _parse_gpu_accounting_log(text: str) -> dict:
    batch_re = re.compile(r"--- TI Batch ([0-9]+)/([0-9]+) ---")
    equil_re = re.compile(r"Equilibrating ([0-9]+) steps\.\.\.")
    prod_re = re.compile(r"Production ([0-9]+) steps\.\.\.")
    gpu_re = re.compile(
        r"TI_DEBUG_GPU iS=([0-9]+) ia=([0-9]+) sign=([0-9eE+\-.]+) force_proj=([0-9eE+\-.]+) "
        r"fe=\(([0-9eE+\-.]+),([0-9eE+\-.]+),([0-9eE+\-.]+)\) cK=\(([0-9eE+\-.]+),([0-9eE+\-.]+),([0-9eE+\-.]+)\)"
    )
    host_avg_re = re.compile(
        r"TI_DEBUG_HOST_AVG batch=([0-9]+) isys=([0-9]+) lambda=([0-9eE+\-.]+) "
        r"avg=\(([0-9eE+\-.]+),([0-9eE+\-.]+),([0-9eE+\-.]+),([0-9eE+\-.]+)\)"
    )
    host_de_re = re.compile(
        r"TI_DEBUG_HOST_DE batch=([0-9]+) isys=([0-9]+) lambda=([0-9eE+\-.]+) "
        r"force_sum=([0-9eE+\-.]+) scale=([0-9eE+\-.]+) dE=([0-9eE+\-.]+)"
    )
    host_stats_re = re.compile(
        r"TI_DEBUG_HOST_STATS batch=([0-9]+) isys=([0-9]+) lambda=([0-9eE+\-.]+) "
        r"nprod=([0-9]+) mean_sq=([0-9eE+\-.]+) var=([0-9eE+\-.]+) sd=([0-9eE+\-.]+) sem=([0-9eE+\-.]+)"
    )
    ti_step_re = re.compile(
        r"TI_DEBUG_TI_STEP batch=([0-9]+) isys=([0-9]+) lambda=([0-9eE+\-.]+) dLambda=([0-9eE+\-.]+) "
        r"prev_dE=([0-9eE+\-.]+) prev_sigma=([0-9eE+\-.]+) in_dE=([0-9eE+\-.]+) in_sigma=([0-9eE+\-.]+) "
        r"cumFE_before=([0-9eE+\-.]+) cumVar_before=([0-9eE+\-.]+) outF=([0-9eE+\-.]+) outFerr=([0-9eE+\-.]+) "
        r"cumFE_after=([0-9eE+\-.]+) cumVar_after=([0-9eE+\-.]+)"
    )
    current_batch = None
    current_phase = None
    out = {"gpu": {}, "host_avg": {}, "host_de": {}, "host_stats": {}, "ti_step": {}}
    for line in text.splitlines():
        m = batch_re.search(line)
        if m:
            current_batch = int(m.group(1))
            current_phase = None
            continue
        m = equil_re.search(line)
        if m:
            current_phase = "equilibration"
            continue
        m = prod_re.search(line)
        if m:
            current_phase = "production"
            continue
        m = gpu_re.search(line)
        if m:
            assert current_batch is not None, f"GPU debug line without current batch: {line}"
            out["gpu"].setdefault(current_batch, []).append(
                {
                    "isys": int(m.group(1)),
                    "ia": int(m.group(2)),
                    "sign": float(m.group(3)),
                    "force_proj": float(m.group(4)),
                    "fe": tuple(float(m.group(i)) for i in range(5, 8)),
                    "cK": tuple(float(m.group(i)) for i in range(8, 11)),
                    "phase": current_phase,
                }
            )
            continue
        m = host_avg_re.search(line)
        if m:
            out["host_avg"][int(m.group(1))] = {
                "isys": int(m.group(2)),
                "lambda": float(m.group(3)),
                "x": float(m.group(4)),
                "y": float(m.group(5)),
                "z": float(m.group(6)),
                "w": float(m.group(7)),
            }
            continue
        m = host_de_re.search(line)
        if m:
            out["host_de"][int(m.group(1))] = {
                "isys": int(m.group(2)),
                "lambda": float(m.group(3)),
                "force_sum": float(m.group(4)),
                "scale": float(m.group(5)),
                "dE": float(m.group(6)),
            }
            continue
        m = host_stats_re.search(line)
        if m:
            out["host_stats"][int(m.group(1))] = {
                "isys": int(m.group(2)),
                "lambda": float(m.group(3)),
                "nprod": int(m.group(4)),
                "mean_sq": float(m.group(5)),
                "var": float(m.group(6)),
                "sd": float(m.group(7)),
                "sem": float(m.group(8)),
            }
            continue
        m = ti_step_re.search(line)
        if m:
            out["ti_step"][int(m.group(1))] = {
                "isys": int(m.group(2)),
                "lambda": float(m.group(3)),
                "dLambda": float(m.group(4)),
                "prev_dE": float(m.group(5)),
                "prev_sigma": float(m.group(6)),
                "in_dE": float(m.group(7)),
                "in_sigma": float(m.group(8)),
                "cumFE_before": float(m.group(9)),
                "cumVar_before": float(m.group(10)),
                "outF": float(m.group(11)),
                "outFerr": float(m.group(12)),
                "cumFE_after": float(m.group(13)),
                "cumVar_after": float(m.group(14)),
            }
    return out


def test_TI_step_python_reference_matches_trapezoid_and_error_propagation():
    series = _ti_step_series([(0.0, 1.0, 0.2), (0.5, 3.0, 0.4), (1.0, 5.0, 0.6)])
    assert series[0] == (0.0, 0.0, 0.0)
    assert math.isclose(series[1][1], 1.0, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(series[2][1], 3.0, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(series[1][2], math.sqrt(0.0125), rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(series[2][2], math.sqrt(0.045), rel_tol=0.0, abs_tol=1e-15)


def test_bond_only_gaussian_chain_reference_curve_is_self_consistent_for_entropic_spring20():
    mean_l2 = _bond_only_gaussian_chain_mean_l2(temperature=300.0)
    assert math.isclose(mean_l2, 1.198**2, rel_tol=2.0e-3, abs_tol=0.0)

    eps = 1.0e-6
    expected_refs = {0.25: 0.00852277075773617, 0.5: 0.011363694343648227, 0.75: 0.014204617929560284}
    for lam, distance in ((0.25, 1.5), (0.5, 2.0), (0.75, 2.5)):
        fd = (
            _bond_only_gaussian_chain_reference_free_energy(distance + eps, temperature=300.0)
            - _bond_only_gaussian_chain_reference_free_energy(distance - eps, temperature=300.0)
        ) / (2.0 * eps)
        dedl = _bond_only_gaussian_chain_reference_dedlambda(distance, dr_dlambda=2.0, temperature=300.0)
        assert math.isclose(dedl, expected_refs[lam], rel_tol=0.0, abs_tol=5.0e-15)
        assert math.isclose(dedl, 2.0 * fd, rel_tol=0.0, abs_tol=5e-12)


def test_computeFreeEnergy_TI_zero_path_on_entropic_spring20_is_exactly_zero():
    result = _run_case("zero_path", nsys=4)
    rows = result["rows"]
    assert len(rows) == 4
    assert abs(result["free_energy"]) < 1e-12
    assert max(abs(row[1]) for row in rows) < 1e-12
    assert max(abs(row[3]) for row in rows) < 1e-12
    assert max(abs(row[10] - 1.0) for row in rows) < 1e-12


@pytest.mark.parametrize("nsys", [2, 4])
def test_computeFreeEnergy_TI_stretch_is_reproducible_for_fixed_nSys(nsys):
    result1 = _run_case("stretch", nsys=nsys)
    result2 = _run_case("stretch", nsys=nsys)
    assert result1 == result2


def test_computeFreeEnergy_TI_zero_path_is_batching_invariant_across_nSys():
    result4 = _run_case("zero_path", nsys=4)
    result2 = _run_case("zero_path", nsys=2)
    assert result4 == result2


def test_computeFreeEnergy_TI_long_T0_dEdlambda_matches_local_finite_difference_reference():
    full = _run_case("stretch_long", nsys=4, temperature=0.0)
    fd13 = _run_case("fd_lambda_1_3_long", nsys=2, temperature=0.0)
    fd23 = _run_case("fd_lambda_2_3_long", nsys=2, temperature=0.0)
    half_width = 1.0 / 18.0
    dE = [row[1] for row in full["rows"]]
    fd_ref_13 = fd13["free_energy"] / (2.0 * half_width)
    fd_ref_23 = fd23["free_energy"] / (2.0 * half_width)
    assert all(math.isfinite(v) and v > 0.0 for v in dE)
    assert math.isclose(dE[1], fd_ref_13, rel_tol=0.35, abs_tol=0.07)
    assert math.isclose(dE[2], fd_ref_23, rel_tol=0.35, abs_tol=0.07)


def test_computeFreeEnergy_TI_gpu_accounting_matches_host_accumulator_and_ti_input(tmp_path: Path):
    log_path = tmp_path / "gpu_accounting.log"
    result = _run_case("gpu_accounting", nsys=1, temperature=0.0, stdout_path=log_path)
    parsed = _parse_gpu_accounting_log(log_path.read_text())
    assert set(parsed["gpu"]) == {1, 2}
    assert set(parsed["host_avg"]) == {1, 2}
    assert set(parsed["host_de"]) == {1, 2}

    for batch in (1, 2):
        gpu_entries = parsed["gpu"][batch]
        assert len(gpu_entries) == 2
        pos = next(entry for entry in gpu_entries if entry["sign"] > 0.0)
        neg = next(entry for entry in gpu_entries if entry["sign"] < 0.0)
        for entry in gpu_entries:
            proj_from_dot = sum(f * c for f, c in zip(entry["fe"], entry["cK"]))
            assert math.isclose(entry["force_proj"], proj_from_dot, rel_tol=1e-6, abs_tol=1e-6)

        host_avg = parsed["host_avg"][batch]
        host_de = parsed["host_de"][batch]
        assert math.isclose(host_avg["z"], pos["force_proj"], rel_tol=1e-6, abs_tol=1e-6)
        assert math.isclose(host_avg["w"], neg["force_proj"], rel_tol=1e-6, abs_tol=1e-6)

        force_sum = pos["force_proj"] + neg["force_proj"]
        assert math.isclose(host_de["force_sum"], force_sum, rel_tol=1e-6, abs_tol=1e-6)
        assert math.isclose(host_de["dE"], -host_de["force_sum"] * host_de["scale"], rel_tol=1e-12, abs_tol=1e-12)
        assert math.isclose(result["rows"][batch - 1][1], host_de["dE"], rel_tol=1e-7, abs_tol=1e-7)


def test_computeFreeEnergy_TI_equilibrated_gpu_accounting_accumulates_the_correct_instantaneous_forces(tmp_path: Path):
    log_path = tmp_path / "gpu_accounting_equilibrated.log"
    result = _run_case("gpu_accounting_equilibrated", nsys=1, temperature=300.0, stdout_path=log_path)
    parsed = _parse_gpu_accounting_log(log_path.read_text())
    assert set(parsed["gpu"]) == {1, 2}
    assert set(parsed["host_avg"]) == {1, 2}
    assert set(parsed["host_de"]) == {1, 2}

    for batch in (1, 2):
        gpu_entries = parsed["gpu"][batch]
        prod_entries = [entry for entry in gpu_entries if entry["phase"] == "production"]
        equil_entries = [entry for entry in gpu_entries if entry["phase"] == "equilibration"]
        assert len(prod_entries) == 40
        assert len(equil_entries) >= 4000

        pos_entries = [entry for entry in prod_entries if entry["sign"] > 0.0]
        neg_entries = [entry for entry in prod_entries if entry["sign"] < 0.0]
        assert len(pos_entries) == 20
        assert len(neg_entries) == 20

        for entry in prod_entries:
            proj_from_dot = sum(f * c for f, c in zip(entry["fe"], entry["cK"]))
            assert math.isclose(entry["force_proj"], proj_from_dot, rel_tol=1e-6, abs_tol=1e-6)

        pos_sum = sum(entry["force_proj"] for entry in pos_entries)
        neg_sum = sum(entry["force_proj"] for entry in neg_entries)
        host_avg = parsed["host_avg"][batch]
        host_de = parsed["host_de"][batch]

        # TI_DEBUG_GPU prints with limited precision, so summed per-step projections can differ
        # slightly from the host accumulator / host_de values printed with higher precision.
        assert math.isclose(host_avg["z"], pos_sum, rel_tol=1e-6, abs_tol=5e-5)
        assert math.isclose(host_avg["w"], neg_sum, rel_tol=1e-6, abs_tol=5e-5)

        force_sum = pos_sum + neg_sum
        assert math.isclose(host_de["force_sum"], force_sum, rel_tol=1e-6, abs_tol=5e-5)
        assert math.isclose(host_de["dE"], -host_de["force_sum"] * host_de["scale"], rel_tol=1e-12, abs_tol=1e-12)
        assert math.isclose(result["rows"][batch - 1][1], host_de["dE"], rel_tol=1e-7, abs_tol=1e-7)


def test_computeFreeEnergy_TI_equilibrated_host_estimator_log_matches_averageForces_conversion_and_ti_step(tmp_path: Path):
    log_path = tmp_path / "gpu_accounting_equilibrated_host_estimator.log"
    result = _run_case("gpu_accounting_equilibrated", nsys=1, temperature=300.0, stdout_path=log_path)
    parsed = _parse_gpu_accounting_log(log_path.read_text())
    assert set(parsed["host_avg"]) == {1, 2}
    assert set(parsed["host_de"]) == {1, 2}
    assert set(parsed["host_stats"]) == {1, 2}
    assert set(parsed["ti_step"]) == {1, 2}

    rows = result["rows"]
    prev_dE = 0.0
    prev_sem = 0.0
    cumulative_fe = 0.0
    cumulative_var = 0.0

    for batch in (1, 2):
        avg = parsed["host_avg"][batch]
        de = parsed["host_de"][batch]
        stats = parsed["host_stats"][batch]
        step = parsed["ti_step"][batch]

        force_sum_f32 = float(np.float32(avg["z"]) + np.float32(avg["w"]))
        assert math.isclose(de["force_sum"], force_sum_f32, rel_tol=0.0, abs_tol=1e-12)
        assert stats["nprod"] == 20
        assert math.isclose(de["scale"], 1.0 / stats["nprod"], rel_tol=0.0, abs_tol=1e-15)
        assert math.isclose(de["dE"], -de["force_sum"] * de["scale"], rel_tol=0.0, abs_tol=1e-12)

        mean_sq = avg["x"] * de["scale"]
        var = mean_sq - de["dE"] * de["dE"]
        sd = math.sqrt(abs(var))
        sem = sd / math.sqrt(stats["nprod"])

        assert math.isclose(stats["mean_sq"], mean_sq, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(stats["var"], var, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(stats["sd"], sd, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(stats["sem"], sem, rel_tol=0.0, abs_tol=1e-12)

        lam = de["lambda"]
        d_lambda = 1.0
        assert math.isclose(step["lambda"], lam, rel_tol=0.0, abs_tol=1e-15)
        assert math.isclose(step["dLambda"], d_lambda, rel_tol=0.0, abs_tol=1e-15)
        assert math.isclose(step["prev_dE"], prev_dE, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(step["prev_sigma"], prev_sem, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(step["in_dE"], de["dE"], rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(step["in_sigma"], stats["sem"], rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(step["cumFE_before"], cumulative_fe, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(step["cumVar_before"], cumulative_var, rel_tol=0.0, abs_tol=1e-12)

        if lam > 1e-9:
            cumulative_fe += 0.5 * (de["dE"] + prev_dE) * d_lambda
            cumulative_var += 0.25 * d_lambda * d_lambda * (stats["sem"] * stats["sem"] + prev_sem * prev_sem)

        assert math.isclose(step["outF"], cumulative_fe, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(step["outFerr"], math.sqrt(cumulative_var), rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(step["cumFE_after"], cumulative_fe, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(step["cumVar_after"], cumulative_var, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(rows[batch - 1][1], de["dE"], rel_tol=0.0, abs_tol=1e-7)
        assert math.isclose(rows[batch - 1][3], step["outF"], rel_tol=0.0, abs_tol=1e-7)
        assert math.isclose(rows[batch - 1][4], step["outFerr"], rel_tol=0.0, abs_tol=1e-7)

        prev_dE = de["dE"]
        prev_sem = stats["sem"]


def test_TI_snapshot_fd_matches_instantaneous_gpu_integrand_after_equilibration(tmp_path: Path):
    log_path = tmp_path / "ti_snapshot_fd_integrand.log"
    result = _run_case("ti_snapshot_fd_integrand", nsys=1, temperature=300.0, stdout_path=log_path)
    lines = log_path.read_text().splitlines()
    fd_re = re.compile(
        r"TI_DEBUG_SNAPSHOT_FD sample=([0-9]+) lambda=([0-9eE+\-.]+) eps=([0-9eE+\-.]+) model=([A-Za-z0-9_]+) "
        r"Eplus=([0-9eE+\-.]+) Eminus=([0-9eE+\-.]+) fd_dU_dlambda=([0-9eE+\-.]+) "
        r"proj0=([0-9eE+\-.]+) proj1=([0-9eE+\-.]+) inst_sum=([0-9eE+\-.]+) inst_dE=([0-9eE+\-.]+)"
    )
    parsed = []
    for line in lines:
        m = fd_re.search(line)
        if not m:
            continue
        parsed.append(
            {
                "sample": int(m.group(1)),
                "lambda": float(m.group(2)),
                "eps": float(m.group(3)),
                "model": m.group(4),
                "e_plus": float(m.group(5)),
                "e_minus": float(m.group(6)),
                "fd_dudl": float(m.group(7)),
                "proj0": float(m.group(8)),
                "proj1": float(m.group(9)),
                "inst_sum": float(m.group(10)),
                "inst_dE": float(m.group(11)),
            }
        )

    assert len(parsed) == len(result["samples"]) == 5
    for got, js in zip(parsed, result["samples"]):
        assert got["sample"] == js["sample"]
        assert math.isclose(got["lambda"], js["lambda"], rel_tol=0.0, abs_tol=1e-15)
        assert math.isclose(got["eps"], js["eps"], rel_tol=0.0, abs_tol=1e-15)
        assert got["model"] == js["model"] == "hardcoded_bond"
        assert math.isclose(got["proj0"] + got["proj1"], got["inst_sum"], rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(got["inst_dE"], -got["inst_sum"], rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(got["fd_dudl"], js["fd_dudl"], rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(got["inst_dE"], js["inst_dE"], rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(got["fd_dudl"], got["inst_dE"], rel_tol=2e-4, abs_tol=5e-4)

    assert result["max_abs_diff"] < 5e-4


def test_computeFreeEnergy_TI_equilibrated_reverse_path_flips_free_energy_and_interior_dEdlambda_sign():
    forward = _run_case("stretch_equilibrated_estimator", nsys=5, temperature=300.0)
    reverse = _run_case("stretch_reverse_equilibrated_estimator", nsys=5, temperature=300.0)
    assert math.isclose(forward["free_energy"], -reverse["free_energy"], rel_tol=0.08, abs_tol=0.01)

    forward_by_lambda = _rows_by_lambda(forward["rows"])
    reverse_by_lambda = _rows_by_lambda(reverse["rows"])
    for lam in (0.25, 0.5, 0.75):
        assert math.isclose(
            forward_by_lambda[lam][1],
            -reverse_by_lambda[1.0 - lam][1],
            rel_tol=0.08,
            abs_tol=0.01,
        )


def test_computeFreeEnergy_TI_equilibrated_is_additive_over_path_segments():
    full = _run_case("stretch_equilibrated_estimator", nsys=5, temperature=300.0)
    first_half = _run_case("stretch_first_half_equilibrated_estimator", nsys=5, temperature=300.0)
    second_half = _run_case("stretch_second_half_equilibrated_estimator", nsys=5, temperature=300.0)

    midpoint_cumulative_fe = _rows_by_lambda(full["rows"])[0.5][3]
    assert math.isclose(midpoint_cumulative_fe, first_half["free_energy"], rel_tol=0.08, abs_tol=0.01)
    assert math.isclose(
        full["free_energy"],
        first_half["free_energy"] + second_half["free_energy"],
        rel_tol=0.08,
        abs_tol=0.01,
    )


@pytest.mark.xfail(reason="Equilibrated TI dE/dlambda does not match local finite-difference free-energy slopes", strict=True)
def test_computeFreeEnergy_TI_equilibrated_local_finite_difference_matches_reported_dEdlambda():
    full = _run_case("stretch_equilibrated_estimator", nsys=5, temperature=300.0)
    fd_quarter = _run_case("fd_lambda_1_4_equilibrated", nsys=2, temperature=300.0)
    fd_threequarter = _run_case("fd_lambda_3_4_equilibrated", nsys=2, temperature=300.0)
    half_width = 1.0 / 16.0
    by_lambda = _rows_by_lambda(full["rows"])
    quarter_ref = fd_quarter["free_energy"] / (2.0 * half_width)
    threequarter_ref = fd_threequarter["free_energy"] / (2.0 * half_width)

    print(f"\n[lambda=0.25] Reported dE/dl: {by_lambda[0.25][1]:.6f}, Expected: {quarter_ref:.6f}, Diff: {abs(by_lambda[0.25][1] - quarter_ref):.6f}")
    print(f"[lambda=0.75] Reported dE/dl: {by_lambda[0.75][1]:.6f}, Expected: {threequarter_ref:.6f}, Diff: {abs(by_lambda[0.75][1] - threequarter_ref):.6f}")

    assert math.isclose(by_lambda[0.25][1], quarter_ref, rel_tol=0.12, abs_tol=0.015)
    assert math.isclose(by_lambda[0.75][1], threequarter_ref, rel_tol=0.12, abs_tol=0.015)


def test_computeFreeEnergy_TI_equilibrated_profile_is_invariant_to_batching():
    single_batch = _run_case("stretch_equilibrated_estimator", nsys=5, temperature=300.0)
    many_batches = _run_case("stretch_equilibrated_estimator", nsys=1, temperature=300.0)
    dE_single = [row[1] for row in single_batch["rows"]]
    dE_batched = [row[1] for row in many_batches["rows"]]

    assert math.isclose(single_batch["free_energy"], many_batches["free_energy"], rel_tol=0.03, abs_tol=0.002)
    assert max(abs(a - b) for a, b in zip(dE_single, dE_batched)) < 0.1


def test_computeFreeEnergy_TI_static_gpu_force_matches_hardcoded_bond_model_and_does_not_accumulate():
    result = _run_case("bond_force_static", nsys=1, temperature=0.0)
    assert result["max_force_norm"] > 1.0
    assert result["max_abs_ref_gpu1"] < 2.0e-5
    assert result["max_abs_gpu1_gpu2"] < 1.0e-12
    assert result["max_pi_force"] < 1.0e-12
    assert max(abs(x) for x in result["sum_force"]) < 1.0e-6
    assert math.isclose(result["first_atom_ref"][0], 14.670737266540527, rel_tol=0.0, abs_tol=1.0e-6)
    assert math.isclose(result["first_atom_gpu"][0], 14.670741081237793, rel_tol=0.0, abs_tol=2.0e-5)


def test_computeFreeEnergy_TI_dynamic_gpu_step_applies_hardcoded_bond_forces_to_atoms():
    result = _run_case("bond_force_dynamic", nsys=1, temperature=0.0)
    assert math.isclose(result["max_force_norm"], 10.0, rel_tol=0.0, abs_tol=1.0e-7)
    assert result["max_force_err_from_velocity"] < 2.0e-5
    assert result["max_vel_err"] < 1.0e-6
    assert result["max_pos_err"] < 1.0e-9
    assert math.isclose(result["first_atom_ref_force"][0], 9.999822616577148, rel_tol=0.0, abs_tol=1.0e-6)
    assert math.isclose(result["first_atom_gpu_velocity"][0], 0.49999114871025085, rel_tol=0.0, abs_tol=2.0e-6)
    assert math.isclose(result["first_atom_gpu_position"][0], -3.9193038940429688, rel_tol=0.0, abs_tol=2.0e-6)


def test_TI_single_lambda_normalization_matches_force_sum_over_nProdSteps():
    result = _run_case("ti_single_lambda_normalization", nsys=1, temperature=300.0)
    avg = result["average_forces"]
    force_sum_f32 = float(np.float32(avg[2]) + np.float32(avg[3]))
    assert result["nProdSteps"] == 4000
    assert math.isclose(result["force_sum"], force_sum_f32, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(result["scale"], 1.0 / result["nProdSteps"], rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(result["dE"], -force_sum_f32 / result["nProdSteps"], rel_tol=0.0, abs_tol=1e-12)


def test_TI_single_lambda_block_split_gives_same_dE_as_single_run():
    single = _run_case("ti_single_lambda_block_single", nsys=1, temperature=300.0)
    split = _run_case("ti_single_lambda_block_split", nsys=1, temperature=300.0)
    assert single["nProdSteps"] == split["nProdSteps"] == 4000
    assert math.isclose(single["force_sum"], split["force_sum"], rel_tol=1e-9, abs_tol=1e-9)
    assert math.isclose(single["dE"], split["dE"], rel_tol=1e-9, abs_tol=1e-12)
    assert math.isclose(single["mean_sq"], split["mean_sq"], rel_tol=1e-9, abs_tol=1e-9)


def test_TI_batch_boundary_state_matches_fresh_initialization_for_same_lambda():
    result = _run_case("ti_batch_boundary_state", nsys=1, temperature=300.0)
    assert result["max_pos_diff"] < 1e-9
    assert result["max_vel_diff"] < 1e-9


def test_TI_rejects_or_handles_nMDsteps_not_divisible_by_nLambda():
    result = _run_case("ti_nondivisible_nmdsteps_guard", nsys=1, temperature=300.0)
    assert result["is_nan"] is True
    assert result["file_exists"] is False


def test_computeFreeEnergy_TI_T0_matches_in_the_first_batch_across_nSys():
    result4 = _run_case("stretch", nsys=4, temperature=0.0)
    result2 = _run_case("stretch", nsys=2, temperature=0.0)
    dE4 = [row[1] for row in result4["rows"]]
    dE2 = [row[1] for row in result2["rows"]]
    assert dE4[:2] == dE2[:2]
    assert [row[3] for row in result4["rows"][:2]] == [row[3] for row in result2["rows"][:2]]


def test_computeFreeEnergy_TI_is_batching_invariant_across_nSys():
    result4 = _run_case("stretch", nsys=4)
    result2 = _run_case("stretch", nsys=2)
    assert math.isclose(result4["free_energy"], result2["free_energy"], rel_tol=0.05, abs_tol=0.5)
    dE4 = [row[1] for row in result4["rows"]]
    dE2 = [row[1] for row in result2["rows"]]
    assert max(abs(a - b) for a, b in zip(dE4, dE2)) < 3.0


def test_computeFreeEnergy_TI_T0_second_batch_should_match_across_nSys():
    result4 = _run_case("stretch", nsys=4, temperature=0.0)
    result2 = _run_case("stretch", nsys=2, temperature=0.0)
    dE4 = [row[1] for row in result4["rows"]]
    dE2 = [row[1] for row in result2["rows"]]
    assert dE4[2:] == dE2[2:]
    assert result4["free_energy"] == result2["free_energy"]


def test_computeFreeEnergy_TI_second_batch_should_match_fresh_second_half_run():
    full = _run_case("stretch", nsys=2, temperature=0.0)
    second_half = _run_case("stretch_second_half", nsys=2, temperature=0.0)
    full_second_batch_contribution = full["rows"][3][3] - full["rows"][2][3]
    assert math.isclose(full_second_batch_contribution, second_half["free_energy"], rel_tol=1e-5, abs_tol=1e-4)


# @pytest.mark.xfail(reason="Equilibrated TI dE/dlambda does not match the 3D Gaussian-chain reference of the current bond-only model", strict=True)
def test_computeFreeEnergy_TI_equilibrated_multi_lambda_profile_matches_bond_only_gaussian_chain_reference():
    result = _run_case("stretch_equilibrated_analytic", nsys=5, temperature=300.0)
    rows = result["rows"]
    assert len(rows) == 5

    expected_lambdas = [0.0, 0.25, 0.5, 0.75, 1.0]
    assert [row[0] for row in rows] == expected_lambdas

    dr_dlambda = rows[-1][10] - rows[0][10]
    assert math.isclose(dr_dlambda, 2.0, rel_tol=0.0, abs_tol=1e-12)

    ref_total = _bond_only_gaussian_chain_reference_free_energy(rows[-1][10], temperature=300.0) - _bond_only_gaussian_chain_reference_free_energy(rows[0][10], temperature=300.0)
    for row in rows[1:4]:
        lam = row[0]
        dE = row[1]
        sigma = row[2]
        distance = row[10]
        ref_from_formula = _bond_only_gaussian_chain_reference_dedlambda(distance, dr_dlambda, temperature=300.0)
        print(f"\n[lambda={lam}] dE/dl: {dE:.6f}, Formula: {ref_from_formula:.6f}, Diff: {abs(dE - ref_from_formula):.6f}, Max allowed: {max(6.0 * sigma, 0.01):.6f}")
        assert sigma > 0.0
        assert math.isclose(
            ref_from_formula,
            {0.25: 0.00852277075773617, 0.5: 0.011363694343648227, 0.75: 0.014204617929560284}[lam],
            rel_tol=0.0,
            abs_tol=5.0e-7,
        )
        assert abs(dE - ref_from_formula) <= max(6.0 * sigma, 0.01)

    print(f"\n[Total Free Energy] Actual: {result['free_energy']:.6f}, Formula: {ref_total:.6f}, Diff: {abs(result['free_energy'] - ref_total):.6f}")
    assert math.isclose(ref_total, 0.011363694343648227, rel_tol=0.0, abs_tol=5.0e-7)
    assert abs(result["free_energy"] - ref_total) <= max(6.0 * rows[-1][4], 0.01)
