import ctypes
import ctypes.util
import json
import subprocess
import sys
from pathlib import Path

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
LOCAL_DATA_DIR = Path(__file__).resolve().parent / "data"


CASE_SCRIPT = r'''
import json, os, sys, numpy as np
repo_root = sys.argv[1]
mode = sys.argv[2]
sys.path.insert(0, repo_root)
os.chdir(os.path.join(repo_root, "cpp"))
from pyBall import MMFF_multi as mm

DT = np.float32(0.05)
L0 = np.float32(1.198)
K = np.float32(80.0)
FLIMIT = np.float32(10.0)
KB = np.float32(8.617333262145e-5)

def init_mm(T=-1.0, gamma=-1.0, nExplore=0, nRelax=0, nSys=1):
    return init_mm_xyz(
        xyz_name=os.path.join(repo_root, "tests", "test_relax_multi", "data", "N2"),
        T=T,
        gamma=gamma,
        nExplore=nExplore,
        nRelax=nRelax,
        nSys=nSys,
    )

def init_mm_xyz(xyz_name, T=-1.0, gamma=-1.0, nExplore=0, nRelax=0, nSys=1):
    base = os.path.join(repo_root, "cpp", "common_resources")
    mm.setVerbosity(-1, 0)
    mm.init(
        nSys_=nSys,
        xyz_name=xyz_name,
        sElementTypes=os.path.join(base, "ElementTypes.dat"),
        sAtomTypes=os.path.join(base, "AtomTypes.dat"),
        sBondTypes=os.path.join(base, "BondTypes.dat"),
        sAngleTypes=os.path.join(base, "AngleTypes.dat"),
        bMMFF=True,
        bEpairs=False,
        nPBC=(0, 0, 0),
        gridnPBC=(0, 0, 0),
        T=T,
        gamma=gamma,
        nExplore=nExplore,
        nRelax=nRelax,
    )
    mm.getBuffs()

def pair_force(apos):
    d = apos[1] - apos[0]
    l = np.float32(np.linalg.norm(d))
    f = d / l * ((l - L0) * K)
    fn = np.float32(np.linalg.norm(f))
    if fn > FLIMIT:
        f *= FLIMIT / fn
    return np.array([f, -f], dtype=np.float32)

def step_relax(apos, vel):
    force = pair_force(apos)
    vel = vel + force * DT
    apos = apos + vel * DT
    return apos, vel

def step_langevin_T0(apos, vel, gamma):
    force = pair_force(apos) - vel * np.float32(gamma)
    vel = vel + force * DT
    apos = apos + vel * DT
    return apos, vel

def mdstep():
    mm.MDloop(perframe=1, Ftol=0.0, iParalel=3, perVF=1)
    mm.download(bForces=True, bVel=True)

def atom_state():
    return (
        np.array(mm.gpu_atoms[0, :2, :3], copy=True),
        np.array(mm.gpu_avel[0, :2, :3], copy=True),
        np.array(mm.gpu_aforces[0, :2, :3], copy=True),
    )

def instantaneous_temperature_from_velocities(velocities):
    # In the current update kernel the velocity is integrated directly with unit mass,
    # therefore <v^2> per Cartesian DOF should approach kB*T.
    return float(np.mean(velocities * velocities) / KB)

def measure_langevin_temperature(xyz_name, target_T, nsys, gamma, ntherm, nsample):
    init_mm_xyz(xyz_name=xyz_name, T=target_T, gamma=gamma, nExplore=ntherm + nsample + 50, nRelax=0, nSys=nsys)
    natoms = int(mm.natoms)
    mdstep()  # warm-up step that enters exploring mode and sets TDrive/MDpars

    thermal_traj = []
    for _ in range(ntherm):
        mdstep()
        vel = np.array(mm.gpu_avel[:, :natoms, :3], copy=False)
        thermal_traj.append(instantaneous_temperature_from_velocities(vel))

    sample_traj = []
    for _ in range(nsample):
        mdstep()
        vel = np.array(mm.gpu_avel[:, :natoms, :3], copy=False)
        sample_traj.append(instantaneous_temperature_from_velocities(vel))

    return {
        "target_temperature": float(target_T),
        "natoms": natoms,
        "nsys": nsys,
        "ntherm": ntherm,
        "nsample": nsample,
        "thermal_start_mean": float(np.mean(thermal_traj[:10])),
        "thermal_end_mean": float(np.mean(thermal_traj[-10:])),
        "sample_mean_temperature": float(np.mean(sample_traj)),
        "sample_std_temperature": float(np.std(sample_traj)),
    }

def run_integrate():
    init_mm()
    mdstep()  # warm-up MDpars
    apos0 = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32)
    vel0 = np.zeros((2, 3), dtype=np.float32)
    mm.gpu_atoms[0, :2, :3] = apos0
    mm.gpu_avel[0, :, :3] = 0.0
    mm.upload(bParams=False, bForces=False, bVel=True)
    apos_ref, vel_ref = step_relax(apos0.copy(), vel0.copy())
    mdstep()
    apos, vel, frc = atom_state()
    return {
        "max_pos_err": float(np.max(np.abs(apos - apos_ref))),
        "max_vel_err": float(np.max(np.abs(vel - vel_ref))),
        "max_force_after_download": float(np.max(np.abs(frc))),
        "moved": bool(np.max(np.abs(apos - apos0)) > 1e-6),
    }

def run_langevin():
    init_mm(T=0.0, gamma=2.0, nExplore=10, nRelax=0)
    mdstep()  # enter host-side exploring state, populate TDrive/MDpars
    apos0 = np.array([[-0.599, 0.0, 0.0], [0.599, 0.0, 0.0]], dtype=np.float32)
    vel0 = np.array([[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]], dtype=np.float32)
    mm.gpu_atoms[0, :2, :3] = apos0
    mm.gpu_avel[0, :, :3] = 0.0
    mm.gpu_avel[0, :2, :3] = vel0
    mm.upload(bParams=False, bForces=False, bVel=True)
    apos_ref, vel_ref = step_langevin_T0(apos0.copy(), vel0.copy(), gamma=2.0)
    mdstep()
    apos, vel, frc = atom_state()
    return {
        "max_pos_err": float(np.max(np.abs(apos - apos_ref))),
        "max_vel_err": float(np.max(np.abs(vel - vel_ref))),
        "max_force_after_download": float(np.max(np.abs(frc))),
    }

def run_relax_path():
    init_mm()
    mdstep()  # warm-up MDpars
    apos_ref = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32)
    vel_ref = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32)
    mm.gpu_atoms[0, :2, :3] = apos_ref
    mm.gpu_avel[0, :, :3] = 0.0
    mm.gpu_avel[0, :2, :3] = vel_ref
    mm.upload(bParams=False, bForces=False, bVel=True)
    max_pos_err = max_vel_err = 0.0
    for _ in range(4):
        apos_ref, vel_ref = step_relax(apos_ref, vel_ref)
        mdstep()
        apos, vel, _ = atom_state()
        max_pos_err = max(max_pos_err, float(np.max(np.abs(apos - apos_ref))))
        max_vel_err = max(max_vel_err, float(np.max(np.abs(vel - vel_ref))))
    return {"max_pos_err": max_pos_err, "max_vel_err": max_vel_err}

def run_langevin_temperature():
    return measure_langevin_temperature(
        xyz_name=os.path.join(repo_root, "tests", "test_relax_multi", "data", "N2"),
        target_T=300.0,
        nsys=16,
        gamma=0.5,
        ntherm=120,
        nsample=200,
    )

def run_langevin_temperature_propandiol(target_T):
    return measure_langevin_temperature(
        xyz_name=os.path.join(repo_root, "cpp", "common_resources", "xyz", "propandiol"),
        target_T=target_T,
        nsys=8,
        gamma=0.5,
        ntherm=200,
        nsample=250,
    )

try:
    if mode == "langevin_temperature_propandiol":
        out = run_langevin_temperature_propandiol(float(sys.argv[3]))
    else:
        out = {
            "integrate": run_integrate,
            "langevin": run_langevin,
            "relax_path": run_relax_path,
            "langevin_temperature": run_langevin_temperature,
        }[mode]()
    print("JSON_RESULT=" + json.dumps(out, sort_keys=True))
finally:
    try: mm.clear()
    except Exception: pass
'''


def _run_case(mode: str, *args: object) -> dict:
    proc = subprocess.run(
        [sys.executable, "-c", CASE_SCRIPT, str(REPO_ROOT), mode, *[str(a) for a in args]],
        cwd=CPP_DIR,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, (
        f"updateAtomsMMFFf4 subprocess failed for {mode}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )
    for line in reversed(proc.stdout.splitlines()):
        if line.startswith("JSON_RESULT="):
            return json.loads(line.split("=", 1)[1])
    raise AssertionError(f"JSON result missing for {mode}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")


def test_updateAtomsMMFFf4_integrates_motion_and_downloads_host_mirrors():
    result = _run_case("integrate")
    assert result["moved"]
    assert result["max_pos_err"] < 1e-6
    assert result["max_vel_err"] < 1e-6
    assert result["max_force_after_download"] < 1e-7


def test_updateAtomsMMFFf4_applies_deterministic_langevin_damping_at_T0():
    result = _run_case("langevin")
    assert result["max_pos_err"] < 1e-6
    assert result["max_vel_err"] < 1e-6
    assert result["max_force_after_download"] < 1e-7


def test_updateAtomsMMFFf4_current_relax_path_matches_neutral_fire_host_settings():
    result = _run_case("relax_path")
    assert result["max_pos_err"] < 1e-6
    assert result["max_vel_err"] < 1e-6


def test_updateAtomsMMFFf4_langevin_thermalizes_to_the_target_temperature():
    result = _run_case("langevin_temperature")
    target = result["target_temperature"]
    assert result["thermal_end_mean"] > 0.8 * target
    assert 0.85 * target < result["sample_mean_temperature"] < 1.20 * target
    assert result["sample_std_temperature"] > 0.0


def test_updateAtomsMMFFf4_langevin_scales_temperature_on_propandiol():
    targets = [100.0, 300.0, 600.0]
    means = []

    for target in targets:
        result = _run_case("langevin_temperature_propandiol", target)
        means.append(result["sample_mean_temperature"])
        assert result["natoms"] > 2
        assert result["thermal_end_mean"] > 0.8 * target
        assert 0.85 * target < result["sample_mean_temperature"] < 1.25 * target
        assert result["sample_std_temperature"] > 0.0

    assert means[0] < means[1] < means[2]
