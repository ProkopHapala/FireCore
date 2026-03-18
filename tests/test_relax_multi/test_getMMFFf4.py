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
DATA_DIR = CPP_DIR / "common_resources"
LOCAL_DATA_DIR = Path(__file__).resolve().parent / "data"

CASE_SCRIPT = r'''
import json
import os
import sys
import numpy as np

repo_root = sys.argv[1]
xyz_base  = sys.argv[2]
sys.path.insert(0, repo_root)
os.chdir(os.path.join(repo_root, "cpp"))

from pyBall import MMFF_multi as mm


KERNEL_BOND_L0 = np.float32(1.198)
KERNEL_BOND_K  = np.float32(80.0)


def eval_bond_python(h):
    l2 = np.float32(np.dot(h, h))
    l  = np.float32(np.sqrt(l2))
    h  = h / l
    fr = (l - KERNEL_BOND_L0) * KERNEL_BOND_K
    return h * fr


def eval_getMMFFf4_python_reference(mm):
    natoms = int(mm.natoms)
    nnode  = int(mm.nnode)
    nvecs  = int(mm.nvecs)

    apos      = np.array(mm.gpu_atoms[0, :, :3], dtype=np.float32, copy=True)
    neighs    = np.array(mm.gpu_neighs[0, :natoms, :], dtype=np.int32, copy=True)
    bkneighs  = np.array(mm.gpu_bkNeighs[0, :nvecs, :], dtype=np.int32, copy=True)
    aforces   = np.zeros((nvecs, 3), dtype=np.float32)

    valid_bk = bkneighs[bkneighs >= 0]
    nneigh_force = int(valid_bk.max()) + 1 if valid_bk.size else 0
    neigh_force  = np.zeros((nneigh_force, 3), dtype=np.float32)

    for iG in range(nnode):
        pa = apos[iG]
        fa = np.zeros(3, dtype=np.float32)
        for i, ing in enumerate(neighs[iG]):
            ing = int(ing)
            if ing < 0:
                break
            if iG < ing:
                f1 = eval_bond_python(apos[ing] - pa)
                fa += f1
                neigh_force[iG * 4 + i] -= f1
        aforces[iG] = fa

    for iv in range(nvecs):
        fe = aforces[iv]
        for ib in bkneighs[iv]:
            ib = int(ib)
            if ib >= 0:
                fe = fe + neigh_force[ib]
        aforces[iv] = fe

    return aforces

base = os.path.join(repo_root, "cpp", "common_resources")
mm.setVerbosity(-1, 0)
mm.init(
    nSys_=1,
    xyz_name=xyz_base,
    sElementTypes=os.path.join(base, "ElementTypes.dat"),
    sAtomTypes=os.path.join(base, "AtomTypes.dat"),
    sBondTypes=os.path.join(base, "BondTypes.dat"),
    sAngleTypes=os.path.join(base, "AngleTypes.dat"),
    bMMFF=True,
    bEpairs=False,
    nPBC=(0, 0, 0),
    gridnPBC=(0, 0, 0),
    T=-1,
    gamma=-1,
)
mm.getBuffs()
pos_before = np.array(mm.gpu_atoms[0, :, :3], copy=True)

py_forces = eval_getMMFFf4_python_reference(mm)

mm.eval_getMMFFf4_ocl()
mm.download(bForces=False, bVel=False)

pos_after = np.array(mm.gpu_atoms[0, :, :3], copy=True)
gpu_forces = np.array(mm.gpu_aforces[0, :, :3], copy=True)

natoms = int(mm.natoms)
npi    = int(mm.npi)
py_atoms = py_forces[:natoms]
py_pi    = py_forces[natoms:natoms+npi]
gpu_atoms = gpu_forces[:natoms]
gpu_pi    = gpu_forces[natoms:natoms+npi]

result = {
    "natoms": natoms,
    "npi": npi,
    "max_abs_atom_error": float(np.max(np.abs(py_atoms - gpu_atoms))),
    "max_abs_pi_error": float(np.max(np.abs(py_pi - gpu_pi))) if npi else 0.0,
    "max_abs_position_shift": float(np.max(np.abs(pos_after - pos_before))),
    "max_abs_atom_force": float(np.max(np.abs(py_atoms))),
    "max_abs_pi_force": float(np.max(np.abs(py_pi))) if npi else 0.0,
}
print("JSON_RESULT=" + json.dumps(result, sort_keys=True))
'''


def _xyz_base(path: Path) -> str:
    return str(path.with_suffix("")) if path.suffix == ".xyz" else str(path)


def _run_case(xyz_path: Path) -> dict:
    proc = subprocess.run(
        [sys.executable, "-c", CASE_SCRIPT, str(REPO_ROOT), _xyz_base(xyz_path)],
        cwd=CPP_DIR,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, (
        f"MMFF_multi subprocess failed for {xyz_path}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )
    for line in reversed(proc.stdout.splitlines()):
        if line.startswith("JSON_RESULT="):
            return json.loads(line.split("=", 1)[1])
    raise AssertionError(f"JSON result missing for {xyz_path}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")


def test_getMMFFf4_matches_python_reference_for_bond_angle_and_torsion_cases():
    cases = {
        "bond_only_N2": LOCAL_DATA_DIR / "N2.xyz",
        "bond_angle_H2O": DATA_DIR / "xyz" / "H2O.xyz",
        "torsion_propandiol": DATA_DIR / "xyz" / "propandiol.xyz",
    }

    failures = []
    for name, xyz_path in cases.items():
        result = _run_case(xyz_path)
        if result["max_abs_atom_error"] >= 1e-4:
            failures.append(f"{name}: atom force error {result['max_abs_atom_error']}")
        if result["max_abs_pi_error"] >= 1e-4:
            failures.append(f"{name}: pi force error {result['max_abs_pi_error']}")
        if result["max_abs_position_shift"] >= 1e-12:
            failures.append(f"{name}: unexpected position shift {result['max_abs_position_shift']}")

    assert not failures, "getMMFFf4 mismatches Python reference:\n- " + "\n- ".join(failures)
