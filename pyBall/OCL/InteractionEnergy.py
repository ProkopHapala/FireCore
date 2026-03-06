"""
InteractionEnergy.py - Headless PyOpenCL module for molecule-substrate interaction energy scanning.
No GUI dependency. Uses OpenCLBase for GPU setup and kernel management.
Uses MMparams (AtomTypes.dat) for proper MMFF REQ parameter initialization.

Usage:
    from pyBall.OCL.InteractionEnergy import InteractionScanner
    scanner = InteractionScanner()
    scanner.load_molecule_xyz('PTCDA.xyz')
    scanner.load_substrate_xyz('NaCl_8x8_L3.xyz')
    results = scanner.evaluate(transforms)
"""

import os
import math
import numpy as np
import pyopencl as cl
from .OpenCLBase import OpenCLBase
from . import ScanUtils
from . import MMparams

# Coulomb constant: e^2/(4*pi*eps0) in eV*Angstrom
COULOMB_CONST = 14.3996

# Path to common_resources for AtomTypes.dat, ElementTypes.dat
_BASE_PATH = os.path.dirname(os.path.abspath(__file__))
_DATA_PATH = os.path.join(_BASE_PATH, "../../cpp/common_resources/")


def _half_step_from_coords(cs):
    cu = np.unique(np.round(np.asarray(cs, dtype=np.float64), 8))
    if len(cu) < 2:
        return 0.0
    ds = np.diff(np.sort(cu))
    ds = ds[ds > 1e-8]
    if len(ds) == 0:
        return 0.0
    return 0.5 * float(ds.min())


def _rect_sheet_F(X, Y, Z):
    R = np.sqrt(X*X + Y*Y + Z*Z)
    t1 = X * np.log(np.maximum(Y + R, 1e-12))
    t2 = Y * np.log(np.maximum(X + R, 1e-12))
    t3 = Z * np.arctan2(X * Y, Z * R + 1e-12)
    return t1 + t2 - t3


def rect_sheet_potential(xs, ys, zs, sigma, xmin, xmax, ymin, ymax):
    x0 = xs - xmin
    x1 = xs - xmax
    y0 = ys - ymin
    y1 = ys - ymax
    return sigma * (_rect_sheet_F(x0, y0, zs) - _rect_sheet_F(x1, y0, zs) - _rect_sheet_F(x0, y1, zs) + _rect_sheet_F(x1, y1, zs))


def rect_dipole_potential(xs, ys, zs, Px, Py, Pz, xmin, xmax, ymin, ymax):
    sum_omega = 0.0
    sum_logy  = 0.0
    sum_logx  = 0.0
    for xc, sx in ((xmin, 1.0), (xmax, -1.0)):
        X = xs - xc
        for yc, sy in ((ymin, 1.0), (ymax, -1.0)):
            Y = ys - yc
            s = sx * sy
            R = np.sqrt(X*X + Y*Y + zs*zs)
            sum_omega += s * np.arctan2(X * Y, zs * R + 1e-12)
            sum_logy  += s * np.log(np.maximum(Y + R, 1e-12))
            sum_logx  += s * np.log(np.maximum(X + R, 1e-12))
    return Pz * sum_omega - Px * sum_logy - Py * sum_logx


def rect_quadrupole_potential(xs, ys, zs, Qxx, Qxy, Qyy, Qxz, Qyz, Qzz, xmin, xmax, ymin, ymax, eps=1e-3):
    I0   = rect_sheet_potential(xs, ys, zs, 1.0, xmin, xmax, ymin, ymax)
    Ixp  = rect_sheet_potential(xs + eps, ys, zs, 1.0, xmin, xmax, ymin, ymax)
    Ixm  = rect_sheet_potential(xs - eps, ys, zs, 1.0, xmin, xmax, ymin, ymax)
    Iyp  = rect_sheet_potential(xs, ys + eps, zs, 1.0, xmin, xmax, ymin, ymax)
    Iym  = rect_sheet_potential(xs, ys - eps, zs, 1.0, xmin, xmax, ymin, ymax)
    Izp  = rect_sheet_potential(xs, ys, zs + eps, 1.0, xmin, xmax, ymin, ymax)
    Izm  = rect_sheet_potential(xs, ys, zs - eps, 1.0, xmin, xmax, ymin, ymax)
    Ixyp = rect_sheet_potential(xs + eps, ys + eps, zs, 1.0, xmin, xmax, ymin, ymax)
    Ixym = rect_sheet_potential(xs + eps, ys - eps, zs, 1.0, xmin, xmax, ymin, ymax)
    Imyp = rect_sheet_potential(xs - eps, ys + eps, zs, 1.0, xmin, xmax, ymin, ymax)
    Imym = rect_sheet_potential(xs - eps, ys - eps, zs, 1.0, xmin, xmax, ymin, ymax)
    Ixzp = rect_sheet_potential(xs + eps, ys, zs + eps, 1.0, xmin, xmax, ymin, ymax)
    Ixzm = rect_sheet_potential(xs + eps, ys, zs - eps, 1.0, xmin, xmax, ymin, ymax)
    Imzp = rect_sheet_potential(xs - eps, ys, zs + eps, 1.0, xmin, xmax, ymin, ymax)
    Imzm = rect_sheet_potential(xs - eps, ys, zs - eps, 1.0, xmin, xmax, ymin, ymax)
    Iyzp = rect_sheet_potential(xs, ys + eps, zs + eps, 1.0, xmin, xmax, ymin, ymax)
    Iyzm = rect_sheet_potential(xs, ys + eps, zs - eps, 1.0, xmin, xmax, ymin, ymax)
    Imzzp = rect_sheet_potential(xs, ys - eps, zs + eps, 1.0, xmin, xmax, ymin, ymax)
    Imzzm = rect_sheet_potential(xs, ys - eps, zs - eps, 1.0, xmin, xmax, ymin, ymax)
    dxx = (Ixp - 2.0 * I0 + Ixm) / (eps * eps)
    dyy = (Iyp - 2.0 * I0 + Iym) / (eps * eps)
    dzz = (Izp - 2.0 * I0 + Izm) / (eps * eps)
    dxy = (Ixyp - Ixym - Imyp + Imym) / (4.0 * eps * eps)
    dxz = (Ixzp - Ixzm - Imzp + Imzm) / (4.0 * eps * eps)
    dyz = (Iyzp - Iyzm - Imzzp + Imzzm) / (4.0 * eps * eps)
    return (Qxx * dxx + Qyy * dyy + Qzz * dzz + 2.0 * Qxy * dxy + 2.0 * Qxz * dxz + 2.0 * Qyz * dyz) / 6.0


def make_REQs_from_enames(enames, qs, atom_types, type_map=None):
    """Build MMFF REQ array (R, sqrt(E), Q, H) from element/type names and charges.
    enames:     list of element names or atom type names
    qs:         array of charges (from xyz 4th column)
    atom_types: dict from MMparams.read_atom_types
    type_map:   optional dict mapping element name -> atom type name (e.g. {'C': 'C_R'})
    Returns: (N,4) float32 array in MMFF REQ layout"""
    N = len(enames)
    REQs = np.zeros((N, 4), dtype=np.float32)
    tm = type_map or {}
    alias = {
        'C.ar': 'C_R', 'C_ar': 'C_R', 'N.ar': 'N_R', 'N_ar': 'N_R', 'O.ar': 'O_R', 'O_ar': 'O_R',
        'C.2': 'C_2', 'C_2': 'C_2', 'N.2': 'N_2', 'N_2': 'N_2', 'O.2': 'O_2', 'O_2': 'O_2',
        'C.1': 'C_1', 'C_1': 'C_1', 'N.1': 'N_1', 'N_1': 'N_1',
        'C.3': 'C_3', 'C_3': 'C_3', 'N.3': 'N_3', 'N_3': 'N_3', 'O.3': 'O_3', 'O_3': 'O_3',
    }
    for i in range(N):
        raw = enames[i]
        aname = tm.get(raw, alias.get(raw, raw))
        if aname not in atom_types:
            raise KeyError(f"Atom type '{aname}' (from '{raw}') not found in AtomTypes.dat for atom {i}")
        at = atom_types[aname]
        REQs[i, 0] = float(at.RvdW)
        REQs[i, 1] = math.sqrt(float(at.EvdW))
        REQs[i, 2] = float(qs[i])
        REQs[i, 3] = float(at.Hb)
    return REQs


def load_xyz_with_REQs(fname, atom_types=None, type_map=None):
    """Load XYZ file and build MMFF REQ parameters from AtomTypes.dat.
    Returns: positions (N,3) float64, REQs (N,4) float32, enames list, Zs array"""
    from pyBall import atomicUtils as au
    apos, Zs, enames, qs, comment = au.load_xyz(fname=fname, bReadN=True)
    lvec = None
    if isinstance(comment, str):
        s = comment.strip()
        if s.startswith('#'): s = s[1:].lstrip()
        if 'lvs' in s:
            i = s.find('lvs')
            parts = s[i+3:].split()
            if len(parts) >= 9:
                nums = [float(x) for x in parts[:9]]
                lvec = np.array(nums, dtype=np.float32).reshape(3, 3)
    if atom_types is None:
        etypes = MMparams.read_element_types(os.path.join(_DATA_PATH, 'ElementTypes.dat'))
        atom_types = MMparams.read_atom_types(os.path.join(_DATA_PATH, 'AtomTypes.dat'), etypes)
    REQs = make_REQs_from_enames(enames, qs, atom_types, type_map=type_map)
    return apos, REQs, enames, Zs, lvec


class InteractionScanner(OpenCLBase):
    """GPU-accelerated molecule-substrate interaction energy scanner.
    Uses MMFF REQ parameters (R, sqrt(E), Q, H) from AtomTypes.dat.
    
    Supports:
    - LJ or Morse van der Waals interactions (switchable)
    - Coulomb electrostatics (switchable)
    - H-bond corrections (switchable)
    - Rigid scan and constrained relaxation modes
    - All scan types from ScanUtils
    """

    def __init__(self, nloc=32, device_index=0):
        super().__init__(nloc=nloc, device_index=device_index)
        cl_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cl')
        kernel_path = os.path.join(cl_dir, 'InteractionEnergy.cl')
        self.load_program(kernel_path=kernel_path, bPrint=False)
        self.krn_evaluate   = self.prg.evaluate_interaction
        self.krn_relax      = self.prg.relax_constrained
        # Load atom types from AtomTypes.dat for parameter assignment
        etypes = MMparams.read_element_types(os.path.join(_DATA_PATH, 'ElementTypes.dat'))
        self.atom_types = MMparams.read_atom_types(os.path.join(_DATA_PATH, 'AtomTypes.dat'), etypes)
        # Default physics settings
        self.enable_LJ      = True
        self.enable_Coulomb  = True
        self.enable_HBond    = False
        self.enable_Morse    = False
        self.Coulomb_const   = np.float32(COULOMB_CONST)
        self.Morse_alpha     = np.float32(1.8)
        # Macro-potential subtraction (continuum embedding)
        self.enable_macro    = False
        self.macro_P         = np.zeros(4, dtype=np.float32)   # (Px,Py,Pz,0) polarization [e/Ang]
        self.macro_AB        = np.zeros(4, dtype=np.float32)   # (Ax,By,0,0) half-sizes [Ang]
        self._macro_layers   = None
        self._macro_bounds   = None
        # Relaxation defaults
        self.spring_k        = np.float32(5.0)
        self.relax_dt        = np.float32(0.005)
        self.relax_nsteps    = 100
        # State
        self._mol_buf     = None
        self._mol_req_buf = None
        self._sub_buf     = None
        self._sub_req_buf = None
        self.nPBC = np.array([1,1,0], dtype=np.int32)
        self.lvec = None
        self.wrap_PBC = False
        self._nmol = 0
        self._nsub = 0
        # Keep host copies for GUI / inspection
        self.mol_apos   = None
        self.mol_REQs   = None
        self.mol_enames = None
        self.sub_apos   = None
        self.sub_REQs   = None
        self.sub_enames = None

    # ======== Data loading ========

    def set_molecule(self, positions, REQs):
        """Upload molecule atoms and MMFF REQ params (R, sqrt(E), Q, H) to GPU."""
        N = len(positions)
        self._nmol = N
        self.mol_apos = np.array(positions, dtype=np.float64)
        self.mol_REQs = np.array(REQs, dtype=np.float32)
        pos4 = np.zeros((N, 4), dtype=np.float32)
        pos4[:, :3] = positions
        self._mol_buf     = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=pos4)
        self._mol_req_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.ascontiguousarray(REQs, dtype=np.float32))

    def set_substrate(self, positions, REQs):
        """Upload substrate atoms and MMFF REQ params (R, sqrt(E), Q, H) to GPU."""
        N = len(positions)
        self._nsub = N
        self.sub_apos = np.array(positions, dtype=np.float64)
        self.sub_REQs = np.array(REQs, dtype=np.float32)
        pos4 = np.zeros((N, 4), dtype=np.float32)
        pos4[:, :3] = positions
        self._sub_buf     = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=pos4)
        self._sub_req_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.ascontiguousarray(REQs, dtype=np.float32))

    def load_molecule_xyz(self, fname, type_map=None):
        """Load molecule from XYZ file with auto REQ assignment from AtomTypes.dat."""
        apos, REQs, enames, Zs, lvec = load_xyz_with_REQs(fname, self.atom_types, type_map=type_map)
        self.mol_enames = enames
        self.set_molecule(apos, REQs)
        return apos, REQs, enames

    def load_substrate_xyz(self, fname, type_map=None):
        """Load substrate from XYZ file with auto REQ assignment from AtomTypes.dat."""
        apos, REQs, enames, Zs, lvec = load_xyz_with_REQs(fname, self.atom_types, type_map=type_map)
        self.sub_enames = enames
        self.lvec = lvec
        self.set_substrate(apos, REQs)
        self._update_macro_from_substrate()
        return apos, REQs, enames

    def _update_macro_from_substrate(self):
        if self.lvec is None: return
        if self.sub_apos is None or self.sub_REQs is None: return
        a = np.array(self.lvec[0], dtype=np.float64)
        b = np.array(self.lvec[1], dtype=np.float64)
        A = float(np.linalg.norm(np.cross(a, b)))
        if A < 1e-12: raise ValueError(f"Invalid surface cell area A={A} from a={a} b={b}")
        qs = self.sub_REQs[:, 2].astype(np.float64)
        ps = self.sub_apos.astype(np.float64)
        mu = (qs[:, None] * ps).sum(axis=0)  # [e*Ang]
        P = mu / A                           # [e/Ang]
        self.macro_P[:] = (P[0], P[1], P[2], 0.0)
        Ax = (float(self.nPBC[0]) + 0.5) * float(np.linalg.norm(a))
        By = (float(self.nPBC[1]) + 0.5) * float(np.linalg.norm(b))
        self.macro_AB[:] = (Ax, By, 0.0, 0.0)
        xmin0, xmax0 = float(self.sub_apos[:, 0].min()), float(self.sub_apos[:, 0].max())
        ymin0, ymax0 = float(self.sub_apos[:, 1].min()), float(self.sub_apos[:, 1].max())
        hx = _half_step_from_coords(self.sub_apos[:, 0])
        hy = _half_step_from_coords(self.sub_apos[:, 1])
        La = float(np.linalg.norm(a))
        Lb = float(np.linalg.norm(b))
        xmin = xmin0 - hx - float(self.nPBC[0]) * La
        xmax = xmax0 + hx + float(self.nPBC[0]) * La
        ymin = ymin0 - hy - float(self.nPBC[1]) * Lb
        ymax = ymax0 + hy + float(self.nPBC[1]) * Lb
        self._macro_bounds = (xmin, xmax, ymin, ymax)
        zs = self.sub_apos[:, 2].astype(np.float64)
        qs = self.sub_REQs[:, 2].astype(np.float64)
        zuniq = []
        qsum = []
        izs = np.argsort(zs)
        tol = 1e-4
        for i in izs:
            z = zs[i]
            q = qs[i]
            if not zuniq or abs(z - zuniq[-1]) > tol:
                zuniq.append(z)
                qsum.append(q)
            else:
                qsum[-1] += q
        sigmas = np.array(qsum, dtype=np.float64) / A
        mus = []
        quads = []
        cx = 0.5 * (xmin + xmax)
        cy = 0.5 * (ymin + ymax)
        for z in zuniq:
            m = np.abs(zs - z) < tol
            dx = ps[m, 0] - cx
            dy = ps[m, 1] - cy
            dz = ps[m, 2] - z
            qq = qs[m]
            mus.append(np.array([np.sum(qq * dx), np.sum(qq * dy), np.sum(qq * dz)]))
            qzz = np.sum(qq * (2.0 * dz * dz - dx * dx - dy * dy))
            qxx = np.sum(qq * (2.0 * dx * dx - dy * dy - dz * dz))
            qyy = np.sum(qq * (2.0 * dy * dy - dx * dx - dz * dz))
            qxy = np.sum(qq * (3.0 * dx * dy))
            qxz = np.sum(qq * (3.0 * dx * dz))
            qyz = np.sum(qq * (3.0 * dy * dz))
            quads.append((qxx / A, qxy / A, qyy / A, qxz / A, qyz / A, qzz / A))
        self._macro_layers = (np.array(zuniq, dtype=np.float64), sigmas, np.array(mus, dtype=np.float64) / A, np.array(quads, dtype=np.float64))

    def _apply_macro_correction(self, transforms, out, relaxed_pos=None):
        if (not self.enable_macro) or (self._macro_layers is None) or (self._macro_bounds is None):
            return out
        zls, sigmas, Pmus, Qdens = self._macro_layers
        xmin, xmax, ymin, ymax = self._macro_bounds
        qs = self.mol_REQs[:, 2].astype(np.float64)
        if relaxed_pos is None:
            T = np.ascontiguousarray(transforms, dtype=np.float64).reshape(-1, 3, 4)
            p0 = self.mol_apos.astype(np.float64)
            xyz = np.empty((len(T), len(p0), 3), dtype=np.float64)
            xyz[:, :, 0] = p0[None, :, 0] * T[:, None, 0, 0] + p0[None, :, 1] * T[:, None, 0, 1] + p0[None, :, 2] * T[:, None, 0, 2] + T[:, None, 0, 3]
            xyz[:, :, 1] = p0[None, :, 0] * T[:, None, 1, 0] + p0[None, :, 1] * T[:, None, 1, 1] + p0[None, :, 2] * T[:, None, 1, 2] + T[:, None, 1, 3]
            xyz[:, :, 2] = p0[None, :, 0] * T[:, None, 2, 0] + p0[None, :, 1] * T[:, None, 2, 1] + p0[None, :, 2] * T[:, None, 2, 2] + T[:, None, 2, 3]
        else:
            xyz = relaxed_pos.astype(np.float64)
        dE = np.zeros(xyz.shape[0], dtype=np.float64)
        xs = xyz[:, :, 0]
        ys = xyz[:, :, 1]
        for zl, sigma, Pmu, Qd in zip(zls, sigmas, Pmus, Qdens):
            dz = xyz[:, :, 2] - zl
            phi = rect_sheet_potential(xs, ys, dz, sigma, xmin, xmax, ymin, ymax)
            phi += rect_dipole_potential(xs, ys, dz, Pmu[0], Pmu[1], Pmu[2], xmin, xmax, ymin, ymax)
            phi += rect_quadrupole_potential(xs, ys, dz, Qd[0], Qd[1], Qd[2], Qd[3], Qd[4], Qd[5], xmin, xmax, ymin, ymax)
            dE -= self.Coulomb_const * np.sum(phi * qs[None, :], axis=1)
        out['Coulomb'] = out['Coulomb'].astype(np.float64) + dE
        out['total'] = out['total'].astype(np.float64) + dE
        return out

    # ======== Evaluation ========

    def _wrap_transforms_PBC(self, transforms):
        """Wrap the xy translation of packed transforms back into the primary cell.

        With a finite image sum (limited nPBC), this makes energies approximately invariant under
        translation by lattice vectors because the evaluated image window stays centered.
        """
        if (not self.wrap_PBC) or (self.lvec is None):
            return transforms
        a = np.array(self.lvec[0], dtype=np.float64)
        b = np.array(self.lvec[1], dtype=np.float64)
        M = np.array([[a[0], b[0]], [a[1], b[1]]], dtype=np.float64)
        det = float(np.linalg.det(M))
        if abs(det) < 1e-12:
            raise ValueError(f"Degenerate lattice vectors for PBC wrap: det={det} a={a} b={b}")
        invM = np.linalg.inv(M)
        T = np.array(transforms, copy=True, dtype=np.float32).reshape(-1, 3, 4)
        txy = T[:, 0:2, 3].astype(np.float64)          # (N,2)
        frac = (invM @ txy.T).T                        # (N,2)
        frac -= np.round(frac)                         # wrap to [-0.5,0.5)
        txy2 = (M @ frac.T).T
        T[:, 0, 3] = txy2[:, 0]
        T[:, 1, 3] = txy2[:, 1]
        return T.reshape(-1, 12)

    def evaluate(self, transforms):
        """Evaluate interaction energy for a batch of rigid-body transforms.
        transforms: (N, 12) float32 array from ScanUtils.pack_transforms
        Returns: dict with 'total', 'LJ', 'Coulomb', 'HBond' arrays of shape (N,)"""
        assert self._mol_buf is not None, "Call set_molecule() or load_molecule_xyz() first"
        assert self._sub_buf is not None, "Call set_substrate() or load_substrate_xyz() first"
        nconf = len(transforms)
        transforms = self._wrap_transforms_PBC(transforms)
        transforms = np.ascontiguousarray(transforms, dtype=np.float32)
        trans_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=transforms)
        sz = nconf * 4  # float32 bytes
        res_total = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=sz)
        res_lj    = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=sz)
        res_coul  = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=sz)
        res_hb    = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=sz)
        local_mem = cl.LocalMemory(self.nloc * 4)
        global_size = (nconf * self.nloc,)
        local_size  = (self.nloc,)
        import pyopencl.array as cl_array
        if self.lvec is None:
            la = cl_array.vec.make_float4(100.0, 0.0, 0.0, 0.0)
            lb = cl_array.vec.make_float4(0.0, 100.0, 0.0, 0.0)
            lc = cl_array.vec.make_float4(0.0, 0.0, 100.0, 0.0)
        else:
            la = cl_array.vec.make_float4(float(self.lvec[0, 0]), float(self.lvec[0, 1]), float(self.lvec[0, 2]), 0.0)
            lb = cl_array.vec.make_float4(float(self.lvec[1, 0]), float(self.lvec[1, 1]), float(self.lvec[1, 2]), 0.0)
            lc = cl_array.vec.make_float4(float(self.lvec[2, 0]), float(self.lvec[2, 1]), float(self.lvec[2, 2]), 0.0)
        npbc = cl_array.vec.make_int3(int(self.nPBC[0]), int(self.nPBC[1]), int(self.nPBC[2]))
        macro_P  = cl_array.vec.make_float4(0.0, 0.0, 0.0, 0.0)
        macro_AB = cl_array.vec.make_float4(0.0, 0.0, 0.0, 0.0)
        self.krn_evaluate(
            self.queue, global_size, local_size,
            self._mol_buf, self._mol_req_buf, np.int32(self._nmol),
            self._sub_buf, self._sub_req_buf, np.int32(self._nsub),
            npbc, la, lb, lc,
            trans_buf,
            np.int32(int(self.enable_LJ)), np.int32(int(self.enable_Coulomb)),
            np.int32(int(self.enable_HBond)), np.int32(int(self.enable_Morse)),
            np.int32(0), macro_P, macro_AB,
            np.float32(self.Coulomb_const), np.float32(self.Morse_alpha),
            local_mem,
            res_total, res_lj, res_coul, res_hb,
        ).wait()
        out = {}
        for name, buf in [('total', res_total), ('LJ', res_lj), ('Coulomb', res_coul), ('HBond', res_hb)]:
            arr = np.empty(nconf, dtype=np.float32)
            cl.enqueue_copy(self.queue, arr, buf)
            out[name] = arr
        return self._apply_macro_correction(transforms, out)

    def evaluate_relaxed(self, transforms):
        """Evaluate with constrained relaxation.
        Returns: dict with energy components + 'relaxed_pos' (nconf, nmol, 3)."""
        assert self._mol_buf is not None, "Call set_molecule() or load_molecule_xyz() first"
        assert self._sub_buf is not None, "Call set_substrate() or load_substrate_xyz() first"
        nconf = len(transforms)
        transforms = self._wrap_transforms_PBC(transforms)
        transforms = np.ascontiguousarray(transforms, dtype=np.float32)
        trans_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=transforms)
        sz = nconf * 4
        res_total = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=sz)
        res_lj    = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=sz)
        res_coul  = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=sz)
        res_hb    = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=sz)
        relax_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=nconf * self._nmol * 4 * 4)
        local_mem = cl.LocalMemory(self.nloc * 4)
        global_size = (nconf * self.nloc,)
        local_size  = (self.nloc,)
        import pyopencl.array as cl_array
        if self.lvec is None:
            la = cl_array.vec.make_float4(100.0, 0.0, 0.0, 0.0)
            lb = cl_array.vec.make_float4(0.0, 100.0, 0.0, 0.0)
            lc = cl_array.vec.make_float4(0.0, 0.0, 100.0, 0.0)
        else:
            la = cl_array.vec.make_float4(float(self.lvec[0, 0]), float(self.lvec[0, 1]), float(self.lvec[0, 2]), 0.0)
            lb = cl_array.vec.make_float4(float(self.lvec[1, 0]), float(self.lvec[1, 1]), float(self.lvec[1, 2]), 0.0)
            lc = cl_array.vec.make_float4(float(self.lvec[2, 0]), float(self.lvec[2, 1]), float(self.lvec[2, 2]), 0.0)
        npbc = cl_array.vec.make_int3(int(self.nPBC[0]), int(self.nPBC[1]), int(self.nPBC[2]))
        macro_P  = cl_array.vec.make_float4(0.0, 0.0, 0.0, 0.0)
        macro_AB = cl_array.vec.make_float4(0.0, 0.0, 0.0, 0.0)
        self.krn_relax(
            self.queue, global_size, local_size,
            self._mol_buf, self._mol_req_buf, np.int32(self._nmol),
            self._sub_buf, self._sub_req_buf, np.int32(self._nsub),
            npbc, la, lb, lc,
            trans_buf,
            np.int32(int(self.enable_LJ)), np.int32(int(self.enable_Coulomb)),
            np.int32(int(self.enable_HBond)), np.int32(int(self.enable_Morse)),
            np.int32(0), macro_P, macro_AB,
            np.float32(self.Coulomb_const), np.float32(self.Morse_alpha),
            self.spring_k, self.relax_dt, np.int32(self.relax_nsteps),
            local_mem,
            res_total, res_lj, res_coul, res_hb,
            relax_buf,
        ).wait()
        out = {}
        for name, buf in [('total', res_total), ('LJ', res_lj), ('Coulomb', res_coul), ('HBond', res_hb)]:
            arr = np.empty(nconf, dtype=np.float32)
            cl.enqueue_copy(self.queue, arr, buf)
            out[name] = arr
        rpos = np.empty((nconf, self._nmol, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, rpos, relax_buf)
        out['relaxed_pos'] = rpos[:, :, :3]
        return self._apply_macro_correction(transforms, out, relaxed_pos=out['relaxed_pos'])

    def evaluate_single(self, pos=(0,0,3), R=None):
        """Evaluate energy for a single translation (identity rotation or given R).
        Returns dict with scalar energies."""
        Rm = R if R is not None else np.eye(3)
        t = ScanUtils.pack_transforms([Rm], [np.array(pos, dtype=np.float64)])
        res = self.evaluate(t)
        return {k: float(v[0]) for k, v in res.items()}

    # ======== Convenience scan methods ========

    def scan_z(self, pos_xy=(0,0), z_range=(1.5, 8.0), nz=50, R=None, relax=False):
        """1D z-approach scan."""
        transforms, info = ScanUtils.scan_z_approach(pos_xy, z_range, R=R, nz=nz)
        results = self.evaluate_relaxed(transforms) if relax else self.evaluate(transforms)
        info['transforms'] = transforms
        results['scan_info'] = info
        return results

    def scan_lateral(self, z=3.0, x_range=(0,8), y_range=(0,8), nx=50, ny=50, R=None, relax=False):
        """2D lateral (x,y) scan at fixed z."""
        transforms, info = ScanUtils.scan_lateral_2d(z, x_range, y_range, R=R, nx=nx, ny=ny)
        results = self.evaluate_relaxed(transforms) if relax else self.evaluate(transforms)
        info['transforms'] = transforms
        results['scan_info'] = info
        return results

    def scan_xz(self, y=0.0, x_range=(0,8), z_range=(1.5,8), nx=50, nz=50, R=None, relax=False):
        """2D vertical slice (x,z) scan."""
        transforms, info = ScanUtils.scan_xz_slice(y, x_range, z_range, R=R, nx=nx, nz=nz)
        results = self.evaluate_relaxed(transforms) if relax else self.evaluate(transforms)
        info['transforms'] = transforms
        results['scan_info'] = info
        return results

    def scan_rotation(self, pos=(0,0,3), axis=(0,0,1), angle_range=(0, 2*np.pi), nrot=36, relax=False):
        """1D rotation scan at fixed position."""
        transforms, info = ScanUtils.scan_rotation_1d(pos, axis, angle_range, nrot=nrot)
        results = self.evaluate_relaxed(transforms) if relax else self.evaluate(transforms)
        info['transforms'] = transforms
        results['scan_info'] = info
        return results

    def scan_rot_z(self, pos_xy=(0,0), z_range=(1.5,8), axis=(0,0,1), angle_range=(0,2*np.pi), nz=30, nrot=36, relax=False):
        """2D rotation-vs-z scan."""
        transforms, info = ScanUtils.scan_rotation_z_2d(pos_xy, z_range, axis, angle_range, nz=nz, nrot=nrot)
        results = self.evaluate_relaxed(transforms) if relax else self.evaluate(transforms)
        info['transforms'] = transforms
        results['scan_info'] = info
        return results

    def scan_slerp(self, q0, q1, t0, t1, npts=50, relax=False):
        """1D SLERP path scan."""
        transforms, info = ScanUtils.scan_slerp_path(q0, q1, t0, t1, npts=npts)
        results = self.evaluate_relaxed(transforms) if relax else self.evaluate(transforms)
        info['transforms'] = transforms
        results['scan_info'] = info
        return results

    def scan_random(self, pos_center, pos_spread, nsamples=1000, seed=42, relax=False):
        """Monte Carlo random sampling."""
        transforms, info = ScanUtils.scan_monte_carlo(pos_center, pos_spread, nsamples=nsamples, seed=seed)
        results = self.evaluate_relaxed(transforms) if relax else self.evaluate(transforms)
        info['transforms'] = transforms
        results['scan_info'] = info
        return results

    def scan_custom(self, dof_specs, relax=False):
        """Multi-DOF Cartesian product scan."""
        transforms, info = ScanUtils.scan_multi_dof(dof_specs)
        results = self.evaluate_relaxed(transforms) if relax else self.evaluate(transforms)
        info['transforms'] = transforms
        results['scan_info'] = info
        return results
