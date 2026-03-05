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
        macro_P  = cl_array.vec.make_float4(float(self.macro_P[0]), float(self.macro_P[1]), float(self.macro_P[2]), 0.0)
        macro_AB = cl_array.vec.make_float4(float(self.macro_AB[0]), float(self.macro_AB[1]), 0.0, 0.0)
        self.krn_evaluate(
            self.queue, global_size, local_size,
            self._mol_buf, self._mol_req_buf, np.int32(self._nmol),
            self._sub_buf, self._sub_req_buf, np.int32(self._nsub),
            npbc, la, lb, lc,
            trans_buf,
            np.int32(int(self.enable_LJ)), np.int32(int(self.enable_Coulomb)),
            np.int32(int(self.enable_HBond)), np.int32(int(self.enable_Morse)),
            np.int32(int(self.enable_macro)), macro_P, macro_AB,
            np.float32(self.Coulomb_const), np.float32(self.Morse_alpha),
            local_mem,
            res_total, res_lj, res_coul, res_hb,
        ).wait()
        out = {}
        for name, buf in [('total', res_total), ('LJ', res_lj), ('Coulomb', res_coul), ('HBond', res_hb)]:
            arr = np.empty(nconf, dtype=np.float32)
            cl.enqueue_copy(self.queue, arr, buf)
            out[name] = arr
        return out

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
        macro_P  = cl_array.vec.make_float4(float(self.macro_P[0]), float(self.macro_P[1]), float(self.macro_P[2]), 0.0)
        macro_AB = cl_array.vec.make_float4(float(self.macro_AB[0]), float(self.macro_AB[1]), 0.0, 0.0)
        self.krn_relax(
            self.queue, global_size, local_size,
            self._mol_buf, self._mol_req_buf, np.int32(self._nmol),
            self._sub_buf, self._sub_req_buf, np.int32(self._nsub),
            npbc, la, lb, lc,
            trans_buf,
            np.int32(int(self.enable_LJ)), np.int32(int(self.enable_Coulomb)),
            np.int32(int(self.enable_HBond)), np.int32(int(self.enable_Morse)),
            np.int32(int(self.enable_macro)), macro_P, macro_AB,
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
        return out

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
