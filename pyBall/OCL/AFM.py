"""
AFMulator: OpenCL-Accelerated AFM Simulation
============================================

Purpose:
--------
AFMulator simulates atomic force microscopy (AFM) imaging using OpenCL GPU acceleration.
It computes tip-sample interactions (van der Waals, electrostatic, Pauli repulsion) and
relaxes the tip to generate constant-force or constant-height images. Used for STM/AFM
image prediction and interpretation.

Major Functionality:
-------------------
1. Force Field Computation
   - make_forcefield(): Compute LJ/Morse force field on 3D grid
   - realloc_forcefield_buffers(): Reallocate atom and parameter buffers
   - Supports both Lennard-Jones and Morse potentials
   - Includes electrostatic convolution with tip charge distribution
   - Dispersion correction via C6/R^6 integration

2. Tip Relaxation
   - run_scan(): Relax tip at each scan point to constant-force equilibrium
   - realloc_scan_buffers(): Reallocate scan point and FE buffers
   - relaxStrokesTilted(): GPU kernel for tip relaxation with tilt
   - get_raw_FE(): Get raw force/energy without relaxation
   - Supports CO-functionalized tip with quadrupole charge distribution

3. Molecule Loading
   - load_molecule(): Load molecular structure from XYZ file
   - assign_params(): Assign LJ/Morse parameters from ElementTypes.dat
   - Combination rules for tip-sample parameters

4. Scan Management
   - Generate scan grid over sample surface
   - Compute force/energy at each scan point
   - Output AFM images (constant-force or constant-height)

Optimization Policy:
-------------------
- Load kernels once during initialization (load_program in __init__)
- Use persistent GPU buffers via realloc_forcefield_buffers and realloc_scan_buffers
- Add bAlloc guards in calling functions (make_forcefield, run_scan) to skip allocation
- Default bAlloc=True for safety, set False for hot paths with fixed buffer sizes
- See OpenCLBase.py for full policy details
"""

import numpy as np
import pyopencl as cl
import os, sys
from pyBall.globals import debug_print
from .OpenCLBase import OpenCLBase

COULOMB_CONST = 14.3996448915  # [eV*Ang/e^2]


def _bytes_to_gb(nbytes):
    return nbytes / (1024.0**3)

class AFMulator(OpenCLBase):
    """
    PyOpenCL AFM simulator (Phase 1: LJ/Morse + point charges).
    Reuses evalLJC_QZs_toImg / evalMorseC_QZs_toImg + relaxStrokesTilted from relax.cl.
    Phase 2: solve_QEq() numpy QEq charge equilibration.
    Coordinate convention: molecule is shifted so grid occupies [0,L] in all dims,
    making dot(pos, dinv) in [0,1] for normalized image sampling.
    """
    # Default CO-tip parameters (from OCL_PP.h)
    DEFAULT_tipA       = np.array([ 1., 0., 0., 0.], dtype=np.float32)
    DEFAULT_tipB       = np.array([ 0., 1., 0., 0.], dtype=np.float32)
    DEFAULT_tipC       = np.array([ 0., 0., 1.,-0.1], dtype=np.float32)  # .w=dtip
    DEFAULT_stiffness  = np.array([-0.03,-0.03,-0.03,-1.0], dtype=np.float32)
    DEFAULT_dpos0      = np.array([ 0., 0.,-4.0, 4.0], dtype=np.float32)
    DEFAULT_relax_pars = np.array([ 0.5, 0.1, 0.02, 0.5], dtype=np.float32)
    DEFAULT_surfFF     = np.array([ 0., 0., 0., 0.], dtype=np.float32)
    DEFAULT_tipQs      = np.array([ 0.,-0.1, 0.1, 0.], dtype=np.float32)
    DEFAULT_tipQZs     = np.array([ 0., 1.8, 3.6, 0.], dtype=np.float32)


    def __init__(self, cl_src_dir=None, use_morse=False, nloc=32, use_fire=True):
        super().__init__(nloc=nloc, preferred_vendor='nvidia', bPrint=True)
        self.use_morse = use_morse
        self._vram_bytes = 0
        if cl_src_dir is None:
            d = os.path.dirname(os.path.abspath(__file__))
            cl_src_dir = os.path.realpath(os.path.join(d,'..','..','cpp','common_resources','cl'))
        self.cl_src_dir = cl_src_dir
        print(f"AFMulator: cl_src_dir={cl_src_dir}")
        dev = self.ctx.devices[0]
        self._max_alloc = dev.get_info(cl.device_info.MAX_MEM_ALLOC_SIZE)
        self._global_mem = dev.get_info(cl.device_info.GLOBAL_MEM_SIZE)
        print(f"AFMulator: device max_alloc={_bytes_to_gb(self._max_alloc):.3f} GB global_mem={_bytes_to_gb(self._global_mem):.3f} GB")
        relax_cl = os.path.join(cl_src_dir, 'relax.cl')
        print(f"AFMulator: compiling {relax_cl}")
        # Build options: -DOPT_FIRE=0 for damped velocity (matches CPU), -DOPT_FIRE=1 for FIRE
        build_options = ['-D', f'OPT_FIRE={1 if use_fire else 0}']
        self.load_program(kernel_path=relax_cl, build_options=build_options)
        print(f"AFMulator: relax.cl compiled OK (use_fire={use_fire})")
        # State
        self.mol = self.elem_types = None
        self.atoms_arr = self.cLJs_arr = None
        self.atoms_cl  = self.cLJs_cl  = None
        self.img_FF = self.n = self.p0 = self.L = None
        self.dA = self.dB = self.dC = None
        self.dinvA = self.dinvB = self.dinvC = None
        self.mol_shift = np.zeros(3, dtype=np.float32)
        # Tip params (mutable)
        self.tipA       = self.DEFAULT_tipA.copy()
        self.tipB       = self.DEFAULT_tipB.copy()
        self.tipC       = self.DEFAULT_tipC.copy()
        self.stiffness  = self.DEFAULT_stiffness.copy()
        self.dpos0      = self.DEFAULT_dpos0.copy()
        self.relax_pars = self.DEFAULT_relax_pars.copy()
        self.surfFF     = self.DEFAULT_surfFF.copy()
        self.tipQs      = self.DEFAULT_tipQs.copy()
        self.tipQZs     = self.DEFAULT_tipQZs.copy()
        # FDBM scan state (set by setup_fdbm_grid)
        self.img_FF_fdbm  = None   # cl.Image holding F_total (Fx,Fy,Fz,E)
        self.fdbm_dinvA   = None   # float4 for interpFE normalized coords
        self.fdbm_dinvB   = None
        self.fdbm_dinvC   = None
        self.fdbm_origin  = None   # (3,) grid origin [Ang]
        self.fdbm_step    = None   # grid spacing [Ang]
        self.fdbm_shape   = None   # (nx,ny,nz)

    # ── buffer management ─────────────────────────────────────────────────────

    def realloc_forcefield_buffers(self, na):
        """(Re-)allocate persistent GPU buffers for make_forcefield via try_make_buffers."""
        sz_f = 4
        clJ_cols = 4 if self.use_morse else 2
        buffs = {
            "atoms": sz_f * 4 * na,       # float4 per atom (pos + charge)
            "cLJs":  sz_f * clJ_cols * na, # LJ or Morse params per atom
        }
        self.try_make_buffers(buffs, suffix="_cl")

    def realloc_scan_buffers(self, n_scan, nz):
        """(Re-)allocate persistent GPU buffers for run_scan / get_raw_FE via try_make_buffers."""
        sz_f = 4
        buffs = {
            "scan_pts": sz_f * 4 * n_scan,         # float4 per scan point
            "scan_FEs": sz_f * 4 * n_scan * nz,    # float4 per (scan_point, z_step)
        }
        self.try_make_buffers(buffs, suffix="_cl")

    def setup_fdbm_grid(self, F_total, origin, step):
        """
        Upload FDBM force-field image to GPU and precompute lattice vectors.
        Call once per new force field; then call scan_fdbm() repeatedly.

        Args:
            F_total: (nx,ny,nz,4) float32  columns=(Fx,Fy,Fz,E), F = -grad(E_total)
            origin:  (3,) float  grid origin [Ang]
            step:    float  grid spacing [Ang]
        """
        nx, ny, nz = F_total.shape[:3]
        Lx, Ly, Lz = nx*step, ny*step, nz*step
        # Reorder (nx,ny,nz,4) -> (nz,ny,nx,4) for OpenCL 3D image layout
        F_img = np.ascontiguousarray(F_total.transpose(2, 1, 0, 3), dtype=np.float32)
        mf  = cl.mem_flags
        fmt = cl.ImageFormat(cl.channel_order.RGBA, cl.channel_type.FLOAT)
        self.img_FF_fdbm = cl.Image(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, fmt, shape=(nx, ny, nz), hostbuf=F_img)
        # Coordinate transform: coord = (pos - origin) / L = pos/L - origin/L
        # dinvA/B/C dot pos gives normalized coord 0..1 within grid
        self.fdbm_dinvA  = np.array([1./Lx, 0.,    0.,    -origin[0]/Lx], dtype=np.float32)
        self.fdbm_dinvB  = np.array([0.,    1./Ly, 0.,    -origin[1]/Ly], dtype=np.float32)
        self.fdbm_dinvC  = np.array([0.,    0.,    1./Lz, -origin[2]/Lz], dtype=np.float32)
        self.fdbm_origin = np.asarray(origin, dtype=np.float32)
        self.fdbm_step   = float(step)
        self.fdbm_shape  = (nx, ny, nz)
        print(f"AFMulator.setup_fdbm_grid: grid={nx}x{ny}x{nz}  step={step:.3f}  L=({Lx:.1f},{Ly:.1f},{Lz:.1f}) Ang")

    def scan_fdbm(self, scan_xs, scan_ys, probe_heights, mol_z=0.0,
                  ppm_mode=True, use_fire=True,
                  K_LAT=0.0031, K_RAD=0.1248, bond_length=3.0,
                  stiffness=None, dpos0=None, relax_pars=None):
        """
        GPU probe-particle relaxation over a 2D scan using relaxStrokes kernel.
        Requires setup_fdbm_grid() to be called first.

        Two tip models (ppm_mode switch):
          ppm_mode=True  (default, physical PPM):
            Probe on a radial bond of length L below tip-apex.
            dpos0 = (0, 0, -bond_length, bond_length)
            stiffness = (-K_LAT, -K_LAT, -K_LAT, -K_RAD)
            Scan heights are tip-apex heights; probe sits ~bond_length below.
            Typical: K_LAT=0.0031 eV/Ang^2 (0.5 N/m), K_RAD=0.1248 eV/Ang^2 (20 N/m), L=3 Ang.

          ppm_mode=False (linear harmonic):
            Pure 3D harmonic spring, no radial bond (stiffness.w=0).
            dpos0 = (0, 0, 0, 0)  probe at tip position.
            stiffness = (-K_LAT, -K_LAT, -K_LAT, 0)  isotropic lateral spring.
            Scan heights are probe heights directly.

        Integrator modes (use_fire switch):
          use_fire=True  (default): Use FIRE optimizer (fast convergence)
          use_fire=False: Use damped velocity (matches CPU pp_relax_2d behavior)
            CPU uses: v = 0.8*v + 0.3*f, pos += v*0.3
            GPU equivalent: damp=0.2, dt=0.3 (for 50 steps N_RELAX_STEP_MAX)

        Args:
            scan_xs:       (nx_s,) x scan positions [Ang, world coords]
            scan_ys:       (ny_s,) y scan positions [Ang]
            probe_heights: (nz_s,) tip-apex heights above mol_z [Ang]
            mol_z:         float  z of molecule top [Ang]
            ppm_mode:      bool  True=PPM radial bond, False=simple harmonic
            use_fire:      bool  True=FIRE integrator, False=damped velocity
            K_LAT:         float  lateral stiffness [eV/Ang^2]  (ppm_mode=True)
            K_RAD:         float  radial bond stiffness [eV/Ang^2]  (ppm_mode=True)
            bond_length:   float  CO bond length [Ang]  (ppm_mode=True)
            stiffness:     (4,) float4 override; None -> auto from mode
            dpos0:         (4,) float4 override; None -> auto from mode
            relax_pars:    (4,) float4 (dt,damp,tmin,tmax); None -> auto from mode/use_fire

        Returns:
            FEs_relax: (nx_s, ny_s, nz_s, 4) float32  relaxed (Fx,Fy,Fz,E) at probe pos
                       iz=0 corresponds to probe_heights[0] (lowest z)
        """
        assert self.img_FF_fdbm is not None, "call setup_fdbm_grid() first"
        nx_s, ny_s, nz_s = len(scan_xs), len(scan_ys), len(probe_heights)
        n_scan = nx_s * ny_s
        origin = self.fdbm_origin

        # --- tip model parameters ---
        if ppm_mode:
            # Physical PPM: probe hangs bond_length below tip-apex.
            # probe_heights are the desired probe positions above mol_z.
            # Tip-apex must be bond_length higher so probe lands at probe_heights.
            stiffness_v  = np.array([-K_LAT, -K_LAT, -K_LAT, -K_RAD], dtype=np.float32) if stiffness is None else np.array(stiffness, dtype=np.float32)
            dpos0_v      = np.array([0., 0., -bond_length, bond_length], dtype=np.float32) if dpos0 is None else np.array(dpos0, dtype=np.float32)
            z_start = float(np.max(probe_heights)) + mol_z + bond_length  # tip-apex = probe_z + L
        else:
            # Linear mode: pure 3D harmonic spring, no radial bond (stiffness.w=0)
            # dpos0.z=-0.001 avoids r=0 in tipForce (prevents 0/0=NaN)
            stiffness_v  = np.array([-K_LAT, -K_LAT, -K_LAT, 0.], dtype=np.float32) if stiffness is None else np.array(stiffness, dtype=np.float32)
            dpos0_v      = np.array([0., 0., -0.001, 0.], dtype=np.float32) if dpos0 is None else np.array(dpos0, dtype=np.float32)
            z_start = float(np.max(probe_heights)) + mol_z

        # --- integrator parameters ---
        if relax_pars is not None:
            relax_pars_v = np.array(relax_pars, dtype=np.float32)
        elif use_fire:
            relax_pars_v = self.relax_pars.copy()  # FIRE: default (dt=0.5, damp=0.1)
        else:
            # Damped velocity: match CPU pp_relax_2d (v=0.8*v+0.3*f, pos+=v*0.3)
            # OpenCL: v*=(1-damp), v+=f*dt, pos+=v*dt -> damp=0.2, dt=0.3
            relax_pars_v = np.array([0.3, 0.2, 0.03, 0.3], dtype=np.float32)

        # --- build scan point buffer (grid-local coordinates) ---
        # Tip starts at highest z, steps down
        heights_desc = np.sort(probe_heights)[::-1]   # descending
        dh = float(heights_desc[0] - heights_desc[-1]) / max(nz_s - 1, 1) if nz_s > 1 else 0.
        dTip = np.array([0., 0., -dh, 0.], dtype=np.float32)

        # vectorised: outer product of scan_xs, scan_ys -> (nx_s*ny_s, 4)
        # pts must be world coordinates (Ang); interpFE applies the full dinv transform internally
        gx, gy = np.meshgrid(scan_xs, scan_ys, indexing='ij')
        pts = np.zeros((n_scan, 4), dtype=np.float32)
        pts[:, 0] = gx.ravel()
        pts[:, 1] = gy.ravel()
        pts[:, 2] = z_start

        # --- allocate/reuse GPU scan buffers ---
        self.realloc_scan_buffers(n_scan, nz_s)
        self.toGPU_(self.scan_pts_cl, pts)

        # --- launch relaxStrokes ---
        self.prg.relaxStrokes(
            self.queue, (max(n_scan, 1),), (1,),
            self.img_FF_fdbm,
            self.scan_pts_cl,
            self.scan_FEs_cl,
            self.fdbm_dinvA, self.fdbm_dinvB, self.fdbm_dinvC,
            dTip,
            stiffness_v, dpos0_v, relax_pars_v,
            np.int32(nz_s)
        )
        self.queue.finish()

        # --- download and reshape ---
        FEs_h = np.zeros((n_scan * nz_s, 4), dtype=np.float32)
        self.fromGPU_(self.scan_FEs_cl, FEs_h)
        self.queue.finish()

        # relaxStrokes: FEs[gid*nz+iz] where iz=0=highest z; flip to iz=0=lowest
        FEs_relax = FEs_h.reshape(nx_s, ny_s, nz_s, 4)[:, :, ::-1, :]
        Fz = FEs_relax[:,:,:,2]
        print(f"  AFMulator.scan_fdbm: Fz min={Fz.min():.4f}  max={Fz.max():.4f}  mean={Fz.mean():.4f} eV/Ang  ppm_mode={ppm_mode}")
        return FEs_relax

    def scan_fdbm_2d(self, scan_xs, scan_ys, probe_heights, mol_z=0.0, K_LAT=0.5, dt=0.3, damp=0.2):
        """
        2D lateral-only relaxation using relaxStrokes2D kernel.
        Exactly matches CPU pp_relax_2d: z fixed per slice, only x,y relax with damped MD.
        Damped update: v *= (1-damp); v += F*dt; pos += v*dt
        Probe resets to anchor at start of each height slice (same as CPU).

        Args:
            scan_xs:       (nx_s,) x scan positions [Ang, world coords]
            scan_ys:       (ny_s,) y scan positions [Ang]
            probe_heights: (nz_s,) probe heights above mol_z, ascending [Ang]
            mol_z:         float  z of molecule top [Ang]
            K_LAT:         float  lateral stiffness [eV/Ang^2]
            dt:            float  time step
            damp:          float  damping coefficient (v *= 1-damp)

        Returns:
            FEs_relax: (nx_s, ny_s, nz_s, 4) float32  (Fx,Fy,Fz,E) at relaxed pos
                       iz=0 = probe_heights[0] (lowest z)
        """
        assert self.img_FF_fdbm is not None, "call setup_fdbm_grid() first"
        nx_s, ny_s, nz_s = len(scan_xs), len(scan_ys), len(probe_heights)
        n_scan = nx_s * ny_s

        # heights descending: iz=0 = highest, iz=nz-1 = lowest
        heights_desc = np.sort(probe_heights)[::-1]
        dh = float(heights_desc[0] - heights_desc[-1]) / max(nz_s - 1, 1) if nz_s > 1 else 0.
        z_start = float(heights_desc[0]) + mol_z   # world z of first (highest) slice

        gx, gy = np.meshgrid(scan_xs, scan_ys, indexing='ij')
        pts = np.zeros((n_scan, 4), dtype=np.float32)
        pts[:, 0] = gx.ravel()
        pts[:, 1] = gy.ravel()
        pts[:, 2] = z_start

        print( "!!!!! scan_fdbm_2d() !!!!!" )

        self.realloc_scan_buffers(n_scan, nz_s)
        self.toGPU_(self.scan_pts_cl, pts)
        self.prg.relaxStrokes2D(
            self.queue, (n_scan,), (1,),
            self.img_FF_fdbm,
            self.scan_pts_cl,
            self.scan_FEs_cl,
            self.fdbm_dinvA, self.fdbm_dinvB, self.fdbm_dinvC,
            np.float32(K_LAT), np.float32(dh),
            np.float32(dt),    np.float32(damp),
            np.int32(nz_s)
        )
        self.queue.finish()
        FEs_h = np.zeros((n_scan * nz_s, 4), dtype=np.float32)
        self.fromGPU_(self.scan_FEs_cl, FEs_h)
        self.queue.finish()
        # flip iz so iz=0 = lowest z (matches CPU pp_relax_2d order)
        FEs_relax = FEs_h.reshape(nx_s, ny_s, nz_s, 4)[:, :, ::-1, :]
        Fz = FEs_relax[:,:,:,2]
        print(f"  AFMulator.scan_fdbm_2d: Fz min={Fz.min():.4f}  max={Fz.max():.4f}  mean={Fz.mean():.4f} eV/Ang")
        return FEs_relax

    def realloc_dispersion_buffers(self, natoms):
        """(Re-)allocate persistent GPU buffers for dispersion computation via try_make_buffers."""
        sz_f = 4
        buffs = {
            "disp_atoms": sz_f * 4 * natoms,      # float4 per atom (pos + charge)
            "disp_C6": sz_f * 2 * natoms,          # float2 per atom (C6_eff, 0)
        }
        self.try_make_buffers(buffs, suffix="_cl")

    # ── molecule loading ──────────────────────────────────────────────────────

    def load_molecule(self, xyz_path):
        sys.path.insert(0, os.path.join(os.path.dirname(__file__),'..','..'))
        from pyBall.AtomicSystem import AtomicSystem
        self.mol = AtomicSystem(fname=xyz_path)
        print(f"AFMulator.load_molecule: {self.mol.natoms} atoms  enames[:5]={self.mol.enames[:5]}")
        return self.mol

    def assign_params(self, params_path=None, tip_R=1.452, tip_E=0.0006808, tip_alpha=-1.8):
        """
        Assign LJ (c6,c12) or Morse (R0,E0,alpha) params from ElementTypes.dat.
        Combination rules: R0_ij = tip_R + R_sample,  E0_ij = sqrt(tip_E * E_sample).
        LJ: c6 = 2*E0*R0^6, c12 = E0*R0^12  (matches getLJ in relax.cl).
        """
        assert self.mol is not None, "call load_molecule() first"
        if params_path is None:
            params_path = os.path.realpath(
                os.path.join(self.cl_src_dir,'..','..','ElementTypes.dat'))
        from .MMparams import read_element_types
        et = read_element_types(params_path)
        self.elem_types = et
        na = self.mol.natoms
        atoms = np.zeros((na,4), dtype=np.float32)
        atoms[:,:3] = self.mol.apos[:na,:3].astype(np.float32)
        if self.mol.qs is not None:
            atoms[:,3] = self.mol.qs[:na].astype(np.float32)
        if self.use_morse:
            cMs = np.zeros((na,4), dtype=np.float32)
            for ia,ename in enumerate(self.mol.enames[:na]):
                e = et.get(ename)
                assert e is not None, f"Element '{ename}' not in ElementTypes.dat"
                R0 = tip_R + e.RvdW
                E0 = np.sqrt(abs(tip_E * e.EvdW))
                cMs[ia] = [R0, E0, tip_alpha, 0.]
            self.cLJs_arr = cMs
        else:
            cLJs = np.zeros((na,2), dtype=np.float32)
            for ia,ename in enumerate(self.mol.enames[:na]):
                e = et.get(ename)
                assert e is not None, f"Element '{ename}' not in ElementTypes.dat"
                R0 = tip_R + e.RvdW
                E0 = np.sqrt(abs(tip_E * e.EvdW))
                cLJs[ia] = [2.*E0*R0**6, E0*R0**12]
            self.cLJs_arr = cLJs
        self.atoms_arr = atoms
        print(f"AFMulator.assign_params: {'Morse' if self.use_morse else 'LJ'} for {na} atoms")
        return atoms, self.cLJs_arr

    # ── grid setup ────────────────────────────────────────────────────────────

    def setup_grid(self, n=(100,100,60), L=None, margin=3.0, z_top=12.0):
        """
        Set up 3D force-field grid. Molecule is shifted so grid occupies [0,L] in all dims,
        making dot(pos, dinv) in [0,1] for normalized image sampling.
        n: (nx,ny,nz)  L: (Lx,Ly,Lz) Ang  margin: xy pad  z_top: extra z above molecule.
        """
        assert self.atoms_arr is not None, "call assign_params() first"
        apos = self.atoms_arr[:,:3]
        mn, mx = apos.min(axis=0), apos.max(axis=0)
        if L is None:
            L = np.array([(mx[0]-mn[0])+2*margin,
                          (mx[1]-mn[1])+2*margin,
                          (mx[2]-mn[2])+margin+z_top], dtype=np.float32)
        else:
            L = np.array(L, dtype=np.float32)
        self.n = np.array(n, dtype=np.int32)
        nx,ny,nz = n
        # Shift atoms so grid origin is at (0,0,0)
        p0_raw = np.array([mn[0]-margin, mn[1]-margin, mn[2]-margin/2], dtype=np.float32)
        self.mol_shift = -p0_raw
        self.atoms_arr[:,:3] += self.mol_shift
        self.p0    = np.zeros(4, dtype=np.float32)          # grid origin = (0,0,0)
        self.dA    = np.array([L[0]/nx, 0., 0., 0.], dtype=np.float32)  # per-voxel step
        self.dB    = np.array([0., L[1]/ny, 0., 0.], dtype=np.float32)
        self.dC    = np.array([0., 0., L[2]/nz, 0.], dtype=np.float32)
        self.dinvA = np.array([1./L[0], 0., 0., 0.], dtype=np.float32)  # for sampler_1
        self.dinvB = np.array([0., 1./L[1], 0., 0.], dtype=np.float32)
        self.dinvC = np.array([0., 0., 1./L[2], 0.], dtype=np.float32)
        self.L = L
        grid_bytes = int(nx*ny*nz*4*4)  # float4 image
        if not hasattr(self, '_vram_bytes'): self._vram_bytes = 0
        warn = "" if grid_bytes <= getattr(self, '_max_alloc', grid_bytes+1) else " [exceeds device max_alloc!]"
        print(f"AFMulator.setup_grid: n={n} L={L} mol_shift={self.mol_shift} | FF image ~{_bytes_to_gb(grid_bytes):.3f} GB{warn}")
        self._grid_bytes = grid_bytes
        return L

    def setup_grid_lvec(self, n=(100, 100, 60), margin_z=2.0, z_top=12.0):
        """Set up primitive-cell force-field grid aligned with mol.lvec.

        - x,y are along mol.lvec[0], mol.lvec[1] (skew cells supported)
        - z is Cartesian along mol.lvec[2] if present, otherwise +z axis
        - grid uses normalized image sampling with periodic wrap (sampler uses repeat/mirrored-repeat)

        This expects atoms_arr to already be in the same coordinate frame as mol.lvec.
        """
        assert self.atoms_arr is not None, "call assign_params() first"
        assert self.mol is not None and getattr(self.mol, 'lvec', None) is not None, \
            "setup_grid_lvec requires mol.lvec (periodic system)"

        nx, ny, nz = [int(x) for x in n]
        if nx <= 1 or ny <= 1 or nz <= 1:
            raise ValueError(f"setup_grid_lvec: invalid n={n}")

        a = self.mol.lvec[0, :3].astype(np.float32)
        b = self.mol.lvec[1, :3].astype(np.float32)
        c = self.mol.lvec[2, :3].astype(np.float32) if self.mol.lvec.shape[0] > 2 else np.array([0., 0., 1.], dtype=np.float32)

        # z extent from atoms + padding
        apos = self.atoms_arr[:, :3]
        zmin = float(apos[:, 2].min())
        zmax = float(apos[:, 2].max())
        Lz = (zmax - zmin) + float(margin_z) + float(z_top)
        if Lz <= 0:
            raise ValueError(f"setup_grid_lvec: invalid Lz={Lz}")

        # Shift atoms so that bottom is near z=margin_z*0.5, keep x,y in cell frame (no bbox shift)
        p0z = zmin - 0.5*float(margin_z)
        self.mol_shift = np.array([0., 0., -p0z], dtype=np.float32)
        self.atoms_arr[:, :3] += self.mol_shift[None, :]

        self.n = np.array([nx, ny, nz], dtype=np.int32)
        self.p0 = np.zeros(4, dtype=np.float32)
        self.dA = np.array([a[0]/nx, a[1]/nx, a[2]/nx, 0.], dtype=np.float32)
        self.dB = np.array([b[0]/ny, b[1]/ny, b[2]/ny, 0.], dtype=np.float32)
        self.dC = np.array([0., 0., Lz/nz, 0.], dtype=np.float32)

        # Inverse mapping for normalized coordinates u,v,w in [0,1]
        # pos = a*u + b*v + (0,0,Lz*w)
        M = np.stack([a, b, np.array([0., 0., Lz], dtype=np.float32)], axis=1).astype(np.float64)  # 3x3
        det = float(np.linalg.det(M))
        if abs(det) < 1e-12:
            raise ValueError(f"setup_grid_lvec: singular cell matrix det={det} for a={a} b={b} Lz={Lz}")
        invMT = np.linalg.inv(M).T.astype(np.float32)
        self.dinvA = np.array([invMT[0, 0], invMT[0, 1], invMT[0, 2], 0.], dtype=np.float32)
        self.dinvB = np.array([invMT[1, 0], invMT[1, 1], invMT[1, 2], 0.], dtype=np.float32)
        self.dinvC = np.array([invMT[2, 0], invMT[2, 1], invMT[2, 2], 0.], dtype=np.float32)

        self.L = np.array([np.linalg.norm(a), np.linalg.norm(b), Lz], dtype=np.float32)
        grid_bytes = int(nx * ny * nz * 4 * 4)
        warn = "" if grid_bytes <= getattr(self, '_max_alloc', grid_bytes+1) else " [exceeds device max_alloc!]"
        print(f"AFMulator.setup_grid_lvec: n=({nx},{ny},{nz}) a={a} b={b} Lz={Lz:.3f} mol_shift={self.mol_shift} | FF image ~{_bytes_to_gb(grid_bytes):.3f} GB{warn}")
        self._grid_bytes = grid_bytes
        return self.L

    # ── force field ───────────────────────────────────────────────────────────

    def make_forcefield(self, bAlloc=True):
        """Upload atom data to GPU and run evalLJC_QZs_toImg (or evalMorseC_QZs_toImg)."""
        assert self.n is not None, "call setup_grid() first"
        nx,ny,nz = self.n
        na = len(self.atoms_arr)
        nMax = int(nx*ny*nz)
        mf = cl.mem_flags
        img_bytes = nx*ny*nz*4*4
        print(f"AFMulator.make_forcefield: na={na} grid=({nx},{ny},{nz}) img={_bytes_to_gb(img_bytes):.3f} GB")
        if img_bytes > self._max_alloc:
            raise MemoryError(f"FF image needs {img_bytes} bytes ({_bytes_to_gb(img_bytes):.3f} GB) > device max_alloc {_bytes_to_gb(self._max_alloc):.3f} GB")
        # 3D image dimension guard
        max_w = self.ctx.devices[0].get_info(cl.device_info.IMAGE3D_MAX_WIDTH)
        max_h = self.ctx.devices[0].get_info(cl.device_info.IMAGE3D_MAX_HEIGHT)
        max_d = self.ctx.devices[0].get_info(cl.device_info.IMAGE3D_MAX_DEPTH)
        if nx > max_w or ny > max_h or nz > max_d:
            raise MemoryError(f"FF image dims {nx}x{ny}x{nz} exceed device limits {max_w}x{max_h}x{max_d}; reduce spacing/margin/z_top")
        # Persistent atom/cLJ buffers
        if bAlloc:
            self.realloc_forcefield_buffers(na)
        self.toGPU_(self.atoms_cl, self.atoms_arr)
        self.toGPU_(self.cLJs_cl,  self.cLJs_arr)
        # img_FF is a cl.Image (not a Buffer) -- reallocate only when shape changes
        fmt = cl.ImageFormat(cl.channel_order.RGBA, cl.channel_type.FLOAT)
        if self.img_FF is None or getattr(self.img_FF, 'shape', None) != (nx, ny, nz):
            try:
                self.img_FF = cl.Image(self.ctx, mf.READ_WRITE, fmt, shape=(nx,ny,nz))
            except Exception as e:
                print(f"AFMulator.make_forcefield: cl.Image allocation failed for shape {nx}x{ny}x{nz} : {e}")
                raise
        gs  = (self._roundup(nMax, self.nloc),)
        ls  = (self.nloc,)
        nGrid = np.array([nx,ny,nz,nMax], dtype=np.int32)
        print(f"AFMulator.make_forcefield: na={na} nGrid={nGrid[:3]} gs={gs[0]}")
        kernel_name = "evalMorseC_QZs_toImg" if self.use_morse else "evalLJC_QZs_toImg"
        try:
            getattr(self.prg, kernel_name)(
                self.queue, gs, ls,
                np.int32(na),
                self.atoms_cl, self.cLJs_cl,
                self.img_FF,
                nGrid, self.p0,
                self.dA, self.dB, self.dC,
                self.tipQs, self.tipQZs
            )
            self.queue.finish()
        except Exception as e:
            print(f"AFMulator.make_forcefield: kernel {kernel_name} failed: {e}")
            raise
        print("AFMulator.make_forcefield: done")

    # ── scan ──────────────────────────────────────────────────────────────────

    def run_scan(self, nxy=(50,50), nz=60, dtip=-0.1,
                 scan_p0=None, scan_da=None, scan_db=None, bAlloc=True):
        """
        Run relaxStrokesTilted for (nx_scan × ny_scan) scan points.
        Returns FEs (nx_scan, ny_scan, nz, 4) and pts (nx_scan, ny_scan, 3).
        scan_p0/da/db: scan grid in kernel-space (post mol_shift); None = auto from bounding box.
        """
        assert self.img_FF is not None, "call make_forcefield() first"
        nx_s, ny_s = nxy
        n_scan = nx_s * ny_s
        if scan_p0 is None:
            apos = self.atoms_arr[:,:3]
            mn, mx = apos.min(axis=0), apos.max(axis=0)
            x0 = mn[0] + (mx[0]-mn[0])*0.05
            y0 = mn[1] + (mx[1]-mn[1])*0.05
            # tip start = mol_top + clearance + |dpos0_z|  so probe starts clearance above mol
            clearance = 5.0  # Ang above molecule surface where probe begins
            z0 = float(mx[2]) + clearance + abs(float(self.dpos0[2]))
            scan_p0 = np.array([x0, y0, z0], dtype=np.float32)
            scan_da = np.array([(mx[0]-mn[0])*0.9/max(nx_s-1,1), 0., 0.], dtype=np.float32)
            scan_db = np.array([0., (mx[1]-mn[1])*0.9/max(ny_s-1,1), 0.], dtype=np.float32)
        print(f"AFMulator.run_scan: nxy={nxy} nz={nz} dtip={dtip}")
        print(f"  scan_p0={scan_p0}  da={scan_da}  db={scan_db}")
        pts = np.zeros((n_scan,4), dtype=np.float32)
        k = 0
        for ix in range(nx_s):
            for iy in range(ny_s):
                pts[k,:3] = scan_p0 + scan_da*ix + scan_db*iy
                k += 1
        FEs_bytes = n_scan*nz*4*4
        if FEs_bytes > self._max_alloc:
            raise MemoryError(f"run_scan: FEs buffer needs {FEs_bytes} bytes ({_bytes_to_gb(FEs_bytes):.3f} GB) > device max_alloc {_bytes_to_gb(self._max_alloc):.3f} GB")
        if bAlloc:
            self.realloc_scan_buffers(n_scan, nz)
        self.toGPU_(self.scan_pts_cl, pts)
        FEs_h  = np.zeros((n_scan*nz, 4), dtype=np.float32)
        tipC = self.tipC.copy(); tipC[3] = np.float32(dtip)
        gs = (self._roundup(n_scan, 1),)
        ls = (1,)
        self.prg.relaxStrokesTilted(
            self.queue, gs, ls,
            self.img_FF, self.scan_pts_cl, self.scan_FEs_cl,
            self.dinvA, self.dinvB, self.dinvC,
            self.tipA, self.tipB, tipC,
            self.stiffness, self.dpos0, self.relax_pars, self.surfFF,
            np.int32(nz)
        )
        self.queue.finish()
        self.fromGPU_(self.scan_FEs_cl, FEs_h)
        self.queue.finish()
        FEs = FEs_h.reshape(nx_s, ny_s, nz, 4)
        print(f"AFMulator.run_scan: done FEs.shape={FEs.shape}")
        return FEs, pts[:,:3].reshape(nx_s, ny_s, 3)

    # ── raw FF (no PP relaxation) ─────────────────────────────────────────────

    def get_raw_FE(self, nxy=(60,60), nz=21, dtip=-0.2,
                   scan_p0=None, scan_da=None, scan_db=None, bAlloc=True):
        """
        Sample force field WITHOUT probe-particle relaxation using getFEinStrokesTilted.
        Probe moves at fixed offset dpos0 below tip, no spring/FIRE minimisation.
        Returns raw_FEs (nx,ny,nz,4) and pts (nx,ny,3).  Same scan grid as run_scan.
        """
        assert self.img_FF is not None, "call make_forcefield() first"
        nx_s, ny_s = nxy
        n_scan = nx_s * ny_s
        if scan_p0 is None:
            apos = self.atoms_arr[:,:3]
            mn, mx = apos.min(axis=0), apos.max(axis=0)
            x0 = mn[0] + (mx[0]-mn[0])*0.05
            y0 = mn[1] + (mx[1]-mn[1])*0.05
            clearance = 5.0
            z0 = float(mx[2]) + clearance + abs(float(self.dpos0[2]))
            scan_p0 = np.array([x0, y0, z0], dtype=np.float32)
            scan_da = np.array([(mx[0]-mn[0])*0.9/max(nx_s-1,1), 0., 0.], dtype=np.float32)
            scan_db = np.array([0., (mx[1]-mn[1])*0.9/max(ny_s-1,1), 0.], dtype=np.float32)
        print(f"AFMulator.get_raw_FE: nxy={nxy} nz={nz} dtip={dtip}")
        print(f"  scan_p0={scan_p0}  da={scan_da}  db={scan_db}")
        pts = np.zeros((n_scan,4), dtype=np.float32)
        k = 0
        for ix in range(nx_s):
            for iy in range(ny_s):
                pts[k,:3] = scan_p0 + scan_da*ix + scan_db*iy
                k += 1
        FEs_bytes = n_scan*nz*4*4
        if FEs_bytes > self._max_alloc:
            raise MemoryError(f"get_raw_FE: FEs buffer needs {FEs_bytes} bytes ({_bytes_to_gb(FEs_bytes):.3f} GB) > device max_alloc {_bytes_to_gb(self._max_alloc):.3f} GB")
        if bAlloc:
            self.realloc_scan_buffers(n_scan, nz)
        self.toGPU_(self.scan_pts_cl, pts)
        FEs_h  = np.zeros((n_scan*nz, 4), dtype=np.float32)
        dTip = np.array([0., 0., dtip, 0.], dtype=np.float32)
        gs = (self._roundup(n_scan, 1),)
        ls = (1,)
        self.prg.getFEinStrokesTilted(
            self.queue, gs, ls,
            self.img_FF, self.scan_pts_cl, self.scan_FEs_cl,
            self.dinvA, self.dinvB, self.dinvC,
            self.tipA, self.tipB, self.tipC,
            dTip, self.dpos0,
            np.int32(nz)
        )
        self.queue.finish()
        self.fromGPU_(self.scan_FEs_cl, FEs_h)
        self.queue.finish()
        FEs = FEs_h.reshape(nx_s, ny_s, nz, 4)
        print(f"AFMulator.get_raw_FE: done FEs.shape={FEs.shape}")
        return FEs, pts[:,:3].reshape(nx_s, ny_s, 3)

    @property
    def mol_z(self):
        """Molecule plane z in kernel-space (after mol_shift). Useful for computing probe heights."""
        if self.atoms_arr is None:
            return 0.0
        return float(self.atoms_arr[:, 2].max())

    # ── frequency shift ───────────────────────────────────────────────────────

    def get_df(self, FEs, k_tip=1800.0):
        """
        Frequency shift from force traces (Sader-Jarvis simplified):
          df(x,y) ≈ -mean_z(dFz/dz) / (2*k_tip)
        FEs: (nx,ny,nz,4); returns df (nx,ny).
        """
        Fz = FEs[:,:,:,2]
        dFzdz = np.gradient(Fz, axis=2)
        return -dFzdz.mean(axis=2) / (2.0 * k_tip)

    # ── Phase 2: QEq ─────────────────────────────────────────────────────────

    def solve_QEq(self, Q_total=0.0):
        """
        CPU numpy Charge Equilibration (QEq).
        Builds screened Coulomb matrix J_ij = COULOMB_CONST / sqrt(r_ij^2 + Ra_ij^2),
        J_ii = Ehard_i, solves [J,1;1^T,0][q;lam]=[-chi;Q_total].
        Updates atoms_arr[:,3] and mol.qs with computed charges.
        """
        assert self.mol is not None and self.elem_types is not None, \
            "call load_molecule() + assign_params() first"
        na = self.mol.natoms
        apos = self.atoms_arr[:,:3]
        et = self.elem_types
        chi = np.zeros(na); hard = np.zeros(na); Ra = np.zeros(na)
        for ia, ename in enumerate(self.mol.enames[:na]):
            e = et.get(ename)
            assert e is not None and e.bQEq, f"No QEq params for element '{ename}'"
            chi[ia]  = -e.Eaff       # electronegativity
            hard[ia] =  e.Ehard      # chemical hardness
            Ra[ia]   =  e.Ra         # screening radius
        dr = apos[:,None,:] - apos[None,:,:]        # (na,na,3)
        r  = np.linalg.norm(dr, axis=2)             # (na,na)
        Ra_ij = 0.5*(Ra[:,None]+Ra[None,:])
        J = COULOMB_CONST / np.sqrt(r**2 + Ra_ij**2)
        np.fill_diagonal(J, hard)
        A = np.zeros((na+1,na+1)); b = np.zeros(na+1)
        A[:na,:na] = J; A[:na,na] = 1.; A[na,:na] = 1.
        b[:na] = -chi; b[na] = Q_total
        sol = np.linalg.solve(A, b)
        qs = sol[:na].astype(np.float32)
        print(f"AFMulator.solve_QEq: sum={qs.sum():.4f} (target {Q_total}) "
              f"range=[{qs.min():.3f},{qs.max():.3f}]")
        self.atoms_arr[:,3] = qs
        if self.mol.qs is not None:
            self.mol.qs[:na] = qs
        # Invalidate cached GPU buffers (need re-upload after charge change)
        self.atoms_cl = None
        return qs

    # ── helpers ───────────────────────────────────────────────────────────────

    def compute_dispersion_grid_cl(self, atomPos, atomTypes, origin, step, ngrid, C6_atom_dict=None, C6_CO=30.0, RA=1.5, bAlloc=True, return_grads=True):
        """
        Compute C6/r^6 dispersion energy grid using OpenCL GPU acceleration.

        Args:
            atomPos: (natoms, 3) atom positions in Angstrom
            atomTypes: (natoms,) atomic numbers
            origin: (3,) grid origin
            step: grid spacing in Angstrom
            ngrid: (3,) grid dimensions
            C6_atom_dict: dict mapping Z to C6 coefficients (default: {1:6.5, 6:24.0, 7:20.0, 8:15.0})
            C6_CO: C6 coefficient for CO tip
            RA: damping radius in Angstrom
            bAlloc: if True, allocate GPU buffers (default True)
            return_grads: if True, compute and return gradients (default True for backward compat)

        Returns:
            E_vdw: (nx, ny, nz) dispersion energy field
            grads: (nx, ny, nz, 3) dispersion gradients (if return_grads=True)
        """
        if C6_atom_dict is None:
            C6_atom_dict = {1: 6.5, 6: 24.0, 7: 20.0, 8: 15.0}

        nx, ny, nz = [int(i) for i in ngrid[:3]]
        natoms = len(atomPos)

        # Prepare atom data
        atoms_arr = np.zeros((natoms, 4), dtype=np.float32)
        atoms_arr[:, :3] = np.asarray(atomPos, dtype=np.float32)
        atoms_arr[:, 3] = 0.0  # charge not used for dispersion

        # Prepare C6 coefficients: C6_eff = sqrt(C6_atom * C6_CO)
        C6_atom = np.array([C6_atom_dict.get(z, 1.0) for z in atomTypes], dtype=np.float32)
        C6_eff = np.sqrt(C6_atom * C6_CO)
        C6_params = np.zeros((natoms, 2), dtype=np.float32)
        C6_params[:, 0] = C6_eff
        C6_params[:, 1] = 0.0

        # Allocate buffers
        if bAlloc:
            self.realloc_dispersion_buffers(natoms)

        # Upload data to GPU
        self.toGPU_(self.disp_atoms_cl, atoms_arr)
        self.toGPU_(self.disp_C6_cl, C6_params)

        # Prepare grid spec
        origin = np.asarray(origin, dtype=np.float32)
        step = np.asarray(step, dtype=np.float32)
        grid_p0 = np.array([origin[0], origin[1], origin[2], 0.0], dtype=np.float32)
        grid_dA = np.array([step, 0.0, 0.0, 0.0], dtype=np.float32)
        grid_dB = np.array([0.0, step, 0.0, 0.0], dtype=np.float32)
        grid_dC = np.array([0.0, 0.0, step, 0.0], dtype=np.float32)
        nMax = int(nx * ny * nz)
        nGrid = np.array([nx, ny, nz, nMax], dtype=np.int32)
        R2damp = RA * RA

        # Allocate output image
        mf = cl.mem_flags
        fmt = cl.ImageFormat(cl.channel_order.RGBA, cl.channel_type.FLOAT)
        img_disp = cl.Image(self.ctx, mf.READ_WRITE, fmt, shape=(nx, ny, nz))

        # Launch kernel
        gs = (nx * ny * nz,)
        ls = (32,)
        self.prg.evalDispersion_toImg(
            self.queue, gs, ls,
            np.int32(natoms),
            self.disp_atoms_cl,
            self.disp_C6_cl,
            img_disp,
            nGrid,
            grid_p0,
            grid_dA,
            grid_dB,
            grid_dC,
            np.float32(R2damp)
        )

        # Read back energy grid
        # OpenCL 3D images have z varying fastest in memory, so we need shape (nz, ny, nx, 4)
        E_vdw_cl = np.zeros((nz, ny, nx, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, E_vdw_cl, img_disp, origin=(0, 0, 0), region=(nx, ny, nz))
        # Transpose from (nz, ny, nx) to (nx, ny, nz) to match expected layout
        E_vdw = E_vdw_cl.transpose(2, 1, 0, 3)  # (nz, ny, nx, 4) -> (nx, ny, nz, 4)
        E_vdw = E_vdw[:, :, :, 3]  # energy is in .w component

        # Compute gradients using GPU (only if requested)
        if return_grads:
            print("  Computing vdW gradients using GPU...")
            grads_cl = self.compute_gradient_cl(E_vdw, step, bAlloc=False)
            # Extract gradients (negate forces): grad = -F
            grads_E_vdw = -grads_cl[..., :3]
            return E_vdw, grads_E_vdw

        return E_vdw

    def compute_gradient_cl(self, E_field, step, bAlloc=True):
        """
        Compute gradient of scalar field on GPU using central differences.
        
        Args:
            E_field: (nx, ny, nz) numpy array, scalar energy field
            step: grid spacing
            bAlloc: allocate GPU buffers if True
        
        Returns:
            grads: (nx, ny, nz, 4) numpy array, (Fx, Fy, Fz, E) where F = -grad(E)
        """
        import numpy as np
        import pyopencl as cl

        nx, ny, nz = E_field.shape
        print(f"AFMulator.compute_gradient_cl: grid={nx}x{ny}x{nz}, step={step}")
        print(f"  Input E_field range: [{E_field.min():.6f}, {E_field.max():.6f}]")

        # Create input image (copy E_field to GPU)
        # Note: OpenCL images have z varying fastest, so we need shape (nz, ny, nx, 4)
        E_cl = np.zeros((nz, ny, nx, 4), dtype=np.float32)
        E_cl[:, :, :, 0] = E_field.transpose(2, 1, 0)  # (nx,ny,nz) -> (nz,ny,nx)
        
        mf = cl.mem_flags
        fmt = cl.ImageFormat(cl.channel_order.RGBA, cl.channel_type.FLOAT)
        
        if bAlloc or not hasattr(self, 'img_E_in') or getattr(self, 'img_E_in', None) is None:
            self.img_E_in = cl.Image(self.ctx, mf.READ_ONLY, fmt, shape=(nx, ny, nz))
        if bAlloc or not hasattr(self, 'img_F_out') or getattr(self, 'img_F_out', None) is None:
            self.img_F_out = cl.Image(self.ctx, mf.WRITE_ONLY, fmt, shape=(nx, ny, nz))

        # Upload input field
        origin = (0, 0, 0)
        region = (nx, ny, nz)
        cl.enqueue_copy(self.queue, self.img_E_in, E_cl, origin=origin, region=region)

        # Launch kernel
        gs = (self._roundup(nx, 8), self._roundup(ny, 8), self._roundup(nz, 4))
        ls = (8, 8, 4)
        self.prg.gradient_central_diff(self.queue, gs, ls, self.img_E_in, self.img_F_out, np.float32(step))

        # Read back results
        # Output is (Fx, Fy, Fz, E) in .xyzw
        F_cl = np.zeros((nz, ny, nx, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, F_cl, self.img_F_out, origin=origin, region=region)
        
        # Transpose back to (nx, ny, nz, 4)
        grads = F_cl.transpose(2, 1, 0, 3)
        
        print(f"AFMulator.compute_gradient_cl: done. grads shape={grads.shape}, range=[{grads.min():.4f},{grads.max():.4f}]")
        print(f"  Output Fx range: [{grads[...,0].min():.6f}, {grads[...,0].max():.6f}]")
        print(f"  Output Fy range: [{grads[...,1].min():.6f}, {grads[...,1].max():.6f}]")
        print(f"  Output Fz range: [{grads[...,2].min():.6f}, {grads[...,2].max():.6f}]")
        return grads

    def compute_gradient_fft_cl(self, E_field, step, bAlloc=True, bDebug=True, debug_dir='./fft_debug'):
        """
        Compute gradient of scalar field on GPU using FFT spectral differentiation.
        
        Uses gpyFFT for efficient spectral differentiation:
        ∂f/∂x = iFFT(ik_x * FFT(f))
        
        Args:
            E_field: (nx, ny, nz) numpy array, scalar energy field
            step: grid spacing
            bAlloc: allocate GPU buffers if True
            bDebug: enable debugging plots and prints
            debug_dir: directory to save debug plots
        
        Returns:
            grads: (nx, ny, nz, 4) numpy array, (Fx, Fy, Fz, E) where F = -grad(E)
        """
        import numpy as np
        import pyopencl as cl
        import pyopencl.array as cl_array
        from . import clUtils as clu
        import matplotlib.pyplot as plt
        import os
        
        nx, ny, nz = E_field.shape
        nxyz = nx * ny * nz
        print(f"AFMulator.compute_gradient_fft_cl: grid={nx}x{ny}x{nz}, step={step}")
        print(f"  Input E_field range: [{E_field.min():.6f}, {E_field.max():.6f}]")
        
        # Debug: Check input field for non-zero values
        idx_max = np.unravel_index(np.argmax(np.abs(E_field)), E_field.shape)
        print(f"  Input max abs at: {idx_max}, value: {E_field[idx_max]:.6f}")
        
        # Check grid size is FFT-friendly (clFFT requires 2,3,5,7 only)
        def check_fft_friendly(n):
            factors = []
            d = 2
            while d * d <= n:
                while n % d == 0:
                    factors.append(d)
                    n //= d
                d += 1
            if n > 1:
                factors.append(n)
            return all(p <= 7 for p in factors)
        
        for dim, name in [(nx, 'nx'), (ny, 'ny'), (nz, 'nz')]:
            if not check_fft_friendly(dim):
                print(f"WARNING: {name}={dim} not FFT-friendly (factors must be 2,3,5,7 only)")
        
        # Load gpyFFT
        clu.try_load_clFFT()
        FFT = clu.FFT
        
        # FFT requires complex buffers with reversed shape (nz, ny, nx)
        shape_fft = (nz, ny, nx)
        
        # Create or reuse buffers
        if bAlloc or not hasattr(self, 'fft_field_buf') or self.fft_field_buf is None:
            self.fft_field_buf = cl_array.Array(self.queue, shape_fft, dtype=np.complex64)
            self.fft_kx_buf = cl_array.Array(self.queue, shape_fft, dtype=np.complex64)
            self.fft_ky_buf = cl_array.Array(self.queue, shape_fft, dtype=np.complex64)
            self.fft_kz_buf = cl_array.Array(self.queue, shape_fft, dtype=np.complex64)
        
        # Upload real field to complex buffer (imaginary part = 0)
        # Note: transpose from (nx,ny,nz) to (nz,ny,nx) for FFT
        E_transposed = E_field.transpose(2, 1, 0).astype(np.float32)
        
        # Debug: Check transposed data
        print(f"  E_transposed shape: {E_transposed.shape}, range: [{E_transposed.min():.4f}, {E_transposed.max():.4f}]")
        
        # Create complex data on host and upload
        E_complex = np.zeros(shape_fft, dtype=np.complex64)
        E_complex.real = E_transposed
        self.fft_field_buf.set(E_complex)
        
        # Debug: Verify buffer contents
        buf_check = self.fft_field_buf.get()
        print(f"  Buffer after upload - shape: {buf_check.shape}, range: [{buf_check.real.min():.4f}, {buf_check.real.max():.4f}]")
        
        # Setup FFT plans
        fft_forward = FFT(self.ctx, self.queue, self.fft_field_buf, axes=(0, 1, 2))
        
        # Compute FFT
        event, = fft_forward.enqueue()
        event.wait()
        
        # Get FFT result for debugging
        fft_result = self.fft_field_buf.get()
        print(f"  FFT output shape: {fft_result.shape}, dtype: {fft_result.dtype}")
        print(f"  FFT output range: [{fft_result.real.min():.4e}, {fft_result.real.max():.4e}] + i[{fft_result.imag.min():.4e}, {fft_result.imag.max():.4e}]")
        
        # Find max abs value position in FFT
        fft_abs = np.abs(fft_result)
        idx_fft_max = np.unravel_index(np.argmax(fft_abs), fft_abs.shape)
        print(f"  FFT max abs at: {idx_fft_max}, value: {fft_abs[idx_fft_max]:.6e}")
        
        # Check if FFT has non-zero elements
        non_zero_fft = np.sum(fft_abs > 1e-10)
        print(f"  FFT non-zero elements (|val|>1e-10): {non_zero_fft} / {nxyz}")
        
        # Compute wave numbers for spectral differentiation
        # k = 2π * [0, 1, ..., N/2-1, -N/2, ..., -1] / (N * step)
        kx = 2.0 * np.pi * np.fft.fftfreq(nx, step)
        ky = 2.0 * np.pi * np.fft.fftfreq(ny, step)
        kz = 2.0 * np.pi * np.fft.fftfreq(nz, step)
        
        # Create 3D k-space grids
        # IMPORTANT: FFT output is (nz, ny, nx), so k-grids must match this ordering
        KZ, KY, KX = np.meshgrid(kz, ky, kx, indexing='ij')
        KX = KX.astype(np.complex64)
        KY = KY.astype(np.complex64)
        KZ = KZ.astype(np.complex64)
        
        print(f"  KX range: [{KX.min():.4f}, {KX.max():.4f}]")
        print(f"  KY range: [{KY.min():.4f}, {KY.max():.4f}]")
        print(f"  KZ range: [{KZ.min():.4f}, {KZ.max():.4f}]")
        
        # Multiply by ik in k-space (spectral differentiation)
        # ∂f/∂x → i * k_x * FFT(f) / nxyz (normalization for inverse FFT)
        # Note: Force = -gradient, so we multiply by -i*k
        fft_data = self.fft_field_buf.get()
        
        # Apply spectral derivative with proper normalization
        # The normalization 1/nxyz is needed because gpyFFT doesn't normalize
        kx_deriv = fft_data * (-1j * KX)  # / nxyz
        ky_deriv = fft_data * (-1j * KY)  # / nxyz
        kz_deriv = fft_data * (-1j * KZ)  # / nxyz
        
        print(f"  kx_deriv range: [{np.abs(kx_deriv).min():.4e}, {np.abs(kx_deriv).max():.4e}]")
        
        # Upload back to GPU
        self.fft_kx_buf.set(kx_deriv)
        self.fft_ky_buf.set(ky_deriv)
        self.fft_kz_buf.set(kz_deriv)
        
        # Create output buffers for inverse FFT
        grad_x_complex = cl_array.Array(self.queue, shape_fft, dtype=np.complex64)
        grad_y_complex = cl_array.Array(self.queue, shape_fft, dtype=np.complex64)
        grad_z_complex = cl_array.Array(self.queue, shape_fft, dtype=np.complex64)
        
        # Copy k-space data to output buffers
        grad_x_complex.set(self.fft_kx_buf.get())
        grad_y_complex.set(self.fft_ky_buf.get())
        grad_z_complex.set(self.fft_kz_buf.get())
        
        # Inverse FFTs
        ifft_x = FFT(self.ctx, self.queue, grad_x_complex, axes=(0, 1, 2))
        ifft_y = FFT(self.ctx, self.queue, grad_y_complex, axes=(0, 1, 2))
        ifft_z = FFT(self.ctx, self.queue, grad_z_complex, axes=(0, 1, 2))
        
        event, = ifft_x.enqueue()
        event.wait()
        event, = ifft_y.enqueue()
        event.wait()
        event, = ifft_z.enqueue()
        event.wait()
        
        # Get inverse FFT results
        grad_x_arr = grad_x_complex.get()
        grad_y_arr = grad_y_complex.get()
        grad_z_arr = grad_z_complex.get()
        
        print(f"  After IFFT - grad_x range: [{grad_x_arr.real.min():.4e}, {grad_x_arr.real.max():.4e}]")
        
        # Normalize by 1/nxyz (gpyFFT forward transform is unnormalized)
        grad_x_arr = grad_x_arr / nxyz
        grad_y_arr = grad_y_arr / nxyz
        grad_z_arr = grad_z_arr / nxyz
        
        print(f"  After norm - grad_x range: [{grad_x_arr.real.min():.4e}, {grad_x_arr.real.max():.4e}]")
        
        # Extract real parts and transpose back to (nx, ny, nz)
        grad_x = grad_x_arr.real.transpose(2, 1, 0)
        grad_y = grad_y_arr.real.transpose(2, 1, 0)
        grad_z = grad_z_arr.real.transpose(2, 1, 0)
        
        print(f"  Final grad_x range: [{grad_x.min():.4e}, {grad_x.max():.4e}]")
        
        # Debug plots
        if bDebug:
            os.makedirs(debug_dir, exist_ok=True)
            
            # Plot FFT magnitude (middle slice)
            mid_z = nz // 2
            mid_y = ny // 2
            mid_x = nx // 2
            
            fig, axes = plt.subplots(2, 3, figsize=(15, 10))
            
            # Input field
            im0 = axes[0, 0].imshow(E_field[:, :, mid_z], cmap='seismic', origin='lower')
            axes[0, 0].set_title(f'Input E_field (z={mid_z})')
            plt.colorbar(im0, ax=axes[0, 0])
            
            # FFT magnitude (log scale)
            fft_mag = np.abs(fft_result)
            im1 = axes[0, 1].imshow(np.log10(fft_mag[:, :, mid_x] + 1e-20), cmap='viridis', origin='lower')
            axes[0, 1].set_title(f'FFT magnitude log10 (x={mid_x})')
            plt.colorbar(im1, ax=axes[0, 1])
            
            # k-space derivative magnitude
            im2 = axes[0, 2].imshow(np.log10(np.abs(kx_deriv[:, :, mid_x]) + 1e-20), cmap='viridis', origin='lower')
            axes[0, 2].set_title(f'k-space d/dx log10 (x={mid_x})')
            plt.colorbar(im2, ax=axes[0, 2])
            
            # Output gradients
            im3 = axes[1, 0].imshow(grad_x[:, :, mid_z], cmap='seismic', origin='lower')
            axes[1, 0].set_title(f'grad_x (z={mid_z})')
            plt.colorbar(im3, ax=axes[1, 0])
            
            im4 = axes[1, 1].imshow(grad_y[:, :, mid_z], cmap='seismic', origin='lower')
            axes[1, 1].set_title(f'grad_y (z={mid_z})')
            plt.colorbar(im4, ax=axes[1, 1])
            
            im5 = axes[1, 2].imshow(grad_z[:, :, mid_z], cmap='seismic', origin='lower')
            axes[1, 2].set_title(f'grad_z (z={mid_z})')
            plt.colorbar(im5, ax=axes[1, 2])
            
            plt.tight_layout()
            plt.savefig(os.path.join(debug_dir, 'fft_gradient_debug.png'), dpi=150)
            plt.close()
            print(f"  Saved debug plot to {debug_dir}/fft_gradient_debug.png")
        
        # Stack into output format (Fx, Fy, Fz, E)
        grads = np.stack([grad_x, grad_y, grad_z, E_field], axis=-1)
        
        print(f"AFMulator.compute_gradient_fft_cl: done. grads shape={grads.shape}, range=[{grads.min():.4f},{grads.max():.4f}]")
        print(f"  Output Fx range: [{grads[...,0].min():.6f}, {grads[...,0].max():.6f}]")
        print(f"  Output Fy range: [{grads[...,1].min():.6f}, {grads[...,1].max():.6f}]")
        print(f"  Output Fz range: [{grads[...,2].min():.6f}, {grads[...,2].max():.6f}]")
        return grads

    @staticmethod
    def _roundup(n, loc):
        loc = max(loc, 1)
        return int(np.ceil(n/loc)*loc)


# ═══════════════════════════════════════════════════════════════════════════════
# Common utilities
# ═══════════════════════════════════════════════════════════════════════════════

ELEM_Z    = {'H':1, 'C':6, 'N':7, 'O':8}
Z_TO_ZVAL = {1:1, 6:4, 7:5, 8:6, 16:6}
RCUT_DEFAULT = {1:2.3, 6:2.6, 7:2.6, 8:2.5}   # Å (.wf rcutoff is Bohr; converted + margin)

def compute_df(Fz, dz):
    """df = -dFz/dz  (frequency-shift proxy).  Fz shape (nx,ny,nz)."""
    return -np.gradient(Fz, abs(dz), axis=2)

def fft_poisson(rho, step):
    """FFT Poisson solver: V(r) from charge density rho(r) on uniform grid with spacing `step`."""
    nx, ny, nz = rho.shape
    rho_k = np.fft.fftn(rho)
    kx = 2*np.pi * np.fft.fftfreq(nx, d=step)
    ky = 2*np.pi * np.fft.fftfreq(ny, d=step)
    kz = 2*np.pi * np.fft.fftfreq(nz, d=step)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    k2 = KX**2 + KY**2 + KZ**2;  k2[0,0,0] = 1.0
    V_k = 4.0*np.pi*COULOMB_CONST*rho_k / k2;  V_k[0,0,0] = 0.0
    return np.real(np.fft.ifftn(V_k)).astype(np.float32)

def build_gaussian_tip(grid_shape, step, sigma):
    """Build normalized Gaussian tip density kernel with FFT wrap-around. Returns (nx,ny,nz) float64."""
    nx, ny, nz = grid_shape
    rx = np.where(np.arange(nx) <= nx//2, np.arange(nx), np.arange(nx) - nx).astype(np.float64) * step
    ry = np.where(np.arange(ny) <= ny//2, np.arange(ny), np.arange(ny) - ny).astype(np.float64) * step
    rz = np.where(np.arange(nz) <= nz//2, np.arange(nz), np.arange(nz) - nz).astype(np.float64) * step
    RX, RY, RZ = np.meshgrid(rx, ry, rz, indexing='ij')
    R2 = RX**2 + RY**2 + RZ**2
    kernel = (1.0 / (2*np.pi*sigma**2)**1.5) * np.exp(-R2 / (2*sigma**2))
    dV = step**3
    tip_int = kernel.sum() * dV
    print(f"  Tip kernel: integral={tip_int:.6f} (should be ~1.0), max={kernel.max():.4f} e/Å³")
    return kernel

def _interp3(g, coords):
    """Trilinear interpolation of 3D grid g at fractional coords (N,3). Returns (N,) float32."""
    from scipy.ndimage import map_coordinates
    return _mc(g, coords.T, order=1, mode='nearest').reshape(coords.shape[:-1]).astype(np.float32)

def get_tip_kernel(rho_t):
    """
    Prepare tip density kernel for convolution.
    Reverse the tip density (for correlation) and roll to center it at index (0,0,0).
    """
    nx, ny, nz = rho_t.shape
    rho_t_rev = rho_t[::-1, ::-1, ::-1]
    return np.roll(np.roll(np.roll(rho_t_rev, -(nx//2), axis=0), -(ny//2), axis=1), -(nz//2), axis=2)

def get_pauli_convolution(rho_s, rho_t, dV, A_pauli=16.0, beta_pauli=1.0):
    """
    Compute Pauli repulsion field using FFT convolution.
    E_pauli = A_pauli * dV * (rho_s**beta * rho_t**beta)
    """
    if beta_pauli != 1.0:
        rho_s = rho_s**beta_pauli
        rho_t = rho_t**beta_pauli
    
    rho_s_k = np.fft.fftn(rho_s.astype(np.float64))
    kernel  = get_tip_kernel(rho_t)
    rho_t_k = np.fft.fftn(kernel.astype(np.float64))
    
    conv = np.real(np.fft.ifftn(rho_s_k * rho_t_k))
    return (A_pauli * dV * conv).astype(np.float32)


# ═══════════════════════════════════════════════════════════════════════════════
# Fireball SCF + density projection
# ═══════════════════════════════════════════════════════════════════════════════

def run_fireball_scf(xyz_path, fdata_dir, nscf=200, verbosity=0):
    """
    Run Fireball SCF on molecule from xyz_path.
    Returns dict with keys: atomTypes, atomPos, rho_sparse, neighs, q_mulliken, natoms.
    NOTE: changes cwd to tests/pyFireball temporarily (Fireball requirement).
    """
    import pyBall.FireCore as fc
    _root = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))
    _fball_cwd = os.path.join(_root, 'tests', 'pyFireball')
    orig_cwd = os.getcwd()
    os.chdir(_fball_cwd)
    # Ensure Fdata symlink
    fdata_local = os.path.join(_fball_cwd, 'Fdata')
    if not os.path.exists(os.path.join(fdata_local, 'info.dat')):
        if os.path.lexists(fdata_local):
            os.unlink(fdata_local)
        os.symlink(fdata_dir, fdata_local)
        print(f"Created Fdata symlink → {fdata_dir}")
    # Parse xyz
    with open(xyz_path) as f:
        lines = f.readlines()
    natoms = int(lines[0])
    atomTypes = []; atomPos = []
    for line in lines[2:2+natoms]:
        p = line.split()
        atomTypes.append(ELEM_Z[p[0]])
        atomPos.append([float(p[1]), float(p[2]), float(p[3])])
    atomTypes = np.array(atomTypes, dtype=np.int32)
    atomPos   = np.array(atomPos,   dtype=np.float64)
    print(f"Fireball SCF: {natoms} atoms  z=[{atomPos[:,2].min():.2f},{atomPos[:,2].max():.2f}]")
    fc.setVerbosity(verbosity)
    fc.preinit()
    fc.init(atomTypes, atomPos)
    fc.SCF(atomPos, nmax_scf=nscf)
    dims   = fc.get_HS_dims()
    neighs = fc.get_HS_neighs(dims)
    neighs = fc.get_rho_sparse(dims, data=neighs)
    rho_sparse = neighs.rho
    print(f"rho_sparse shape={rho_sparse.shape}  |max|={np.abs(rho_sparse).max():.4f}")
    # Mulliken charges
    charges2d = np.zeros((natoms, 1), dtype=np.float64)
    fc.getCharges(charges2d)
    q_pop = charges2d[:,0]
    z_val = np.array([Z_TO_ZVAL.get(int(z), int(z)) for z in atomTypes], dtype=float)
    q_mulliken = (q_pop - z_val).astype(np.float32)
    print(f"Mulliken charges: sum={q_mulliken.sum():.4f}  range=[{q_mulliken.min():.3f},{q_mulliken.max():.3f}]")
    os.chdir(orig_cwd)
    return dict(atomTypes=atomTypes, atomPos=atomPos, rho_sparse=rho_sparse,
                neighs=neighs, q_mulliken=q_mulliken, natoms=natoms)


def setup_density_grid(atomPos, step=0.1, margin=4.0, z_extra=6.0, block=8):
    """
    Setup 3D grid spec for density projection.
    Grid dimensions are rounded up to multiples of block for GPU performance.
    Margin is increased to accommodate rounding, ensuring at least specified margin around molecule.
    Returns (grid_spec, origin, ngrid).
    """
    pos_min_raw = atomPos.min(axis=0) - margin
    pos_max_raw = atomPos.max(axis=0) + np.array([margin, margin, margin + z_extra])
    span_raw = pos_max_raw - pos_min_raw
    ngrid_raw = np.ceil(span_raw / step)
    ngrid = (np.ceil(ngrid_raw / block).astype(int) * block)
    total_span = ngrid * step
    # Center grid on molecule, extending equally on both sides to accommodate rounding
    origin = (0.5*(pos_min_raw + pos_max_raw) - 0.5*total_span).astype(np.float32)
    grid_spec = {
        'origin': origin,
        'dA': [step, 0., 0.], 'dB': [0., step, 0.], 'dC': [0., 0., step],
        'ngrid': ngrid.astype(int),
    }
    return grid_spec, origin, ngrid


def _onsite_occ(Z):
    """Neutral-atom valence occupations for (s, py, pz, px) kernel convention."""
    if Z == 1:  return np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float32)
    if Z == 6:  return np.array([2.0, 2.0/3, 2.0/3, 2.0/3], dtype=np.float32)
    if Z == 7:  return np.array([2.0, 1.0, 1.0, 1.0], dtype=np.float32)
    if Z == 8:  return np.array([2.0, 4.0/3, 4.0/3, 4.0/3], dtype=np.float32)
    return np.array([float(Z_TO_ZVAL.get(int(Z), int(Z))), 0.0, 0.0, 0.0], dtype=np.float32)


def build_neutral_atom_rho(atomTypes, neighs, natoms):
    """Build neutral-atom (promolecule) on-site density matrix matching rho_sparse shape."""
    rho_na = np.zeros_like(neighs.rho, dtype=np.float32)
    neigh_j = neighs.neigh_j.reshape(natoms, -1)
    for i in range(natoms):
        slots = np.where(neigh_j[i] == (i+1))[0]
        if len(slots) == 0:
            raise RuntimeError(f"No self-neighbor slot for atom i={i}")
        iself = int(slots[0])
        occ = _onsite_occ(int(atomTypes[i]))
        rho_na[i, iself, :, :] = 0.0
        for k in range(4):
            rho_na[i, iself, k, k] = occ[k]
    return rho_na


def project_density_grids(rho_sparse, rho_na_sparse, neighs, atomTypes, atomPos, grid_spec,
                          fdata_basis_dir, step=0.15, proj_kwargs=None, verbosity=0):
    """
    Project SCF and neutral-atom densities onto grid.
    Returns (rho_grid, rho_na_grid, rho_diff, projector).
    proj_kwargs: dict passed to GridProjector and .project() (e.g. use_tiled, debug flags).
    """
    from pyBall.FireballOCL import Grid as ocl_grid
    if proj_kwargs is None:
        proj_kwargs = {}
    # Extract GridProjector constructor kwargs
    gp_kw = {k: proj_kwargs[k] for k in
             ['debug_early_exit','debug_clear_only','debug_return0',
              'debug_read_task','debug_read_grid'] if k in proj_kwargs}
    projector = ocl_grid.GridProjector(fdata_dir=fdata_basis_dir, verbosity=verbosity, **gp_kw)
    projector.load_basis(sorted(set(atomTypes.tolist())))
    atoms_dict = {
        'pos':  atomPos,
        'Rcut': np.array([RCUT_DEFAULT.get(int(z), 4.5) for z in atomTypes]),
        'type': atomTypes,
    }
    use_tiled = proj_kwargs.get('use_tiled', True)
    nMaxAtom  = proj_kwargs.get('nMaxAtom', 64)
    tasks     = proj_kwargs.get('tasks', None)
    print("Projecting SCF density to grid...")
    rho_grid = projector.project(rho_sparse, neighs, atoms_dict, grid_spec,
                                  tasks=tasks, nMaxAtom=nMaxAtom, use_tiled=use_tiled)
    dV = step**3
    print(f"rho_grid shape={rho_grid.shape} range=[{rho_grid.min():.5f},{rho_grid.max():.5f}]")
    print(f"  Integrated electrons = {rho_grid.sum()*dV:.2f}")
    print("Projecting neutral-atom density...")
    rho_na_grid = projector.project(rho_na_sparse, neighs, atoms_dict, grid_spec,
                                     tasks=tasks, nMaxAtom=nMaxAtom, use_tiled=use_tiled)
    rho_diff = (rho_grid - rho_na_grid).astype(np.float32)
    print(f"rho_na_grid: range=[{rho_na_grid.min():.5f},{rho_na_grid.max():.5f}] integral={rho_na_grid.sum()*dV:.2f} e")
    print(f"rho_diff: range=[{rho_diff.min():.5f},{rho_diff.max():.5f}] integral={rho_diff.sum()*dV:.2f} e")
    return rho_grid, rho_na_grid, rho_diff, projector


# ═══════════════════════════════════════════════════════════════════════════════
# FDBM force field computation
# ═══════════════════════════════════════════════════════════════════════════════

def compute_fdbm_forcefield(rho_grid, V_ES, origin, step, atomPos,
                             scan_xs, scan_ys, probe_heights, mol_z=0.0,
                             A_pauli=16.0, C6_CO=30.0, q_CO=-0.05, sigma_tip=0.7,
                             force_model='fdbm', es_model='poisson', q_mulliken=None,
                             rho_smooth=None):
    """
    Compute FDBM/gradient forces on scan grid.
    Returns dict with keys: FEs_raw, Fz_raw, F_pauli, F_vdw, F_es,
        grads_E_Pauli, grads_E_ES (or grads_rho, grads_VES for gradient model),
        E_Pauli_field, E_ES_field (only for fdbm model).
    """
    from scipy.ndimage import gaussian_filter
    nx_d, ny_d, nz_d = rho_grid.shape
    dV = step**3
    nx_s = len(scan_xs);  ny_s = len(scan_ys);  nz = len(probe_heights)
    natoms = len(atomPos)
    # Scan positions
    XX, YY, ZZ = np.meshgrid(scan_xs, scan_ys, probe_heights + mol_z, indexing='ij')
    flat_pos = np.stack([XX.ravel(), YY.ravel(), ZZ.ravel()], axis=1)
    grid_c = (flat_pos - origin) / step
    # Smooth density (for gradient model fallback)
    if rho_smooth is None:
        sigma_elec = 0.8 / step
        rho_smooth = gaussian_filter(rho_grid.astype(np.float64), sigma=sigma_elec).astype(np.float32)
    result = {}
    if force_model == 'fdbm':
        print(f"  Tip model: Gaussian sigma={sigma_tip:.2f} Å")
        rho_tip_kernel = build_gaussian_tip((nx_d, ny_d, nz_d), step, sigma_tip)
        print("  FFT convolution: Pauli (density overlap) and ES (tip in Hartree)...")
        fft_rho_s   = np.fft.rfftn(rho_grid.astype(np.float64))
        fft_rho_tip = np.fft.rfftn(rho_tip_kernel)
        fft_VES     = np.fft.rfftn(V_ES.astype(np.float64))
        E_Pauli_field = (A_pauli * dV * np.fft.irfftn(fft_rho_s * fft_rho_tip, s=(nx_d, ny_d, nz_d))).astype(np.float32)
        E_ES_field    = (q_CO    * dV * np.fft.irfftn(fft_VES   * fft_rho_tip, s=(nx_d, ny_d, nz_d))).astype(np.float32)
        print(f"  E_Pauli field: range=[{E_Pauli_field.min():.4f},{E_Pauli_field.max():.4f}] eV")
        print(f"  E_ES field:    range=[{E_ES_field.min():.4f},{E_ES_field.max():.4f}] eV")
        print("  Computing total energy field and gradient (optimized: one gradient computation)...")
        # Compute total energy and gradient once
        E_total = E_Pauli_field + E_ES_field
        print(f"  E_total range: [{E_total.min():.4f},{E_total.max():.4f}] eV")
        # GPU gradient computation: compute_gradient_cl returns (Fx,Fy,Fz,E) where F = -grad(E)
        # So we need to negate to get gradients: grad = -F
        grads_cl_total = self.compute_gradient_cl(E_total, step, bAlloc=True)
        # Extract gradients (negate forces): grad = -F
        grads_total = [-grads_cl_total[..., a].astype(np.float32) for a in range(3)]
        result['E_Pauli_field'] = E_Pauli_field
        result['E_ES_field']    = E_ES_field
        result['grads_E_Pauli'] = grads_total  # Total gradient contains both Pauli and ES
        result['grads_E_ES']    = grads_total  # Same total gradient (for backward compat)
        result['grads_total']   = grads_total  # New: explicit total gradient
    else:
        print("  Computing gradient-based Pauli force (old model) using GPU (optimized: one gradient computation)...")
        # Compute total potential: A_pauli * rho_smooth + q_CO * V_ES
        E_total = A_pauli * rho_smooth + q_CO * V_ES
        print(f"  E_total range: [{E_total.min():.4f},{E_total.max():.4f}] eV")
        # Compute gradient of total potential once
        grads_cl_total = self.compute_gradient_cl(E_total, step, bAlloc=True)
        # Extract gradients (negate forces): grad = -F
        grads_total = [-grads_cl_total[..., a].astype(np.float32) for a in range(3)]
        result['grads_rho'] = grads_total  # Total gradient (for backward compat)
        result['grads_VES'] = grads_total  # Total gradient (for backward compat)
        result['grads_total'] = grads_total  # New: explicit total gradient
    # Interpolate forces at scan positions
    print("  Interpolating forces at scan positions...")
    F_pauli = np.zeros((flat_pos.shape[0], 3), dtype=np.float32)
    F_vdw   = np.zeros((flat_pos.shape[0], 3), dtype=np.float32)
    F_es    = np.zeros((flat_pos.shape[0], 3), dtype=np.float32)
    if force_model == 'fdbm':
        # Use total gradient for both Pauli and ES (gradient of E_total = gradient of E_Pauli + E_ES)
        for a in range(3):
            F_pauli[:,a] = -_interp3(grads_total[a], grid_c)
        # ES is included in total gradient, set to zero for accounting
        # In practice, F_pauli now contains both Pauli and ES forces
    else:
        # Use total gradient for gradient model
        for a in range(3):
            F_pauli[:,a] = -_interp3(grads_total[a], grid_c)
        if es_model != 'poisson':
            # For non-poisson ES model, still use direct atom loop (ES not in gradient)
            RA2_ES = 2.0**2
            for ia in range(natoms):
                dr = flat_pos - atomPos[ia]; r2 = np.sum(dr**2, axis=1); r3_es = (r2 + RA2_ES)**1.5
                for a in range(3):
                    F_es[:,a] += COULOMB_CONST * q_mulliken[ia] * q_CO * dr[:,a] / r3_es
    # London C6/r^6 atom sum
    print("  Atom loop: London vdW...")
    RA2_VDW = 1.5**2
    for ia in range(natoms):
        dr = flat_pos - atomPos[ia]; r2 = np.sum(dr**2, axis=1); r8 = (r2 + RA2_VDW)**4
        for a in range(3):
            F_vdw[:,a] -= C6_CO * 6.0 * dr[:,a] / r8
    F_fdbm = np.zeros((flat_pos.shape[0], 4), dtype=np.float32)
    F_fdbm[:,:3] = F_pauli + F_vdw + F_es
    FEs_raw = F_fdbm.reshape(nx_s, ny_s, nz, 4)
    Fz_raw  = FEs_raw[:,:,:,2]
    print(f"Fz_raw: min={Fz_raw.min():.4f}  max={Fz_raw.max():.4f}  mean={Fz_raw.mean():.4f} eV/Å")
    result.update(FEs_raw=FEs_raw, Fz_raw=Fz_raw,
                  F_pauli=F_pauli.reshape(nx_s,ny_s,nz,3),
                  F_vdw=F_vdw.reshape(nx_s,ny_s,nz,3),
                  F_es=F_es.reshape(nx_s,ny_s,nz,3),
                  rho_smooth=rho_smooth)
    return result


def make_fdbm_force_func(origin, step, grads_list_pauli, grads_list_es, atomPos,
                          C6_CO=30.0, force_model='fdbm',
                          A_pauli=16.0, q_CO=-0.05,
                          grads_rho=None, grads_VES=None, es_model='poisson', q_mulliken=None):
    """
    Return a callable force_func(positions) → (N,3) float32
    for use in pp_relax_2d.
    """
    RA2_VDW = 1.5**2
    natoms = len(atomPos)
    def force_func(positions):
        gc = (positions - origin) / step
        F = np.zeros((len(positions), 3), dtype=np.float64)
        if force_model == 'fdbm':
            for a in range(3):
                F[:,a] += -_interp3(grads_list_pauli[a], gc).astype(np.float64)
                F[:,a] += -_interp3(grads_list_es[a],    gc).astype(np.float64)
        else:
            for a in range(3):
                F[:,a] += -A_pauli * _interp3(grads_rho[a], gc).astype(np.float64)
            if es_model == 'poisson':
                for a in range(3):
                    F[:,a] += -q_CO * _interp3(grads_VES[a], gc).astype(np.float64)
        for ia in range(natoms):
            dr = positions - atomPos[ia]; r2 = np.sum(dr**2, axis=1); r8 = (r2 + RA2_VDW)**4
            for a in range(3):
                F[:,a] -= C6_CO * 6.0 * dr[:,a] / r8
        return F.astype(np.float32)
    return force_func


def pp_relax_2d(force_func, scan_xs, scan_ys, probe_heights, mol_z=0.0,
                K_LAT=0.03, N_RELAX=30, PP_MARGIN=2.0, step=0.15):
    """
    2D lateral PP relaxation per height slice.
    force_func(positions_Nx3) → (N,3) forces.
    Returns (FEs_relax, tip_disp) where:
        FEs_relax: (nx_s, ny_s, nz, 4) forces at relaxed positions
        tip_disp: dict with 'dx' (nx_s, ny_s, nz) and 'dy' (nx_s, ny_s, nz) displacement arrays
    """
    from scipy.ndimage import map_coordinates as _mc
    nx_s = len(scan_xs); ny_s = len(scan_ys); nz = len(probe_heights)
    XX2, YY2 = np.meshgrid(scan_xs, scan_ys, indexing='ij')
    # Force-map grid with lateral margin
    pp_x0 = scan_xs[0] - PP_MARGIN;  pp_x1 = scan_xs[-1] + PP_MARGIN
    pp_y0 = scan_ys[0] - PP_MARGIN;  pp_y1 = scan_ys[-1] + PP_MARGIN
    pp_nx = int(np.ceil((pp_x1 - pp_x0) / step)) + 1
    pp_ny = int(np.ceil((pp_y1 - pp_y0) / step)) + 1
    pp_xs = pp_x0 + np.arange(pp_nx) * step
    pp_ys = pp_y0 + np.arange(pp_ny) * step
    PP_X, PP_Y = np.meshgrid(pp_xs, pp_ys, indexing='ij')
    def _interp2d(F2d, px, py):
        cx = np.clip((px - pp_x0) / step, 0, F2d.shape[0]-1.001)
        cy = np.clip((py - pp_y0) / step, 0, F2d.shape[1]-1.001)
        return _mc(F2d, [cx.ravel(), cy.ravel()], order=1, mode='nearest').reshape(px.shape).astype(np.float32)
    FEs_relax = np.zeros((nx_s, ny_s, nz, 4), dtype=np.float32)
    tip_disp = {'dx': np.zeros((nx_s, ny_s, nz), dtype=np.float32),
                'dy': np.zeros((nx_s, ny_s, nz), dtype=np.float32)}
    for iz in range(nz):
        probe_z = probe_heights[iz] + mol_z
        PP_Z = np.full_like(PP_X, probe_z)
        pp_flat = np.stack([PP_X.ravel(), PP_Y.ravel(), PP_Z.ravel()], axis=1)
        debug_print(1, f"  iz={iz:2d} h={probe_heights[iz]:.1f} Å  probe_z={probe_z:.2f} Å  building {pp_nx}x{pp_ny} force map...")
        FF = force_func(pp_flat)
        FF_x = FF[:,0].reshape(pp_nx, pp_ny)
        FF_y = FF[:,1].reshape(pp_nx, pp_ny)
        FF_z = FF[:,2].reshape(pp_nx, pp_ny)
        probe_x = XX2.astype(np.float32).copy()
        probe_y = YY2.astype(np.float32).copy()
        vx = np.zeros_like(probe_x); vy = np.zeros_like(probe_y)
        for _ in range(N_RELAX):
            Fx_s = _interp2d(FF_x, probe_x, probe_y) - K_LAT*(probe_x - XX2)
            Fy_s = _interp2d(FF_y, probe_x, probe_y) - K_LAT*(probe_y - YY2)
            vx = 0.8*vx + 0.3*Fx_s;  probe_x += vx*0.3
            vy = 0.8*vy + 0.3*Fy_s;  probe_y += vy*0.3
        # Store displacement (final - initial position)
        tip_disp['dx'][:,:,iz] = probe_x - XX2
        tip_disp['dy'][:,:,iz] = probe_y - YY2
        FEs_relax[:,:,iz,0] = _interp2d(FF_x, probe_x, probe_y)
        FEs_relax[:,:,iz,1] = _interp2d(FF_y, probe_x, probe_y)
        FEs_relax[:,:,iz,2] = _interp2d(FF_z, probe_x, probe_y)
    Fz_relax = FEs_relax[:,:,:,2]
    debug_print(1, f"Fz_relax: min={Fz_relax.min():.4f}  max={Fz_relax.max():.4f}  mean={Fz_relax.mean():.4f} eV/Å")
    return FEs_relax, tip_disp


def pp_relax_2d_cl(afmulator, F_total, origin, step,
                   scan_xs, scan_ys, probe_heights, mol_z=0.0,
                   K_LAT=0.0031, dpos0=None, stiffness=None, relax_pars=None,
                   ppm_mode=True):
    """
    Backward-compat wrapper: calls AFMulator.setup_fdbm_grid + scan_fdbm.
    Prefer calling those methods directly for repeated scans (avoids re-uploading image).

    ppm_mode=True:  physical PPM with radial CO bond (L=3 Ang, K_RAD=0.1248 eV/Ang^2).
    ppm_mode=False: simplified harmonic, probe pinned at scan heights.
    """
    afmulator.setup_fdbm_grid(F_total, origin, step)
    FEs_relax = afmulator.scan_fdbm(
        scan_xs, scan_ys, probe_heights, mol_z=mol_z,
        ppm_mode=ppm_mode, K_LAT=K_LAT,
        stiffness=stiffness, dpos0=dpos0, relax_pars=relax_pars
    )
    nx_s, ny_s, nz_s = FEs_relax.shape[:3]
    tip_disp = {'dx': np.zeros((nx_s, ny_s, nz_s), dtype=np.float32),
                'dy': np.zeros((nx_s, ny_s, nz_s), dtype=np.float32)}
    return FEs_relax, tip_disp


# ═══════════════════════════════════════════════════════════════════════════════
# FDBM AFM field computation helpers
# ═══════════════════════════════════════════════════════════════════════════════

# Pauli parameters fitted against DFTB z-scan with raw overlap (A=1,beta=1 convolution)
# E_pauli = A_pauli * overlap_raw^beta_pauli
PAULI_FITTED_DEFAULTS = {
    'mio-1-1': {'A': 787.22, 'beta': 1.2371},  # fitted for pentacene atom 0, raw overlap (A=1,beta=1 convolution)
    '3ob-3-1': {'A': 509.28, 'beta': 1.0586},  # fitted for pentacene atom 0, raw overlap (A=1,beta=1 convolution)
}

def compute_pauli_overlap(rho_grid, rho_tip_total, step, tip_rolled=False):
    """Compute raw Pauli overlap via FFT cross-correlation (A=1, beta=1).

    overlap(R) = dV * IFFT(FFT(rho_grid) * conj(FFT(rho_tip)))
    This is the pure density-density overlap integral at each tip position R.
    Clipped to [1e-30, inf] to allow safe power-law scaling.

    If tip_rolled=True the tip already has O at index (0,0,0).
    Returns overlap_raw (nx,ny,nz) float32.
    """
    dV = step**3
    overlap_raw = dV * np.real(np.fft.ifftn(np.fft.fftn(rho_grid) * np.conj(np.fft.fftn(rho_tip_total)))).astype(np.float32)
    return np.clip(overlap_raw, 1e-30, None)


def scale_pauli_field(overlap_raw, step, A_pauli, beta_pauli, return_grads=True):
    """Scale raw overlap into energy field: E_pauli = A_pauli * overlap^beta_pauli.

    Args:
        overlap_raw: (nx, ny, nz) raw overlap field
        step: grid spacing
        A_pauli: Pauli amplitude parameter
        beta_pauli: Pauli exponent parameter
        return_grads: if True, compute and return gradients (default True for backward compat)

    Returns:
        E_pauli: (nx, ny, nz) Pauli energy field
        grads: (nx, ny, nz, 3) Pauli gradients (if return_grads=True)
    """
    E_pauli = A_pauli * (overlap_raw ** beta_pauli)
    if return_grads:
        grads = np.stack([np.gradient(E_pauli, step, axis=i) for i in range(3)], axis=-1)
        return E_pauli, grads
    return E_pauli


def compute_pauli_field(rho_grid, rho_tip_total, step, A_pauli=1.0, beta_pauli=1.0, tip_rolled=False):
    """Compute Pauli repulsion field via FFT cross-correlation.

    Convenience wrapper: compute_pauli_overlap then scale_pauli_field.
    Default A=1, beta=1 returns raw overlap as energy field.
    Returns (E_pauli_field, grads_E_pauli) where grads shape is (nx,ny,nz,3).
    """
    overlap_raw = compute_pauli_overlap(rho_grid, rho_tip_total, step, tip_rolled=tip_rolled)
    return scale_pauli_field(overlap_raw, step, A_pauli, beta_pauli)

def compute_es_conv_field(V_ES, rho_tip_delta, step, tip_rolled=False, return_grads=True):
    """Convolve electrostatic potential with tip delta-density.

    Args:
        V_ES: (nx, ny, nz) electrostatic potential
        rho_tip_delta: (nx, ny, nz) tip delta density
        step: grid spacing
        tip_rolled: if True, use rolled tip (default False)
        return_grads: if True, compute and return gradients (default True for backward compat)

    Returns:
        E_es: (nx, ny, nz) electrostatic energy field
        grads: (nx, ny, nz, 3) electrostatic gradients (if return_grads=True)
    """
    dV = step**3
    nx_t, ny_t, nz_t = rho_tip_delta.shape
    flipped = rho_tip_delta[::-1,::-1,::-1]
    if tip_rolled:
        tip_kernel = flipped
    else:
        tip_kernel = np.roll(np.roll(np.roll(flipped, -(nx_t//2), axis=0), -(ny_t//2), axis=1), -(nz_t//2), axis=2)
    E_es = dV * np.real(np.fft.ifftn(np.fft.fftn(V_ES) * np.fft.fftn(tip_kernel))).astype(np.float32)
    if return_grads:
        grads = np.stack([np.gradient(E_es, step, axis=i) for i in range(3)], axis=-1)
        return E_es, grads
    return E_es

def compute_vdw_field(atomPos, atomTypes, origin, step, ngrid, C6_table=None, C6_CO=30.0, RA=1.5):
    """Compute vdW dispersion field C6/r^6 on grid.

    Args:
        atomPos: (natoms, 3) positions in Angstrom
        atomTypes: list/array of atomic Z numbers or element symbols
        origin: (3,) grid origin
        step: grid spacing
        ngrid: (3,) grid dimensions
        C6_table: dict mapping Z -> C6 value. Default: {1:6.5, 6:24.0, 7:20.0, 8:15.0}
        C6_CO: C6 value for CO tip
        RA: damping radius in Angstrom

    Returns (E_vdw, grads_vdw).
    """
    if C6_table is None:
        C6_table = {1: 6.5, 6: 24.0, 7: 20.0, 8: 15.0}
    nx, ny, nz = [int(i) for i in ngrid[:3]]
    xs = origin[0] + np.arange(nx)*step
    ys = origin[1] + np.arange(ny)*step
    zs = origin[2] + np.arange(nz)*step
    XX, YY, ZZ = np.meshgrid(xs, ys, zs, indexing='ij')
    E_vdw = np.zeros((nx, ny, nz), dtype=np.float32)
    RA2 = RA**2
    for ia in range(len(atomPos)):
        Z = C6_table.get(atomTypes[ia], 10.0)
        E_vdw -= np.sqrt(Z * C6_CO) / ((XX-atomPos[ia,0])**2 + (YY-atomPos[ia,1])**2 + (ZZ-atomPos[ia,2])**2 + RA2)**3
    grads = np.stack([np.gradient(E_vdw, step, axis=i) for i in range(3)], axis=-1)
    return E_vdw, grads

def read_grid_spec_from_log(log_path):
    """Read grid origin and ngrid from FDBM step1_density/log.txt.

    Args:
        log_path: path to log.txt file

    Returns:
        (origin, ngrid, step) as numpy arrays, or (None, None, None) if not found
    """
    if not os.path.exists(log_path):
        return None, None, None
    with open(log_path, 'r') as f:
        for line in f:
            if 'Grid:' in line and 'origin=' in line and 'ngrid=' in line:
                parts = line.split('ngrid=')
                origin_str = parts[0].split('origin=')[1].strip()
                ngrid_str = parts[1].strip()
                origin = np.array([float(x) for x in origin_str.strip('[]').split()])
                ngrid = np.array([int(x) for x in ngrid_str.strip('[]').split()])
                # step is usually after ngrid, e.g. "ngrid=[152 88 96] step=0.15"
                step = 0.15
                if 'step=' in line:
                    step_parts = line.split('step=')
                    if len(step_parts) > 1:
                        try:
                            step = float(step_parts[1].split()[0])
                        except ValueError:
                            pass
                return origin, ngrid, step
    return None, None, None

# ═══════════════════════════════════════════════════════════════════════════════
# Single-atom projection test helper
# ═══════════════════════════════════════════════════════════════════════════════

def project_single_atom(Z, rho_4x4, step, margin, fdata_basis_dir, use_tiled=True, verbosity=0):
    """
    Project density for a single atom at origin with given 4x4 density matrix block.
    Returns (rho_grid, grid_spec, integral, projector).
    """
    from pyBall.FireballOCL import Grid as ocl_grid
    block = 8
    n = int(np.ceil(2*margin / step / block)) * block
    origin = np.array([-n*step/2, -n*step/2, -n*step/2], dtype=np.float32)
    grid_spec = {
        'origin': origin,
        'dA': [step, 0., 0.], 'dB': [0., step, 0.], 'dC': [0., 0., step],
        'ngrid': np.array([n, n, n], dtype=np.int32),
    }
    atomPos   = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
    atomTypes = np.array([Z], dtype=np.int32)
    atoms_dict = {
        'pos': atomPos,
        'Rcut': np.array([RCUT_DEFAULT.get(Z, 4.5)], dtype=np.float64),
        'type': atomTypes,
    }
    neigh_max  = 1; numorb_max = 4
    rho_sparse = np.zeros((1, neigh_max, numorb_max, numorb_max), dtype=np.float32)
    rho_sparse[0, 0, :, :] = rho_4x4
    neigh_j = np.zeros((1, neigh_max), dtype=np.int32)
    neigh_j[0, 0] = 1
    class _FakeNeighs:
        pass
    neighs = _FakeNeighs()
    neighs.rho = rho_sparse; neighs.neigh_j = neigh_j.ravel()
    neighs.neigh_max = neigh_max; neighs.numorb_max = numorb_max
    projector = ocl_grid.GridProjector(fdata_dir=fdata_basis_dir, verbosity=verbosity)
    projector.load_basis([Z])
    print(f"  [project_single_atom] Z={Z} rho_diag={np.diag(rho_4x4)} grid={n}^3 step={step}")
    rho_grid = projector.project(rho_sparse, neighs, atoms_dict, grid_spec, nMaxAtom=64, use_tiled=use_tiled)
    dV = step**3
    integral = float(rho_grid.sum() * dV)
    print(f"  rho_grid shape={rho_grid.shape}  range=[{rho_grid.min():.6f}, {rho_grid.max():.6f}]")
    print(f"  Integrated electrons = {integral:.6f}")
    return rho_grid, grid_spec, integral, projector




# Singleton AFMulator instance for dispersion computation
_dispersion_afmulator = None

def compute_dispersion_grid(atomPos, atomTypes, origin, step, ngrid, C6_atom_dict=None, C6_CO=30.0, RA=1.5, use_opencl=True, return_grads=True):
    """
    Compute C6/r^6 dispersion energy grid for CO tip.

    Args:
        atomPos: (natoms, 3) atom positions in Angstrom
        atomTypes: (natoms,) atomic numbers
        origin: (3,) grid origin
        step: grid spacing in Angstrom
        ngrid: (3,) grid dimensions
        C6_atom_dict: dict mapping Z to C6 coefficients (default: {1:6.5, 6:24.0, 7:20.0, 8:15.0})
        C6_CO: C6 coefficient for CO tip
        RA: damping radius in Angstrom
        use_opencl: if True, use OpenCL GPU acceleration (default True)
        return_grads: if True, compute and return gradients (default True for backward compat)

    Returns:
        E_vdw: (nx, ny, nz) dispersion energy field
        grads: (nx, ny, nz, 3) dispersion gradients (if return_grads=True)
    """
    if C6_atom_dict is None:
        C6_atom_dict = {1: 6.5, 6: 24.0, 7: 20.0, 8: 15.0}

    # Use OpenCL version if requested
    if use_opencl:
        global _dispersion_afmulator
        # Create singleton AFMulator instance on first call
        if _dispersion_afmulator is None:
            _dispersion_afmulator = AFMulator(use_morse=False, nloc=32)

        # Allocate buffers on first call (check if attribute exists)
        if not hasattr(_dispersion_afmulator, 'disp_atoms_cl') or _dispersion_afmulator.disp_atoms_cl is None:
            _dispersion_afmulator.realloc_dispersion_buffers(len(atomPos))

        if return_grads:
            return _dispersion_afmulator.compute_dispersion_grid_cl(
                atomPos, atomTypes, origin, step, ngrid,
                C6_atom_dict=C6_atom_dict, C6_CO=C6_CO, RA=RA, bAlloc=False, return_grads=True
            )
        else:
            # Compute energy only (no gradient)
            return _dispersion_afmulator.compute_dispersion_grid_cl(
                atomPos, atomTypes, origin, step, ngrid,
                C6_atom_dict=C6_atom_dict, C6_CO=C6_CO, RA=RA, bAlloc=False, return_grads=False
            )

    # Original Python implementation (backup/reference)
    nx, ny, nz = [int(i) for i in ngrid[:3]]
    xs = origin[0] + np.arange(nx)*step
    ys = origin[1] + np.arange(ny)*step
    zs = origin[2] + np.arange(nz)*step
    XX, YY, ZZ = np.meshgrid(xs, ys, zs, indexing='ij')
    E_vdw = np.zeros((nx, ny, nz), dtype=np.float32)

    C6_atom = np.array([C6_atom_dict.get(z, 1.0) for z in atomTypes])
    RA2 = RA**2
    for ia in range(len(atomPos)):
        r2 = (XX-atomPos[ia,0])**2 + (YY-atomPos[ia,1])**2 + (ZZ-atomPos[ia,2])**2
        E_vdw -= np.sqrt(C6_atom[ia]*C6_CO) / (r2 + RA2)**3

    if return_grads:
        grads_E_vdw = np.stack([np.gradient(E_vdw, step, axis=i) for i in range(3)], axis=-1)
        return E_vdw, grads_E_vdw
    return E_vdw

