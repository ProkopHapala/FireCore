import os
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyopencl as cl
import importlib.util

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.MMFF import MMFF as MMFF_pyocl
from pyBall.OCL.MolecularDynamics import MolecularDynamics
from pyBall.tests import ocl_GridFF_new as gff_ocl


DEFAULT_MOL_TYPE_MAP = {'C': 'C_R', 'O': 'O_2', 'H': 'H'}


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def write_xyz_frame(fout, enames, pos3, comment=""):
    nat = pos3.shape[0]
    fout.write(f"{nat}\n")
    fout.write(f"{comment}\n")
    for i in range(nat):
        fout.write(f"{enames[i]:3s} {pos3[i,0]:12.6f} {pos3[i,1]:12.6f} {pos3[i,2]:12.6f}\n")


def load_gridff_array(path):
    arr = np.load(path)
    if arr.ndim != 4:
        raise ValueError(f"GridFF must be 4D, got {arr.shape} from {path}")
    if arr.shape[3] == 3:
        arr4 = np.zeros(arr.shape[:3] + (4,), dtype=np.float32)
        arr4[:, :, :, :3] = arr.astype(np.float32)
        return arr4
    if arr.shape[3] == 4:
        return np.ascontiguousarray(arr.astype(np.float32))
    raise ValueError(f"GridFF channels must be 3 or 4, got {arr.shape} from {path}")


def find_generated_gridff(workdir, src_xyz):
    bname = os.path.splitext(os.path.basename(src_xyz))[0]
    cands = [
        os.path.join(workdir, 'data', bname, 'Bspline_PLQd.npy'),
        os.path.join(workdir, 'data', bname, 'Bspline_PLQd_ocl.npy'),
    ]
    for p in cands:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"Generated GridFF not found, tried {cands}")


def _load_caf2_rect_module():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    mod_path = os.path.join(root, 'scripts', 'caf2_rect_surface.py')
    spec = importlib.util.spec_from_file_location('caf2_rect_surface', mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def ensure_rectangular_caf2_xyz(src_xyz, out_dir, nx=2, nz=1, layers=2, layers_from='top'):
    src_xyz = os.path.abspath(src_xyz)
    bname = os.path.basename(src_xyz)
    if 'CaF2' not in bname:
        return src_xyz
    if ('generated_rect' in src_xyz) or ('_rect_' in os.path.basename(src_xyz)):
        return src_xyz
    mod = _load_caf2_rect_module()
    ensure_dir(out_dir)
    stem = os.path.splitext(os.path.basename(src_xyz))[0]
    prefix = f'{stem}_rect_nx{int(nx)}_nz{int(nz)}_L{int(layers)}_{layers_from}'
    xyz_path = os.path.join(out_dir, prefix + '.xyz')
    if os.path.exists(xyz_path):
        return xyz_path
    sys0 = AtomicSystem(fname=src_xyz, bPreinit=False)
    base, lvec_rect0, info = mod.build_rectangular_base(sys0, coeff_search=4)
    layer_sel = mod.select_layers(base['apos'], base['enames'], base['atypes'], base['qs'], n_keep=layers, keep_from=layers_from, ztol=0.12, gap_scale=1.5)
    vac0 = float(np.linalg.norm(sys0.lvec[2]) - (np.max(sys0.apos[:, 2]) - np.min(sys0.apos[:, 2])))
    vacuum = max(0.0, vac0)
    apos_trim, lvec_trim, thickness = mod.make_output_lvec(lvec_rect0, layer_sel['apos'], vacuum=vacuum)
    apos_out, enames_out, atypes_out, qs_out, lvec_out = mod.replicate_rectangular(apos_trim, layer_sel['enames'], layer_sel['atypes'], layer_sel['qs'], lvec_trim, nx=nx, ny=nz)
    qs_out = mod.assign_caf2_charges(enames_out, 1.0)
    sys_out = AtomicSystem(apos=apos_out, atypes=atypes_out, enames=enames_out, lvec=lvec_out, qs=qs_out, bPreinit=False)
    sys_out.saveXYZ(xyz_path, blvec=True, bQs=True)
    mod.plot_structure(apos_out, enames_out, lvec_out, os.path.join(out_dir, prefix + '.png'), title=prefix)
    return xyz_path


def ensure_gridff_file(src_xyz, out_dir, gridff_path=None, dg=(23.175/225.0, 20.070/200.0, 48.472/382.0), save_name='double3', job='PLQ', use_CG=True, nmaxiter=1000, nPerStep=25, damp=0.15, sigma=0.0, alpha_morse=1.5, atom_types_path=None, element_types_path=None):
    if gridff_path is not None and os.path.exists(gridff_path):
        return os.path.abspath(gridff_path)
    ensure_dir(out_dir)
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if atom_types_path is None:
        atom_types_path = os.path.join(root, 'cpp', 'common_resources', 'AtomTypes.dat')
    if element_types_path is None:
        element_types_path = os.path.join(root, 'cpp', 'common_resources', 'ElementTypes.dat')
    cwd = os.getcwd()
    try:
        os.chdir(out_dir)
        gff_ocl.test_gridFF_ocl(
            fname=os.path.abspath(src_xyz),
            Atom_Types_name=os.path.abspath(atom_types_path),
            Element_Types_name=os.path.abspath(element_types_path),
            job=job,
            save_name=save_name,
            use_CG=bool(use_CG),
            nmaxiter=int(nmaxiter),
            nPerStep=int(nPerStep),
            damp=float(damp),
            save_fig=False,
            dg=tuple(float(x) for x in dg),
            sigma=float(sigma),
            alpha_morse=float(alpha_morse),
        )
        return find_generated_gridff(out_dir, src_xyz)
    finally:
        os.chdir(cwd)


def prepare_molecule(mol_path, type_map=None, z_shift=0.0, find_bonds=True):
    mol = AtomicSystem(mol_path)
    mol.apos = np.asarray(mol.apos, dtype=np.float32)
    mol.enames = np.array([str(e).split('_')[0] for e in mol.enames], dtype=object)
    if z_shift != 0.0:
        mol.apos[:, 2] += np.float32(z_shift)
    if find_bonds and getattr(mol, 'bonds', None) is None:
        mol.findBonds()
    mol.neighs()
    tmap = dict(DEFAULT_MOL_TYPE_MAP)
    if type_map is not None:
        tmap.update(type_map)
    mol.atom_types_mmff = [tmap.get(e, e) for e in mol.enames]
    return mol


def find_corner_oxygen_indices(mol):
    en = np.asarray(mol.enames)
    ios = np.where(en == 'O')[0]
    if len(ios) < 2:
        raise ValueError('Need at least two oxygen atoms in molecule to define anchor/opposite oxygen')
    xy_sum = mol.apos[ios, 0] + mol.apos[ios, 1]
    i_anchor = int(ios[np.argmin(xy_sum)])
    i_opposite = int(ios[np.argmax(xy_sum)])
    return i_anchor, i_opposite


def download_state(md, natoms, nvecs):
    apos_buf = np.empty((nvecs, 4), dtype=np.float32)
    af_buf = np.empty((nvecs, 4), dtype=np.float32)
    cl.enqueue_copy(md.queue, apos_buf, md.buffer_dict['apos'])
    cl.enqueue_copy(md.queue, af_buf, md.buffer_dict['aforce'])
    md.queue.finish()
    return apos_buf[:natoms, :3].copy(), af_buf[:natoms, :3].copy(), af_buf[:natoms, 3].copy()


class GridFFRelaxedScan:
    def __init__(self, mol_path, sub_xyz_path, gridff_path=None, out_dir='out_gridff_relaxed_scan', mol_type_map=None, grid_p0=(0.0, 0.0, 0.0), grid_step=(23.175/225.0, 20.070/200.0, 48.472/382.0), grid_alpha=1.5, grid_sigma=0.0, dt=0.02, damp=0.01, Fconv=1e-4, nstep_max=2000, out_stride=100, anchor_k=2000.0):
        self.mol_path = os.path.abspath(mol_path)
        self.sub_xyz_path = os.path.abspath(sub_xyz_path)
        self.gridff_path = None if gridff_path is None else os.path.abspath(gridff_path)
        self.out_dir = os.path.abspath(out_dir)
        self.grid_p0 = np.array(grid_p0, dtype=np.float32)
        self.grid_step = np.array(grid_step, dtype=np.float32)
        self.grid_alpha = float(grid_alpha)
        self.grid_sigma = float(grid_sigma)
        self.dt = float(dt)
        self.damp = float(damp)
        self.Fconv = float(Fconv)
        self.nstep_max = int(nstep_max)
        self.out_stride = int(out_stride)
        self.anchor_k = float(anchor_k)
        self.mol_type_map = dict(DEFAULT_MOL_TYPE_MAP)
        if mol_type_map is not None:
            self.mol_type_map.update(mol_type_map)
        self.sub = None
        self.mol = None
        self.mm = None
        self.md = None
        self.gridff = None
        self.anchor_idx = None
        self.opposite_idx = None
        self.surface_top_z = None

    def prepare(self, anchor_z_above=4.0, lateral_shift=(0.0, 0.0), generate_grid=True):
        ensure_dir(self.out_dir)
        rect_dir = os.path.join(self.out_dir, 'generated_rect')
        self.sub_xyz_path = ensure_rectangular_caf2_xyz(self.sub_xyz_path, rect_dir, nx=2, nz=1, layers=2, layers_from='top')
        self.sub = AtomicSystem(self.sub_xyz_path, bPreinit=False)
        self.surface_top_z = float(np.max(self.sub.apos[:, 2]))
        self.gridff_path = ensure_gridff_file(
            src_xyz=self.sub_xyz_path,
            out_dir=os.path.join(self.out_dir, 'gridff_cache'),
            gridff_path=self.gridff_path if not generate_grid else self.gridff_path,
            dg=tuple(float(x) for x in self.grid_step),
            sigma=self.grid_sigma,
            alpha_morse=self.grid_alpha,
        ) if (generate_grid or self.gridff_path is not None) else self.gridff_path
        if self.gridff_path is None or (not os.path.exists(self.gridff_path)):
            raise FileNotFoundError('GridFF path not available; provide --gridff or allow generation')
        self.gridff = load_gridff_array(self.gridff_path)
        z_shift = self.surface_top_z + float(anchor_z_above)
        self.mol = prepare_molecule(self.mol_path, type_map=self.mol_type_map, z_shift=z_shift)
        self.mol.apos[:, 0] += np.float32(lateral_shift[0])
        self.mol.apos[:, 1] += np.float32(lateral_shift[1])
        self.anchor_idx, self.opposite_idx = find_corner_oxygen_indices(self.mol)
        self.mm = MMFF_pyocl(bTorsion=False, verbosity=0)
        self.mm.capping_atoms = set()
        self.mm.reorder_nodes_first = False
        self.mm.nPBC = (0, 0, 0)
        self.mm.lvec = np.eye(3, dtype=np.float32) * 100.0
        self.mm.toMMFFsp3_loc(mol=self.mol, atom_types=self.mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
        if getattr(self.mol, 'qs', None) is not None:
            qs = np.asarray(self.mol.qs, dtype=np.float32)
            if len(qs) == self.mm.natoms:
                self.mm.REQs[:, 2] = qs
        self.mm.make_back_neighs(b_cap_neighs=False)
        self.mm.excl = self.mm._make_excl_1_2_3(self.mm.neighs, neighCell=self.mm.neighCell, npbc=self.mm.npbc, EXCL_MAX=16)
        self.md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0', enable_nonbond=True)
        self.md.realloc(self.mm, nSystems=1)
        self.md.setup_kernels()
        self.md.pack_system(0, self.mm)
        self.md.initGridFF(
            grid_shape=tuple(int(x) for x in self.gridff.shape[:3]),
            bspline_data=np.ascontiguousarray(self.gridff),
            grid_p0=tuple(float(x) for x in self.grid_p0),
            grid_step=tuple(float(x) for x in self.grid_step),
            use_texture=False,
            r_damp=0.0,
            alpha_morse=self.grid_alpha,
            bKernels=True,
        )
        if getattr(self.md, 'kernel_args_getNonBond_GridFF_Bspline_ex2', None) is None:
            raise RuntimeError('GridFF ex2 kernel is not initialized')
        self.zero_velocities()
        return self

    def zero_velocities(self):
        self.md.toGPU('avel', np.zeros((self.mm.nvecs, 4), dtype=np.float32))

    def get_anchor_position(self):
        pos, _, _ = download_state(self.md, self.mm.natoms, self.mm.nvecs)
        return pos[self.anchor_idx].copy()

    def set_positions(self, pos3, zero_vel=True):
        apos_buf = np.zeros((self.mm.nvecs, 4), dtype=np.float32)
        apos_buf[:self.mm.natoms, :3] = np.asarray(pos3, dtype=np.float32)
        self.md.toGPU('apos', apos_buf)
        if zero_vel:
            self.zero_velocities()

    def set_constraint(self, target_pos):
        nat = self.mm.natoms
        constr = np.zeros((nat, 4), dtype=np.float32)
        constrK = np.zeros((nat, 4), dtype=np.float32)
        constr[self.anchor_idx, :3] = np.asarray(target_pos, dtype=np.float32)
        constr[self.anchor_idx, 3] = 1.0
        constrK[self.anchor_idx, 0] = self.anchor_k
        constrK[self.anchor_idx, 1] = self.anchor_k
        constrK[self.anchor_idx, 2] = self.anchor_k
        self.md.toGPU('constr', constr)
        self.md.toGPU('constrK', constrK)
        vfac = np.float32(1.0 - self.damp)
        self.md.toGPU('MDparams', np.array([self.dt, 1e+6, vfac, 0.0], dtype=np.float32), byte_offset=0)

    def relax_to_constraint(self, target_pos, verbose=False, batch_size=100):
        self.set_constraint(target_pos)
        hist = []
        pos = None
        F = None
        istep = 0
        while istep < self.nstep_max:
            # Run batch of OpenCL iterations
            batch_end = min(istep + batch_size, self.nstep_max)
            for _ in range(batch_size):
                self.md.run_cleanForceMMFFf4()
                self.md.run_getNonBond_GridFF_Bspline_ex2()
                self.md.run_getMMFFf4()
                self.md.run_updateAtomsMMFFf4()
                istep += 1
            
            # Download and check after batch
            want = (istep == batch_size) or (istep == self.nstep_max) or ((istep % self.out_stride) == 0)
            if want:
                pos, F, E = download_state(self.md, self.mm.natoms, self.mm.nvecs)
                fmag = np.sqrt(np.sum(F * F, axis=1))
                fmask = np.ones(len(fmag), dtype=bool)
                fmask[self.anchor_idx] = False
                fmax = float(fmag[fmask].max()) if np.any(fmask) else float(fmag.max())
                hist.append((istep, fmax))
                if verbose:
                    print(f'relax step {istep:6d} fmax={fmax:12.6e}')
                if not np.isfinite(pos).all():
                    raise RuntimeError(f'NaN/Inf positions at relax step {istep}')
                if fmax < self.Fconv:
                    return pos, F, E, istep, hist
        if pos is None:
            pos, F, E = download_state(self.md, self.mm.natoms, self.mm.nvecs)
        return pos, F, E, self.nstep_max, hist

    def make_linear_anchor_path(self, x0, x1, y, z, nscan):
        xs = np.linspace(float(x0), float(x1), int(nscan) + 1)
        path = np.zeros((len(xs), 3), dtype=np.float32)
        path[:, 0] = xs
        path[:, 1] = float(y)
        path[:, 2] = float(z)
        return path

    def run_path_scan(self, anchor_path, xyz_path=None, warm_start=True, verbose=False):
        anchor_path = np.asarray(anchor_path, dtype=np.float32)
        nscan = len(anchor_path)
        A_scan = np.zeros((nscan, self.mm.natoms, 3), dtype=np.float32)
        F_scan = np.zeros((nscan, self.mm.natoms, 3), dtype=np.float32)
        opposite_path = np.zeros((nscan, 3), dtype=np.float32)
        anchor_real = np.zeros((nscan, 3), dtype=np.float32)
        iters = np.zeros(nscan, dtype=np.int32)
        fmaxs = np.zeros(nscan, dtype=np.float32)
        if xyz_path is None:
            xyz_path = os.path.join(self.out_dir, 'scan_relaxed.xyz')
        fout = open(xyz_path, 'w')
        try:
            if not warm_start:
                pos0 = self.mol.apos.copy()
            for i, target in enumerate(anchor_path):
                if (not warm_start) and i > 0:
                    shift = target - anchor_path[0]
                    pos = pos0.copy()
                    pos[:, :3] += shift[None, :3]
                    self.set_positions(pos, zero_vel=True)
                pos, F, E, niter, hist = self.relax_to_constraint(target, verbose=verbose)
                fmag = np.sqrt(np.sum(F * F, axis=1))
                fmask = np.ones(len(fmag), dtype=bool)
                fmask[self.anchor_idx] = False
                A_scan[i] = pos
                F_scan[i] = F
                opposite_path[i] = pos[self.opposite_idx]
                anchor_real[i] = pos[self.anchor_idx]
                iters[i] = niter
                fmaxs[i] = float(fmag[fmask].max()) if np.any(fmask) else float(fmag.max())
                comment = f'# scan={i} target=({target[0]:.6f},{target[1]:.6f},{target[2]:.6f}) anchor=({anchor_real[i,0]:.6f},{anchor_real[i,1]:.6f},{anchor_real[i,2]:.6f}) opposite=({opposite_path[i,0]:.6f},{opposite_path[i,1]:.6f},{opposite_path[i,2]:.6f}) niter={niter} fmax={fmaxs[i]:.6e}'
                write_xyz_frame(fout, self.mol.enames, pos, comment=comment)
            return {
                'A_scan': A_scan,
                'F_scan': F_scan,
                'anchor_target_path': anchor_path,
                'anchor_real_path': anchor_real,
                'opposite_path': opposite_path,
                'iters': iters,
                'fmaxs': fmaxs,
                'anchor_idx': int(self.anchor_idx),
                'opposite_idx': int(self.opposite_idx),
                'xyz_path': xyz_path,
            }
        finally:
            fout.close()

    def plot_gridff_xz(self, out_png, iy=None, plq=(1.0, 1.0, 0.0), vmax=0.1, test_req=None, test_particle='H', reference_z=19.5):
        dat = self.get_gridff_xz_slice(iy=iy, plq=plq, reference_z=reference_z)
        fig, axes = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True)
        names = ['VPaul', 'VLond', 'VCoul', 'VTotal']
        fields = [dat['VPaul'], dat['VLond'], dat['VCoul'], dat['VTotal']]
        
        # Test particle info with REQ parameters
        if test_req:
            test_info = (
                f"Test particle: {test_particle}\n"
                f"R={test_req['R']:.2f} A\n"
                f"E={test_req['E']:.3f} eV\n"
                f"Q={test_req['Q']:+.2f} e\n"
                f"PLQ: P={plq[0]:.2f} L={plq[1]:.2f} Q={plq[2]:.2f}"
            )
        else:
            test_info = f"Test particle: {test_particle}\nPLQ: P={plq[0]:.2f} L={plq[1]:.2f} Q={plq[2]:.2f}"
        
        for ax, name, field in zip(axes, names, fields):
            im = ax.pcolormesh(dat['x_edges'], dat['z_edges'], field.T, shading='auto', cmap='bwr', vmin=-vmax, vmax=vmax)
            ax.set_title(f'{name} xz iy={dat["iy"]}')
            ax.set_xlabel('x [A]')
            ax.set_ylabel('z [A]')
            ax.set_aspect('equal')
            fig.colorbar(im, ax=ax, fraction=0.046)
            # Add horizontal line for XY slice height
            ax.axhline(y=self.plot_xy_z if hasattr(self, 'plot_xy_z') else 3.0, color='green', linestyle='--', alpha=0.7, linewidth=1)
        
        fig.text(0.02, 0.98, test_info, transform=fig.transFigure, fontsize=8, verticalalignment='top', family='monospace')
        fig.text(0.02, 0.02, f"Green line: XY slice at z={self.plot_xy_z if hasattr(self, 'plot_xy_z') else 3.0:.1f}A", transform=fig.transFigure, fontsize=8, verticalalignment='bottom', family='monospace')
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        return out_png

    def plot_gridff_xy(self, out_png, z=7.0, plq=(1.0, 1.0, 0.0), vmax=0.1, test_req=None, test_particle='H', reference_z=19.5):
        """Plot XY slice of GridFF at specified height with Coulomb reference subtraction"""
        self.plot_xy_z = z  # Store for crosshair reference
        dat = self.get_gridff_xy_slice(z=z, plq=plq, reference_z=reference_z)
        fig, axes = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True)
        names = ['VPaul', 'VLond', 'VCoul', 'VTotal']
        fields = [dat['VPaul'], dat['VLond'], dat['VCoul'], dat['VTotal']]
        
        # Test particle info
        if test_req:
            test_info = (
                f"Test particle: {test_particle}\n"
                f"R={test_req['R']:.2f} A\n"
                f"E={test_req['E']:.3f} eV\n"
                f"Q={test_req['Q']:+.2f} e\n"
                f"PLQ: P={plq[0]:.2f} L={plq[1]:.2f} Q={plq[2]:.2f}"
            )
        else:
            test_info = f"Test particle: {test_particle}\nPLQ: P={plq[0]:.2f} L={plq[1]:.2f} Q={plq[2]:.2f}"
        
        for ax, name, field in zip(axes, names, fields):
            im = ax.pcolormesh(dat['x_edges'], dat['y_edges'], field.T, shading='auto', cmap='bwr', vmin=-vmax, vmax=vmax)
            ax.set_title(f'{name} xy z={z:.1f}A')
            ax.set_xlabel('x [A]')
            ax.set_ylabel('y [A]')
            ax.set_aspect('equal')
            fig.colorbar(im, ax=ax, fraction=0.046)
            # Add vertical line for XZ slice position
            ax.axvline(x=dat.get('xz_y', self.grid_p0[1]), color='red', linestyle='--', alpha=0.7, linewidth=1)
        
        fig.text(0.02, 0.98, test_info, transform=fig.transFigure, fontsize=8, verticalalignment='top', family='monospace')
        fig.text(0.02, 0.02, f"Red line: XZ slice at y={dat.get('xz_y', self.grid_p0[1]):.1f}A", 
                transform=fig.transFigure, fontsize=7, color='red', verticalalignment='bottom')
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        return out_png

    def plot_gridff_xz_arbitrary(self, out_png, point1, point2, plq=(1.0, 1.0, 0.0), vmax=0.1, test_req=None, test_particle='H', n_points=100):
        """Plot GridFF XZ slice along arbitrary line defined by two points (y,z) coordinates"""
        # Convert points to numpy arrays
        p1 = np.array(point1, dtype=float)
        p2 = np.array(point2, dtype=float)
        
        # Generate points along the line
        t = np.linspace(0, 1, n_points)
        line_y = p1[0] + t * (p2[0] - p1[0])  # y coordinates
        line_z = p1[1] + t * (p2[1] - p1[1])  # z coordinates
        
        # Get x range from grid
        x_range = self.grid_p0[0] + np.arange(self.gridff.shape[0]) * self.grid_step[0]
        
        # Create 2D array for the slice
        VPaul = np.zeros((len(x_range), n_points))
        VLond = np.zeros((len(x_range), n_points))
        VCoul = np.zeros((len(x_range), n_points))
        
        # Sample GridFF along the arbitrary line
        for i, (y, z) in enumerate(zip(line_y, line_z)):
            # Find nearest grid indices
            iy = int(np.clip(np.round((y - self.grid_p0[1]) / self.grid_step[1]), 0, self.gridff.shape[1] - 1))
            iz = int(np.clip(np.round((z - self.grid_p0[2]) / self.grid_step[2]), 0, self.gridff.shape[2] - 1))
            
            # Extract slice at this (y,z) point
            VPaul[:, i] = self.gridff[:, iy, iz, 0]
            VLond[:, i] = self.gridff[:, iy, iz, 1]
            VCoul[:, i] = self.gridff[:, iy, iz, 2]
        
        # Apply PLQ scaling
        VTotal = VPaul * plq[0] + VLond * plq[1] + VCoul * plq[2]
        
        # Create meshgrid for plotting
        X, Y_line = np.meshgrid(x_range, t)
        
        # Plot
        fig, axes = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True)
        names = ['VPaul', 'VLond', 'VCoul', 'VTotal']
        fields = [VPaul.T, VLond.T, VCoul.T, VTotal.T]
        
        # Test particle info
        if test_req:
            test_info = (
                f"Test particle: {test_particle}\n"
                f"R={test_req['R']:.2f} A\n"
                f"E={test_req['E']:.3f} eV\n"
                f"Q={test_req['Q']:+.2f} e\n"
                f"PLQ: P={plq[0]:.2f} L={plq[1]:.2f} Q={plq[2]:.2f}\n"
                f"Slice: ({p1[0]:.1f},{p1[1]:.1f})→({p2[0]:.1f},{p2[1]:.1f})"
            )
        else:
            test_info = f"Test particle: {test_particle}\nPLQ: P={plq[0]:.2f} L={plq[1]:.2f} Q={plq[2]:.2f}"
        
        for ax, name, field in zip(axes, names, fields):
            im = ax.pcolormesh(X, Y_line, field, shading='auto', cmap='bwr', vmin=-vmax, vmax=vmax)
            ax.set_title(f'{name} arbitrary slice')
            ax.set_xlabel('x [A]')
            ax.set_ylabel('parametric position')
            ax.set_aspect('equal')
            fig.colorbar(im, ax=ax, fraction=0.046)
        
        fig.text(0.02, 0.98, test_info, transform=fig.transFigure, fontsize=8, verticalalignment='top', family='monospace')
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        return out_png

    def get_gridff_xy_slice(self, z=3.0, plq=(1.0, 1.0, 0.0), reference_z=20.0):
        """Get XY slice of GridFF at specified height with Coulomb reference subtraction"""
        plq = np.asarray(plq, dtype=float)
        iz = int(np.clip(np.round((z - self.grid_p0[2]) / self.grid_step[2]), 0, self.gridff.shape[2] - 1))
        iz_ref = int(np.clip(np.round((reference_z - self.grid_p0[2]) / self.grid_step[2]), 0, self.gridff.shape[2] - 1))
        x_edges = self.grid_p0[0] + np.arange(self.gridff.shape[0] + 1) * self.grid_step[0]
        y_edges = self.grid_p0[1] + np.arange(self.gridff.shape[1] + 1) * self.grid_step[1]
        VPaul = self.gridff[:, :, iz, 0].copy()
        VLond = self.gridff[:, :, iz, 1].copy()
        VCoul = self.gridff[:, :, iz, 2].copy()
        VCoul_ref = self.gridff[:, :, iz_ref, 2].copy()  # Reference at large distance
        VCoul_corrected = VCoul - VCoul_ref  # Subtract reference
        VTotal = VPaul * plq[0] + VLond * plq[1] + VCoul_corrected * plq[2]
        return {
            'iz': int(iz),
            'x_edges': x_edges,
            'y_edges': y_edges,
            'VPaul': VPaul,
            'VLond': VLond,
            'VCoul': VCoul_corrected,  # Return corrected Coulomb
            'VTotal': VTotal,
            'plq': np.asarray(plq, dtype=np.float32),
            'xz_y': self.get_anchor_position()[1] if hasattr(self, 'anchor_idx') else self.grid_p0[1],
        }

    def get_gridff_xz_slice(self, iy=-1, plq=(1.0, 1.0, 0.0), reference_z=20.0):
        """Get XZ slice of GridFF at specified Y index with Coulomb reference subtraction"""
        plq = np.asarray(plq, dtype=float)
        if iy is None:
            ay = self.get_anchor_position()[1]
            iy = int(np.clip(np.round((ay - self.grid_p0[1]) / self.grid_step[1]), 0, self.gridff.shape[1] - 1))
        iy = np.clip(iy, 0, self.gridff.shape[1] - 1)
        x_edges = self.grid_p0[0] + np.arange(self.gridff.shape[0] + 1) * self.grid_step[0]
        z_edges = self.grid_p0[2] + np.arange(self.gridff.shape[2] + 1) * self.grid_step[2]
        VPaul = self.gridff[:, iy, :, 0].copy()
        VLond = self.gridff[:, iy, :, 1].copy()
        VCoul = self.gridff[:, iy, :, 2].copy()
        
        # Reference Coulomb potential at large distance (z=20A) to remove surface dipole offset
        iz_ref = int(np.clip(np.round((reference_z - self.grid_p0[2]) / self.grid_step[2]), 0, self.gridff.shape[2] - 1))
        VCoul_ref = VCoul[:, iz_ref:iz_ref+1]  # Reference slice
        VCoul_corrected = VCoul - VCoul_ref  # Subtract reference
        
        VTotal = VPaul * plq[0] + VLond * plq[1] + VCoul_corrected * plq[2]
        return {
            'iy': int(iy),
            'x_edges': x_edges,
            'z_edges': z_edges,
            'VPaul': VPaul,
            'VLond': VLond,
            'VCoul': VCoul_corrected,  # Return corrected Coulomb
            'VTotal': VTotal,
            'plq': np.asarray(plq, dtype=np.float32),
            'xy_x': self.get_anchor_position()[0] if hasattr(self, 'anchor_idx') else self.grid_p0[0],
        }

    def step_relaxation(self, target_pos, nsub=1, verbose=False):
        self.set_constraint(target_pos)
        pos = None
        F = None
        E = None
        for _ in range(int(nsub)):
            self.md.run_cleanForceMMFFf4()
            self.md.run_getNonBond_GridFF_Bspline_ex2()
            self.md.run_getMMFFf4()
            self.md.run_updateAtomsMMFFf4()
        pos, F, E = download_state(self.md, self.mm.natoms, self.mm.nvecs)
        fmag = np.sqrt(np.sum(F * F, axis=1))
        fmask = np.ones(len(fmag), dtype=bool)
        fmask[self.anchor_idx] = False
        fmax = float(fmag[fmask].max()) if np.any(fmask) else float(fmag.max())
        if verbose:
            print(f'step_relaxation() nsub={nsub} fmax={fmax:12.6e}')
        return {
            'pos': pos,
            'force': F,
            'energy_terms': E,
            'fmax': fmax,
            'anchor_idx': int(self.anchor_idx),
            'opposite_idx': int(self.opposite_idx),
            'anchor_target': np.asarray(target_pos, dtype=np.float32).copy(),
            'anchor_real': pos[self.anchor_idx].copy(),
            'opposite_real': pos[self.opposite_idx].copy(),
        }

    def plot_opposite_trajectory(self, scan_result, out_png):
        opp = np.asarray(scan_result['opposite_path'])
        anc = np.asarray(scan_result['anchor_real_path'])
        tgt = np.asarray(scan_result['anchor_target_path'])
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
        axes[0].plot(opp[:, 0], opp[:, 2], '-o', ms=3, lw=1.5, label='opposite O')
        axes[0].plot(anc[:, 0], anc[:, 2], '-o', ms=3, lw=1.0, label='anchor real')
        axes[0].plot(tgt[:, 0], tgt[:, 2], '--', lw=1.5, label='anchor target')
        axes[0].set_xlabel('x [A]')
        axes[0].set_ylabel('z [A]')
        axes[0].grid(True)
        axes[0].legend()
        axes[1].plot(opp[:, 0], opp[:, 1], '-o', ms=3, lw=1.5, label='opposite O')
        axes[1].plot(anc[:, 0], anc[:, 1], '-o', ms=3, lw=1.0, label='anchor real')
        axes[1].plot(tgt[:, 0], tgt[:, 1], '--', lw=1.5, label='anchor target')
        axes[1].set_xlabel('x [A]')
        axes[1].set_ylabel('y [A]')
        axes[1].grid(True)
        axes[1].legend()
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        return out_png

    def save_scan_npz(self, scan_result, out_npz):
        np.savez(out_npz, **scan_result)
        return out_npz
