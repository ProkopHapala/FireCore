import os
import math
import ast
import struct
import numpy as np

from pyBall.OCL import MMparams
from pyBall.OCL.InteractionEnergy import InteractionScanner, load_xyz_with_REQs
from pyBall.OCL.MolecularDynamics import MolecularDynamics, FOLDED_BASIS_MAX
from pyBall.tests.utils import getPLQH

_BASE_PATH = os.path.dirname(os.path.abspath(__file__))
_DATA_PATH = os.path.join(_BASE_PATH, '../cpp/common_resources')


def _load_atom_types():
    etypes = MMparams.read_element_types(os.path.join(_DATA_PATH, 'ElementTypes.dat'))
    return MMparams.read_atom_types(os.path.join(_DATA_PATH, 'AtomTypes.dat'), etypes)


_ATOM_TYPES = None


def get_atom_types():
    global _ATOM_TYPES
    if _ATOM_TYPES is None:
        _ATOM_TYPES = _load_atom_types()
    return _ATOM_TYPES


class ProbeParticle:
    def __init__(self, name='H', charge=0.2, alpha=2.0, E0=None, R0=None, H=0.0):
        ats = get_atom_types()
        if name not in ats:
            raise KeyError(f"ProbeParticle: atom type '{name}' not found in AtomTypes.dat")
        at = ats[name]
        self.name = str(name)
        self.charge = float(charge)
        self.alpha = float(alpha)
        self.R0 = float(at.RvdW if R0 is None else R0)
        self.E0 = float(at.EvdW if E0 is None else E0)
        self.H = float(H)

    def req(self):
        return np.array([[self.R0, math.sqrt(self.E0), self.charge, self.H]], dtype=np.float32)

    def plqh(self):
        return np.asarray(getPLQH(R0=self.R0, E0=self.E0, a=self.alpha, Q=self.charge, H=self.H), dtype=np.float32)

    def as_dict(self):
        return {'name': self.name, 'charge': self.charge, 'alpha': self.alpha, 'R0': self.R0, 'E0': self.E0, 'H': self.H}


def infer_gridff_xyz(gridff_path):
    gridff_path = os.path.abspath(gridff_path)
    base = os.path.basename(os.path.dirname(gridff_path))
    cand = os.path.join(_DATA_PATH, 'xyz', base + '.xyz')
    if os.path.exists(cand):
        return cand
    raise FileNotFoundError(f"infer_gridff_xyz(): cannot infer xyz for grid '{gridff_path}', tried '{cand}'. Pass xyz_path explicitly")


def load_npy_relaxed(path):
    try:
        return np.load(path)
    except ValueError as e:
        msg = str(e)
        if 'Cannot parse header' not in msg:
            raise
    with open(path, 'rb') as f:
        magic = f.read(6)
        if magic != b'\x93NUMPY':
            raise ValueError(f"load_npy_relaxed(): '{path}' is not a .npy file")
        ver = struct.unpack('BB', f.read(2))
        if ver[0] == 1:
            hlen = struct.unpack('<H', f.read(2))[0]
        elif ver[0] in (2, 3):
            hlen = struct.unpack('<I', f.read(4))[0]
        else:
            raise ValueError(f"load_npy_relaxed(): unsupported npy version {ver} in '{path}'")
        header = f.read(hlen).decode('latin1').strip()
        if not header.endswith('}'):
            i = header.rfind('}')
            if i < 0:
                raise ValueError(f"load_npy_relaxed(): malformed header in '{path}': {header!r}")
            header = header[:i + 1]
        meta = ast.literal_eval(header)
        shape = tuple(int(x) for x in meta['shape'])
        dt = np.dtype(meta['descr'])
        count = int(np.prod(shape))
        arr = np.fromfile(f, dtype=dt, count=count)
        if arr.size != count:
            raise ValueError(f"load_npy_relaxed(): expected {count} items in '{path}', got {arr.size}")
        arr = arr.reshape(shape, order='F' if bool(meta.get('fortran_order', False)) else 'C')
        return arr


def _identity_transforms(points):
    pts = np.asarray(points, dtype=np.float32).reshape(-1, 3)
    T = np.zeros((len(pts), 3, 4), dtype=np.float32)
    T[:, 0, 0] = 1.0
    T[:, 1, 1] = 1.0
    T[:, 2, 2] = 1.0
    T[:, :, 3] = pts
    return T.reshape(-1, 12)


def _component_plqh(component=None):
    key = 'total' if component is None else str(component).lower()
    if key == 'total':
        return np.array([1.0, 1.0, 1.0, 0.0], dtype=np.float32)
    if key == 'nonel':
        return np.array([1.0, 1.0, 0.0, 0.0], dtype=np.float32)
    if key == 'pauli':
        return np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float32)
    if key == 'london':
        return np.array([0.0, 1.0, 0.0, 0.0], dtype=np.float32)
    if key == 'coulomb':
        return np.array([0.0, 0.0, 1.0, 0.0], dtype=np.float32)
    raise ValueError(f"_component_plqh(): unsupported component '{component}'")


class GridFFProbeBackend:
    component_names = ('pauli', 'london', 'coulomb', 'nonel', 'total')

    def __init__(self, gridff_path, probe, xyz_path=None, z0=None, chunk_size=8192, use_gpu=True, type_map=None):
        self.gridff_path = os.path.abspath(gridff_path)
        self.xyz_path = os.path.abspath(xyz_path or infer_gridff_xyz(self.gridff_path))
        self.probe = probe
        self.type_map = dict(type_map or {})
        self.chunk_size = int(chunk_size)
        self.use_gpu = bool(use_gpu)
        self.grid = load_npy_relaxed(self.gridff_path)
        if self.grid.ndim != 4 or self.grid.shape[-1] < 3:
            raise ValueError(f"GridFFProbeBackend: expected PLQ grid shape (nx,ny,nz,>=3), got {self.grid.shape} from {self.gridff_path}")
        apos, _, enames, _, lvec = load_xyz_with_REQs(self.xyz_path, atom_types=get_atom_types(), type_map=self.type_map)
        if lvec is None:
            raise ValueError(f"GridFFProbeBackend: xyz '{self.xyz_path}' must contain lattice vectors")
        self.sub_apos = np.asarray(apos, dtype=np.float64)
        self.sub_enames = list(enames)
        self.lvec = np.asarray(lvec, dtype=np.float64)
        self.nx, self.ny, self.nz = self.grid.shape[:3]
        Lx = float(np.linalg.norm(self.lvec[0]))
        Ly = float(np.linalg.norm(self.lvec[1]))
        Lz = float(np.linalg.norm(self.lvec[2]))
        self.dg = np.array([Lx / self.nx, Ly / self.ny, Lz / self.nz], dtype=np.float64)
        self.z_top = float(np.max(self.sub_apos[:, 2]))
        z0_eff = float(z0) if z0 is not None else (0.0 if self.z_top >= 0.0 else (self.z_top - 3.0 * self.dg[2]))
        self.g0 = np.array([-0.5 * Lx, -0.5 * Ly, z0_eff], dtype=np.float64)
        self.plqh = probe.plqh().astype(np.float64)
        self.z_mode = 'height_above_top'
        self.md = None
        if self.use_gpu:
            grid4 = np.zeros(self.grid.shape[:3] + (4,), dtype=np.float32)
            grid4[..., :min(4, self.grid.shape[3])] = np.asarray(self.grid[..., :min(4, self.grid.shape[3])], dtype=np.float32)
            self.md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0')
            self.md.init_rigid_molecule_batch(np.zeros((1, 3), dtype=np.float32), probe.req(), nSystems=self.chunk_size)
            self.md.initGridFF(grid_shape=(self.nx, self.ny, self.nz), bspline_data=grid4, grid_p0=self.g0, grid_step=self.dg, use_texture=False, r_damp=0.0, alpha_morse=probe.alpha, bKernels=True)

    def make_world_zs(self, z_range, nz):
        hs = np.linspace(float(z_range[1]), float(z_range[0]), int(nz), dtype=np.float64)
        return self.z_top + hs

    def z_world_to_report(self, z_world):
        return np.asarray(z_world, dtype=np.float64) - self.z_top

    def _interp_channel(self, pts, ich):
        p = np.asarray(pts, dtype=np.float64).reshape(-1, 3)
        gx = (p[:, 0] - self.g0[0]) / self.dg[0]
        gy = (p[:, 1] - self.g0[1]) / self.dg[1]
        gz = (p[:, 2] - self.g0[2]) / self.dg[2]
        ix0 = np.floor(gx).astype(np.int64)
        iy0 = np.floor(gy).astype(np.int64)
        iz0 = np.floor(gz).astype(np.int64)
        tx = gx - ix0
        ty = gy - iy0
        tz = gz - iz0
        valid = (gz >= 0.0) & (gz <= (self.nz - 1))
        ix0 %= self.nx
        iy0 %= self.ny
        ix1 = (ix0 + 1) % self.nx
        iy1 = (iy0 + 1) % self.ny
        iz0 = np.clip(iz0, 0, self.nz - 2)
        iz1 = iz0 + 1
        g = self.grid[..., ich]
        c000 = g[ix0, iy0, iz0]
        c100 = g[ix1, iy0, iz0]
        c010 = g[ix0, iy1, iz0]
        c110 = g[ix1, iy1, iz0]
        c001 = g[ix0, iy0, iz1]
        c101 = g[ix1, iy0, iz1]
        c011 = g[ix0, iy1, iz1]
        c111 = g[ix1, iy1, iz1]
        c00 = c000 * (1.0 - tx) + c100 * tx
        c10 = c010 * (1.0 - tx) + c110 * tx
        c01 = c001 * (1.0 - tx) + c101 * tx
        c11 = c011 * (1.0 - tx) + c111 * tx
        c0 = c00 * (1.0 - ty) + c10 * ty
        c1 = c01 * (1.0 - ty) + c11 * ty
        out = c0 * (1.0 - tz) + c1 * tz
        out = np.asarray(out, dtype=np.float64)
        out[~valid] = np.nan
        return out

    def eval_components(self, points, components=('total', 'coulomb')):
        p = np.asarray(points, dtype=np.float64).reshape(-1, 3)
        out = {}
        need = set(str(c).lower() for c in components)
        if self.md is not None:
            if {'pauli', 'nonel', 'total'} & need:
                out['pauli'] = np.asarray(self.md.eval_gridff_points(p, PLQH=np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float32), chunk_size=self.chunk_size)['total'], dtype=np.float64) * self.plqh[0]
            if {'london', 'nonel', 'total'} & need:
                out['london'] = np.asarray(self.md.eval_gridff_points(p, PLQH=np.array([0.0, 1.0, 0.0, 0.0], dtype=np.float32), chunk_size=self.chunk_size)['total'], dtype=np.float64) * self.plqh[1]
            if {'coulomb', 'total'} & need:
                out['coulomb'] = np.asarray(self.md.eval_gridff_points(p, PLQH=np.array([0.0, 0.0, 1.0, 0.0], dtype=np.float32), chunk_size=self.chunk_size)['total'], dtype=np.float64) * self.plqh[2]
        else:
            if {'pauli', 'nonel', 'total'} & need:
                out['pauli'] = self._interp_channel(p, 0) * self.plqh[0]
            if {'london', 'nonel', 'total'} & need:
                out['london'] = self._interp_channel(p, 1) * self.plqh[1]
            if {'coulomb', 'total'} & need:
                out['coulomb'] = self._interp_channel(p, 2) * self.plqh[2]
        if 'nonel' in need or 'total' in need:
            out['nonel'] = out.get('pauli', 0.0) + out.get('london', 0.0)
        if 'total' in need:
            out['total'] = out.get('nonel', 0.0) + out.get('coulomb', 0.0)
        return out

    def eval_surface_map_gpu(self, x_range, y_range, z_range, nx=64, ny=64, nz=80, selector='total', color_component='coulomb', mode='threshold', threshold=0.0):
        if self.md is None:
            raise ValueError('GridFFProbeBackend.eval_surface_map_gpu(): GPU backend not initialized')
        sel = str(selector).lower()
        csel = str(color_component).lower()
        res = self.md.eval_surface_iso_gridff(
            probe_req=self.probe.req()[0],
            sel_PLQH=_component_plqh(sel) * self.probe.plqh(),
            col_PLQH=_component_plqh(csel) * self.probe.plqh(),
            x_range=x_range,
            y_range=y_range,
            z_range=(self.z_top + float(z_range[0]), self.z_top + float(z_range[1])),
            nx=nx, ny=ny, nz=nz, mode=mode, threshold=threshold, z_top=self.z_top,
        )
        return res


class XYZProbeBackend:
    component_names = ('total', 'nonel', 'coulomb', 'hbond', 'pauli', 'london')

    def __init__(self, substrate_xyz, probe, type_map=None, backend='reference', nPBC=(4, 4, 0), chunk_size=8192):
        self.substrate_xyz = os.path.abspath(substrate_xyz)
        self.probe = probe
        self.type_map = dict(type_map or {})
        self.backend = str(backend)
        self.nPBC = tuple(int(x) for x in nPBC)
        self.chunk_size = int(chunk_size)
        self.sub_apos, self.sub_REQs, self.sub_enames, _, self.lvec = load_xyz_with_REQs(self.substrate_xyz, atom_types=get_atom_types(), type_map=self.type_map)
        if self.lvec is None:
            raise ValueError(f"XYZProbeBackend: substrate '{self.substrate_xyz}' must contain lattice vectors")
        self.z_top = float(np.max(self.sub_apos[:, 2]))
        self.scanner = None
        self.md = None
        self._folded_ready = False
        self._bMacro = self._detect_rect_macro_ok()
        if self.backend == 'reference':
            self.scanner = InteractionScanner(nloc=32)
            self.scanner.set_molecule(np.zeros((1, 3), dtype=np.float32), probe.req())
            self.scanner.load_substrate_xyz(self.substrate_xyz, type_map=self.type_map)
            self.scanner.nPBC[:] = self.nPBC
            self.scanner.enable_macro = self._bMacro
            self.scanner.wrap_PBC = True
            self.scanner._update_macro_from_substrate()
        elif self.backend in ('fast_gpu', 'folded_gpu'):
            self.md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0')
            self.md.init_rigid_molecule_batch(np.zeros((1, 3), dtype=np.float32), probe.req(), nSystems=self.chunk_size)
            self.md.set_surface(self.substrate_xyz, nPBC=self.nPBC, pos0=(0.0, 0.0, 0.0), alpha_morse=probe.alpha, r_damp=0.0, bMacro=self._bMacro, type_map=self.type_map)
        else:
            raise ValueError(f"XYZProbeBackend: unsupported backend '{self.backend}', expected reference|fast_gpu|folded_gpu")
        self.z_mode = 'height_above_top'

    def _detect_rect_macro_ok(self):
        zs = np.asarray(self.sub_apos[:, 2], dtype=np.float64)
        if zs.size == 0:
            return False
        zuniq = []
        tol = 1e-4
        for z in np.sort(zs):
            if (not zuniq) or (abs(z - zuniq[-1]) > tol):
                zuniq.append(z)
        return len(zuniq) <= 3

    def make_world_zs(self, z_range, nz):
        hs = np.linspace(float(z_range[1]), float(z_range[0]), int(nz), dtype=np.float64)
        return self.z_top + hs

    def z_world_to_report(self, z_world):
        return np.asarray(z_world, dtype=np.float64) - self.z_top

    def _eval_reference_component(self, points, component):
        comp = str(component).lower()
        self.scanner.enable_Morse = False
        self.scanner.enable_LJ = False
        self.scanner.enable_Coulomb = False
        self.scanner.enable_HBond = False
        if comp == 'coulomb':
            self.scanner.enable_Coulomb = True
        elif comp in ('nonel', 'lj', 'morse'):
            self.scanner.enable_LJ = True
        elif comp == 'hbond':
            self.scanner.enable_HBond = True
        elif comp == 'total':
            self.scanner.enable_LJ = True
            self.scanner.enable_Coulomb = True
            self.scanner.enable_HBond = False
        else:
            raise ValueError(f"XYZProbeBackend(reference): unsupported component '{component}', expected total|nonel|coulomb|hbond")
        res = self.scanner.evaluate(_identity_transforms(points))
        if comp in ('nonel', 'lj', 'morse'):
            return np.asarray(res['LJ'], dtype=np.float64)
        if comp == 'total':
            return np.asarray(res['total'], dtype=np.float64)
        if comp == 'coulomb':
            return np.asarray(res['Coulomb'], dtype=np.float64)
        if comp == 'hbond':
            return np.asarray(res['HBond'], dtype=np.float64)
        raise ValueError(f"XYZProbeBackend(reference): internal unsupported component '{component}'")

    def _ensure_folded_ready(self, components):
        if self.backend != 'folded_gpu':
            return
        if self._folded_ready:
            return
        comp_keys = []
        for c in components:
            cc = str(c).lower()
            if cc == 'nonel':
                comp_keys += ['pauli', 'london']
            elif cc in ('total', 'pauli', 'london', 'coulomb'):
                comp_keys.append(cc)
        comp_keys = tuple(sorted(set(comp_keys)))
        if not comp_keys:
            comp_keys = ('total',)
        nu, nv = 4, 4
        nz = max(2, min(4, FOLDED_BASIS_MAX // max(1, nu * nv)))
        z_top = float(np.max(self.sub_apos[:, 2]))
        z_range = (z_top + 0.5, z_top + 8.0)
        self.md.fit_folded_surface_basis(self.substrate_xyz, type_map=self.type_map, nPBC=self.nPBC, z_range=z_range, nu=nu, nv=nv, nz=nz, nxy=24, nz_samp=32, r_damp=0.0, alpha_morse=self.probe.alpha, bMacro=True, components=comp_keys)
        self._folded_ready = True

    def eval_components(self, points, components=('total', 'coulomb')):
        comps = [str(c).lower() for c in components]
        pts = np.asarray(points, dtype=np.float32).reshape(-1, 3)
        out = {}
        if self.backend == 'reference':
            need = set(comps)
            for comp in sorted(need):
                if comp == 'nonel':
                    out['nonel'] = self._eval_reference_component(pts, 'nonel')
                else:
                    out[comp] = self._eval_reference_component(pts, comp)
            return out
        T = _identity_transforms(pts)
        if self.backend == 'fast_gpu':
            need = []
            for comp in comps:
                if comp == 'nonel':
                    need += ['pauli', 'london']
                elif comp in ('total', 'pauli', 'london', 'coulomb'):
                    need.append(comp)
                else:
                    raise ValueError(f"XYZProbeBackend(fast_gpu): unsupported component '{comp}', expected total|nonel|pauli|london|coulomb")
            res = self.md.eval_rigid_getSurfMorse_components(T, chunk_size=self.chunk_size, components=tuple(sorted(set(need))))
            for k in sorted(set(need)):
                out[k] = np.asarray(res[k], dtype=np.float64)
            if 'nonel' in comps:
                out['nonel'] = out.get('pauli', 0.0) + out.get('london', 0.0)
            return out
        if self.backend == 'folded_gpu':
            self._ensure_folded_ready(comps)
            uniq = []
            for comp in comps:
                if comp == 'nonel':
                    uniq += ['pauli', 'london']
                elif comp in ('total', 'pauli', 'london', 'coulomb'):
                    uniq.append(comp)
                else:
                    raise ValueError(f"XYZProbeBackend(folded_gpu): unsupported component '{comp}', expected total|nonel|pauli|london|coulomb")
            for comp in sorted(set(uniq)):
                out[comp] = np.asarray(self.md.eval_rigid_getSurfFolded(T, chunk_size=self.chunk_size, component=comp)['total'], dtype=np.float64)
            if 'nonel' in comps:
                out['nonel'] = out.get('pauli', 0.0) + out.get('london', 0.0)
            return out
        raise ValueError(f"XYZProbeBackend: invalid backend '{self.backend}'")

    def eval_surface_map_gpu(self, x_range, y_range, z_range, nx=64, ny=64, nz=80, selector='total', color_component='coulomb', mode='threshold', threshold=0.0):
        if self.backend == 'reference':
            raise ValueError('XYZProbeBackend.eval_surface_map_gpu(): reference backend is debug-only and has no production GPU surface path')
        sel = str(selector).lower()
        csel = str(color_component).lower()
        if self.backend == 'fast_gpu':
            xs = np.linspace(float(x_range[0]), float(x_range[1]), int(nx), dtype=np.float32)
            ys = np.linspace(float(y_range[0]), float(y_range[1]), int(ny), dtype=np.float32)
            zs = self.make_world_zs(z_range, nz).astype(np.float32)
            XX, YY, ZZ = np.meshgrid(xs, ys, zs, indexing='ij')
            pts = np.stack([XX, YY, ZZ], axis=-1).reshape(-1, 3)
            comps = [sel] if sel == csel else [sel, csel]
            vals = self.eval_components(pts, components=tuple(comps))
            vsel = np.asarray(vals[sel], dtype=np.float64).reshape(len(xs), len(ys), len(zs))
            if mode == 'threshold':
                Z = np.full((len(xs), len(ys)), np.nan, dtype=np.float64)
                ok = np.zeros((len(xs), len(ys)), dtype=bool)
                for iy in range(len(ys)):
                    zh, okm = _refine_threshold_rows(zs.astype(np.float64), vsel[:, iy, :], threshold)
                    Z[:, iy] = zh
                    ok[:, iy] = okm
            elif mode == 'first_minimum':
                Z = np.full((len(xs), len(ys)), np.nan, dtype=np.float64)
                ok = np.zeros((len(xs), len(ys)), dtype=bool)
                zr = zs[::-1].astype(np.float64)
                for iy in range(len(ys)):
                    zh, okm = _refine_first_min_rows(zr, vsel[:, iy, ::-1])
                    Z[:, iy] = zh
                    ok[:, iy] = okm
            else:
                raise ValueError(f"XYZProbeBackend.eval_surface_map_gpu(): unsupported mode '{mode}'")
            points = np.zeros((len(xs), len(ys), 3), dtype=np.float64)
            points[:, :, 0] = xs[:, None]
            points[:, :, 1] = ys[None, :]
            points[:, :, 2] = Z
            good = ok & np.isfinite(Z)
            color = np.full((len(xs), len(ys)), np.nan, dtype=np.float64)
            if np.any(good):
                pgood = points[good]
                color[good] = np.asarray(self.eval_components(pgood, components=(csel,))[csel], dtype=np.float64)
            return {'points_world': points, 'ok_mask': good, 'z_report': self.z_world_to_report(Z), 'color': color}
        if self.backend == 'folded_gpu':
            raise ValueError('XYZProbeBackend.eval_surface_map_gpu(): folded_gpu production surface map is not implemented for the new 2D GPU iso path')
        raise ValueError(f"XYZProbeBackend.eval_surface_map_gpu(): invalid backend '{self.backend}'")


def build_backend(source='xyz', substrate_xyz=None, gridff_path=None, probe=None, xyz_type_map=None, backend='reference', xyz_path=None, nPBC=(4, 4, 0), chunk_size=8192, z0=None):
    if probe is None:
        raise ValueError('build_backend(): probe must be provided')
    src = str(source).lower()
    if src == 'gridff':
        if gridff_path is None:
            raise ValueError('build_backend(source=gridff): gridff_path is required')
        return GridFFProbeBackend(gridff_path=gridff_path, probe=probe, xyz_path=xyz_path, z0=z0, chunk_size=chunk_size, use_gpu=True, type_map=xyz_type_map)
    if src == 'xyz':
        if substrate_xyz is None:
            raise ValueError('build_backend(source=xyz): substrate_xyz is required')
        return XYZProbeBackend(substrate_xyz=substrate_xyz, probe=probe, type_map=xyz_type_map, backend=backend, nPBC=nPBC, chunk_size=chunk_size)
    raise ValueError(f"build_backend(): unsupported source '{source}', expected gridff|xyz")


def refine_threshold_crossing(zs, vals, threshold):
    z = np.asarray(zs, dtype=np.float64)
    v = np.asarray(vals, dtype=np.float64) - float(threshold)
    for i in range(len(z) - 1):
        if (v[i] <= 0.0 and v[i + 1] >= 0.0) or (v[i] >= 0.0 and v[i + 1] <= 0.0):
            dv = v[i + 1] - v[i]
            t = 0.5 if abs(dv) < 1e-16 else (-v[i] / dv)
            t = min(1.0, max(0.0, t))
            return float(z[i] + t * (z[i + 1] - z[i])), i
    return np.nan, -1


def refine_first_minimum(zs, vals):
    z = np.asarray(zs, dtype=np.float64)
    v = np.asarray(vals, dtype=np.float64)
    for i in range(1, len(z) - 1):
        if (v[i] <= v[i - 1]) and (v[i] <= v[i + 1]) and ((v[i] < v[i - 1]) or (v[i] < v[i + 1])):
            z0, z1, z2 = z[i - 1], z[i], z[i + 1]
            v0, v1, v2 = v[i - 1], v[i], v[i + 1]
            den = (z0 - z1) * (z0 - z2) * (z1 - z2)
            if abs(den) < 1e-16:
                return float(z1), i
            A = (z2 * (v1 - v0) + z1 * (v0 - v2) + z0 * (v2 - v1)) / den
            B = (z2 * z2 * (v0 - v1) + z1 * z1 * (v2 - v0) + z0 * z0 * (v1 - v2)) / den
            if abs(A) < 1e-16:
                return float(z1), i
            zm = -B / (2.0 * A)
            if zm < min(z0, z2) or zm > max(z0, z2):
                zm = z1
            return float(zm), i
    return np.nan, -1


def _refine_threshold_rows(zs, vals, thr):
    s = np.asarray(vals, dtype=np.float64) - float(thr)
    nrow, nz_ = s.shape
    idx = np.full(nrow, -1, dtype=np.int64)
    okm = np.zeros(nrow, dtype=bool)
    for i in range(nz_ - 1):
        hit = (~okm) & np.isfinite(s[:, i]) & np.isfinite(s[:, i + 1]) & (((s[:, i] <= 0.0) & (s[:, i + 1] >= 0.0)) | ((s[:, i] >= 0.0) & (s[:, i + 1] <= 0.0)))
        idx[hit] = i
        okm[hit] = True
    zh = np.full(nrow, np.nan, dtype=np.float64)
    for m in np.where(okm)[0]:
        i = idx[m]
        dv = s[m, i + 1] - s[m, i]
        t = 0.5 if abs(dv) < 1e-16 else (-s[m, i] / dv)
        t = min(1.0, max(0.0, t))
        zh[m] = float(zs[i] + t * (zs[i + 1] - zs[i]))
    return zh, okm


def _refine_first_min_rows(zs, vals):
    z = np.asarray(zs, dtype=np.float64)
    vals = np.asarray(vals, dtype=np.float64)
    nrow, nz_ = vals.shape
    zh = np.full(nrow, np.nan, dtype=np.float64)
    okm = np.zeros(nrow, dtype=bool)
    if nz_ < 3:
        return zh, okm
    vm = vals[:, 1:-1]
    v0 = vals[:, :-2]
    v2 = vals[:, 2:]
    cand = np.isfinite(vm) & np.isfinite(v0) & np.isfinite(v2) & (vm <= v0) & (vm <= v2) & ((vm < v0) | (vm < v2))
    anyc = np.any(cand, axis=1)
    if not np.any(anyc):
        return zh, okm
    idx = np.argmax(cand, axis=1) + 1
    for m in np.where(anyc)[0]:
        i = int(idx[m])
        z0, z1, z2 = z[i - 1], z[i], z[i + 1]
        v0_, v1_, v2_ = vals[m, i - 1], vals[m, i], vals[m, i + 1]
        den = (z0 - z1) * (z0 - z2) * (z1 - z2)
        if abs(den) < 1e-16:
            zh[m] = float(z1)
        else:
            A = (z2 * (v1_ - v0_) + z1 * (v0_ - v2_) + z0 * (v2_ - v1_)) / den
            B = (z2 * z2 * (v0_ - v1_) + z1 * z1 * (v2_ - v0_) + z0 * z0 * (v1_ - v2_)) / den
            zm = float(z1) if abs(A) < 1e-16 else float(-B / (2.0 * A))
            if (zm < min(z0, z2)) or (zm > max(z0, z2)):
                zm = float(z1)
            zh[m] = zm
        okm[m] = True
    return zh, okm


class SurfaceMapResult:
    def __init__(self, xs, ys, zs, color, selector, color_component, ok_mask, diagnostics, points, zs_world=None):
        self.xs = np.asarray(xs, dtype=np.float64)
        self.ys = np.asarray(ys, dtype=np.float64)
        self.zs = np.asarray(zs, dtype=np.float64)
        self.zs_world = np.asarray(self.zs if zs_world is None else zs_world, dtype=np.float64)
        self.color = np.asarray(color, dtype=np.float64)
        self.selector = str(selector)
        self.color_component = str(color_component)
        self.ok_mask = np.asarray(ok_mask, dtype=bool)
        self.diagnostics = diagnostics
        self.points = np.asarray(points, dtype=np.float64)


def make_surface_backend(source='xyz', substrate_xyz=None, gridff_path=None, probe=None, xyz_type_map=None, backend='reference', xyz_path=None, nPBC=(4, 4, 0), chunk_size=8192, z0=None):
    return build_backend(source=source, substrate_xyz=substrate_xyz, gridff_path=gridff_path, probe=probe, xyz_type_map=xyz_type_map, backend=backend, xyz_path=xyz_path, nPBC=nPBC, chunk_size=chunk_size, z0=z0)


def build_surface_map(source='xyz', substrate_xyz=None, gridff_path=None, probe=None, xyz_type_map=None, backend='reference', xyz_path=None, nPBC=(4, 4, 0), chunk_size=8192, z0=None, x_range=(0.0, 4.0), y_range=(0.0, 4.0), z_range=(0.5, 8.0), nx=64, ny=64, nz=80, selector='total', color_component='coulomb', mode='threshold', threshold=0.0, bFailOnMiss=False):
    bk = make_surface_backend(source=source, substrate_xyz=substrate_xyz, gridff_path=gridff_path, probe=probe, xyz_type_map=xyz_type_map, backend=backend, xyz_path=xyz_path, nPBC=nPBC, chunk_size=chunk_size, z0=z0)
    return sample_surface_map(bk, x_range=x_range, y_range=y_range, z_range=z_range, nx=nx, ny=ny, nz=nz, selector=selector, color_component=color_component, mode=mode, threshold=threshold, bFailOnMiss=bFailOnMiss)



def sample_surface_map(backend, x_range, y_range, z_range, nx=64, ny=64, nz=80, selector='total', color_component='coulomb', mode='threshold', threshold=0.0, bFailOnMiss=False):
    if hasattr(backend, 'eval_surface_map_gpu') and (getattr(backend, 'backend', None) != 'reference'):
        res = backend.eval_surface_map_gpu(x_range=x_range, y_range=y_range, z_range=z_range, nx=nx, ny=ny, nz=nz, selector=selector, color_component=color_component, mode=mode, threshold=threshold)
        ok = np.asarray(res['ok_mask'], dtype=bool)
        if bFailOnMiss and np.any(~ok):
            bad = np.argwhere(~ok)
            ib = bad[0]
            raise RuntimeError(f"sample_surface_map(): failed at {len(bad)}/{ok.size} xy points; first_fail=({int(ib[0])},{int(ib[1])})")
        pts = np.asarray(res['points_world'], dtype=np.float64)
        zs_world = np.asarray(pts[:, :, 2], dtype=np.float64)
        zs_report = np.asarray(res['z_report'], dtype=np.float64)
        color = np.asarray(res['color'], dtype=np.float64)
        xs = np.linspace(float(x_range[0]), float(x_range[1]), int(nx), dtype=np.float64)
        ys = np.linspace(float(y_range[0]), float(y_range[1]), int(ny), dtype=np.float64)
        fail = [(int(ix), int(iy), float(xs[ix]), float(ys[iy])) for ix, iy in np.argwhere(~ok)]
        return SurfaceMapResult(xs=xs, ys=ys, zs=zs_report, zs_world=zs_world, color=color, selector=str(selector).lower(), color_component=str(color_component).lower(), ok_mask=ok, diagnostics={'fail_points': fail, 'mode': mode, 'threshold': float(threshold), 'z_range': tuple(map(float, z_range)), 'z_mode': getattr(backend, 'z_mode', 'world_absolute'), 'n_fail': len(fail), 'n_total': int(ok.size), 'path': 'gpu_2d_iso'}, points=pts)

    # DEBUG / REFERENCE ONLY: deprecated host-driven fallback kept for validation, not production.
    xs = np.linspace(float(x_range[0]), float(x_range[1]), int(nx), dtype=np.float64)
    ys = np.linspace(float(y_range[0]), float(y_range[1]), int(ny), dtype=np.float64)
    z_desc_world = backend.make_world_zs(z_range, nz) if hasattr(backend, 'make_world_zs') else np.linspace(float(z_range[1]), float(z_range[0]), int(nz), dtype=np.float64)
    sel = str(selector).lower()
    csel = str(color_component).lower()
    need = [sel]
    if csel not in need:
        need.append(csel)
    Z = np.full((len(xs), len(ys)), np.nan, dtype=np.float64)
    C = np.full((len(xs), len(ys)), np.nan, dtype=np.float64)
    ok = np.zeros((len(xs), len(ys)), dtype=bool)
    fail = []
    eval_points = np.zeros((len(xs), len(ys), 3), dtype=np.float64)
    def _refine_threshold_rows(zs, vals, thr):
        s = vals - float(thr)
        nrow, nz_ = s.shape
        idx = np.full(nrow, -1, dtype=np.int64)
        okm = np.zeros(nrow, dtype=bool)
        for i in range(nz_ - 1):
            hit = (~okm) & np.isfinite(s[:, i]) & np.isfinite(s[:, i + 1]) & (((s[:, i] <= 0.0) & (s[:, i + 1] >= 0.0)) | ((s[:, i] >= 0.0) & (s[:, i + 1] <= 0.0)))
            idx[hit] = i
            okm[hit] = True
        zh = np.full(nrow, np.nan, dtype=np.float64)
        for m in np.where(okm)[0]:
            i = idx[m]
            dv = s[m, i + 1] - s[m, i]
            t = 0.5 if abs(dv) < 1e-16 else (-s[m, i] / dv)
            t = min(1.0, max(0.0, t))
            zh[m] = float(zs[i] + t * (zs[i + 1] - zs[i]))
        return zh, okm

    def _refine_first_min_rows(zs, vals):
        nrow, nz_ = vals.shape
        zh = np.full(nrow, np.nan, dtype=np.float64)
        okm = np.zeros(nrow, dtype=bool)
        if nz_ < 3:
            return zh, okm
        vm = vals[:, 1:-1]
        v0 = vals[:, :-2]
        v2 = vals[:, 2:]
        cand = np.isfinite(vm) & np.isfinite(v0) & np.isfinite(v2) & (vm <= v0) & (vm <= v2) & ((vm < v0) | (vm < v2))
        anyc = np.any(cand, axis=1)
        if not np.any(anyc):
            return zh, okm
        idx = np.argmax(cand, axis=1) + 1
        for m in np.where(anyc)[0]:
            i = int(idx[m])
            z0, z1, z2 = zs[i - 1], zs[i], zs[i + 1]
            v0_, v1_, v2_ = vals[m, i - 1], vals[m, i], vals[m, i + 1]
            den = (z0 - z1) * (z0 - z2) * (z1 - z2)
            if abs(den) < 1e-16:
                zh[m] = float(z1)
            else:
                A = (z2 * (v1_ - v0_) + z1 * (v0_ - v2_) + z0 * (v2_ - v1_)) / den
                B = (z2 * z2 * (v0_ - v1_) + z1 * z1 * (v2_ - v0_) + z0 * z0 * (v1_ - v2_)) / den
                zm = float(z1) if abs(A) < 1e-16 else float(-B / (2.0 * A))
                if (zm < min(z0, z2)) or (zm > max(z0, z2)):
                    zm = float(z1)
                zh[m] = zm
            okm[m] = True
        return zh, okm

    for ix, x in enumerate(xs):
        pts = np.zeros((len(ys) * len(z_desc_world), 3), dtype=np.float64)
        pts[:, 0] = x
        pts[:, 1] = np.repeat(ys, len(z_desc_world))
        pts[:, 2] = np.tile(z_desc_world, len(ys))
        vals = backend.eval_components(pts, components=need)
        vsel = np.asarray(vals[sel], dtype=np.float64).reshape(len(ys), len(z_desc_world))
        if mode == 'threshold':
            z_hits, ok_row = _refine_threshold_rows(z_desc_world, vsel, threshold)
        elif mode == 'first_minimum':
            z_hits, ok_row = _refine_first_min_rows(z_desc_world[::-1], vsel[:, ::-1])
        else:
            raise ValueError(f"sample_surface_map(): unsupported mode '{mode}', expected threshold|first_minimum")
        good = np.where(ok_row & np.isfinite(z_hits))[0]
        bad = np.where(~(ok_row & np.isfinite(z_hits)))[0]
        for iy in bad:
            fail.append((ix, int(iy), float(x), float(ys[iy])))
        if len(good) == 0:
            continue
        pgood = np.zeros((len(good), 3), dtype=np.float64)
        pgood[:, 0] = x
        pgood[:, 1] = ys[good]
        pgood[:, 2] = z_hits[good]
        cvals = np.asarray(backend.eval_components(pgood, components=(csel,))[csel], dtype=np.float64)
        z_rep = backend.z_world_to_report(z_hits[good]) if hasattr(backend, 'z_world_to_report') else z_hits[good]
        Z[ix, good] = np.asarray(z_rep, dtype=np.float64)
        C[ix, good] = cvals
        ok[ix, good] = True
        eval_points[ix, good, :] = pgood
    if bFailOnMiss and len(fail) > 0:
        raise RuntimeError(f"sample_surface_map(): failed at {len(fail)}/{len(xs)*len(ys)} xy points; first_fail={fail[0]}")
    return SurfaceMapResult(xs=xs, ys=ys, zs=Z, zs_world=eval_points[:, :, 2], color=C, selector=sel, color_component=csel, ok_mask=ok, diagnostics={'fail_points': fail, 'mode': mode, 'threshold': float(threshold), 'z_range': tuple(map(float, z_range)), 'z_mode': getattr(backend, 'z_mode', 'world_absolute'), 'n_fail': len(fail), 'n_total': int(len(xs) * len(ys)), 'path': 'debug_host_fallback'}, points=eval_points)


def surface_normals(xs, ys, zs):
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    z = np.asarray(zs, dtype=np.float64)
    if z.shape != (len(x), len(y)):
        raise ValueError(f"surface_normals(): z.shape={z.shape} expected ({len(x)},{len(y)})")
    dzdx = np.gradient(z, x, axis=0, edge_order=1)
    dzdy = np.gradient(z, y, axis=1, edge_order=1)
    n = np.zeros(z.shape + (3,), dtype=np.float64)
    n[..., 0] = -dzdx
    n[..., 1] = -dzdy
    n[..., 2] = 1.0
    nn = np.linalg.norm(n, axis=2)
    nn[nn <= 1e-16] = 1.0
    n /= nn[:, :, None]
    return n
