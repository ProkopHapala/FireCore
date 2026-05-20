import os
import numpy as np
import pyBall.atomicUtils as au
from pyBall.OCL import AFM as afm
from pyBall.OCL import AFM_utils as afm_utils
from pyBall.config_utils import get_config, get_path, get_dftb_basis_path

class ModularAFMPipeline:
    """
    Decoupled, stage-based modular pipeline for AFM and STM simulations.
    Saves intermediate results to disk allowing fast, independent stage execution.

    Supports multiple quantum chemistry backends:
    - 'dftb': DFTB+ with Slater-type orbitals (GPU-accelerated projection)
    - 'pyscf': pySCF with Gaussian-type orbitals (CPU-based evaluation)
    """
    def __init__(self, xyz_file, output_dir, basis='mio-1-1', slako_prefix='mio-1-1',
                 work_dir=None, step=0.1, margin=4.0, z_extra=6.0,
                 scan_range=3.0, scan_step=0.1, height_range=(2.8, 3.6), height_step=0.1,
                 co_tip_dir=None,
                 atomPos=None, enames=None,  # Optional: inject geometry directly instead of xyz_file
                 backend='dftb', pyscf_params=None):  # Backend selection and parameters
        self.xyz_file = xyz_file
        self.output_dir = output_dir
        self._injected_atomPos = atomPos   # If provided, skip loading xyz_file
        self._injected_enames = enames
        self.basis = basis
        self.slako_prefix = slako_prefix
        self.co_tip_dir = co_tip_dir
        self.backend = backend.lower()
        self.pyscf_params = pyscf_params or {'method': 'RHF', 'basis': 'sto-3g', 'xc': None}
        
        self.work_dir = work_dir or os.path.join(output_dir, 'dftb_work')
        self.step = step
        self.margin = margin
        self.z_extra = z_extra
        self.scan_range = scan_range
        self.scan_step = scan_step
        self.height_range = height_range
        self.height_step = height_step
        
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.work_dir, exist_ok=True)
        
        # Cache file paths
        self.cache_stage1 = os.path.join(self.output_dir, 'cache_stage1_scf.npz')
        self.cache_stage2 = os.path.join(self.output_dir, 'cache_stage2_grids.npz')
        self.cache_stage3 = os.path.join(self.output_dir, 'cache_stage3_potentials.npz')
        self.cache_stage4 = os.path.join(self.output_dir, 'cache_stage4_relax.npz')
        
        # Grid parameters
        self.origin = None
        self.ngrid = None
        self.grid_spec = None
        self.scan_xs = None
        self.scan_ys = None
        self.heights = None
        
        # DFTB structures (only used for backend='dftb')
        self.atomPos = None
        self.atomTypes = None
        self.enames = None
        self.projector = None
        self.atoms_dict = None
        self.norb_per_atom = None
        self.orb_offsets = None

        # pySCF structures (only used for backend='pyscf')
        self._pyscf_data = None  # Cache for mol, mf, dm from pySCF
        
        # Load molecule and scan grid parameters
        self._init_geometry_and_grids()

    def _init_geometry_and_grids(self):
        """Load molecular structure and define grid parameters."""
        ELEM_Z = {'H':1,'C':6,'N':7,'O':8,'P':15,'S':16,'Br':35,'I':53}
        inv_z = {v:k for k,v in ELEM_Z.items()}
        
        if self._injected_atomPos is not None and self._injected_enames is not None:
            print(f"\n[ModularPipeline] Using injected geometry ({len(self._injected_atomPos)} atoms)")
            self.atomPos   = np.array(self._injected_atomPos, dtype=np.float64)
            self.enames    = list(self._injected_enames)
            self.atomTypes = np.array([ELEM_Z.get(e, 6) for e in self.enames], dtype=np.int32)
        else:
            print(f"\n[ModularPipeline] Loading molecule from {self.xyz_file}")
            pos, _, names, _, _ = au.load_xyz(self.xyz_file)
            self.atomPos  = np.array(pos, dtype=np.float64)
            self.atomTypes = np.array([ELEM_Z.get(e, 6) for e in names], dtype=np.int32)
            self.enames = [inv_z.get(int(z), 'C') for z in self.atomTypes]
        print(f"  {len(self.atomPos)} atoms loaded.")
        
        # Scan grid coordinates
        x_min = self.atomPos[:,0].min() - self.scan_range
        x_max = self.atomPos[:,0].max() + self.scan_range
        y_min = self.atomPos[:,1].min() - self.scan_range
        y_max = self.atomPos[:,1].max() + self.scan_range
        scan_points_x = int(np.ceil((x_max - x_min) / self.scan_step))
        scan_points_y = int(np.ceil((y_max - y_min) / self.scan_step))
        self.scan_xs = np.linspace(x_min, x_max, scan_points_x)
        self.scan_ys = np.linspace(y_min, y_max, scan_points_y)
        self.heights  = np.arange(self.height_range[0], self.height_range[1], self.height_step)
        
        # Setup projector / backend-specific initialization
        if self.backend == 'dftb':
            self._init_dftb_backend()
        elif self.backend == 'pyscf':
            self._init_pyscf_backend()
        else:
            raise ValueError(f"Unknown backend: {backend}. Use 'dftb' or 'pyscf'.")

    def _init_dftb_backend(self):
        """Initialize DFTB backend: setup GridProjector with STO basis."""
        from pyBall import dftb_utils as du
        from pyBall.DFTB.DFTBplusParser import parse_wfc_hsd, convert_wfc_to_species_list_ang
        from pyBall.DFTB import Grid_dftb as dg

        basis_name = self.basis
        if basis_name == 'mio-1-1':
            self.slako_prefix = du.SK_PATHS.get('mio-1-1', self.slako_prefix)
        elif basis_name == '3ob-3-1':
            self.slako_prefix = du.SK_PATHS.get('3ob-3-1', self.slako_prefix)

        basis_name = self.slako_prefix.rstrip('/').split('/')[-1] if '/' in self.slako_prefix else self.slako_prefix
        if not basis_name:
            basis_name = '3ob-3-1'

        # Use config system to find basis file
        basis_hsd_path = get_dftb_basis_path(basis_name)
        if basis_hsd_path is None:
            # Fallback to old hardcoded path
            _ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
            basis_hsd_path = os.path.join(_ROOT, 'pyBall', 'DFTB', 'data', f'wfc.{basis_name}.hsd')
        
        if os.path.exists(basis_hsd_path):
            basis_data = parse_wfc_hsd(basis_hsd_path)
            basis_ang = convert_wfc_to_species_list_ang(basis_data, resolution_bohr=0.04)
            self.norb_per_atom, self.orb_offsets, max_l = afm_utils.build_orbital_layout(basis_data, self.enames)
            max_shells = 3 if max_l >= 2 else 2

            coords_bohr = self.atomPos * 1.8897259886
            species_per_atom = list(range(len(self.enames)))
            dftb_data = {
                'coords_bohr': coords_bohr,
                'species_per_atom': species_per_atom,
                'species_names': self.enames
            }
            self.projector, self.atoms_dict = dg.setup_gridprojector_from_dftb(dftb_data, basis_ang, verbosity=0, max_shells=max_shells)
            print(f"[ModularPipeline] DFTB backend initialized with {len(self.enames)} atoms")
        else:
            print(f"[ModularPipeline] WARNING: Basis file not found: {basis_hsd_path}")

    def _init_pyscf_backend(self):
        """Initialize pySCF backend: no projector setup (CPU-based evaluation)."""
        print(f"[ModularPipeline] pySCF backend initialized (method={self.pyscf_params.get('method', 'RHF')}, basis={self.pyscf_params.get('basis', 'sto-3g')})")
        # No projector needed for pySCF - density is computed directly on grid
        self.projector = None
        self.atoms_dict = None
        self.norb_per_atom = None
        self.orb_offsets = None

    def stage1_scf(self, force_recompute=False):
        """Stage 1: SCF computation (DFTB or pySCF depending on backend)."""
        if not force_recompute and os.path.exists(self.cache_stage1):
            print(f"\n[ModularPipeline] Loading Stage 1 (SCF) from cache...")
            data = np.load(self.cache_stage1, allow_pickle=True)
            if self.backend == 'dftb':
                return data['dm_dense'], data['eigvecs'], data['eigvals']
            else:  # pySCF
                self._pyscf_data = {k: data[k] for k in data.keys() if k.startswith('mol_') or k in ['dm', 'eigvecs', 'eigvals']}
                return data.get('dm'), data['eigvecs'], data['eigvals']

        print(f"\n[ModularPipeline] Running Stage 1 (SCF) with backend='{self.backend}'...")

        if self.backend == 'dftb':
            return self._stage1_scf_dftb()
        else:  # pySCF
            return self._stage1_scf_pyscf()

    def _stage1_scf_dftb(self):
        """DFTB backend: run DFTBcore SCF and extract density matrix."""
        from pyBall.DFTB.DFTBcore import DFTBcore
        from pyBall.DFTB.DFTBplusParser import parse_wfc_hsd
        from pyBall import dftb_utils as du
        import shutil

        basis_name = self.slako_prefix.rstrip('/').split('/')[-1] if '/' in self.slako_prefix else self.slako_prefix
        if not basis_name:
            basis_name = '3ob-3-1'
        
        # Use config system to find basis file
        basis_hsd_path = get_dftb_basis_path(basis_name)
        if basis_hsd_path is None:
            # Fallback to old hardcoded path
            _ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
            basis_hsd_path = os.path.join(_ROOT, 'pyBall', 'DFTB', 'data', f'wfc.{basis_name}.hsd')

        sk_dir = du.SK_PATHS.get(basis_name, os.path.join(os.environ.get('DFTB_SK_PATH', ''), basis_name))

        # Set up DFTBcore directory and input
        xyz_path = os.path.join(self.work_dir, 'geom.xyz')
        hsd_path = os.path.join(self.work_dir, 'dftb_in.hsd')
        au.save_xyz(xyz_path, self.enames, self.atomPos)

        basis_data = parse_wfc_hsd(basis_hsd_path)
        species = sorted(set(self.enames))
        max_am_map = {0: 's', 1: 'p', 2: 'd'}
        max_ang_lines = []
        for elem in species:
            elem_data = basis_data[elem]
            max_l = max(orb['AngularMomentum'] for orb in elem_data['orbitals'])
            max_ang_lines.append(f'    {elem} = "{max_am_map[max_l]}"')
        max_ang_str = '\n'.join(max_ang_lines)

        with open(hsd_path, 'w') as f:
            f.write(f'''Geometry = xyzFormat {{
  <<< "geom.xyz"
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SCCTolerance = 1e-7
  MaxSCCIterations = 200
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{sk_dir}/"
    Separator = "-"
    Suffix = ".skf"
    LowerCaseTypeName = No
  }}
  MaxAngularMomentum = {{
{max_ang_str}
  }}
}}
''')
        for i, elem1 in enumerate(species):
            for elem2 in species[i:]:
                for sk_file in [f"{elem1}-{elem2}.skf", f"{elem2}-{elem1}.skf"]:
                    src = os.path.join(sk_dir, sk_file)
                    if os.path.exists(src):
                        shutil.copy(src, self.work_dir)

        old_cwd = os.getcwd()
        try:
            os.chdir(self.work_dir)
            dftb = DFTBcore()
            dftb.init('dftb_in.hsd')
            dftb.enable_matrix_collection(dm=True, h=False, s=False)
            dftb.run_scf()
            dm_dense = dftb.get_dm_dense()
            eigvecs, eigvals = dftb.get_eigvecs_dense()
            dftb.finalize()
        finally:
            os.chdir(old_cwd)

        np.savez_compressed(self.cache_stage1, dm_dense=dm_dense, eigvecs=eigvecs, eigvals=eigvals)
        print(f"  Stage 1 (DFTB) complete and cached.")
        return dm_dense, eigvecs, eigvals

    def _stage1_scf_pyscf(self):
        """pySCF backend: run pySCF SCF and compute density on grid (combined Stage 1+2)."""
        # For pySCF, we combine Stage 1 and Stage 2 (density is computed on grid directly)
        # This avoids caching dm_dense which is specific to STO basis representation
        result = afm_utils.get_density_from_pyscf(
            self.atomPos, self.atomTypes,
            step=self.step, margin=self.margin, z_extra=self.z_extra,
            **self.pyscf_params
        )

        # Cache the result for Stage 2 to load
        # Note: pySCF doesn't have dm_dense in the same sense as DFTB
        np.savez_compressed(self.cache_stage1,
                           dm=result['dm'],  # Density matrix in AO basis
                           eigvecs=result['eigvecs'],
                           eigvals=result['eigvals'],
                           rho_scf=result['rho_scf'],  # Pre-computed on grid
                           rho_na=result['rho_na'],
                           rho_diff=result['rho_diff'],
                           origin=result['origin'],
                           ngrid=result['ngrid'],
                           grid_spec=result['grid_spec'])

        self._pyscf_data = {'mol': result['mol'], 'mf': result['mf'], 'dm': result['dm']}

        print(f"  Stage 1 (pySCF) complete and cached.")
        # Return None for dm_dense since pySCF doesn't use dense STO projection
        return None, result['eigvecs'], result['eigvals']

    def stage2_project(self, dm_dense, force_recompute=False):
        """Stage 2: Grid density projection (SCF, Neutral Atom, and Diff).

        For DFTB backend: uses GPU projector with dm_dense.
        For pySCF backend: densities already computed in Stage 1, just loads from cache.
        """
        # For pySCF backend, densities were computed in Stage 1
        if self.backend == 'pyscf':
            if not force_recompute and os.path.exists(self.cache_stage1):
                print(f"\n[ModularPipeline] Loading Stage 2 (pySCF densities) from Stage 1 cache...")
                data = np.load(self.cache_stage1, allow_pickle=True)
                self.origin = data['origin']
                self.ngrid = data['ngrid']
                rho_scf, rho_na, rho_diff = data['rho_scf'], data['rho_na'], data['rho_diff']
                # Reconstruct grid_spec
                if 'grid_spec' in data:
                    self.grid_spec = data['grid_spec'].item() if isinstance(data['grid_spec'], np.ndarray) else data['grid_spec']
                else:
                    self.grid_spec = {
                        'origin': self.origin,
                        'dA': np.array([self.step, 0.0, 0.0], dtype=np.float32),
                        'dB': np.array([0.0, self.step, 0.0], dtype=np.float32),
                        'dC': np.array([0.0, 0.0, self.step], dtype=np.float32),
                        'ngrid': self.ngrid,
                    }
                z_profile = rho_scf.sum(axis=(0, 1))
                iz_max = int(np.argmax(z_profile))
                print(f"  [Stage2 pySCF] rho_scf: shape={rho_scf.shape} range=[{rho_scf.min():.4e},{rho_scf.max():.4e}] sum={rho_scf.sum():.4e}")
                print(f"  [Stage2 pySCF] density z-peak at iz={iz_max}, z={float(self.origin[2]) + iz_max*self.step:.3f} A")
                return rho_scf, rho_na, rho_diff
            elif force_recompute:
                # If force_recompute, we need to re-run Stage 1
                print(f"\n[ModularPipeline] force_recompute=True for pySCF, re-running Stage 1...")
                _, _, _ = self.stage1_scf(force_recompute=True)
                data = np.load(self.cache_stage1, allow_pickle=True)
                self.origin = data['origin']
                self.ngrid = data['ngrid']
                rho_scf, rho_na, rho_diff = data['rho_scf'], data['rho_na'], data['rho_diff']
                if 'grid_spec' in data:
                    self.grid_spec = data['grid_spec'].item() if isinstance(data['grid_spec'], np.ndarray) else data['grid_spec']
                else:
                    self.grid_spec = {
                        'origin': self.origin,
                        'dA': np.array([self.step, 0.0, 0.0], dtype=np.float32),
                        'dB': np.array([0.0, self.step, 0.0], dtype=np.float32),
                        'dC': np.array([0.0, 0.0, self.step], dtype=np.float32),
                        'ngrid': self.ngrid,
                    }
                return rho_scf, rho_na, rho_diff
            else:
                raise RuntimeError("pySCF Stage 1 must be run before Stage 2")

        # DFTB backend: standard GPU projection
        if not force_recompute and os.path.exists(self.cache_stage2):
            print(f"\n[ModularPipeline] Loading Stage 2 (density grids) from cache...")
            data = np.load(self.cache_stage2)
            self.origin = data['origin']
            self.ngrid = data['ngrid']
            rho_scf, rho_na, rho_diff = data['rho_scf'], data['rho_na'], data['rho_diff']
            # Reconstruct grid_spec
            self.grid_spec = {
                'origin': self.origin,
                'dA': np.array([self.step, 0.0, 0.0], dtype=np.float32),
                'dB': np.array([0.0, self.step, 0.0], dtype=np.float32),
                'dC': np.array([0.0, 0.0, self.step], dtype=np.float32),
                'ngrid': self.ngrid,
            }
            z_profile = rho_scf.sum(axis=(0, 1))
            iz_max = int(np.argmax(z_profile))
            print(f"  [Stage2 cache] rho_scf: shape={rho_scf.shape} range=[{rho_scf.min():.4e},{rho_scf.max():.4e}] sum={rho_scf.sum():.4e}")
            print(f"  [Stage2 cache] density z-peak at iz={iz_max}, z={float(self.origin[2]) + iz_max*self.step:.3f} A")
            return rho_scf, rho_na, rho_diff

        print(f"\n[ModularPipeline] Projecting Stage 2 (density grids)...")
        from pyBall.DFTB.DFTBplusParser import parse_wfc_hsd, convert_wfc_to_species_list_ang
        from pyBall.DFTB import Grid_dftb as dg

        # Grid parameters setup
        padding = self.margin + self.z_extra
        x_min, x_max = self.atomPos[:,0].min() - self.margin, self.atomPos[:,0].max() + self.margin
        y_min, y_max = self.atomPos[:,1].min() - self.margin, self.atomPos[:,1].max() + self.margin
        z_min, z_max = self.atomPos[:,2].min() - self.margin, self.atomPos[:,2].max() + padding

        origin = np.array([x_min, y_min, z_min], dtype=np.float32)
        ngrid = np.ceil(np.array([x_max - x_min, y_max - y_min, z_max - z_min]) / self.step).astype(np.int32)
        # Round to nearest multiple of 8
        ngrid = ((ngrid + 7) // 8) * 8

        self.origin = origin
        self.ngrid = ngrid
        self.grid_spec = {
            'origin': origin,
            'dA': np.array([self.step, 0.0, 0.0], dtype=np.float32),
            'dB': np.array([0.0, self.step, 0.0], dtype=np.float32),
            'dC': np.array([0.0, 0.0, self.step], dtype=np.float32),
            'ngrid': ngrid,
        }

        basis_name = self.slako_prefix.rstrip('/').split('/')[-1] if '/' in self.slako_prefix else self.slako_prefix
        if not basis_name:
            basis_name = '3ob-3-1'
        
        # Use config system to find basis file
        basis_hsd_path = get_dftb_basis_path(basis_name)
        if basis_hsd_path is None:
            # Fallback to old hardcoded path
            _ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
            basis_hsd_path = os.path.join(_ROOT, 'pyBall', 'DFTB', 'data', f'wfc.{basis_name}.hsd')

        basis_data = parse_wfc_hsd(basis_hsd_path)
        basis_ang = convert_wfc_to_species_list_ang(basis_data, resolution_bohr=0.04)

        rho_scf = self.projector.project_density_dense(dm_dense.astype(np.float32), self.norb_per_atom, self.orb_offsets, self.atoms_dict, self.grid_spec)
        print(f"  [Stage2] rho_scf: shape={rho_scf.shape} range=[{rho_scf.min():.4e},{rho_scf.max():.4e}] sum={rho_scf.sum():.4e}")

        coords_bohr = self.atomPos * 1.8897259886
        species_per_atom = list(range(len(self.enames)))
        geo = {
            'natoms': len(self.enames),
            'species_per_atom': species_per_atom,
            'species_names': self.enames,
            'coords_bohr': coords_bohr
        }
        rho_na = dg.project_neutral_density(geo, self.projector, self.atoms_dict, self.grid_spec, basis_ang)
        print(f"  [Stage2] rho_na:  shape={rho_na.shape} range=[{rho_na.min():.4e},{rho_na.max():.4e}] sum={rho_na.sum():.4e}")
        rho_diff = (rho_scf - rho_na).astype(np.float32)
        print(f"  [Stage2] rho_diff:range=[{rho_diff.min():.4e},{rho_diff.max():.4e}] sum={rho_diff.sum():.4e}")
        # z-profile to show where density lives
        z_profile = rho_scf.sum(axis=(0,1))
        iz_max = int(np.argmax(z_profile))
        print(f"  [Stage2] density z-peak at iz={iz_max}, z={float(self.origin[2]) + iz_max*self.step:.3f} A (mol z range [{self.atomPos[:,2].min():.2f},{self.atomPos[:,2].max():.2f}])")

        np.savez_compressed(self.cache_stage2, rho_scf=rho_scf, rho_na=rho_na, rho_diff=rho_diff, origin=self.origin, ngrid=self.ngrid)
        print(f"  Stage 2 complete and cached.")
        return rho_scf, rho_na, rho_diff

    def stage3_potentials(self, rho_scf, rho_na, rho_diff, force_recompute=False,
                          pauli_params=None, vdw_params={'C6_CO': 30.0}):
        """Stage 3: Poisson Electrostatic, Pauli Repulsion, Dispersion, and Total Field (F_total) computation."""
        if not force_recompute and os.path.exists(self.cache_stage3):
            print(f"\n[ModularPipeline] Loading Stage 3 (potentials) from cache...")
            data = np.load(self.cache_stage3)
            return data['V_ES'], data['E_pauli_field'], data['E_ES_field'], data['E_vdw'], data['F_total']
            
        print(f"\n[ModularPipeline] Computing Stage 3 (FDBM potentials)...")
        
        # Set default Pauli params based on backend
        if pauli_params is None:
            if self.backend == 'pyscf':
                from pyBall.OCL import AFM as afm_mod
                pyscf_basis_key = f"pyscf_{self.pyscf_params.get('basis', 'sto-3g')}"
                if pyscf_basis_key in afm_mod.PAULI_FITTED_DEFAULTS:
                    pauli_params = afm_mod.PAULI_FITTED_DEFAULTS[pyscf_basis_key]
                    print(f"  Using pySCF Pauli defaults for {pyscf_basis_key}: A={pauli_params['A']:.2f}, beta={pauli_params['beta']:.2f}")
                else:
                    print(f"  WARNING: No fitted Pauli params found for {pyscf_basis_key}, using default")
                    pauli_params = {'A': 1.15, 'beta': 0.36}  # Fitted for sto-3g
            else:
                from pyBall.OCL import AFM as afm_mod
                if self.basis in afm_mod.PAULI_FITTED_DEFAULTS:
                    pauli_params = afm_mod.PAULI_FITTED_DEFAULTS[self.basis]
                    print(f"  Using DFTB Pauli defaults for {self.basis}: A={pauli_params['A']:.2f}, beta={pauli_params['beta']:.2f}")
                else:
                    print(f"  WARNING: No fitted Pauli params found for {self.basis}, using default")
                    pauli_params = {'A': 787.22, 'beta': 1.2371}
        
        # Step 2: Electrostatics
        V_ES = afm.fft_poisson(rho_diff, self.step)
        
        # Load CO tip density
        target_shape = tuple(self.ngrid)
        fdata_dir = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tests', 'pyFireball', 'Fdata'))
        fdata_basis = os.path.join(fdata_dir, 'basis')
        
        # Check local co_tip_dir first (like full pipeline)
        if self.co_tip_dir is not None and os.path.isdir(self.co_tip_dir):
            print(f"  Loading precomputed CO tip from {self.co_tip_dir}...")
            co_rho_total_raw = np.load(os.path.join(self.co_tip_dir, 'co_rho_total.npy'))
            co_rho_delta_raw = np.load(os.path.join(self.co_tip_dir, 'co_rho_delta.npy'))
            print(f"  Raw CO tip shape: {co_rho_total_raw.shape}")
        else:
            # Check global cache
            cached = afm_utils._get_cached_co_tip(self.step, self.margin, fdata_dir, fdata_basis)
            if cached is not None:
                print(f"  Loading cached CO tip (step={self.step}, margin={self.margin})...")
                co_rho_total_raw, co_rho_delta_raw = cached
                print(f"  Raw CO tip shape: {co_rho_total_raw.shape}")
            else:
                # Compute on-the-fly
                print(f"  Computing CO tip on-the-fly (step={self.step})...")
                co_tip_work = os.path.join(self.output_dir, 'co_tip_work')
                os.makedirs(co_tip_work, exist_ok=True)
                co_grid_spec, co_ngrid, co_origin = afm_utils._compute_co_tip_grid(step=self.step, margin=self.margin)
                print(f"  CO grid: ngrid={co_ngrid}, origin={co_origin}")
                afm_utils._call_compute_co_tip_script(co_tip_work, co_grid_spec, self.step, 100, fdata_dir, fdata_basis)
                co_rho_total_raw = np.load(os.path.join(co_tip_work, 'co_rho_total.npy'))
                co_rho_delta_raw = np.load(os.path.join(co_tip_work, 'co_rho_delta.npy'))
                # Save to global cache for reuse
                afm_utils._save_cached_co_tip(co_rho_total_raw, co_rho_delta_raw, self.step, self.margin, fdata_dir, fdata_basis)
                print(f"  CO tip saved to global cache")
            
        co_rho_total = afm_utils._pad_and_roll_co_tip(co_rho_total_raw, target_shape)
        co_rho_delta = afm_utils._pad_and_roll_co_tip(co_rho_delta_raw, target_shape)
        
        # Step 3: Pauli repulsion
        overlap_raw = afm.compute_pauli_overlap(rho_scf, co_rho_total, self.step, tip_rolled=True)
        A_pauli = pauli_params.get('A', 787.22)
        beta_pauli = pauli_params.get('beta', 1.2371)
        E_pauli_field = afm.scale_pauli_field(overlap_raw, self.step, A_pauli, beta_pauli, return_grads=False)
        
        # Step 4: Electrostatic convolution
        E_ES_field = afm.compute_es_conv_field(V_ES, co_rho_delta, self.step, tip_rolled=True, return_grads=False)
        
        # Step 5: Dispersion
        E_vdw = afm.compute_dispersion_grid(
            self.atomPos, self.atomTypes, self.origin, self.step, self.ngrid,
            C6_CO=vdw_params['C6_CO'], return_grads=False
        )
        
        # Total force field
        E_total = E_pauli_field + E_ES_field + E_vdw
        
        # Use GPU for gradient computation
        afmulator = afm.AFMulator(use_morse=False, nloc=32)
        F_total = afmulator.compute_gradient_cl(E_total, self.step, bAlloc=True)
        
        np.savez_compressed(self.cache_stage3, V_ES=V_ES, E_pauli_field=E_pauli_field,
                            E_ES_field=E_ES_field, E_vdw=E_vdw, F_total=F_total)
        print(f"  Stage 3 complete and cached.")
        return V_ES, E_pauli_field, E_ES_field, E_vdw, F_total

    def stage4_relax(self, F_total, force_recompute=False, relax_params={'K_LAT': 0.5}, ppm_mode=True):
        """Stage 4: Probe-particle MD relaxation (yielding AFM signal and tip displacements)."""
        if not force_recompute and os.path.exists(self.cache_stage4):
            print(f"\n[ModularPipeline] Loading Stage 4 (relaxation) from cache...")
            data = np.load(self.cache_stage4)
            tip_disp = {'dx': data['tip_disp_dx'], 'dy': data['tip_disp_dy'], 'dz': data['tip_disp_dz']}
            return data['df'], tip_disp, data['FEs_relax']
            
        print(f"\n[ModularPipeline] Running Stage 4 (probe relaxation)...")
        afmulator = afm.AFMulator(use_morse=False, nloc=32)
        
        df, tip_disp = afm_utils.compose_and_relax_total(
            F_total,
            self.scan_xs, self.scan_ys, self.heights,
            self.origin, self.step, self.atomPos, K_LAT=relax_params['K_LAT'],
            use_gpu_relax=True, ppm_mode=ppm_mode, afmulator=afmulator
        )
        
        # Output FEs_relax from composition
        mol_z = float(self.atomPos[:,2].max())
        if ppm_mode:
            relax_pars_ppm = [0.1, 0.1, 0.03, 0.1]
            FEs_relax, _ = afmulator.scan_fdbm(
                self.scan_xs, self.scan_ys, self.heights, mol_z=mol_z,
                K_LAT=relax_params['K_LAT'], relax_pars=relax_pars_ppm
            )
        else:
            FEs_relax, _ = afmulator.scan_fdbm_2d(
                self.scan_xs, self.scan_ys, self.heights, mol_z=mol_z,
                K_LAT=relax_params['K_LAT']
            )
            
        np.savez_compressed(self.cache_stage4, df=df, tip_disp_dx=tip_disp['dx'],
                            tip_disp_dy=tip_disp['dy'], tip_disp_dz=tip_disp['dz'], FEs_relax=FEs_relax)
        print(f"  Stage 4 complete and cached.")
        return df, tip_disp, FEs_relax

    def stage5_stm(self, eigvecs, eigvals, lumo_offsets=[1, 2, 3], mo_indices=None,
                  field='ldos', use_exp_basis=True, exp_beta=1.0, exp_r0=3.0):
        """Stage 5: Standard STM projection on height slices."""
        if self.backend == 'pyscf':
            raise NotImplementedError(
                "STM imaging (Stage 5) is not yet supported with pySCF backend. "
                "The current Phase 1 implementation computes density directly on grid "
                "but doesn't export MOs for GPU projection. Use DFTB backend for STM, "
                "or implement Phase 2 (GTO GPU projection) for pySCF STM support."
            )
        print(f"\n[ModularPipeline] Running Stage 5 (Standard STM)...")
        return afm_utils.compute_stm(
            self.projector, eigvecs, eigvals, self.scan_xs, self.scan_ys, self.heights,
            self.norb_per_atom, self.orb_offsets, self.atoms_dict,
            lumo_offsets=lumo_offsets, mo_indices=mo_indices, field=field,
            use_exp_basis=use_exp_basis, exp_beta=exp_beta, exp_r0=exp_r0
        )

    def stage6_br_stm(self, eigvecs, eigvals, tip_disp, lumo_offsets=[1, 2, 3],
                      mo_indices=None, field='ldos', use_exp_basis=True, exp_beta=1.0, exp_r0=3.0):
        """Stage 6: Bond-Resolved STM (STM at AFM-relaxed tip positions)."""
        if self.backend == 'pyscf':
            raise NotImplementedError(
                "Bond-Resolved STM (Stage 6) is not yet supported with pySCF backend. "
                "The current Phase 1 implementation computes density directly on grid "
                "but doesn't export MOs for GPU projection. Use DFTB backend for STM, "
                "or implement Phase 2 (GTO GPU projection) for pySCF STM support."
            )
        print(f"\n[ModularPipeline] Running Stage 6 (Bond-Resolved STM)...")
        return afm_utils.compute_bond_resolved_stm(
            self.projector, eigvecs, eigvals, self.scan_xs, self.scan_ys, self.heights,
            tip_disp, self.norb_per_atom, self.orb_offsets, self.atoms_dict,
            lumo_offsets=lumo_offsets, mo_indices=mo_indices, field=field,
            use_exp_basis=use_exp_basis, exp_beta=exp_beta, exp_r0=exp_r0
        )
