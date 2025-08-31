# fitting_driver.py
from operator import iand
import sys
import os
import re
import numpy as np
import pyopencl as cl
import time
from sympy.ntheory import is_abundant
from .OpenCLBase import OpenCLBase

def checkSizeAndStop( var, bStop=True, name="" ):
    if len(var) == 0:
        print( name+".shape: ", var.shape, "is empty" )
        if bStop: 
            print(" => exit()")
            exit()

class FittingDriver(OpenCLBase):

    default_REQH=(1.7, 0.1, 0.0, 0.0)
    alphaMorse = 1.5

    # def __init__(self, platform_index=0, device_index=0):
    #     super().__init__(platform_index=platform_index, device_index=device_index)
    #     self.atom_type_map = {}
    #     self.atom_type_names = []
    #     self.dof_definitions = []
    #     self.tREQHs_base = None

    def __init__(self, nloc=32, perBatch=10, verbose=0, use_type_charges=False, bCompile=False, serial_mode=False):
        print("FittingDriver.__init__(): bCompile = ", bCompile)
        # Initialize the base class
        super().__init__(nloc=nloc, device_index=0)
        
        if bCompile:
            # Load the OpenCL program
            base_path = os.path.dirname(os.path.abspath(__file__))
            rel_path = "../../cpp/common_resources/cl/FitREQ.cl"
            if not self.load_program(rel_path=rel_path, base_path=base_path, bPrint=False):
                exit(1)

        self.atom_type_map = {}
        self.atom_type_names = []
        self.dof_definitions = []
        self.tREQHs_base = None  
        self.verbose = int(verbose)
        # Runtime charge source switch: 0=use per-atom atoms.w, 1=use type-based REQH.z
        self.use_type_charges = int(bool(use_type_charges))
        # Runtime kernel selection: False=parallel kernel, True=serial kernel
        self.serial_mode = bool(serial_mode)
        # Regularization control
        self.regularize_enabled = True
        self.host_regParams_base = None

        # --- Lightweight user hack hooks (non-invasive) ---
        # Type-level edits applied after reading AtomTypes.dat and before DOFs:
        #   - type_set   : list of (name, comp, value) where comp in {'R','E','Q','H'}; for 'E' value is EvdW and will be sqrt()'ed for storage
        #   - type_scale : list of (name, comp, factor); for 'E' factor scales EvdW, we apply sqrt(factor) to stored sqrt(E)
        #   - pair_like_scale : list of (nameA, nameB, factor) -> convenience: scales both A and B 'E' by sqrt(factor)
        # Per-atom charge overrides applied just before upload:
        #   - charge_overrides : dict {atom_index: charge}
        self.type_set = []
        self.type_scale = []
        self.pair_like_scale = []
        self.charge_overrides = {}
        print("FittingDriver.__init__(): END")

    @staticmethod
    def _comp_index(comp):
        if isinstance(comp, int): return comp
        m = {'R':0,'E':1,'Q':2,'H':3}
        return m[comp]

    def apply_type_overrides(self):
        """Apply simple set/scale edits to `self.tREQHs_base`.
        Semantics:
          - type_set: ('X','E',EvdW_abs) stores sqrt(EvdW_abs) in column 1; other comps stored verbatim
          - type_scale: ('X','E',fE) scales EvdW by factor fE -> multiply stored sqrt(E) by sqrt(fE);
                        other comps multiply stored value by factor
          - pair_like_scale: ('A','B',f) convenience -> as if scaling EvdW for A and B by factor f (again via sqrt on stored value)
        """
        if self.tREQHs_base is None: return
        # Pair-like -> expand into type_scale on 'E' for each participant
        if self.pair_like_scale:
            for a,b,f in self.pair_like_scale:
                sf = float(np.sqrt(float(f)))
                self.type_scale.append((a,'E',sf))
                self.type_scale.append((b,'E',sf))
            # do not keep growing on repeated calls
            self.pair_like_scale = []
        # Apply type_set
        for name, comp, val in list(self.type_set):
            if name not in self.atom_type_map: continue
            i = self.atom_type_map[name]
            ci = self._comp_index(comp)
            x = float(val)
            if ci == 1: x = float(np.sqrt(x))  # store sqrt(E)
            self.tREQHs_base[i, ci] = x
        # Apply type_scale
        for name, comp, fac in list(self.type_scale):
            if name not in self.atom_type_map: continue
            i = self.atom_type_map[name]
            ci = self._comp_index(comp)
            f = float(fac)
            if ci == 1:  # scale EvdW -> scale stored sqrt(E) by sqrt(f)
                self.tREQHs_base[i, ci] *= float(np.sqrt(f))
            else:
                self.tREQHs_base[i, ci] *= f

    def apply_charge_overrides(self):
        """Override per-atom charges in `self.host_atoms[:,3]` from `self.charge_overrides`.
        Charges are taken from atoms.w in kernels; REQH.z is only used for electron-pair subtraction.
        """
        if not self.charge_overrides: return
        for i, q in self.charge_overrides.items():
            if 0 <= int(i) < self.host_atoms.shape[0]:
                self.host_atoms[int(i), 3] = float(q)

    def dump_used_type_params(self, out=None):
        """Print a compact table of used atom types and their REQH just before upload.
        Order matches GPU buffer `tREQHs` (i indexes `self.atom_type_names`).
        E printed is EvdW (stored internally as sqrt(E)).
        """
        if out is None: out = sys.stdout
        if not hasattr(self, 'tREQHs_base') or self.tREQHs_base is None: return
        n = len(self.atom_type_names)
        print(f"# Used types (n={n}) — order == GPU tREQHs", file=out)
        print(f"#   id  name         R[Å]          E[eV]            Q              H", file=out)
        for i, name in enumerate(self.atom_type_names):
            R, sE, Q, H = self.tREQHs_base[i]
            E = float(sE)*float(sE)
            print(f"{i:4d}  {name:>4s}  {R:12.6f}  {E:12.6f}  {Q:12.6f}  {H:12.6f}", file=out)

    def dump_dof_regularization(self, out=None):
        """Print DOF regularization constraints loaded from `dof_file` for used types only."""
        if out is None: out = sys.stdout
        if not hasattr(self, 'dof_definitions') or not self.dof_definitions: return
        print(f"# DOF regularization (used types only)", file=out)
        print(f"#   id  name  comp  min        max        xlo        xhi        Klo        Khi        K0         xstart", file=out)
        for i, d in enumerate(self.dof_definitions):
            ti = self.atom_type_map.get(d['typename'])
            if ti is None: continue
            comp = d['comp'] if isinstance(d['comp'], int) else int(d['comp'])
            # Map comp index to letters for readability
            compL = {0:'R',1:'E',2:'Q',3:'H'}.get(comp, '?')
            print(f"{i:4d}  {d['typename']:>4s}   {compL}   {d['min']:10.6f} {d['max']:10.6f} {d['xlo']:10.6f} {d['xhi']:10.6f} {d['Klo']:10.6f} {d['Khi']:10.6f} {d['K0']:10.6f} {d['xstart']:10.6f}", file=out)

    def set_charge_for_local_index(self, i_local, q, frag='both'):
        """Convenience: set charge of i-th atom in every sample.
        Args:
          i_local: 0-based index within fragment (per sample)
          q:       charge value to assign
          frag:    'frag1', 'frag2', or 'both'
        Side-effect: populates self.charge_overrides and writes into self.host_atoms.
        Call after load_data() and before init_and_upload*().
        """
        if not hasattr(self, 'host_ranges'):
            raise RuntimeError("set_charge_for_local_index() requires load_data() first")
        q = float(q)
        for (i0, j0, ni, nj) in self.host_ranges:
            if frag in ('frag1','both') and 0 <= i_local < ni:
                idx = i0 + i_local
                self.charge_overrides[idx] = q
                self.host_atoms[idx,3] = q
            if frag in ('frag2','both') and 0 <= i_local < nj:
                idx = j0 + i_local
                self.charge_overrides[idx] = q
                self.host_atoms[idx,3] = q

    def set_charges_for_sample_atoms(self, pairs):
        """Batch-set charges for sample-local indices (ignoring fragments).
        pairs: list of (ia, q), ia is 0-based within each sample. Keeps charge_overrides for logging/replay.
        """
        if not hasattr(self, 'host_ranges'):
            raise RuntimeError("set_charges_for_sample_atoms() requires load_data() first")
        if not pairs: return
        if not hasattr(self, '_i0') or not hasattr(self, '_nat'):
            hr = np.asarray(self.host_ranges)
            self._i0  = hr[:, 0].astype(np.int64)
            self._nat = (hr[:, 2] + hr[:, 3]).astype(np.int64)
        i0, nat = self._i0, self._nat
        idx_blocks = [i0[ia < nat] + int(ia) for ia, _ in pairs]
        if not idx_blocks: return
        idx = np.concatenate(idx_blocks)
        qv = np.concatenate([np.full(b.shape, float(q), np.float32) for (_, q), b in zip(pairs, idx_blocks)])
        self.host_atoms[idx, 3] = qv
        self.charge_overrides.update({int(i): float(v) for i, v in zip(idx.tolist(), qv.tolist())})

    def set_charge_for_sample_atom(self, i_local, q):
        """Vectorized wrapper: set charge of a single sample-local index across all samples.
        Equivalent to set_charges_for_sample_atoms([(i_local, q)]).
        """
        self.set_charges_for_sample_atoms([(i_local, q)])

    def load_atom_types(self, atom_types_file):
        """
        Loads default non-bonded parameters from an AtomTypes.dat file.
        Stores them in a dictionary for later lookup.
        """
        if self.verbose>0:
            print(f"Loading base atom type parameters from {atom_types_file}...")
        self.base_params = {}
        with open(atom_types_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split()
                if len(parts) < 11: continue # Skip malformed lines

                type_name = parts[0]
                # REQH components from the file
                # RvdW is column 10, EvdW is column 11 (1-based index)
                # Qbase is column 12, Hb is column 13
                RvdW = float(parts[9])
                EvdW = float(parts[10])
                Qbase = float(parts[11]) if len(parts) > 11 else 0.0
                Hb = float(parts[12]) if len(parts) > 12 else 0.0
                
                self.base_params[type_name] = {'R': RvdW, 'E': EvdW, 'Q': Qbase, 'H': Hb}
        if self.verbose>0:
            print(f"Loaded base parameters for {len(self.base_params)} atom types.")

    def load_data(self, xyz_file):
        """Loads molecular geometries from a concatenated XYZ file.
        Supports comment formats:
          - "# E = <eV> eV"
          - "... n0 <int> Etot <eV> x0 <r> z <deg> <name>"
        If n0 is missing, defaults to natoms//2.
        """
        if self.verbose>0:
            print(f"Loading data from {xyz_file}...")
        atoms_list, atypes_list, ranges_list = [], [], []
        ia = 0

        type_names = []
        ErefW = []
        system_params = []

        # Regex for robust parsing
        reE  = re.compile(r"(?:#\s*E\s*=\s*|\bEtot\s+)([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")
        reN0 = re.compile(r"\bn0\s+(\d+)")
        reX  = re.compile(r"\bx0\s+([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")
        reZ  = re.compile(r"\bz\s+([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)")

        with open(xyz_file, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                line_s = line.strip()
                if line_s == '' or not line_s.isdigit():
                    # Skip stray lines until a natoms line is found
                    continue
                natoms = int(line_s)
                comment = f.readline()
                ws = comment.split()
                if self.verbose>1:
                    print(f"load_data() isys: {len(ranges_list)} comment: ", ws)

                mE = reE.search(comment)
                Etot = float(mE.group(1)) if mE else float('nan')
                mN = reN0.search(comment)
                n0 = int(mN.group(1)) if mN else natoms // 2
                mX = reX.search(comment); mZ = reZ.search(comment)
                x = float(mX.group(1)) if mX else float('nan')
                z = float(mZ.group(1)) if mZ else float('nan')
                molname = ws[-1] if len(ws) > 0 and not ws[-1].startswith('eV') else ''

                system_params.append((molname, (x, z), Etot))
                ranges_list.append((ia, ia + n0, n0, natoms - n0))
                ErefW.append(Etot)

                for _ in range(natoms):
                    parts = f.readline().strip().split()
                    if len(parts) < 4:
                        continue
                    type_name = parts[0]
                    if type_name not in self.atom_type_map:
                        self.atom_type_map[type_name] = len(self.atom_type_names)
                        self.atom_type_names.append(type_name)
                    type_names.append(type_name)
                    atypes_list.append(self.atom_type_map[type_name])
                    # Expect up to 4 floats; pad if missing to keep float4 layout
                    floats = [float(p) for p in parts[1:5]]
                    if len(floats) < 4:
                        floats += [0.0] * (4 - len(floats))
                    atoms_list.append(floats)
                ia += natoms

        # Pack reference energies and weights into float2 array: [Eref, W]
        # Default weight W=1.0 unless provided elsewhere
        self.host_ErefW = np.zeros((len(ErefW), 2), dtype=np.float32)
        self.host_ErefW[:, 0] = np.asarray(ErefW, dtype=np.float32)
        self.host_ErefW[:, 1] = 1.0
        self.atom_types_set = set(self.atom_type_names)
        if self.verbose>0:
            print("load_data() atom_types_set: ", self.atom_types_set)
        self.host_atoms    = np.array(atoms_list,  dtype=np.float32)
        self.host_atypes   = np.array(atypes_list, dtype=np.int32)
        self.host_ranges   = np.array(ranges_list, dtype=np.int32)

        host_ranges2      = self.host_ranges.copy()
        host_ranges2[:,:] = self.host_ranges[:,[1,0,3,2]] # flip i0,j0 and ni,nj
        self.host_ranges2 = host_ranges2

        print("self.host_ranges:  ", self.host_ranges)
        print("self.host_ranges2: ", self.host_ranges2)

        self.n_atoms_total = len(self.host_atoms)
        self.n_samples     = len(self.host_ranges)
        self.host_ieps     = np.full((self.n_atoms_total, 2), -1, dtype=np.int32)
        self.system_params = system_params
        if self.verbose>0:
            print(f"Loaded {self.n_samples} samples, {self.n_atoms_total} atoms total.")

        if self.verbose>2:
            print("Atoms:")
            for i in range(len(self.host_atoms)):
                print(f"{i}: {type_names[i]} {self.host_atypes[i]} {self.host_atoms[i]}")

            print("Atom types:")
            for i in range(len(self.atom_type_names)):
                if self.verbose>2:
                    print(f"{i}: {self.atom_type_names[i]} {self.atom_type_map[self.atom_type_names[i]]} ")

    def load_dofs(self, dof_file):
        """Loads the definitions of the degrees of freedom (DOFs)."""
        if self.verbose>0:
            print(f"Loading DOFs from {dof_file}...")
        with open(dof_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.strip().split()
                self.dof_definitions.append({
                    'typename': parts[0],
                    'comp': int(parts[1]),
                    'min': float(parts[2]),
                    'max': float(parts[3]),
                    'xlo': float(parts[4]),
                    'xhi': float(parts[5]),
                    'Klo': float(parts[6]),
                    'Khi': float(parts[7]),
                    'K0': float(parts[8]),
                    'xstart': float(parts[9])
                })
        self.n_dofs = len(self.dof_definitions)
        if self.verbose>0:
            print(f"Loaded {self.n_dofs} degrees of freedom.")

    def prepare_host_data(self, REQH=default_REQH):
        """
        Prepares all host-side numpy arrays needed for the kernels.
        This version safely handles DOFs for atom types not present in the input data.
        """
        # --- Build initial REQH parameter matrix ---
        n_types = len(self.atom_type_names)
        self.tREQHs_base = np.tile(np.array(REQH, dtype=np.float32), (n_types, 1))


        # Prepare host-side data for fitted parameter overrides
        n_types = len(self.atom_type_names)
        self.host_typToREQ = np.full((n_types, 4), -1, dtype=np.int32)
        
        # Set up mapping from atom types to DOF indices
        for dof_idx, dof in enumerate(self.dof_definitions):
            if dof['typename'] in self.atom_type_map:
                type_idx = self.atom_type_map[dof['typename']]
                comp = dof['comp']
                self.host_typToREQ[type_idx, comp] = dof_idx

        # Create inverse mapping: DOF index -> (atom_type, component)
        self.host_DOFtoTypeComp = np.full((self.n_dofs, 2), -1, dtype=np.int32)
        for dof_idx, dof in enumerate(self.dof_definitions):
            if dof['typename'] in self.atom_type_map:
                type_idx = self.atom_type_map[dof['typename']]
                comp = dof['comp']
                self.host_DOFtoTypeComp[dof_idx] = [type_idx, comp]

        # --- Build initial REQH parameter matrix ---
        n_types_in_data = len(self.atom_type_names)
        self.tREQHs_base = np.zeros((n_types_in_data, 4), dtype=np.float32)

        # 1. Initialize from AtomTypes.dat defaults
        for i, type_name in enumerate(self.atom_type_names):
            # Look up the base parameters for the type name found in the XYZ data
            params = self.base_params.get(type_name)
            if params:
                self.tREQHs_base[i, 0] = params['R']
                # Store sqrt(EvdW) to realize geometric mixing Eij = sqrt(Eii*Ejj) via product in kernel
                self.tREQHs_base[i, 1] = np.sqrt(params['E'])
                self.tREQHs_base[i, 2] = params['Q']
                self.tREQHs_base[i, 3] = params['H']
            else:
                print(f"Warning: Atom type '{type_name}' from XYZ data not found in AtomTypes.dat. Using zeros.")

        # 1b. Apply simple user overrides (optional)
        self.apply_type_overrides()

        # 2. Override with 'xstart' values for the parameters being optimized (DOFs)
        for dof in self.dof_definitions:
            type_idx = self.atom_type_map.get(dof['typename'])
            if type_idx is None:
                # This safeguard is still important
                continue
            # Override the default value with the starting value for the optimization
            x = dof['xstart']
            if dof['comp'] == 1: x = np.sqrt(x)  # comp 1 = EvdW -> store sqrt(E)
            self.tREQHs_base[type_idx, dof['comp']] = x

        # This loop populates the initial values from the 'xstart' column
        for dof in self.dof_definitions:
            # SAFEGUARD 1: Check if the atom type from the DOF file exists in our XYZ data.
            type_idx = self.atom_type_map.get(dof['typename'])
            
            # If get() returns None, the type was not found. Silently skip it here.
            if type_idx is None:
                continue
                
            x = dof['xstart']
            if dof['comp'] == 1: x = np.sqrt(x)  # comp 1 = EvdW -> store sqrt(E)
            self.tREQHs_base[type_idx, dof['comp']] = x

        # --- Build mapping arrays for the assembly kernel ---
        dof_to_atom_list = []
        dof_coeffs_list = []
        dof_nis_list = []

        # IMPORTANT: The derivative kernel writes dEdREQs only for fragment-1 atoms (i in [i0, i0+ni)).
        # Therefore, assemble DOF mappings must reference only fragment-1 atoms; otherwise contributions
        # from fragment-2 would read zeros and can null the assembled forces. Build a boolean mask of
        # fragment-1 atoms across all samples and use it to filter indices below.
        mask_frag1 = np.zeros(self.n_atoms_total, dtype=bool)
        for (i0, j0, ni, nj) in self.host_ranges: mask_frag1[i0:i0+ni] = True
        if self.verbose > 1:
            n_f1 = int(mask_frag1.sum())
            print(f"prepare_host_data(): restricting DOF mapping to frag-1 atoms only (count={n_f1}/{self.n_atoms_total})")

        for dof in self.dof_definitions:
            # SAFEGUARD 2: Perform the same check for the main assembly logic.
            # This is where the original crash occurred.
            type_idx = self.atom_type_map.get(dof['typename'])

            if type_idx is None:
                continue

            # Restrict to FRAGMENT-1 atoms only, matching derivative coverage.
            atom_indices = np.where((self.host_atypes == type_idx) & mask_frag1)[0]
            if len(atom_indices) == 0:
                continue

            # Define the coefficient vector for this DOF.
            coeff_vec = np.zeros(4, dtype=np.float32)
            coeff_vec[dof['comp']] = 1.0

            # For each atom affected by this DOF, add its index and coefficient.
            start_idx = len(dof_to_atom_list)
            dof_to_atom_list.extend(atom_indices)
            dof_coeffs_list.extend([coeff_vec] * len(atom_indices))
            dof_nis_list.append((start_idx, len(atom_indices)))

        if len(dof_to_atom_list) == 0:
            print("prepare_host_data() ERROR: len(dof_to_atom_list)==0 => Exiting.")
            # Diagnostics: show mismatch between XYZ atom types and DOF selection types
            try:
                types_xyz = sorted(set(self.atom_type_names))
            except Exception:
                types_xyz = []
            types_dof = [d.get('typename', '') for d in self.dof_definitions]
            print("# Atom types found in XYZ ({}): {}".format(len(types_xyz), ", ".join(types_xyz)))
            print("# Types requested by DOFs ({}): {}".format(len(types_dof), ", ".join(types_dof)))
            try:
                missing = sorted(set(types_dof) - set(types_xyz))
                if missing:
                    print("# DOF types not present in XYZ ({}): {}".format(len(missing), ", ".join(missing)))
            except Exception:
                pass
            exit(1)
            

        self.host_dof_nis     = np.array(dof_nis_list,     dtype=np.int32)
        self.host_dof_to_atom = np.array(dof_to_atom_list, dtype=np.int32)
        self.host_dof_coeffs  = np.array(dof_coeffs_list,  dtype=np.float32)
        
        # Convert reference energies to numpy array
        self.host_ErefW = np.array(self.host_ErefW, dtype=np.float32)
        
        # Initialize derivatives buffer
        self.host_dEdREQs = np.zeros((self.n_atoms_total, 4), dtype=np.float32)

    def init_and_upload(self):
        """Allocates GPU buffers for the final, most efficient workflow."""
        self.prepare_host_data()
        if self.verbose>0:
            self.dump_used_type_params()
        
        # Prepare regularization parameters for upload as a (n_dofs, 8) float array
        try:            
            self.host_regParams = np.array([
                ( d['min'], d['max'], d['xlo'], d['xhi'], d['Klo'], d['Khi'], d['K0'], d['xstart']) 
                for d in self.dof_definitions], dtype=np.float32).reshape(self.n_dofs, 8)
        except Exception as e:
            print("self.dof_definitions:")
            for d in self.dof_definitions: print(d)
            raise e

        # Keep a pristine copy to allow toggling regularization on/off
        self.host_regParams_base = self.host_regParams.copy()
        if not self.regularize_enabled:
            rp = self.host_regParams.copy(); rp[:, 4:7] = 0.0  # zero Klo,Khi,K0
            self.host_regParams = rp

        # checkSizeAndStop(self.host_dof_nis,     name="host_dof_nis")
        # checkSizeAndStop(self.host_dof_to_atom, name="host_dof_to_atom")
        # checkSizeAndStop(self.host_dof_coeffs,  name="host_dof_coeffs")
        # checkSizeAndStop(self.host_regParams,   name="host_regParams")
        # checkSizeAndStop(self.host_ranges,      name="host_ranges")
        # checkSizeAndStop(self.tREQHs_base,      name="tREQHs_base")
        # checkSizeAndStop(self.host_atypes,      name="host_atypes")
        # checkSizeAndStop(self.host_ieps,        name="host_ieps")
        # checkSizeAndStop(self.host_atoms,       name="host_atoms")
        # checkSizeAndStop(self.host_ErefW,       name="host_ErefW")
        # checkSizeAndStop(self.host_dEdREQs,     name="host_dEdREQs")

        buffs_ = {
            "ranges":           self.host_ranges,
            "ranges2":          self.host_ranges2, # ranges with swapped fragments 
            "tREQHs":           self.tREQHs_base,
            "atypes":           self.host_atypes,
            "ieps":             self.host_ieps,
            "atoms":            self.host_atoms,
            "ErefW":            self.host_ErefW,
            "dEdREQs":          self.host_dEdREQs,
            "DOFnis":           self.host_dof_nis,
            "DOFtoAtom":        self.host_dof_to_atom,
            "DOFcofefs":        self.host_dof_coeffs,
            "regParams":        self.host_regParams,
            "DOFtoTypeComp":    self.host_DOFtoTypeComp
        }
        for k,v in buffs_.items(): checkSizeAndStop( v, name=k)
        buffs = { k: v.nbytes for k,v in buffs_.items() }
        buffs.update( { "fDOFs":self.n_dofs*4, "DOFs":self.n_dofs*4,   })
        self.try_make_buffers(buffs)

        self.try_make_buffers({
            "Emols": self.n_samples*4,
            "Jmols": self.n_samples*4,
        })
        
        # Upload all static data to the GPU
        # Apply any per-atom charge overrides before upload
        self.apply_charge_overrides()
        self.toGPU_(self.ranges_buff,    self.host_ranges)
        self.toGPU_(self.ranges2_buff,   self.host_ranges2)
        self.toGPU_(self.atypes_buff,    self.host_atypes)
        self.toGPU_(self.ieps_buff,      self.host_ieps)
        self.toGPU_(self.atoms_buff,     self.host_atoms)
        self.toGPU_(self.tREQHs_buff,    self.tREQHs_base)
        self.toGPU_(self.ErefW_buff,     self.host_ErefW)
        self.toGPU_(self.DOFnis_buff,    self.host_dof_nis)
        self.toGPU_(self.DOFtoAtom_buff, self.host_dof_to_atom)
        self.toGPU_(self.DOFcofefs_buff, self.host_dof_coeffs)
        self.toGPU_(self.regParams_buff, self.host_regParams)
        self.toGPU_(self.DOFtoTypeComp_buff, self.host_DOFtoTypeComp)
        
        # Defer binding kernel arguments until after compile_with_model() injects the model.
        # compile_with_model() will call set_kernel_args() once the program is (re)built.
        if self.verbose>0:
            print("GPU buffers allocated and static data uploaded. Waiting for compile_with_model() to bind kernels.")

    def set_regularization_enabled(self, enabled=True):
        """Enable/disable regularization by zeroing/restoring stiffness terms.

        Effects:
        - Updates self.host_regParams (Klo,Khi,K0 columns) and uploads to GPU if buffer exists.
        - Stores an immutable copy in self.host_regParams_base on first call (if not already set).
        """
        self.regularize_enabled = bool(enabled)
        if self.host_regParams is None:
            return
        if self.host_regParams_base is None:
            self.host_regParams_base = self.host_regParams.copy()
        if self.regularize_enabled:
            self.host_regParams = self.host_regParams_base.copy()
        else:
            rp = self.host_regParams_base.copy(); rp[:, 4:7] = 0.0
            self.host_regParams = rp
        if hasattr(self, 'regParams_buff'):
            self.toGPU_(self.regParams_buff, self.host_regParams)
        if self.verbose>0:
            print(f"Regularization {'ENABLED' if self.regularize_enabled else 'DISABLED'}")

    def disable_regularization(self):
        self.set_regularization_enabled(False)

    def enable_regularization(self):
        self.set_regularization_enabled(True)

    def setup_derivative_kernel(self):
        """Setup and configure both derivative kernel variants (parallel and serial)."""
        if not hasattr(self, 'kernelheaders'):
            raise RuntimeError("Templated kernels not available. Call compile_with_model(macros=...) before binding args.")
            
        # Initialize both kernel variants
        self.kernel_deriv = self.prg.evalSampleDerivatives_template
        self.kernel_deriv_serial = self.prg.evalSampleDerivatives_template_serial
        
        # Common argument binding for both kernels
        globParams = np.array([self.alphaMorse, 0.0, 0.0, 0.0], dtype=np.float32)            
        self.deriv_args1 = (
            self.ranges_buff, 
            self.tREQHs_buff, 
            self.atypes_buff, 
            self.ieps_buff,
            self.atoms_buff, 
            self.dEdREQs_buff, 
            self.ErefW_buff,
            self.Emols_buff, 
            self.Jmols_buff,
            np.int32(self.use_type_charges), 
            globParams
        )

        self.deriv_args2 = (
            self.ranges2_buff, 
            self.tREQHs_buff, 
            self.atypes_buff, 
            self.ieps_buff,
            self.atoms_buff, 
            self.dEdREQs_buff, 
            self.ErefW_buff,
            self.Emols_buff, 
            self.Jmols_buff,
            np.int32(self.use_type_charges), 
            globParams
        )
        
        #self.kernel_deriv.set_args(*args)
        #self.kernel_deriv_serial.set_args(*args)

    def setup_energy_kernel(self):
        """Setup and configure the energy evaluation kernel.
        
        Sets self.energy_kern to evalSampleEnergy_template and binds its arguments.
        """
        if 'evalSampleEnergy_template' not in getattr(self, 'kernelheaders', {}):
            raise RuntimeError("Energy kernel not available. Call compile_energy_with_model() first.")
        self.energy_kern = self.prg.evalSampleEnergy_template
        globParams = np.array( [self.alphaMorse,0.0,0.0,0.0], dtype=np.float32)
        self.energy_kern.set_args(
            self.ranges_buff, 
            self.tREQHs_buff, 
            self.atypes_buff, 
            self.ieps_buff,
            self.atoms_buff, 
            self.Emols_buff,
            np.int32(self.use_type_charges), globParams
        )

    def setup_assembly_kernel(self):
        """Setup and configure the assembly and regularization kernel.
        
        Sets self.assemble_kern to assembleAndRegularize and binds its arguments.
        """
        self.assemble_kern = self.prg.assembleAndRegularize
        self.assemble_kern.set_args(
            np.int32(self.n_dofs),
            self.fDOFs_buff,
            self.DOFnis_buff,
            self.DOFtoAtom_buff,
            self.DOFcofefs_buff,
            self.dEdREQs_buff,
            self.DOFs_buff,
            self.regParams_buff,
            self.tREQHs_buff,
            self.DOFtoTypeComp_buff
        )

    def set_kernel_args(self):
        """Bind arguments for templated derivative and assembly kernels. No fallback."""
        if 'evalSampleDerivatives_template' not in getattr(self, 'kernelheaders', {}):
            raise RuntimeError("Templated kernels not available. Call compile_with_model(macros=...) before binding args.")

        # Derivative kernel (templated)
        self.setup_derivative_kernel()

        # Assembly + regularization kernel
        self.setup_assembly_kernel()
        
        # Energy kernel
        self.setup_energy_kernel()
        
        if self.verbose>0:
            print("Templated kernel arguments bound.")

    def compile_with_model(self, macros=None, output_path=None, bPrint=False):
        """
        Build the OpenCL program with a model-specific pair accumulation snippet injected
        into both templated kernels (parallel and serial variants).
        """
        if macros is None:
            macros = {}
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path = "../../cpp/common_resources/cl/FitREQ.cl"
        source_path = os.path.abspath(os.path.join(base_path, rel_path))

        # Preprocess source with macro substitution for the template markers
        pre_src = self.preprocess_opencl_source(
            source_path,
            substitutions={ 'macros': macros },
            output_path=output_path,
            bPrint=bPrint
        )

        # Build program from preprocessed source and refresh kernel headers
        self.prg = cl.Program(self.ctx, pre_src).build()
        self.kernelheaders = self.extract_kernel_headers(pre_src)
        self.use_template = True
        
        # Verify both kernel variants are available
        required_kernels = {
            'evalSampleDerivatives_template',
            'evalSampleDerivatives_template_serial',
            'assembleAndRegularize'
        }
        available = set(self.kernelheaders.keys())
        if not required_kernels.issubset(available):
            missing = required_kernels - available
            raise RuntimeError(f"Missing required kernels: {missing}")
            
        # Bind kernel args now that the program is built
        self.set_kernel_args()
        if bPrint:
            print(f"compile_with_model(): kernels available: {list(self.kernelheaders.keys())}")

    def compile_energy_with_model(self, macros=None, output_path=None, bPrint=False):
        """
        Build the OpenCL program with a model-specific pair ENERGY snippet injected
        into the templated kernel `evalSampleEnergy_template` at marker
            //<<<MODEL_PAIR_ENERGY

        Pass snippet via macros={'MODEL_PAIR_ENERGY': code}.
        """
        if macros is None:
            macros = {}
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path = "../../cpp/common_resources/cl/FitREQ.cl"
        source_path = os.path.abspath(os.path.join(base_path, rel_path))

        pre_src = self.preprocess_opencl_source(
            source_path,
            substitutions={ 'macros': macros },
            output_path=output_path,
            bPrint=bPrint
        )
        self.prg = cl.Program(self.ctx, pre_src).build()
        self.kernelheaders = self.extract_kernel_headers(pre_src)
        if bPrint:
            print(f"compile_energy_with_model(): kernels available: {list(self.kernelheaders.keys())}")
        # No args set here; energy-only path uses dedicated setters

    def prepare_host_data_energy_only(self):
        """Prepare only what's needed for energy evaluation (no DOFs, no derivatives)."""
        # Build REQH table from known atom types and base_params
        n_types_in_data = len(self.atom_type_names)
        self.tREQHs_base = np.zeros((n_types_in_data, 4), dtype=np.float32)
        for i, type_name in enumerate(self.atom_type_names):
            params = self.base_params.get(type_name)
            if params:
                self.tREQHs_base[i, 0] = params['R']
                # Store sqrt(EvdW) to realize geometric mixing via product in kernel
                self.tREQHs_base[i, 1] = np.sqrt(params['E'])
                self.tREQHs_base[i, 2] = params['Q']
                self.tREQHs_base[i, 3] = params['H']
            else:
                print(f"Warning: Atom type '{type_name}' from XYZ data not found in AtomTypes.dat. Using zeros.")
        # Optional user overrides
        self.apply_type_overrides()

    def init_and_upload_energy_only(self):
        """Allocate and upload minimal buffers for energy evaluation."""
        self.prepare_host_data_energy_only()
        if self.verbose>0:
            self.dump_used_type_params()
        buffs = {
            "ranges": self.host_ranges.nbytes,
            "tREQHs": self.tREQHs_base.nbytes,
            "atypes": self.host_atypes.nbytes,
            "ieps":   self.host_ieps.nbytes,
            "atoms":  self.host_atoms.nbytes,
            "Emols":  self.n_samples*4,
        }
        self.try_make_buffers(buffs)
        self.toGPU_(self.ranges_buff, self.host_ranges)
        self.toGPU_(self.atypes_buff, self.host_atypes)
        self.toGPU_(self.ieps_buff,   self.host_ieps)
        # Apply any per-atom charge overrides before upload
        self.apply_charge_overrides()
        self.toGPU_(self.atoms_buff,  self.host_atoms)
        self.toGPU_(self.tREQHs_buff, self.tREQHs_base)

    def set_charge_source(self, use_type_charges=False):
        """Set runtime charge source switch and rebind kernel args next time they're set.
        This affects both derivative and energy templated kernels.
        """
        self.use_type_charges = int(bool(use_type_charges))
        # If kernels already exist, rebind with updated arg
        if hasattr(self, 'kernel_deriv') and hasattr(self, 'ranges_buff'):  self._set_eval_ranges(self.ranges_buff)
        if hasattr(self, 'energy_kern') and hasattr(self, 'ranges_buff'):   self.setup_energy_kernel()

    def evaluate_energies(self, workgroup_size=32):
        """Run energy kernel and return per-sample energies as np.float32 array."""
        # Lazily allocate Emols buffer and bind args if needed (when using full init path)
        if not hasattr(self, 'Emols_buff'):
            self.try_make_buffers({"Emols": self.n_samples*4})
            self.setup_energy_kernel()
        cl.enqueue_nd_range_kernel(self.queue, self.energy_kern, (self.n_samples*workgroup_size,), (workgroup_size,))
        Emols = self.fromGPU_(self.Emols_buff, shape=(self.n_samples,), dtype='f4')
        self.queue.finish()
        return Emols

    def update_tREQHs_from_dofs(self, dofs_vec):
        """Return a fresh tREQHs array updated from DOF vector.
        DOF order follows `self.dof_definitions`. Comp 1 (EvdW) is stored as sqrt(E).
        """
        T = np.array(self.tREQHs_base, copy=True)
        for i, dof in enumerate(self.dof_definitions):
            ti = self.atom_type_map.get(dof['typename'])
            if ti is None: continue
            ci = int(dof['comp'])
            x  = float(dofs_vec[i])
            if ci == 1: x = float(np.sqrt(max(x, 0.0)))
            T[ti, ci] = x
        return T.astype(np.float32, copy=False)

    def getErrorDerivs(self, dofs_vec, niter=1, bDownload=True, bBothSides=False):
        """Calculate both objective (J) and its derivatives (g) for given DOF vector.
        
        This implementation properly sequences kernel execution:
        1. First runs assembleAndRegularize to update tREQs from DOFs
        2. Runs evalSampleDerivatives (parallel or serial) to compute sample derivatives
        3. Runs assembleAndRegularize again to assemble final derivatives
        
        Returns:
            tuple: (J, g) where J is the objective value and g is the gradient vector
        """
        # Convert input to float32 and upload DOFs
        dofs_vec_f32 = np.asarray(dofs_vec, dtype=np.float32)
        self.toGPU_(self.DOFs_buff, dofs_vec_f32)
        deriv_kern     = self.kernel_deriv_serial if self.serial_mode else self.kernel_deriv
        nloc           = 1                        if self.serial_mode else 32
        
        cl.    enqueue_nd_range_kernel( self.queue, self.assemble_kern, (self.n_dofs * 128,),     (128,)  )  # First pass: update tREQs from DOFs
        for itr in range(niter):
            #cl.enqueue_nd_range_kernel( self.queue,      deriv_kern,    (self.n_samples * nloc,), (nloc,) )  # Run derivative kernel (parallel or serial)
            #deriv_kern                (self.queue, (self.n_samples * nloc,), (nloc,), *self.deriv_args1)
            deriv_kern                (self.queue, (self.n_samples * nloc,), (nloc,), *self.deriv_args2)
            self.queue.finish()

            if bBothSides:  deriv_kern(self.queue, (self.n_samples * nloc,), (nloc,), *self.deriv_args2)
            self.queue.finish()
            cl.enqueue_nd_range_kernel( self.queue, self.assemble_kern, (self.n_dofs * 128,),     (128,)  )  # Second pass: assemble final derivatives
        
        if bDownload:
            # Download results
            self.queue.finish()
            dEdREQs = self.fromGPU_(self.dEdREQs_buff, self.host_dEdREQs )
            fDOFs   = self.fromGPU_(self.fDOFs_buff,  shape=(self.n_dofs,))
            Emols   = self.fromGPU_(self.Emols_buff, shape=(self.n_samples,))
            Jmols   = self.fromGPU_(self.Jmols_buff, shape=(self.n_samples,))

            print("dEdREQs:\n", dEdREQs);
            print("fDOFs:", fDOFs);
            print("Emols:", Emols);
            print("Jmols:", Jmols);

            J = np.sum(Jmols)            
            return J, -fDOFs  # Return objective and negative forces (gradient)
        return None, None