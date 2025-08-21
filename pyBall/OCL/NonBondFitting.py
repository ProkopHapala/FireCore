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

class FittingDriver(OpenCLBase):

    default_REQH=(1.7, 0.1, 0.0, 0.0)
    alphaMorse = 1.5

    # def __init__(self, platform_index=0, device_index=0):
    #     super().__init__(platform_index=platform_index, device_index=device_index)
    #     self.atom_type_map = {}
    #     self.atom_type_names = []
    #     self.dof_definitions = []
    #     self.tREQHs_base = None

    def __init__(self, nloc=32, perBatch=10, verbose=0):
        # Initialize the base class
        super().__init__(nloc=nloc, device_index=0)
        
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
                # We assume Q (charge) and H (H-bond) defaults are 0 for now,
                # as they are typically derived or fitted.
                RvdW = float(parts[9])
                EvdW = float(parts[10])
                
                self.base_params[type_name] = {'R': RvdW, 'E': EvdW, 'Q': 0.0, 'H': 0.0}
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

        for dof in self.dof_definitions:
            # SAFEGUARD 2: Perform the same check for the main assembly logic.
            # This is where the original crash occurred.
            type_idx = self.atom_type_map.get(dof['typename'])

            if type_idx is None:
                # The type name from dofSelection.dat was not found in any of the
                # loaded XYZ files. Print a warning and skip this DOF.
                print(f"Warning: DOF typename '{dof['typename']}' not found in input data. Skipping.")
                continue # Move to the next DOF in the list

            # Find all atoms across all samples that this DOF applies to
            atom_indices = np.where(self.host_atypes == type_idx)[0]
            
            # It's also possible the type exists in principle, but no atoms of that
            # type are in this specific dataset. In that case, we also skip.
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

        self.host_dof_nis = np.array(dof_nis_list, dtype=np.int32)
        self.host_dof_to_atom = np.array(dof_to_atom_list, dtype=np.int32)
        self.host_dof_coeffs = np.array(dof_coeffs_list, dtype=np.float32)
        
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

        # Define all required buffers
        buffs = {
            "ranges":           self.host_ranges.nbytes,
            "tREQHs":           self.tREQHs_base.nbytes,
            "atypes":           self.host_atypes.nbytes,
            "ieps":             self.host_ieps.nbytes,
            "atoms":            self.host_atoms.nbytes,
            "ErefW":            self.host_ErefW.nbytes,
            "dEdREQs":          self.host_atoms.nbytes,
            "fDOFs":            self.n_dofs*4,
            "DOFnis":           self.host_dof_nis.nbytes,
            "DOFtoAtom":        self.host_dof_to_atom.nbytes,
            "DOFcofefs":        self.host_dof_coeffs.nbytes,
            "DOFs":             self.n_dofs*4,
            "regParams":        self.host_regParams.nbytes
        }
        self.try_make_buffers(buffs)
        
        # Upload all static data to the GPU
        # Apply any per-atom charge overrides before upload
        self.apply_charge_overrides()
        self.toGPU_(self.ranges_buff,    self.host_ranges)
        self.toGPU_(self.atypes_buff,    self.host_atypes)
        self.toGPU_(self.ieps_buff,      self.host_ieps)
        self.toGPU_(self.atoms_buff,     self.host_atoms)
        self.toGPU_(self.tREQHs_buff,    self.tREQHs_base)
        self.toGPU_(self.ErefW_buff,     self.host_ErefW)
        self.toGPU_(self.DOFnis_buff,    self.host_dof_nis)
        self.toGPU_(self.DOFtoAtom_buff, self.host_dof_to_atom)
        self.toGPU_(self.DOFcofefs_buff, self.host_dof_coeffs)
        self.toGPU_(self.regParams_buff, self.host_regParams)
        
        self.set_kernel_args()
        if self.verbose>0:
            print("GPU buffers allocated and static data uploaded.")

    def set_kernel_args(self):
        """Sets arguments for the final two kernels."""
        
        # __kernel void evalSampleDerivatives(
        #     //const int4 ns,                  
        #     __global int4*    ranges,    // 1: [nsample]  (i0,ni, j0,nj) star and number of atoms in fragments 1,2 
        #     __global float4*  tREQHs,    // 2: [ntypes]   non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
        #     __global int*     atypes,    // 3: [natomTot] atom types
        #     //__global int*     hosts,   // [natomTot] index of host atoms (it electron pair), -1 if it is not electron pair
        #     __global int2*    ieps,      // 4: [natomTot] {iep1,iep2} index of electron pair type ( used to subtract  charge)
        #     __global float4*  atoms,     // 5: [natomTot] positions and REPS charge of each atom (x,y,z,Q) 
        #     __global float4*  dEdREQs,   // 6: [natomTot] output derivatives of type REQH parameters
        #     __global float2*  ErefW      // 7: [nsamp]   {E,W} reference energies and weights for each sample (molecule,system)
        # ){
        # Prefer original kernel by default; use templated kernel only when explicitly compiled
        use_templated = getattr(self, 'use_template', False) and \
                        ('evalSampleDerivatives_template' in getattr(self, 'kernelheaders', {}))
        if use_templated:
            self.eval_kern = self.prg.evalSampleDerivatives_template
            globParams = np.array([self.alphaMorse, 0.0, 0.0, 0.0], dtype=np.float32)
            self.eval_kern.set_args(
                self.ranges_buff, self.tREQHs_buff, self.atypes_buff, self.ieps_buff,
                self.atoms_buff, self.dEdREQs_buff, self.ErefW_buff, globParams
            )
        else:
            self.eval_kern = self.prg.evalSampleDerivatives
            self.eval_kern.set_args(
                self.ranges_buff, self.tREQHs_buff, self.atypes_buff, self.ieps_buff,
                self.atoms_buff, self.dEdREQs_buff, self.ErefW_buff
            )

        # __kernel void assembleAndRegularize(
        #     int nDOFs,                       // 1: number of DOFs
        #     __global float*   fDOFs,         // 2: [nDOFs]    derivatives of REQH parameters
        #     __global int2*    DOFnis,        // 3: [nDOFs]    (i0,ni) star and number of atoms in fragments 1,2    
        #     __global int*     DOFtoAtom,     // 4: [nInds]    list of atom indexes relevant for each DOF (non-uniform ranges indexed by DOFnis)
        #     __global float4*  DOFcofefs,     // 5: [nInds]    factors for update of each DOF from the atom dEdREQH parameters   fDOFi = dot( DOFcofefs[i], dEdREQH[i] )
        #     __global float4*  dEdREQs,       // 6: [nAtomTot] derivatives of REQH parameters
        #     __global const float*   DOFs,     // 7: [nDOFs]    current values of DOFs (need only for regularization)
        #     __global const float8*  regParams // 8: [nDOFs]   {min,max, xlo,xhi, Klo,Khi, K0,x0}
        # ){
        # This kernel assembles physical forces and adds regularization
        self.assemble_kern = self.prg.assembleAndRegularize
        self.assemble_kern.set_args(
            np.int32(self.n_dofs),
            self.fDOFs_buff,
            self.DOFnis_buff,
            self.DOFtoAtom_buff,
            self.DOFcofefs_buff,
            self.dEdREQs_buff,
            self.DOFs_buff,
            self.regParams_buff
        )
        if self.verbose>0:
            print("Kernel arguments pre-set.")

    def compile_with_model(self, macros=None, output_path=None, bPrint=False):
        """
        Build the OpenCL program with a model-specific pair accumulation snippet injected
        into the templated kernel `evalSampleDerivatives_template`.

        The snippet should be passed via the 'macros' dict with key 'MODEL_PAIR_ACCUMULATION',
        e.g. macros={ 'MODEL_PAIR_ACCUMULATION': "/* code updating fREQi and Ei */" }.

        Variables in scope for the snippet inside the i–j loop:
          atomi, atomj (float4); REQi, REQj (float4); dij (float3); r (float); ir (float)
          fREQi (float4 accumulator of derivatives for atom i); Ei (float energy accumulator)
        """
        if macros is None:
            macros = {}
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path = "../../cpp/common_resources/cl/FitREQ.cl"
        source_path = os.path.abspath(os.path.join(base_path, rel_path))

        # Preprocess source with macro substitution for the template marker
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
        if bPrint:
            print(f"compile_with_model(): kernels available: {list(self.kernelheaders.keys())}")
        # Reset args to ensure we bind to the (possibly new) eval kernel object
        # Only if buffers are already allocated; otherwise let init_and_upload() call set_kernel_args()
        needed = [
            'ranges_buff','tREQHs_buff','atypes_buff','ieps_buff',
            'atoms_buff','dEdREQs_buff','ErefW_buff',
            'fDOFs_buff','DOFnis_buff','DOFtoAtom_buff','DOFcofefs_buff','DOFs_buff','regParams_buff'
        ]
        if all(hasattr(self, n) for n in needed):
            self.set_kernel_args()
        else:
            if bPrint:
                print("compile_with_model(): Buffers not yet allocated; kernel args will be set later.")

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

    def set_energy_kernel_args(self):
        """Bind args for evalSampleEnergy_template."""
        if 'evalSampleEnergy_template' not in getattr(self, 'kernelheaders', {}):
            raise RuntimeError("Energy kernel not available. Call compile_energy_with_model() first.")
        self.energy_kern = self.prg.evalSampleEnergy_template
        globParams = np.array( [self.alphaMorse,0.0,0.0,0.0], dtype=np.float32)
        # Ensure Emols buffer exists for the full init path as well (not only energy-only path)
        if not hasattr(self, 'Emols_buff'):
            self.try_make_buffers({"Emols": self.n_samples*4})
        self.energy_kern.set_args(
            self.ranges_buff, self.tREQHs_buff, self.atypes_buff, self.ieps_buff,
            self.atoms_buff, self.Emols_buff, globParams
        )

    def evaluate_energies(self, workgroup_size=32):
        """Run energy kernel and return per-sample energies as np.float32 array."""
        # Lazily allocate Emols buffer and bind args if needed (when using full init path)
        if not hasattr(self, 'Emols_buff'):
            self.try_make_buffers({"Emols": self.n_samples*4})
            self.set_energy_kernel_args()
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

    def get_forces(self, dofs_vec, bCheck=False):
        """Calculates the total force (physical + regularization) for a given DOF vector."""
        dofs_vec_f32 = dofs_vec.astype(np.float32)
        # 1) Upload DOFs (used by regularization)
        self.toGPU_(self.DOFs_buff, dofs_vec_f32)
        # 2) Refresh and upload tREQHs (used by per-sample eval)
        current_tREQHs = self.update_tREQHs_from_dofs(dofs_vec_f32)
        self.toGPU_(self.tREQHs_buff, current_tREQHs)
        # 3) Kernels
        workgroup_size_eval = 32
        workgroup_size_assemble = 128
        t0 = time.perf_counter() if self.verbose>0 else 0.0
        cl.enqueue_nd_range_kernel(self.queue, self.eval_kern, (self.n_samples * workgroup_size_eval,), (workgroup_size_eval,))
        if self.verbose>0:
            self.queue.finish(); t1 = time.perf_counter()
        cl.enqueue_nd_range_kernel(self.queue, self.assemble_kern, (self.n_dofs * workgroup_size_assemble,), (workgroup_size_assemble,))
        if self.verbose>0:
            self.queue.finish(); t2 = time.perf_counter()
        # 4) Download
        forces = self.fromGPU_(self.fDOFs_buff, shape=(self.n_dofs,))
        self.queue.finish()
        if self.verbose>0:
            t3 = time.perf_counter()
            print(f"get_forces(): dt_eval={t1-t0:.4e}s dt_assemble={t2-t1:.4e}s dt_download={t3-t2:.4e}s")
        if bCheck:
            vmin = float(np.min(forces)); vmax = float(np.max(forces))
            if (vmin!=vmin) or (vmax!=vmax) or (vmin==0) or (vmax==0):
                print("ERROR: get_forces() forces are not valid:", vmin, vmax)
                exit(1)
            if self.verbose>0:
                print("get_forces(): Forces are OK: min", vmin, "max", vmax)
        return forces

    def evaluate_objective(self, dofs_vec):
        """Scalar objective consistent with derivative kernel + regularization.

        J(dofs) = 0.5 * sum_s W_s * (Emol_s(dofs) - Eref_s)^2  -  E_reg(dofs)

        where E_reg uses the same parameters and clamping as the OpenCL
        assembleAndRegularize kernel:
            x  = clamp(x, [min,max])
            E_reg_i = 0.5*K0*(x-x0)^2 + 0.5*Klo*max(0, xlo-x)^2 + 0.5*Khi*max(0, x-xhi)^2

        Note: the minus sign on E_reg matches the current kernel which subtracts
        the regularization gradient contributions. This ensures analytic and
        numerical derivatives agree for the derivative test.
        """
        # Upload DOFs (for regularization kernel path) and update type params
        x = np.asarray(dofs_vec, dtype=np.float32)
        if hasattr(self, 'DOFs_buff'):
            self.toGPU_(self.DOFs_buff, x)
        T = self.update_tREQHs_from_dofs(x)
        self.toGPU_(self.tREQHs_buff, T)

        # Per-sample energies from energy kernel
        Em = self.evaluate_energies().astype(np.float64)
        EW = self.host_ErefW
        if EW.ndim == 1:
            Eref = EW.astype(np.float64)
            Wv   = np.ones_like(Em, dtype=np.float64)
        else:
            Eref = EW[:, 0].astype(np.float64)
            Wv   = EW[:, 1].astype(np.float64)

        d = Em - Eref
        J_phys = 0.5 * np.dot(Wv, d*d)

        # Regularization energy to mirror assembleAndRegularize (with clamping)
        J_reg = 0.0
        if hasattr(self, 'host_regParams') and self.host_regParams is not None and self.host_regParams.size > 0:
            p = self.host_regParams.astype(np.float64)
            xmin, xmax = p[:, 0], p[:, 1]
            xlo,  xhi  = p[:, 2], p[:, 3]
            Klo,  Khi  = p[:, 4], p[:, 5]
            K0,   x0   = p[:, 6], p[:, 7]
            xc = np.clip(x.astype(np.float64), xmin, xmax)
            J_reg = 0.5*np.sum(K0*(xc - x0)**2) \
                  + 0.5*np.sum(Klo*np.maximum(0.0, xlo - xc)**2) \
                  + 0.5*np.sum(Khi*np.maximum(0.0, xc - xhi)**2)

        # Match kernel sign convention (kernel subtracts reg gradient)
        return float(J_phys - J_reg)

    def evaluate_grad_fd(self, x, eps=1e-4, scheme='central'):
        """Finite-difference gradient of objective using energy kernel."""
        x0 = np.asarray(x, dtype=np.float32)
        n = x0.size
        g = np.zeros_like(x0)
        if scheme != 'central':
            scheme = 'central'
        for i in range(n):
            xp = x0.copy(); xp[i] += eps
            xm = x0.copy(); xm[i] -= eps
            Jp = self.evaluate_objective(xp)
            Jm = self.evaluate_objective(xm)
            g[i] = (Jp - Jm) / (2.0*eps)
        return g

def run_derivative_test(
    model_name='MODEL_MorseQ_PAIR',
    energy_model_name='ENERGY_MorseQ_PAIR',
    atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
    xyz_file="/home/prokop/git/FireCore/tests/tFitREQ/HHalogens/porcessed/HF-A1_HF-D1.xyz",
    dof_file="/home/prokop/git/FireCore/tests/tFitREQ/dofSelection_MorseSR.dat",
    eps=1e-4,
    verbose=0,
):
    """Compile templated derivative and energy kernels, then compare analytical
    gradient from get_forces() with finite-difference gradient from evaluate_objective()."""
    print("\n--- Running derivative test ---")
    drv = FittingDriver(verbose=verbose)
    drv.load_atom_types(atom_types_file)
    drv.load_data(xyz_file)
    drv.load_dofs(dof_file)
    drv.init_and_upload()

    # Inject both derivative and energy snippets in one build
    base_path = os.path.dirname(os.path.abspath(__file__))
    forces_path = os.path.abspath(os.path.join(base_path, "../../cpp/common_resources/cl/Forces.cl"))
    macro_der = extract_macro_block(forces_path, model_name)
    macro_en  = extract_macro_block(forces_path, energy_model_name)
    macros = { 'MODEL_PAIR_ACCUMULATION': macro_der, 'MODEL_PAIR_ENERGY': macro_en }
    drv.compile_with_model(macros=macros, bPrint=False)
    # Bind energy kernel args (derivative kernel args are bound in set_kernel_args() during init/compile)
    drv.set_energy_kernel_args()

    x0 = np.array([d['xstart'] for d in drv.dof_definitions], dtype=np.float32)
    g_an = drv.get_forces(x0)
    g_fd = drv.evaluate_grad_fd(x0, eps=eps)

    diff = g_an - g_fd
    nrm2 = float(np.linalg.norm(diff))
    linf = float(np.max(np.abs(diff))) if diff.size>0 else 0.0
    rel  = nrm2 / (float(np.linalg.norm(g_fd)) + 1e-12)
    print(f"Derivative test: L2={nrm2:.3e} Linf={linf:.3e} Rel={rel:.3e}")
    if verbose>0:
        print("g_an[:8]", g_an[:8])
        print("g_fd[:8]", g_fd[:8])
    return g_an, g_fd

def optimizer_FIRE(driver, initial_dofs, max_steps=1000, dt_start=0.01, fmax=1e-4):
    """
    Performs optimization using the FIRE (Fast Inertial Relaxation Engine) algorithm.
    """
    print("\n--- Starting FIRE Optimization ---")
    # FIRE parameters
    N_min=5; finc=1.1; fdec=0.5; alpha_start=0.1; fa=0.99
    
    dofs = initial_dofs.copy()
    vels = np.zeros_like(dofs)
    alpha = alpha_start
    dt = dt_start
    steps_since_negative_P = 0

    for i_step in range(max_steps):
        forces = driver.get_forces(dofs)
        force_norm = np.linalg.norm(forces)

        if i_step % 10 == 0:
            print(f"FIRE Step {i_step:04d} | Force Norm = {force_norm:.6e} | dt = {dt:.4e} | alpha = {alpha:.4f}")

        if force_norm < fmax:
            print(f"FIRE converged in {i_step} steps. Final force norm: {force_norm:.6e}")
            break

        P = np.dot(forces, vels)

        vels_norm = np.linalg.norm(vels)
        if P > 0:
            steps_since_negative_P += 1
            if steps_since_negative_P > N_min:
                dt = min(dt * finc, 0.1) # Max dt
                alpha *= fa
            # MD step
            vels = (1.0 - alpha) * vels + alpha * forces / force_norm * vels_norm
        else:
            steps_since_negative_P = 0
            dt *= fdec
            alpha = alpha_start
            vels[:] = 0.0 # Stop and reset
        
        dofs += vels * dt
    
    if i_step == max_steps - 1:
        print("FIRE finished due to max steps.")

    return dofs

def extract_macro_block(file_path, macro_name):
    """Extracts a macro code block from Forces.cl delimited by
    a line '//>>>macro <macro_name>' followed by a brace-balanced block.
    Returns the block including the surrounding braces, suitable for injection.
    """
    with open(file_path, 'r') as f:
        s = f.read()
    tag = f"//>>>macro {macro_name}"
    i = s.find(tag)
    if i < 0:
        raise ValueError(f"Macro tag not found: {tag}")
    j = s.find('{', i)
    if j < 0:
        raise ValueError(f"Opening '{{' not found for macro: {macro_name}")
    depth = 0
    k = j
    while k < len(s):
        c = s[k]
        if c == '{': depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                k += 1
                break
        k += 1
    if depth != 0:
        raise ValueError(f"Unbalanced braces while parsing macro: {macro_name}")
    return s[j:k]

def setup_driver(
    model_name='ENERGY_MorseQ_PAIR',
    atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
    verbose=0,
):
    # Create driver, load atom types, and compile the energy kernel template.
    # Do NOT upload buffers here; that requires XYZ data loaded first.
    drv = FittingDriver(verbose=verbose)
    drv.load_atom_types(atom_types_file)
    base_path   = os.path.dirname(os.path.abspath(__file__))
    forces_path = os.path.abspath(os.path.join(base_path, "../../cpp/common_resources/cl/Forces.cl"))
    macro_code  = extract_macro_block(forces_path, model_name)
    macros      = { 'MODEL_PAIR_ENERGY': macro_code }
    drv.compile_energy_with_model(macros=macros, bPrint=False)
    return drv


def run_energy_imshow(model_name='ENERGY_MorseQ_PAIR',
                      xyz_file="/home/prokop/git/FireCore/tests/tFitREQ/HHalogens/HBr-D1_HBr-A1.xyz",
                      atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
                      kcal=False, sym=True, bColorbar=True,
                      verbose=0, show=True, save=None, cmap='bwr', lines=True,
                      rmax=8.0, drv=None):
    """Evaluate energy-only kernel on a packed XYZ and show 3-panel imshow:
    Reference, Model, Difference. Uses helper functions from tests/tFitREQ/split_scan_imshow_new.py
    to parse headers, detect rows, reshape to grid, and plot.
    """
    import matplotlib.pyplot as plt
    from tests.tFitREQ.split_scan_imshow_new import (
        parse_headers_ra, read_scan_atomicutils, parse_xyz_blocks,
        compute_ra_vec, detect_rows_by_r, reshape_to_grid,
        plot_imshow, compute_shift_from_grid
    )

    # 1) Extract reference energies and geometry (r, a) from headers
    Eh, Rh, Ah = parse_headers_ra(xyz_file)
    if Eh.size == 0:
        raise RuntimeError(f"No reference energies parsed from headers: {xyz_file}")

    title=os.path.basename(xyz_file)

    if drv is None:
        drv = setup_driver(model_name=model_name, atom_types_file=atom_types_file, verbose=verbose)
    # Allow user to call drv.load_data() beforehand to set per-atom hacks
    if not hasattr(drv, 'host_ranges'):
        drv.load_data(xyz_file)

    # 2) Evaluate model energies in the same order
    # Now that data is loaded, allocate/upload minimal buffers and bind kernel args
    drv.init_and_upload_energy_only()
    drv.set_energy_kernel_args()


    Em = drv.evaluate_energies()
    if verbose>0: print(f"Model energies: min={np.nanmin(Em):.6e} max={np.nanmax(Em):.6e}")

    if Em.shape[0] != Eh.shape[0]:
        print(f"Warning: Em(n={Em.shape[0]}) != Eh(n={Eh.shape[0]}). Proceeding with min length.")
        n = min(Em.shape[0], Eh.shape[0])
        Em = Em[:n]
        Eh = Eh[:n]
        Rh = Rh[:n]
        Ah = Ah[:n]

    # 3) Derive geometry r/a from atomic positions with header overrides; then row detection
    Es_geo, Ps = read_scan_atomicutils(xyz_file)
    if Es_geo.size == 0:
        # Fallback to local parser for positions if atomicUtils isn't available
        _, _, Ps = parse_xyz_blocks(xyz_file, natoms=None)
    r, a = compute_ra_vec(Ps, signed=True)
    if Rh.size == r.size and np.any(np.isfinite(Rh)):
        r = np.where(np.isfinite(Rh), Rh, r)
    if Ah.size == a.size and np.any(np.isfinite(Ah)):
        a = np.where(np.isfinite(Ah), Ah, a)
    rows, _ = detect_rows_by_r(r)
    Vref, Rg, Arow, rv = reshape_to_grid(Eh, r, a, rows)
    Vmod, _,   _,   _  = reshape_to_grid(Em, r, a, rows)
    if verbose>0: print(f"Grids: Vref shape={Vref.shape}, Vmod shape={Vmod.shape}")

    # 4) Reference alignment: subtract a common reference (asymptotic min at max r col)
    sref = compute_shift_from_grid(Vref)
    print("sref", sref)
    Vref -= sref    # shift just reference, not model!
    Vdif  = Vmod - Vref

    Vref_min = np.nanmin(Vref)
    Vmod_min = np.nanmin(Vmod)

    if verbose>0:
        print("Vref[:,0]", Vref[:,0])
        print("Vref[0,:]", Vref[0,:])
        print("Vmod[:,0]", Vmod[:,0])
        print("Vmod[0,:]", Vmod[0,:])
        print("Vref.shape , min, max", Vref.shape, np.nanmin(Vref), np.nanmax(Vref))
        print("Vmod.shape , min, max", Vmod.shape, np.nanmin(Vmod), np.nanmax(Vmod))

    # 5) Plot 3 panels via helper
    fig, axs = plot_energy_panels( Vref, Vmod, Vdif, rv, Arow, model_name, sym=sym, cmap=cmap, bColorbar=bColorbar, title=title )

    if save is not None:
        try:
            os.makedirs(os.path.dirname(save), exist_ok=True)
        except Exception:
            pass
        plt.savefig(save, dpi=150)
        if verbose:
            print(f"Saved figure to: {save}")
    # 6) Optional line profiles via helper (each dataset computes its own y-limits etc.)
    if lines:
        # Draw ref and model profiles on the same axes; internals computed per-call
        ax = plot_energy_profile(Vref, rv, Arow, [':', ':', ':', ':'], [1.5, 1.5, 1.5, 1.5], ['radial', 'radial', 'angular', 'min over r per angle'], marker='.', rmax=rmax, kcal=kcal, sym=sym, name='ref')
        ax = plot_energy_profile(Vmod, rv, Arow, ['-', '-', '-', '-'], [0.5, 0.5, 0.5, 0.5], ['radial', 'radial', 'angular', 'min over r per angle'], marker='.', rmax=rmax, kcal=kcal, sym=sym, name='mod', ax=ax)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return Vref, Vmod, Vdif, rv, Arow

# ----------------------------
# Modular plotting helpers
# ----------------------------

def plot_energy_panels(Vref, Vmod, Vdif, rv, Arow, model_name, title=None, sym=True, cmap='bwr', bColorbar=True):
    """Plot 3-panel imshow (Reference, Model, Difference).
    Returns (emin_ref, vmax_ref, emin_mod, vmax_mod, fig, axs).
    Limits are in displayed units (kcal if requested).
    """
    import matplotlib.pyplot as plt
    from tests.tFitREQ.split_scan_imshow_new import plot_imshow
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    plot_imshow(Vref, rv, Arow, title='Reference',           ax=axs[0], bColorbar=bColorbar, cmap=cmap, bSym=sym, bByMin=True)
    plot_imshow(Vmod, rv, Arow, title=f'Model {model_name}', ax=axs[1], bColorbar=bColorbar, cmap=cmap, bSym=sym, bByMin=True)
    plot_imshow(Vdif, rv, Arow, title='Difference',          ax=axs[2], bColorbar=bColorbar, cmap=cmap, bSym=sym, bByMin=True)
    if title is not None: fig.suptitle(title)
    plt.tight_layout()
    return fig, axs


def plot_energy_profile(
    V, rv, Arow,
    ls, lw, labels,
    marker='.', ms=3.0, rmax=None,
    kcal=False, sym=True,
    row_colors=('k','r'), ang_color='g', min_color='b',
    name='', ylims=None, ax=None
):
    """Plot 4 concise profiles for a single dataset V; computes helpers internally.
    Curves mapping:
      0-> radial @ start angle, 1-> radial @ mid angle,
      2-> angular @ radius of V-min, 3-> per-angle minima.
    """
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(7, 4))

    fac = 23.060548 if kcal else 1.0
    # Rows to plot
    ny = V.shape[0] if hasattr(V, 'shape') and len(V.shape) > 0 else 0
    rows = sorted({i for i in (0, ny//2) if 0 <= i < ny})
    # Angle mapping and sorting (0..π)
    theta = np.radians(Arow)
    theta = (theta + (np.pi/2.0)) % np.pi
    srt = np.argsort(theta)
    theta_s = theta[srt]
    # Radius column at global minimum of THIS dataset
    try:
        _, ix_ang = np.unravel_index(np.nanargmin(V), V.shape)
    except Exception:
        ix_ang = None
    # Symmetric y-limits per dataset if requested and not overridden
    if (ylims is None) and sym:
        vmin=np.nanmin(V)
        ylims = (vmin, -vmin)

    rr = rv.astype(float)
    # Radial rows
    for j, irow in enumerate(rows):
        ee = (V[irow, :] * fac).astype(float)
        m = np.isfinite(rr) & np.isfinite(ee)
        if rmax is not None:
            m &= (rr <= float(rmax))
        if not np.any(m):
            continue
        yy = ee[m]
        col = row_colors[j if j < len(row_colors) else -1]
        ax.plot(rr[m], yy, color=col, ls=ls[j], lw=lw[j], marker=marker, ms=ms,
                label=f"{name} {labels[j]} @{Arow[irow]:.0f}°")
    # Angular @ ix_ang
    if ix_ang is not None:
        e_ang = (V[:, int(ix_ang)] * fac)
        yy = e_ang[srt]
        lab_r = rv[int(ix_ang)] if rv is not None else np.nan
        ax.plot(theta_s, yy, color=ang_color, ls=ls[2], lw=lw[2], marker=marker, ms=ms,
                label=f"{name} {labels[2]} @ r≈{lab_r:.2f}")
    # Min over r per angle
    try:
        e_min = np.nanmin(V, axis=1) * fac
        yy = e_min[srt]
        ax.plot(theta_s, yy, color=min_color, ls=ls[3], lw=lw[3], marker=marker,
                label=f"{name} {labels[3]}")
    except Exception:
        pass
    if ylims is not None:
        ax.set_ylim(float(ylims[0]), float(ylims[1]))
    ax.grid(True, ls=':')
    ax.set_xlabel('r [Å] / θ [rad]')
    ax.set_ylabel('E [kcal/mol]' if kcal else 'E [eV]')
    ax.legend(loc='best', fontsize=8)
    return ax



# def plot_energy_profiles(Vref, Vmod, rv, Arow, rmax=8.0, title=None, ylims=None, ax=None):
#     """Overlay 4 profiles for ref and mod by calling a single drawer per dataset:
#     - Radial at first and middle angle rows
#     - Angular at radius index of global minimum of REFERENCE (shared ix)
#     - Per-angle minima over r
#     """
#     ny = Vref.shape[0]
#     if ny <= 0:
#         return
#     # Rows to show (unique and valid)
#     rows = sorted({i for i in (0, ny//2) if 0 <= i < ny})
#     # Shared angle mapping and sorting (0..π)
#     theta = np.radians(Arow)
#     theta = (theta + (np.pi/2.0)) % np.pi
#     srt = np.argsort(theta)
#     theta_s = theta[srt]
#     # Shared angular radius column from REFERENCE minimum
#     try:
#         _, ix_ref = np.unravel_index(np.nanargmin(Vref), Vref.shape)
#     except Exception:
#         ix_ref = None

#     if ax is None:
#         fig, ax = plt.subplots(1, 1, figsize=(7, 4))

#     plot_energy_profile(ax, Vref, rv, Arow, rows, srt, theta_s, ix_ref, [':', ':', ':', ':'], [1.5, 1.5, 1.5, 1.5], ['radial', 'radial', 'angular', 'min over r per angle'], marker='.', rmax=rmax, ylims=ylims, fac=fac, name='ref')
#     plot_energy_profile(ax, Vmod, rv, Arow, rows, srt, theta_s, ix_ref, ['-', '-', '-', '-'], [0.5, 0.5, 0.5, 0.5], ['radial', 'radial', 'angular', 'min over r per angle'], marker='.', rmax=rmax, ylims=ylims, fac=fac, name='mod')

#     #suptitle
#     if title is not None: fig.suptitle(title)
#     ax.grid(True, ls=':')
#     ax.set_xlabel('r [Å] / θ [rad]')
#     ax.set_ylabel('E [kcal/mol]' if kcal else 'E [eV]')
#     ax.legend(loc='best', fontsize=8)
#     ax.set_title('Profiles: ref (:) vs model (-)')
#     plt.tight_layout()

#     return fig, ax

def run_templated_example(model_name='MODEL_MorseQ_PAIR',
                          atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
                          xyz_file="/home/prokop/git/FireCore/tests/tFitREQ/input_example.xyz",
                          dof_file="/home/prokop/git/FireCore/tests/tFitREQ/dofSelection_MorseSR.dat"):
    """Demo: compile and run the templated kernel with a selected model macro.
    Keeps the original baseline runnable separately.
    """
    print("\n--- Running templated fitting demo ---")
    drv = FittingDriver()
    drv.load_atom_types(atom_types_file)
    drv.load_data(xyz_file)
    drv.load_dofs(dof_file)
    drv.init_and_upload()

    # Load macro snippet from Forces.cl and compile templated kernel
    base_path = os.path.dirname(os.path.abspath(__file__))
    forces_path = os.path.abspath(os.path.join(base_path, "../../cpp/common_resources/cl/Forces.cl"))
    macro_code = extract_macro_block(forces_path, model_name)
    macros = { 'MODEL_PAIR_ACCUMULATION': macro_code }
    drv.compile_with_model(macros=macros, bPrint=True)

    x0 = np.array([d['xstart'] for d in drv.dof_definitions])
    x_final = optimizer_FIRE(drv, x0)

    print("\n" + "="*40)
    print(f"Templated model: {model_name}")
    print("Final Optimized DOF values:")
    print(x_final)
    return x_final

def run_energy_example(model_name='ENERGY_MorseQ_PAIR',
                       atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
                       xyz_file="/home/prokop/git/FireCore/tests/tFitREQ/HHalogens/HBr-D1_HBr-A1.xyz"):
    """Compile and run the energy-only templated kernel on provided XYZ file.
    model_name: one of ENERGY_LJQH2_PAIR, ENERGY_LJr8QH2_PAIR, ENERGY_MorseQ_PAIR
    """
    print("\n--- Running energy-only evaluation ---")
    drv = FittingDriver()
    drv.load_atom_types(atom_types_file)
    drv.load_data(xyz_file)

    # Build and upload minimal buffers
    drv.init_and_upload_energy_only()

    # Load energy-only macro and compile energy template
    base_path = os.path.dirname(os.path.abspath(__file__))
    forces_path = os.path.abspath(os.path.join(base_path, "../../cpp/common_resources/cl/Forces.cl"))
    macro_code = extract_macro_block(forces_path, model_name)
    macros = { 'MODEL_PAIR_ENERGY': macro_code }
    drv.compile_energy_with_model(macros=macros, bPrint=True)
    drv.set_energy_kernel_args()

    Emols = drv.evaluate_energies()
    print("\n" + "="*40)
    print(f"Energy model: {model_name}")
    for i, E in enumerate(Emols):
        ref = drv.host_ErefW[i] if hasattr(drv, 'host_ErefW') and len(drv.host_ErefW)>i else np.nan
        print(f"Sample {i}: E = {E:.6f}  (Eref = {ref:.6f})")
    return Emols

def run_baseline_example():
    """Wraps the existing baseline example into a function."""
    print("\n--- Running baseline example ---")
    driver = FittingDriver()
    driver.load_atom_types("/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat")
    driver.load_data("/home/prokop/git/FireCore/tests/tFitREQ/input_example.xyz")
    driver.load_dofs("/home/prokop/git/FireCore/tests/tFitREQ/dofSelection_MorseSR.dat")
    driver.init_and_upload()

    initial_dofs = np.array([d['xstart'] for d in driver.dof_definitions])
    final_dofs_fire = optimizer_FIRE(driver, initial_dofs)
    print("\n" + "="*40)
    print("Final Optimized DOF values:")
    print(final_dofs_fire)
    return final_dofs_fire

# --- Example Usage ---
if __name__ == '__main__':
    # run like
    # python -u -m pyBall.OCL.NonBondFitting

    # Default: baseline example (kept intact)
    #run_baseline_example()

    # To run the templated-kernel demo with a selected model macro, uncomment one:
    #run_templated_example('MODEL_LJQH2_PAIR')
    # run_templated_example('MODEL_LJr8QH2_PAIR')
    #run_templated_example('MODEL_MorseQ_PAIR')

    run_energy_imshow('ENERGY_MorseQ_PAIR')
    #print("="*40)