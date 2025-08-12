# fitting_driver.py
from operator import iand
import sys
import os
import re
import numpy as np
import pyopencl as cl
from sympy.ntheory import is_abundant
from .OpenCLBase import OpenCLBase

class FittingDriver(OpenCLBase):

    default_REQH=(1.7, 0.1, 0.0, 0.0)

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

        self.host_ErefW = np.array(ErefW, dtype=np.float32)
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
                self.tREQHs_base[i, 1] = params['E']
                self.tREQHs_base[i, 2] = params['Q']
                self.tREQHs_base[i, 3] = params['H']
            else:
                print(f"Warning: Atom type '{type_name}' from XYZ data not found in AtomTypes.dat. Using zeros.")

        # 2. Override with 'xstart' values for the parameters being optimized (DOFs)
        for dof in self.dof_definitions:
            type_idx = self.atom_type_map.get(dof['typename'])
            if type_idx is None:
                # This safeguard is still important
                continue
            # Override the default value with the starting value for the optimization
            self.tREQHs_base[type_idx, dof['comp']] = dof['xstart']

        # This loop populates the initial values from the 'xstart' column
        for dof in self.dof_definitions:
            # SAFEGUARD 1: Check if the atom type from the DOF file exists in our XYZ data.
            type_idx = self.atom_type_map.get(dof['typename'])
            
            # If get() returns None, the type was not found. Silently skip it here.
            if type_idx is None:
                continue
                
            self.tREQHs_base[type_idx, dof['comp']] = dof['xstart']

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
        if getattr(self, 'use_template', False) and 'evalSampleDerivatives_template' in getattr(self, 'kernelheaders', {}):
            self.eval_kern = self.prg.evalSampleDerivatives_template
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
                self.tREQHs_base[i, 1] = params['E']
                self.tREQHs_base[i, 2] = params['Q']
                self.tREQHs_base[i, 3] = params['H']
            else:
                print(f"Warning: Atom type '{type_name}' from XYZ data not found in AtomTypes.dat. Using zeros.")

    def init_and_upload_energy_only(self):
        """Allocate and upload minimal buffers for energy evaluation."""
        self.prepare_host_data_energy_only()
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
        self.toGPU_(self.atoms_buff,  self.host_atoms)
        self.toGPU_(self.tREQHs_buff, self.tREQHs_base)

    def set_energy_kernel_args(self):
        """Bind args for evalSampleEnergy_template."""
        if 'evalSampleEnergy_template' not in getattr(self, 'kernelheaders', {}):
            raise RuntimeError("Energy kernel not available. Call compile_energy_with_model() first.")
        self.energy_kern = self.prg.evalSampleEnergy_template
        self.energy_kern.set_args(
            self.ranges_buff, self.tREQHs_buff, self.atypes_buff, self.ieps_buff,
            self.atoms_buff, self.Emols_buff
        )

    def evaluate_energies(self, workgroup_size=32):
        """Run energy kernel and return per-sample energies as np.float32 array."""
        cl.enqueue_nd_range_kernel(self.queue, self.energy_kern, (self.n_samples*workgroup_size,), (workgroup_size,))
        Emols = self.fromGPU_(self.Emols_buff, shape=(self.n_samples,), dtype='f4')
        self.queue.finish()
        return Emols

    def get_forces(self, dofs_vec, bCheck=True):
        """Calculates the total force (physical + regularization) for a given DOF vector."""
        dofs_vec_f32 = dofs_vec.astype(np.float32)
        
        # 1. Upload current DOF values (needed by assembleAndRegularize)
        self.toGPU_(self.DOFs_buff, dofs_vec_f32)
        
        # 2. Update the main REQH parameter table (needed by evalAndScalePerSample)
        #current_tREQHs = self.update_tREQHs_from_dofs(dofs_vec_f32)
        #self.toGPU_(self.tREQHs_buff, current_tREQHs)

        # --- Launch Kernels ---
        workgroup_size_eval = 32
        workgroup_size_assemble = 128
        
        # 3. Launch evaluation kernel
        cl.enqueue_nd_range_kernel(self.queue, self.eval_kern,      (self.n_samples * workgroup_size_eval,), (workgroup_size_eval,))
        
        # 4. Launch assembly and regularization kernel
        cl.enqueue_nd_range_kernel(self.queue, self.assemble_kern,  (self.n_dofs * workgroup_size_assemble,), (workgroup_size_assemble,))

        # 5. Download final forces from the GPU
        forces = self.fromGPU_(self.fDOFs_buff, shape=(self.n_dofs,))
        self.queue.finish()

        if bCheck:
            vmin=np.min(forces); 
            vmax=np.max(forces);
            # check NaN, inf etc
            if(vmin!=vmin or vmax!=vmax or vmin==0 or vmax==0):
                print("ERROR: get_forces() forces are not valid:", vmin, vmax)
                exit(1)
            print("get_forces(): Forces are OK: min", vmin, "max", vmax)
            exit(0)

        return forces

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

def run_energy_imshow(model_name='ENERGY_MorseQ_PAIR',
                      xyz_file="/home/prokop/git/FireCore/tests/tFitREQ/HHalogens/HBr-D1_HBr-A1.xyz",
                      atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
                      kcal=False, sym=True, Emin=2.0, bColorbar=True,
                      verbose=0, show=True, save=None):
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

    # 2) Evaluate model energies in the same order
    drv = FittingDriver(verbose=verbose)
    drv.load_atom_types(atom_types_file)
    drv.load_data(xyz_file)
    drv.init_and_upload_energy_only()
    base_path = os.path.dirname(os.path.abspath(__file__))
    forces_path = os.path.abspath(os.path.join(base_path, "../../cpp/common_resources/cl/Forces.cl"))
    macro_code = extract_macro_block(forces_path, model_name)
    macros = { 'MODEL_PAIR_ENERGY': macro_code }
    drv.compile_energy_with_model(macros=macros, bPrint=False)
    drv.set_energy_kernel_args()
    Em = drv.evaluate_energies()
    if verbose:
        print(f"Model energies: min={np.nanmin(Em):.6e} max={np.nanmax(Em):.6e}")

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
    if verbose:
        print(f"Grids: Vref shape={Vref.shape}, Vmod shape={Vmod.shape}")

    # 4) Reference alignment: subtract a common reference (asymptotic min at max r col)
    sref = compute_shift_from_grid(Vref)
    Vref_ = Vref - sref
    Vmod_ = Vmod #- sref
    Vdif  = Vmod_ - Vref_


    print("Vref_[:,0]", Vref_[:,0])
    print("Vref_[0,:]", Vref_[0,:])
    print("Vmod_[:,0]", Vmod_[:,0])
    print("Vmod_[0,:]", Vmod_[0,:])

    print("Vref_.shape , min, max", Vref_.shape, Vref_.min(), Vref_.max())
    print("Vmod_.shape , min, max", Vmod_.shape, Vmod_.min(), Vmod_.max())

    # 5) Plot 3 panels
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    vmax = None
    emin = -abs(Emin) if sym else None
    vmax = +abs(Emin) if sym else None
    # Reference
    plot_imshow(Vref_, rv, Arow, emin=emin, vmax=vmax, title='Reference', kcal=kcal, ax=axs[0], bColorbar=bColorbar)
    # Model
    plot_imshow(Vmod_, rv, Arow, emin=emin, vmax=vmax, title=f'Model {model_name}', kcal=kcal, ax=axs[1], bColorbar=bColorbar)
    # Difference (auto symmetric around 0)
    vmax_d = np.nanmax(np.abs(Vdif))
    plot_imshow(Vdif, rv, Arow, emin=-vmax_d, vmax=+vmax_d, title='Difference', kcal=kcal, ax=axs[2], bColorbar=bColorbar)
    plt.tight_layout()
    if save is not None:
        try:
            os.makedirs(os.path.dirname(save), exist_ok=True)
        except Exception:
            pass
        plt.savefig(save, dpi=150)
        if verbose:
            print(f"Saved figure to: {save}")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return Vref_, Vmod_, Vdif, rv, Arow

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