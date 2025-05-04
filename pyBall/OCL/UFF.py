import sys
import os
import re # Import regular expressions for parsing
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
# import pyopencl.cltypes as cltypes
import time

# Assuming clUtils exists and has roundup_global_size, get_cl_info
# If not, define simple versions:
class clUtils:
    @staticmethod
    def roundup_global_size(gs, ls):
        return (gs + ls - 1) // ls * ls
    @staticmethod
    def get_cl_info(device):
        print(f"Device: {device.name}")
        print(f"  Type: {cl.device_type.to_string(device.type)}")
        print(f"  Max Compute Units: {device.max_compute_units}")
        print(f"  Max Work Group Size: {device.max_work_group_size}")
        # Add more info as needed

clu = clUtils()

# Assuming MMFF class exists (placeholder)
# class MMFF: pass

class UFF_CL:
    """
    PyOpenCL interface for running UFF calculations on GPU,
    following the structure of the provided C++ reference.
    Extracts kernel headers automatically from the source file.
    """

    # Remove the hardcoded kernelheaders dictionary

    def __init__(self, nloc=32, cl_platform_index=0, cl_device_index=0, kernel_path=None):
        """
        Initializes the OpenCL context, queue, compiles the UFF kernels,
        and extracts kernel headers.

        Args:
            nloc (int): Default local work size (1D).
            cl_platform_index (int): Index of the OpenCL platform to use.
            cl_device_index (int): Index of the OpenCL device to use.
            kernel_path (str, optional): Explicit path to the OpenCL kernel file.
                                         If None, attempts to find it automatically.
        """
        self.nloc = nloc
        self.ctx = cl.create_some_context(answers=[0])
        self.queue = cl.CommandQueue(self.ctx)
        
        # Print device info
        clu.get_cl_info(self.ctx.devices[0])

        # Load and compile the OpenCL program - try multiple possible locations
        # Define potential paths where the kernel file might be located
        base_path = os.path.abspath(__file__)
        rel_path = "../../../cpp/common_resources/cl/UFF.cl"
        # form relative path to absolute
        kernel_path = os.path.abspath(os.path.join(base_path, rel_path))
        kernel_found = False
        if os.path.exists(kernel_path):
            with open(kernel_path, 'r') as f:
                self.prg = cl.Program(self.ctx, f.read()).build()
            kernel_found = True
            print(f"Successfully loaded kernel from: {kernel_path}")
        else:
            print(f"Kernel file not found at: {kernel_path}")
            exit(1)

        print(f"Loading UFF kernels from: {kernel_path}")
        with open(kernel_path, 'r') as f:
            self.kernel_source = f.read()

        # --- Extract Kernel Headers ---
        self.kernelheaders = self._extract_kernel_headers(self.kernel_source)
        if not self.kernelheaders:
             raise RuntimeError(f"Could not extract any kernel headers from {kernel_path}")
        print(f"Extracted headers for kernels: {list(self.kernelheaders.keys())}")

        # --- Build Program ---
        try:
            self.prg = cl.Program(self.ctx, self.kernel_source).build()
            print("UFF OpenCL program built successfully.")
        except cl.RuntimeError as e:
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("!!! OpenCL Build Error !!!")
            print(e)
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            # Optionally save the source code that failed to build
            # with open("failed_build_source.cl", "w") as f:
            #     f.write(self.kernel_source)
            # print("Failed source code saved to failed_build_source.cl")
            raise e # Re-raise the exception to stop execution


        # --- Get Kernel Functions ---
        # Use extracted names to get kernel objects
        self.kernels = {}
        for name in self.kernelheaders.keys():
             try:
                 self.kernels[name] = getattr(self.prg, name)
             except AttributeError:
                 print(f"Warning: Kernel function '{name}' found in header but not in compiled program.")
                 # Handle this case - maybe remove from kernelheaders?
        # Ensure all expected kernels are present
        expected_kernels = [
            "evalBondsAndHNeigh_UFF", "evalAngles_UFF", "evalDihedrals_UFF",
            "evalInversions_UFF", "assembleForces_UFF"
        ]
        for k_name in expected_kernels:
            if k_name not in self.kernels:
                raise RuntimeError(f"Expected kernel '{k_name}' not found or failed to compile.")


        # Placeholders - these will be set in realloc
        self.nSystems = 0
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0
        self.ninversions = 0
        self.npbc = 0
        self.nf = 0
        self.i0dih = 0
        self.i0inv = 0
        self.i0ang = 0

        self.buffer_dict = {}
        self.kernel_args = {}
        self.kernel_params = {}
        self.local_size = (nloc,) # Default local size

        # Define precision used (matches OpenCL code)
        self.np_float = np.float32 # Use np.float64 if using double in OpenCL
        self.np_int   = np.int32

    def _extract_kernel_headers(self, source_code):
        """
        Extracts kernel function signatures from OpenCL source code.

        Args:
            source_code (str): The entire OpenCL source code as a string.

        Returns:
            dict: A dictionary where keys are kernel names and values are
                  the full signature string (from '__kernel'/'kernel' up to
                  the closing parenthesis of the arguments).
        """
        headers = {}
        # Regex to find kernel definitions (handles __kernel or kernel, return type, name, and captures up to opening parenthesis)
        # It accounts for potential whitespace and basic pointer usage in return types (like 'kernel void* func(')
        kernel_pattern = re.compile(r'(?:__)?kernel\s+[\w\s\*]+\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*\(')

        current_pos = 0
        while True:
            match = kernel_pattern.search(source_code, current_pos)
            if not match:
                break

            kernel_name = match.group(1)
            kernel_start_pos = match.start()
            open_paren_pos = match.end() - 1 # Position of the opening '('

            # Find the matching closing parenthesis, handling nesting
            paren_level = 1
            close_paren_pos = -1
            search_start = open_paren_pos + 1
            in_block_comment = False
            in_line_comment = False

            for i in range(search_start, len(source_code)):
                char = source_code[i]
                prev_char = source_code[i-1] if i > search_start else ''
                next_char = source_code[i+1] if i < len(source_code)-1 else ''

                # Handle comments - basic version
                if char == '/' and next_char == '*' and not in_line_comment:
                    in_block_comment = True
                elif char == '*' and next_char == '/' and in_block_comment:
                    in_block_comment = False
                    continue # Skip the '/' character of the closing comment marker
                elif char == '/' and prev_char == '/' and not in_block_comment:
                     in_line_comment = True
                elif char == '\n' and in_line_comment:
                     in_line_comment = False

                if in_block_comment or in_line_comment:
                    continue

                # Handle parenthesis nesting
                if char == '(':
                    paren_level += 1
                elif char == ')':
                    paren_level -= 1
                    if paren_level == 0:
                        close_paren_pos = i
                        break

            if close_paren_pos == -1:
                print(f"Warning: Could not find matching ')' for kernel '{kernel_name}' starting at position {kernel_start_pos}. Skipping.")
                current_pos = open_paren_pos + 1 # Move search past the likely broken kernel start
            else:
                # Extract the full signature string
                header_signature = source_code[kernel_start_pos : close_paren_pos + 1]
                headers[kernel_name] = header_signature
                # print(f"Extracted header for: {kernel_name}") # Debug
                current_pos = close_paren_pos + 1 # Continue search after this kernel

        return headers

    # --- Methods for buffer allocation, data upload, kernel setup, execution, download ---
    # (Keep the methods: _check_buf, realloc_buffers, set_a2f_map_size,
    # upload_topology_params, upload_positions, clear_forces, _parse_kernel_header,
    # _init_kernel_params, update_kernel_param, _generate_kernel_args,
    # setup_kernel_args, run_eval_step, get_forces, get_positions, get_energies
    # exactly as they were in the previous version, they rely on self.kernelheaders
    # being populated correctly by the new _extract_kernel_headers method)
    # ... (paste the definitions of these methods here from the previous response) ...
    # --- Start Paste ---
    def _check_buf(self, name, required_size, flags=cl.mem_flags.READ_WRITE):
        """ Helper to create or resize a buffer if needed. """
        current_buf = self.buffer_dict.get(name)
        if current_buf is None or current_buf.size < required_size:
            if current_buf: current_buf.release() # Release old buffer if resizing
            if required_size > 0:
                print(f"Allocating buffer '{name}' with size {required_size} bytes")
                self.buffer_dict[name] = cl.Buffer(self.ctx, flags, size=required_size)
            else:
                 print(f"Warning: Buffer '{name}' has zero size, skipping allocation.")
                 self.buffer_dict[name] = None # Handle zero-size case
        # Ensure the buffer exists if size > 0
        elif required_size == 0 and current_buf is not None:
             # If size is now 0, release the buffer
             print(f"Releasing buffer '{name}' as required size is 0.")
             current_buf.release()
             self.buffer_dict[name] = None
        elif self.buffer_dict.get(name) is None and required_size > 0:
             # This case shouldn't happen if the initial check works, but as safety:
             print(f"Re-Allocating buffer '{name}' with size {required_size} bytes")
             self.buffer_dict[name] = cl.Buffer(self.ctx, flags, size=required_size)


    def realloc_buffers(self, natoms, nbonds, nangles, ndihedrals, ninversions, npbc, nSystems=1):
        """
        Allocates or reallocates all necessary OpenCL buffers for UFF calculations
        based on the provided topology sizes for a single system. Assumes all systems are identical.
        """
        print(f"Reallocating buffers for {nSystems} system(s):")
        print(f"  natoms={natoms}, nbonds={nbonds}, nangles={nangles}, ndihedrals={ndihedrals}, ninversions={ninversions}, npbc={npbc}")

        self.nSystems = nSystems
        self.natoms = natoms
        self.nbonds = nbonds
        self.nangles = nangles
        self.ndihedrals = ndihedrals
        self.ninversions = ninversions
        self.npbc = npbc # Number of PBC shift vectors per system

        # Calculate fint size and offsets based on C++ logic (adjusted for angle/dih/inv only)
        # nf = ndihedrals*4 + ninversions*4 + nangles*3 # Size needed for storing pieces
        # Assuming fint stores float4 {fx,fy,fz,E_contrib}
        f4size = 4 * self.np_float().itemsize
        self.i0dih = 0
        self.i0inv = self.i0dih + self.ndihedrals * 4
        self.i0ang = self.i0inv + self.ninversions * 4
        nf_total_slots = self.i0ang + self.nangles * 3 # Total number of float4 slots in fint
        self.nf = nf_total_slots # Store total slots for potential use

        print(f"  Calculated fint offsets: i0dih={self.i0dih}, i0inv={self.i0inv}, i0ang={self.i0ang}")
        print(f"  Total fint slots (float4): {nf_total_slots}")

        # Define sizes for buffers (per system sizes multiplied by nSystems)
        natoms_tot = nSystems * natoms
        nbonds_tot = nSystems * nbonds
        nangles_tot = nSystems * nangles
        ndihedrals_tot = nSystems * ndihedrals
        ninversions_tot = nSystems * ninversions
        npbc_tot = nSystems * npbc
        nf_buf_size = nSystems * nf_total_slots * f4size # Total byte size for fint buffer

        # --- Allocate Buffers (using _check_buf for potential resizing) ---
        mf = cl.mem_flags
        self._check_buf('apos',       natoms_tot * f4size, mf.READ_WRITE) # Positions
        self._check_buf('fapos',      natoms_tot * f4size, mf.READ_WRITE) # Forces (bonds accumulate here, then others)
        self._check_buf('hneigh',     natoms_tot * 4 * f4size, mf.READ_WRITE) # H-vectors {nx,ny,nz,1/L}

        self._check_buf('neighs',     natoms_tot * 4 * self.np_int().itemsize, mf.READ_ONLY) # Neighbor indices
        self._check_buf('neighCell',  natoms_tot * 4 * self.np_int().itemsize, mf.READ_ONLY) # Neighbor cell indices
        self._check_buf('pbc_shifts', npbc_tot * f4size, mf.READ_ONLY)    # PBC shift vectors

        # UFF Parameters
        self._check_buf('REQs',       natoms_tot * f4size, mf.READ_ONLY)    # Non-bond parameters {R0,E0,Q,_}
        self._check_buf('bonParams',  nbonds_tot * 2 * self.np_float().itemsize, mf.READ_ONLY) # Bond params {K, l0} (float2)
        self._check_buf('angParams1', nangles_tot * f4size, mf.READ_ONLY)   # Angle params {K,c0,c1,c2}
        self._check_buf('angParams2_w',nangles_tot * self.np_float().itemsize, mf.READ_ONLY) # Angle param {c3}
        self._check_buf('dihParams',  ndihedrals_tot * 3 * self.np_float().itemsize, mf.READ_ONLY)# Dihedral params {V,d,n} (float3)
        self._check_buf('invParams',  ninversions_tot * f4size, mf.READ_ONLY)  # Inversion params {K,c0,c1,c2}

        # Precomputed Topology Indices
        self._check_buf('neighBs',    natoms_tot * 4 * self.np_int().itemsize, mf.READ_ONLY) # Map (atom, neigh_slot) -> bond_idx
        self._check_buf('angAtoms',   nangles_tot * 3 * self.np_int().itemsize, mf.READ_ONLY) # Angle atom indices {i,j,k}
        self._check_buf('angNgs',     nangles_tot * 2 * self.np_int().itemsize, mf.READ_ONLY) # Angle h-neigh indices {ji, kj} (int2)
        self._check_buf('dihAtoms',   ndihedrals_tot * 4 * self.np_int().itemsize, mf.READ_ONLY) # Dihedral atom indices {i,j,k,l}
        self._check_buf('dihNgs',     ndihedrals_tot * 3 * self.np_int().itemsize, mf.READ_ONLY) # Dihedral h-neigh indices {ji,kj,lk} (int3)
        self._check_buf('invAtoms',   ninversions_tot * 4 * self.np_int().itemsize, mf.READ_ONLY) # Inversion atom indices {i,j,k,l} (i=center)
        self._check_buf('invNgs',     ninversions_tot * 3 * self.np_int().itemsize, mf.READ_ONLY) # Inversion h-neigh indices {ji,ki,li} (int3)

        # Assembly Map (a2f)
        # Sizes depend on the total number of references, must be provided by host
        # Allocate placeholder - size needs to be determined and buffer filled on host
        self._check_buf('a2f_offsets', natoms_tot * self.np_int().itemsize, mf.READ_ONLY)
        self._check_buf('a2f_counts',  natoms_tot * self.np_int().itemsize, mf.READ_ONLY)
        # Size of a2f_indices = sum(a2f_counts) * itemsize
        # self._check_buf('a2f_indices', total_a2f_refs * self.np_int().itemsize, mf.READ_ONLY)
        # ---> Allocation of a2f_indices requires knowing total_a2f_refs beforehand <---
        # ---> Add a separate method or argument to set this size <---
        self.total_a2f_refs = 0 # Placeholder, needs to be set

        # Temporary Force Storage
        self._check_buf('fint',       nf_buf_size, mf.READ_WRITE)       # Stores angle/dih/inv force pieces

        # Optional Energy Buffers
        self._check_buf('Eb_contrib', natoms_tot * self.np_float().itemsize, mf.READ_WRITE) # Per-atom bond energy
        self._check_buf('Ea_contrib', nangles_tot * self.np_float().itemsize, mf.READ_WRITE)# Per-angle energy
        self._check_buf('Ed_contrib', ndihedrals_tot * self.np_float().itemsize, mf.READ_WRITE)# Per-dihedral energy
        self._check_buf('Ei_contrib', ninversions_tot * self.np_float().itemsize, mf.READ_WRITE)# Per-inversion energy

        print("Buffer allocation/resizing complete.")
        self.setup_kernel_args() # Setup args after buffers are ready

    def set_a2f_map_size(self, total_a2f_refs_per_system):
        """ Sets the size for the a2f_indices buffer. Must be called after realloc_buffers. """
        self.total_a2f_refs = self.nSystems * total_a2f_refs_per_system
        print(f"Setting a2f_indices buffer size for {self.total_a2f_refs} total references.")
        self._check_buf('a2f_indices', self.total_a2f_refs * self.np_int().itemsize, cl.mem_flags.READ_ONLY)


    def upload_topology_params(self, uff_data, iSys=0):
        """
        Uploads static topology and parameter data for one system to the GPU buffers.
        'uff_data' should be an object or dictionary containing NumPy arrays for:
        neighs, neighCell, pbc_shifts, REQs, bonParams, angParams1, angParams2_w,
        dihParams, invParams, neighBs, angAtoms, angNgs, dihAtoms, dihNgs,
        invAtoms, invNgs, a2f_offsets, a2f_counts, a2f_indices.
        """
        print(f"Uploading topology and parameters for system {iSys}...")

        # Calculate byte offsets for this system
        off_natoms_i4  = iSys * self.natoms * 4 * self.np_int().itemsize
        off_natoms_f4  = iSys * self.natoms * 4 * self.np_float().itemsize
        off_nbonds_f2  = iSys * self.nbonds * 2 * self.np_float().itemsize
        off_nangles_f4 = iSys * self.nangles * 4 * self.np_float().itemsize
        off_nangles_f1 = iSys * self.nangles * self.np_float().itemsize
        off_nangles_i3 = iSys * self.nangles * 3 * self.np_int().itemsize
        off_nangles_i2 = iSys * self.nangles * 2 * self.np_int().itemsize
        off_ndihedrals_f3 = iSys * self.ndihedrals * 3 * self.np_float().itemsize
        off_ndihedrals_i4 = iSys * self.ndihedrals * 4 * self.np_int().itemsize
        off_ndihedrals_i3 = iSys * self.ndihedrals * 3 * self.np_int().itemsize
        off_ninversions_f4= iSys * self.ninversions * 4 * self.np_float().itemsize
        off_ninversions_i4= iSys * self.ninversions * 4 * self.np_int().itemsize
        off_ninversions_i3= iSys * self.ninversions * 3 * self.np_int().itemsize
        off_npbc_f4    = iSys * self.npbc * 4 * self.np_float().itemsize
        off_a2f_offs   = iSys * self.natoms * self.np_int().itemsize
        off_a2f_counts = iSys * self.natoms * self.np_int().itemsize
        # Offset for a2f_indices needs the count from previous systems
        a2f_refs_before = iSys * (self.total_a2f_refs // self.nSystems) if self.nSystems > 0 else 0 # Assumes equal systems
        off_a2f_indices = a2f_refs_before * self.np_int().itemsize

        # Helper to upload if data exists
        def _upload_if_present(key, buf_name, offset, dtype):
            data = getattr(uff_data, key, None)
            buf = self.buffer_dict.get(buf_name)
            if buf is None:
                # print(f"  Skipping upload for non-allocated buffer: {buf_name}")
                return
            if data is not None:
                if not isinstance(data, np.ndarray): data = np.array(data, dtype=dtype)
                if data.size * data.dtype.itemsize > 0: # Check if array has size > 0 bytes
                    # print(f"  Uploading {key} to {buf_name} at offset {offset}, size {data.nbytes}") # Debug
                    cl.enqueue_copy(self.queue, buf, data.astype(dtype), device_offset=offset, is_blocking=False)
                else: print(f"  Skipping upload for empty array: {key}")
            else: print(f"  Warning: Data key '{key}' not found in uff_data.")

        # Upload data
        _upload_if_present('neighs',    'neighs',    off_natoms_i4, self.np_int)
        _upload_if_present('neighCell', 'neighCell', off_natoms_i4, self.np_int)
        _upload_if_present('pbc_shifts','pbc_shifts',off_npbc_f4,   self.np_float)
        _upload_if_present('REQs',      'REQs',      off_natoms_f4, self.np_float)
        _upload_if_present('bonParams', 'bonParams', off_nbonds_f2, self.np_float)
        _upload_if_present('angParams1','angParams1',off_nangles_f4,self.np_float)
        _upload_if_present('angParams2_w','angParams2_w',off_nangles_f1,self.np_float)
        _upload_if_present('dihParams', 'dihParams', off_ndihedrals_f3, self.np_float)
        _upload_if_present('invParams', 'invParams', off_ninversions_f4, self.np_float)

        _upload_if_present('neighBs',   'neighBs',   off_natoms_i4, self.np_int)
        _upload_if_present('angAtoms',  'angAtoms',  off_nangles_i3, self.np_int)
        _upload_if_present('angNgs',    'angNgs',    off_nangles_i2, self.np_int)
        _upload_if_present('dihAtoms',  'dihAtoms',  off_ndihedrals_i4, self.np_int)
        _upload_if_present('dihNgs',    'dihNgs',    off_ndihedrals_i3, self.np_int)
        _upload_if_present('invAtoms',  'invAtoms',  off_ninversions_i4, self.np_int)
        _upload_if_present('invNgs',    'invNgs',    off_ninversions_i3, self.np_int)

        _upload_if_present('a2f_offsets','a2f_offsets', off_a2f_offs, self.np_int)
        _upload_if_present('a2f_counts', 'a2f_counts',  off_a2f_counts, self.np_int)
        # Only upload a2f_indices if the buffer was allocated
        if self.buffer_dict.get('a2f_indices'):
             _upload_if_present('a2f_indices','a2f_indices', off_a2f_indices, self.np_int)

        self.queue.finish() # Wait for uploads to complete
        print(f"Finished uploading topology and parameters for system {iSys}.")

    def upload_positions(self, apos_all, iSys=None):
        """ Uploads atom positions for one or all systems. """
        buf = self.buffer_dict.get('apos')
        if buf is None: raise RuntimeError("Position buffer 'apos' not allocated.")

        if iSys is None: # Upload all systems
            expected_shape = (self.nSystems, self.natoms, 4)
            if not isinstance(apos_all, np.ndarray): apos_all = np.array(apos_all, dtype=self.np_float)
            if apos_all.shape != expected_shape:
                 raise ValueError(f"Shape mismatch for apos_all. Expected {expected_shape}, got {apos_all.shape}")
            cl.enqueue_copy(self.queue, buf, apos_all.astype(self.np_float))
        else: # Upload single system
             expected_shape = (self.natoms, 4)
             if not isinstance(apos_all, np.ndarray): apos_all = np.array(apos_all, dtype=self.np_float)
             if apos_all.shape != expected_shape:
                 raise ValueError(f"Shape mismatch for single system apos. Expected {expected_shape}, got {apos_all.shape}")
             offset = iSys * self.natoms * 4 * self.np_float().itemsize
             cl.enqueue_copy(self.queue, buf, apos_all.astype(self.np_float), device_offset=offset)

    def clear_forces(self):
        """ Sets the fapos buffer to zero using enqueue_fill_buffer. """
        buf = self.buffer_dict.get('fapos')
        if buf:
            print("Clearing force buffer (fapos)...")
            try:
                # Pattern needs to match the data type (float4 -> 4 floats)
                pattern = np.zeros(4, dtype=self.np_float) # Or just np.float32(0)? Check docs
                self.queue.enqueue_fill_buffer(buf, pattern, 0, buf.size) # Fill with zeros
                self.queue.finish()
            except Exception as e:
                # Fallback if fill_buffer fails or is not supported well
                print(f"enqueue_fill_buffer failed ({e}), falling back to enqueue_copy.")
                zero_forces = np.zeros((self.nSystems * self.natoms, 4), dtype=self.np_float)
                cl.enqueue_copy(self.queue, buf, zero_forces)
                self.queue.finish()
        else:
             print("Warning: Force buffer 'fapos' not allocated, cannot clear.")


    def _parse_kernel_header(self, header_string):
        """ Parses kernel header to extract arguments (names and types: 0=buffer, 1=scalar). """
        # --- This method remains the same as in the previous version ---
        param_block = header_string[header_string.find('(') + 1:]
        end_idx = param_block.rfind(')')
        if end_idx != -1: param_block = param_block[:end_idx]
        param_lines = [line.strip() for line in param_block.split('\n')]
        param_lines = [line for line in param_lines if line and not line.strip().startswith('//')]

        args = []
        for line in param_lines:
            if '//' in line: line = line.split('//')[0].strip() # Remove inline comments
            if not line: continue
            if line.endswith(','): line = line[:-1].strip()

            parts = line.split()
            if not parts: continue

            if '__global' in line:
                name = parts[-1].replace('*', '').strip(',;')
                args.append((name, 0)) # Type 0 for buffer
            elif 'const' in line:
                # Find the last word, assume it's the name
                name = parts[-1].strip(',;')
                # Basic check if it might be a type instead (e.g., "const int")
                if name.lower() in ['int', 'float', 'int4', 'float4', 'float3', 'float2', 'int2', 'int3']:
                     # This case likely means the name was on the previous part if split weirdly
                     # Or it's a very simple const type definition (less common for kernels)
                     # Let's assume name is the last word for simplicity here.
                     pass # Keep the assumption
                args.append((name, 1)) # Type 1 for scalar/const struct
        return args

    def _init_kernel_params(self):
        """ Initialize dictionary of standard SCALAR kernel parameters. """
        # --- This method remains the same as in the previous version ---
        # Values that are constant or depend only on topology/offsets
        self.kernel_params = {
            # Topology sizes per system
            'natoms':      self.np_int(self.natoms),
            'nangles':     self.np_int(self.nangles),
            'ndihedrals':  self.np_int(self.ndihedrals),
            'ninversions': self.np_int(self.ninversions),
            'npbc':        self.np_int(self.npbc),
            # Offsets in fint
            'i0dih':       self.np_int(self.i0dih),
            'i0inv':       self.np_int(self.i0inv),
            'i0ang':       self.np_int(self.i0ang),
            # 'i0bon':       self.np_int(0), # Unused if bonds write to fapos
            # Common scalar simulation parameters (can be updated later)
            'bSubtractBondNonBond': self.np_int(0), # Default: off
            'bSubtractAngleNonBond':self.np_int(0), # Default: off
            'SubNBTorsionFactor':   self.np_float(0.0), # Default: off
            'Rdamp':                self.np_float(0.1), # Example default
            'FmaxNonBonded':        self.np_float(10.0),# Example default
        }

    def update_kernel_param(self, name, value):
        """ Update a scalar kernel parameter. """
        # --- This method remains the same as in the previous version ---
        if name in self.kernel_params:
            # Update with correct numpy type
            current_type = type(self.kernel_params[name])
            self.kernel_params[name] = current_type(value)
            # Re-generate kernel args that might use this parameter
            self.setup_kernel_args() # Regenerate all args for simplicity
            print(f"Updated kernel parameter '{name}' to {value}")
        else:
            print(f"Warning: Kernel parameter '{name}' not found.")


    def _generate_kernel_args(self, kernel_name):
        """ Generate argument list for a kernel using buffer_dict and kernel_params. """
        # --- This method remains the same as in the previous version ---
        if not hasattr(self, 'kernel_params'): self._init_kernel_params()

        header = self.kernelheaders.get(kernel_name)
        if not header: raise ValueError(f"Header not found for kernel: {kernel_name}")

        arg_specs = self._parse_kernel_header(header)
        args = []
        # print(f"Generating args for {kernel_name}: {arg_specs}") # Debug print
        for name, typ in arg_specs:
            if typ == 0: # Buffer
                buf = self.buffer_dict.get(name)
                if buf is None and name not in ['Eb_contrib', 'Ea_contrib', 'Ed_contrib', 'Ei_contrib']: # Optional energy bufs
                    raise ValueError(f"Buffer '{name}' not allocated or is None for kernel '{kernel_name}'")
                elif buf is None and name in ['Eb_contrib', 'Ea_contrib', 'Ed_contrib', 'Ei_contrib']:
                    # Pass a null pointer if buffer is optional and not allocated
                    # This requires the kernel to handle null pointers for these optional args
                    # PyOpenCL might require explicit handling for None -> cl_mem(0)
                    # Easier: Ensure these buffers are always allocated (even if size 0?) or require them?
                    # For now, let's assume they must be allocated if the kernel expects them.
                    # --> Modify _check_buf to handle size 0 better or remove optional aspect.
                    # --> Let's assume they *are* allocated (even if size 0) or raise error.
                    # --> Correction: Check if kernel actually *uses* the buffer name before erroring
                    if name in header: # Simple check if the name appears in the header text
                         raise ValueError(f"Optional buffer '{name}' is None but required by kernel '{kernel_name}' header.")
                    else:
                         # If name isn't even in the header, skip (unlikely scenario)
                         continue
                args.append(buf)
            elif typ == 1: # Scalar/Const
                param = self.kernel_params.get(name)
                if param is None: raise ValueError(f"Scalar parameter '{name}' not set for kernel '{kernel_name}'")
                args.append(param)
        # print(f"Args for {kernel_name}: {args}") # Debug print
        return args

    def setup_kernel_args(self):
        """ Prepares argument lists for all defined UFF kernels. """
        # --- This method remains the same as in the previous version ---
        self._init_kernel_params() # Ensure params are initialized
        self.kernel_args = {}
        for name in self.kernelheaders.keys():
            # Check if the kernel actually exists in the compiled program
            if name in self.kernels:
                try:
                    self.kernel_args[name] = self._generate_kernel_args(name)
                    print(f"Successfully generated arguments for kernel: {name}")
                except ValueError as e:
                    print(f"Error generating arguments for kernel {name}: {e}")
                    # Decide whether to raise the error or continue
                    # raise e # Option: Stop if any kernel fails setup
            else:
                print(f"Skipping argument generation for kernel '{name}' (not found in compiled program).")
        print("Kernel arguments setup complete.")


    def run_eval_step(self, bClearForce=True):
        """ Executes one step of UFF evaluation kernels. """
        # --- This method remains the same as in the previous version ---
        if not self.kernel_args:
             print("Kernel arguments not set up. Call setup_kernel_args() first.")
             return

        # Ensure kernel_args were generated for all expected kernels
        expected_kargs = [
             "evalBondsAndHNeigh_UFF", "evalAngles_UFF", "evalDihedrals_UFF",
             "evalInversions_UFF", "assembleForces_UFF"
        ]
        for k_name in expected_kargs:
            if k_name not in self.kernel_args:
                 raise RuntimeError(f"Kernel arguments for '{k_name}' are missing. Cannot run evaluation step.")

        natoms_glob = clu.roundup_global_size(self.natoms * self.nSystems, self.nloc)
        nangles_glob = clu.roundup_global_size(self.nangles * self.nSystems, self.nloc)
        ndihedrals_glob = clu.roundup_global_size(self.ndihedrals * self.nSystems, self.nloc)
        ninversions_glob = clu.roundup_global_size(self.ninversions * self.nSystems, self.nloc)

        gs_bonds = (natoms_glob,)       # Atom-centric for bonds/hneigh
        gs_angs = (nangles_glob,)      # Interaction-centric
        gs_dihs = (ndihedrals_glob,)   # Interaction-centric
        gs_invs = (ninversions_glob,)  # Interaction-centric
        gs_assem = (natoms_glob,)      # Atom-centric

        ls = (self.nloc,) # 1D local size

        q = self.queue

        if bClearForce:
            self.clear_forces() # Zero out fapos before accumulating bond forces

        # 1. Bonds and H-Neigh (accumulates bond force into fapos)
        self.kernels['evalBondsAndHNeigh_UFF'](q, gs_bonds, ls, *self.kernel_args['evalBondsAndHNeigh_UFF'])

        # 2. Angles (writes pieces to fint)
        self.kernels['evalAngles_UFF'](q, gs_angs, ls, *self.kernel_args['evalAngles_UFF'])

        # 3. Dihedrals (writes pieces to fint)
        self.kernels['evalDihedrals_UFF'](q, gs_dihs, ls, *self.kernel_args['evalDihedrals_UFF'])

        # 4. Inversions (writes pieces to fint)
        self.kernels['evalInversions_UFF'](q, gs_invs, ls, *self.kernel_args['evalInversions_UFF'])

        # 5. Assemble (adds fint pieces to fapos which contains bond forces)
        self.kernels['assembleForces_UFF'](q, gs_assem, ls, *self.kernel_args['assembleForces_UFF'])

        # q.finish() # Optional: wait for completion here if needed immediately

    def get_forces(self, iSys=None):
        """ Downloads forces for one or all systems. """
        # --- This method remains the same as in the previous version ---
        buf = self.buffer_dict.get('fapos')
        if buf is None: return None
        fapos_np = np.empty((self.nSystems, self.natoms, 4), dtype=self.np_float)
        cl.enqueue_copy(self.queue, fapos_np, buf).wait()
        if iSys is None:
            return fapos_np
        else:
            return fapos_np[iSys]

    def get_positions(self, iSys=None):
        """ Downloads positions for one or all systems. """
        # --- This method remains the same as in the previous version ---
        buf = self.buffer_dict.get('apos')
        if buf is None: return None
        apos_np = np.empty((self.nSystems, self.natoms, 4), dtype=self.np_float)
        cl.enqueue_copy(self.queue, apos_np, buf).wait()
        if iSys is None:
            return apos_np
        else:
            return apos_np[iSys]

    def get_energies(self):
        """ Downloads optional energy contribution buffers. """
        # --- This method remains the same as in the previous version ---
        energies = {}
        buf_map = {
            'Eb_contrib': (self.natoms, 'Eb_per_atom'),
            'Ea_contrib': (self.nangles, 'Ea_per_angle'),
            'Ed_contrib': (self.ndihedrals, 'Ed_per_dihedral'),
            'Ei_contrib': (self.ninversions, 'Ei_per_inversion'),
        }
        for buf_name, (count_per_sys, out_key) in buf_map.items():
            buf = self.buffer_dict.get(buf_name)
            if buf:
                 total_count = self.nSystems * count_per_sys
                 if total_count > 0:
                     e_np = np.empty(total_count, dtype=self.np_float)
                     cl.enqueue_copy(self.queue, e_np, buf).wait()
                     energies[out_key] = e_np.reshape(self.nSystems, count_per_sys)
                 else:
                     energies[out_key] = np.zeros((self.nSystems, 0), dtype=self.np_float) # Handle zero counts

        return energies
    # --- End Paste ---


# ======================================================
# Example Usage (Remains the same as previous version)
# ======================================================
if __name__ == "__main__":
    # --- 1. Setup Phase ---
    # Create UFF_CL instance
    # You might need to provide the kernel path explicitly if auto-detection fails
    # uff_cl = UFF_CL(nloc=64, kernel_path="/path/to/your/uff_kernels.cl")
    uff_cl = UFF_CL(nloc=64)

    # Define system parameters (replace with actual data loading)
    nSystems    = 1
    natoms      = 100
    nbonds      = 99
    nangles     = 150
    ndihedrals  = 200
    ninversions = 10
    npbc        = 27 # Number of shift vectors (e.g., for 3x3x3 cell images, ignoring (0,0,0))

    # Allocate GPU buffers
    uff_cl.realloc_buffers(natoms, nbonds, nangles, ndihedrals, ninversions, npbc, nSystems)

    # --- Host Data Preparation (CRUCIAL) ---
    print("TODO: Prepare host data (NumPy arrays) for topology and parameters.")
    # Using the same MockUFFData as before for demonstration
    class MockUFFData:
        def __init__(self, cl_instance):
            nSys = cl_instance.nSystems
            na = cl_instance.natoms
            nb = cl_instance.nbonds
            ng = cl_instance.nangles
            nd = cl_instance.ndihedrals
            ni = cl_instance.ninversions
            npbc_ = cl_instance.npbc
            nf = cl_instance.nf # Total fint slots

            self.natoms = na; self.nbonds = nb; self.nangles = ng
            self.ndihedrals = nd; self.ninversions = ni; self.npbc = npbc_

            self.apos       = np.random.rand(na, 4).astype(cl_instance.np_float) * 10; self.apos[:,3]=1.0
            self.neighs     = np.random.randint(-1, na, size=(na, 4), dtype=cl_instance.np_int)
            self.neighCell  = np.zeros((na, 4), dtype=cl_instance.np_int)
            self.pbc_shifts = np.zeros((npbc_, 4), dtype=cl_instance.np_float)
            self.REQs       = np.random.rand(na, 4).astype(cl_instance.np_float)
            self.bonParams  = np.random.rand(nb, 2).astype(cl_instance.np_float)
            self.angParams1 = np.random.rand(ng, 4).astype(cl_instance.np_float)
            self.angParams2_w= np.random.rand(ng).astype(cl_instance.np_float)
            self.dihParams  = np.random.rand(nd, 3).astype(cl_instance.np_float)
            self.invParams  = np.random.rand(ni, 4).astype(cl_instance.np_float)
            self.neighBs    = np.random.randint(-1, nb, size=(na, 4), dtype=cl_instance.np_int)
            self.angAtoms   = np.random.randint(0, na, size=(ng * 3), dtype=cl_instance.np_int)
            self.angNgs     = np.zeros((ng, 2), dtype=cl_instance.np_int) # Placeholder
            self.dihAtoms   = np.random.randint(0, na, size=(nd * 4), dtype=cl_instance.np_int)
            self.dihNgs     = np.zeros((nd * 3), dtype=cl_instance.np_int) # Placeholder
            self.invAtoms   = np.random.randint(0, na, size=(ni * 4), dtype=cl_instance.np_int)
            self.invNgs     = np.zeros((ni * 3), dtype=cl_instance.np_int) # Placeholder
            # Ensure topology indices in *Ngs point within valid hneigh range [0, natoms*4 - 1]
            # This needs careful generation based on actual topology
            max_h_idx = na*4
            if ng > 0: self.angNgs = np.random.randint(0, max_h_idx, size=(ng, 2), dtype=cl_instance.np_int)
            if nd > 0: self.dihNgs = np.random.randint(0, max_h_idx, size=(nd * 3), dtype=cl_instance.np_int) # Flattened
            if ni > 0: self.invNgs = np.random.randint(0, max_h_idx, size=(ni * 3), dtype=cl_instance.np_int) # Flattened

            # Assembly Map (a2f) - Example structure
            self.a2f_counts = np.random.randint(1, 10, size=na, dtype=cl_instance.np_int)
            self.a2f_offsets= np.cumsum(np.concatenate(([0], self.a2f_counts[:-1]))).astype(cl_instance.np_int)
            self.total_a2f_refs = int(self.a2f_counts.sum()) # Ensure it's standard int
            # Indices must point within the valid fint range [0, nf - 1]
            # nf = i0ang + nangles*3
            max_fint_idx = cl_instance.i0ang + ng*3
            if self.total_a2f_refs > 0 and max_fint_idx > 0:
                 self.a2f_indices= np.random.randint(0, max_fint_idx, size=self.total_a2f_refs, dtype=cl_instance.np_int)
            elif self.total_a2f_refs > 0 and max_fint_idx <= 0:
                 print("Warning: Cannot generate a2f_indices as max_fint_idx <= 0")
                 self.a2f_indices = np.array([], dtype=cl_instance.np_int)
            else:
                 self.a2f_indices= np.array([], dtype=cl_instance.np_int)

    # Create mock data
    mock_data = MockUFFData(uff_cl)
    if mock_data.total_a2f_refs <= 0:
        print("Warning: Mock data generated 0 total a2f references. Assembly map will be empty.")

    # Set a2f map size (important!)
    uff_cl.set_a2f_map_size(mock_data.total_a2f_refs) # Use refs per system

    # --- 2. Upload Static Data ---
    uff_cl.upload_topology_params(mock_data, iSys=0)

    # --- 3. Simulation/Evaluation Loop ---
    apos_host = mock_data.apos
    uff_cl.upload_positions(apos_host, iSys=0)
    uff_cl.update_kernel_param('bSubtractBondNonBond', 0) # Turn off NB subtraction

    print("\nRunning one UFF evaluation step...")
    start_time = time.time()
    uff_cl.run_eval_step(bClearForce=True)
    uff_cl.queue.finish()
    end_time = time.time()
    print(f"UFF step execution time: {end_time - start_time:.6f} seconds")

    # --- 4. Download Results ---
    print("\nDownloading results...")
    final_forces = uff_cl.get_forces(iSys=0)
    energies = uff_cl.get_energies()

    if final_forces is not None:
        print("Forces (first 5 atoms):\n", final_forces[:5])
    if energies:
        print("Energies:")
        for k, v in energies.items():
            if v.size > 0:
                 print(f"  {k}: shape={v.shape}, sum={np.sum(v)}")
            else:
                 print(f"  {k}: shape={v.shape} (Empty)")


    print("\nPyOpenCL UFF Example Finished.")