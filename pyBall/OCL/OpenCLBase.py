"""
OpenCLBase: PyOpenCL Utility Base Class
========================================

Purpose:
--------
OpenCLBase is a base class for PyOpenCL-accelerated computations in FireCore.
It provides common infrastructure for:
- OpenCL device selection and context management
- Kernel program loading and caching
- GPU buffer allocation and reuse
- Data transfer between host and device
- Kernel launch utilities

Major Functionality:
-------------------
1. Device Management
   - Select GPU device by vendor or index
   - Create OpenCL context and command queue
   - Print device information for debugging

2. Kernel Program Management
   - Load and compile OpenCL kernel files (.cl)
   - Cache compiled programs to avoid recompilation
   - Extract kernel headers for introspection
   - Support build options for conditional compilation

3. Buffer Management
   - try_make_buffers(): Allocate buffers with size checking
   - try_make_buff(): Allocate single buffer with size checking
   - toGPU_(): Copy host data to GPU buffer
   - fromGPU_(): Copy GPU buffer to host data
   - buffer_dict: Centralized buffer storage for reuse

4. Utility Functions
   - _flat32(): Convert array to C-contiguous float32
   - _int32(): Convert array to C-contiguous int32
   - extract_kernel_headers(): Parse kernel signatures from source

Design Pattern:
--------------
Subclasses inherit from OpenCLBase and implement domain-specific kernels
(e.g., AFM, Grid projection, Hubbard solver). The base class handles all
boilerplate for device management, buffer allocation, and data transfer.

PyOpenCL Resource Management and Optimization Policy
==================================================

Core Principle:
--------------
Kernel calls are expected to be invoked many times with similar or only slightly changed input.
Python operations (memory allocation, dict creation, kernel compilation) are costly overhead that must be minimized.

Optimization Strategies:
-----------------------

1. Kernel Caching (load_program)
   - Compile kernels once during initialization, not per-call
   - Cache compiled program in self.prg
   - Skip compilation if program already loaded
   - Print messages only when actually compiling, not when using cached program

2. Persistent Buffer Management (try_make_buffers, realloc_* methods)
   - Allocate buffers once and reuse across calls
   - Use try_make_buffers() which checks if buffer exists and has correct size
   - Only reallocate if size changes or buffer doesn't exist
   - Store buffers in self.buffer_dict for centralized management

3. bTryAllocate Guards (calling functions)
   - Add bTryAllocate=True parameter to run functions (default True for safety)
   - Guard dict creation and buffer allocation with "if bTryAllocate:"
   - When False, skip dict creation and try_make_buffers call entirely
   - Allows zero-allocation hot paths for repeated calls with same buffer sizes
   - Pattern: guard in calling function, NOT in realloc_* helper function

4. Separation of Initialization vs Runtime
   - Load kernels during __init__ or load_basis_sto, not in per-call methods
   - Allocate buffers during first call or when sizes change, not every call
   - Use bAlloc/bTryAllocate flags to control when reallocation happens

Example Pattern:
--------------
def my_kernel_function(self, data, bTryAllocate=True):
    # Guard buffer allocation - skip dict creation when not needed
    if bTryAllocate:
        buffs = {"input": sz, "output": sz}
        self.try_make_buffers(buffs, suffix="_buff")
    
    # Data transfer (always needed)
    self.toGPU_(self.input_buff, data)
    
    # Kernel launch
    self.prg.my_kernel(self.queue, gs, ls, self.input_buff, self.output_buff)
    
    # Read back (always needed)
    return self.fromGPU_(self.output_buff)

Usage:
------
# First call: allocates buffers
result = my_kernel_function(data, bTryAllocate=True)

# Subsequent calls with same buffer sizes: skips allocation
result = my_kernel_function(new_data, bTryAllocate=False)
"""

import os
import re
import numpy as np
import pyopencl as cl
from . import clUtils as clu


def print_devices(platforms=None ):
    if platforms is None:
        platforms = cl.get_platforms()
    for i, platform in enumerate(platforms):
        print(f"Platform {i}: {platform.name}")
        devices = platform.get_devices()
        for j, device in enumerate(devices):
            print(f"  Device {j}: {device.name}")

def select_device(platforms=None, preferred_vendor='nvidia', bPrint=False, device_index=0, return_queue=False):
    if platforms is None:
        platforms = cl.get_platforms()
    if bPrint:
        print_devices(platforms)
    preferred_devices = []
    for platform in platforms:
        for device in platform.get_devices():
            name = device.name.lower()
            vendor = device.vendor.lower()
            if (preferred_vendor is not None) and ((preferred_vendor.lower() in name) or (preferred_vendor.lower() in vendor)):
                preferred_devices.append((platform, device))
    if preferred_devices:
        platform, device = preferred_devices[0]
        ctx = cl.Context([device])
        if bPrint:
            print(f"Selected {preferred_vendor} device: {device.name}")
    else:
        devices = []
        for platform in platforms:
            for device in platform.get_devices():
                devices.append(device)
        if len(devices) > 0:
            i = max(0, min(int(device_index), len(devices) - 1))
            device = devices[i]
            ctx = cl.Context([device])
            if bPrint:
                print(f"Selected fallback device {i}: {device.name}")
        else:
            if bPrint:
                print(f"Selected default device {device_index}")
            ctx = cl.create_some_context(answers=[device_index])
            device = ctx.devices[0]
    queue = cl.CommandQueue(ctx)
    if bPrint and len(ctx.devices) > 0:
        clu.get_cl_info(ctx.devices[0])
    if return_queue:
        return ctx, queue
    return ctx

class OpenCLBase:
    """
    Base class for OpenCL applications providing common functionality.
    
    This class handles:
    - OpenCL context and queue initialization
    - Kernel loading and compilation
    - Kernel header extraction
    - Buffer management
    - Common utility functions for OpenCL operations
    """
    
    @staticmethod
    def _flat32(arr):
        return np.require(arr, dtype=np.float32, requirements=('C',)) # .view(np.float32)   .ravel()

    @staticmethod
    def _int32(arr):
        return np.require(arr, dtype=np.int32, requirements=('C',)) # .view(np.int32).ravel()


    def __init__(self, nloc=32, device_index=0, preferred_vendor='nvidia', bPrint=True):
        """
        Initialize the OpenCL environment.
        
        Args:
            nloc (int): Local work group size
            device_index (int): Index of the device to use (default: 0)
        """
        self.nloc = nloc
        self.ctx, self.queue = select_device(preferred_vendor=preferred_vendor, bPrint=bPrint, device_index=device_index, return_queue=True)
            
        self.buffer_dict = {}
        self.kernelheaders = {}
        self.prg = None
    
    def load_program(self, kernel_path=None, rel_path=None, base_path=None, bPrint=False, bMakeHeaders=True, build_options=None):
        """
        Load and compile an OpenCL program.
        Only compiles if program is not already loaded.

        Args:
            kernel_path (str): Absolute path to the kernel file
            rel_path (str): Relative path to the kernel file from base_path
            base_path (str): Base directory path (defaults to this file's directory)
            bPrint (bool): Print debug info
            bMakeHeaders (bool): Extract kernel headers
            build_options (str): OpenCL build options

        Returns:
            bool: True if successful, False otherwise
        """
        # Return early if program is already loaded
        if self.prg is not None:
            print(f"load_program({kernel_path!r} - already loaded)")
            return True

        if kernel_path is None and rel_path is not None:
            if base_path is None:
                base_path = os.path.dirname(os.path.abspath(__file__))
            kernel_path = os.path.abspath(os.path.join(base_path, rel_path))
        if not os.path.exists(kernel_path):
            print(f"OpenCLBase::load_program() ERROR: Kernel file not found at: {kernel_path}")
            return False
        with open(kernel_path, 'r') as f:
            try:
                kernel_source = f.read()
                print(f"[OpenCLBase] Starting compilation for {kernel_path}...")
                if build_options is None:
                    self.prg = cl.Program(self.ctx, kernel_source).build()
                else:
                    self.prg = cl.Program(self.ctx, kernel_source).build(options=build_options)
                print(f"[OpenCLBase] Compilation complete for {kernel_path}")
                # Extract kernel headers automatically
                if bMakeHeaders:
                    self.kernelheaders = self.extract_kernel_headers(kernel_source)
                    if bPrint:
                        for kernel_name, kernel_header in self.kernelheaders.items(): print(f"OpenCLBase::extract_kernel_headers() Kernel name:: {kernel_name} \n {kernel_header}")
                        print(f"Extracted headers for kernels: {list(self.kernelheaders.keys())}")
            except Exception as e:
                print(f"OpenCLBase::load_program() ERROR: Failed to build kernel: {kernel_path}")
                #print(f"Kernel source:\n {kernel_source}")
                print(f"Error: {str(e)}")
                raise e
        if bPrint:
            print(f"OpenCLBase::load_program() Successfully loaded kernel from: {kernel_path}")
        return True
    
    def extract_kernel_headers(self, source_code):
        """
        Extract kernel headers from OpenCL source code.
        
        Args:
            source_code (str): The OpenCL source code as a string
            
        Returns:
            dict: A dictionary where keys are kernel names (without parentheses)
                  and values are the full header string
        """
        headers = {}

        # Robust extraction: for each "__kernel void name(" find the matching ')' of the parameter list
        import re
        for m in re.finditer(r'__kernel\s+void\s+([A-Za-z0-9_]+)\s*\(', source_code):
            kname = m.group(1)
            i0 = m.start()
            # Find the '(' that starts the parameter list
            ip = source_code.find('(', m.end() - 1)
            if ip < 0:
                continue
            level = 1
            i = ip + 1
            while i < len(source_code) and level > 0:
                c = source_code[i]
                if c == '(':
                    level += 1
                elif c == ')':
                    level -= 1
                i += 1
            if level != 0:
                continue
            # Include lines from kernel start up to the line containing the closing ')'
            i1 = i
            # Expand to end of line for readability (optional)
            while i1 < len(source_code) and source_code[i1] != '\n':
                i1 += 1
            headers[kname] = source_code[i0:i1]

        return headers
    
    def create_buffer(self, name, size, flags=cl.mem_flags.READ_WRITE):
        """
        Create an OpenCL buffer and store it in the buffer dictionary.
        
        Args:
            name (str): Name of the buffer
            size (int): Size of the buffer in bytes
            flags (cl.mem_flags): Memory flags for the buffer
            
        Returns:
            cl.Buffer: The created buffer
        """
        if size <= 0:
            raise ValueError(f"Invalid buffer size for {name}: {size}")
            
        buffer = cl.Buffer(self.ctx, flags, size=size)
        self.buffer_dict[name] = buffer
        return buffer
    
    def check_buf(self, name, required_size, flags=cl.mem_flags.READ_WRITE):
        """ Helper to create or resize a buffer if needed. """
        current_buf = self.buffer_dict.get(name)
        if current_buf is None or current_buf.size < required_size:
            if current_buf: current_buf.release() # Release old buffer if resizing
            if required_size > 0:
                print(f"OpenCLBase::check_buf() Allocating buffer '{name}' with size {required_size} bytes")
                self.buffer_dict[name] = cl.Buffer(self.ctx, flags, size=required_size)
            else:
                print(f"OpenCLBase::check_buf() Warning: Buffer '{name}' has zero size, skipping allocation.")
                self.buffer_dict[name] = None # Handle zero-size case
        # Ensure the buffer exists if size > 0
        elif required_size == 0 and current_buf is not None:
            # If size is now 0, release the buffer
            print(f"Releasing buffer '{name}' as required size is 0.")
            current_buf.release()
            self.buffer_dict[name] = None
        elif self.buffer_dict.get(name) is None and required_size > 0:
            # This case shouldn't happen if the initial check works, but as safety:
            print(f"OpenCLBase::check_buf() Re-Allocating buffer '{name}' with size {required_size} bytes")
            self.buffer_dict[name] = cl.Buffer(self.ctx, flags, size=required_size)

    def try_make_buff( self, buff_name, sz):
        if hasattr(self,buff_name): 
            buff = getattr(self, buff_name)
            if not ( buff is None or buff.size != sz ):
                return buff, False
        #print( "try_make_buff(",buff_name,") reallocate to [bytes]: ", sz )
        buff = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, sz )
        setattr(self, buff_name, buff )
        return buff, True

    def try_buff(self, name, names, sz, suffix="_buff" ):
        if name in names: self.try_make_buff( name+suffix , sz)

    def try_make_buffers( self, buffs, suffix="_buff" ):
        for name, sz in buffs.items():
            #final_name = alias if alias else buffname
            buff_name = name + suffix
            # if getattr(self, buff_name, None) is None or getattr(self, buff_name).size != sz:
            #     self.try_make_buff(buff_name, sz)
            buff,new = self.try_make_buff( buff_name, sz)
            # Ensure the buffer is also in the buffer_dict for toGPU/fromGPU
            if new:  # Only add if it was newly created
                self.buffer_dict[buff_name] = buff

    def toGPU_(self, buf, host_data, byte_offset=0 ):
        cl.enqueue_copy(self.queue, buf, host_data, device_offset=byte_offset)

    def fromGPU_(self, buf, host_data=None, byte_offset=0, shape=None, dtype='f4' ):
        if host_data is None: host_data = np.empty(shape, dtype=dtype)
        cl.enqueue_copy(self.queue, host_data, buf, device_offset=byte_offset)
        return host_data
        
    def toGPU(self, buf_name, host_data, byte_offset=0):
        """
        Upload data to a GPU buffer.
        
        Args:
            buf_name (str): Name of the buffer in the buffer dictionary
            host_data (numpy.ndarray): Data to upload
            byte_offset (int): Offset in bytes
        """
        cl.enqueue_copy(self.queue, self.buffer_dict[buf_name], host_data, device_offset=byte_offset)
    
    def fromGPU(self, buf_name, host_data, byte_offset=0):
        """
        Download data from a GPU buffer.
        
        Args:
            buf_name (str): Name of the buffer in the buffer dictionary
            host_data (numpy.ndarray): Array to store the downloaded data
            byte_offset (int): Offset in bytes
        """
        cl.enqueue_copy(self.queue, host_data, self.buffer_dict[buf_name], device_offset=byte_offset)
    
    def download_buf(self, name_or_buf, dtype=None):
        """Minimal helper to download a full OpenCL buffer to a flat numpy array.

        - Accepts either a buffer name (string key into buffer_dict) or a cl.Buffer
        - Infers dtype from name to avoid boilerplate (ints for index-like buffers, float otherwise)
        - Returns a 1D numpy array; caller can reshape as needed
        """
        # Resolve buffer object
        if isinstance(name_or_buf, str):
            buf = self.buffer_dict.get(name_or_buf)
            name = name_or_buf
        else:
            buf = name_or_buf
            name = None
        if buf is None:
            return None

        # Infer dtype if not provided
        if dtype is None:
            int_markers = (
                'Atoms','Ngs','neigh','atype','indices','counts','offsets','Cell','a2f_'
            )
            if isinstance(name_or_buf, str) and any(m in name_or_buf for m in int_markers):
                dtype = np.int32
            else:
                dtype = np.float32

        # Allocate host array and copy
        n = buf.size // np.dtype(dtype).itemsize
        host = np.empty(n, dtype=dtype)
        cl.enqueue_copy(self.queue, host, buf)
        return host
    
    def bufflist(self, names):
        """
        Get a list of buffers by name.
        
        Args:
            names (list): List of buffer names
            
        Returns:
            list: List of buffer objects
        """
        return [self.buffer_dict[name] for name in names]
    
    def roundUpGlobalSize(self, global_size):
        """
        Round up the global work size to a multiple of the local work size.
        
        Args:
            global_size (int): The global work size
            
        Returns:
            int: Rounded up global work size
        """
        return (global_size + self.nloc - 1) // self.nloc * self.nloc
    
    def parse_kernel_header(self, header_string):
        """
        Parse a kernel header to extract buffer and parameter information.
        Improved version that properly handles comments and multi-line declarations.
        """
        # Extract parameter block between the first '(' and its matching ')'.
        # Do NOT use rfind(')') because header_string may accidentally include
        # pieces of the kernel body (e.g., printf format strings containing ')').
        i0 = header_string.find('(')
        if i0 < 0:
            return []
        level = 1
        i1 = i0 + 1
        while i1 < len(header_string) and level > 0:
            c = header_string[i1]
            if c == '(':
                level += 1
            elif c == ')':
                level -= 1
            i1 += 1
        if level != 0:
            raise ValueError("parse_kernel_header: unmatched parentheses in kernel header")
        param_block = header_string[i0 + 1:i1 - 1]
        
        # Split into lines and clean them
        param_lines = []
        for line in param_block.split('\n'):
            line = line.strip()
            # Skip empty lines and full-line comments
            if not line or line.startswith('//'):
                continue
            # Remove inline comments
            if '//' in line:
                line = line.split('//')[0].strip()
            if line:  # If anything remains after cleaning
                param_lines.append(line)
        
        # Join lines and split by commas to get individual parameters
        params = []
        current_param = ''
        for line in param_lines:
            current_param += ' ' + line if current_param else line
            if line.endswith(','):
                params.append(current_param[:-1].strip())  # Remove trailing comma
                current_param = ''
        if current_param:  # Add last parameter if no trailing comma
            params.append(current_param.strip())
        
        # Extract parameter names
        args = []
        for param in params:
            if not param:
                continue
            
            # Handle image parameters (__read_only image*_t)
            if '__read_only' in param and 'image' in param:
                parts = param.split()
                # Find the image type (image3d_t, image2d_t, etc)
                for i, part in enumerate(parts):
                    if part.startswith('image'):
                        param_name = parts[i+1].replace(',', '').strip()
                        args.append((param_name, 0))  # 2 indicates image type
                        break
                continue
                
            # Handle buffer parameters (__global)
            if '__global' in param:
                parts = param.split()
                param_name = parts[-1].replace('*', '').strip()
                args.append((param_name, 0))  # 0 indicates buffer
                continue
                
            # Handle other parameters (constants)
            parts = param.split()
            param_name = parts[-1].replace(',', '').strip()
            args.append((param_name, 1))  # 1 indicates constant
        
        return args
    
    def generate_kernel_args(self, kname, overrides=None, bPrint=False):
        """
        Generate argument list for a kernel based on its header definition.
        
        Args:
            kname (str): Kernel name (without parentheses)
        
        Returns:
            list: List of arguments for the kernel call
        """
        if not hasattr(self, 'kernel_params'):
            raise AttributeError("kernel_params dictionary not initialized")
            
        if kname not in self.kernelheaders:
            print(f"OpenCLBase::generate_kernel_args() Kernel '{kname}' not found in kernel headers")
            print("Available kernels:", list(self.kernelheaders.keys()))
            raise KeyError(f"Kernel '{kname}' not found in kernel headers")
            
        kernel_header = self.kernelheaders[kname]
        args_names    = self.parse_kernel_header(kernel_header)

        if bPrint:
            print(f"OpenCLBase::generate_kernel_args() Kernel '{kname}' header:")
            print(kernel_header)
            print(f"OpenCLBase::generate_kernel_args() Kernel '{kname}' args_names:")
            for i,arg in enumerate(args_names): print("    ", i, arg)

        args = []
        overrides = overrides or {}

        try:
            for aname, typ in args_names:
                if typ == 0:
                    # Some kernels may have scalar args mis-classified as buffers by header parsing.
                    # Prefer buffers when present; otherwise fall back to kernel_params.
                    if aname in self.buffer_dict:
                        args.append(self.buffer_dict[aname])
                    elif aname in overrides:
                        args.append(overrides[aname])
                    else:
                        args.append(self.kernel_params[aname])
                else:
                    if aname in overrides:
                        args.append(overrides[aname])
                    else:
                        args.append(self.kernel_params[aname])
        except KeyError as e:
            print ( "kernel_header ", kernel_header )
            print(f"OpenCLBase::generate_kernel_args() KeyError: {e}")
            raise

        if bPrint:
            print(f"OpenCLBase::generate_kernel_args() Kernel '{kname}' args:")
            for i,arg in enumerate(args): print("    ", i, arg)

        return args

    def preprocess_opencl_source(self, source_path, substitutions=None, output_path=None, bPrint=False):
        """
        Preprocess OpenCL source code with file/function/macro substitutions.
        
        Args:
            source_path (str): Path to the source .cl file
            substitutions (dict): Dictionary with:
                - 'files': dict mapping file markers to file paths
                - 'functions': dict mapping function markers to function names
                - 'macros': dict mapping macro markers to macro definitions
        
        Returns:
            str: Preprocessed source code
        """
        print(f"OpenCLBase::preprocess_opencl_source() source_path: {source_path}")
        if substitutions is None:
            substitutions = {}
            
        with open(source_path, 'r') as f:
            source = f.read()
        
        # Process file inclusions
        if 'files' in substitutions:
            for marker, file_path in substitutions['files'].items():
                if os.path.exists(file_path):
                    with open(file_path, 'r') as f:
                        file_content = f.read()
                    source = source.replace(f"//<<<file {marker}", file_content)
        
        # Process function substitutions
        if 'functions' in substitutions:
            for marker, func_name in substitutions['functions'].items():
                source = source.replace(f"//<<<function {marker}", func_name)
        
        # Process macro substitutions
        if 'macros' in substitutions:
            # Replace only exact sentinel lines to avoid replacing occurrences inside comments/docs
            # We match lines whose stripped content is exactly `//<<<{marker}`
            lines = source.split('\n')
            for idx, line in enumerate(lines):
                stripped = line.strip()
                for marker, macro_def in substitutions['macros'].items():
                    if stripped == f"//<<<{marker}":
                        # Insert macro definition in place of the marker line
                        # Ensure we keep surrounding code structure; macro_def may contain newlines
                        lines[idx] = macro_def
                        break
            source = '\n'.join(lines)

        # if bPrint:
        #     print("================================") 
        #     print(f"preprocess_opencl_source() preprocess_opencl_source saved to :\n", output_path)
        #     print(f"preprocess_opencl_source() preprocess_opencl_source :\n", source)
        #     print("================================") 
        
        # Save to file if output_path specified
        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'w') as f: f.write(source)
            if bPrint: print(f"preprocess_opencl_source() preprocessed source {source_path} saved to: {output_path}")

        return source

    def parse_cl_lib(self, cl_path):
        """
        Parse an OpenCL source file to extract functions and macros marked by //>>> sections.
        Sections should start with one of:
          //>>>function <name and optional signature>
          //>>>macro    <name>

        Returns a dict with:
          { 'functions': { name: body_str }, 'macros': { name: body_str } }
        """
        functions = {}
        macros = {}

        with open(cl_path, 'r') as f:
            content = f.read()

        lines = content.split('\n')
        current_kind = None   # 'function' | 'macro'
        current_name = None
        current_content = []

        for line in lines:
            raw = line
            line = line.strip()

            if line.startswith('//>>>'):
                # Save previous section
                if current_kind and current_name and current_content:
                    body = '\n'.join(current_content)
                    if current_kind == 'function':
                        functions[current_name] = body
                    elif current_kind == 'macro':
                        macros[current_name] = body

                # Start new section
                header = line[len('//>>>'):].strip()
                # Expect format: "function NAME..." or "macro NAME"
                if header.startswith('function'):
                    current_kind = 'function'
                    name_part = header[len('function'):].strip()
                    # Use the token before first '(' as canonical function key
                    fname = name_part.split('(')[0].strip() if '(' in name_part else name_part
                    current_name = fname
                elif header.startswith('macro'):
                    current_kind = 'macro'
                    name_part = header[len('macro'):].strip()
                    # Macro key is the next token
                    mname = name_part.split()[0] if name_part else ''
                    current_name = mname
                else:
                    # Unknown header kind; treat as macro with whole header as name
                    current_kind = 'macro'
                    current_name = header
                current_content = []
            else:
                if current_kind:
                    current_content.append(raw)

        # Save last section
        if current_kind and current_name and current_content:
            body = '\n'.join(current_content)
            if current_kind == 'function':
                functions[current_name] = body
            elif current_kind == 'macro':
                macros[current_name] = body

        return {'functions': functions, 'macros': macros}

    def parse_forces_cl(self, forces_path):
        """
        Parse Forces.cl file to extract functions and macros.
        
        Args:
            forces_path (str): Path to Forces.cl file
            
        Returns:
            dict: Dictionary with 'functions' and 'macros' extracted from Forces.cl
        """
        # Keep for backward compatibility; delegate to generic parser
        return self.parse_cl_lib(forces_path)

    # def load_preprocessed_program(self, source_path, substitutions=None, output_path=None):
    #     """
    #     Load and compile preprocessed OpenCL source with substitutions.
        
    #     Args:
    #         source_path: Path to the original .cl file
    #         substitutions: Dictionary with file/function/macro substitutions
    #         output_path: Optional path to save the preprocessed file
    #     """
    #     if substitutions is None:
    #         substitutions = {}
    #     # Preprocess the source
    #     preprocessed_source = self.preprocess_opencl_source(source_path, substitutions)
    #     # Compile the preprocessed source
    #     self.prg = cl.Program(self.ctx, preprocessed_source).build()
    #     self.kernelheaders = self.extract_kernel_headers(preprocessed_source)
    #     return self.prg