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

def select_device(platforms=None, preferred_vendor='nvidia', bPrint=False, device_index=0):
    if platforms is None:
        platforms = cl.get_platforms()
    if bPrint:
        print_devices(platforms)
    # Try to find preferred vendor device
    preferred_devices = []
    for platform in platforms:
        for device in platform.get_devices():
            if preferred_vendor.lower() in device.name.lower():
                preferred_devices.append((platform, device))
    if preferred_devices:
        platform, device = preferred_devices[0]
        ctx = cl.Context([device])
        if bPrint:
            print(f"Selected {preferred_vendor} device: {device.name}")
    else:
        # Fall back to default behavior
        if bPrint:
            print(f"Selected default device {device_index}")
        ctx = cl.create_some_context(answers=[device_index])
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
    
    def __init__(self, nloc=32, device_index=0):
        """
        Initialize the OpenCL environment.
        
        Args:
            nloc (int): Local work group size
            device_index (int): Index of the device to use (default: 0)
        """
        self.nloc = nloc
        self.ctx = select_device(preferred_vendor='nvidia', bPrint=True)
        clu.get_cl_info(self.ctx.devices[0])
        self.queue = cl.CommandQueue(self.ctx)
            
        self.buffer_dict = {}
        self.kernelheaders = {}
        self.prg = None
    
    def load_program(self, kernel_path=None, rel_path=None, base_path=None, bPrint=False):
        """
        Load and compile an OpenCL program.
        
        Args:
            kernel_path (str): Absolute path to the kernel file
            rel_path (str): Relative path to the kernel file from base_path
            base_path (str): Base directory path (defaults to this file's directory)
            
        Returns:
            bool: True if successful, False otherwise
        """
        if kernel_path is None and rel_path is not None:
            if base_path is None:
                base_path = os.path.dirname(os.path.abspath(__file__))
            kernel_path = os.path.abspath(os.path.join(base_path, rel_path))
        
        if not os.path.exists(kernel_path):
            print(f"OpenCLBase::load_program() ERROR: Kernel file not found at: {kernel_path}")
            return False
            
        with open(kernel_path, 'r') as f:
            kernel_source = f.read()
            self.prg = cl.Program(self.ctx, kernel_source).build()
            # Extract kernel headers automatically
            self.kernelheaders = self.extract_kernel_headers(kernel_source)

            if bPrint:
                for kernel_name, kernel_header in self.kernelheaders.items():
                    print(f"OpenCLBase::extract_kernel_headers() Kernel name:: {kernel_name} \n {kernel_header}")
                
        if bPrint:
            print(f"OpenCLBase::load_program() Successfully loaded kernel from: {kernel_path}")
            print(f"Extracted headers for kernels: {list(self.kernelheaders.keys())}")
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
        
        # Split into lines and process line by line
        lines = source_code.split('\n')
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            # Skip empty lines and commented lines
            if not line or line.startswith('//'):
                i += 1
                continue
                
            # Check for kernel definition
            if line.startswith('__kernel'):
                kernel_start = i
                # Get clean kernel name (split on whitespace and take third element)
                kernel_name = line.split()[2].split('(')[0]  # Remove any trailing parenthesis
                
                # Find opening parenthesis
                while '(' not in line and i < len(lines):
                    i += 1
                    line = lines[i].strip()
                    
                # Find closing parenthesis
                paren_level = 1
                i += 1
                while i < len(lines) and paren_level > 0:
                    line = lines[i].strip()
                    if not line.startswith('//'):  # Skip comment lines
                        paren_level += line.count('(')
                        paren_level -= line.count(')')
                    i += 1
                    
                # Extract full header
                header = '\n'.join(lines[kernel_start:i])
                headers[kernel_name] = header
            else:
                i += 1
                
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
        # Extract parameter block (everything between parentheses)
        param_block = header_string[header_string.find('(') + 1:header_string.rfind(')')]
        
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
            # Handle buffer parameters (__global)
            if '__global' in param:
                parts = param.split()
                param_name = parts[-1].replace('*', '').strip()
                args.append((param_name, 0))
            # Handle direct parameters (const)
            elif 'const' in param:
                parts = param.split()
                if len(parts) >= 3:
                    param_name = parts[2].strip(',;')
                    args.append((param_name, 1))
        
        return args
    
    def generate_kernel_args(self, kname):
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
        args_names   = self.parse_kernel_header(kernel_header)
        
        args = []
        try:
            for aname, typ in args_names:
                if typ == 0:
                    args.append(self.buffer_dict[aname])
                else:
                    args.append(self.kernel_params[aname])
        except KeyError as e:
            print ( "kernel_header ", kernel_header )
            print(f"OpenCLBase::generate_kernel_args() KeyError: {e}")
            raise
                
        return args