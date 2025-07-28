import os
import re
import numpy as np
import pyopencl as cl

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

class OCLSystem:
    """
    A system for managing OpenCL context, programs, kernels, and buffers.
    """
    def __init__(self, nloc=32, device_index=0, preferred_vendor='nvidia'):
        self.nloc = nloc
        self.ctx = select_device(preferred_vendor=preferred_vendor, bPrint=True, device_index=device_index)
        self.queue = cl.CommandQueue(self.ctx)
        self.programs = {}
        self.buffers = {}
        self.kernels = {} # Cache for kernel objects
        self.kernel_headers = {}
        self.kernel_params = {}

    def load_program(self, name, kernel_filepath):
        """
        Load and compile an OpenCL program from a file.
        Args:
            name (str): A unique name for this program.
            kernel_filepath (str): Absolute path to the kernel file.
        """
        print(f"OCLSystem::load_program() Loading kernel '{name}' from: {kernel_filepath}")
        if not os.path.exists(kernel_filepath):
            print(f"OCLSystem::load_program() ERROR: Kernel file not found at: {kernel_filepath}")
            return False
        with open(kernel_filepath, 'r') as f:
            kernel_source = f.read()
            try:
                prg = cl.Program(self.ctx, kernel_source).build()
                self.programs[name] = prg
                self.kernel_headers[name] = self._extract_kernel_headers(kernel_source)
                print(f"OCLSystem::load_program() Successfully loaded kernel '{name}' from: {kernel_filepath}")
                return True
            except cl.LogicError as e:
                print(f"OCLSystem::load_program() OpenCL compilation error for '{name}': {e}")
                return False

    def _extract_kernel_headers(self, source_code):
        """
        Extract kernel headers from OpenCL source code.
        """
        headers = {}
        lines = source_code.split('\n')
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('//'):
                i += 1
                continue
            if line.startswith('__kernel'):
                kernel_name = line.split()[2].split('(')[0]
                
                kernel_signature_lines = []
                # Collect lines until the opening curly brace '{' is found
                while i < len(lines):
                    current_line = lines[i]
                    kernel_signature_lines.append(current_line)
                    if '{' in current_line and not current_line.strip().startswith('//'): # Ensure it's not a commented out brace
                        break
                    i += 1
                
                # Join the collected lines and remove the part after the first '{'
                full_signature_with_body_start = "\n".join(kernel_signature_lines)
                kernel_signature = full_signature_with_body_start.split('{', 1)[0].strip() # Get only the part before the first '{'
                headers[kernel_name] = kernel_signature
            i += 1 # Move to the next line
        return headers

    def create_buffer(self, name, size, flags=cl.mem_flags.READ_WRITE, hostbuf=None):
        """
        Create an OpenCL buffer and store it.
        """
        if size <= 0:
            # print(f"Warning: Buffer '{name}' has zero size, skipping allocation.")
            self.buffers[name] = None # Store None for zero-sized buffers
            return None
        
        # Release old buffer if it exists and needs resizing
        if name in self.buffers and self.buffers[name] is not None:
            if self.buffers[name].size != size:
                print(f"OCLSystem::create_buffer() Re-allocating buffer '{name}' with new size {size} bytes (old size {self.buffers[name].size} bytes)")
                self.buffers[name].release()
                self.buffers[name] = cl.Buffer(self.ctx, flags, size=size, hostbuf=hostbuf)
            # else: buffer already exists with correct size, do nothing
        else:
            print(f"OCLSystem::create_buffer() Allocating buffer '{name}' with size {size} bytes")
            self.buffers[name] = cl.Buffer(self.ctx, flags, size=size, hostbuf=hostbuf)
        return self.buffers[name]

    def toGPU(self, buffer_name, host_data, byte_offset=0):
        """
        Upload data to a GPU buffer.
        """
        if self.buffers.get(buffer_name) is None:
            raise ValueError(f"Buffer '{buffer_name}' does not exist or is zero-sized.")
        cl.enqueue_copy(self.queue, self.buffers[buffer_name], host_data, device_offset=byte_offset)

    def fromGPU(self, buffer_name, host_data, byte_offset=0):
        """
        Download data from a GPU buffer.
        """
        if self.buffers.get(buffer_name) is None:
            raise ValueError(f"Buffer '{buffer_name}' does not exist or is zero-sized.")
        cl.enqueue_copy(self.queue, host_data, self.buffers[buffer_name], device_offset=byte_offset)

    def get_buffer(self, name):
        """
        Get a buffer by name.
        """
        return self.buffers.get(name)

    def round_up_global_size(self, global_size):
        """
        Round up the global work size to a multiple of the local work size.
        """
        return (global_size + self.nloc - 1) // self.nloc * self.nloc

    def _parse_kernel_header(self, header_string):
        """
        Parse a kernel header to extract buffer and parameter information.
        """
        param_block = header_string[header_string.find('(') + 1:header_string.rfind(')')]
        # Consolidate multi-line parameter block into a single line for easier parsing
        single_line_params = ' '.join(param_block.splitlines()).strip()
        # Split the consolidated string into individual parameter declarations
        params = [p.strip() for p in single_line_params.split(',') if p.strip()]
        
        args = []
        for param in params:
            if not param:
                continue
            # Determine argument type based on keywords
            if '__global' in param:
                # This is a buffer argument
                param_name = param.split('*')[-1].replace(',', '').strip()
                args.append((param_name, 'buffer'))
            elif 'image' in param and ('__read_only' in param or '__write_only' in param):
                # This is an image argument
                param_name = param.split()[-1].replace(',', '').strip()
                args.append((param_name, 'image'))
            else:
                # Assume it's a scalar argument
                param_name = param.split()[-1].replace(',', '').strip()
                args.append((param_name, 'scalar'))
        return args

    def get_kernel_args(self, program_name, kernel_name):
        """
        Generate argument list for a kernel based on its header definition and current buffers/parameters.
        """
        if program_name not in self.programs:
            raise ValueError(f"Program '{program_name}' not loaded.")
        if kernel_name not in self.kernel_headers[program_name]:
            raise ValueError(f"Kernel '{kernel_name}' not found in program '{program_name}'.")

        kernel_header = self.kernel_headers[program_name][kernel_name]
        arg_info = self._parse_kernel_header(kernel_header)

        args = []
        for name, arg_type in arg_info:
            if arg_type == 'buffer':
                if name not in self.buffers or self.buffers[name] is None:
                    raise ValueError(f"Required buffer '{name}' for kernel '{kernel_name}' is not allocated or is zero-sized.")
                args.append(self.buffers[name])
            elif arg_type == 'scalar':
                if name not in self.kernel_params:
                    raise ValueError(f"Required scalar parameter '{name}' for kernel '{kernel_name}' is not set.")
                args.append(self.kernel_params[name])
            elif arg_type == 'image':
                # Handle image arguments (not fully implemented in current scope, placeholder)
                if name not in self.buffers or self.buffers[name] is None:
                    raise ValueError(f"Required image '{name}' for kernel '{kernel_name}' is not allocated or is zero-sized.")
                args.append(self.buffers[name]) # Assuming images are also stored as buffers for now
            else:
                raise ValueError(f"Unknown argument type '{arg_type}' for argument '{name}'.")
        return args

    def set_kernel_param(self, name, value):
        """
        Set a scalar parameter for kernels.
        """
        self.kernel_params[name] = value

    def get_kernel_param(self, name):
        """
        Get a scalar parameter by name.
        """
        return self.kernel_params.get(name)

    def execute_kernel_with_args(self, kernel_obj, global_size, local_size, args):
        """
        Execute an OpenCL kernel with pre-prepared arguments.
        Args:
            kernel_obj (pyopencl.Kernel): The actual PyOpenCL kernel object.
            global_size (tuple or int): Global work size.
            local_size (tuple or int, optional): Local work size. Defaults to None.
            args (list): List of arguments for the kernel.
        """
        kernel_obj(self.queue, global_size, local_size, *args)
        self.queue.finish()

    def execute_kernel(self, program_name, kernel_name, global_size, local_size=None, args=None):
        """
        Execute an OpenCL kernel.
        Args:
            program_name (str): Name of the loaded program.
            kernel_name (str): Name of the kernel function.
            global_size (tuple or int): Global work size.
            local_size (tuple or int, optional): Local work size. Defaults to None.
            args (list, optional): List of arguments for the kernel. If None, arguments are generated automatically.
        """
        if program_name not in self.programs:
            raise ValueError(f"Program '{program_name}' not loaded.")
        
        prg = self.programs[program_name]
        # Retrieve kernel from cache or create if not exists
        if kernel_name not in self.kernels:
            self.kernels[kernel_name] = getattr(prg, kernel_name)
        kernel = self.kernels[kernel_name]

        if args is None:
            args = self.get_kernel_args(program_name, kernel_name)

        kernel(self.queue, global_size, local_size, *args)
        self.queue.finish()

    def clear_buffers(self):
        """
        Clear all allocated buffers and kernel parameters.
        """
        self.buffers.clear()
        self.kernel_params.clear()
        self.kernels.clear()
        self.programs.clear()
