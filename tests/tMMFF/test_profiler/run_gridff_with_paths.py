#!/usr/bin/env python3
"""
Wrapper script for running GridFF with correct paths

This script sets the correct paths for the GridFF OpenCL kernels
and then runs the specified test script.
"""

import os
import sys
import time

# Get the project root directory
try:
    script_path = __file__
except NameError:
    script_path = sys.argv[0]

project_root = os.path.abspath(os.path.join(os.path.dirname(script_path), '../../..'))

# Add the project root to Python path if not already there
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Set environment variables for OpenCL kernel paths
os.environ['FIRECORE_ROOT'] = project_root

# Print environment information
print(f"FireCore Project Root: {project_root}")
print(f"Python Path: {sys.path}")

# Import the test module
from pyBall.tests import ocl_GridFF_new as gff

# Patch the GridFF_cl class to use the correct path for OpenCL kernels
original_init = gff.GridFF_cl.__init__

def patched_init(self, *args, **kwargs):
    # Get the current directory
    current_dir = os.getcwd()
    print(f"GridFF_cl() called from path= {current_dir}")
    
    # Call the original init
    original_init(self, *args, **kwargs)
    
    # Now patch the initialization method that loads the OpenCL kernel
    original_init_cl = self.init_cl
    
    def patched_init_cl(device=None, platform_id=None):
        # Set up the OpenCL context and queue
        if hasattr(self, 'ctx') and self.ctx is None:
            self.ctx, self.queue = gff.clu.get_nvidia_device()
        
        # Use absolute path for the kernel file
        kernel_path = '/home/indranil/git/FireCore/cpp/common_resources/cl/GridFF.cl'
        
        # Verify the kernel file exists
        if not os.path.exists(kernel_path):
            print(f"Error: Kernel file not found at {kernel_path}")
            return False
        
        print(f"Loading OpenCL kernel from {kernel_path}")
        
        # Compile the OpenCL program
        try:
            with open(kernel_path, 'r') as f:
                self.prg = self.cl.Program(self.ctx, f.read()).build()
            return True
        except Exception as e:
            print(f"Error compiling OpenCL program: {e}")
            return False
    
    # Replace the init_cl method
    self.init_cl = patched_init_cl

# Apply the patch
gff.GridFF_cl.__init__ = patched_init

# Run the test function
if __name__ == "__main__":
    # Set default parameters
    job = os.environ.get('GRIDFF_JOB', 'PLQ')
    dataset = os.environ.get('GRIDFF_DATASET', 'NaCl_1x1_L1')
    
    # Print test information
    print(f"Running GridFF test with job={job} and dataset={dataset}")
    
    # Start timing
    start_time = time.time()
    
    # Run the test
    gff.test_gridFF_ocl(job=job, fname=f"./data/xyz/{dataset}.xyz")
    
    # End timing
    end_time = time.time()
    print(f"Test completed in {end_time - start_time:.2f} seconds")
