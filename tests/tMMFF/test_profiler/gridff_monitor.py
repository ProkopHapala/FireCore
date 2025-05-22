#!/usr/bin/env python3
"""
GridFF Monitor

This module provides a way to monitor GridFF execution and capture
performance metrics without modifying the original code.
"""

import os
import sys
import time
import json
import datetime
import argparse
import traceback
import builtins
import numpy as np

# Store the original open function
original_open = builtins.open

# Store the original numpy save function
original_save = np.save

# Function to patch file operations
def patch_file_operations():
    # Get the absolute path to the project root
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
    
    # Define paths to important files
    kernel_path = os.path.join(project_root, 'cpp/common_resources/cl/GridFF.cl')
    element_types_path = os.path.join(project_root, 'cpp/common_resources/ElementTypes.dat')
    
    # Create output directory for redirected saves
    output_dir = os.path.join(os.path.dirname(__file__), 'profile_results')
    os.makedirs(output_dir, exist_ok=True)
    
    # Define the patched open function
    def patched_open(file, *args, **kwargs):
        # Check if this is an attempt to open the GridFF.cl file
        if isinstance(file, str) and file.endswith('GridFF.cl'):
            print(f"Redirecting GridFF.cl path from '{file}' to '{kernel_path}'")
            return original_open(kernel_path, *args, **kwargs)
        # Check if this is an attempt to open the ElementTypes.dat file
        elif isinstance(file, str) and file.endswith('ElementTypes.dat'):
            print(f"Redirecting ElementTypes.dat path from '{file}' to '{element_types_path}'")
            return original_open(element_types_path, *args, **kwargs)
        # Check if this is an attempt to write to a file that doesn't exist
        elif isinstance(file, str) and len(args) > 0 and 'w' in args[0] and not os.path.exists(os.path.dirname(file)):
            # Redirect to our output directory
            new_path = os.path.join(output_dir, os.path.basename(file))
            print(f"Redirecting write from '{file}' to '{new_path}'")
            return original_open(new_path, *args, **kwargs)
        # For all other files, use the original open function
        return original_open(file, *args, **kwargs)
    
    # Define the patched numpy save function
    def patched_save(file, arr, *args, **kwargs):
        if isinstance(file, str) and not os.path.exists(os.path.dirname(file)):
            # Redirect to our output directory
            new_path = os.path.join(output_dir, os.path.basename(file))
            print(f"Redirecting numpy save from '{file}' to '{new_path}'")
            return original_save(new_path, arr, *args, **kwargs)
        return original_save(file, arr, *args, **kwargs)
    
    # Replace the built-in open function
    builtins.open = patched_open
    
    # Replace the numpy save function
    np.save = patched_save
    
    print("File operations patched for GridFF")

# Function to set up OpenCL profiling environment variables
def setup_opencl_profiling():
    # Basic OpenCL profiling
    os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
    os.environ['OPENCL_PROFILE'] = '1'
    
    # NVIDIA-specific profiling
    os.environ['CUDA_PROFILE'] = '1'
    os.environ['CUDA_PROFILE_CSV'] = '1'
    os.environ['CUDA_PROFILE_LOG'] = 'cuda_profile_%p.log'
    
    # AMD-specific profiling
    os.environ['AMD_OCL_BUILD_OPTIONS_APPEND'] = '-cl-opt-disable'
    os.environ['AMD_OPENCL_PROFILER_ENABLE'] = '1'
    
    # Intel-specific profiling
    os.environ['INTEL_OPENCL_PROFILE'] = '1'
    
    print("OpenCL profiling environment variables set")

# Function to run GridFF with monitoring
def run_gridff_with_monitoring(job="PLQ", dataset="NaCl_1x1_L1", output_dir="profile_results"):
    # Apply patches
    patch_file_operations()
    setup_opencl_profiling()
    
    # Add FireCore project root to Python path
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate timestamp for output files
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Output files
    summary_file = os.path.join(output_dir, f"profile_{job}_{dataset}_{timestamp}_summary.txt")
    
    # Record start time
    start_time = time.time()
    
    # Run GridFF test
    try:
        print(f"\nRunning GridFF test with job={job} and dataset={dataset}...")
        
        # Find the correct path for the dataset file
        dataset_path = os.path.join(project_root, 'cpp/common_resources/xyz', f"{dataset}.xyz")
        
        # Verify the dataset file exists
        if not os.path.exists(dataset_path):
            print(f"Error: Dataset file not found at {dataset_path}")
            return 1
        
        print(f"Using dataset file: {dataset_path}")
        
        # Import the GridFF test module
        from pyBall.tests import ocl_GridFF_new as gff
        
        # Create a try-except block to catch any errors during the test
        try:
            # Run the test function
            gff.test_gridFF_ocl(job=job, fname=dataset_path)
        except Exception as e:
            print(f"Warning: GridFF test encountered an error: {e}")
            print("Continuing with monitoring results...")
        
        # Record end time
        end_time = time.time()
        total_time = end_time - start_time
        
        # Write summary file
        with open(summary_file, 'w') as f:
            f.write(f"GridFF Execution Summary\n")
            f.write(f"======================\n\n")
            f.write(f"Job: {job}\n")
            f.write(f"Dataset: {dataset}\n")
            f.write(f"Timestamp: {timestamp}\n\n")
            f.write(f"Total execution time: {total_time:.2f} seconds\n\n")
            
            # Add environment variables
            f.write(f"OpenCL Environment Variables:\n")
            for var in ['PYOPENCL_COMPILER_OUTPUT', 'OPENCL_PROFILE', 'CUDA_PROFILE', 
                       'CUDA_PROFILE_CSV', 'AMD_OCL_BUILD_OPTIONS_APPEND', 'INTEL_OPENCL_PROFILE']:
                f.write(f"  {var}={os.environ.get(var, 'Not set')}\n")
        
        print(f"\n{'='*50}")
        print(f"Monitoring complete for GridFF")
        print(f"Job: {job}, Dataset: {dataset}")
        print(f"Total execution time: {total_time:.2f} seconds\n")
        print(f"Summary saved to: {summary_file}")
        
        return 0
    
    except Exception as e:
        print(f"Error running GridFF test: {e}")
        traceback.print_exc()
        return 1

# Main function
def main():
    parser = argparse.ArgumentParser(description='GridFF Monitor')
    parser.add_argument('--job', default='PLQ', help='Job type for GridFF (PLQ, Morse, etc.)')
    parser.add_argument('--dataset', default='NaCl_1x1_L1', help='Dataset name for GridFF')
    parser.add_argument('--output', default='profile_results', help='Output directory for monitoring results')
    
    args = parser.parse_args()
    
    return run_gridff_with_monitoring(args.job, args.dataset, args.output)

if __name__ == '__main__':
    sys.exit(main())
