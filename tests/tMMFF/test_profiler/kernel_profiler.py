#!/usr/bin/env python3
"""
OpenCL Kernel Profiler

This module provides a way to profile OpenCL kernel execution times
without modifying any of the original code. It patches PyOpenCL
to track kernel execution times and redirects file paths as needed.
"""

import os
import sys
import time
import json
import datetime
import argparse
import traceback
import importlib
import builtins
import numpy as np
import pyopencl as cl

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
        elif isinstance(file, str) and 'w' in args and not os.path.exists(os.path.dirname(file)):
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

# Global variable to store kernel execution times
kernel_times = {}

# Function to patch PyOpenCL for kernel profiling
def patch_pyopencl():
    global kernel_times
    
    # Store the original CommandQueue constructor
    original_queue_init = cl.CommandQueue.__init__
    
    # Create patched CommandQueue constructor
    def patched_queue_init(self, context, device=None, properties=None):
        # Add profiling flag to properties
        if properties is None:
            properties = cl.command_queue_properties.PROFILING_ENABLE
        else:
            properties |= cl.command_queue_properties.PROFILING_ENABLE
        
        # Call original constructor
        original_queue_init(self, context, device, properties)
    
    # Apply the patch
    cl.CommandQueue.__init__ = patched_queue_init
    
    # Patch the enqueue_nd_range_kernel function
    original_enqueue = cl.enqueue_nd_range_kernel
    
    def patched_enqueue(queue, kernel, global_size, local_size=None, *args, **kwargs):
        # Execute the kernel and get the event
        event = original_enqueue(queue, kernel, global_size, local_size, *args, **kwargs)
        
        try:
            # Wait for the event to complete
            event.wait()
            
            # Calculate execution time
            start_time = event.profile.start
            end_time = event.profile.end
            execution_time = (end_time - start_time) * 1e-9  # Convert nanoseconds to seconds
            
            # Get kernel name
            kernel_name = kernel.function_name
            
            # Create a unique identifier for this kernel invocation
            kernel_id = f"{kernel_name}:{global_size}:{local_size}"
            
            # Update kernel statistics
            if kernel_id not in kernel_times:
                kernel_times[kernel_id] = {
                    'name': kernel_name,
                    'global_size': str(global_size),
                    'local_size': str(local_size),
                    'calls': 0,
                    'total_time': 0.0,
                    'min_time': float('inf'),
                    'max_time': 0.0
                }
            
            kernel_times[kernel_id]['calls'] += 1
            kernel_times[kernel_id]['total_time'] += execution_time
            kernel_times[kernel_id]['min_time'] = min(kernel_times[kernel_id]['min_time'], execution_time)
            kernel_times[kernel_id]['max_time'] = max(kernel_times[kernel_id]['max_time'], execution_time)
            
            # Print kernel execution information
            print(f"Kernel: {kernel_name}, Time: {execution_time*1000:.3f} ms, Global: {global_size}, Local: {local_size}")
        except Exception as e:
            print(f"Error capturing kernel profiling data: {e}")
        
        return event
    
    # Apply the patch
    cl.enqueue_nd_range_kernel = patched_enqueue
    
    print("PyOpenCL patched for kernel profiling")
    return kernel_times

# Function to save kernel profiling results
def save_kernel_results(kernel_times, output_dir, job, dataset):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate timestamp for output files
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Output files
    csv_file = os.path.join(output_dir, f"profile_{job}_{dataset}_{timestamp}_kernels.csv")
    json_file = os.path.join(output_dir, f"profile_{job}_{dataset}_{timestamp}_kernels.json")
    summary_file = os.path.join(output_dir, f"profile_{job}_{dataset}_{timestamp}_summary.txt")
    
    # Check if we have any kernel data
    if not kernel_times:
        print("Warning: No kernel execution data captured!")
        with open(summary_file, 'w') as f:
            f.write(f"GridFF OpenCL Kernel Profiling Summary\n")
            f.write(f"====================================\n\n")
            f.write(f"Job: {job}\n")
            f.write(f"Dataset: {dataset}\n")
            f.write(f"Timestamp: {timestamp}\n\n")
            f.write("No kernel execution data captured.\n")
            f.write("\nPossible reasons:\n")
            f.write("1. PyOpenCL profiling patch was not applied correctly\n")
            f.write("2. No OpenCL kernels were executed\n")
            f.write("3. The OpenCL implementation doesn't support profiling\n")
        return csv_file, json_file, summary_file
    
    # Write CSV file
    with open(csv_file, 'w') as f:
        f.write("Kernel,GlobalSize,LocalSize,Calls,TotalTime,AvgTime,MinTime,MaxTime\n")
        for kernel_id, stats in sorted(kernel_times.items(), key=lambda x: x[1]['total_time'], reverse=True):
            avg_time = stats['total_time'] / stats['calls'] if stats['calls'] > 0 else 0
            f.write(f"{stats['name']},{stats['global_size']},{stats['local_size']},{stats['calls']},{stats['total_time']},{avg_time},{stats['min_time']},{stats['max_time']}\n")
    
    # Write JSON file
    with open(json_file, 'w') as f:
        json.dump(kernel_times, f, indent=2)
    
    # Write summary file
    with open(summary_file, 'w') as f:
        f.write(f"GridFF OpenCL Kernel Profiling Summary\n")
        f.write(f"====================================\n\n")
        f.write(f"Job: {job}\n")
        f.write(f"Dataset: {dataset}\n")
        f.write(f"Timestamp: {timestamp}\n\n")
        
        f.write(f"{'Kernel':<30} {'Global Size':<15} {'Local Size':<15} {'Calls':<8} {'Total Time (s)':<15} {'Avg Time (s)':<15} {'Min Time (s)':<15} {'Max Time (s)':<15}\n")
        f.write("-" * 130 + "\n")
        
        # Sort kernels by total time
        for kernel_id, stats in sorted(kernel_times.items(), key=lambda x: x[1]['total_time'], reverse=True):
            avg_time = stats['total_time'] / stats['calls'] if stats['calls'] > 0 else 0
            f.write(f"{stats['name']:<30} {stats['global_size']:<15} {stats['local_size']:<15} {stats['calls']:<8} {stats['total_time']:<15.6f} {avg_time:<15.6f} {stats['min_time']:<15.6f} {stats['max_time']:<15.6f}\n")
    
    return csv_file, json_file, summary_file
    
    return csv_file, json_file, summary_file

# Function to run GridFF with profiling
def run_gridff_with_profiling(job="PLQ", dataset="NaCl_1x1_L1", output_dir="profile_results"):
    # Apply patches
    patch_file_operations()
    patch_pyopencl()
    
    # Add FireCore project root to Python path
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    
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
            print("Continuing with profiling results...")
        
        # Record end time
        end_time = time.time()
        total_time = end_time - start_time
        
        # Save kernel profiling results
        global kernel_times
        csv_file, json_file, summary_file = save_kernel_results(kernel_times, output_dir, job, dataset)
        
        print(f"\n{'='*50}")
        print(f"Profiling complete for GridFF")
        print(f"Job: {job}, Dataset: {dataset}")
        print(f"Total execution time: {total_time:.2f} seconds\n")
        print(f"Kernel profiling results saved to:")
        print(f"  CSV: {csv_file}")
        print(f"  JSON: {json_file}")
        print(f"  Summary: {summary_file}")
        
        return 0
    
    except Exception as e:
        print(f"Error running GridFF test: {e}")
        traceback.print_exc()
        return 1

# Main function
def main():
    parser = argparse.ArgumentParser(description='OpenCL Kernel Profiler')
    parser.add_argument('--job', default='PLQ', help='Job type for GridFF (PLQ, Morse, etc.)')
    parser.add_argument('--dataset', default='NaCl_1x1_L1', help='Dataset name for GridFF')
    parser.add_argument('--output', default='profile_results', help='Output directory for profiling results')
    
    args = parser.parse_args()
    
    return run_gridff_with_profiling(args.job, args.dataset, args.output)

if __name__ == '__main__':
    sys.exit(main())
