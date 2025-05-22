#!/usr/bin/env python3
"""
GridFF Path Profiler

This module provides a way to profile GridFF OpenCL code by patching
the file open operation to correctly locate the GridFF.cl kernel file.
It doesn't modify any of the original code.
"""

import os
import sys
import time
import datetime
import argparse
import subprocess
import builtins
import traceback
import importlib
import pyopencl as cl

# Store the original open function
original_open = builtins.open

# Function to patch the open function
def patch_open():
    # Get the absolute path to the project root
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
    
    # Define paths to important files
    kernel_path = os.path.join(project_root, 'cpp/common_resources/cl/GridFF.cl')
    element_types_path = os.path.join(project_root, 'cpp/common_resources/ElementTypes.dat')
    
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
        # For all other files, use the original open function
        return original_open(file, *args, **kwargs)
    
    # Replace the built-in open function
    builtins.open = patched_open
    print("File open function patched for GridFF.cl")

# Function to patch PyOpenCL for kernel profiling
def patch_pyopencl():
    # Track kernel execution times
    kernel_times = {}
    
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
    
    # Patch enqueue_nd_range_kernel to track execution times
    original_enqueue = cl.enqueue_nd_range_kernel
    
    def patched_enqueue(queue, kernel, global_size, local_size=None, *args, **kwargs):
        # Add wait_for to kwargs if not present
        if 'wait_for' not in kwargs:
            kwargs['wait_for'] = None
        
        # Execute the kernel and get the event
        event = original_enqueue(queue, kernel, global_size, local_size, *args, **kwargs)
        
        # Wait for the event to complete
        event.wait()
        
        try:
            # Calculate execution time
            start_time = event.profile.start
            end_time = event.profile.end
            execution_time = (end_time - start_time) * 1e-9  # Convert nanoseconds to seconds
            
            # Update kernel statistics
            kernel_name = kernel.function_name
            if kernel_name not in kernel_times:
                kernel_times[kernel_name] = {'calls': 0, 'total_time': 0.0, 'min_time': float('inf'), 'max_time': 0.0}
            
            kernel_times[kernel_name]['calls'] += 1
            kernel_times[kernel_name]['total_time'] += execution_time
            kernel_times[kernel_name]['min_time'] = min(kernel_times[kernel_name]['min_time'], execution_time)
            kernel_times[kernel_name]['max_time'] = max(kernel_times[kernel_name]['max_time'], execution_time)
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
    
    # Write CSV file
    with open(csv_file, 'w') as f:
        f.write("Kernel,Calls,TotalTime,AvgTime,MinTime,MaxTime\n")
        for kernel, stats in sorted(kernel_times.items(), key=lambda x: x[1]['total_time'], reverse=True):
            avg_time = stats['total_time'] / stats['calls'] if stats['calls'] > 0 else 0
            f.write(f"{kernel},{stats['calls']},{stats['total_time']},{avg_time},{stats['min_time']},{stats['max_time']}\n")
    
    # Write JSON file
    import json
    with open(json_file, 'w') as f:
        json.dump(kernel_times, f, indent=2)
    
    # Write summary file
    with open(summary_file, 'w') as f:
        f.write(f"GridFF OpenCL Kernel Profiling Summary\n")
        f.write(f"====================================\n\n")
        f.write(f"Job: {job}\n")
        f.write(f"Dataset: {dataset}\n")
        f.write(f"Timestamp: {timestamp}\n\n")
        
        if kernel_times:
            f.write(f"{'Kernel':<40} {'Calls':<8} {'Total Time (s)':<15} {'Avg Time (s)':<15} {'Min Time (s)':<15} {'Max Time (s)':<15}\n")
            f.write("-" * 110 + "\n")
            
            # Sort kernels by total time
            for kernel, stats in sorted(kernel_times.items(), key=lambda x: x[1]['total_time'], reverse=True):
                avg_time = stats['total_time'] / stats['calls'] if stats['calls'] > 0 else 0
                f.write(f"{kernel:<40} {stats['calls']:<8} {stats['total_time']:<15.6f} {avg_time:<15.6f} {stats['min_time']:<15.6f} {stats['max_time']:<15.6f}\n")
        else:
            f.write("No kernel execution data captured.\n")
    
    return csv_file, json_file, summary_file

# Function to run GridFF with profiling
def run_gridff_with_profiling(job="PLQ", dataset="NaCl_1x1_L1", output_dir="profile_results"):
    # Apply patches
    patch_open()
    kernel_times = patch_pyopencl()
    
    # Add FireCore project root to Python path
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    
    # Record start time
    start_time = time.time()
    
    # Run GridFF test
    try:
        print(f"\nRunning GridFF test with job={job} and dataset={dataset}...")
        
        # Import the GridFF test module
        from pyBall.tests import ocl_GridFF_new as gff
        
        # Find the correct path for the dataset file
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
        dataset_path = os.path.join(project_root, 'cpp/common_resources/xyz', f"{dataset}.xyz")
        
        # Verify the dataset file exists
        if not os.path.exists(dataset_path):
            print(f"Error: Dataset file not found at {dataset_path}")
            return 1
        
        print(f"Using dataset file: {dataset_path}")
        
        # Run the test function
        gff.test_gridFF_ocl(job=job, fname=dataset_path)
        
        # Record end time
        end_time = time.time()
        total_time = end_time - start_time
        
        # Save kernel profiling results
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
    parser = argparse.ArgumentParser(description='GridFF Path Profiler')
    parser.add_argument('--job', default='PLQ', help='Job type for GridFF (PLQ, Morse, etc.)')
    parser.add_argument('--dataset', default='NaCl_1x1_L1', help='Dataset name for GridFF')
    parser.add_argument('--output', default='profile_results', help='Output directory for profiling results')
    
    args = parser.parse_args()
    
    return run_gridff_with_profiling(args.job, args.dataset, args.output)

if __name__ == '__main__':
    sys.exit(main())
