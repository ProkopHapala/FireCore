#!/usr/bin/env python3
"""
Direct GridFF Profiler

This script directly runs the GridFF test with profiling,
using a monkeypatched version of the GridFF_cl class to load the kernel file
from an absolute path and track kernel execution times.
"""

import os
import sys
import time
import datetime
import argparse
import traceback
import csv
import json
from collections import defaultdict

# Get the project root directory
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))

# Add the project root to Python path
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Global variables to store profiling data
kernel_stats = defaultdict(lambda: {'calls': 0, 'total_time': 0, 'avg_time': 0, 'min_time': float('inf'), 'max_time': 0})

# Function to save kernel profiling results
def save_kernel_stats(output_dir="profile_results", job="PLQ", dataset="NaCl_1x1_L1"):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate timestamp for output files
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Output files
    kernels_csv = os.path.join(output_dir, f"profile_{job}_{dataset}_{timestamp}_kernels.csv")
    kernels_json = os.path.join(output_dir, f"profile_{job}_{dataset}_{timestamp}_kernels.json")
    summary_txt = os.path.join(output_dir, f"profile_{job}_{dataset}_{timestamp}_summary.txt")
    
    # Save kernel stats to CSV
    with open(kernels_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Kernel', 'Calls', 'Total Time (ms)', 'Avg Time (ms)', 'Min Time (ms)', 'Max Time (ms)'])
        
        for kernel_name, stats in sorted(kernel_stats.items(), key=lambda x: x[1]['total_time'], reverse=True):
            writer.writerow([
                kernel_name,
                stats['calls'],
                f"{stats['total_time']:.3f}",
                f"{stats['avg_time']:.3f}",
                f"{stats['min_time']:.3f}",
                f"{stats['max_time']:.3f}"
            ])
    
    # Save kernel stats to JSON
    with open(kernels_json, 'w') as f:
        json.dump(kernel_stats, f, indent=2)
    
    # Save summary to text file
    with open(summary_txt, 'w') as f:
        f.write(f"GridFF OpenCL Kernel Profiling Summary\n")
        f.write(f"===================================\n\n")
        f.write(f"Job: {job}\n")
        f.write(f"Dataset: {dataset}\n")
        f.write(f"Timestamp: {timestamp}\n\n")
        
        f.write(f"{'Kernel':<40} {'Calls':<8} {'Total (ms)':<12} {'Avg (ms)':<12} {'Min (ms)':<12} {'Max (ms)':<12}\n")
        f.write("-" * 100 + "\n")
        
        for kernel_name, stats in sorted(kernel_stats.items(), key=lambda x: x[1]['total_time'], reverse=True):
            f.write(f"{kernel_name:<40} {stats['calls']:<8} {stats['total_time']:<12.3f} {stats['avg_time']:<12.3f} "
                   f"{stats['min_time']:<12.3f} {stats['max_time']:<12.3f}\n")
    
    print(f"\nKernel profiling results saved to:")
    print(f"  CSV: {kernels_csv}")
    print(f"  JSON: {kernels_json}")
    print(f"  Summary: {summary_txt}")
    
    return kernels_csv, kernels_json, summary_txt

# Main function
def main():
    parser = argparse.ArgumentParser(description='Direct GridFF OpenCL Profiler')
    parser.add_argument('--job', default='PLQ', help='Job type (PLQ, Morse, etc.)')
    parser.add_argument('--dataset', default='NaCl_1x1_L1', help='Dataset name')
    parser.add_argument('--output', default='profile_results', help='Output directory for profiling results')
    
    args = parser.parse_args()
    
    # Set environment variables for OpenCL profiling
    os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
    os.environ['OPENCL_PROFILE'] = '1'
    os.environ['CUDA_PROFILE'] = '1'
    os.environ['CUDA_PROFILE_CSV'] = '1'
    
    # Import PyOpenCL and patch it for profiling
    import pyopencl as cl
    
    # Save original enqueue_nd_range_kernel function
    original_enqueue = cl.enqueue_nd_range_kernel
    
    # Create patched function for kernel profiling
    def profiled_enqueue(queue, kernel, global_work_size, local_work_size=None, *args, **kwargs):
        # Get kernel name
        kernel_name = kernel.function_name if hasattr(kernel, 'function_name') else str(kernel)
        
        # Check if profiling is enabled for this queue
        has_profiling = False
        try:
            if hasattr(queue, 'properties'):
                has_profiling = bool(queue.properties & cl.command_queue_properties.PROFILING_ENABLE)
        except Exception:
            pass
        
        # Execute kernel and get event
        start_time = time.time()
        event = original_enqueue(queue, kernel, global_work_size, local_work_size, *args, **kwargs)
        
        # Wait for kernel to complete
        event.wait()
        
        # Calculate CPU execution time
        cpu_time_ms = (time.time() - start_time) * 1000
        
        # Try to get GPU time if profiling is enabled
        gpu_time_ms = 0
        try:
            if has_profiling:
                gpu_start = event.get_profiling_info(cl.profiling_info.START)
                gpu_end = event.get_profiling_info(cl.profiling_info.END)
                gpu_time_ms = (gpu_end - gpu_start) / 1e6  # Convert ns to ms
        except Exception:
            # Silently fail - profiling info might not be available
            pass
        
        # Use GPU time if available, otherwise use CPU time
        exec_time_ms = gpu_time_ms if gpu_time_ms > 0 else cpu_time_ms
        
        # Record kernel execution time
        kernel_stats[kernel_name]['calls'] += 1
        kernel_stats[kernel_name]['total_time'] += exec_time_ms
        kernel_stats[kernel_name]['avg_time'] = kernel_stats[kernel_name]['total_time'] / kernel_stats[kernel_name]['calls']
        kernel_stats[kernel_name]['min_time'] = min(kernel_stats[kernel_name]['min_time'], exec_time_ms)
        kernel_stats[kernel_name]['max_time'] = max(kernel_stats[kernel_name]['max_time'], exec_time_ms)
        
        # Print kernel execution time
        time_source = "GPU" if gpu_time_ms > 0 else "CPU"
        print(f"Kernel {kernel_name}: {exec_time_ms:.3f} ms ({time_source} time)")
        
        return event
    
    # Apply patch to PyOpenCL
    cl.enqueue_nd_range_kernel = profiled_enqueue
    print("PyOpenCL patched for kernel profiling")
    
    # Now import and patch GridFF_cl
    from pyBall.OCL.GridFF import GridFF_cl
    
    # Save original __init__ method
    original_init = GridFF_cl.__init__
    
    # Create patched __init__ method to override the kernel file path
    def patched_init(self, *args, **kwargs):
        # Call original init
        original_init(self, *args, **kwargs)
        
        # Define a method to patch the OpenCL program loading
        def load_opencl_program(self):
            # Use absolute path for the kernel file
            kernel_path = os.path.join(project_root, 'cpp/common_resources/cl/GridFF.cl')
            
            # Verify the kernel file exists
            if not os.path.exists(kernel_path):
                print(f"Error: Kernel file not found at {kernel_path}")
                return False
            
            print(f"Loading OpenCL kernel from {kernel_path}")
            
            # Make sure we have a context and queue
            if not hasattr(self, 'ctx') or self.ctx is None:
                from pyBall.OCL import clUtils
                self.ctx, self.queue = clUtils.get_nvidia_device()
            
            # Compile the OpenCL program
            try:
                with open(kernel_path, 'r') as f:
                    self.prg = cl.Program(self.ctx, f.read()).build()
                return True
            except Exception as e:
                print(f"Error compiling OpenCL program: {e}")
                traceback.print_exc()
                return False
        
        # Add the method to the instance
        self.load_opencl_program = lambda: load_opencl_program(self)
    
    # Apply patch to GridFF_cl
    GridFF_cl.__init__ = patched_init
    print(f"GridFF_cl patched to use kernel file from {os.path.join(project_root, 'cpp/common_resources/cl/GridFF.cl')}")
    
    # Record start time
    start_time = time.time()
    
    try:
        # Import the test module
        print(f"\nRunning GridFF test with job={args.job} and dataset={args.dataset}...")
        from pyBall.tests import ocl_GridFF_new as gff
        
        # Create an instance of GridFF_cl and load the OpenCL program
        gridff = gff.GridFF_cl()
        print("Loading OpenCL program...")
        gridff.load_opencl_program()
        
        # Run the test function
        gff.test_gridFF_ocl(job=args.job, fname=f"./data/xyz/{args.dataset}.xyz")
    
    except Exception as e:
        print(f"Error running GridFF test: {e}")
        traceback.print_exc()
    
    finally:
        # Calculate total time
        end_time = time.time()
        total_time = end_time - start_time
        
        print(f"\n{'='*50}")
        print(f"Profiling complete for GridFF")
        print(f"Job: {args.job}, Dataset: {args.dataset}")
        print(f"Total execution time: {total_time:.2f} seconds")
        
        # Save kernel profiling results
        save_kernel_stats(args.output, args.job, args.dataset)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
