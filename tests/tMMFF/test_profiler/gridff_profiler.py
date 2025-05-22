#!/usr/bin/env python3
"""
GridFF Profiler

This module provides a specialized profiler for GridFF OpenCL code.
It patches the GridFF_cl class to properly load the OpenCL kernel file
and tracks kernel execution times.
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

# Global variables to store profiling data
kernel_stats = defaultdict(lambda: {'calls': 0, 'total_time': 0, 'avg_time': 0, 'min_time': float('inf'), 'max_time': 0})

# Function to patch PyOpenCL for kernel profiling
def patch_pyopencl():
    try:
        import pyopencl as cl
        
        # Save original enqueue_nd_range_kernel function
        original_enqueue = cl.enqueue_nd_range_kernel
        
        # Create patched function
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
            
            # If profiling is not enabled, create a new queue with profiling enabled
            if not has_profiling:
                try:
                    ctx = queue.context
                    new_queue = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
                    queue = new_queue
                    has_profiling = True
                    print(f"Created new queue with profiling enabled for kernel {kernel_name}")
                except Exception as e:
                    print(f"Warning: Could not create queue with profiling enabled: {e}")
            
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
        
        # Apply patch
        cl.enqueue_nd_range_kernel = profiled_enqueue
        
        print("PyOpenCL patched for kernel profiling")
        return True
    except ImportError:
        print("PyOpenCL not available, skipping kernel profiling")
        return False
    except Exception as e:
        print(f"Error patching PyOpenCL: {e}")
        traceback.print_exc()
        return False

# Function to patch GridFF_cl to use the correct kernel path
def patch_gridff(kernel_path=None):
    try:
        # Find the absolute path to the GridFF.cl kernel file
        if kernel_path is None:
            # Try to find the kernel file
            project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
            kernel_path = os.path.join(project_root, 'cpp/common_resources/cl/GridFF.cl')
            
            if not os.path.exists(kernel_path):
                print(f"Warning: Kernel file not found at {kernel_path}")
                # Try to find the kernel file in other locations
                possible_paths = [
                    os.path.join(project_root, 'cpp/common_resources/cl/GridFF.cl'),
                    os.path.join(project_root, 'pyBall/OCL/cl/GridFF.cl'),
                    '/home/indranil/git/FireCore/cpp/common_resources/cl/GridFF.cl'
                ]
                
                for path in possible_paths:
                    if os.path.exists(path):
                        kernel_path = path
                        print(f"Found kernel file at {kernel_path}")
                        break
        
        # Import GridFF_cl
        from pyBall.OCL.GridFF import GridFF_cl
        
        # Save the original __init__ method
        original_init = GridFF_cl.__init__
        
        # Create a patched __init__ method
        def patched_init(self, *args, **kwargs):
            # Call the original __init__
            original_init(self, *args, **kwargs)
            
            # Save the original init_cl method
            original_init_cl = self.init_cl
            
            # Create a patched init_cl method
            def patched_init_cl(device=None, platform_id=None):
                # Call the original init_cl to set up the context and queue
                original_init_cl(device, platform_id)
                
                # Now override the OpenCL program loading
                try:
                    print(f"Loading OpenCL kernel from {kernel_path}")
                    with open(kernel_path, 'r') as f:
                        self.prg = self.cl.Program(self.ctx, f.read()).build()
                    return True
                except Exception as e:
                    print(f"Error compiling OpenCL program: {e}")
                    traceback.print_exc()
                    return False
            
            # Replace the init_cl method
            self.init_cl = patched_init_cl
        
        # Apply the patch
        GridFF_cl.__init__ = patched_init
        
        print(f"GridFF_cl patched to use kernel file: {kernel_path}")
        return True
    except ImportError:
        print("Could not import pyBall.OCL.GridFF, skipping GridFF patching")
        return False
    except Exception as e:
        print(f"Error patching GridFF_cl: {e}")
        traceback.print_exc()
        return False

# Function to save kernel profiling results
def save_kernel_stats(output_dir="profile_results"):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate timestamp for output files
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Output files
    kernels_csv = os.path.join(output_dir, f"profile_{timestamp}_kernels.csv")
    kernels_json = os.path.join(output_dir, f"profile_{timestamp}_kernels.json")
    summary_txt = os.path.join(output_dir, f"profile_{timestamp}_kernel_summary.txt")
    
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
        f.write(f"OpenCL Kernel Profiling Summary\n")
        f.write(f"=============================\n\n")
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

# Function to show OpenCL device information
def show_opencl_info():
    try:
        import pyopencl as cl
        print("\nOpenCL Information:")
        print("-" * 50)
        
        # Get platform and device information
        platforms = cl.get_platforms()
        for i, platform in enumerate(platforms):
            print(f"Platform {i}: {platform.name} ({platform.version})")
            
            devices = platform.get_devices()
            for j, device in enumerate(devices):
                print(f"  Device {j}: {device.name} (Type: {cl.device_type.to_string(device.type)})")
                print(f"    Compute Units: {device.max_compute_units}")
                print(f"    Global Memory: {device.global_mem_size / (1024**2):.2f} MB")
                print(f"    Local Memory: {device.local_mem_size / 1024:.2f} KB")
                print(f"    Max Work Group Size: {device.max_work_group_size}")
                if hasattr(device, 'max_clock_frequency'):
                    print(f"    Max Clock Frequency: {device.max_clock_frequency} MHz")
    except ImportError:
        print("PyOpenCL not available, skipping OpenCL information")
    except Exception as e:
        print(f"Error getting OpenCL information: {e}")

# Function to enable OpenCL profiling via environment variables
def enable_opencl_profiling():
    # Set environment variables for OpenCL profiling
    os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
    os.environ['OPENCL_PROFILE'] = '1'
    
    # For NVIDIA
    os.environ['CUDA_PROFILE'] = '1'
    os.environ['CUDA_PROFILE_CSV'] = '1'
    
    # For AMD
    os.environ['AMD_OCL_BUILD_OPTIONS_APPEND'] = '-cl-opt-disable'
    
    # For Intel
    os.environ['INTEL_OPENCL_PROFILE'] = '1'
    
    print("OpenCL profiling enabled via environment variables")

# Main function
def main():
    parser = argparse.ArgumentParser(description='GridFF OpenCL Profiler')
    parser.add_argument('--job', default='PLQ', help='Job type (PLQ, Morse, etc.)')
    parser.add_argument('--dataset', default='NaCl_1x1_L1', help='Dataset name')
    parser.add_argument('--output', default='profile_results', help='Output directory for profiling results')
    parser.add_argument('--kernel-path', help='Path to GridFF.cl kernel file')
    
    args = parser.parse_args()
    
    # Enable OpenCL profiling
    enable_opencl_profiling()
    
    # Show OpenCL information
    show_opencl_info()
    
    # Patch PyOpenCL
    if not patch_pyopencl():
        print("Failed to patch PyOpenCL, kernel profiling will not be available")
        return 1
    
    # Patch GridFF_cl
    if not patch_gridff(args.kernel_path):
        print("Failed to patch GridFF_cl, using original implementation")
    
    # Add project root to Python path
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
        print(f"Added project root to Python path: {project_root}")
    
    # Set environment variables for the test
    os.environ['GRIDFF_JOB'] = args.job
    os.environ['GRIDFF_DATASET'] = args.dataset
    
    # Record start time
    start_time = time.time()
    
    try:
        # Import the test module
        print(f"\nRunning GridFF test with job={args.job} and dataset={args.dataset}...")
        from pyBall.tests import ocl_GridFF_new as gff
        
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
        save_kernel_stats(args.output)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
