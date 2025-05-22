#!/usr/bin/env python3
"""
OpenCL Kernel Profiler

This module provides a way to profile OpenCL kernel execution times
by patching the PyOpenCL library.
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
    parser = argparse.ArgumentParser(description='OpenCL Kernel Profiler')
    parser.add_argument('script', help='Python script to profile')
    parser.add_argument('args', nargs='*', help='Arguments to pass to the script')
    parser.add_argument('--output', default='profile_results', help='Output directory for profiling results')
    
    args = parser.parse_args()
    
    # Enable OpenCL profiling
    enable_opencl_profiling()
    
    # Show OpenCL information
    show_opencl_info()
    
    # Patch PyOpenCL
    if not patch_pyopencl():
        print("Failed to patch PyOpenCL, kernel profiling will not be available")
        return 1
    
    # Add script's directory to Python path
    script_dir = os.path.dirname(os.path.abspath(args.script))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    
    # Add project root to Python path (assuming FireCore structure)
    project_root = os.path.abspath(os.path.join(script_dir, '..', '..', '..'))
    if os.path.exists(os.path.join(project_root, 'pyBall')):
        print(f"Adding project root to Python path: {project_root}")
        if project_root not in sys.path:
            sys.path.insert(0, project_root)
    
    # Save original sys.argv and set it to the target script's args
    original_argv = sys.argv
    sys.argv = [args.script] + args.args
    
    # Record start time
    start_time = time.time()
    
    try:
        # Execute the script
        print(f"\nRunning {args.script} with kernel profiling...")
        
        with open(args.script, 'rb') as f:
            code = compile(f.read(), args.script, 'exec')
            exec(code, {'__name__': '__main__'})
    
    except Exception as e:
        print(f"Error running script: {e}")
        traceback.print_exc()
    
    finally:
        # Calculate total time
        end_time = time.time()
        total_time = end_time - start_time
        
        print(f"\n{'='*50}")
        print(f"Profiling complete for {args.script}")
        print(f"Total execution time: {total_time:.2f} seconds")
        
        # Restore original sys.argv
        sys.argv = original_argv
        
        # Save kernel profiling results
        save_kernel_stats(args.output)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
