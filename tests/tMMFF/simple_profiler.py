#!/usr/bin/env python

import sys
import os
import numpy as np
import time
import pyopencl as cl

sys.path.append("../../")
from pyBall.tests import ocl_GridFF_new as gff

# Setup profiling command queue
def enable_opencl_profiling():
    """Replace the command queue with a profiling-enabled one"""
    ctx = gff.clgff.ctx
    gff.clgff.queue = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
    print("\n===== OpenCL Profiling Enabled =====\n")

# Dictionary to store kernel timing information
kernel_timings = {}
def profile_kernel(name, evt, details=None):
    """Record kernel execution timing"""
    evt.wait()
    start = evt.get_profiling_info(cl.profiling_info.START)
    end = evt.get_profiling_info(cl.profiling_info.END)
    duration_ms = (end - start) * 1e-6  # Convert nanoseconds to milliseconds
    
    if name not in kernel_timings:
        kernel_timings[name] = []
    kernel_timings[name].append(duration_ms)
    
    # Print immediate feedback
    print(f"Kernel: {name:<30} Time: {duration_ms:>10.3f} ms{' - ' + details if details else ''}")
    return duration_ms

# Wrap key kernel functions for profiling
def patch_kernels():
    # Get list of available kernels
    available_kernels = dir(gff.clgff.prg)
    originals = {}
    
    # Common kernel names we want to profile if they exist
    kernel_names = [
        "make_MorseFF", 
        "set",
        "project_atoms_on_grid_quintic_pbc", 
        "poissonW", 
        "poissonW_old",
        "laplace_real_pbc", 
        "slabPotential_zyx"
    ]
    
    # Only patch kernels that exist
    for kernel_name in kernel_names:
        if kernel_name in available_kernels:
            # Store original kernel
            original_kernel = getattr(gff.clgff.prg, kernel_name)
            originals[kernel_name] = original_kernel
            
            # Create wrapper function to profile this kernel
            def make_wrapper(kernel_name, original_kernel):
                def wrapper(*args, **kwargs):
                    # Special handling for iteration-based kernels
                    details = ""
                    if kernel_name == "laplace_real_pbc" and len(args) > 7:
                        details = f"Iteration {args[7]}"
                    
                    # Call the original kernel and profile it
                    evt = original_kernel(*args, **kwargs)
                    profile_kernel(kernel_name, evt, details)
                    return evt
                return wrapper
            
            # Set the wrapped function
            setattr(gff.clgff.prg, kernel_name, make_wrapper(kernel_name, original_kernel))
            print(f"Profiling enabled for kernel: {kernel_name}")
    
    return originals

def restore_kernels(originals):
    # Restore all patched kernels dynamically
    for kernel_name, original_function in originals.items():
        if hasattr(gff.clgff.prg, kernel_name):
            setattr(gff.clgff.prg, kernel_name, original_function)
            print(f"Restored original kernel: {kernel_name}")

# Wrap memory transfer operations
def track_transfer_time():
    # Keep track of data movement
    transfer_stats = {
        "host_to_device": {"count": 0, "bytes": 0, "time_ms": 0},
        "device_to_host": {"count": 0, "bytes": 0, "time_ms": 0}
    }
    
    # Store original enqueue_copy
    original_enqueue_copy = cl.enqueue_copy
    
    def wrapped_enqueue_copy(queue, dest, src, *args, **kwargs):
        # Determine transfer direction
        is_buffer_src = isinstance(src, cl.Buffer)
        direction = "device_to_host" if is_buffer_src else "host_to_device"
        
        # Try to determine transfer size
        size = 0
        if hasattr(src, 'nbytes'):
            size = src.nbytes
        elif hasattr(dest, 'nbytes'):
            size = dest.nbytes
        
        # Profile the transfer
        start_time = time.time()
        evt = original_enqueue_copy(queue, dest, src, *args, **kwargs)
        
        if queue.properties & cl.command_queue_properties.PROFILING_ENABLE:
            evt.wait()  # Wait for the transfer to complete
            start = evt.get_profiling_info(cl.profiling_info.START)
            end = evt.get_profiling_info(cl.profiling_info.END)
            duration_ms = (end - start) * 1e-6  # Convert ns to ms
            
            # Update statistics
            transfer_stats[direction]["count"] += 1
            transfer_stats[direction]["bytes"] += size
            transfer_stats[direction]["time_ms"] += duration_ms
            
            # Print immediate feedback for large transfers
            if size > 1024*1024:  # Only print for transfers > 1MB
                size_mb = size / (1024*1024)
                print(f"Transfer: {direction:<15} Size: {size_mb:>8.2f} MB Time: {duration_ms:>8.2f} ms Bandwidth: {size_mb/(duration_ms/1000):>8.2f} MB/s")
        
        return evt
    
    # Apply our wrapper
    cl.enqueue_copy = wrapped_enqueue_copy
    
    return transfer_stats, original_enqueue_copy

def restore_transfer(original_enqueue_copy):
    cl.enqueue_copy = original_enqueue_copy

# Print summary statistics
def print_summary():
    if not kernel_timings:
        print("No kernel executions recorded.")
        return
        
    print("\n===== KERNEL EXECUTION SUMMARY =====\n")
    print(f"{'Kernel':<30} {'Calls':>6} {'Total (ms)':>12} {'Average (ms)':>14} {'Min (ms)':>10} {'Max (ms)':>10}")
    print("-" * 90)
    
    # Calculate and print statistics for each kernel
    total_time = 0
    kernels_by_total_time = sorted(kernel_timings.items(), key=lambda x: sum(x[1]), reverse=True)
    
    for name, times in kernels_by_total_time:
        num_calls = len(times)
        total = sum(times)
        avg = total / num_calls
        min_time = min(times)
        max_time = max(times)
        
        print(f"{name:<30} {num_calls:>6} {total:>12.3f} {avg:>14.3f} {min_time:>10.3f} {max_time:>10.3f}")
        total_time += total
    
    print("-" * 90)
    print(f"{'Total':<30} {sum(len(times) for times in kernel_timings.values()):>6} {total_time:>12.3f}")
    
    # Memory transfer summary
    if hasattr(cl, '_transfer_stats'):
        h2d = cl._transfer_stats["host_to_device"]
        d2h = cl._transfer_stats["device_to_host"]
        
        print("\n===== MEMORY TRANSFER SUMMARY =====\n")
        print(f"Direction       Count     Size (MB)   Time (ms)   Bandwidth (MB/s)")
        print("-" * 70)
        
        if h2d["count"] > 0:
            h2d_mb = h2d["bytes"] / (1024*1024)
            h2d_bw = h2d_mb / (h2d["time_ms"]/1000) if h2d["time_ms"] > 0 else 0
            print(f"Host to Device {h2d['count']:>6} {h2d_mb:>12.2f} {h2d['time_ms']:>12.2f} {h2d_bw:>12.2f}")
            
        if d2h["count"] > 0:
            d2h_mb = d2h["bytes"] / (1024*1024)
            d2h_bw = d2h_mb / (d2h["time_ms"]/1000) if d2h["time_ms"] > 0 else 0
            print(f"Device to Host {d2h['count']:>6} {d2h_mb:>12.2f} {d2h['time_ms']:>12.2f} {d2h_bw:>12.2f}")

# Main function to run profiling
def profile_gridff(fname="data/xyz/NaCl_1x1_L1.xyz", job="PLQ"):
    # Enable OpenCL profiling
    enable_opencl_profiling()
    
    # Patch kernels
    originals = patch_kernels()
    
    # Setup transfer tracking
    cl._transfer_stats, original_enqueue_copy = track_transfer_time()
    
    try:
        # Start timing
        overall_start = time.time()
        
        # Run the actual code
        print(f"Running: {fname} with job={job}")
        gff.test_gridFF_ocl(
            fname=fname,
            Element_Types_name="./data/ElementTypes.dat",
            save_name="profiled",
            job=job
        )
        
        # Print overall time
        overall_time = time.time() - overall_start
        print(f"\n===== TOTAL EXECUTION TIME: {overall_time:.4f}s =====\n")
        
        # Print detailed summary
        print_summary()
        
    finally:
        # Restore original functions
        restore_kernels(originals)
        restore_transfer(original_enqueue_copy)

# Main entry point
if __name__ == "__main__":
    # Get command line arguments
    if len(sys.argv) > 1:
        name = sys.argv[1]
    else:
        name = "NaCl_1x1_L1"
        
    if len(sys.argv) > 2:
        job = sys.argv[2]
    else:
        job = "PLQ"
    
    profile_gridff(f"data/xyz/{name}.xyz", job)
