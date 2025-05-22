#!/usr/bin/env python

import sys
import os
import numpy as np
import time
import pyopencl as cl
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall.tests import ocl_GridFF_new as gff

# Dictionary to store kernel timings
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

# Patch the test_gridFF_ocl function to avoid visualization errors
def patch_test_function():
    # Store the original function
    original_test_function = gff.test_gridFF_ocl
    
    # Create a wrapper that skips the problematic visualization
    def wrapped_test_function(*args, **kwargs):
        # Run the original function but capture any IndexError during plotting
        try:
            result = original_test_function(*args, **kwargs)
            return result
        except IndexError as e:
            print(f"\nVisualization error caught: {e}")
            print("Skipping visualization, continuing with profiling...")
            # We still want to return some data for analysis
            return None
    
    # Replace the function
    gff.test_gridFF_ocl = wrapped_test_function
    return original_test_function

# Setup OpenCL profiling
def enable_profiling():
    # Patch kernels to enable profiling
    ctx = gff.clgff.ctx
    props = cl.command_queue_properties.PROFILING_ENABLE
    gff.clgff.queue = cl.CommandQueue(ctx, properties=props)
    print("\n===== OpenCL Profiling Enabled =====\n")
    
    # Get available kernels
    kernel_names = [name for name in dir(gff.clgff.prg) 
                   if not name.startswith('__') and 
                   callable(getattr(gff.clgff.prg, name))]
    
    # Store original kernels
    original_kernels = {}
    
    # Patch each kernel with profiling
    for name in kernel_names:
        # Keep reference to original kernel
        original_kernel = getattr(gff.clgff.prg, name)
        original_kernels[name] = original_kernel
        
        # Create profiling wrapper
        def make_wrapper(name, original):
            def wrapper(*args, **kwargs):
                # Execute and profile
                evt = original(*args, **kwargs)
                profile_kernel(name, evt)
                return evt
            return wrapper
        
        # Apply the wrapper
        setattr(gff.clgff.prg, name, make_wrapper(name, original_kernel))
        print(f"Profiling enabled for kernel: {name}")
    
    return original_kernels

def restore_profiling(original_kernels):
    # Restore original kernels
    for name, original in original_kernels.items():
        setattr(gff.clgff.prg, name, original)
    print(f"Restored {len(original_kernels)} original kernels")

# Print profiling summary
def print_summary():
    if not kernel_timings:
        print("No kernel executions recorded")
        return
    
    print("\n===== OpenCL KERNEL EXECUTION SUMMARY =====\n")
    print(f"{'Kernel':<40} {'Calls':>6} {'Total (ms)':>12} {'Average (ms)':>14} {'Min (ms)':>10} {'Max (ms)':>10}")
    print("-" * 100)
    
    # Sort kernels by total execution time (most expensive first)
    sorted_kernels = sorted(kernel_timings.items(), 
                            key=lambda x: sum(x[1]), 
                            reverse=True)
    
    overall_total = 0
    for name, times in sorted_kernels:
        calls = len(times)
        total = sum(times)
        avg = total / calls
        min_time = min(times)
        max_time = max(times)
        
        print(f"{name:<40} {calls:>6} {total:>12.3f} {avg:>14.3f} {min_time:>10.3f} {max_time:>10.3f}")
        overall_total += total
    
    print("-" * 100)
    print(f"{'TOTAL':<40} {sum(len(times) for times in kernel_timings.values()):>6} {overall_total:>12.3f}\n")
    
    # Create a bar chart showing kernel execution times
    plt.figure(figsize=(12, 8))
    
    # Prepare data for plotting
    kernel_names = [k for k, _ in sorted_kernels]
    total_times = [sum(times) for _, times in sorted_kernels]
    avg_times = [sum(times)/len(times) for _, times in sorted_kernels]
    calls = [len(times) for _, times in sorted_kernels]
    
    # Plot total times
    plt.subplot(1, 2, 1)
    plt.barh(range(len(kernel_names)), total_times, color='skyblue')
    plt.yticks(range(len(kernel_names)), kernel_names)
    plt.xlabel('Total Execution Time (ms)')
    plt.title('OpenCL Kernel Total Execution Times')
    plt.tight_layout()
    
    # Plot average times
    plt.subplot(1, 2, 2)
    plt.barh(range(len(kernel_names)), avg_times, color='lightgreen')
    plt.yticks(range(len(kernel_names)), [f"{name} ({calls[i]} calls)" for i, name in enumerate(kernel_names)])
    plt.xlabel('Average Execution Time (ms)')
    plt.title('OpenCL Kernel Average Execution Times')
    
    plt.tight_layout()
    plt.savefig('kernel_profile.png', dpi=150)
    print(f"Profiling visualization saved to kernel_profile.png")

# Main function
def profile_gridff(fname="data/xyz/NaCl_1x1_L1.xyz", job="PLQ"):
    # Patch the visualization function
    original_test = patch_test_function()
    
    try:
        # Enable OpenCL profiling
        original_kernels = enable_profiling()
        
        # Run with time measurement
        start_time = time.time()
        print(f"\nRunning GridFF with: {fname}, job={job}\n")  
        
        try:
            # Run the test
            gff.test_gridFF_ocl(
                fname=fname,
                Element_Types_name="./data/ElementTypes.dat",
                save_name="profiled",
                job=job
            )
        except Exception as e:
            print(f"\nError during execution: {e}")
        
        # Total execution time
        duration = time.time() - start_time
        print(f"\n===== TOTAL EXECUTION TIME: {duration:.2f} seconds =====\n")
        
        # Print and visualize summary
        print_summary()
        
    finally:
        # Restore original functions
        restore_profiling(original_kernels)
        gff.test_gridFF_ocl = original_test

# Main entry point
if __name__ == "__main__":
    # Parse command line arguments
    if len(sys.argv) > 1:
        name = sys.argv[1]
    else:
        name = "NaCl_1x1_L1"
        
    if len(sys.argv) > 2:
        job = sys.argv[2]
    else:
        job = "PLQ"
    
    profile_gridff(f"data/xyz/{name}.xyz", job)
