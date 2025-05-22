#!/usr/bin/env python3
"""
Comprehensive Grid Generation Profiler

This script provides detailed profiling for the grid generation process,
tracking both CPU and GPU operations across the execution flow.
"""

import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt
import pyopencl as cl
import csv

sys.path.append("../../")

# Import the necessary modules
from pyBall.tests import ocl_GridFF_new as gff
from pyBall.OCL import GridFF

# Dictionary to store timing information
timings = {}
kernel_timings = {}

# Function to start timing
def start_timing(name):
    timings[name] = {'start': time.time()}

# Function to stop timing
def stop_timing(name):
    if name in timings:
        timings[name]['end'] = time.time()
        timings[name]['elapsed'] = (timings[name]['end'] - timings[name]['start']) * 1000  # ms
        return timings[name]['elapsed']
    return 0

# Function to record kernel execution time
def record_kernel_time(name, elapsed_ms):
    if name not in kernel_timings:
        kernel_timings[name] = []
    kernel_timings[name].append(elapsed_ms)

# Create a subclass of GridFF_cl with profiling capabilities
class ProfilingGridFF_cl(GridFF.GridFF_cl):
    def __init__(self, *args, **kwargs):
        # Initialize the parent class
        super().__init__(*args, **kwargs)
        
        # Enable profiling for OpenCL queue
        if self.ctx is not None and hasattr(self.ctx, 'devices') and len(self.ctx.devices) > 0:
            # Create a new queue with profiling enabled
            self.queue = cl.CommandQueue(self.ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
            
        # Dictionary to store kernel events for profiling
        self.kernel_events = {}
        
        # Track memory usage
        self.memory_usage = {
            'buffers': {},
            'total_allocated': 0,
            'peak_allocated': 0,
            'allocation_history': []
        }
    
    # Track buffer allocation
    def track_buffer_allocation(self, name, size_bytes):
        self.memory_usage['buffers'][name] = size_bytes
        self.memory_usage['total_allocated'] += size_bytes
        if self.memory_usage['total_allocated'] > self.memory_usage['peak_allocated']:
            self.memory_usage['peak_allocated'] = self.memory_usage['total_allocated']
        
        # Record allocation history with timestamp
        self.memory_usage['allocation_history'].append({
            'time': time.time(),
            'operation': f'Allocate {name}',
            'size_bytes': size_bytes,
            'total_allocated': self.memory_usage['total_allocated']
        })
    
    # Track buffer deallocation
    def track_buffer_deallocation(self, name):
        if name in self.memory_usage['buffers']:
            size_bytes = self.memory_usage['buffers'][name]
            self.memory_usage['total_allocated'] -= size_bytes
            del self.memory_usage['buffers'][name]
            
            # Record deallocation history
            self.memory_usage['allocation_history'].append({
                'time': time.time(),
                'operation': f'Deallocate {name}',
                'size_bytes': -size_bytes,
                'total_allocated': self.memory_usage['total_allocated']
            })
    
    # Override buffer creation to track memory
    def try_make_buff(self, name, size_bytes, *args, **kwargs):
        # Call the original method
        result = super().try_make_buff(name, size_bytes, *args, **kwargs)
        
        # Track the allocation
        self.track_buffer_allocation(name, size_bytes)
        
        return result
    
    # Profile kernel execution
    def profile_kernel(self, kernel_name, kernel_func, *args, **kwargs):
        # Start CPU timing
        start_time = time.time()
        
        # Execute the kernel and get the event
        event = kernel_func(*args, **kwargs)
        
        # Wait for the kernel to complete
        if event and hasattr(event, 'wait'):
            event.wait()
            
            # Calculate GPU execution time in milliseconds
            if hasattr(event, 'profile'):
                start_time_ns = event.profile.start
                end_time_ns = event.profile.end
                execution_time_ms = (end_time_ns - start_time_ns) / 1e6
                
                # Record the kernel execution time
                record_kernel_time(kernel_name, execution_time_ms)
                
                # Store the event for later analysis
                if kernel_name not in self.kernel_events:
                    self.kernel_events[kernel_name] = []
                self.kernel_events[kernel_name].append(event)
        
        # Calculate CPU overhead (includes kernel execution + Python overhead)
        cpu_time_ms = (time.time() - start_time) * 1000
        record_kernel_time(f"{kernel_name}_cpu", cpu_time_ms)
        
        return event
    
    # Override methods to add profiling
    def makeCoulombEwald(self, atoms, bOld=False):
        method_name = 'makeCoulombEwald'
        start_timing(method_name)
        
        # Store original methods to restore later
        original_methods = {}
        
        # Intercept kernel calls by monkey patching the kernel execution
        if hasattr(self.prg, 'slabPotential'):
            original_methods['slabPotential'] = self.prg.slabPotential
            
            def profiled_slabPotential(*args, **kwargs):
                return self.profile_kernel('slabPotential', original_methods['slabPotential'], *args, **kwargs)
            
            self.prg.slabPotential = profiled_slabPotential
        
        try:
            result = super().makeCoulombEwald(atoms, bOld)
        finally:
            # Restore original methods
            for name, method in original_methods.items():
                setattr(self.prg, name, method)
        
        stop_timing(method_name)
        return result
    
    def make_MorseFF(self, atoms, REQs, nPBC=(4, 4, 0), dg=(0.1, 0.1, 0.1), ng=None, 
                    lvec=[[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]], 
                    g0=(0.0, 0.0, 0.0), GFFParams=(0.1, 1.5, 0.0, 0.0), bTime=True, bReturn=True):
        method_name = 'make_MorseFF'
        start_timing(method_name)
        
        # Store original methods to restore later
        original_methods = {}
        
        # Intercept kernel calls
        if hasattr(self.prg, 'evalMorse'):
            original_methods['evalMorse'] = self.prg.evalMorse
            
            def profiled_evalMorse(*args, **kwargs):
                return self.profile_kernel('evalMorse', original_methods['evalMorse'], *args, **kwargs)
            
            self.prg.evalMorse = profiled_evalMorse
        
        try:
            result = super().make_MorseFF(atoms, REQs, nPBC, dg, ng, lvec, g0, GFFParams, bTime, bReturn)
        finally:
            # Restore original methods
            for name, method in original_methods.items():
                setattr(self.prg, name, method)
        
        stop_timing(method_name)
        return result
    
    def make_Coulomb_points(self, atoms, ps, nPBC=(25, 25, 0), lvec=None, Ls=None, 
                           GFFParams=(0.1, 1.5, 0.0, 0.0), bTime=True, bReturn=True):
        method_name = 'make_Coulomb_points'
        start_timing(method_name)
        
        # Store original methods to restore later
        original_methods = {}
        
        # Intercept kernel calls
        if hasattr(self.prg, 'evalCoulomb'):
            original_methods['evalCoulomb'] = self.prg.evalCoulomb
            
            def profiled_evalCoulomb(*args, **kwargs):
                return self.profile_kernel('evalCoulomb', original_methods['evalCoulomb'], *args, **kwargs)
            
            self.prg.evalCoulomb = profiled_evalCoulomb
        
        try:
            result = super().make_Coulomb_points(atoms, ps, nPBC, lvec, Ls, GFFParams, bTime, bReturn)
        finally:
            # Restore original methods
            for name, method in original_methods.items():
                setattr(self.prg, name, method)
        
        stop_timing(method_name)
        return result
    
    def set_grid(self, gsh):
        method_name = 'set_grid'
        start_timing(method_name)
        result = super().set_grid(gsh)
        stop_timing(method_name)
        return result
    
    def sample3D(self, ps, buff):
        method_name = 'sample3D'
        start_timing(method_name)
        
        # Store original methods to restore later
        original_methods = {}
        
        # Intercept kernel calls
        if hasattr(self.prg, 'sample3D'):
            original_methods['sample3D'] = self.prg.sample3D
            
            def profiled_sample3D(*args, **kwargs):
                return self.profile_kernel('sample3D_kernel', original_methods['sample3D'], *args, **kwargs)
            
            self.prg.sample3D = profiled_sample3D
        
        try:
            result = super().sample3D(ps, buff)
        finally:
            # Restore original methods
            for name, method in original_methods.items():
                setattr(self.prg, name, method)
        
        stop_timing(method_name)
        return result
    
    def sample3D_grid(self, V_buff, samp_grid, bReturn=True):
        method_name = 'sample3D_grid'
        start_timing(method_name)
        
        # Store original methods to restore later
        original_methods = {}
        
        # Intercept kernel calls
        if hasattr(self.prg, 'sample3D_grid'):
            original_methods['sample3D_grid'] = self.prg.sample3D_grid
            
            def profiled_sample3D_grid(*args, **kwargs):
                return self.profile_kernel('sample3D_grid_kernel', original_methods['sample3D_grid'], *args, **kwargs)
            
            self.prg.sample3D_grid = profiled_sample3D_grid
        
        try:
            result = super().sample3D_grid(V_buff, samp_grid, bReturn)
        finally:
            # Restore original methods
            for name, method in original_methods.items():
                setattr(self.prg, name, method)
        
        stop_timing(method_name)
        return result
    
    def poisson(self, bReturn=True, sh=None, dV=None):
        method_name = 'poisson'
        start_timing(method_name)
        
        # Store original methods to restore later
        original_methods = {}
        
        # Intercept kernel calls
        kernels_to_profile = ['initFFT', 'forwardFFT', 'solvePoisson', 'inverseFFT']
        for kernel_name in kernels_to_profile:
            if hasattr(self.prg, kernel_name):
                original_methods[kernel_name] = getattr(self.prg, kernel_name)
                
                # Create a closure to capture the current kernel_name
                def make_profiled_kernel(k_name, original_method):
                    def profiled_kernel(*args, **kwargs):
                        return self.profile_kernel(k_name, original_method, *args, **kwargs)
                    return profiled_kernel
                
                setattr(self.prg, kernel_name, make_profiled_kernel(kernel_name, original_methods[kernel_name]))
        
        try:
            result = super().poisson(bReturn, sh, dV)
        finally:
            # Restore original methods
            for name, method in original_methods.items():
                setattr(self.prg, name, method)
        
        stop_timing(method_name)
        return result

# Function to profile GPU device information
def profile_gpu_info():
    print("\n===== GPU Information =====")
    
    # Get OpenCL platforms and devices
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

# Monkey patch the visualization code in ocl_GridFF_new.py to avoid index errors
def safe_test_gridFF_ocl(*args, **kwargs):
    # Store the original function
    original_func = gff.test_gridFF_ocl
    
    # Define a wrapper function to catch index errors
    def wrapper(*args, **kwargs):
        # Save the original imshow function
        original_imshow = plt.imshow
        
        # Define a safe imshow function that checks array bounds
        def safe_imshow(arr, **kwargs):
            try:
                return original_imshow(arr, **kwargs)
            except IndexError as e:
                print(f"Warning: IndexError in visualization: {e}")
                # Return a dummy imshow with a small valid array
                return original_imshow(np.zeros((2, 2)), **kwargs)
        
        # Replace plt.imshow with our safe version
        plt.imshow = safe_imshow
        
        try:
            # Add a flag to disable visualization if needed
            if 'no_vis' in kwargs and kwargs['no_vis']:
                # Store the original show function
                original_show = plt.show
                # Replace with a dummy function
                plt.show = lambda: None
                
                # Remove the no_vis kwarg before passing to the original function
                del kwargs['no_vis']
            
            # Call the original function
            result = original_func(*args, **kwargs)
            return result
        finally:
            # Restore the original functions
            plt.imshow = original_imshow
            if 'original_show' in locals():
                plt.show = original_show
    
    return wrapper

# Function to run profiling on a specific job
def run_profiling(job="PLQ", name="NaCl_1x1_L1", iterations=1):
    print(f"\n===== Profiling {job} job on {name} dataset ({iterations} iterations) =====")
    
    # Create a profiling instance of GridFF_cl
    profiling_clgff = ProfilingGridFF_cl()
    
    # Save the original GridFF_cl instance and test_gridFF_ocl function
    original_clgff = gff.clgff
    original_test_func = gff.test_gridFF_ocl
    
    # Replace with our profiling instance and patched test function
    gff.clgff = profiling_clgff
    gff.test_gridFF_ocl = safe_test_gridFF_ocl()
    
    try:
        # Run multiple iterations for more accurate profiling
        for i in range(iterations):
            if iterations > 1:
                print(f"\nIteration {i+1}/{iterations}...")
            
            # Profile the overall execution
            start_timing(f"total_{job}")
            
            # Profile data loading
            start_timing("data_loading")
            atoms = gff.make_atoms_arrays(
                fname=f"data/xyz/{name}.xyz",
                Element_Types_name="./data/ElementTypes.dat"
            )
            stop_timing("data_loading")
            
            # Run the test with visualization disabled for profiling
            start_timing(f"test_{job}")
            try:
                gff.test_gridFF_ocl(
                    fname=f"data/xyz/{name}.xyz",
                    Element_Types_name="./data/ElementTypes.dat",
                    job=job,
                    save_name="profile_run",
                    no_vis=True  # Disable visualization for profiling
                )
            except Exception as e:
                print(f"Warning: Error during test execution: {e}")
                print("Continuing with profiling...")
            stop_timing(f"test_{job}")
            
            stop_timing(f"total_{job}")
    finally:
        # Restore the original GridFF_cl instance and test function
        gff.clgff = original_clgff
        gff.test_gridFF_ocl = original_test_func
    
    # Include memory usage information if available
    memory_usage = None
    if hasattr(profiling_clgff, 'memory_usage'):
        memory_usage = profiling_clgff.memory_usage
    
    # Return the timing and memory results
    return {
        'timings': timings,
        'kernel_timings': kernel_timings,
        'memory_usage': memory_usage
    }

# Function to display profiling results
def display_profiling_results(results, job, name):
    timings = results['timings']
    kernel_timings = results['kernel_timings']
    memory_usage = results.get('memory_usage', None)
    
    print("\n===== Profiling Results =====")
    
    # Sort timings by elapsed time
    sorted_timings = sorted(
        [(k, v['elapsed']) for k, v in timings.items() if 'elapsed' in v],
        key=lambda x: x[1],
        reverse=True
    )
    
    # Print overall timings
    print("\nOverall Timings:")
    print(f"{'Operation':<40} {'Time (ms)':<15} {'% of Total':<15}")
    print("-" * 70)
    
    total_time = next((t for op, t in sorted_timings if op == f"total_{job}"), None)
    if total_time is None and sorted_timings:
        total_time = sorted_timings[0][1]  # Use the longest operation as total if total not found
    
    for operation, time_ms in sorted_timings:
        percent = (time_ms / total_time * 100) if total_time else 0
        print(f"{operation:<40} {time_ms:<15.2f} {percent:<15.2f}%")
    
    # Print kernel timings if available
    if kernel_timings:
        print("\nGPU Kernel Timings:")
        print(f"{'Kernel':<40} {'Calls':<8} {'Total (ms)':<15} {'Avg (ms)':<15} {'% of GPU':<15}")
        print("-" * 95)
        
        # Calculate total GPU time (excluding CPU overhead)
        gpu_kernels = [k for k in kernel_timings.keys() if not k.endswith('_cpu')]
        total_gpu_time = sum([sum(kernel_timings[k]) for k in gpu_kernels]) if gpu_kernels else 0
        
        sorted_kernels = sorted(
            [(k, len(v), sum(v), sum(v)/len(v) if v else 0) for k, v in kernel_timings.items()],
            key=lambda x: x[2],  # Sort by total time
            reverse=True
        )
        
        for kernel, calls, kernel_time, avg_time in sorted_kernels:
            # Skip CPU overhead entries for this display
            if kernel.endswith('_cpu'):
                continue
                
            percent = (kernel_time / total_gpu_time * 100) if total_gpu_time > 0 else 0
            print(f"{kernel:<40} {calls:<8} {kernel_time:<15.2f} {avg_time:<15.2f} {percent:<15.2f}%")
        
        # Print CPU vs GPU time comparison
        print("\nCPU vs GPU Time:")
        print(f"{'Operation':<40} {'GPU Time (ms)':<15} {'CPU Time (ms)':<15} {'Overhead (ms)':<15} {'Overhead %':<15}")
        print("-" * 100)
        
        # Group kernel timings by operation (removing _cpu suffix)
        gpu_ops = {}
        cpu_ops = {}
        
        for k, v in kernel_timings.items():
            if k.endswith('_cpu'):
                base_name = k[:-4]  # Remove _cpu suffix
                cpu_ops[base_name] = sum(v)
            else:
                gpu_ops[k] = sum(v)
        
        # Display CPU vs GPU time for each operation
        for op in sorted(set(gpu_ops.keys()).union(set(cpu_ops.keys()))):
            gpu_time = gpu_ops.get(op, 0)
            cpu_time = cpu_ops.get(op, 0)
            overhead = cpu_time - gpu_time
            overhead_percent = (overhead / gpu_time * 100) if gpu_time > 0 else 0
            
            print(f"{op:<40} {gpu_time:<15.2f} {cpu_time:<15.2f} {overhead:<15.2f} {overhead_percent:<15.2f}%")
    
    # Print memory usage if available
    if memory_usage:
        print("\nMemory Usage:")
        print(f"Peak GPU memory allocated: {memory_usage['peak_allocated'] / (1024*1024):.2f} MB")
        
        # Print largest buffers
        if memory_usage['buffers']:
            print("\nLargest Buffers:")
            print(f"{'Buffer':<30} {'Size (MB)':<15}")
            print("-" * 45)
            
            sorted_buffers = sorted(
                [(name, size) for name, size in memory_usage['buffers'].items()],
                key=lambda x: x[1],
                reverse=True
            )
            
            for name, size in sorted_buffers[:10]:  # Show top 10 buffers
                print(f"{name:<30} {size / (1024*1024):<15.2f}")
    
    # Save results to CSV
    with open(f"profile_{job}_{name}.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Category', 'Operation', 'Time (ms)', '% of Total'])
        
        # Write overall timings
        for operation, time_ms in sorted_timings:
            percent = (time_ms / total_time * 100) if total_time else 0
            writer.writerow(['Overall', operation, f"{time_ms:.2f}", f"{percent:.2f}"])
        
        # Write kernel timings
        if kernel_timings:
            for kernel, calls, total_time, avg_time in sorted_kernels:
                writer.writerow(['Kernel', kernel, f"{total_time:.2f}", f"{calls}"])
        
        # Write memory usage
        if memory_usage and memory_usage['buffers']:
            for name, size in sorted_buffers:
                writer.writerow(['Memory', name, f"{size / (1024*1024):.2f} MB", ""])
    
    print(f"\nResults saved to profile_{job}_{name}.csv")
    
    # Create visualizations
    plt.figure(figsize=(15, 15))
    
    # Plot overall timings
    plt.subplot(3, 1, 1)
    
    # Filter out very small timings and the total time for better visualization
    significant_timings = [(op, t) for op, t in sorted_timings 
                          if t > total_time * 0.01 and not op.startswith('total_')]
    
    operations = [op for op, _ in significant_timings]
    times = [t for _, t in significant_timings]
    
    if operations:  # Only create plot if we have data
        y_pos = range(len(operations))
        plt.barh(y_pos, times)
        plt.yticks(y_pos, operations)
        plt.xlabel('Time (ms)')
        plt.title(f'Operation Execution Times - {job} job on {name} dataset')
        
        # Add time labels
        for i, v in enumerate(times):
            plt.text(v + 0.1, i, f"{v:.2f} ms", va='center')
    
    # Plot kernel timings if available
    if kernel_timings:
        plt.subplot(3, 1, 2)
        
        # Get top 10 kernels by time (excluding CPU overhead)
        gpu_kernels = [(k, c, t, a) for k, c, t, a in sorted_kernels if not k.endswith('_cpu')][:10]
        
        kernel_names = [k for k, _, _, _ in gpu_kernels]
        kernel_times = [t for _, _, t, _ in gpu_kernels]
        
        if kernel_names:  # Only create plot if we have data
            y_pos = range(len(kernel_names))
            plt.barh(y_pos, kernel_times)
            plt.yticks(y_pos, kernel_names)
            plt.xlabel('Time (ms)')
            plt.title('GPU Kernel Execution Times')
            
            # Add time labels
            for i, v in enumerate(kernel_times):
                plt.text(v + 0.1, i, f"{v:.2f} ms", va='center')
        
        # Plot CPU vs GPU comparison
        plt.subplot(3, 1, 3)
        
        # Prepare data for CPU vs GPU comparison
        ops = []
        gpu_times = []
        cpu_times = []
        
        # Find operations with both CPU and GPU measurements
        for op in sorted(set(gpu_ops.keys())):
            if op in gpu_ops and op + '_cpu' in kernel_timings:
                ops.append(op)
                gpu_times.append(gpu_ops[op])
                cpu_times.append(cpu_ops.get(op, 0))
        
        if ops:  # Only create plot if we have data
            x = range(len(ops))
            width = 0.35
            
            plt.bar(x, gpu_times, width, label='GPU Time')
            plt.bar([i + width for i in x], cpu_times, width, label='CPU Time (incl. overhead)')
            
            plt.xlabel('Operation')
            plt.ylabel('Time (ms)')
            plt.title('CPU vs GPU Execution Time Comparison')
            plt.xticks([i + width/2 for i in x], ops, rotation=45, ha='right')
            plt.legend()
    
    plt.tight_layout()
    plt.savefig(f"profile_{job}_{name}.png")
    print(f"Visualization saved to profile_{job}_{name}.png")
    
    # Create a separate memory usage plot if available
    if memory_usage and memory_usage['allocation_history']:
        plt.figure(figsize=(12, 6))
        
        # Extract time and memory usage from history
        times = []
        memory_values = []
        start_time = memory_usage['allocation_history'][0]['time']
        
        for entry in memory_usage['allocation_history']:
            times.append(entry['time'] - start_time)  # Time since start
            memory_values.append(entry['total_allocated'] / (1024*1024))  # MB
        
        plt.plot(times, memory_values)
        plt.xlabel('Time (s)')
        plt.ylabel('Memory Usage (MB)')
        plt.title('GPU Memory Usage Over Time')
        plt.grid(True)
        
        # Add peak memory annotation
        peak_memory = memory_usage['peak_allocated'] / (1024*1024)
        plt.axhline(y=peak_memory, color='r', linestyle='--')
        plt.text(times[-1] * 0.5, peak_memory * 1.05, f'Peak: {peak_memory:.2f} MB', color='r')
        
        plt.tight_layout()
        plt.savefig(f"memory_profile_{job}_{name}.png")
        print(f"Memory usage visualization saved to memory_profile_{job}_{name}.png")
    
    # Return the figure for further customization if needed
    return plt.gcf()

# Main function
def main():
    # Parse command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Profile Grid Generation Process')
    parser.add_argument('--job', type=str, default='PLQ', help='Job type (PLQ, Morse, Ewald)')
    parser.add_argument('--name', type=str, default='NaCl_1x1_L1', help='Dataset name')
    parser.add_argument('--iterations', type=int, default=1, help='Number of iterations for averaging')
    args = parser.parse_args()
    
    # Profile GPU information
    profile_gpu_info()
    
    # Run profiling
    results = run_profiling(args.job, args.name, args.iterations)
    
    # Display results
    display_profiling_results(results, args.job, args.name)
    
    # Show plots
    plt.show()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error during profiling: {e}")
        import traceback
        traceback.print_exc()
