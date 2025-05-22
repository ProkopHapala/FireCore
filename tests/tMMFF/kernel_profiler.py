#!/usr/bin/env python

import sys
import os
import numpy as np
import time
import pyopencl as cl
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from functools import wraps
from collections import defaultdict

sys.path.append("../../")
import pyBall.tests.ocl_GridFF_new as gff
from pyBall.OCL.GridFF import GridFF_cl

# Global variables to store profiling results
kernel_timings = defaultdict(list)
function_timings = defaultdict(list)

# Function to profile kernel execution time
def profile_event(name, event):
    """Profile a PyOpenCL event"""
    event.wait()
    start = event.get_profiling_info(cl.profiling_info.START)
    end = event.get_profiling_info(cl.profiling_info.END)
    duration_ms = (end - start) * 1e-6  # nanoseconds to milliseconds
    
    # Store timing
    kernel_timings[name].append(duration_ms)
    return duration_ms

# Function to monkey-patch GridFF_cl methods with profiling
def patch_gridff():
    """Patch critical GridFF_cl methods with profiling"""
    original_methods = {}
    
    # Store original methods and create wrapped versions
    # 1. project_atoms_on_grid_quintic_pbc - key kernel for charge distribution
    original_methods['_project_atoms_on_grid_quintic_pbc'] = GridFF_cl._project_atoms_on_grid_quintic_pbc
    
    @wraps(GridFF_cl._project_atoms_on_grid_quintic_pbc)
    def project_atoms_profiled(self, sz_glob, sz_loc, na, ns):
        # Start timing
        start_time = time.time()
        
        # Ensure we use a profiling-enabled queue
        original_queue = self.queue
        profile_queue = cl.CommandQueue(self.ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
        self.queue = profile_queue
        
        # Execute the 'set' kernel
        print("\n=== Profiling project_atoms_on_grid kernel ===\n")
        nxyz2 = np.int32(ns[0]*ns[1]*ns[2] * 2)
        set_event = self.prg.set(self.queue, sz_glob, sz_loc, nxyz2, self.Qgrid_buff, np.float32(0.0))
        set_time = profile_event('set', set_event)
        print(f"Kernel: {'set':<20} Time: {set_time:>10.3f} ms")
        
        # Execute the projection kernel
        project_event = self.prg.project_atoms_on_grid_quintic_pbc(
            self.queue, sz_glob, sz_loc,
            na, self.atoms_buff, self.Qgrid_buff,
            ns, self.gcl.g0, self.gcl.dg
        )
        project_time = profile_event('project_atoms_on_grid_quintic_pbc', project_event)
        print(f"Kernel: {'project_atoms':<20} Time: {project_time:>10.3f} ms")
        
        # Restore original queue
        self.queue = original_queue
        
        # Record total function time
        func_time = (time.time() - start_time) * 1000  # seconds to ms
        function_timings['project_atoms'].append(func_time)
        print(f"Total function time: {func_time:.3f} ms (includes Python overhead)")
    
    # 2. poisson - key kernel for solving the Poisson equation
    original_methods['poisson'] = GridFF_cl.poisson
    
    @wraps(GridFF_cl.poisson)
    def poisson_profiled(self, bReturn=True, sh=None, dV=None):
        start_time = time.time()
        print("\n=== Profiling poisson solver ===\n")
        
        # Create a profiling queue
        original_queue = self.queue
        profile_queue = cl.CommandQueue(self.ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
        self.queue = profile_queue
        
        # Execute the method with profiling
        try:
            result = original_methods['poisson'](self, bReturn=bReturn, sh=sh, dV=dV)
        except Exception as e:
            print(f"Error in poisson: {e}")
            self.queue = original_queue
            raise e
        
        # Restore original queue
        self.queue = original_queue
        
        # Record total function time
        func_time = (time.time() - start_time) * 1000  # seconds to ms
        function_timings['poisson'].append(func_time)
        print(f"Total poisson time: {func_time:.3f} ms (includes FFT and Python overhead)")
        
        return result
    
    # 3. laplace_real_loop_inert - solver for the Poisson equation in real space
    original_methods['laplace_real_loop_inert'] = GridFF_cl.laplace_real_loop_inert
    
    @wraps(GridFF_cl.laplace_real_loop_inert)
    def laplace_profiled(self, niter=16, cSOR=0.0, cV=0.6, bReturn=False, sh=None):
        start_time = time.time()
        print("\n=== Profiling laplace_real_loop_inert ===\n")
        
        # Create a profiling queue
        original_queue = self.queue
        profile_queue = cl.CommandQueue(self.ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
        self.queue = profile_queue
        
        # Execute with profiling
        try:
            result = original_methods['laplace_real_loop_inert'](self, niter=niter, cSOR=cSOR, cV=cV, bReturn=bReturn, sh=sh)
        except Exception as e:
            print(f"Error in laplace_real_loop_inert: {e}")
            self.queue = original_queue
            raise e
        
        # Restore original queue
        self.queue = original_queue
        
        # Record total function time
        func_time = (time.time() - start_time) * 1000  # seconds to ms
        function_timings['laplace_real_loop'].append(func_time)
        print(f"Total laplace solver time: {func_time:.3f} ms")
        
        return result
    
    # 4. make_MorseFF - very important kernel for force field calculation
    original_methods['make_MorseFF'] = GridFF_cl.make_MorseFF
    
    @wraps(GridFF_cl.make_MorseFF)
    def make_morse_profiled(self, *args, **kwargs):
        start_time = time.time()
        print("\n=== Profiling make_MorseFF ===\n")
        
        # Create a profiling queue
        original_queue = self.queue
        profile_queue = cl.CommandQueue(self.ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
        self.queue = profile_queue
        
        # Execute with profiling
        try:
            result = original_methods['make_MorseFF'](self, *args, **kwargs)
        except Exception as e:
            print(f"Error in make_MorseFF: {e}")
            self.queue = original_queue
            raise e
        
        # Restore original queue
        self.queue = original_queue
        
        # Record total function time
        func_time = (time.time() - start_time) * 1000  # seconds to ms
        function_timings['make_MorseFF'].append(func_time)
        print(f"Total make_MorseFF time: {func_time:.3f} ms")
        
        return result
    
    # Patch the methods
    GridFF_cl._project_atoms_on_grid_quintic_pbc = project_atoms_profiled
    GridFF_cl.poisson = poisson_profiled
    GridFF_cl.laplace_real_loop_inert = laplace_profiled
    GridFF_cl.make_MorseFF = make_morse_profiled
    
    return original_methods

# Function to restore original methods
def restore_gridff(original_methods):
    GridFF_cl._project_atoms_on_grid_quintic_pbc = original_methods['_project_atoms_on_grid_quintic_pbc']
    GridFF_cl.poisson = original_methods['poisson']
    GridFF_cl.laplace_real_loop_inert = original_methods['laplace_real_loop_inert']
    GridFF_cl.make_MorseFF = original_methods['make_MorseFF']
    print("Original methods restored")

# Function to generate profiling visualizations
def generate_visualizations():
    if not kernel_timings and not function_timings:
        print("No profiling data collected")
        return
    
    # Create a figure for kernel timings
    if kernel_timings:
        plt.figure(figsize=(10, 6))
        kernel_names = list(kernel_timings.keys())
        avg_times = [np.mean(times) for times in kernel_timings.values()]
        call_counts = [len(times) for times in kernel_timings.values()]
        
        # Sort by average execution time
        sorted_indices = np.argsort(avg_times)[::-1]  # Descending order
        sorted_names = [kernel_names[i] for i in sorted_indices]
        sorted_times = [avg_times[i] for i in sorted_indices]
        sorted_counts = [call_counts[i] for i in sorted_indices]
        
        plt.barh(range(len(sorted_names)), sorted_times, color='skyblue')
        plt.yticks(range(len(sorted_names)), [f"{name} ({count} calls)" for name, count in zip(sorted_names, sorted_counts)])
        plt.xlabel('Average Execution Time (ms)')
        plt.title('OpenCL Kernel Performance')
        plt.tight_layout()
        plt.savefig('kernel_timings.png', dpi=150)
        print("Kernel timing visualization saved to kernel_timings.png")
    
    # Create a figure for function timings
    if function_timings:
        plt.figure(figsize=(10, 6))
        func_names = list(function_timings.keys())
        avg_times = [np.mean(times) for times in function_timings.values()]
        call_counts = [len(times) for times in function_timings.values()]
        
        # Sort by average execution time
        sorted_indices = np.argsort(avg_times)[::-1]  # Descending order
        sorted_names = [func_names[i] for i in sorted_indices]
        sorted_times = [avg_times[i] for i in sorted_indices]
        sorted_counts = [call_counts[i] for i in sorted_indices]
        
        plt.barh(range(len(sorted_names)), sorted_times, color='lightgreen')
        plt.yticks(range(len(sorted_names)), [f"{name} ({count} calls)" for name, count in zip(sorted_names, sorted_counts)])
        plt.xlabel('Average Execution Time (ms)')
        plt.title('Function Performance (Including Python Overhead)')
        plt.tight_layout()
        plt.savefig('function_timings.png', dpi=150)
        print("Function timing visualization saved to function_timings.png")

# Function to print detailed profiling results
def print_profiling_results():
    if not kernel_timings and not function_timings:
        print("No profiling data collected")
        return
    
    print("\n===== DETAILED PROFILING RESULTS =====\n")
    
    # Kernel timings
    if kernel_timings:
        print("KERNEL TIMINGS (OpenCL execution only):\n")
        print(f"{'Kernel':<30} {'Calls':>6} {'Avg (ms)':>10} {'Min (ms)':>10} {'Max (ms)':>10} {'Total (ms)':>12}")
        print("-" * 80)
        
        # Sort kernels by average execution time
        sorted_kernels = sorted(kernel_timings.items(), key=lambda x: np.mean(x[1]), reverse=True)
        
        for name, times in sorted_kernels:
            calls = len(times)
            avg = np.mean(times)
            min_time = np.min(times)
            max_time = np.max(times)
            total = np.sum(times)
            print(f"{name:<30} {calls:>6} {avg:>10.3f} {min_time:>10.3f} {max_time:>10.3f} {total:>12.3f}")
    
    # Function timings
    if function_timings:
        print("\nFUNCTION TIMINGS (Including Python overhead):\n")
        print(f"{'Function':<30} {'Calls':>6} {'Avg (ms)':>10} {'Min (ms)':>10} {'Max (ms)':>10} {'Total (ms)':>12}")
        print("-" * 80)
        
        # Sort functions by average execution time
        sorted_funcs = sorted(function_timings.items(), key=lambda x: np.mean(x[1]), reverse=True)
        
        for name, times in sorted_funcs:
            calls = len(times)
            avg = np.mean(times)
            min_time = np.min(times)
            max_time = np.max(times)
            total = np.sum(times)
            print(f"{name:<30} {calls:>6} {avg:>10.3f} {min_time:>10.3f} {max_time:>10.3f} {total:>12.3f}")

# Main profiling function
def profile_gridff(fname="data/xyz/NaCl_1x1_L1.xyz", job="PLQ"):
    # Apply profiling patches
    original_methods = patch_gridff()
    
    try:
        # Record total runtime
        start_time = time.time()
        
        # Run the test
        print(f"\nRunning GridFF with: {fname}, job={job}\n")
        try:
            gff.test_gridFF_ocl(
                fname=fname,
                Element_Types_name="./data/ElementTypes.dat",
                save_name="profiled",
                job=job
            )
        except Exception as e:
            print(f"Error during test execution: {e}")
        
        # Print timing
        elapsed = time.time() - start_time
        print(f"\n===== TOTAL EXECUTION TIME: {elapsed:.2f} seconds =====\n")
        
        # Print and visualize results
        print_profiling_results()
        generate_visualizations()
        
    finally:
        # Restore original methods
        restore_gridff(original_methods)

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
