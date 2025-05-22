#!/usr/bin/env python

import sys
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import pyopencl as cl
import pandas as pd

sys.path.append("../../")
from pyBall.tests import ocl_GridFF_new as gff

# Create profiling result storage
kernel_times = {}
kernel_calls = {}

class ProfilingMonitor:
    def __init__(self):
        self.kernel_times = {}
        self.kernel_calls = {}
        self.memory_transfers = {}
        self.start_time = None
        self.overall_times = {}
        
    def start_timer(self, name):
        self.start_time = time.time()
        
    def stop_timer(self, name):
        if self.start_time is not None:
            elapsed = time.time() - self.start_time
            if name not in self.overall_times:
                self.overall_times[name] = 0
            self.overall_times[name] += elapsed
            self.start_time = None
    
    def record_kernel(self, name, event):
        # Wait for event to complete
        event.wait()
        
        # Record profiling info
        start = event.get_profiling_info(cl.profiling_info.START)
        end = event.get_profiling_info(cl.profiling_info.END)
        duration_ns = end - start
        duration_ms = duration_ns / 1000000  # Convert to milliseconds
        
        if name not in self.kernel_times:
            self.kernel_times[name] = []
            self.kernel_calls[name] = 0
        
        self.kernel_times[name].append(duration_ms)
        self.kernel_calls[name] += 1
        
    def record_memory_transfer(self, direction, size_bytes, event):
        # Wait for transfer to complete
        event.wait()
        
        # Record profiling info
        start = event.get_profiling_info(cl.profiling_info.START)
        end = event.get_profiling_info(cl.profiling_info.END)
        duration_ns = end - start
        duration_ms = duration_ns / 1000000  # Convert to milliseconds
        bandwidth = size_bytes / duration_ms * 1000 / (1024*1024)  # MB/s
        
        if direction not in self.memory_transfers:
            self.memory_transfers[direction] = {'sizes': [], 'times': [], 'bandwidths': []}
            
        self.memory_transfers[direction]['sizes'].append(size_bytes)
        self.memory_transfers[direction]['times'].append(duration_ms)
        self.memory_transfers[direction]['bandwidths'].append(bandwidth)
    
    def print_stats(self):
        print("\n===== OPENCL KERNEL PROFILING =====\n")
        
        # Print kernel statistics
        kernel_df = pd.DataFrame({
            'Kernel': list(self.kernel_times.keys()),
            'Calls': [self.kernel_calls[k] for k in self.kernel_times.keys()],
            'Total Time (ms)': [sum(self.kernel_times[k]) for k in self.kernel_times.keys()],
            'Average Time (ms)': [np.mean(self.kernel_times[k]) for k in self.kernel_times.keys()],
            'Min Time (ms)': [np.min(self.kernel_times[k]) for k in self.kernel_times.keys()],
            'Max Time (ms)': [np.max(self.kernel_times[k]) for k in self.kernel_times.keys()]
        })
        kernel_df = kernel_df.sort_values('Total Time (ms)', ascending=False)
        print(kernel_df.to_string(index=False))
        
        # Print memory transfer statistics
        if self.memory_transfers:
            print("\n===== MEMORY TRANSFER PROFILING =====\n")
            for direction, data in self.memory_transfers.items():
                transfer_df = pd.DataFrame({
                    'Size (MB)': [s/(1024*1024) for s in data['sizes']],
                    'Time (ms)': data['times'],
                    'Bandwidth (MB/s)': data['bandwidths']
                })
                print(f"\n{direction} - {len(data['sizes'])} transfers:")
                print(f"  Total: {sum(data['sizes'])/(1024*1024):.2f} MB in {sum(data['times']):.2f} ms")
                print(f"  Average Bandwidth: {np.mean(data['bandwidths']):.2f} MB/s")
        
        # Print overall time statistics
        if self.overall_times:
            print("\n===== OVERALL TIME PROFILING =====\n")
            overall_df = pd.DataFrame({
                'Operation': list(self.overall_times.keys()),
                'Time (s)': list(self.overall_times.values())
            })
            overall_df = overall_df.sort_values('Time (s)', ascending=False)
            print(overall_df.to_string(index=False))
    
    def plot_kernel_times(self):
        # Create the plot
        plt.figure(figsize=(12, 8))
        
        # Sort kernels by total execution time
        sorted_kernels = sorted(self.kernel_times.keys(), 
                              key=lambda k: sum(self.kernel_times[k]), 
                              reverse=True)
        
        # Prepare data
        total_times = [sum(self.kernel_times[k]) for k in sorted_kernels]
        avg_times = [np.mean(self.kernel_times[k]) for k in sorted_kernels]
        calls = [self.kernel_calls[k] for k in sorted_kernels]
        
        # Create subplots
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14))
        
        # Plot total times
        ax1.barh(sorted_kernels, total_times)
        ax1.set_title('Total Kernel Execution Time (ms)')
        ax1.set_xlabel('Time (ms)')
        
        # Plot average times
        ax2.barh(sorted_kernels, avg_times)
        ax2.set_title('Average Kernel Execution Time (ms)')
        ax2.set_xlabel('Time (ms)')
        
        # Plot number of calls
        ax3.barh(sorted_kernels, calls)
        ax3.set_title('Number of Kernel Calls')
        ax3.set_xlabel('Calls')
        
        plt.tight_layout()
        plt.savefig('kernel_profile.png')
        print("Kernel profile plot saved to 'kernel_profile.png'")

# Create the profiling monitor
profiler = ProfilingMonitor()

# Store original methods to patch
original_methods = {}

def patch_gridff_cl():
    # Store original methods
    original_methods['make_MorseFF'] = gff.clgff.make_MorseFF
    original_methods['fit3D'] = gff.clgff.fit3D
    original_methods['poisson_old'] = gff.clgff.poisson_old
    original_methods['project_atoms_on_grid_quintic_pbc'] = gff.clgff.project_atoms_on_grid_quintic_pbc
    original_methods['laplace_real_loop_inert'] = gff.clgff.laplace_real_loop_inert
    original_methods['slabPotential'] = gff.clgff.slabPotential
    
    # Create profiled versions of the methods
    def profiled_make_MorseFF(self, atoms, REQs, **kwargs):
        # Track all timing for this function
        profiler.start_timer('make_MorseFF_total')
        
        # Extract common parameters with defaults
        Params = kwargs.get('Params', (0.00001, 1.5, 0.0, 0.0))
        nPBC = kwargs.get('nPBC', (0, 0, 0))
        # Handle either Ls or lvec (depending on which version is called)
        Ls = kwargs.get('Ls', (20.0, 20.0, 20.0))
        lvec = kwargs.get('lvec', None)
        g0 = kwargs.get('g0', None)  # Some calls include g0
        bDamp = kwargs.get('bDamp', True)
        bReturn = kwargs.get('bReturn', True)
        
        # Create a profiling-enabled queue
        original_queue = self.queue
        try:
            self.queue = cl.CommandQueue(self.ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
            
            # Prepare the call
            na = len(atoms)
            if isinstance(nPBC, tuple): nPBC = np.array(nPBC, dtype=np.int32)
            if isinstance(Ls, tuple): Ls = np.array(Ls, dtype=np.float32)
            if isinstance(Params, tuple): Params = np.array(Params, dtype=np.float32)
            buff_names = {'atoms', 'REQs', 'V_Paul', 'V_Lond'}
            
            # Profile buffer setup
            profiler.start_timer('buffer_setup')
            nxyz = self.gcl.nxyz
            self.try_make_buffs(buff_names, na, nxyz)
            profiler.stop_timer('buffer_setup')
            
            # Profile data transfer to device
            profiler.start_timer('data_transfer_h2d')
            atoms_np = np.array(atoms, dtype=np.float32)
            reqs_np = np.array(REQs, dtype=np.float32)
            
            # Profile the memory transfers
            atoms_event = cl.enqueue_copy(self.queue, self.atoms_buff, atoms_np, is_blocking=False)
            profiler.record_memory_transfer('Host to Device', atoms_np.nbytes, atoms_event)
            
            reqs_event = cl.enqueue_copy(self.queue, self.REQs_buff, reqs_np, is_blocking=False)
            profiler.record_memory_transfer('Host to Device', reqs_np.nbytes, reqs_event)
            profiler.stop_timer('data_transfer_h2d')
            
            # Profile kernel execution
            profiler.start_timer('kernel_execution_MorseFF')
            nL = self.nloc
            nG = self.gcl.nxyz
            
            # Call the kernel with appropriate arguments
            kernel_args = [
                self.queue, (nG,), (nL,),
                np.int32(na),
                self.atoms_buff,
                self.REQs_buff,
                Params,
                self.V_Paul_buff,
                self.V_Lond_buff,
                nPBC,
                Ls if lvec is None else lvec,  # Use lvec if provided, otherwise Ls
                np.int32(bDamp)
            ]
            
            # Execute the kernel and profile it
            kernel_event = self.prg.make_MorseFF(*kernel_args)
            profiler.record_kernel('make_MorseFF', kernel_event)
            profiler.stop_timer('kernel_execution_MorseFF')
            
            # Handle result retrieval if requested
            result = None
            if bReturn:
                profiler.start_timer('data_transfer_d2h')
                V_Paul = np.zeros(self.gsh.ns[::-1], dtype=np.float32)
                V_Lond = np.zeros(self.gsh.ns[::-1], dtype=np.float32)
                
                # Start timing (for compatibility with original code)
                t0 = time.time()
                prep_t = 0
                
                # Transfer results back to host with profiling
                paul_event = cl.enqueue_copy(self.queue, V_Paul, self.V_Paul_buff, is_blocking=False)
                profiler.record_memory_transfer('Device to Host', V_Paul.nbytes, paul_event)
                
                lond_event = cl.enqueue_copy(self.queue, V_Lond, self.V_Lond_buff, is_blocking=False)
                profiler.record_memory_transfer('Device to Host', V_Lond.nbytes, lond_event)
                
                result = (V_Paul, V_Lond)
                profiler.stop_timer('data_transfer_d2h')
                
                # Print timing information (like the original method)
                dt = time.time() - t0
                nabc = na * nG * np.prod(nPBC)
                nGOPs = nabc * 1e-9
                speed = nGOPs / dt
                print(f"GridFF_cl::make_MorseFF() time[s]  {dt}   preparation[s]   {prep_t} [s] nGOPs:  {nGOPs}  speed[GOPs/s]:  {speed}  na,nxyz,npbc  {na} {nG} {np.prod(nPBC)}")
        finally:
            # Always restore the original queue
            self.queue = original_queue
        
        profiler.stop_timer('make_MorseFF_total')
        return result
    
    # Patch the methods
    gff.clgff.make_MorseFF = profiled_make_MorseFF
    
    # Setup profiling command queue
    original_queue = gff.clgff.queue
    gff.clgff.queue = cl.CommandQueue(gff.clgff.ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)

def restore_gridff_cl():
    # Restore original methods
    for name, method in original_methods.items():
        setattr(gff.clgff, name, method)

def profile_gridff(fname="data/xyz/NaCl_1x1_L1.xyz", job="PLQ"):
    # Apply the patches
    patch_gridff_cl()
    
    # Run the code with profiling
    try:
        start_time = time.time()
        
        # Execute the code
        gff.test_gridFF_ocl(
            fname=fname,
            Element_Types_name="./data/ElementTypes.dat",
            save_name="profiled",
            job=job
        )
        
        total_time = time.time() - start_time
        print(f"\n===== TOTAL EXECUTION TIME: {total_time:.4f}s =====\n")
        
        # Print profiling statistics
        profiler.print_stats()
        
        # Generate profiling visualization
        profiler.plot_kernel_times()
        
    finally:
        # Restore original methods
        restore_gridff_cl()

if __name__ == "__main__":
    # Get command line arguments if provided
    if len(sys.argv) > 1:
        name = sys.argv[1]
    else:
        name = "NaCl_1x1_L1"
        
    if len(sys.argv) > 2:
        job = sys.argv[2]
    else:
        job = "PLQ"
    
    print(f"\n===== PROFILING ENABLED: {name}, {job} =====\n")
    profile_gridff(f"data/xyz/{name}.xyz", job)
