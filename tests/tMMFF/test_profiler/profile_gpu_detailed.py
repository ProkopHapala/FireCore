import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import pyopencl as cl

sys.path.append("../../")

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

# Start overall timing
start_timing('total')

# Import the necessary modules
start_timing('import')
from pyBall.tests import ocl_GridFF_new as gff
stop_timing('import')

# Profile OpenCL initialization and device selection
start_timing('opencl_init')
try:
    # Get OpenCL platforms and devices
    platforms = cl.get_platforms()
    nvidia_platform = None
    nvidia_device = None
    
    # Look for NVIDIA platform
    for platform in platforms:
        if 'NVIDIA' in platform.name:
            nvidia_platform = platform
            break
    
    # If NVIDIA platform found, get NVIDIA device
    if nvidia_platform:
        devices = nvidia_platform.get_devices()
        if devices:
            nvidia_device = devices[0]
            print(f"Selected device: {nvidia_device.name}")
            print(f"Device Name: {nvidia_device.name}")
            print(f"Max Compute Units: {nvidia_device.max_compute_units}")
            print(f"Max Work Group Size: {nvidia_device.max_work_group_size}")
            print(f"Global Memory Size: {nvidia_device.global_mem_size / (1024**2)} MB")
            print(f"Local Memory Size: {nvidia_device.local_mem_size / 1024} KB")
            print(f"Max Clock Frequency: {nvidia_device.max_clock_frequency} MHz")
            
            # Create context and queue with profiling enabled
            ctx = cl.Context([nvidia_device])
            queue = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
            
            # Try to set the context and queue in the GridFF_cl instance
            if hasattr(gff, 'clgff') and gff.clgff is not None:
                gff.clgff.ctx = ctx
                gff.clgff.queue = queue
    else:
        # If no NVIDIA platform found, use the first available platform and device
        if platforms:
            devices = platforms[0].get_devices()
            if devices:
                print(f"OpenCL Platform: {platforms[0].name}")
                print(f"OpenCL Device: {devices[0].name}")
                print(f"Compute Units: {devices[0].max_compute_units}")
                
                # Create context and queue with profiling enabled
                ctx = cl.Context([devices[0]])
                queue = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
                
                # Try to set the context and queue in the GridFF_cl instance
                if hasattr(gff, 'clgff') and gff.clgff is not None:
                    gff.clgff.ctx = ctx
                    gff.clgff.queue = queue
except Exception as e:
    print(f"Error initializing OpenCL: {e}")
stop_timing('opencl_init')

# Profile data loading and preparation
start_timing('data_preparation')
name = "NaCl_1x1_L1"
job = "PLQ"
try:
    atoms = gff.make_atoms_arrays(
        fname=f"data/xyz/{name}.xyz",
        Element_Types_name="./data/ElementTypes.dat"
    )
    print(f"Loaded {len(atoms)} atoms from {name}")
except Exception as e:
    print(f"Error loading data: {e}")
stop_timing('data_preparation')

# Monkey patch the GridFF_cl methods to add profiling
start_timing('monkey_patching')
if hasattr(gff, 'clgff') and gff.clgff is not None:
    # Save original methods
    original_enqueue_kernel = gff.clgff.enqueue_kernel
    original_makeCoulombEwald = gff.clgff.makeCoulombEwald
    original_make_Coulomb_points = gff.clgff.make_Coulomb_points
    
    # Add profiling to enqueue_kernel
    def profiled_enqueue_kernel(kernel, global_size, local_size=None, args=None):
        kernel_name = kernel.function_name
        start_timing(f"kernel_{kernel_name}")
        event = original_enqueue_kernel(kernel, global_size, local_size, args)
        
        # Get profiling info if available
        if event is not None:
            try:
                event.wait()
                start = event.get_profiling_info(cl.profiling_info.START)
                end = event.get_profiling_info(cl.profiling_info.END)
                gpu_time = (end - start) * 1e-6  # Convert nanoseconds to milliseconds
                record_kernel_time(kernel_name, gpu_time)
            except:
                pass
        
        stop_timing(f"kernel_{kernel_name}")
        return event
    
    # Add profiling to makeCoulombEwald
    def profiled_makeCoulombEwald(xyzq, bOld=False):
        start_timing('GPU_makeCoulombEwald')
        result = original_makeCoulombEwald(xyzq, bOld)
        stop_timing('GPU_makeCoulombEwald')
        return result
    
    # Add profiling to make_Coulomb_points
    def profiled_make_Coulomb_points(atoms, ps, nPBC=(0,0,0), Ls=None, GFFParams=(0.0001, 1.5, 0.0, 0.0), bReturn=False):
        start_timing('GPU_make_Coulomb_points')
        result = original_make_Coulomb_points(atoms, ps, nPBC, Ls, GFFParams, bReturn)
        stop_timing('GPU_make_Coulomb_points')
        return result
    
    # Apply the monkey patches
    gff.clgff.enqueue_kernel = profiled_enqueue_kernel
    gff.clgff.makeCoulombEwald = profiled_makeCoulombEwald
    gff.clgff.make_Coulomb_points = profiled_make_Coulomb_points
stop_timing('monkey_patching')

# Profile the main test function
start_timing('test_execution')
try:
    # Run the test
    gff.test_gridFF_ocl(
        fname=f"data/xyz/{name}.xyz",
        Element_Types_name="./data/ElementTypes.dat",
        job=job
    )
    
    # Restore original methods if patched
    if hasattr(gff, 'clgff') and gff.clgff is not None:
        if 'original_enqueue_kernel' in locals():
            gff.clgff.enqueue_kernel = original_enqueue_kernel
        if 'original_makeCoulombEwald' in locals():
            gff.clgff.makeCoulombEwald = original_makeCoulombEwald
        if 'original_make_Coulomb_points' in locals():
            gff.clgff.make_Coulomb_points = original_make_Coulomb_points
except Exception as e:
    print(f"Error running test: {e}")
stop_timing('test_execution')

# Stop overall timing
stop_timing('total')

# Print profiling results
print("\n===== Profiling Results =====")
for name, data in sorted(timings.items(), key=lambda x: x[1].get('elapsed', 0), reverse=True):
    if 'elapsed' in data:
        print(f"{name}: {data['elapsed']:.2f} ms")

# Print kernel profiling results if available
if kernel_timings:
    print("\n===== GPU Kernel Profiling Results =====")
    for kernel_name, times in sorted(kernel_timings.items(), key=lambda x: sum(x[1]), reverse=True):
        total_time = sum(times)
        avg_time = total_time / len(times) if times else 0
        print(f"{kernel_name}: {total_time:.2f} ms total, {len(times)} calls, {avg_time:.2f} ms avg")

# Create visualization
plt.figure(figsize=(15, 10))

# Extract elapsed times
elapsed_times = {}
for name, data in timings.items():
    if 'elapsed' in data:
        elapsed_times[name] = data['elapsed']

# Sort by elapsed time
sorted_times = sorted(elapsed_times.items(), key=lambda x: x[1], reverse=True)

# Bar chart of execution times (excluding very small times)
significant_times = [(k, v) for k, v in sorted_times if v > 1.0]  # Only show times > 1ms

plt.subplot(2, 1, 1)
operations = [k for k, _ in significant_times]
times = [t for _, t in significant_times]

y_pos = range(len(operations))
plt.barh(y_pos, times)
plt.yticks(y_pos, operations)
plt.xlabel('Time (ms)')
plt.title(f'Operation Execution Times - {job} job on {name} dataset')

# Add time labels
for i, v in enumerate(times):
    plt.text(v + 0.1, i, f"{v:.2f} ms", va='center')

# Plot kernel execution times if available
if kernel_timings:
    plt.subplot(2, 1, 2)
    kernel_names = []
    kernel_times = []
    
    for kernel_name, times in sorted(kernel_timings.items(), key=lambda x: sum(x[1]), reverse=True):
        kernel_names.append(kernel_name)
        kernel_times.append(sum(times))
    
    y_pos = range(len(kernel_names))
    plt.barh(y_pos, kernel_times)
    plt.yticks(y_pos, kernel_names)
    plt.xlabel('Time (ms)')
    plt.title('GPU Kernel Execution Times')
    
    # Add time labels
    for i, v in enumerate(kernel_times):
        plt.text(v + 0.1, i, f"{v:.2f} ms", va='center')

plt.tight_layout()
plt.savefig(f"profile_results_{job}_{name}_detailed.png")
print(f"\nVisualization saved to profile_results_{job}_{name}_detailed.png")

# Show plot
plt.show()
