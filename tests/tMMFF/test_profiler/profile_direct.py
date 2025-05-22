import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")

# Dictionary to store timing information
timings = {}

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

# Start overall timing
start_timing('total')

# Import the necessary modules
start_timing('import')
import pyopencl as cl
from pyBall.tests import ocl_GridFF_new as gff
stop_timing('import')

# Fix for OpenCL initialization
start_timing('opencl_init')
try:
    # Try to initialize OpenCL manually
    platforms = cl.get_platforms()
    if platforms:
        devices = platforms[0].get_devices()
        if devices:
            ctx = cl.Context(devices=devices)
            queue = cl.CommandQueue(ctx)
            
            # Print OpenCL info
            print("\nOpenCL Platform:", platforms[0].name)
            print("OpenCL Device:", devices[0].name)
            print("Compute Units:", devices[0].max_compute_units)
            
            # Try to set the context and queue in the GridFF_cl instance
            if hasattr(gff, 'clgff') and gff.clgff is not None:
                gff.clgff.ctx = ctx
                gff.clgff.queue = queue
                print("Successfully set OpenCL context and queue in GridFF_cl instance")
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

# Profile the main test function
start_timing('test_execution')
try:
    # Monkey patch the GridFF_cl methods to add profiling
    if hasattr(gff, 'clgff') and gff.clgff is not None:
        original_makeCoulombEwald = gff.clgff.makeCoulombEwald
        original_make_Coulomb_points = gff.clgff.make_Coulomb_points
        
        def profiled_makeCoulombEwald(xyzq, bOld=False):
            start_timing('GPU_makeCoulombEwald')
            result = original_makeCoulombEwald(xyzq, bOld)
            stop_timing('GPU_makeCoulombEwald')
            return result
        
        def profiled_make_Coulomb_points(atoms, ps, nPBC=(0,0,0), Ls=None, GFFParams=(0.0001, 1.5, 0.0, 0.0), bReturn=False):
            start_timing('GPU_make_Coulomb_points')
            result = original_make_Coulomb_points(atoms, ps, nPBC, Ls, GFFParams, bReturn)
            stop_timing('GPU_make_Coulomb_points')
            return result
        
        # Apply the monkey patches
        gff.clgff.makeCoulombEwald = profiled_makeCoulombEwald
        gff.clgff.make_Coulomb_points = profiled_make_Coulomb_points
    
    # Run the test
    gff.test_gridFF_ocl(
        fname=f"data/xyz/{name}.xyz",
        Element_Types_name="./data/ElementTypes.dat",
        job=job,
        bPlot=False  # Disable plotting for profiling
    )
    
    # Restore original methods if patched
    if hasattr(gff, 'clgff') and gff.clgff is not None:
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

# Create visualization
plt.figure(figsize=(12, 8))

# Extract elapsed times
elapsed_times = {}
for name, data in timings.items():
    if 'elapsed' in data:
        elapsed_times[name] = data['elapsed']

# Sort by elapsed time
sorted_times = sorted(elapsed_times.items(), key=lambda x: x[1], reverse=True)

# Bar chart of execution times
operations = [k for k, _ in sorted_times]
times = [t for _, t in sorted_times]

y_pos = range(len(operations))
plt.barh(y_pos, times)
plt.yticks(y_pos, operations)
plt.xlabel('Time (ms)')
plt.title(f'Operation Execution Times - {job} job on {name} dataset')

# Add time labels
for i, v in enumerate(times):
    plt.text(v + 0.1, i, f"{v:.2f} ms", va='center')

plt.tight_layout()
plt.savefig(f"profile_results_{job}_{name}.png")
print(f"\nVisualization saved to profile_results_{job}_{name}.png")

# Show plot
plt.show()
