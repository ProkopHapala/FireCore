import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import pyopencl as cl

sys.path.append("../../")

from pyBall.tests import ocl_GridFF_new as gff

# Enable OpenCL profiling
os.environ['PYOPENCL_CTX'] = '0'
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'

# Create a context with profiling enabled
platforms = cl.get_platforms()
devices = platforms[0].get_devices(device_type=cl.device_type.GPU)
ctx = cl.Context(devices=devices)
queue = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)

# Override the default queue in GridFF_cl
gff.clgff.ctx = ctx
gff.clgff.queue = queue

# Function to profile kernel execution
def profile_kernel(name, func, *args, **kwargs):
    # CPU timing
    cpu_start = time.time()
    result = func(*args, **kwargs)
    cpu_end = time.time()
    cpu_time = (cpu_end - cpu_start) * 1000  # Convert to ms
    
    print(f"\n{name} CPU time: {cpu_time:.2f} ms")
    
    # Get OpenCL event timing if available
    if hasattr(gff.clgff, 'last_events') and gff.clgff.last_events:
        for i, event in enumerate(gff.clgff.last_events):
            start = event.get_profiling_info(cl.profiling_info.START)
            end = event.get_profiling_info(cl.profiling_info.END)
            gpu_time = (end - start) * 1e-6  # Convert from ns to ms
            print(f"{name} Kernel {i} GPU time: {gpu_time:.2f} ms")
    
    return result, cpu_time

# Profile different operations
def run_profiling(name="NaCl_1x1_L1", job="PLQ"):
    # Prepare the profiling results container
    profile_results = {
        'operation': [],
        'cpu_time_ms': [],
        'gpu_time_ms': []
    }
    
    print(f"\n===== Profiling GridFF OpenCL with {name} dataset =====\n")
    
    # Modify the GridFF_cl class to track events
    if not hasattr(gff.clgff, 'last_events'):
        gff.clgff.last_events = []
    
    # Add event tracking to the enqueue_kernel method
    original_enqueue = gff.clgff.enqueue_kernel
    
    def enqueue_with_profiling(*args, **kwargs):
        event = original_enqueue(*args, **kwargs)
        gff.clgff.last_events.append(event)
        return event
    
    # Monkey patch the method
    gff.clgff.enqueue_kernel = enqueue_with_profiling
    
    # Run the test with profiling
    try:
        # Profile the main test function
        _, cpu_time = profile_kernel(
            "test_gridFF_ocl",
            gff.test_gridFF_ocl,
            fname=f"data/xyz/{name}.xyz",
            Element_Types_name="./data/ElementTypes.dat",
            job=job,
            bPlot=False  # Disable plotting for profiling
        )
        
        profile_results['operation'].append(f"Full {job} calculation")
        profile_results['cpu_time_ms'].append(cpu_time)
        
        # Reset for individual kernel profiling
        gff.clgff.last_events = []
        
        # Profile individual operations
        # 1. Create a sample system
        atoms = gff.make_atoms_arrays(fname=f"data/xyz/{name}.xyz", 
                                     Element_Types_name="./data/ElementTypes.dat")
        
        # 2. Profile grid setup
        gff.clgff.last_events = []
        _, cpu_time = profile_kernel("Grid setup", gff.clgff.set_grid, 
                                   gff.GridShape(dg=(0.1, 0.1, 0.1), 
                                                Ls=[10.0, 10.0, 10.0], 
                                                g0=(-5.0, -5.0, -5.0)))
        profile_results['operation'].append("Grid setup")
        profile_results['cpu_time_ms'].append(cpu_time)
        
        # 3. Profile Coulomb calculation
        gff.clgff.last_events = []
        xyzq = np.zeros((len(atoms), 4), dtype=np.float32)
        xyzq[:, 0:3] = atoms[:, 0:3]
        xyzq[:, 3] = atoms[:, 3]  # Charges
        _, cpu_time = profile_kernel("Coulomb calculation", gff.clgff.makeCoulombEwald, xyzq)
        profile_results['operation'].append("Coulomb calculation")
        profile_results['cpu_time_ms'].append(cpu_time)
        
        # 4. Profile point evaluation
        gff.clgff.last_events = []
        ps = np.zeros((100, 4), dtype=np.float32)
        ps[:, 0] = np.linspace(-5, 5, 100)
        ps[:, 1] = np.linspace(-5, 5, 100)
        ps[:, 2] = 0.0
        _, cpu_time = profile_kernel("Point evaluation", gff.clgff.make_Coulomb_points, 
                                   atoms, ps, nPBC=(10, 10, 0), 
                                   Ls=[10.0, 10.0, 10.0], 
                                   GFFParams=(0.00001, 1.5, 0.0, 0.0), 
                                   bReturn=True)
        profile_results['operation'].append("Point evaluation")
        profile_results['cpu_time_ms'].append(cpu_time)
        
        # Print summary table
        print("\n===== Profiling Summary =====")
        print(f"{'Operation':<30} {'CPU Time (ms)':<15} {'GPU Time (ms)':<15}")
        print("-" * 60)
        for i, op in enumerate(profile_results['operation']):
            gpu_time = "N/A"
            print(f"{op:<30} {profile_results['cpu_time_ms'][i]:<15.2f} {gpu_time:<15}")
        
    finally:
        # Restore original method
        gff.clgff.enqueue_kernel = original_enqueue

if __name__ == "__main__":
    # Run profiling with default parameters
    run_profiling()
    
    # You can also profile different jobs or datasets
    # run_profiling(name="NaCl_8x8_L3", job="Morse")
    # run_profiling(job="Ewald")
