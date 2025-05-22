import sys
import numpy as np
import os
import time
import matplotlib.pyplot as plt

sys.path.append("../../")

# Set environment variables for OpenCL
os.environ['PYOPENCL_CTX'] = '0'

# Function to measure execution time
def measure_execution_time(func, *args, **kwargs):
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    return result, (end_time - start_time) * 1000  # Convert to milliseconds

# Main profiling function
def profile_gridff_execution(name="NaCl_1x1_L1", job="PLQ", iterations=3):
    from pyBall.tests import ocl_GridFF_new as gff
    
    print(f"\n===== Profiling GridFF OpenCL with {name} dataset, {job} job =====\n")
    
    # Results containers
    timings = {
        'total': [],
        'setup': [],
        'computation': [],
        'teardown': []
    }
    
    # Run multiple iterations
    for i in range(iterations):
        print(f"\nRunning iteration {i+1}/{iterations}...")
        
        # Measure total execution time
        start_total = time.time()
        
        # Setup phase - load data and prepare
        start_setup = time.time()
        if job == "PLQ":
            # Prepare data without running the full test
            atoms = gff.make_atoms_arrays(
                fname=f"data/xyz/{name}.xyz",
                Element_Types_name="./data/ElementTypes.dat"
            )
        elif job == "Ewald":
            # For Ewald test, create a simple system
            d = 0.6
            atoms = np.array([
                [-d, 0.0, 0.0, 1.0],
                [+d, 0.0, 0.0, 1.0],
                [0.0, -d, 0.0, -1.0],
                [0.0, +d, 0.0, -1.0],
            ])
        end_setup = time.time()
        
        # Computation phase - run the actual OpenCL computation
        start_computation = time.time()
        if job == "PLQ":
            # Run the computation part only
            result, _ = measure_execution_time(
                gff.test_gridFF_ocl,
                fname=f"data/xyz/{name}.xyz",
                Element_Types_name="./data/ElementTypes.dat",
                job=job,
                bPlot=False  # Disable plotting for profiling
            )
        elif job == "Ewald":
            # Run Ewald computation
            apos = atoms[:, 0:3]
            qs = atoms[:, 3]
            result, _ = measure_execution_time(
                gff.test_Ewald,
                apos, qs,
                ns=(100, 100, 100),
                dg=(0.10, 0.10, 0.10),
                order=3,
                bPlot=False
            )
        end_computation = time.time()
        
        # Teardown phase - cleanup
        start_teardown = time.time()
        # No explicit teardown needed in this case
        end_teardown = time.time()
        
        # Record timings
        end_total = time.time()
        timings['total'].append((end_total - start_total) * 1000)
        timings['setup'].append((end_setup - start_setup) * 1000)
        timings['computation'].append((end_computation - start_computation) * 1000)
        timings['teardown'].append((end_teardown - start_teardown) * 1000)
    
    # Calculate average timings
    avg_timings = {}
    for phase, times in timings.items():
        avg_timings[phase] = sum(times) / len(times)
    
    # Print results
    print("\n===== Profiling Results =====")
    print(f"Average execution times over {iterations} iterations:")
    print(f"  Total:       {avg_timings['total']:.2f} ms")
    print(f"  Setup:       {avg_timings['setup']:.2f} ms ({avg_timings['setup']/avg_timings['total']*100:.1f}%)")
    print(f"  Computation: {avg_timings['computation']:.2f} ms ({avg_timings['computation']/avg_timings['total']*100:.1f}%)")
    print(f"  Teardown:    {avg_timings['teardown']:.2f} ms ({avg_timings['teardown']/avg_timings['total']*100:.1f}%)")
    
    # Create visualization
    plt.figure(figsize=(10, 6))
    
    # Bar chart of execution phases
    phases = ['Setup', 'Computation', 'Teardown']
    times = [avg_timings['setup'], avg_timings['computation'], avg_timings['teardown']]
    
    plt.bar(phases, times, color=['blue', 'green', 'red'])
    plt.ylabel('Time (ms)')
    plt.title(f'Execution Time Breakdown - {job} job on {name} dataset')
    
    # Add percentage labels
    for i, v in enumerate(times):
        plt.text(i, v + 5, f"{v:.1f} ms\n({v/avg_timings['total']*100:.1f}%)", 
                 ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(f"profile_results_{job}_{name}.png")
    print(f"Visualization saved to profile_results_{job}_{name}.png")
    
    return avg_timings

# Function to profile individual OpenCL kernels
def profile_opencl_kernels():
    import pyopencl as cl
    from pyBall.OCL.GridFF import GridFF_cl
    
    # Get OpenCL platform and device info
    platforms = cl.get_platforms()
    if not platforms:
        print("No OpenCL platforms found!")
        return
    
    print("\nOpenCL Platforms and Devices:")
    for i, platform in enumerate(platforms):
        print(f"  Platform {i}: {platform.name} ({platform.version})")
        devices = platform.get_devices()
        for j, device in enumerate(devices):
            print(f"    Device {j}: {device.name} (Type: {cl.device_type.to_string(device.type)})")
            print(f"      Compute Units: {device.max_compute_units}")
            print(f"      Global Memory: {device.global_mem_size / (1024**2):.2f} MB")
            print(f"      Local Memory: {device.local_mem_size / 1024:.2f} KB")
            print(f"      Max Work Group Size: {device.max_work_group_size}")

# Run the profiling
if __name__ == "__main__":
    try:
        # First, check OpenCL environment
        profile_opencl_kernels()
        
        # Run profiling for PLQ job
        profile_gridff_execution(name="NaCl_1x1_L1", job="PLQ", iterations=3)
        
        # Uncomment to run additional profiling jobs
        # profile_gridff_execution(name="NaCl_1x1_L1", job="Ewald", iterations=3)
        
        # Show plots
        plt.show()
    except Exception as e:
        print(f"Error during profiling: {e}")
        import traceback
        traceback.print_exc()
