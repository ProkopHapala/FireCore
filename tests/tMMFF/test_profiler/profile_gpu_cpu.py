import sys
import numpy as np
import os
import time
import matplotlib.pyplot as plt

sys.path.append("../../")

# Import the modules we need to profile
from pyBall.tests import ocl_GridFF_new as gff

# Dictionary to store timing information
timings = {}

# Function to start timing
def start_timing(name):
    if name not in timings:
        timings[name] = {}
    timings[name]['start'] = time.time()

# Function to stop timing
def stop_timing(name):
    if name in timings and 'start' in timings[name]:
        timings[name]['end'] = time.time()
        timings[name]['elapsed'] = (timings[name]['end'] - timings[name]['start']) * 1000  # ms
        return timings[name]['elapsed']
    return 0

# Function to profile the main test function
def profile_test_gridFF_ocl(name="NaCl_1x1_L1", job="PLQ"):
    # Start overall timing
    start_timing('total')
    
    # Profile data loading and preparation
    start_timing('data_preparation')
    atoms = gff.make_atoms_arrays(
        fname=f"data/xyz/{name}.xyz",
        Element_Types_name="./data/ElementTypes.dat"
    )
    stop_timing('data_preparation')
    
    # Store the original methods we need to profile
    original_makeCoulombEwald = gff.clgff.makeCoulombEwald
    original_make_Coulomb_points = gff.clgff.make_Coulomb_points
    
    # Create profiling wrappers
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
    
    try:
        # Profile the main test function
        start_timing('test_gridFF_ocl')
        gff.test_gridFF_ocl(
            fname=f"data/xyz/{name}.xyz",
            Element_Types_name="./data/ElementTypes.dat",
            job=job,
            bPlot=False  # Disable plotting for profiling
        )
        stop_timing('test_gridFF_ocl')
    finally:
        # Restore original methods
        gff.clgff.makeCoulombEwald = original_makeCoulombEwald
        gff.clgff.make_Coulomb_points = original_make_Coulomb_points
    
    # Stop overall timing
    stop_timing('total')
    
    # Return the timing results
    return timings

# Function to profile the Ewald test
def profile_test_Ewald():
    # Start overall timing
    start_timing('total_Ewald')
    
    # Create a simple system for Ewald test
    d = 0.6
    apos = np.array([
        [-d, 0.0, 0.0],
        [+d, 0.0, 0.0],
        [0.0, -d, 0.0],
        [0.0, +d, 0.0],
    ])
    qs = [+1.0, +1.0, -1.0, -1.0]
    
    # Store the original methods we need to profile
    original_makeCoulombEwald = gff.clgff.makeCoulombEwald
    
    # Create profiling wrapper
    def profiled_makeCoulombEwald(xyzq, bOld=False):
        start_timing('GPU_makeCoulombEwald_Ewald')
        result = original_makeCoulombEwald(xyzq, bOld)
        stop_timing('GPU_makeCoulombEwald_Ewald')
        return result
    
    # Apply the monkey patch
    gff.clgff.makeCoulombEwald = profiled_makeCoulombEwald
    
    try:
        # Profile the Ewald test
        start_timing('test_Ewald')
        gff.test_Ewald(
            apos, qs,
            ns=(100, 100, 100),
            dg=(0.10, 0.10, 0.10),
            order=3,
            bPlot=False
        )
        stop_timing('test_Ewald')
    finally:
        # Restore original method
        gff.clgff.makeCoulombEwald = original_makeCoulombEwald
    
    # Stop overall timing
    stop_timing('total_Ewald')
    
    # Return the timing results
    return timings

# Function to display profiling results
def display_profiling_results(results, title):
    print(f"\n===== {title} =====\n")
    
    # Extract elapsed times
    elapsed_times = {}
    for name, data in results.items():
        if 'elapsed' in data:
            elapsed_times[name] = data['elapsed']
    
    # Sort by elapsed time
    sorted_times = sorted(elapsed_times.items(), key=lambda x: x[1], reverse=True)
    
    # Calculate total time
    if 'total' in elapsed_times:
        total_time = elapsed_times['total']
    else:
        total_time = sum(elapsed_times.values())
    
    # Print results
    print(f"{'Operation':<30} {'Time (ms)':<15} {'% of Total':<15}")
    print("-" * 60)
    for name, time_ms in sorted_times:
        percent = (time_ms / total_time) * 100
        print(f"{name:<30} {time_ms:<15.2f} {percent:<15.2f}%")
    
    # Create visualization
    plt.figure(figsize=(12, 10))
    
    # Bar chart of execution times
    plt.subplot(2, 1, 1)
    operations = [k for k, _ in sorted_times[:8]]  # Top 8 operations
    times = [t for _, t in sorted_times[:8]]
    
    y_pos = range(len(operations))
    plt.barh(y_pos, times)
    plt.yticks(y_pos, operations)
    plt.xlabel('Time (ms)')
    plt.title(f'Top Operations by Execution Time - {title}')
    
    # Pie chart of GPU vs CPU time
    plt.subplot(2, 1, 2)
    gpu_time = sum(t for n, t in sorted_times if n.startswith('GPU_'))
    cpu_time = total_time - gpu_time
    
    labels = ['GPU Operations', 'CPU Operations']
    sizes = [gpu_time, cpu_time]
    
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
    plt.axis('equal')
    plt.title('GPU vs CPU Time Distribution')
    
    plt.tight_layout()
    plt.savefig(f"profile_results_{title.replace(' ', '_')}.png")
    print(f"\nVisualization saved to profile_results_{title.replace(' ', '_')}.png")

# Main function
def main():
    # Profile the PLQ job
    print("\nProfiling PLQ job...")
    results_plq = profile_test_gridFF_ocl(name="NaCl_1x1_L1", job="PLQ")
    display_profiling_results(results_plq, "PLQ Job")
    
    # Profile the Ewald test
    print("\nProfiling Ewald test...")
    results_ewald = profile_test_Ewald()
    display_profiling_results(results_ewald, "Ewald Test")
    
    # Show plots
    plt.show()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error during profiling: {e}")
        import traceback
        traceback.print_exc()
