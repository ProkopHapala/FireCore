import sys
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import csv

sys.path.append("../../")

# Import the profiled GridFF_cl class
from pyBall.OCL.GridFF_profiled import ProfilingGridFF_cl

# Import the original test module
from pyBall.tests import ocl_GridFF_new as gff

# Create a profiling class to track overall execution time
class Profiler:
    def __init__(self):
        self.timings = {}
        self.current_timers = {}
    
    def start(self, name):
        self.current_timers[name] = time.time()
    
    def stop(self, name):
        if name in self.current_timers:
            elapsed = (time.time() - self.current_timers[name]) * 1000  # ms
            if name not in self.timings:
                self.timings[name] = []
            self.timings[name].append(elapsed)
            del self.current_timers[name]
            return elapsed
        return 0
    
    def get_results(self):
        results = {}
        for name, times in self.timings.items():
            results[name] = {
                'count': len(times),
                'total': sum(times),
                'average': sum(times) / len(times) if times else 0,
                'min': min(times) if times else 0,
                'max': max(times) if times else 0
            }
        return results

# Create a global profiler
profiler = Profiler()

def run_profiling(name="NaCl_1x1_L1", job="PLQ", iterations=1):
    print(f"\nRunning profiling for {job} job on {name} dataset ({iterations} iterations)...")
    
    # Create a profiling instance of GridFF_cl
    prof_clgff = ProfilingGridFF_cl()
    
    # Save the original GridFF_cl instance
    original_clgff = gff.clgff
    
    # Replace with our profiling instance
    gff.clgff = prof_clgff
    
    try:
        # Run multiple iterations for more accurate profiling
        for i in range(iterations):
            print(f"\nIteration {i+1}/{iterations}...")
            
            # Profile the overall execution
            profiler.start(f"total_{job}")
            
            # Profile data loading
            profiler.start("data_loading")
            atoms = gff.make_atoms_arrays(
                fname=f"data/xyz/{name}.xyz",
                Element_Types_name="./data/ElementTypes.dat"
            )
            profiler.stop("data_loading")
            
            # Run the test
            profiler.start(f"test_{job}")
            if job == "PLQ":
                gff.test_gridFF_ocl(
                    fname=f"data/xyz/{name}.xyz",
                    Element_Types_name="./data/ElementTypes.dat",
                    job=job,
                    bPlot=False  # Disable plotting for profiling
                )
            elif job == "Ewald":
                # For Ewald test, we need a simple system
                d = 0.6
                apos = np.array([
                    [-d, 0.0, 0.0],
                    [+d, 0.0, 0.0],
                    [0.0, -d, 0.0],
                    [0.0, +d, 0.0],
                ])
                qs = [+1.0, +1.0, -1.0, -1.0]
                gff.test_Ewald(
                    apos, qs,
                    ns=(100, 100, 100),
                    dg=(0.10, 0.10, 0.10),
                    order=3,
                    bPlot=False
                )
            profiler.stop(f"test_{job}")
            
            profiler.stop(f"total_{job}")
    finally:
        # Restore the original GridFF_cl instance
        gff.clgff = original_clgff
    
    # Get profiling results
    overall_results = profiler.get_results()
    gpu_results = prof_clgff.get_profiling_results()
    
    # Combine results
    combined_results = {
        'overall': overall_results,
        'gpu': gpu_results
    }
    
    # Print profiling results
    print("\n===== Overall Profiling Results =====")
    print(f"{'Operation':<30} {'Count':<8} {'Total (ms)':<12} {'Avg (ms)':<10}")
    print("-" * 70)
    
    sorted_overall = sorted(overall_results.items(), key=lambda x: x[1]['total'], reverse=True)
    for name, stats in sorted_overall:
        print(f"{name:<30} {stats['count']:<8} {stats['total']:<12.2f} {stats['average']:<10.2f}")
    
    # Print GPU profiling results
    prof_clgff.print_profiling_results()
    
    # Save results to CSV
    with open(f"profile_{job}_{name}.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Category', 'Operation', 'Count', 'Total (ms)', 'Avg (ms)', 'Min (ms)', 'Max (ms)'])
        
        # Write overall results
        for name, stats in sorted_overall:
            writer.writerow(['Overall', name, stats['count'], f"{stats['total']:.2f}", 
                            f"{stats['average']:.2f}", f"{stats['min']:.2f}", f"{stats['max']:.2f}"])
        
        # Write CPU results
        sorted_cpu = sorted(gpu_results['cpu_timings'].items(), key=lambda x: x[1]['total'], reverse=True)
        for name, stats in sorted_cpu:
            writer.writerow(['CPU', name, stats['count'], f"{stats['total']:.2f}", 
                            f"{stats['average']:.2f}", f"{stats['min']:.2f}", f"{stats['max']:.2f}"])
        
        # Write GPU results
        sorted_gpu = sorted(gpu_results['gpu_timings'].items(), key=lambda x: x[1]['total'], reverse=True)
        for name, stats in sorted_gpu:
            writer.writerow(['GPU', name, stats['count'], f"{stats['total']:.2f}", 
                            f"{stats['average']:.2f}", f"{stats['min']:.2f}", f"{stats['max']:.2f}"])
    
    print(f"\nDetailed results saved to profile_{job}_{name}.csv")
    
    # Create visualizations
    fig = plt.figure(figsize=(15, 12))
    
    # Plot overall timings
    plt.subplot(2, 2, 1)
    names = [name for name, _ in sorted_overall[:8]]
    times = [stats['total'] for _, stats in sorted_overall[:8]]
    
    y_pos = range(len(names))
    plt.barh(y_pos, times)
    plt.yticks(y_pos, names)
    plt.xlabel('Time (ms)')
    plt.title('Top Overall Operations')
    
    # Plot CPU timings
    plt.subplot(2, 2, 2)
    sorted_cpu = sorted(gpu_results['cpu_timings'].items(), key=lambda x: x[1]['total'], reverse=True)[:8]
    names = [name for name, _ in sorted_cpu]
    times = [stats['total'] for _, stats in sorted_cpu]
    
    y_pos = range(len(names))
    plt.barh(y_pos, times)
    plt.yticks(y_pos, names)
    plt.xlabel('Time (ms)')
    plt.title('Top CPU Operations')
    
    # Plot GPU timings if available
    if gpu_results['gpu_timings']:
        plt.subplot(2, 2, 3)
        sorted_gpu = sorted(gpu_results['gpu_timings'].items(), key=lambda x: x[1]['total'], reverse=True)[:8]
        names = [name for name, _ in sorted_gpu]
        times = [stats['total'] for _, stats in sorted_gpu]
        
        y_pos = range(len(names))
        plt.barh(y_pos, times)
        plt.yticks(y_pos, names)
        plt.xlabel('Time (ms)')
        plt.title('Top GPU Kernels')
    
    # Plot CPU vs GPU time distribution
    plt.subplot(2, 2, 4)
    cpu_time = sum(stats['total'] for _, stats in sorted_cpu)
    gpu_time = sum(stats['total'] for _, stats in sorted_gpu) if gpu_results['gpu_timings'] else 0
    
    labels = ['CPU Operations', 'GPU Operations']
    sizes = [cpu_time, gpu_time]
    
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
    plt.axis('equal')
    plt.title('CPU vs GPU Time Distribution')
    
    plt.tight_layout()
    plt.savefig(f"profile_{job}_{name}.png")
    print(f"Visualization saved to profile_{job}_{name}.png")
    
    return combined_results

def main():
    # Run profiling for PLQ job
    run_profiling(name="NaCl_1x1_L1", job="PLQ", iterations=1)
    
    # Uncomment to run Ewald profiling
    # run_profiling(name="NaCl_1x1_L1", job="Ewald", iterations=1)
    
    # Show plots
    plt.show()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error during profiling: {e}")
        import traceback
        traceback.print_exc()
