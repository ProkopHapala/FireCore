import sys
import numpy as np
import os
import time
import matplotlib.pyplot as plt

sys.path.append("../../")

# Set environment variables for OpenCL
os.environ['PYOPENCL_CTX'] = '0'

# Simple timer class for manual profiling
class Timer:
    def __init__(self, name):
        self.name = name
        self.start_time = None
        self.end_time = None
    
    def start(self):
        self.start_time = time.time()
        return self
    
    def stop(self):
        self.end_time = time.time()
        return self
    
    def elapsed_ms(self):
        if self.start_time is None:
            return 0
        end = self.end_time if self.end_time is not None else time.time()
        return (end - self.start_time) * 1000  # Convert to milliseconds

# Profiler class to track multiple timers
class Profiler:
    def __init__(self):
        self.timers = {}
        self.current_timers = {}
    
    def start(self, name):
        if name not in self.current_timers:
            self.current_timers[name] = Timer(name).start()
        return self.current_timers[name]
    
    def stop(self, name):
        if name in self.current_timers:
            timer = self.current_timers[name].stop()
            if name not in self.timers:
                self.timers[name] = []
            self.timers[name].append(timer.elapsed_ms())
            del self.current_timers[name]
    
    def get_average(self, name):
        if name in self.timers and self.timers[name]:
            return sum(self.timers[name]) / len(self.timers[name])
        return 0
    
    def get_total(self, name):
        if name in self.timers and self.timers[name]:
            return sum(self.timers[name])
        return 0
    
    def get_count(self, name):
        if name in self.timers:
            return len(self.timers[name])
        return 0
    
    def get_summary(self):
        summary = {}
        for name in self.timers:
            summary[name] = {
                'count': self.get_count(name),
                'total': self.get_total(name),
                'average': self.get_average(name)
            }
        return summary

# Global profiler instance
profiler = Profiler()

# Function to run a test with profiling
def run_profiled_test(name="NaCl_1x1_L1", job="PLQ"):
    from pyBall.tests import ocl_GridFF_new as gff
    
    # Monkey patch the GridFF_cl class to add profiling
    from pyBall.OCL.GridFF import GridFF_cl
    
    # Save original methods
    original_enqueue_kernel = GridFF_cl.enqueue_kernel
    original_makeCoulombEwald = GridFF_cl.makeCoulombEwald
    original_make_Coulomb_points = GridFF_cl.make_Coulomb_points
    
    # Add profiling to enqueue_kernel
    def profiled_enqueue_kernel(self, kernel, global_size, local_size=None, args=None):
        kernel_name = kernel.function_name
        timer_name = f"kernel:{kernel_name}"
        profiler.start(timer_name)
        result = original_enqueue_kernel(self, kernel, global_size, local_size, args)
        profiler.stop(timer_name)
        return result
    
    # Add profiling to makeCoulombEwald
    def profiled_makeCoulombEwald(self, xyzq, bOld=False):
        profiler.start("makeCoulombEwald")
        result = original_makeCoulombEwald(self, xyzq, bOld)
        profiler.stop("makeCoulombEwald")
        return result
    
    # Add profiling to make_Coulomb_points
    def profiled_make_Coulomb_points(self, atoms, ps, nPBC=(0,0,0), Ls=None, GFFParams=(0.0001, 1.5, 0.0, 0.0), bReturn=False):
        profiler.start("make_Coulomb_points")
        result = original_make_Coulomb_points(self, atoms, ps, nPBC, Ls, GFFParams, bReturn)
        profiler.stop("make_Coulomb_points")
        return result
    
    # Apply the monkey patches
    GridFF_cl.enqueue_kernel = profiled_enqueue_kernel
    GridFF_cl.makeCoulombEwald = profiled_makeCoulombEwald
    GridFF_cl.make_Coulomb_points = profiled_make_Coulomb_points
    
    try:
        # Run the test with overall profiling
        profiler.start("total_execution")
        
        if job == "PLQ":
            gff.test_gridFF_ocl(
                fname=f"data/xyz/{name}.xyz",
                Element_Types_name="./data/ElementTypes.dat",
                job=job,
                bPlot=False  # Disable plotting for profiling
            )
        elif job == "Ewald":
            # For Ewald test, we need to create a simple system
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
        
        profiler.stop("total_execution")
        
    finally:
        # Restore original methods
        GridFF_cl.enqueue_kernel = original_enqueue_kernel
        GridFF_cl.makeCoulombEwald = original_makeCoulombEwald
        GridFF_cl.make_Coulomb_points = original_make_Coulomb_points
    
    # Return profiling results
    return profiler.get_summary()

# Function to display profiling results
def display_profiling_results(results, job, name):
    print(f"\n===== Profiling Results for {job} job on {name} dataset =====\n")
    
    # Sort results by total time
    sorted_results = sorted(results.items(), key=lambda x: x[1]['total'], reverse=True)
    
    # Print total execution time
    if "total_execution" in results:
        total_time = results["total_execution"]["total"]
        print(f"Total execution time: {total_time:.2f} ms")
    else:
        total_time = sum(item[1]['total'] for item in sorted_results)
        print(f"Sum of all profiled operations: {total_time:.2f} ms")
    
    # Print detailed results
    print("\nDetailed timing breakdown:")
    print(f"{'Operation':<40} {'Count':<8} {'Total (ms)':<12} {'Avg (ms)':<10} {'% of Total':<10}")
    print("-" * 80)
    
    for name, stats in sorted_results:
        percent = (stats['total'] / total_time) * 100 if total_time > 0 else 0
        print(f"{name:<40} {stats['count']:<8} {stats['total']:<12.2f} {stats['average']:<10.2f} {percent:<10.2f}%")
    
    # Create visualization
    plt.figure(figsize=(15, 10))
    
    # Pie chart of operation times
    plt.subplot(2, 1, 1)
    labels = [k for k, _ in sorted_results[:5] if not k.startswith("total_")]  # Top 5 operations excluding total
    sizes = [s['total'] for n, s in sorted_results[:5] if not n.startswith("total_")]
    
    if len(sorted_results) > 5:
        others_total = sum(s['total'] for n, s in sorted_results[5:] if not n.startswith("total_"))
        if others_total > 0:
            labels.append('Others')
            sizes.append(others_total)
    
    if sizes:  # Only create pie chart if we have data
        plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
        plt.axis('equal')
        plt.title('Operation Time Distribution')
    
    # Bar chart of operation times
    plt.subplot(2, 1, 2)
    operations = [k for k, _ in sorted_results[:10] if not k.startswith("total_")]  # Top 10 operations excluding total
    times = [s['total'] for n, s in sorted_results[:10] if not n.startswith("total_")]
    
    if operations:  # Only create bar chart if we have data
        y_pos = range(len(operations))
        plt.barh(y_pos, times)
        plt.yticks(y_pos, operations)
        plt.xlabel('Time (ms)')
        plt.title('Top 10 Operation Times')
    
    plt.tight_layout()
    plt.savefig(f"profile_results_{job}_{name}.png")
    print(f"\nVisualization saved to profile_results_{job}_{name}.png")

# Main function
def main():
    # Run profiling for PLQ job
    print("\nRunning profiling for PLQ job...")
    results_plq = run_profiled_test(name="NaCl_1x1_L1", job="PLQ")
    display_profiling_results(results_plq, "PLQ", "NaCl_1x1_L1")
    
    # Uncomment to run Ewald profiling
    # print("\nRunning profiling for Ewald job...")
    # results_ewald = run_profiled_test(name="NaCl_1x1_L1", job="Ewald")
    # display_profiling_results(results_ewald, "Ewald", "NaCl_1x1_L1")
    
    plt.show()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error during profiling: {e}")
        import traceback
        traceback.print_exc()
