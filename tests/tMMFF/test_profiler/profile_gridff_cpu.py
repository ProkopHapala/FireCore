import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import pyopencl as cl
import csv

sys.path.append("../../")

from pyBall.tests import ocl_GridFF_new as gff
from pyBall.OCL.GridFF import GridFF_cl, GridShape

# Enable OpenCL profiling
os.environ['PYOPENCL_CTX'] = '0'  # Use the first available device
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'

# Create a context with profiling enabled
platforms = cl.get_platforms()
devices = platforms[0].get_devices()
ctx = cl.Context(devices=devices)
queue = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)

# Create a profiling GridFF_cl class
class ProfilingGridFF_cl(GridFF_cl):
    def __init__(self):
        super().__init__()
        self.ctx = ctx
        self.queue = queue
        self.kernel_events = {}
        self.kernel_times = {}
        self.total_time = 0
        self.kernel_counts = {}
    
    def enqueue_kernel(self, kernel, global_size, local_size=None, args=None):
        if args is None:
            args = []
        
        # Track kernel invocation
        kernel_name = kernel.function_name
        if kernel_name not in self.kernel_counts:
            self.kernel_counts[kernel_name] = 0
        self.kernel_counts[kernel_name] += 1
        
        # Prepare arguments
        kernel_args = []
        for arg in args:
            kernel_args.append(arg)
        
        # Execute kernel with profiling
        event = kernel(queue, global_size, local_size, *kernel_args)
        
        # Store event for later analysis
        if kernel_name not in self.kernel_events:
            self.kernel_events[kernel_name] = []
        self.kernel_events[kernel_name].append(event)
        
        return event
    
    def calculate_kernel_times(self):
        self.kernel_times = {}
        self.total_time = 0
        
        for kernel_name, events in self.kernel_events.items():
            times = []
            for event in events:
                start = event.get_profiling_info(cl.profiling_info.START)
                end = event.get_profiling_info(cl.profiling_info.END)
                time_ns = end - start
                times.append(time_ns * 1e-6)  # Convert to ms
            
            self.kernel_times[kernel_name] = {
                'count': len(times),
                'total': sum(times),
                'min': min(times) if times else 0,
                'max': max(times) if times else 0,
                'avg': sum(times) / len(times) if times else 0
            }
            self.total_time += sum(times)
    
    def reset_profiling(self):
        self.kernel_events = {}
        self.kernel_times = {}
        self.total_time = 0
        self.kernel_counts = {}

# Create profiling instance
prof_clgff = ProfilingGridFF_cl()

# Function to profile CPU vs OpenCL performance
def profile_cpu_vs_opencl(name="NaCl_1x1_L1", job="PLQ", iterations=3):
    print(f"\n===== Profiling GridFF CPU vs OpenCL with {name} dataset, {job} job =====\n")
    
    # Get device info
    device = devices[0]
    device_info = {
        'name': device.name,
        'type': cl.device_type.to_string(device.type),
        'vendor': device.vendor,
        'version': device.version,
        'driver_version': device.driver_version,
        'compute_units': device.max_compute_units,
        'global_mem': device.global_mem_size / (1024**2),  # MB
        'local_mem': device.local_mem_size / 1024,  # KB
        'max_work_group_size': device.max_work_group_size
    }
    
    print("Device Information:")
    for key, value in device_info.items():
        print(f"  {key}: {value}")
    
    # Results containers
    results = {
        'opencl': {'total': 0, 'kernels': {}},
        'cpu': {'total': 0, 'functions': {}}
    }
    
    # Run multiple iterations for averaging
    for i in range(iterations):
        print(f"\nRunning iteration {i+1}/{iterations}...")
        
        # Reset profiling data
        prof_clgff.reset_profiling()
        
        # Profile OpenCL execution
        print("\nProfiling OpenCL execution...")
        start_time = time.time()
        
        if job == "PLQ":
            gff.test_gridFF_ocl(
                fname=f"data/xyz/{name}.xyz",
                Element_Types_name="./data/ElementTypes.dat",
                job=job,
                bPlot=False,  # Disable plotting for profiling
                clgff=prof_clgff  # Use our profiling instance
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
                bPlot=False,
                clgff=prof_clgff
            )
        
        # Calculate elapsed time
        opencl_time = time.time() - start_time
        
        # Calculate kernel times
        prof_clgff.calculate_kernel_times()
        
        # Store OpenCL results for this iteration
        results['opencl']['total'] += opencl_time
        for kernel_name, stats in prof_clgff.kernel_times.items():
            if kernel_name not in results['opencl']['kernels']:
                results['opencl']['kernels'][kernel_name] = {
                    'count': 0,
                    'total': 0,
                    'gpu_time': 0
                }
            results['opencl']['kernels'][kernel_name]['count'] += stats['count']
            results['opencl']['kernels'][kernel_name]['total'] += stats['total']
            results['opencl']['kernels'][kernel_name]['gpu_time'] += stats['total']
    
    # Calculate averages
    results['opencl']['total'] /= iterations
    for kernel_name in results['opencl']['kernels']:
        results['opencl']['kernels'][kernel_name]['count'] //= iterations
        results['opencl']['kernels'][kernel_name]['total'] /= iterations
        results['opencl']['kernels'][kernel_name]['gpu_time'] /= iterations
    
    # Print results
    print("\n===== Profiling Results =====")
    print(f"Total OpenCL execution time: {results['opencl']['total']*1000:.2f} ms (averaged over {iterations} iterations)")
    print(f"Total GPU kernel time: {prof_clgff.total_time:.2f} ms")
    print(f"CPU overhead: {(results['opencl']['total']*1000 - prof_clgff.total_time):.2f} ms")
    
    print("\nTop 10 OpenCL kernels by execution time:")
    sorted_kernels = sorted(prof_clgff.kernel_times.items(), key=lambda x: x[1]['total'], reverse=True)[:10]
    print(f"{'Kernel':<50} {'Count':<8} {'Total (ms)':<12} {'Avg (ms)':<10} {'% of Total':<10}")
    print("-" * 90)
    for kernel_name, stats in sorted_kernels:
        percent = (stats['total'] / prof_clgff.total_time) * 100
        print(f"{kernel_name:<50} {stats['count']:<8} {stats['total']:<12.2f} {stats['avg']:<10.2f} {percent:<10.2f}%")
    
    # Save results to CSV
    output_file = f"profile_results_{job}_{name}.csv"
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Profiling Results for", f"{job} job on {name} dataset"])
        writer.writerow(["Device", device_info['name']])
        writer.writerow(["Total OpenCL execution time (ms)", f"{results['opencl']['total']*1000:.2f}"])
        writer.writerow(["Total GPU kernel time (ms)", f"{prof_clgff.total_time:.2f}"])
        writer.writerow(["CPU overhead (ms)", f"{(results['opencl']['total']*1000 - prof_clgff.total_time):.2f}"])
        writer.writerow([])
        writer.writerow(["Kernel", "Count", "Total (ms)", "Avg (ms)", "% of Total"])
        for kernel_name, stats in sorted_kernels:
            percent = (stats['total'] / prof_clgff.total_time) * 100
            writer.writerow([kernel_name, stats['count'], f"{stats['total']:.2f}", f"{stats['avg']:.2f}", f"{percent:.2f}%"])
    
    print(f"\nResults saved to {output_file}")
    
    # Create visualization
    plt.figure(figsize=(15, 10))
    
    # Pie chart of kernel times
    plt.subplot(2, 1, 1)
    labels = [k for k, _ in sorted_kernels[:5]]  # Top 5 kernels
    sizes = [s['total'] for _, s in sorted_kernels[:5]]
    if len(sorted_kernels) > 5:
        labels.append('Others')
        sizes.append(sum(s['total'] for _, s in sorted_kernels[5:]))
    
    plt.pie(sizes, labels=[k.split('.')[-1] if '.' in k else k for k in labels], autopct='%1.1f%%', startangle=90)
    plt.axis('equal')
    plt.title('Kernel Execution Time Distribution')
    
    # Bar chart of kernel times
    plt.subplot(2, 1, 2)
    kernels = [k.split('.')[-1] if '.' in k else k for k, _ in sorted_kernels[:10]]  # Top 10 kernels
    times = [s['total'] for _, s in sorted_kernels[:10]]
    
    plt.barh(range(len(kernels)), times, align='center')
    plt.yticks(range(len(kernels)), kernels)
    plt.xlabel('Execution Time (ms)')
    plt.title('Top 10 Kernel Execution Times')
    
    plt.tight_layout()
    plt.savefig(f"profile_results_{job}_{name}.png")
    print(f"Visualization saved to profile_results_{job}_{name}.png")
    
    return results

# Run the profiling
if __name__ == "__main__":
    try:
        # Run profiling for different jobs
        profile_cpu_vs_opencl(name="NaCl_1x1_L1", job="PLQ")
        
        # Uncomment to run additional profiling jobs
        # profile_cpu_vs_opencl(name="NaCl_1x1_L1", job="Ewald")
        # profile_cpu_vs_opencl(name="NaCl_1x1_L1", job="Morse")
        
        # Show plots
        plt.show()
    except Exception as e:
        print(f"Error during profiling: {e}")
        import traceback
        traceback.print_exc()
