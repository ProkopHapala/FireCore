import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import pyopencl as cl
import csv
import argparse
from tabulate import tabulate

sys.path.append("../../")

from pyBall.tests import ocl_GridFF_new as gff
from pyBall.OCL.GridFF import GridFF_cl, GridShape

# Parse command line arguments
parser = argparse.ArgumentParser(description='Profile GridFF OpenCL code')
parser.add_argument('--job', type=str, default='PLQ', help='Job type (PLQ, Morse, Ewald)')
parser.add_argument('--name', type=str, default='NaCl_1x1_L1', help='Dataset name')
parser.add_argument('--output', type=str, default='profile_results', help='Output file prefix')
parser.add_argument('--iterations', type=int, default=5, help='Number of iterations for averaging')
parser.add_argument('--device', type=int, default=0, help='OpenCL device index')
args = parser.parse_args()

# Enable OpenCL profiling
os.environ['PYOPENCL_CTX'] = str(args.device)
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'

# Get platform and device info
platforms = cl.get_platforms()
print(f"\nAvailable platforms:")
for i, platform in enumerate(platforms):
    print(f"  [{i}] {platform.name} ({platform.version})")
    devices = platform.get_devices()
    for j, device in enumerate(devices):
        print(f"      [{j}] {device.name} (Type: {cl.device_type.to_string(device.type)})")

# Create a context with profiling enabled
devices = platforms[0].get_devices(device_type=cl.device_type.GPU)
ctx = cl.Context(devices=devices)
queue = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)

# Create a new GridFF_cl instance with profiling enabled
class ProfilingGridFF_cl(GridFF_cl):
    def __init__(self):
        super().__init__()
        self.ctx = ctx
        self.queue = queue
        self.kernel_events = {}
        self.kernel_times = {}
        self.total_time = 0
    
    def enqueue_kernel(self, kernel, global_size, local_size=None, args=None):
        if args is None:
            args = []
        
        # Prepare arguments
        kernel_args = []
        for arg in args:
            kernel_args.append(arg)
        
        # Execute kernel with profiling
        event = kernel(queue, global_size, local_size, *kernel_args)
        
        # Store event for later analysis
        kernel_name = kernel.function_name
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
                'min': min(times),
                'max': max(times),
                'avg': sum(times) / len(times)
            }
            self.total_time += sum(times)
    
    def reset_profiling(self):
        self.kernel_events = {}
        self.kernel_times = {}
        self.total_time = 0

# Create profiling instance
prof_clgff = ProfilingGridFF_cl()

# Function to run a profiling test
def run_profiling_test(name=args.name, job=args.job, iterations=args.iterations):
    results = []
    
    print(f"\n===== Profiling GridFF OpenCL with {name} dataset, {job} job =====\n")
    
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
    
    # Run multiple iterations for averaging
    iteration_times = []
    for i in range(iterations):
        print(f"\nRunning iteration {i+1}/{iterations}...")
        
        # Reset profiling data
        prof_clgff.reset_profiling()
        
        # Measure total CPU time
        start_time = time.time()
        
        # Run the test
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
        elapsed_time = time.time() - start_time
        
        # Calculate kernel times
        prof_clgff.calculate_kernel_times()
        
        # Store results for this iteration
        iteration_times.append({
            'cpu_time': elapsed_time * 1000,  # ms
            'gpu_time': prof_clgff.total_time,
            'kernel_times': prof_clgff.kernel_times.copy()
        })
    
    # Calculate average times
    avg_cpu_time = sum(it['cpu_time'] for it in iteration_times) / iterations
    avg_gpu_time = sum(it['gpu_time'] for it in iteration_times) / iterations
    
    # Combine kernel statistics across iterations
    kernel_stats = {}
    for it in iteration_times:
        for kernel_name, stats in it['kernel_times'].items():
            if kernel_name not in kernel_stats:
                kernel_stats[kernel_name] = {
                    'count': 0,
                    'total': 0,
                    'times': []
                }
            kernel_stats[kernel_name]['count'] += stats['count']
            kernel_stats[kernel_name]['total'] += stats['total']
            kernel_stats[kernel_name]['times'].extend([stats['avg']] * stats['count'])
    
    # Calculate final statistics
    for kernel_name, stats in kernel_stats.items():
        times = stats['times']
        stats['min'] = min(times)
        stats['max'] = max(times)
        stats['avg'] = stats['total'] / stats['count']
        stats['percent'] = (stats['total'] / (avg_gpu_time * iterations)) * 100
    
    # Sort kernels by total time
    sorted_kernels = sorted(kernel_stats.items(), key=lambda x: x[1]['total'], reverse=True)
    
    # Print results
    print("\n===== Profiling Results =====")
    print(f"Total CPU time: {avg_cpu_time:.2f} ms (averaged over {iterations} iterations)")
    print(f"Total GPU time: {avg_gpu_time:.2f} ms (sum of all kernel executions)")
    print(f"CPU/GPU ratio: {avg_cpu_time / avg_gpu_time:.2f}")
    
    print("\nKernel Execution Times:")
    table_data = []
    for kernel_name, stats in sorted_kernels:
        table_data.append([
            kernel_name,
            stats['count'] // iterations,
            f"{stats['total'] / iterations:.2f}",
            f"{stats['min']:.2f}",
            f"{stats['avg']:.2f}",
            f"{stats['max']:.2f}",
            f"{stats['percent']:.2f}%"
        ])
    
    headers = ["Kernel", "Calls", "Total (ms)", "Min (ms)", "Avg (ms)", "Max (ms)", "% of Total"]
    print(tabulate(table_data, headers=headers, tablefmt="grid"))
    
    # Save results to CSV
    output_file = f"{args.output}_{job}_{name}.csv"
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Profiling Results for", f"{job} job on {name} dataset"])
        writer.writerow(["Device", device_info['name']])
        writer.writerow(["Total CPU time (ms)", f"{avg_cpu_time:.2f}"])
        writer.writerow(["Total GPU time (ms)", f"{avg_gpu_time:.2f}"])
        writer.writerow(["CPU/GPU ratio", f"{avg_cpu_time / avg_gpu_time:.2f}"])
        writer.writerow([])
        writer.writerow(headers)
        for row in table_data:
            writer.writerow(row)
    
    print(f"\nResults saved to {output_file}")
    
    # Create visualization
    plt.figure(figsize=(12, 6))
    
    # Pie chart of kernel times
    plt.subplot(1, 2, 1)
    labels = [k for k, _ in sorted_kernels[:5]]  # Top 5 kernels
    sizes = [s['percent'] for _, s in sorted_kernels[:5]]
    if len(sorted_kernels) > 5:
        labels.append('Others')
        sizes.append(sum(s['percent'] for _, s in sorted_kernels[5:]))
    
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
    plt.axis('equal')
    plt.title('Kernel Execution Time Distribution')
    
    # Bar chart of kernel times
    plt.subplot(1, 2, 2)
    kernels = [k[:20] + '...' if len(k) > 20 else k for k, _ in sorted_kernels[:10]]  # Top 10 kernels
    times = [s['total'] / iterations for _, s in sorted_kernels[:10]]
    
    plt.barh(range(len(kernels)), times, align='center')
    plt.yticks(range(len(kernels)), kernels)
    plt.xlabel('Execution Time (ms)')
    plt.title('Top 10 Kernel Execution Times')
    
    plt.tight_layout()
    plt.savefig(f"{args.output}_{job}_{name}.png")
    print(f"Visualization saved to {args.output}_{job}_{name}.png")
    
    return {
        'avg_cpu_time': avg_cpu_time,
        'avg_gpu_time': avg_gpu_time,
        'kernel_stats': kernel_stats
    }

# Run the profiling
if __name__ == "__main__":
    try:
        # Install tabulate if not present
        import importlib.util
        if importlib.util.find_spec("tabulate") is None:
            print("Installing tabulate package...")
            import subprocess
            subprocess.check_call([sys.executable, "-m", "pip", "install", "tabulate"])
            from tabulate import tabulate
        
        # Run profiling
        results = run_profiling_test()
        
        # Show plot if requested
        if '--show-plot' in sys.argv:
            plt.show()
    except Exception as e:
        print(f"Error during profiling: {e}")
        import traceback
        traceback.print_exc()
