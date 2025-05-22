import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import pyopencl as cl
import argparse
import json
from collections import defaultdict

sys.path.append("../../")

from pyBall.tests import ocl_GridFF_new as gff
from pyBall.OCL.GridFF import GridFF_cl, GridShape

# Parse command line arguments
parser = argparse.ArgumentParser(description='Profile memory usage in GridFF OpenCL code')
parser.add_argument('--job', type=str, default='PLQ', help='Job type (PLQ, Morse, Ewald)')
parser.add_argument('--name', type=str, default='NaCl_1x1_L1', help='Dataset name')
parser.add_argument('--output', type=str, default='memory_profile', help='Output file prefix')
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

# Create a memory-tracking GridFF_cl class
class MemoryTrackingGridFF_cl(GridFF_cl):
    def __init__(self):
        super().__init__()
        self.ctx = ctx
        self.queue = queue
        self.memory_allocations = []
        self.buffer_sizes = {}
        self.kernel_events = {}
        self.kernel_times = {}
        self.total_time = 0
        self.memory_timeline = []
        self.current_memory = 0
        self.peak_memory = 0
        self.kernel_args = defaultdict(list)
    
    def create_buffer(self, size, flags=cl.mem_flags.READ_WRITE, hostbuf=None):
        # Track memory allocation
        buffer_id = len(self.memory_allocations)
        allocation_time = time.time()
        
        # Create the buffer
        if hostbuf is not None:
            buffer = cl.Buffer(self.ctx, flags, hostbuf=hostbuf)
            buffer_size = hostbuf.nbytes
        else:
            buffer = cl.Buffer(self.ctx, flags, size=size)
            buffer_size = size
        
        # Record allocation
        self.memory_allocations.append({
            'id': buffer_id,
            'size': buffer_size,
            'time': allocation_time,
            'type': 'allocation',
            'buffer': buffer
        })
        
        # Update current and peak memory usage
        self.current_memory += buffer_size
        self.peak_memory = max(self.peak_memory, self.current_memory)
        
        # Record in timeline
        self.memory_timeline.append({
            'time': allocation_time,
            'operation': 'allocate',
            'size': buffer_size,
            'current_memory': self.current_memory
        })
        
        # Store buffer size for later reference
        self.buffer_sizes[buffer] = buffer_size
        
        return buffer
    
    def release_buffer(self, buffer):
        if buffer in self.buffer_sizes:
            release_time = time.time()
            buffer_size = self.buffer_sizes[buffer]
            
            # Update current memory usage
            self.current_memory -= buffer_size
            
            # Record in timeline
            self.memory_timeline.append({
                'time': release_time,
                'operation': 'release',
                'size': buffer_size,
                'current_memory': self.current_memory
            })
            
            # Remove from tracking
            del self.buffer_sizes[buffer]
            
            # Release the buffer
            buffer.release()
    
    def enqueue_kernel(self, kernel, global_size, local_size=None, args=None):
        if args is None:
            args = []
        
        # Track kernel arguments
        kernel_name = kernel.function_name
        arg_sizes = []
        for arg in args:
            if isinstance(arg, cl.Buffer) and arg in self.buffer_sizes:
                arg_sizes.append(self.buffer_sizes[arg])
            elif isinstance(arg, np.ndarray):
                arg_sizes.append(arg.nbytes)
            else:
                arg_sizes.append(0)  # Scalar arguments
        
        self.kernel_args[kernel_name].append({
            'global_size': global_size,
            'local_size': local_size,
            'arg_sizes': arg_sizes,
            'time': time.time()
        })
        
        # Execute kernel with profiling
        event = kernel(queue, global_size, local_size, *args)
        
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
                'avg': sum(times) / len(times) if times else 0,
                'times': times
            }
            self.total_time += sum(times)
    
    def reset_tracking(self):
        # Release any remaining buffers
        for buffer in list(self.buffer_sizes.keys()):
            self.release_buffer(buffer)
        
        self.memory_allocations = []
        self.buffer_sizes = {}
        self.kernel_events = {}
        self.kernel_times = {}
        self.total_time = 0
        self.memory_timeline = []
        self.current_memory = 0
        self.peak_memory = 0
        self.kernel_args = defaultdict(list)

# Create memory tracking instance
mem_clgff = MemoryTrackingGridFF_cl()

# Function to run a memory profiling test
def run_memory_profiling(name=args.name, job=args.job):
    print(f"\n===== Memory Profiling GridFF OpenCL with {name} dataset, {job} job =====\n")
    
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
    
    # Reset tracking
    mem_clgff.reset_tracking()
    
    # Run the test
    start_time = time.time()
    
    if job == "PLQ":
        gff.test_gridFF_ocl(
            fname=f"data/xyz/{name}.xyz",
            Element_Types_name="./data/ElementTypes.dat",
            job=job,
            bPlot=False,  # Disable plotting for profiling
            clgff=mem_clgff  # Use our tracking instance
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
            clgff=mem_clgff
        )
    
    # Calculate elapsed time
    elapsed_time = time.time() - start_time
    
    # Calculate kernel times
    mem_clgff.calculate_kernel_times()
    
    # Prepare results
    results = {
        'device_info': device_info,
        'job': job,
        'dataset': name,
        'cpu_time': elapsed_time * 1000,  # ms
        'gpu_time': mem_clgff.total_time,
        'peak_memory': mem_clgff.peak_memory / (1024**2),  # MB
        'memory_timeline': mem_clgff.memory_timeline,
        'kernel_times': mem_clgff.kernel_times,
        'kernel_args': dict(mem_clgff.kernel_args)
    }
    
    # Print results
    print("\n===== Memory Profiling Results =====")
    print(f"Total CPU time: {results['cpu_time']:.2f} ms")
    print(f"Total GPU time: {results['gpu_time']:.2f} ms")
    print(f"Peak memory usage: {results['peak_memory']:.2f} MB")
    
    # Print top memory-consuming kernels
    print("\nTop kernels by memory usage:")
    kernel_memory = []
    for kernel_name, args_list in mem_clgff.kernel_args.items():
        total_arg_size = sum(sum(args['arg_sizes']) for args in args_list)
        kernel_memory.append((kernel_name, total_arg_size / (1024**2)))  # MB
    
    kernel_memory.sort(key=lambda x: x[1], reverse=True)
    for i, (kernel_name, memory) in enumerate(kernel_memory[:10]):
        print(f"{i+1}. {kernel_name}: {memory:.2f} MB")
    
    # Print top time-consuming kernels
    print("\nTop kernels by execution time:")
    kernel_times = [(k, v['total']) for k, v in mem_clgff.kernel_times.items()]
    kernel_times.sort(key=lambda x: x[1], reverse=True)
    for i, (kernel_name, time_ms) in enumerate(kernel_times[:10]):
        print(f"{i+1}. {kernel_name}: {time_ms:.2f} ms")
    
    # Save results to JSON
    output_file = f"{args.output}_{job}_{name}.json"
    with open(output_file, 'w') as f:
        # Convert numpy types to native Python types for JSON serialization
        def convert_to_serializable(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj
        
        # Create a serializable copy of results
        serializable_results = json.loads(json.dumps(results, default=convert_to_serializable))
        json.dump(serializable_results, f, indent=2)
    
    print(f"\nDetailed results saved to {output_file}")
    
    # Create visualizations
    plt.figure(figsize=(15, 10))
    
    # Memory usage timeline
    plt.subplot(2, 2, 1)
    timeline = mem_clgff.memory_timeline
    times = [(t['time'] - start_time) * 1000 for t in timeline]  # Convert to ms from start
    memory = [t['current_memory'] / (1024**2) for t in timeline]  # Convert to MB
    
    plt.plot(times, memory)
    plt.xlabel('Time (ms)')
    plt.ylabel('Memory Usage (MB)')
    plt.title('Memory Usage Timeline')
    plt.grid(True)
    
    # Kernel execution times
    plt.subplot(2, 2, 2)
    top_kernels = [k for k, _ in kernel_times[:10]]
    times = [mem_clgff.kernel_times[k]['total'] for k in top_kernels]
    
    plt.barh(range(len(top_kernels)), times)
    plt.yticks(range(len(top_kernels)), [k[:20] + '...' if len(k) > 20 else k for k in top_kernels])
    plt.xlabel('Execution Time (ms)')
    plt.title('Top 10 Kernels by Execution Time')
    plt.grid(True)
    
    # Memory usage by kernel
    plt.subplot(2, 2, 3)
    top_memory_kernels = [k for k, _ in kernel_memory[:10]]
    memory_sizes = [m for _, m in kernel_memory[:10]]
    
    plt.barh(range(len(top_memory_kernels)), memory_sizes)
    plt.yticks(range(len(top_memory_kernels)), [k[:20] + '...' if len(k) > 20 else k for k in top_memory_kernels])
    plt.xlabel('Memory Usage (MB)')
    plt.title('Top 10 Kernels by Memory Usage')
    plt.grid(True)
    
    # Kernel execution distribution (pie chart)
    plt.subplot(2, 2, 4)
    labels = [k for k, _ in kernel_times[:5]]  # Top 5 kernels
    sizes = [mem_clgff.kernel_times[k]['total'] for k in labels]
    if len(kernel_times) > 5:
        labels.append('Others')
        sizes.append(sum(t for _, t in kernel_times[5:]))
    
    plt.pie(sizes, labels=[k[:15] + '...' if len(k) > 15 else k for k in labels], 
            autopct='%1.1f%%', startangle=90)
    plt.axis('equal')
    plt.title('Kernel Execution Time Distribution')
    
    plt.tight_layout()
    plt.savefig(f"{args.output}_{job}_{name}.png")
    print(f"Visualization saved to {args.output}_{job}_{name}.png")
    
    return results

# Run the profiling
if __name__ == "__main__":
    try:
        # Run memory profiling
        results = run_memory_profiling()
        
        # Show plot if requested
        if '--show-plot' in sys.argv:
            plt.show()
    except Exception as e:
        print(f"Error during profiling: {e}")
        import traceback
        traceback.print_exc()
