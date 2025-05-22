#!/usr/bin/env python3
"""
General-purpose Python profiler with OpenCL support

This profiler can be used to analyze any Python script, with special handling for
OpenCL code. It uses various profiling techniques including cProfile for Python code
and PyOpenCL's event profiling for GPU kernels.

Usage:
    python profiler.py [options] script_to_profile.py [script_args...]

Example:
    python profiler.py --trace --opencl ../run_test_GridFF_ocl_new.py
"""

import argparse
import builtins
import cProfile
import io
import os
import pstats
import sys
import time
import importlib.util
import inspect
import traceback
import functools
import json
import datetime
from pathlib import Path
from collections import defaultdict

# Optional imports - will be checked at runtime
try:
    import pyopencl as cl
    HAVE_OPENCL = True
except ImportError:
    HAVE_OPENCL = False

try:
    import numpy as np
    HAVE_NUMPY = True
except ImportError:
    HAVE_NUMPY = False

try:
    import matplotlib.pyplot as plt
    HAVE_MATPLOTLIB = True
except ImportError:
    HAVE_MATPLOTLIB = False

# Global variables to store profiling data
function_stats = defaultdict(lambda: {'calls': 0, 'total_time': 0, 'own_time': 0})
kernel_stats = defaultdict(lambda: {'calls': 0, 'total_time': 0, 'avg_time': 0})
memory_stats = {'buffers': {}, 'total_allocated': 0, 'peak_allocated': 0, 'history': []}

# Timing decorator for function profiling
def profile_func(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        
        elapsed = end_time - start_time
        func_name = f"{func.__module__}.{func.__qualname__}"
        
        function_stats[func_name]['calls'] += 1
        function_stats[func_name]['total_time'] += elapsed
        function_stats[func_name]['own_time'] += elapsed
        
        return result
    return wrapper

# OpenCL profiling wrapper
def patch_opencl():
    if not HAVE_OPENCL:
        return False
    
    # Save original methods
    original_methods = {}
    
    # Track memory allocation
    def track_buffer_allocation(size_bytes, name="unknown"):
        memory_stats['buffers'][name] = size_bytes
        memory_stats['total_allocated'] += size_bytes
        if memory_stats['total_allocated'] > memory_stats['peak_allocated']:
            memory_stats['peak_allocated'] = memory_stats['total_allocated']
        
        memory_stats['history'].append({
            'time': time.time(),
            'operation': f'Allocate {name}',
            'size_bytes': size_bytes,
            'total_allocated': memory_stats['total_allocated']
        })
    
    # Patch Buffer creation
    original_methods['Buffer'] = cl.Buffer
    def profiled_buffer(context, flags, size=None, hostbuf=None):
        buffer = original_methods['Buffer'](context, flags, size, hostbuf)
        buffer_size = size if size is not None else (hostbuf.nbytes if hasattr(hostbuf, 'nbytes') else 0)
        
        # Get caller info for better naming
        frame = inspect.currentframe().f_back
        caller_info = inspect.getframeinfo(frame)
        buffer_name = f"Buffer_{os.path.basename(caller_info.filename)}:{caller_info.lineno}"
        
        track_buffer_allocation(buffer_size, buffer_name)
        return buffer
    
    cl.Buffer = profiled_buffer
    
    # Patch CommandQueue
    original_methods['CommandQueue'] = cl.CommandQueue
    def profiled_command_queue(context, device=None, properties=None):
        # Always enable profiling
        if properties is None:
            properties = cl.command_queue_properties.PROFILING_ENABLE
        elif not (properties & cl.command_queue_properties.PROFILING_ENABLE):
            properties |= cl.command_queue_properties.PROFILING_ENABLE
        
        return original_methods['CommandQueue'](context, device, properties)
    
    cl.CommandQueue = profiled_command_queue
    
    # In PyOpenCL, we need to patch the actual methods of CommandQueue instances
    # We'll do this by creating a subclass and monkey-patching the module
    
    # First, save the original CommandQueue class
    original_CommandQueue = cl.CommandQueue
    
    # Create a profiling subclass
    class ProfilingCommandQueue(cl.CommandQueue):
        def __init__(self, *args, **kwargs):
            # Ensure profiling is enabled
            if 'properties' not in kwargs:
                kwargs['properties'] = cl.command_queue_properties.PROFILING_ENABLE
            elif not (kwargs['properties'] & cl.command_queue_properties.PROFILING_ENABLE):
                kwargs['properties'] |= cl.command_queue_properties.PROFILING_ENABLE
            
            super().__init__(*args, **kwargs)
        
        def enqueue_nd_range_kernel(self, kernel, global_work_size, local_work_size=None, *args, **kwargs):
            # Get kernel name
            kernel_name = kernel.function_name
            
            # Execute kernel and get event
            start_time = time.time()
            event = super().enqueue_nd_range_kernel(kernel, global_work_size, local_work_size, *args, **kwargs)
            
            # Wait for kernel to complete
            event.wait()
            
            # Calculate GPU execution time
            gpu_start = event.profile.start
            gpu_end = event.profile.end
            gpu_time_ms = (gpu_end - gpu_start) / 1e6  # Convert ns to ms
            
            # Record kernel execution time
            kernel_stats[kernel_name]['calls'] += 1
            kernel_stats[kernel_name]['total_time'] += gpu_time_ms
            kernel_stats[kernel_name]['avg_time'] = (
                kernel_stats[kernel_name]['total_time'] / 
                kernel_stats[kernel_name]['calls']
            )
            
            # Record CPU overhead
            cpu_time_ms = (time.time() - start_time) * 1000
            kernel_stats[f"{kernel_name}_cpu_overhead"]['calls'] += 1
            kernel_stats[f"{kernel_name}_cpu_overhead"]['total_time'] += (cpu_time_ms - gpu_time_ms)
            kernel_stats[f"{kernel_name}_cpu_overhead"]['avg_time'] = (
                kernel_stats[f"{kernel_name}_cpu_overhead"]['total_time'] / 
                kernel_stats[f"{kernel_name}_cpu_overhead"]['calls']
            )
            
            return event
    
    # Replace the CommandQueue class with our profiling version
    cl.CommandQueue = ProfilingCommandQueue
    original_methods['CommandQueue'] = original_CommandQueue
    
    # Save original methods for restoration
    return original_methods

def restore_opencl(original_methods):
    if not HAVE_OPENCL or not original_methods:
        return
    
    for name, method in original_methods.items():
        if name == 'Buffer':
            cl.Buffer = method
        elif name == 'CommandQueue':
            cl.CommandQueue = method

# Function to patch modules for profiling
def patch_module(module_name, module=None):
    if module is None:
        if module_name in sys.modules:
            module = sys.modules[module_name]
        else:
            return False
    
    # Get all functions and methods in the module
    for name, obj in inspect.getmembers(module):
        if inspect.isfunction(obj) and obj.__module__ == module_name:
            # Replace the function with the profiled version
            setattr(module, name, profile_func(obj))
        elif inspect.isclass(obj) and obj.__module__ == module_name:
            # Process class methods
            for method_name, method in inspect.getmembers(obj):
                if inspect.isfunction(method) and not method_name.startswith('__'):
                    try:
                        setattr(obj, method_name, profile_func(method))
                    except (AttributeError, TypeError):
                        # Skip methods that can't be replaced
                        pass
    
    return True

# Function to run the target script with profiling
def run_with_profiling(args):
    # Save original sys.argv and set it to the target script's args
    original_argv = sys.argv.copy()
    sys.argv = [args.script] + args.script_args
    
    # Add script's directory to sys.path
    script_dir = os.path.dirname(os.path.abspath(args.script))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    
    # Set up profiling
    profiler = cProfile.Profile()
    
    # Patch OpenCL if requested
    original_opencl_methods = None
    if args.opencl and HAVE_OPENCL:
        print("Enabling OpenCL profiling...")
        original_opencl_methods = patch_opencl()
    
    # Patch modules for tracing if requested
    if args.trace:
        print("Enabling function tracing...")
        # We'll patch modules as they are imported
        original_import = __import__
        
        def patched_import(name, globals=None, locals=None, fromlist=(), level=0):
            module = original_import(name, globals, locals, fromlist, level)
            if args.trace_all or name.startswith(args.trace_modules):
                patch_module(name, module)
            return module
        
        builtins.__import__ = patched_import
    
    # Start profiling
    start_time = time.time()
    profiler.enable()
    
    try:
        # Load and run the target script
        spec = importlib.util.spec_from_file_location(
            os.path.basename(args.script).split('.')[0],
            args.script
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        
        # If script has a main function, call it
        if hasattr(module, 'main') and callable(module.main):
            module.main()
    except Exception as e:
        print(f"Error running script: {e}")
        traceback.print_exc()
    finally:
        # Stop profiling
        profiler.disable()
        end_time = time.time()
        
        # Restore original imports if needed
        if args.trace:
            import builtins
            builtins.__import__ = original_import
        
        # Restore OpenCL if patched
        if original_opencl_methods:
            restore_opencl(original_opencl_methods)
        
        # Restore original sys.argv
        sys.argv = original_argv
    
    return profiler, end_time - start_time

# Function to save and display profiling results
def save_and_display_results(profiler, total_time, args):
    # Create output directory if it doesn't exist
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    # Base filename for outputs
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_filename = os.path.join(
        output_dir, 
        f"{os.path.basename(args.script).split('.')[0]}_{timestamp}"
    )
    
    # Save cProfile results
    s = io.StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats('cumulative')
    ps.print_stats(args.limit)
    
    with open(f"{base_filename}_profile.txt", 'w') as f:
        f.write(s.getvalue())
    
    # Print summary to console
    print(f"\n{'='*50}")
    print(f"Profiling results for {args.script}")
    print(f"{'='*50}")
    print(f"Total execution time: {total_time:.2f} seconds")
    print(f"\nTop {args.limit} functions by cumulative time:")
    ps.print_stats(args.limit)
    
    # Save detailed function stats if we have any
    if function_stats and args.trace:
        # Sort by total time
        sorted_funcs = sorted(
            function_stats.items(),
            key=lambda x: x[1]['total_time'],
            reverse=True
        )
        
        # Save to CSV
        with open(f"{base_filename}_functions.csv", 'w') as f:
            f.write("Function,Calls,Total Time (s),Own Time (s)\n")
            for func_name, stats in sorted_funcs:
                f.write(f"{func_name},{stats['calls']},{stats['total_time']:.6f},{stats['own_time']:.6f}\n")
        
        # Print top functions
        print(f"\n{'='*50}")
        print(f"Top {min(args.limit, len(sorted_funcs))} functions by total time:")
        print(f"{'Function':<60} {'Calls':<10} {'Total Time (s)':<15} {'Own Time (s)':<15}")
        print("-" * 100)
        
        for func_name, stats in sorted_funcs[:args.limit]:
            print(f"{func_name:<60} {stats['calls']:<10} {stats['total_time']:<15.6f} {stats['own_time']:<15.6f}")
    
    # Save OpenCL kernel stats if we have any
    if kernel_stats and args.opencl:
        # Sort by total time
        sorted_kernels = sorted(
            [(k, v) for k, v in kernel_stats.items() if not k.endswith('_cpu_overhead')],
            key=lambda x: x[1]['total_time'],
            reverse=True
        )
        
        # Save to CSV
        with open(f"{base_filename}_kernels.csv", 'w') as f:
            f.write("Kernel,Calls,Total Time (ms),Avg Time (ms)\n")
            for kernel_name, stats in sorted_kernels:
                f.write(f"{kernel_name},{stats['calls']},{stats['total_time']:.3f},{stats['avg_time']:.3f}\n")
        
        # Print kernel stats
        if sorted_kernels:
            print(f"\n{'='*50}")
            print(f"OpenCL Kernel Execution Times:")
            print(f"{'Kernel':<40} {'Calls':<8} {'Total Time (ms)':<15} {'Avg Time (ms)':<15}")
            print("-" * 80)
            
            for kernel_name, stats in sorted_kernels:
                print(f"{kernel_name:<40} {stats['calls']:<8} {stats['total_time']:<15.3f} {stats['avg_time']:<15.3f}")
            
            # Calculate CPU overhead
            cpu_overhead = []
            for kernel_name, _ in sorted_kernels:
                overhead_key = f"{kernel_name}_cpu_overhead"
                if overhead_key in kernel_stats:
                    cpu_overhead.append((kernel_name, kernel_stats[overhead_key]))
            
            if cpu_overhead:
                print(f"\nCPU Overhead for Kernel Execution:")
                print(f"{'Kernel':<40} {'Overhead (ms)':<15} {'%':<10}")
                print("-" * 70)
                
                for kernel_name, stats in cpu_overhead:
                    kernel_time = kernel_stats[kernel_name]['total_time']
                    overhead = stats['total_time']
                    overhead_pct = (overhead / kernel_time * 100) if kernel_time > 0 else 0
                    print(f"{kernel_name:<40} {overhead:<15.3f} {overhead_pct:<10.2f}%")
    
    # Save memory stats if we have any
    if memory_stats['peak_allocated'] > 0 and args.opencl:
        # Save to JSON
        with open(f"{base_filename}_memory.json", 'w') as f:
            json.dump(memory_stats, f, indent=2)
        
        # Print memory stats
        print(f"\n{'='*50}")
        print(f"OpenCL Memory Usage:")
        print(f"Peak memory allocated: {memory_stats['peak_allocated'] / (1024*1024):.2f} MB")
        
        # Print top buffers
        if memory_stats['buffers']:
            sorted_buffers = sorted(
                memory_stats['buffers'].items(),
                key=lambda x: x[1],
                reverse=True
            )
            
            print(f"\nTop {min(10, len(sorted_buffers))} largest buffers:")
            print(f"{'Buffer':<50} {'Size (MB)':<15}")
            print("-" * 65)
            
            for name, size in sorted_buffers[:10]:
                print(f"{name:<50} {size / (1024*1024):<15.2f}")
    
    # Create visualizations if matplotlib is available
    if HAVE_MATPLOTLIB and (function_stats or kernel_stats):
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Plot function stats if available
        if function_stats and args.trace and len(sorted_funcs) > 0:
            plt.subplot(2, 1, 1)
            
            # Get top N functions
            top_n = min(10, len(sorted_funcs))
            func_names = [f"{name.split('.')[-1]}" for name, _ in sorted_funcs[:top_n]]
            func_times = [stats['total_time'] for _, stats in sorted_funcs[:top_n]]
            
            y_pos = range(len(func_names))
            plt.barh(y_pos, func_times)
            plt.yticks(y_pos, func_names)
            plt.xlabel('Time (s)')
            plt.title('Top Functions by Execution Time')
            
            # Add time labels
            for i, v in enumerate(func_times):
                plt.text(v + 0.01, i, f"{v:.2f}s", va='center')
        
        # Plot kernel stats if available
        if kernel_stats and args.opencl and sorted_kernels:
            plt.subplot(2, 1, 2)
            
            # Get top N kernels
            top_n = min(10, len(sorted_kernels))
            kernel_names = [name for name, _ in sorted_kernels[:top_n]]
            kernel_times = [stats['total_time'] for _, stats in sorted_kernels[:top_n]]
            
            y_pos = range(len(kernel_names))
            plt.barh(y_pos, kernel_times)
            plt.yticks(y_pos, kernel_names)
            plt.xlabel('Time (ms)')
            plt.title('OpenCL Kernel Execution Times')
            
            # Add time labels
            for i, v in enumerate(kernel_times):
                plt.text(v + 0.1, i, f"{v:.2f}ms", va='center')
        
        plt.tight_layout()
        plt.savefig(f"{base_filename}_chart.png")
        print(f"\nVisualization saved to {base_filename}_chart.png")
        
        # Plot memory usage over time if available
        if memory_stats['history'] and args.opencl:
            plt.figure(figsize=(12, 6))
            
            # Extract time and memory usage
            if memory_stats['history']:
                times = []
                memory_values = []
                start_time = memory_stats['history'][0]['time']
                
                for entry in memory_stats['history']:
                    times.append(entry['time'] - start_time)
                    memory_values.append(entry['total_allocated'] / (1024*1024))  # MB
                
                plt.plot(times, memory_values)
                plt.xlabel('Time (s)')
                plt.ylabel('Memory Usage (MB)')
                plt.title('OpenCL Memory Usage Over Time')
                plt.grid(True)
                
                # Add peak memory line
                peak_memory = memory_stats['peak_allocated'] / (1024*1024)
                plt.axhline(y=peak_memory, color='r', linestyle='--')
                plt.text(times[-1] * 0.5, peak_memory * 1.05, f'Peak: {peak_memory:.2f} MB', color='r')
                
                plt.tight_layout()
                plt.savefig(f"{base_filename}_memory.png")
                print(f"Memory usage visualization saved to {base_filename}_memory.png")
    
    print(f"\nDetailed profiling results saved to {output_dir}/")
    print(f"{'='*50}")

def main():
    parser = argparse.ArgumentParser(description='Python profiler with OpenCL support')
    parser.add_argument('script', help='Python script to profile')
    parser.add_argument('script_args', nargs='*', help='Arguments to pass to the script')
    parser.add_argument('--opencl', action='store_true', help='Enable OpenCL profiling')
    parser.add_argument('--trace', action='store_true', help='Enable function tracing')
    parser.add_argument('--trace-all', action='store_true', help='Trace all modules (can be slow)')
    parser.add_argument('--trace-modules', default='', help='Comma-separated list of module prefixes to trace')
    parser.add_argument('--output-dir', default='./profile_results', help='Directory to save results')
    parser.add_argument('--limit', type=int, default=20, help='Limit the number of functions in the output')
    
    args = parser.parse_args()
    
    # Check if the script exists
    if not os.path.isfile(args.script):
        print(f"Error: Script '{args.script}' not found")
        return 1
    
    # Check if OpenCL is available
    if args.opencl and not HAVE_OPENCL:
        print("Warning: OpenCL profiling requested but PyOpenCL is not installed")
    
    # Check if matplotlib is available
    if not HAVE_MATPLOTLIB:
        print("Warning: matplotlib is not installed, visualizations will be skipped")
    
    # Run the script with profiling
    print(f"Profiling {args.script}...")
    profiler, total_time = run_with_profiling(args)
    
    # Save and display results
    save_and_display_results(profiler, total_time, args)
    
    return 0

if __name__ == '__main__':
    import builtins
    sys.exit(main())
