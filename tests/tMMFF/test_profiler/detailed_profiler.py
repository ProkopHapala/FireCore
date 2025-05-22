#!/usr/bin/env python3
"""
Detailed Python profiler with OpenCL support

This profiler provides detailed information about Python function calls and OpenCL kernels
without modifying the target code.

Usage:
    python detailed_profiler.py your_script.py [your_script_args]
"""

import argparse
import atexit
import csv
import inspect
import json
import os
import sys
import time
import datetime
import traceback
from collections import defaultdict

# Global variables to store profiling data
function_stats = defaultdict(lambda: {'calls': 0, 'total_time': 0, 'own_time': 0, 'last_start': 0})
function_stack = []
call_graph = defaultdict(set)

# OpenCL profiling data
kernel_stats = defaultdict(lambda: {'calls': 0, 'total_time': 0, 'avg_time': 0})
memory_stats = {'buffers': {}, 'total_allocated': 0, 'peak_allocated': 0, 'history': []}

# Function to trace Python function calls
def trace_calls(frame, event, arg):
    if event != 'call' and event != 'return':
        return trace_calls
    
    # Get function name
    code = frame.f_code
    func_name = code.co_name
    if func_name == 'write':
        # Ignore write calls to avoid excessive logging
        return trace_calls
    
    # Get module name
    module_name = frame.f_globals.get('__name__', '')
    if module_name.startswith('_') or module_name.startswith('importlib'):
        # Skip internal modules
        return trace_calls
    
    # Get file name
    file_name = code.co_filename
    if '/site-packages/' in file_name or '/dist-packages/' in file_name:
        # Skip library code
        return trace_calls
    
    # Full function name
    full_name = f"{module_name}.{func_name}"
    
    if event == 'call':
        # Record call time
        function_stats[full_name]['last_start'] = time.time()
        function_stats[full_name]['calls'] += 1
        
        # Record call graph
        if function_stack:
            caller = function_stack[-1]
            call_graph[caller].add(full_name)
        
        # Push function onto stack
        function_stack.append(full_name)
    
    elif event == 'return':
        # Pop function from stack
        if function_stack and function_stack[-1] == full_name:
            function_stack.pop()
        
        # Calculate execution time
        if full_name in function_stats and 'last_start' in function_stats[full_name]:
            elapsed = time.time() - function_stats[full_name]['last_start']
            function_stats[full_name]['total_time'] += elapsed
            function_stats[full_name]['own_time'] += elapsed
            
            # Subtract time from parent functions (own time calculation)
            for parent in function_stack:
                if parent in function_stats:
                    function_stats[parent]['own_time'] -= elapsed
    
    return trace_calls

# Function to patch PyOpenCL for kernel profiling
def patch_pyopencl():
    try:
        import pyopencl as cl
        
        # Save original methods for kernel execution and buffer creation
        original_enqueue_kernel = cl.enqueue_nd_range_kernel
        original_buffer = cl.Buffer
        
        # Patch kernel execution
        def profiled_enqueue_kernel(queue, kernel, global_work_size, local_work_size=None, *args, **kwargs):
            # Get kernel name
            kernel_name = kernel.function_name if hasattr(kernel, 'function_name') else str(kernel)
            
            # Make sure profiling is enabled for this queue
            has_profiling = False
            try:
                if hasattr(queue, 'properties'):
                    has_profiling = bool(queue.properties & cl.command_queue_properties.PROFILING_ENABLE)
            except:
                pass  # Ignore errors when checking properties
            
            # Execute kernel and get event
            start_time = time.time()
            event = original_enqueue_kernel(queue, kernel, global_work_size, local_work_size, *args, **kwargs)
            
            # Wait for kernel to complete
            event.wait()
            
            # Calculate execution time
            cpu_time_ms = (time.time() - start_time) * 1000
            
            # Try to get GPU time if profiling is enabled
            gpu_time_ms = 0
            try:
                if has_profiling:
                    gpu_start = event.get_profiling_info(cl.profiling_info.START)
                    gpu_end = event.get_profiling_info(cl.profiling_info.END)
                    gpu_time_ms = (gpu_end - gpu_start) / 1e6  # Convert ns to ms
            except Exception as e:
                # Silently fail - profiling info might not be available
                pass
            
            # Record kernel execution time
            kernel_stats[kernel_name]['calls'] += 1
            kernel_stats[kernel_name]['total_time'] += gpu_time_ms if gpu_time_ms > 0 else cpu_time_ms
            kernel_stats[kernel_name]['avg_time'] = (
                kernel_stats[kernel_name]['total_time'] / 
                kernel_stats[kernel_name]['calls']
            )
            
            # Record CPU overhead if we have GPU time
            if gpu_time_ms > 0:
                overhead_key = f"{kernel_name}_cpu_overhead"
                kernel_stats[overhead_key]['calls'] += 1
                kernel_stats[overhead_key]['total_time'] += (cpu_time_ms - gpu_time_ms)
                kernel_stats[overhead_key]['avg_time'] = (
                    kernel_stats[overhead_key]['total_time'] / 
                    kernel_stats[overhead_key]['calls']
                )
            
            return event
        
        # Create a wrapper class for Buffer
        class ProfiledBuffer(original_buffer):
            def __new__(cls, context, flags, size=None, hostbuf=None):
                # Create the buffer using the original class
                buffer = original_buffer(context, flags, size, hostbuf)
                
                # Get the buffer size
                buffer_size = size if size is not None else (hostbuf.nbytes if hasattr(hostbuf, 'nbytes') else 0)
                
                # Get caller info for better naming
                try:
                    frame = inspect.currentframe().f_back
                    caller_info = inspect.getframeinfo(frame)
                    buffer_name = f"Buffer_{os.path.basename(caller_info.filename)}:{caller_info.lineno}"
                except:
                    buffer_name = f"Buffer_{time.time()}"
                
                # Track memory allocation
                memory_stats['buffers'][buffer_name] = buffer_size
                memory_stats['total_allocated'] += buffer_size
                if memory_stats['total_allocated'] > memory_stats['peak_allocated']:
                    memory_stats['peak_allocated'] = memory_stats['total_allocated']
                
                memory_stats['history'].append({
                    'time': time.time(),
                    'operation': f'Allocate {buffer_name}',
                    'size_bytes': buffer_size,
                    'total_allocated': memory_stats['total_allocated']
                })
                
                return buffer
        
        # Apply patches
        cl.enqueue_nd_range_kernel = profiled_enqueue_kernel
        cl.Buffer = profiled_buffer
        
        print("PyOpenCL patched for profiling")
        return True
    except ImportError:
        print("PyOpenCL not available, skipping OpenCL profiling")
        return False
    except Exception as e:
        print(f"Error patching PyOpenCL: {e}")
        traceback.print_exc()
        return False

# Function to show OpenCL device information
def show_opencl_info():
    try:
        import pyopencl as cl
        print("\nOpenCL Information:")
        print("-" * 50)
        
        # Get platform and device information
        platforms = cl.get_platforms()
        for i, platform in enumerate(platforms):
            print(f"Platform {i}: {platform.name} ({platform.version})")
            
            devices = platform.get_devices()
            for j, device in enumerate(devices):
                print(f"  Device {j}: {device.name} (Type: {cl.device_type.to_string(device.type)})")
                print(f"    Compute Units: {device.max_compute_units}")
                print(f"    Global Memory: {device.global_mem_size / (1024**2):.2f} MB")
                print(f"    Local Memory: {device.local_mem_size / 1024:.2f} KB")
                print(f"    Max Work Group Size: {device.max_work_group_size}")
                if hasattr(device, 'max_clock_frequency'):
                    print(f"    Max Clock Frequency: {device.max_clock_frequency} MHz")
    except ImportError:
        print("PyOpenCL not available, skipping OpenCL information")
    except Exception as e:
        print(f"Error getting OpenCL information: {e}")

# Function to save profiling results
def save_results(output_dir, limit=20, sort_key='time'):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Base filename for outputs
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_filename = os.path.join(output_dir, f"profile_{timestamp}")
    
    # Save function stats
    if function_stats:
        # Determine sort key
        if sort_key == 'time' or sort_key == 'cumulative':
            sort_lambda = lambda x: x[1]['total_time']
        elif sort_key == 'calls':
            sort_lambda = lambda x: x[1]['calls']
        elif sort_key == 'name':
            sort_lambda = lambda x: x[0]
        else:
            sort_lambda = lambda x: x[1]['total_time']
        
        # Sort by the selected key
        sorted_funcs = sorted(
            function_stats.items(),
            key=sort_lambda,
            reverse=(sort_key != 'name')  # Don't reverse if sorting by name
        )
        
        # Save to CSV
        with open(f"{base_filename}_functions.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Function', 'Calls', 'Total Time (s)', 'Own Time (s)', 'Avg Time (s)'])
            for func_name, stats in sorted_funcs:
                avg_time = stats['total_time'] / stats['calls'] if stats['calls'] > 0 else 0
                writer.writerow([
                    func_name, 
                    stats['calls'], 
                    f"{stats['total_time']:.6f}", 
                    f"{stats['own_time']:.6f}",
                    f"{avg_time:.6f}"
                ])
        
        # Save call graph as JSON
        graph_data = {caller: list(callees) for caller, callees in call_graph.items()}
        with open(f"{base_filename}_call_graph.json", 'w') as f:
            json.dump(graph_data, f, indent=2)
            
        # Save a summary text file with all profiling information
        with open(f"{base_filename}_summary.txt", 'w') as f:
            f.write(f"Profiling Summary\n")
            f.write(f"================\n\n")
            
            # Write function stats
            f.write(f"Top {limit} functions by {sort_key}:\n")
            f.write(f"{'Function':<60} {'Calls':<10} {'Total Time (s)':<15} {'Own Time (s)':<15} {'Avg Time (ms)':<15}\n")
            f.write("-" * 115 + "\n")
            
            for func_name, stats in sorted_funcs[:limit]:
                avg_time = (stats['total_time'] / stats['calls'] * 1000) if stats['calls'] > 0 else 0
                f.write(f"{func_name:<60} {stats['calls']:<10} {stats['total_time']:<15.6f} {stats['own_time']:<15.6f} {avg_time:<15.2f}\n")
        
        # Print top functions
        print(f"\n{'='*50}")
        print(f"Top {limit} functions by {sort_key}:")
        print(f"{'Function':<60} {'Calls':<10} {'Total Time (s)':<15} {'Own Time (s)':<15} {'Avg Time (ms)':<15}")
        print("-" * 115)
        
        for func_name, stats in sorted_funcs[:limit]:
            avg_time = (stats['total_time'] / stats['calls'] * 1000) if stats['calls'] > 0 else 0
            print(f"{func_name:<60} {stats['calls']:<10} {stats['total_time']:<15.6f} {stats['own_time']:<15.6f} {avg_time:<15.2f}")
    
    # Save kernel stats
    if kernel_stats:
        # Sort by total time
        sorted_kernels = sorted(
            [(k, v) for k, v in kernel_stats.items() if not k.endswith('_cpu_overhead')],
            key=lambda x: x[1]['total_time'],
            reverse=True
        )
        
        # Save to CSV
        with open(f"{base_filename}_kernels.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Kernel', 'Calls', 'Total Time (ms)', 'Avg Time (ms)'])
            for kernel_name, stats in sorted_kernels:
                writer.writerow([
                    kernel_name, 
                    stats['calls'], 
                    f"{stats['total_time']:.3f}", 
                    f"{stats['avg_time']:.3f}"
                ])
        
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
    
    # Save memory stats
    if memory_stats['peak_allocated'] > 0:
        # Save to JSON
        with open(f"{base_filename}_memory.json", 'w') as f:
            json.dump(memory_stats, f, indent=2)
        
        # Add memory stats to summary file
        with open(f"{base_filename}_summary.txt", 'a') as f:
            f.write(f"\n\nOpenCL Memory Usage:\n")
            f.write(f"Peak memory allocated: {memory_stats['peak_allocated'] / (1024*1024):.2f} MB\n")
        
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
            
            print(f"\nTop 10 largest buffers:")
            print(f"{'Buffer':<50} {'Size (MB)':<15}")
            print("-" * 65)
            
            for name, size in sorted_buffers[:10]:
                print(f"{name:<50} {size / (1024*1024):<15.2f}")
    
    print(f"\nDetailed profiling results saved to {output_dir}/")
    return base_filename

# Function to create visualizations
def create_visualizations(base_filename, limit=15):
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        
        # Create figure for function stats
        if function_stats:
            plt.figure(figsize=(12, 8))
            
            # Sort by total time
            sorted_funcs = sorted(
                function_stats.items(),
                key=lambda x: x[1]['total_time'],
                reverse=True
            )
            
            # Get top functions
            top_n = min(limit, len(sorted_funcs))
            func_names = [name.split('.')[-1] for name, _ in sorted_funcs[:top_n]]
            func_times = [stats['total_time'] for _, stats in sorted_funcs[:top_n]]
            
            # Create horizontal bar chart
            y_pos = np.arange(len(func_names))
            plt.barh(y_pos, func_times)
            plt.yticks(y_pos, func_names)
            plt.xlabel('Time (s)')
            plt.title('Top Functions by Execution Time')
            
            # Add time labels
            for i, v in enumerate(func_times):
                plt.text(v + 0.01, i, f"{v:.2f}s", va='center')
            
            plt.tight_layout()
            plt.savefig(f"{base_filename}_functions.png")
            print(f"Function visualization saved to {base_filename}_functions.png")
        
        # Create figure for kernel stats
        if kernel_stats:
            # Sort by total time
            sorted_kernels = sorted(
                [(k, v) for k, v in kernel_stats.items() if not k.endswith('_cpu_overhead')],
                key=lambda x: x[1]['total_time'],
                reverse=True
            )
            
            # Add kernel stats to summary file
            with open(f"{base_filename}_summary.txt", 'a') as f:
                if sorted_kernels:
                    f.write(f"\n\nOpenCL Kernel Execution Times:\n")
                    f.write(f"{'Kernel':<40} {'Calls':<8} {'Total Time (ms)':<15} {'Avg Time (ms)':<15}\n")
                    f.write("-" * 80 + "\n")
                    
                    for kernel_name, stats in sorted_kernels:
                        f.write(f"{kernel_name:<40} {stats['calls']:<8} {stats['total_time']:<15.3f} {stats['avg_time']:<15.3f}\n")
            
            if sorted_kernels:
                plt.figure(figsize=(12, 8))
                
                # Get top kernels
                top_n = min(limit, len(sorted_kernels))
                kernel_names = [name for name, _ in sorted_kernels[:top_n]]
                kernel_times = [stats['total_time'] for _, stats in sorted_kernels[:top_n]]
                
                # Create horizontal bar chart
                y_pos = np.arange(len(kernel_names))
                plt.barh(y_pos, kernel_times)
                plt.yticks(y_pos, kernel_names)
                plt.xlabel('Time (ms)')
                plt.title('OpenCL Kernel Execution Times')
                
                # Add time labels
                for i, v in enumerate(kernel_times):
                    plt.text(v + 0.1, i, f"{v:.2f}ms", va='center')
                
                plt.tight_layout()
                plt.savefig(f"{base_filename}_kernels.png")
                print(f"Kernel visualization saved to {base_filename}_kernels.png")
        
        # Create memory usage visualization
        if memory_stats['history']:
            plt.figure(figsize=(12, 6))
            
            # Extract time and memory usage
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
    
    except ImportError:
        print("Matplotlib not available, skipping visualizations")
    except Exception as e:
        print(f"Error creating visualizations: {e}")
        traceback.print_exc()

# Main function to run the profiler
def main():
    parser = argparse.ArgumentParser(description='Detailed Python profiler with OpenCL support')
    parser.add_argument('script', help='Python script to profile')
    parser.add_argument('args', nargs='*', help='Arguments to pass to the script')
    parser.add_argument('--output', default='profile_results', help='Output directory for profiling results')
    parser.add_argument('--no-opencl', action='store_true', help='Disable OpenCL profiling')
    parser.add_argument('--no-trace', action='store_true', help='Disable function tracing')
    parser.add_argument('--limit', type=int, default=20, help='Limit the number of functions to show in output')
    parser.add_argument('--sort', choices=['time', 'cumulative', 'calls', 'name'], default='time', 
                        help='Sort order for function statistics')
    
    args = parser.parse_args()
    
    # Set up OpenCL environment variables for profiling
    if not args.no_opencl:
        os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
        os.environ['OPENCL_PROFILE'] = '1'
        # For NVIDIA
        os.environ['CUDA_PROFILE'] = '1'
        os.environ['CUDA_PROFILE_CSV'] = '1'
        # For AMD
        os.environ['AMD_OCL_BUILD_OPTIONS_APPEND'] = '-cl-opt-disable'
        # For Intel
        os.environ['INTEL_OPENCL_PROFILE'] = '1'
    
    # Check if the script exists
    if not os.path.isfile(args.script):
        print(f"Error: Script '{args.script}' not found")
        return 1
    
    # Show OpenCL information
    show_opencl_info()
    
    # Patch PyOpenCL if requested
    if not args.no_opencl:
        patch_pyopencl()
    
    # Set up function tracing if requested
    if not args.no_trace:
        sys.settrace(trace_calls)
    
    # Register exit handler to save results
    atexit.register(lambda: save_results(args.output, args.limit, args.sort))
    
    # Save original sys.argv and set it to the target script's args
    original_argv = sys.argv.copy()
    sys.argv = [args.script] + args.args
    
    # Add script's directory and project root to sys.path
    script_dir = os.path.dirname(os.path.abspath(args.script))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    
    # Add project root (assuming FireCore structure)
    project_root = os.path.abspath(os.path.join(script_dir, '..', '..', '..'))
    if os.path.exists(os.path.join(project_root, 'pyBall')):
        print(f"Adding project root to Python path: {project_root}")
        if project_root not in sys.path:
            sys.path.insert(0, project_root)
    
    # Record start time
    start_time = time.time()
    
    try:
        # Load and run the target script
        print(f"Running {args.script} with profiling...")
        
        # Execute the script
        with open(args.script, 'rb') as f:
            code = compile(f.read(), args.script, 'exec')
            exec(code, {'__name__': '__main__'})
    
    except Exception as e:
        print(f"Error running script: {e}")
        traceback.print_exc()
    
    finally:
        # Calculate total time
        end_time = time.time()
        total_time = end_time - start_time
        
        print(f"\n{'='*50}")
        print(f"Profiling complete for {args.script}")
        print(f"Total execution time: {total_time:.2f} seconds")
        
        # Disable tracing
        sys.settrace(None)
        
        # Restore original sys.argv
        sys.argv = original_argv
        
        # Save results and create visualizations
        base_filename = save_results(args.output, args.limit, args.sort)
        create_visualizations(base_filename)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
