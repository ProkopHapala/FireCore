#!/usr/bin/env python3
"""
Simple Python profiler with OpenCL support

This profiler can be used to analyze any Python script, with special handling for
OpenCL code. It uses cProfile for Python code profiling and a simple decorator-based
approach for tracking function execution times.

Usage:
    python simple_profiler.py [options] script_to_profile.py [script_args...]

Example:
    python simple_profiler.py --opencl ../run_test_GridFF_ocl_new.py
"""

import argparse
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

# Global variables to store profiling data
function_stats = defaultdict(lambda: {'calls': 0, 'total_time': 0})

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
        
        return result
    return wrapper

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
            try:
                setattr(module, name, profile_func(obj))
            except (AttributeError, TypeError):
                # Skip functions that can't be replaced
                pass
        elif inspect.isclass(obj) and obj.__module__ == module_name:
            # Process class methods
            for method_name, method in inspect.getmembers(obj):
                if inspect.isfunction(method) and not method_name.startswith('__'):
                    try:
                        # Try to patch the method
                        if hasattr(obj, method_name) and callable(getattr(obj, method_name)):
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
    if args.opencl:
        print("Enabling OpenCL profiling via environment variables...")
        # Set environment variable for OpenCL profiling
        os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
        os.environ['OPENCL_PROFILE'] = '1'
    
    # Patch modules for tracing if requested
    if args.trace:
        print("Enabling function tracing...")
        # We'll patch modules as they are imported
        import builtins
        original_import = builtins.__import__
        
        def patched_import(name, globals=None, locals=None, fromlist=(), level=0):
            module = original_import(name, globals, locals, fromlist, level)
            if args.trace_all or any(name.startswith(prefix) for prefix in args.trace_modules.split(',') if prefix):
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
            f.write("Function,Calls,Total Time (s)\n")
            for func_name, stats in sorted_funcs:
                f.write(f"{func_name},{stats['calls']},{stats['total_time']:.6f}\n")
        
        # Print top functions
        print(f"\n{'='*50}")
        print(f"Top {min(args.limit, len(sorted_funcs))} functions by total time:")
        print(f"{'Function':<60} {'Calls':<10} {'Total Time (s)':<15}")
        print("-" * 85)
        
        for func_name, stats in sorted_funcs[:args.limit]:
            print(f"{func_name:<60} {stats['calls']:<10} {stats['total_time']:<15.6f}")
    
    # Try to find OpenCL profiling information if available
    if args.opencl:
        try:
            import pyopencl as cl
            print(f"\n{'='*50}")
            print("OpenCL Information:")
            
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
    
    print(f"\nDetailed profiling results saved to {output_dir}/")
    print(f"{'='*50}")

def main():
    parser = argparse.ArgumentParser(description='Simple Python profiler with OpenCL support')
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
    
    # Run the script with profiling
    print(f"Profiling {args.script}...")
    profiler, total_time = run_with_profiling(args)
    
    # Save and display results
    save_and_display_results(profiler, total_time, args)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
