#!/usr/bin/env python3
"""
Plug-and-play Python profiler

This profiler runs any Python command as-is and collects profiling information.
It doesn't modify the target script or interfere with its arguments or execution flow.

Usage:
    python run_profiler.py [profiler_options] -- python your_script.py [your_script_args]

Example:
    python run_profiler.py --cprofile --output profile_results -- python ../run_test_GridFF_ocl_new.py
"""

import argparse
import cProfile
import io
import os
import pstats
import subprocess
import sys
import time
import datetime
from pathlib import Path

def parse_args():
    # First, check if -- is in the arguments
    if '--' not in sys.argv:
        print("Error: Command to profile must be specified after --")
        print("Example: python run_profiler.py --cprofile -- python your_script.py")
        sys.exit(1)
    
    # Split arguments before and after --
    dash_index = sys.argv.index('--')
    profiler_args = sys.argv[1:dash_index]
    command_args = sys.argv[dash_index+1:]
    
    # Check if there's a command after --
    if not command_args:
        print("Error: No command specified after --")
        print("Example: python run_profiler.py --cprofile -- python your_script.py")
        sys.exit(1)
    
    # Parse profiler arguments
    parser = argparse.ArgumentParser(description='Plug-and-play Python profiler')
    parser.add_argument('--cprofile', action='store_true', help='Enable cProfile profiling')
    parser.add_argument('--output', default='profile_results', help='Output directory for profiling results')
    parser.add_argument('--limit', type=int, default=20, help='Limit the number of functions in the output')
    parser.add_argument('--sort', default='cumulative', 
                       choices=['cumulative', 'time', 'calls', 'pcalls', 'name'], 
                       help='Sort order for profiling results')
    parser.add_argument('--opencl-info', action='store_true', help='Show OpenCL device information')
    
    args = parser.parse_args(profiler_args)
    args.command = command_args
    
    return args

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

def run_with_cprofile(command, output_dir, limit, sort_by):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Base filename for outputs
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_filename = os.path.join(
        output_dir, 
        f"profile_{timestamp}"
    )
    
    # Set up profiling
    profiler = cProfile.Profile()
    
    # Start profiling
    print(f"Running command with cProfile: {' '.join(command)}")
    start_time = time.time()
    profiler.enable()
    
    # Run the command
    process = subprocess.Popen(command)
    process.wait()
    
    # Stop profiling
    profiler.disable()
    end_time = time.time()
    
    # Calculate total time
    total_time = end_time - start_time
    
    # Save cProfile results
    s = io.StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats(sort_by)
    ps.print_stats(limit)
    
    with open(f"{base_filename}_profile.txt", 'w') as f:
        f.write(s.getvalue())
    
    # Print summary to console
    print(f"\n{'='*50}")
    print(f"Profiling results for: {' '.join(command)}")
    print(f"{'='*50}")
    print(f"Total execution time: {total_time:.2f} seconds")
    print(f"Exit code: {process.returncode}")
    print(f"\nTop {limit} functions by {sort_by} time:")
    ps.print_stats(limit)
    
    print(f"\nDetailed profiling results saved to {base_filename}_profile.txt")
    print(f"{'='*50}")
    
    return process.returncode

def run_with_time(command):
    # Just run the command and time it
    print(f"Running command with time measurement: {' '.join(command)}")
    start_time = time.time()
    
    # Run the command
    process = subprocess.Popen(command)
    process.wait()
    
    # Calculate and print total time
    end_time = time.time()
    total_time = end_time - start_time
    
    print(f"\n{'='*50}")
    print(f"Time measurement for: {' '.join(command)}")
    print(f"{'='*50}")
    print(f"Total execution time: {total_time:.2f} seconds")
    print(f"Exit code: {process.returncode}")
    
    return process.returncode

def main():
    args = parse_args()
    
    # Show OpenCL information if requested
    if args.opencl_info:
        show_opencl_info()
    
    # Set environment variables for OpenCL profiling
    os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
    os.environ['OPENCL_PROFILE'] = '1'
    
    # Run with profiling
    if args.cprofile:
        return run_with_cprofile(args.command, args.output, args.limit, args.sort)
    else:
        return run_with_time(args.command)

if __name__ == '__main__':
    sys.exit(main())
