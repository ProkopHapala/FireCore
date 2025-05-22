#!/usr/bin/env python3
"""
Wrapper script for profiling GridFF code

This script provides a convenient way to profile GridFF code using the general-purpose profiler.

Usage:
    python profile_grid.py [options]
"""

import argparse
import os
import sys
import subprocess

def main():
    parser = argparse.ArgumentParser(description='Profile GridFF code')
    parser.add_argument('--job', type=str, default='PLQ', help='Job type (PLQ, Morse, Ewald)')
    parser.add_argument('--dataset', type=str, default='NaCl_1x1_L1', help='Dataset name')
    parser.add_argument('--trace-modules', type=str, default='pyBall', 
                       help='Comma-separated list of module prefixes to trace')
    parser.add_argument('--no-opencl', action='store_true', help='Disable OpenCL profiling')
    parser.add_argument('--no-trace', action='store_true', help='Disable function tracing')
    parser.add_argument('--output-dir', default='./profile_results', help='Directory to save results')
    
    args = parser.parse_args()
    
    # Path to the profiler script
    profiler_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'simple_profiler.py')
    
    # Path to the GridFF script
    grid_script = os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 
        '..', 
        'run_test_GridFF_ocl_new.py'
    ))
    
    # Build the command
    cmd = [
        sys.executable,  # Current Python interpreter
        profiler_path,
    ]
    
    # Add profiling options
    if not args.no_opencl:
        cmd.append('--opencl')
    if not args.no_trace:
        cmd.append('--trace')
        cmd.append('--trace-modules')
        cmd.append(args.trace_modules)
    
    cmd.append('--output-dir')
    cmd.append(args.output_dir)
    
    # Add the script to profile
    cmd.append(grid_script)
    
    # Check if the grid script accepts command line arguments
    # For now, we'll modify the environment variables instead of passing command line args
    # since run_test_GridFF_ocl_new.py might not accept command line arguments
    os.environ['FIRECORE_JOB'] = args.job
    os.environ['FIRECORE_DATASET'] = args.dataset
    
    print(f"\nProfiling GridFF with job={args.job}, dataset={args.dataset}")
    
    # Print the command being run
    print(f"Running: {' '.join(cmd)}")
    
    # Run the profiler
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running profiler: {e}")
        return 1
    except KeyboardInterrupt:
        print("Profiling interrupted by user")
        return 130
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
