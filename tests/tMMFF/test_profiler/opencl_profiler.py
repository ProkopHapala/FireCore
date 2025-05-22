#!/usr/bin/env python3
"""
Simple OpenCL Profiler

This module provides a simple, reliable way to profile OpenCL kernels
without modifying your existing code.
"""

import os
import sys
import time
import datetime
import argparse
import subprocess
import traceback

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

# Function to enable OpenCL profiling via environment variables
def enable_opencl_profiling():
    # Set environment variables for OpenCL profiling
    os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
    os.environ['OPENCL_PROFILE'] = '1'
    
    # For NVIDIA
    os.environ['CUDA_PROFILE'] = '1'
    os.environ['CUDA_PROFILE_CSV'] = '1'
    
    # For AMD
    os.environ['AMD_OCL_BUILD_OPTIONS_APPEND'] = '-cl-opt-disable'
    
    # For Intel
    os.environ['INTEL_OPENCL_PROFILE'] = '1'
    
    print("OpenCL profiling enabled via environment variables")

# Function to run a command with profiling
def run_with_profiling(cmd, output_dir="profile_results", cwd=None):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate timestamp for output files
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Output files
    stdout_file = os.path.join(output_dir, f"profile_{timestamp}_stdout.txt")
    stderr_file = os.path.join(output_dir, f"profile_{timestamp}_stderr.txt")
    summary_file = os.path.join(output_dir, f"profile_{timestamp}_summary.txt")
    
    # Enable OpenCL profiling
    enable_opencl_profiling()
    
    # Show OpenCL information
    show_opencl_info()
    
    # Print the command being executed
    print(f"\nRunning command with profiling: {' '.join(cmd)}")
    
    # Record start time
    start_time = time.time()
    
    # Run the command and capture output
    try:
        with open(stdout_file, 'w') as stdout_f, open(stderr_file, 'w') as stderr_f:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                cwd=cwd
            )
            
            # Process output in real-time
            while True:
                stdout_line = process.stdout.readline()
                stderr_line = process.stderr.readline()
                
                if stdout_line == '' and stderr_line == '' and process.poll() is not None:
                    break
                
                if stdout_line:
                    print(stdout_line.rstrip())
                    stdout_f.write(stdout_line)
                    stdout_f.flush()
                
                if stderr_line:
                    print(stderr_line.rstrip(), file=sys.stderr)
                    stderr_f.write(stderr_line)
                    stderr_f.flush()
            
            # Get return code
            return_code = process.poll()
        
        # Record end time
        end_time = time.time()
        total_time = end_time - start_time
        
        # Write summary
        with open(summary_file, 'w') as f:
            f.write(f"OpenCL Profiling Summary\n")
            f.write(f"======================\n\n")
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"Working Directory: {cwd or os.getcwd()}\n")
            f.write(f"Timestamp: {timestamp}\n")
            f.write(f"Exit Code: {return_code}\n")
            f.write(f"Total Execution Time: {total_time:.2f} seconds\n\n")
            
            # Add OpenCL environment variables
            f.write(f"OpenCL Environment Variables:\n")
            for var in ['PYOPENCL_COMPILER_OUTPUT', 'OPENCL_PROFILE', 'CUDA_PROFILE', 
                       'CUDA_PROFILE_CSV', 'AMD_OCL_BUILD_OPTIONS_APPEND', 'INTEL_OPENCL_PROFILE']:
                f.write(f"  {var}={os.environ.get(var, 'Not set')}\n")
        
        print(f"\n{'='*50}")
        print(f"Profiling complete")
        print(f"Total execution time: {total_time:.2f} seconds")
        print(f"Exit code: {return_code}")
        print(f"Output saved to:")
        print(f"  Standard output: {stdout_file}")
        print(f"  Standard error: {stderr_file}")
        print(f"  Summary: {summary_file}")
        
        return return_code, stdout_file, stderr_file, summary_file
    
    except Exception as e:
        print(f"Error running command: {e}")
        traceback.print_exc()
        return 1, None, None, None

# Main function
def main():
    parser = argparse.ArgumentParser(description='Simple OpenCL Profiler')
    parser.add_argument('command', nargs='+', help='Command to profile')
    parser.add_argument('--output', default='profile_results', help='Output directory for profiling results')
    parser.add_argument('--cwd', help='Working directory for the command')
    parser.add_argument('--no-opencl-info', action='store_true', help='Skip OpenCL information display')
    
    args = parser.parse_args()
    
    # Run the command with profiling
    return run_with_profiling(args.command, args.output, args.cwd)[0]

if __name__ == '__main__':
    sys.exit(main())
