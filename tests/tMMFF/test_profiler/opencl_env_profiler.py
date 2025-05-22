#!/usr/bin/env python3
"""
OpenCL Environment Profiler

This module provides a way to profile OpenCL applications by setting
the appropriate environment variables and running the target script.
It doesn't modify any of the original code.
"""

import os
import sys
import time
import datetime
import argparse
import subprocess
import glob
import re

# Function to set up OpenCL profiling environment variables
def setup_opencl_profiling():
    # Basic OpenCL profiling
    os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
    os.environ['OPENCL_PROFILE'] = '1'
    
    # NVIDIA-specific profiling
    os.environ['CUDA_PROFILE'] = '1'
    os.environ['CUDA_PROFILE_CSV'] = '1'
    os.environ['CUDA_PROFILE_LOG'] = 'cuda_profile_%p.log'
    os.environ['CUDA_PROFILE_CONFIG'] = 'cuda_profile_config.txt'
    
    # AMD-specific profiling
    os.environ['AMD_OCL_BUILD_OPTIONS_APPEND'] = '-cl-opt-disable'
    os.environ['AMD_OPENCL_PROFILER_ENABLE'] = '1'
    
    # Intel-specific profiling
    os.environ['INTEL_OPENCL_PROFILE'] = '1'
    
    # Create CUDA profile config file if it doesn't exist
    config_file = 'cuda_profile_config.txt'
    if not os.path.exists(config_file):
        with open(config_file, 'w') as f:
            f.write("gpustarttimestamp\n")
            f.write("gpuendtimestamp\n")
            f.write("method\n")
            f.write("gputime\n")
            f.write("cputime\n")
            f.write("occupancy\n")
    
    print("OpenCL profiling environment variables set")

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
    
    # Set up OpenCL profiling environment variables
    setup_opencl_profiling()
    
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
                cwd=cwd,
                env=os.environ
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
        
        # Process CUDA profile logs if they exist
        cuda_logs = glob.glob("cuda_profile_*.log")
        if cuda_logs:
            print(f"Found CUDA profile logs: {cuda_logs}")
            # Process and summarize CUDA profile logs
            process_cuda_logs(cuda_logs, os.path.join(output_dir, f"profile_{timestamp}_cuda.txt"))
        
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
        import traceback
        traceback.print_exc()
        return 1, None, None, None

# Function to process CUDA profile logs
def process_cuda_logs(log_files, output_file):
    kernel_stats = {}
    
    for log_file in log_files:
        try:
            with open(log_file, 'r') as f:
                lines = f.readlines()
                
                # Parse header to find column indices
                header = lines[0].strip().split(',')
                method_idx = header.index('method') if 'method' in header else -1
                gputime_idx = header.index('gputime') if 'gputime' in header else -1
                
                if method_idx >= 0 and gputime_idx >= 0:
                    # Process each line
                    for line in lines[1:]:
                        cols = line.strip().split(',')
                        if len(cols) > max(method_idx, gputime_idx):
                            method = cols[method_idx]
                            gputime = float(cols[gputime_idx])
                            
                            if method not in kernel_stats:
                                kernel_stats[method] = {'calls': 0, 'total_time': 0.0, 'min_time': float('inf'), 'max_time': 0.0}
                            
                            kernel_stats[method]['calls'] += 1
                            kernel_stats[method]['total_time'] += gputime
                            kernel_stats[method]['min_time'] = min(kernel_stats[method]['min_time'], gputime)
                            kernel_stats[method]['max_time'] = max(kernel_stats[method]['max_time'], gputime)
        except Exception as e:
            print(f"Error processing CUDA log file {log_file}: {e}")
    
    # Write kernel statistics to output file
    with open(output_file, 'w') as f:
        f.write("CUDA Kernel Profiling Summary\n")
        f.write("===========================\n\n")
        
        if kernel_stats:
            f.write(f"{'Kernel':<40} {'Calls':<8} {'Total Time (ms)':<15} {'Avg Time (ms)':<15} {'Min Time (ms)':<15} {'Max Time (ms)':<15}\n")
            f.write("-" * 110 + "\n")
            
            # Sort kernels by total time
            for kernel, stats in sorted(kernel_stats.items(), key=lambda x: x[1]['total_time'], reverse=True):
                avg_time = stats['total_time'] / stats['calls'] if stats['calls'] > 0 else 0
                f.write(f"{kernel:<40} {stats['calls']:<8} {stats['total_time']:<15.3f} {avg_time:<15.3f} {stats['min_time']:<15.3f} {stats['max_time']:<15.3f}\n")
        else:
            f.write("No kernel execution data found in CUDA profile logs.\n")
    
    print(f"CUDA kernel profiling summary saved to: {output_file}")

# Function to run GridFF with profiling
def run_gridff_with_profiling(job="PLQ", dataset="NaCl_1x1_L1", output_dir="profile_results"):
    # Set environment variables for GridFF
    os.environ['GRIDFF_JOB'] = job
    os.environ['GRIDFF_DATASET'] = dataset
    
    # Add FireCore project root to Python path
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    
    # Command to run GridFF test
    cmd = ["python", "-c", f"""import sys; sys.path.insert(0, '{project_root}'); from pyBall.tests import ocl_GridFF_new as gff; gff.test_gridFF_ocl(job='{job}', fname='./data/xyz/{dataset}.xyz')"""]
    
    # Run with profiling
    return run_with_profiling(cmd, output_dir)

# Main function
def main():
    parser = argparse.ArgumentParser(description='OpenCL Environment Profiler')
    parser.add_argument('--job', default='PLQ', help='Job type for GridFF (PLQ, Morse, etc.)')
    parser.add_argument('--dataset', default='NaCl_1x1_L1', help='Dataset name for GridFF')
    parser.add_argument('--output', default='profile_results', help='Output directory for profiling results')
    parser.add_argument('command', nargs='*', help='Command to profile (if not using GridFF)')
    
    args = parser.parse_args()
    
    # If a command is provided, run it with profiling
    if args.command:
        return run_with_profiling(args.command, args.output)[0]
    else:
        # Otherwise, run GridFF with profiling
        return run_gridff_with_profiling(args.job, args.dataset, args.output)[0]

if __name__ == '__main__':
    sys.exit(main())
