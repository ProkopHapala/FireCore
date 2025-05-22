#!/usr/bin/env python3
"""
NVIDIA OpenCL Profiler for GridFF

This script provides profiling capabilities for OpenCL code running on NVIDIA GPUs
by using manual timing instrumentation and environment variables.
"""

import os
import sys
import time
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Set environment variables for OpenCL profiling
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
os.environ['OPENCL_PROFILE'] = '1'  # Enable OpenCL profiling if supported

# Simple timer class for manual profiling
class Timer:
    def __init__(self, name):
        self.name = name
        self.start_time = None
        self.end_time = None
    
    def start(self):
        self.start_time = time.time()
        return self
    
    def stop(self):
        self.end_time = time.time()
        return self
    
    def elapsed_ms(self):
        if self.start_time is None:
            return 0
        end = self.end_time if self.end_time is not None else time.time()
        return (end - self.start_time) * 1000  # Convert to milliseconds

def run_profiled_test(job="PLQ", name="NaCl_1x1_L1", output_dir="./profile_results"):
    """Run a profiled test by executing the original script and capturing timing information"""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a modified version of run_test_GridFF_ocl_new.py with timing code
    profiled_script = f"""
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")

# Start timing
start_time = time.time()

# Set environment variables for OpenCL profiling
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
os.environ['OPENCL_PROFILE'] = '1'  # Enable OpenCL profiling if supported

# Import the test module
from pyBall.tests import ocl_GridFF_new as gff

# Dictionary to store timing information
timings = {{}}

# Function to start timing
def start_timing(name):
    timings[name] = {{'start': time.time()}}

# Function to stop timing
def stop_timing(name):
    if name in timings:
        timings[name]['end'] = time.time()
        timings[name]['elapsed'] = (timings[name]['end'] - timings[name]['start']) * 1000  # ms
        return timings[name]['elapsed']
    return 0

# Profile data loading
start_timing('data_loading')
name = "{name}"
job = "{job}"
stop_timing('data_loading')

# Profile test execution
start_timing('test_execution')
gff.test_gridFF_ocl(
    fname="data/xyz/"+name+".xyz",
    Element_Types_name="./data/ElementTypes.dat", 
    save_name="profile_run", 
    job=job,
    bPlot=False  # Disable plotting for profiling
)
stop_timing('test_execution')

# Calculate total time
total_time = (time.time() - start_time) * 1000  # ms
print(f"\\nTotal execution time: {{total_time:.2f}} ms")

# Print timing results
print("\\n===== Timing Results =====")
for name, data in timings.items():
    if 'elapsed' in data:
        print(f"{{name}}: {{data['elapsed']:.2f}} ms")

# Save timing results to file
with open('{output_dir}/timing_results.txt', 'w') as f:
    f.write(f"Total execution time: {{total_time:.2f}} ms\\n")
    for name, data in timings.items():
        if 'elapsed' in data:
            f.write(f"{{name}}: {{data['elapsed']:.2f}} ms\\n")

# Don't show plots
# plt.show()
"""
    
    # Write the profiled script to a temporary file
    with open(f"{output_dir}/profiled_test.py", 'w') as f:
        f.write(profiled_script)
    
    # Run the profiled script and capture output
    print(f"Running profiled test for {job} job on {name} dataset...")
    
    # Start timing the subprocess execution
    timer = Timer("subprocess_execution").start()
    
    # Run the profiled script
    process = subprocess.Popen(
        ["python", f"{output_dir}/profiled_test.py"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    # Capture output in real-time
    stdout_lines = []
    stderr_lines = []
    
    while True:
        stdout_line = process.stdout.readline()
        stderr_line = process.stderr.readline()
        
        if stdout_line:
            print(stdout_line.strip())
            stdout_lines.append(stdout_line)
        
        if stderr_line:
            print(stderr_line.strip(), file=sys.stderr)
            stderr_lines.append(stderr_line)
        
        if not stdout_line and not stderr_line and process.poll() is not None:
            break
    
    # Get remaining output
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.strip())
        stdout_lines.append(stdout)
    if stderr:
        print(stderr.strip(), file=sys.stderr)
        stderr_lines.append(stderr)
    
    # Stop timing
    timer.stop()
    
    # Save full output
    with open(f"{output_dir}/stdout.log", 'w') as f:
        f.writelines(stdout_lines)
    
    with open(f"{output_dir}/stderr.log", 'w') as f:
        f.writelines(stderr_lines)
    
    # Parse timing results
    timing_results = {}
    try:
        with open(f"{output_dir}/timing_results.txt", 'r') as f:
            for line in f:
                if ':' in line:
                    name, time_str = line.split(':', 1)
                    time_ms = float(time_str.strip().split()[0])
                    timing_results[name.strip()] = time_ms
    except FileNotFoundError:
        print(f"Warning: Timing results file not found. The profiled script may have failed.")
        timing_results["Total execution time"] = timer.elapsed_ms()
    
    # Add subprocess execution time
    timing_results["Subprocess execution"] = timer.elapsed_ms()
    
    return timing_results

def analyze_opencl_info():
    """Analyze OpenCL information on the system"""
    print("\n===== OpenCL Information =====")
    
    try:
        # Create a simple script to print OpenCL information
        opencl_script = """
import pyopencl as cl

print("\\nOpenCL Platforms and Devices:")
platforms = cl.get_platforms()
for i, platform in enumerate(platforms):
    print(f"  Platform {i}: {platform.name} ({platform.version})")
    devices = platform.get_devices()
    for j, device in enumerate(devices):
        print(f"    Device {j}: {device.name} (Type: {cl.device_type.to_string(device.type)})")
        print(f"      Compute Units: {device.max_compute_units}")
        print(f"      Global Memory: {device.global_mem_size / (1024**2):.2f} MB")
        print(f"      Local Memory: {device.local_mem_size / 1024:.2f} KB")
        print(f"      Max Work Group Size: {device.max_work_group_size}")
"""
        
        # Write the script to a temporary file
        with open("opencl_info.py", 'w') as f:
            f.write(opencl_script)
        
        # Run the script
        subprocess.run(["python", "opencl_info.py"], check=True)
        
        # Clean up
        os.remove("opencl_info.py")
    except Exception as e:
        print(f"Error analyzing OpenCL information: {e}")

def visualize_results(timing_results, job, name, output_dir="./profile_results"):
    """Visualize the timing results"""
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Sort results by time
    sorted_results = sorted(timing_results.items(), key=lambda x: x[1], reverse=True)
    
    # Bar chart of execution times
    labels = [item[0] for item in sorted_results]
    times = [item[1] for item in sorted_results]
    
    y_pos = range(len(labels))
    plt.barh(y_pos, times)
    plt.yticks(y_pos, labels)
    plt.xlabel('Time (ms)')
    plt.title(f'Execution Times - {job} job on {name} dataset')
    
    # Add time labels
    for i, v in enumerate(times):
        plt.text(v + 0.1, i, f"{v:.2f} ms", va='center')
    
    # Save figure
    plt.tight_layout()
    plt.savefig(f"{output_dir}/profile_{job}_{name}.png")
    print(f"\nVisualization saved to {output_dir}/profile_{job}_{name}.png")
    
    # Show figure
    plt.show()

def main():
    """Main function"""
    # Parse command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Profile OpenCL code running on NVIDIA GPU')
    parser.add_argument('--job', type=str, default='PLQ', help='Job type (PLQ, Morse, Ewald)')
    parser.add_argument('--name', type=str, default='NaCl_1x1_L1', help='Dataset name')
    parser.add_argument('--output-dir', type=str, default='./profile_results', help='Output directory')
    args = parser.parse_args()
    
    # Analyze OpenCL information
    analyze_opencl_info()
    
    # Run profiled test
    timing_results = run_profiled_test(args.job, args.name, args.output_dir)
    
    # Print timing results
    print("\n===== Timing Results =====")
    for name, time_ms in sorted(timing_results.items(), key=lambda x: x[1], reverse=True):
        print(f"{name}: {time_ms:.2f} ms")
    
    # Visualize results
    visualize_results(timing_results, args.job, args.name, args.output_dir)

if __name__ == "__main__":
    main()
