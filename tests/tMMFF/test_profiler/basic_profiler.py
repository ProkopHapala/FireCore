#!/usr/bin/env python3
"""
Basic Profiler for GridFF OpenCL Code

This script provides basic profiling capabilities for your GridFF OpenCL code
by instrumenting the run_test_GridFF_ocl_new.py script with timing code.
"""

import sys
import os
import time
import subprocess
import matplotlib.pyplot as plt
import numpy as np

def profile_command(cmd, cwd=None, env=None):
    """Run a command and measure its execution time"""
    start_time = time.time()
    
    # Run the command
    process = subprocess.Popen(
        cmd,
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    # Capture output in real-time
    stdout_lines = []
    stderr_lines = []
    
    print(f"Running command: {' '.join(cmd)}")
    
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
    
    # Calculate execution time
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    return {
        'elapsed_time': elapsed_time,
        'stdout': ''.join(stdout_lines),
        'stderr': ''.join(stderr_lines),
        'exit_code': process.returncode
    }

def create_instrumented_script(original_script, instrumented_script, job="PLQ", name="NaCl_1x1_L1"):
    """Create an instrumented version of the original script with timing code"""
    with open(original_script, 'r') as f:
        lines = f.readlines()
    
    # Add timing imports and functions
    instrumented_code = """import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

# Dictionary to store timing information
timings = {}

# Function to start timing
def start_timing(name):
    timings[name] = {'start': time.time()}

# Function to stop timing
def stop_timing(name):
    if name in timings:
        timings[name]['end'] = time.time()
        timings[name]['elapsed'] = (timings[name]['end'] - timings[name]['start']) * 1000  # ms
        return timings[name]['elapsed']
    return 0

# Start overall timing
start_timing('total')

"""
    
    # Find the line that calls test_gridFF_ocl
    test_line_idx = -1
    for i, line in enumerate(lines):
        if "test_gridFF_ocl" in line and "fname=" in line:
            test_line_idx = i
            break
    
    if test_line_idx >= 0:
        # Add the original imports
        for i in range(test_line_idx):
            if "import" in lines[i] or "sys.path.append" in lines[i]:
                instrumented_code += lines[i]
        
        # Add profiling code
        instrumented_code += f"""
# Profile the test execution
start_timing('test_gridFF_ocl')
name="{name}"
gff.test_gridFF_ocl(
    fname="data/xyz/"+name+".xyz",
    Element_Types_name="./data/ElementTypes.dat", 
    save_name="profile_run", 
    job="{job}",
    bPlot=False  # Disable plotting for profiling
)
stop_timing('test_gridFF_ocl')

# Stop overall timing
stop_timing('total')

# Print profiling results
print("\\n===== Profiling Results =====")
for name, data in sorted(timings.items(), key=lambda x: x[1].get('elapsed', 0), reverse=True):
    if 'elapsed' in data:
        print(f"{{name}}: {{data['elapsed']:.2f}} ms")

# Save results to file
with open('profile_results.txt', 'w') as f:
    for name, data in sorted(timings.items(), key=lambda x: x[1].get('elapsed', 0), reverse=True):
        if 'elapsed' in data:
            f.write(f"{{name}}: {{data['elapsed']:.2f}} ms\\n")

# Don't show plots
# plt.show()
"""
    else:
        # If we can't find the test line, just use the original script
        instrumented_code += ''.join(lines)
    
    # Write the instrumented script
    with open(instrumented_script, 'w') as f:
        f.write(instrumented_code)
    
    return instrumented_script

def run_profiling(job="PLQ", name="NaCl_1x1_L1"):
    """Run profiling on the specified job and dataset"""
    # Create instrumented script
    original_script = "run_test_GridFF_ocl_new.py"
    instrumented_script = "run_test_GridFF_ocl_profiled.py"
    
    create_instrumented_script(original_script, instrumented_script, job, name)
    
    # Run the instrumented script
    result = profile_command(["python", instrumented_script])
    
    # Parse profiling results
    profiling_results = {}
    
    try:
        with open('profile_results.txt', 'r') as f:
            for line in f:
                if ':' in line:
                    parts = line.strip().split(':', 1)
                    if len(parts) == 2:
                        name = parts[0].strip()
                        time_str = parts[1].strip().split()[0]
                        try:
                            time_ms = float(time_str)
                            profiling_results[name] = time_ms
                        except ValueError:
                            pass
    except FileNotFoundError:
        # If the file doesn't exist, use the overall execution time
        profiling_results['total'] = result['elapsed_time'] * 1000  # ms
    
    # Add the overall execution time if not already present
    if 'total' not in profiling_results:
        profiling_results['total'] = result['elapsed_time'] * 1000  # ms
    
    return profiling_results

def analyze_system_info():
    """Analyze system information"""
    print("\n===== System Information =====")
    
    # Get CPU info
    try:
        with open('/proc/cpuinfo', 'r') as f:
            cpu_info = f.read()
        
        # Extract CPU model
        for line in cpu_info.split('\n'):
            if 'model name' in line:
                print(f"CPU: {line.split(':', 1)[1].strip()}")
                break
    except:
        print("Could not determine CPU information")
    
    # Get memory info
    try:
        with open('/proc/meminfo', 'r') as f:
            mem_info = f.read()
        
        # Extract total memory
        for line in mem_info.split('\n'):
            if 'MemTotal' in line:
                print(f"Memory: {line.split(':', 1)[1].strip()}")
                break
    except:
        print("Could not determine memory information")
    
    # Check for NVIDIA GPU
    try:
        nvidia_smi = subprocess.run(['nvidia-smi'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if nvidia_smi.returncode == 0:
            print("\nNVIDIA GPU detected:")
            for line in nvidia_smi.stdout.split('\n'):
                if '|' in line and ('NVIDIA' in line or 'GeForce' in line or 'Tesla' in line or 'Quadro' in line):
                    print(f"GPU: {line.strip()}")
    except:
        print("NVIDIA GPU not detected or nvidia-smi not available")
    
    # Check OpenCL
    try:
        # Create a simple script to check OpenCL
        with open('check_opencl.py', 'w') as f:
            f.write("""
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
""")
        
        # Run the script
        subprocess.run(['python', 'check_opencl.py'])
        
        # Clean up
        os.remove('check_opencl.py')
    except:
        print("Could not check OpenCL information")

def visualize_results(profiling_results, job, name):
    """Visualize profiling results"""
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Sort results by time
    sorted_results = sorted(profiling_results.items(), key=lambda x: x[1], reverse=True)
    
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
    plt.savefig(f"profile_{job}_{name}.png")
    print(f"\nVisualization saved to profile_{job}_{name}.png")
    
    # Show figure
    plt.show()

def profile_with_time_command(job="PLQ", name="NaCl_1x1_L1"):
    """Profile using the time command"""
    print("\n===== Profiling with time command =====")
    
    # Create a script that runs the test with the job and name
    with open('run_test_with_params.py', 'w') as f:
        f.write(f"""
import sys
sys.path.append("../../")
from pyBall.tests import ocl_GridFF_new as gff

# Run the test
gff.test_gridFF_ocl(
    fname="data/xyz/{name}.xyz",
    Element_Types_name="./data/ElementTypes.dat", 
    save_name="profile_run", 
    job="{job}",
    bPlot=False  # Disable plotting for profiling
)
""")
    
    # Run the script with time command
    result = profile_command(["time", "python", "run_test_with_params.py"])
    
    # Clean up
    os.remove('run_test_with_params.py')
    
    return result

def main():
    """Main function"""
    # Parse command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Profile GridFF OpenCL code')
    parser.add_argument('--job', type=str, default='PLQ', help='Job type (PLQ, Morse, Ewald)')
    parser.add_argument('--name', type=str, default='NaCl_1x1_L1', help='Dataset name')
    args = parser.parse_args()
    
    # Analyze system information
    analyze_system_info()
    
    # Run profiling
    print(f"\n===== Profiling {args.job} job on {args.name} dataset =====")
    
    # Profile with time command first
    time_result = profile_with_time_command(args.job, args.name)
    
    # Run instrumented profiling
    profiling_results = run_profiling(args.job, args.name)
    
    # Print profiling results
    print("\n===== Profiling Results =====")
    for name, time_ms in sorted(profiling_results.items(), key=lambda x: x[1], reverse=True):
        print(f"{name}: {time_ms:.2f} ms")
    
    # Visualize results
    visualize_results(profiling_results, args.job, args.name)

if __name__ == "__main__":
    main()
