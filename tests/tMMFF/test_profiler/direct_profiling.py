import sys
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import csv

sys.path.append("../../")

# Simple timer class for direct profiling
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

# Global dictionary to store timing information
timings = {}

# Function to start timing
def start_timing(name):
    if name not in timings:
        timings[name] = []
    timings[name].append(Timer(name).start())

# Function to stop timing
def stop_timing(name):
    if name in timings and timings[name]:
        timer = timings[name][-1].stop()
        return timer.elapsed_ms()
    return 0

# Function to run the test with direct profiling
def run_test_with_profiling(name="NaCl_1x1_L1", job="PLQ"):
    # Import the necessary modules
    import subprocess
    
    # Create a modified version of run_test_GridFF_ocl_new.py with profiling code
    with open('run_test_GridFF_ocl_new.py', 'r') as f:
        original_code = f.read()
    
    # Create a profiled version of the script
    profiled_script = f"""import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")

from pyBall.tests import ocl_GridFF_new as gff

# Dictionary to store timing information
timings = {{}}

# Function to start timing
def start_timing(name):
    if name not in timings:
        timings[name] = {{}}
    timings[name]['start'] = time.time()

# Function to stop timing
def stop_timing(name):
    if name in timings and 'start' in timings[name]:
        timings[name]['end'] = time.time()
        timings[name]['elapsed'] = (timings[name]['end'] - timings[name]['start']) * 1000  # ms
        return timings[name]['elapsed']
    return 0

# Profile the whole test
start_timing('total')

# Run the test
name="{name}"
print(f"Running {{job}} job on {{name}} dataset with profiling...")

# Profile data loading
start_timing('data_loading')
# Your actual test code here
gff.test_gridFF_ocl(
    fname="data/xyz/"+name+".xyz",
    Element_Types_name="./data/ElementTypes.dat", 
    save_name="profile_run", 
    job="{job}",
    bPlot=False
)
stop_timing('data_loading')

# Stop total timing
stop_timing('total')

# Print profiling results
print("\n===== Profiling Results =====")
for name, data in timings.items():
    if 'elapsed' in data:
        print(f"{{name}}: {{data['elapsed']:.2f}} ms")

# Save results to file
with open('profile_results.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Operation', 'Time (ms)'])
    for name, data in timings.items():
        if 'elapsed' in data:
            writer.writerow([name, f"{{data['elapsed']:.2f}}"])

print("\nResults saved to profile_results.csv")

# Don't show plots in profiling mode
# plt.show()
"""
    
    # Write the profiled script to a temporary file
    with open('profiled_test.py', 'w') as f:
        f.write(profiled_script)
    
    # Run the profiled script
    print(f"\nRunning {job} job on {name} dataset with profiling...")
    start_timing('subprocess')
    result = subprocess.run(['python', 'profiled_test.py'], 
                           cwd=os.getcwd(),
                           capture_output=True, 
                           text=True)
    stop_timing('subprocess')
    
    # Print the output
    print(result.stdout)
    if result.stderr:
        print("Errors:")
        print(result.stderr)
    
    # Parse the results if available
    try:
        results = {}
        with open('profile_results.csv', 'r') as f:
            reader = csv.reader(f)
            next(reader)  # Skip header
            for row in reader:
                if len(row) >= 2:
                    results[row[0]] = float(row[1])
        return results
    except Exception as e:
        print(f"Error parsing results: {e}")
        return {}

# Function to display profiling results
def display_profiling_results(results, job, name):
    if not results:
        print("No profiling results available.")
        return
    
    print(f"\n===== Profiling Results for {job} job on {name} dataset =====\n")
    
    # Sort results by time
    sorted_results = sorted(results.items(), key=lambda x: x[1], reverse=True)
    
    # Print total execution time
    if "total" in results:
        total_time = results["total"]
        print(f"Total execution time: {total_time:.2f} ms")
    else:
        total_time = sum(results.values())
        print(f"Sum of all profiled operations: {total_time:.2f} ms")
    
    # Print detailed results
    print("\nDetailed timing breakdown:")
    print(f"{'Operation':<30} {'Time (ms)':<15} {'% of Total':<15}")
    print("-" * 60)
    
    for name, time_ms in sorted_results:
        percent = (time_ms / total_time) * 100 if total_time > 0 else 0
        print(f"{name:<30} {time_ms:<15.2f} {percent:<15.2f}%")
    
    # Create visualization
    plt.figure(figsize=(12, 8))
    
    # Pie chart of operation times
    plt.subplot(2, 1, 1)
    labels = [k for k, _ in sorted_results[:5]]  # Top 5 operations
    sizes = [t for _, t in sorted_results[:5]]
    
    if len(sorted_results) > 5:
        others_total = sum(t for _, t in sorted_results[5:])
        if others_total > 0:
            labels.append('Others')
            sizes.append(others_total)
    
    if sizes:  # Only create pie chart if we have data
        plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
        plt.axis('equal')
        plt.title('Operation Time Distribution')
    
    # Bar chart of operation times
    plt.subplot(2, 1, 2)
    operations = [k for k, _ in sorted_results[:10]]  # Top 10 operations
    times = [t for _, t in sorted_results[:10]]
    
    if operations:  # Only create bar chart if we have data
        y_pos = range(len(operations))
        plt.barh(y_pos, times)
        plt.yticks(y_pos, operations)
        plt.xlabel('Time (ms)')
        plt.title('Top 10 Operation Times')
    
    plt.tight_layout()
    plt.savefig(f"profile_results_{job}_{name}.png")
    print(f"\nVisualization saved to profile_results_{job}_{name}.png")

# Function to add profiling to the GridFF.cl file
def add_profiling_to_opencl():
    # Path to the OpenCL file
    cl_file_path = "../../cpp/common_resources/cl/GridFF.cl"
    
    try:
        # Read the original OpenCL file
        with open(cl_file_path, 'r') as f:
            cl_code = f.read()
        
        # Create a backup
        with open(cl_file_path + '.backup', 'w') as f:
            f.write(cl_code)
        
        print(f"Created backup of {cl_file_path} as {cl_file_path}.backup")
        
        # We can't directly modify the OpenCL file for profiling as it would require
        # significant changes to the kernel code and host code. Instead, we'll focus
        # on host-side profiling.
        
        print("Note: Direct OpenCL kernel profiling requires modifications to both")
        print("the kernel code and host code. We'll focus on host-side profiling.")
        
        return True
    except Exception as e:
        print(f"Error backing up OpenCL file: {e}")
        return False

# Function to add profiling to the GridFF_cl class
def add_profiling_to_gridff_cl():
    # Path to the GridFF.py file
    py_file_path = "../../pyBall/OCL/GridFF.py"
    
    try:
        # Read the original Python file
        with open(py_file_path, 'r') as f:
            py_code = f.read()
        
        # Create a backup
        with open(py_file_path + '.backup', 'w') as f:
            f.write(py_code)
        
        print(f"Created backup of {py_file_path} as {py_file_path}.backup")
        
        # We can't directly modify the Python file for profiling as it would require
        # significant changes to the class implementation. Instead, we'll focus
        # on external profiling.
        
        print("Note: Direct Python class profiling requires modifications to the")
        print("class implementation. We'll focus on external profiling.")
        
        return True
    except Exception as e:
        print(f"Error backing up Python file: {e}")
        return False

# Main function
def main():
    # Backup files before profiling
    add_profiling_to_opencl()
    add_profiling_to_gridff_cl()
    
    # Run profiling for PLQ job
    print("\nRunning profiling for PLQ job...")
    results_plq = run_test_with_profiling(name="NaCl_1x1_L1", job="PLQ")
    display_profiling_results(results_plq, "PLQ", "NaCl_1x1_L1")
    
    # Uncomment to run Ewald profiling
    # print("\nRunning profiling for Ewald job...")
    # results_ewald = run_test_with_profiling(name="NaCl_1x1_L1", job="Ewald")
    # display_profiling_results(results_ewald, "Ewald", "NaCl_1x1_L1")
    
    plt.show()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error during profiling: {e}")
        import traceback
        traceback.print_exc()
