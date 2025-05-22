import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")

from pyBall.tests import ocl_GridFF_new as gff

# Dictionary to store timing information
timings = {}

# Function to start timing
def start_timing(name):
    if name not in timings:
        timings[name] = {}
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
name="NaCl_1x1_L1"
print(f"Running {job} job on {name} dataset with profiling...")

# Profile data loading
start_timing('data_loading')
# Your actual test code here
gff.test_gridFF_ocl(
    fname="data/xyz/"+name+".xyz",
    Element_Types_name="./data/ElementTypes.dat", 
    save_name="profile_run", 
    job="PLQ",
    bPlot=False
)
stop_timing('data_loading')

# Stop total timing
stop_timing('total')

# Print profiling results
print("
===== Profiling Results =====")
for name, data in timings.items():
    if 'elapsed' in data:
        print(f"{name}: {data['elapsed']:.2f} ms")

# Save results to file
with open('profile_results.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Operation', 'Time (ms)'])
    for name, data in timings.items():
        if 'elapsed' in data:
            writer.writerow([name, f"{data['elapsed']:.2f}"])

print("
Results saved to profile_results.csv")

# Don't show plots in profiling mode
# plt.show()
