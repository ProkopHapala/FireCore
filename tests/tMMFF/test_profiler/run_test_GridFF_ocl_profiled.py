import sys
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

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
sys.path.append("../../")
from pyBall.tests import ocl_GridFF_new as gff
import pyopencl as cl

# Profile the test execution
start_timing('test_gridFF_ocl')
name="NaCl_1x1_L1"
gff.test_gridFF_ocl(
    fname="data/xyz/"+name+".xyz",
    Element_Types_name="./data/ElementTypes.dat", 
    save_name="profile_run", 
    job="PLQ",
    bPlot=False  # Disable plotting for profiling
)
stop_timing('test_gridFF_ocl')

# Stop overall timing
stop_timing('total')

# Print profiling results
print("\n===== Profiling Results =====")
for name, data in sorted(timings.items(), key=lambda x: x[1].get('elapsed', 0), reverse=True):
    if 'elapsed' in data:
        print(f"{name}: {data['elapsed']:.2f} ms")

# Save results to file
with open('profile_results.txt', 'w') as f:
    for name, data in sorted(timings.items(), key=lambda x: x[1].get('elapsed', 0), reverse=True):
        if 'elapsed' in data:
            f.write(f"{name}: {data['elapsed']:.2f} ms\n")

# Don't show plots
# plt.show()
