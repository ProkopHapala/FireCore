#!/usr/bin/env python3
"""
Run FireCore test with enhanced OpenCL profiling
"""

import sys
import os
import importlib.util

# Add FireCore root directory to Python path
sys.path.append("../../")

# Import the GridFF_cl class
from pyBall.OCL import GridFF

# Run the test script
print("Running FireCore test script...")
spec = importlib.util.spec_from_file_location("test_script", "run_test_GridFF_ocl_new.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

# Save profiling results
print("Saving profiling results...")
output_dir = os.path.join(os.getcwd(), "firecore_enhanced_profile_results")
os.makedirs(output_dir, exist_ok=True)
GridFF.save_profiling_results(output_dir)

print("Profiling completed successfully")
