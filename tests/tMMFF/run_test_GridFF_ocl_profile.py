#!/usr/bin/env python3
"""
Modified version of run_test_GridFF_ocl_new.py for profiling
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add FireCore root directory to Python path
firecore_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, firecore_root)
print(f"Added FireCore root directory to Python path: {firecore_root}")

# Now import the module
from pyBall.tests import ocl_GridFF_new as gff

import pyopencl as cl

# Set the name of the structure to use
name = "NaCl_1x1_L1"

# Set the correct paths
xyz_path = os.path.join(firecore_root, "cpp/common_resources/xyz", f"{name}.xyz")
element_types_path = os.path.join(firecore_root, "cpp/common_resources/ElementTypes.dat")

# Create output directory
output_dir = os.path.join(os.getcwd(), "profile_output")
os.makedirs(output_dir, exist_ok=True)

print(f"XYZ file path: {xyz_path}")
print(f"Element types file path: {element_types_path}")
print(f"Output directory: {output_dir}")

# Monkey patch the path calculation in the module
original_basename = os.path.basename

def patched_basename(path):
    # If it's the XYZ file, return a special name
    if xyz_path in path:
        return "profile_output"
    # Otherwise, use the original function
    return original_basename(path)

# Replace the function
os.path.basename = patched_basename

# Run the GridFF OpenCL test
print(f"Running GridFF OpenCL test with structure: {name}")
gff.test_gridFF_ocl(
    fname=xyz_path,
    Element_Types_name=element_types_path,
    save_name="profile_output",
    job="PLQ"
)

print("GridFF OpenCL test completed successfully")
