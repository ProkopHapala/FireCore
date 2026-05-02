#!/usr/bin/env python
"""
Test DFTB+ Python API with compiled library.
Based on documentation in DFTB_docs.md Level 2: Native Python API
"""

import sys
import os
import numpy as np

# Try to import the DFTB+ Python API
try:
    from dftbplus import DftbPlus
    print("✓ Successfully imported dftbplus.DftbPlus")
except ImportError as e:
    print(f"✗ Failed to import dftbplus: {e}")
    print("Make sure you ran: source ~/.bashrc")
    sys.exit(1)

# Library path from your compilation
LIB_PATH = '/home/prokophapala/git_SW/dftbplus/_build/src/dftbp/libdftbplus.so'

if not os.path.exists(LIB_PATH):
    print(f"✗ Library not found: {LIB_PATH}")
    sys.exit(1)

print(f"✓ Library found: {LIB_PATH}")

# Create minimal dftb_in.hsd for testing
hsd_content = """Geometry = xyzFormat {
  <<< "h2o.xyz"
}

ParserOptions {
  ParserVersion = 15
}

Hamiltonian = DFTB {
  Scc = Yes
  MaxAngularMomentum {
    O = "p"
    H = "s"
  }
  SlaterKosterFiles = Type2FileNames {
    Prefix = "/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/"
    Separator = "-"
    Suffix = ".skf"
  }
}

Options {
  WriteResultsTag = Yes
}
"""

# Create XYZ file
xyz_content = """3
H2O test
O      0.0000000000   0.0000000000   0.0000000000
H      0.9580000000   0.0000000000   0.0000000000
H     -0.2400000000   0.9270000000   0.0000000000
"""

with open('h2o.xyz', 'w') as f:
    f.write(xyz_content)
print("✓ Created h2o.xyz")

with open('dftb_in.hsd', 'w') as f:
    f.write(hsd_content)
print("✓ Created dftb_in.hsd")

# Simple test: Initialize calculator
try:
    calc = DftbPlus(libpath=LIB_PATH)
    print("✓ DftbPlus calculator initialized")
    
    # Check available methods
    print("\n--- Available methods ---")
    methods = [m for m in dir(calc) if not m.startswith('_')]
    for m in sorted(methods):
        print(f"  - {m}")
        
except Exception as e:
    print(f"✗ Failed to initialize calculator: {e}")
    sys.exit(1)

calc.close()
print("\n✓ Calculator initialized successfully!")
print("\nNote: The calculator reads dftb_in.hsd from current directory by default")
