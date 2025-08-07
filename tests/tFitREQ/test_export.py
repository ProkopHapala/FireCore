#!/usr/bin/python3
import sys, os
import numpy as np
# Initialize import paths
sys.path.clear()
#sys.path.append("/home/prokop/git/FireCore-fitREQH")
sys.path.append("../../")
from pyBall import FitREQ as fit
from pyBall.FitREQ import exportAllSamplesToXYZ

ref_path = "/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/FitREQ/DFT_2D/"
#ref_path = "/home/prokophapala/Desktop/CARBSIS/DFT_ref_2D/"
#name = "H2O-D1_H2O-A1"
#name = "H2O-D1_H2O-A1"
name = "HCOOH-D1_HCOOH-A1"

donors    = [
'H2O-D1', 
'NH3-D1'
'CH2NH-D1',
'HCOOH-D1', 
'HCONH2-D1', 
#'C4H3NO2-D1', 
#'C4H5N-D1', 
] 
acceptors = [
'H2O-A1', 
'NH3-A1',
'CH2O-A1', 
#'CH2NH-A1',    # This makes it crash
'HCOOH-A1', 
'HCOOH-A2', 
#'HCONH2-A1', 
#'C4H3NO2-A1', 
#'C5H5N-A1', 
]

#ref_dirs = fit.combine_fragments( donors, acceptors, path=ref_path, ext=".xyz" )  ;print( "ref_dirs:\n", ref_dirs )
ref_dirs = fit.combine_fragments( donors, acceptors, path=ref_path )               ;print( "ref_dirs:\n", ref_dirs )


marks    = fit.concatenate_xyz_files( directories=ref_dirs, base_path=ref_path, fname='all.xyz', output_file='all.xyz' )
#marks    = fit.concatenate_xyz_files_flat( names=ref_dirs, base_path=ref_path, output_file='all.xyz' )

print( "marks:\n", marks )

#fname = "input_single.xyz"
#fname = ref_path +"/"+ name + "/all.xyz"
#fname = ref_path +"/"+ "/concatenated_all.xyz"
fname = 'all.xyz'
#fname = "input_2CH2NH.xyz"








fname = "input_example.xyz"
print(f"Reading data from {fname}")
Erefs, x0s = fit.read_xyz_data(fname)
print(f"Read {len(Erefs)} energy entries")


# --- load types using default resources
fit.loadTypes()


# # --- load samples into C++
# nsamp = fit.loadXYZ(fname)
# print(f"Loaded {nsamp} samples into C++")

# --- export all samples via C++
output_fname = "test_export_output.xyz"
print(f"Exporting all samples to {output_fname}")
#fit.exportAllSamplesToXYZ(output_fname)

exportAllSamplesToXYZ(output_fname)


print(f"Exported all samples to {output_fname}")

print("TEST PASSED")

