#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 09:11:41 2023
@author: jestilj1
"""
import numpy as np
import sys
import os
from math import sqrt
# Check if the correct number of command line arguments are provided
if len(sys.argv) != 2:
    print("Usage: python uniqueSymmetry_fcc111_1.2.py /path/to/filename.dat")
    sys.exit(1)
# Get the full path from the command line argument
full_path = sys.argv[1]
# Extract the base filename without the extension
filename = os.path.splitext(os.path.basename(full_path))[0]
# load data #
data = np.loadtxt(full_path, usecols=(2,3,4,5,6,7,8)) # first two coordinates are iter and npts
# Define tolerance values for the symmetries #
t_trans = 0.4
t_rot = 25
t_energy = 0.01
# Translation vectors based on Cu lattice (7D so they match the variable array) 
a_1 = np.array([(3.632/sqrt(2)), 0, 0, 0, 0, 0, 0])
a_2 = np.array([((3.632/sqrt(2))/2), 2.224, 0, 0, 0, 0, 0]) 
a_3 = np.array([-((3.632/sqrt(2))/2), 2.224, 0, 0, 0, 0, 0])
a_4 = a_2 + a_3
a_5 = 2*a_2
x2h = round(2.224, 4) # half of unit cell in y-direction
# Define symmetries by means of unit cell lengths #
x1s = round(0.000, 4) # x low bound
x2s = round(0.000, 4) # y low bound
x1h = round(1.284, 4) # half of unit cell in x-direction
x2h = round(2.224, 4) # half of unit cell in y-direction
x1e = round(2.568, 4) # x high bound
x2e = round(4.448, 4) # y high bound
    
# Rounded variables for comparison-purposes #
x2sym1 = round((x2e/3), 4)
x2sym2 = round((x2e*(5/6)), 4)
x2sym3 = round(((x2e/6)+x2h), 4)
x2sym4 = round(((x2e/6)), 4)

# Helper function to filter out identical points #
def filter_identical(data, tol):
    unique_coords = []
    for coord in data:
        is_identical = False
        for unique_coord in unique_coords:
            if np.allclose(coord, unique_coord, atol=tol):
                is_identical = True
                break
        if not is_identical:
            unique_coords.append(coord)
    unique_data = np.array(unique_coords)
    return unique_data

# Function to place data points within periodic bounds if outside #
def return_periodic(data_input):
    data = np.copy(data_input)
    for i in range(len(data)):
        if data[i,0] > x1e:
            data[i,0] -= x1e
        elif data[i,0] < x1s:
            data[i,0] += x1e
        if data[i,1] > x2e:
            data[i,1] -= x2e
        elif data[i,1] < x2s:
            data[i,1] += x2e
        if data[i,3] > 180.0:
            data[i,3] -= 360.0
        elif data[i,3] < -180.0:
            data[i,3] += 360.0
        if data[i,4] > 180.0:
            data[i,4] -= 360.0
        elif data[i,4] < -180.0:
            data[i,4] += 360.0
        if data[i,5] > 180.0:
            data[i,5] -= 360.0
        elif data[i,5] < -180.0:
            data[i,5] += 360.0
    return data 
                 
## Do symmetry operations on the remaining unique data points ##
# Translations #
def translational_symmetry(array_input):
    array_1 = np.copy(array_input) + -a_1
    array_2 = np.copy(array_input) + -a_2
    array_3 = np.copy(array_input) + -a_3
    array_4 = np.copy(array_input) + -a_4
    array_5 = np.copy(array_input) + -a_5
    array_6 = np.copy(array_input) + a_1
    array_7 = np.copy(array_input) + a_2
    array_8 = np.copy(array_input) + a_3
    array_9 = np.copy(array_input) + a_4
    array_10 = np.copy(array_input) + a_5
      
    return array_1, array_2, array_3, array_4, array_5, array_6, array_7, array_8, array_9, array_10
            
# Rotations #
def rotational_symmetry(array_input, tol):
    array1 = np.copy(array_input)
    array2 = np.copy(array_input)
    if np.isclose(array1[0], x1s, atol=tol) and np.isclose(array1[1], x2s, atol=tol):  # First layer
        if array1[5] <= -60.0:
            array1[5] += 240.0
            array2[5] += 120.0
        elif array1[5] >= 60.0:
            array1[5] -= 240.0
            array2[5] -= 120.0
        elif array1[5] < 60.0 and array1[5] > -60.0:
            array1[5] -= 120.0
            array2[5] += 120.0            
    elif np.isclose(array1[0], x1h, atol=tol) and np.isclose(array1[1], x2h, atol=tol):        
        if array1[5] <= -60.0:       
            array1[5] += 240.0
            array2[5] += 120.0
        elif array1[5] >= 60.0:
            array1[5] -= 240.0
            array2[5] -= 120.0
        elif array1[5] < 60.0 and array1[5] > -60.0: 
            array1[5] -= 120.0
            array2[5] += 120.0           
    elif np.isclose(array1[0], x1s, atol=tol) and np.isclose(array1[1], x2e, atol=tol):
        if array1[5] <= -60.0:       
            array1[5] += 240.0
            array2[5] += 120.0
        elif array1[5] >= 60.0:
            array1[5] -= 240.0
            array2[5] -= 120.0
        elif array1[5] < 60.0 and array1[5] > -60.0:
            array1[5] -= 120.0
            array2[5] += 120.0 
    elif np.isclose(array1[0], x1e, atol=tol) and np.isclose(array1[1], x2e, atol=tol):        
        if array1[5] <= -60.0:       
            array1[5] += 240.0
            array2[5] += 120.0
        elif array1[5] >= 60.0:
            array1[5] -= 240.0
            array2[5] -= 120.0
        elif array1[5] < 60.0 and array1[5] > -60.0: 
            array1[5] -= 120.0
            array2[5] += 120.0 
    elif np.isclose(array1[0], x1e, atol=tol) and np.isclose(array1[1], x2s, atol=tol):        
        if array1[5] <= -60.0:       
            array1[5] += 240.0
            array2[5] += 120.0
        elif array1[5] >= 60.0:
            array1[5] -= 240.0
            array2[5] -= 120.0
        elif array1[5] < 60.0 and array1[5] > -60.0: 
            array1[5] -= 120.0
            array2[5] += 120.0 
    elif np.isclose(array1[0], x1s, atol=tol) and np.isclose(array1[1], x2sym1, atol=tol): # Second layer
        if array1[5] <= -60.0:       
            array1[5] += 240.0
            array2[5] += 120.0
        elif array1[5] >= 60.0:
            array1[5] -= 240.0
            array2[5] -= 120.0
        elif array1[5] < 60.0 and array1[5] > -60.0: 
            array1[5] -= 120.0
            array2[5] += 120.0 
    elif np.isclose(array1[0], x1h, atol=tol) and np.isclose(array1[1], x2sym2, atol=tol):        
        if array1[5] <= -60.0:       
            array1[5] += 240.0
            array2[5] += 120.0
        elif array1[5] >= 60.0:
            array1[5] -= 240.0
            array2[5] -= 120.0
        elif array1[5] < 60.0 and array1[5] > -60.0: 
            array1[5] -= 120.0
            array2[5] += 120.0 
    elif np.isclose(array1[0], x1s, atol=tol) and np.isclose(array1[1], x2sym3, atol=tol): # Third layer
        if array1[5] <= -60.0:       
            array1[5] += 240.0
            array2[5] += 120.0
        elif array1[5] >= 60.0:
            array1[5] -= 240.0
            array2[5] -= 120.0
        elif array1[5] < 60.0 and array1[5] > -60.0: 
            array1[5] -= 120.0
            array2[5] += 120.0 
    elif np.isclose(array1[0], x1h, atol=tol) and np.isclose(array1[1], x2sym4, atol=tol):        
        if array1[5] <= -60.0:       
            array1[5] += 240.0
            array2[5] += 120.0
        elif array1[5] >= 60.0:
            array1[5] -= 240.0
            array2[5] -= 120.0
        elif array1[5] < 60.0 and array1[5] > -60.0: 
            array1[5] -= 120.0
            array2[5] += 120.0 
    return array1, array2

## Check for matches between translationally and rotationally equivalent points in the set of unique data points and remove ##
data=np.round(data, 4)
#unique_data = filter_identical(data, t_trans) # uncomment to see effect of filtering before symmetry checks
periodic_unique_data = return_periodic(data) # checking if any variables have ended up beyond the periodic bounds and returns them inside.
new_unique_data = filter_identical(periodic_unique_data, t_trans) # removing duplicate points within t_trans tolerance
to_remove = np.zeros(len(new_unique_data), dtype=bool)

for i in range(len(new_unique_data)):
    trans_sym_1, trans_sym_2, trans_sym_3, trans_sym_4, trans_sym_5, trans_sym_6, trans_sym_7, trans_sym_8, trans_sym_9, trans_sym_10 = translational_symmetry(new_unique_data[i,:])
    rot_sym1, rot_sym2 =  rotational_symmetry(new_unique_data[i,:], 0.1)
    trans_rot_sym1, trans_rot_sym2, trans_rot_sym3, trans_rot_sym4, trans_rot_sym5, trans_rot_sym6, trans_rot_sym7, trans_rot_sym8, trans_rot_sym9, trans_rot_sym10 = translational_symmetry(rot_sym1)  
    trans_rot_sym11, trans_rot_sym12, trans_rot_sym13, trans_rot_sym14, trans_rot_sym15, trans_rot_sym16, trans_rot_sym17, trans_rot_sym18, trans_rot_sym19, trans_rot_sym20 = translational_symmetry(rot_sym2)
  
    for j in range(i+1,len(new_unique_data)):
        if np.allclose(new_unique_data[i,:3], new_unique_data[j,:3], atol=t_trans) and np.allclose(new_unique_data[i,3:], new_unique_data[j,3:], atol=t_rot): #remove duplicates based on rotations
            to_remove[j] = True
        
        if np.allclose(new_unique_data[j,:3], trans_sym_1[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_sym_1[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_sym_2[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_sym_2[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_sym_3[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_sym_3[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_sym_4[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_sym_4[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_sym_5[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_sym_5[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_sym_6[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_sym_6[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_sym_7[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_sym_7[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_sym_8[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_sym_8[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_sym_9[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_sym_9[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_sym_10[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_sym_10[3:], atol=t_rot):
            to_remove[j] = True

        if np.allclose(new_unique_data[j,:3], rot_sym1[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], rot_sym1[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], rot_sym2[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], rot_sym2[3:], atol=t_rot):
            to_remove[j] = True

        if np.allclose(new_unique_data[j,:3], trans_rot_sym1[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym1[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym2[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym2[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym3[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym3[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym4[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym4[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym5[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym5[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym6[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym6[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym7[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym7[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym8[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym8[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym9[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym9[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym10[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym10[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym11[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym11[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym12[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym12[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym13[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym13[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym14[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym14[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym15[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym15[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym16[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym16[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym17[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym17[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym18[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym18[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym19[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym19[3:], atol=t_rot):
            to_remove[j] = True
        if np.allclose(new_unique_data[j,:3], trans_rot_sym20[:3], atol=t_trans) and np.allclose(new_unique_data[j,3:], trans_rot_sym20[3:], atol=t_rot):
            to_remove[j] = True
            
new_unique_data = np.delete(new_unique_data, np.where(to_remove), axis=0)  
np.savetxt(filename + "_unique.dat", new_unique_data, fmt='%15.7f')