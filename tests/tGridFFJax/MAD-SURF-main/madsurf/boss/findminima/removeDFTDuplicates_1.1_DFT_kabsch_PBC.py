import os
import re
import shutil
import numpy as np
from ase import io
from ase import Atoms
from scipy.spatial.transform import Rotation 
from ase.io import read, write
from ase.io.aims import read_aims, write_aims
from ase.geometry import distance
# Set the base directory
base_dir = '.'
# Get a sorted list of all the locmin_N folders
folders = sorted([f for f in os.listdir(base_dir) if f.startswith('locmin_')])
# Define a list to contain identical structures
identical_structures = []
def kabsch_rmsd(A, B, pbc=None):
    """
    Calculate the root-mean-square deviation (RMSD) between two structures using Kabsch algorithm.
    """
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    A_centered = A - centroid_A
    B_centered = B - centroid_B
    H = np.dot(A_centered.T, B_centered)
    U, S, Vt = np.linalg.svd(H)
    V = Vt.T
    
    R = np.dot(U, Vt)

    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = np.dot(U, Vt)
    
    # Apply rotation matrix to B
    B_rotated = np.dot(B_centered, R.T) + centroid_A
    
    # Apply periodic boundary conditions if provided
    if pbc is not None:
        B_rotated = np.mod(B_rotated, pbc)

    rmsd = np.sqrt(np.mean((A_centered - np.dot(B_centered, R.T))**2))

    return rmsd, B_rotated

unique_file = os.path.join(base_dir, 'unique_entries_sorted.dat')
if not os.path.isfile(unique_file):
    # Open the output file
    output_file = open('identical_pairs.dat', 'w')
    
    # Loop over all pairs of structures in the locmin_N folders
    for i in range(len(folders)):
        for j in range(i + 1, len(folders)):
            folder1 = folders[i]
            folder2 = folders[j]
    
            # Define the path to the structure files
            file_path1 = os.path.join(base_dir, folder1, 'geometry_{}_opt.in'.format(folder1.split('_')[-1]))
            file_path2 = os.path.join(base_dir, folder2, 'geometry_{}_opt.in'.format(folder2.split('_')[-1]))
    
            try:
                # Load the structures
                atoms1 = read(file_path1)
                atoms2 = read(file_path2)
                positions_A = atoms1.positions
                positions_B = atoms2.positions
                pbc_A = atoms1.get_pbc()
                
                # Get the distance (RMSD) between the structure pair
                rmsd, B_rotated = kabsch_rmsd(positions_A, positions_B, pbc_A)
                
                
                #Save aligned structure for inspection (applies only for inspecting the last structures checked by the script)
                #atoms_rotated = atoms2.copy()
                #atoms_rotated.positions = B_rotated 
                #write('geometry_{}_opt_aligned.in'.format(folder2), atoms_rotated)
                
                if rmsd < 0.1:
                    output_file.write('RMSD between structures {} and {}: {:.4f} Angstroms\n'.format(folder1, folder2, rmsd))
             
                # Print the distance value
                print('RMSD between structures {} and {}: {:.4f} Angstroms'.format(folder1, folder2, rmsd))
    
            except FileNotFoundError:
                # Skip the current pair of structures if the file is not found
                print('Skipping pair {} and {}: geometry file not found'.format(folder1, folder2))
                continue
    
    output_file.close()
    
    ## Checking if a structure is identical with several other structures ##
    # read the output file and extract locmin entries
    with open("identical_pairs.dat", "r") as f:
        lines = f.readlines()
    
    locmin_entries = set()
    for line in lines:
        locmin_matches = re.findall(r"locmin_\d+", line)
        if locmin_matches:
            locmin_entries.update(locmin_matches)
    
    # create dictionary with each locmin entry as a key
    # and a list of identical entries as the value
    identical_entries = {}
    for locmin_entry in locmin_entries:
        identical_entries[locmin_entry] = [locmin_entry]
    
    # iterate over each line again and update the dictionary accordingly
    for line in lines:
        locmin_matches = re.findall(r"locmin_\d+", line)
        if locmin_matches and len(locmin_matches) == 2:
            entry1, entry2 = locmin_matches
            if entry1 in identical_entries and entry2 in identical_entries:
                identical_entries[entry1].append(entry2)
                identical_entries[entry2].append(entry1)
    
    # convert dictionary to list of sets with unique entries
    unique_identical_entries = []
    for entry_set in identical_entries.values():
        if set(entry_set) not in unique_identical_entries:
            unique_identical_entries.append(set(entry_set))
    
    # sort the unique identical entries by the smallest locmin_N entry in each set
    unique_identical_entries.sort(key=lambda entry_set: int(sorted(entry_set)[0].split("_")[1]))
    
    # write the unique identical entries to a file
    with open("identical_locmin.dat", "w") as outfile:
        for entry_set in unique_identical_entries:
            outfile.write(' '.join(sorted(entry_set)) + "\n")
    
    # extract all locmin_N entries from the identical_locmin file
    locmin_entries = set()
    with open("identical_locmin.dat", "r") as f:
        lines = f.readlines()
    for line in lines:
        locmin_matches = re.findall(r"locmin_\d+", line)
        if locmin_matches:
            locmin_entries.update(locmin_matches)
    
    # create a list of locmin_N entries that are missing from identical_entries
    missing_entries = []
    with open("Final_energies.dat", "r") as f:
        for line in f:
            N, energy,_ = line.split()
            locmin_entry = f"locmin_{N}"
            if locmin_entry not in locmin_entries:
                missing_entries.append(locmin_entry)
    
    # sort the missing entries by N
    missing_entries.sort(key=lambda x: int(x.split("_")[1]))
    
    # combine unique identical entries and missing entries into a list of unique entries
    unique_entries = []
    for entry_set in unique_identical_entries:
        first_entry = sorted(entry_set, key=lambda x: int(x.split("_")[1]))[0]
        unique_entries.append(first_entry)
    for missing_entry in missing_entries:
        unique_entries.append(missing_entry)
    
    # sort the unique entries by N
    unique_entries_sorted = sorted(unique_entries, key=lambda x: int(x.split("_")[1]))
    
    # write the sorted unique entries to a file
    with open("unique_entries_sorted.dat", "w") as outfile:
        for entry in unique_entries_sorted:
            outfile.write(entry + "\n")

    # check for folders that are not in unique_entries_sorted.dat
    unique_entries = []
    with open("unique_entries_sorted.dat", "r") as f:
        for line in f:
            unique_entries.append(line.strip())
    
    duplicates_dir = "duplicates"
    if not os.path.exists(duplicates_dir):
        os.makedirs(duplicates_dir)
    
    for folder in os.listdir("."):
        if folder.startswith("locmin_") and folder not in unique_entries:
            shutil.move(folder, os.path.join(duplicates_dir, folder))

else:
    # check for folders that are not in unique_entries_sorted.dat
    unique_entries = []
    with open("unique_entries_sorted.dat", "r") as f:
        for line in f:
            unique_entries.append(line.strip())
    
    duplicates_dir = "duplicates"
    if not os.path.exists(duplicates_dir):
        os.makedirs(duplicates_dir)
    
    for folder in os.listdir("."):
        if folder.startswith("locmin_") and folder not in unique_entries:
            shutil.move(folder, os.path.join(duplicates_dir, folder))