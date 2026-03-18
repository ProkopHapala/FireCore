import os
import re
import numpy as np
from math import sqrt 
import shutil
from ase import io
from ase import atoms
from ase.build import molecule
from ase.io import read
from ase.geometry import distance
from ase.visualize import view
# Set the base directory
base_dir = '.'
# Translation vectors based on Cu lattice 
a_1 = np.array([(3.632/sqrt(2)), 0, 0])
a_2 = np.array([((3.632/sqrt(2))/2), 2.224, 0]) 
a_3 = np.array([-((3.632/sqrt(2))/2), 2.224, 0])
a_4 = a_2 + a_3
a_5 = 2*a_2
# Get a sorted list of all the locmin_N folders
folders = sorted([f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f)) and f.startswith('locmin_')])
# Define a list to contain identical structures
equivalent_structures = []
unique_file = os.path.join(base_dir, 'equivalent_input_locmin.dat')
if not os.path.isfile(unique_file):
    # Open the output file
    output_file = open('equivalent_input_pairs.dat', 'w')
    
    # Loop over all pairs of structures in the locmin_N folders
    for i in range(len(folders)):
        for j in range(i + 1, len(folders)):
            folder1 = folders[i]
            folder2 = folders[j]
    
            # Define the path to the structure files
            file_path1 = os.path.join(base_dir, folder1, 'geometry.in'.format(folder1.split('_')[-1]))
            file_path2 = os.path.join(base_dir, folder2, 'geometry.in'.format(folder2.split('_')[-1]))
    
            try:
                    # Load the structures
                atoms1 = read(file_path1)
                atoms2 = read(file_path2)
                del atoms1[[atom.index for atom in atoms1 if atom.symbol=='Cu']]
                del atoms2[[atom.index for atom in atoms2 if atom.symbol=='Cu']]
                
                # Use the translate method of ASE here, this will avoid false positives:
                for vector in [-a_1, -a_2, -a_3, -a_4, -a_5, a_1, a_2, a_3, a_4, a_5]:
                    translated_atoms1 = atoms2.copy()
                    translated_atoms1.translate(vector)
                    if np.allclose(atoms1.positions, translated_atoms1.positions, atol=0.1):
                        output_file.write('Structures {} and {} are translationally equivalent\n'.format(folder1, folder2))
                    # Check for rotational symmetry for translated adsorbates
                    for angle in [120, 240]:
                        rotated_atoms1 = translated_atoms1.copy()
                        rotated_atoms1.rotate(angle, 'z', center='COP')
                        if np.allclose(atoms1.positions, rotated_atoms1.positions, atol=1.0):
                            output_file.write('Structures {} and the translated {} are rotationally equivalent at {} degrees\n'.format(folder1, folder2, angle))  

                # Check for rotational symmetry for untranslated adsorbates
                for angle in [120, 240]:
                    rotated_atoms2 = atoms2.copy()
                    rotated_atoms2.rotate(angle, 'z', center='COP')
                    if np.allclose(atoms1.positions, rotated_atoms2.positions, atol=1.0):
                        output_file.write('Structures {} and {} are rotationally equivalent at {} degrees\n'.format(folder1, folder2, angle))

                ## Check for translational symmetry (alternative method)
                # dist = np.linalg.norm(atoms1.positions - atoms2.positions)
                # translated_dist1 = np.linalg.norm((atoms1.positions) - (atoms1.positions + a))
                # translated_dist2 = np.linalg.norm((atoms1.positions) - (atoms1.positions + a_2))
                # if abs(dist - translated_dist1) < 0.1:
                #     output_file.write('Structures {} and {} are translationally equivalent\n'.format(folder1, folder2)) 
                # if abs(dist - translated_dist2) < 0.1:
                #     output_file.write('Structures {} and {} are translationally equivalent\n'.format(folder1, folder2))  
                
            except FileNotFoundError:
                # Skip the current pair of structures if the file is not found
                print('Skipping pair {} and {}: geometry file not found'.format(folder1, folder2))
                continue
    
    output_file.close()        
    
    ## Checking if a structure is identical with several other structures ##
    # read the output file and extract locmin entries
    def find(x, parents):
        if parents[x] == x:
            return x
        parents[x] = find(parents[x], parents)
        return parents[x]
    
    def union(x, y, parents):
        px, py = find(x, parents), find(y, parents)
        if px != py:
            parents[py] = px
    
    equivalent_pairs_file = open("equivalent_input_pairs.dat", "r")
    
    parents = {}
    for line in equivalent_pairs_file:
        structure_names = re.findall(r'locmin_\d+', line)
        if structure_names[0] not in parents:
            parents[structure_names[0]] = structure_names[0]
        for name in structure_names[1:]:
            if name not in parents:
                parents[name] = name
            union(structure_names[0], name, parents)
    
    equivalent_pairs_file.close()
    
    equivalent_structures = {}
    for name in parents:
        parent = find(name, parents)
        if parent not in equivalent_structures:
            equivalent_structures[parent] = set()
        equivalent_structures[parent].add(name)
    
    # sort the unique identical entries by the smallest locmin_N entry in each set
    equivalent_structures_list = list(map(list, equivalent_structures.values()))
    sorted_list = sorted(equivalent_structures_list, key=lambda x: int(x[0].split('_')[1]))
    
    # write the unique equivalent entries to a file
    with open("equivalent_input_locmin.dat", "w") as outfile:#a_2 = np.array[2 * x for x in a_1]

# vectors going diagonally
        for equivalent_set in sorted_list:
            outfile.write(' '.join(equivalent_set) + "\n")
    
    # Get a sorted list of all the locmin_N folders
    folders = sorted([f for f in os.listdir(base_dir) if f.startswith('locmin_')]) 
    
    # Create the "equivalents" directory if it does not exist
    if not os.path.exists('input_equivalents'):
        os.makedirs('input_equivalents')
    
    # Flatten the list of equivalent folders into a single list
    with open('equivalent_input_locmin.dat', 'r') as f:
        equivalent_set_print = [folder for line in f for folder in line.strip().split()[1:]]
    
    # Loop through the directories and move those that do not contain N values in the first column of the file
    for f in folders:
        if f.startswith('locmin_') and os.path.isdir(f):
            if f in equivalent_set_print:
                shutil.move(f, 'input_equivalents')

else:
    # Get a sorted list of all the locmin_N folders
    folders = sorted([f for f in os.listdir(base_dir) if f.startswith('locmin_')]) 
    
    # Create the "equivalents" directory if it does not exist
    if not os.path.exists('input_equivalents'):
        os.makedirs('input_equivalents')
    
    # Flatten the list of equivalent folders into a single list
    with open('equivalent_input_locmin.dat', 'r') as f:
        equivalent_set_print = [folder for line in f for folder in line.strip().split()[1:]]
    
    # Loop through the directories and move those that do not contain N values in the first column of the file
    for f in folders:
        if f.startswith('locmin_') and os.path.isdir(f):
            if f in equivalent_set_print:
                shutil.move(f, 'input_equivalents')