from ase.io import read, write
from sys import argv
import os

# Load the extxyz file
traj = read(argv[1], index=':')

# Define the indices of the structures you want to save in the subset
indices_to_save = list(range(0, len(traj), 10))  # Adjust this list as per your requirements

# Create a new trajectory with the subset of structures
test_set_traj = [traj[i] for i in indices_to_save]
train_set_traj = [traj[i] for i in range(len(traj)) if i not in indices_to_save]

# Output file name
output_file_test = 'test_set.extxyz'
output_file_train = 'train_set.extxyz'

# Check if the append option is provided as a command line argument
if len(argv) > 2 and argv[2].lower() == '--append':
    # Append to the existing file
    write(output_file_test, test_set_traj, append=True)
    write(output_file_train, train_set_traj, append=True)
else:
    # Write a new file (overwrite if it already exists)
    write(output_file_test, test_set_traj)
    write(output_file_train, train_set_traj)
