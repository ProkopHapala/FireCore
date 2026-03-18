import os
import shutil
import subprocess
import numpy as np

import numpy as np
import matplotlib.pyplot as plt

from dscribe.descriptors import SOAP
from dscribe.kernels import REMatchKernel

from ase import Atoms
from ase.io import read, write
from ase.io.trajectory import Trajectory

from sklearn.preprocessing import normalize, StandardScaler
from sklearn.decomposition import PCA

def extract_mode_number(filename):
    # Extract mode number from the filename, assuming it follows the format "modified_traj_mode{mode_number}.traj"
    return int(filename.split('_')[-1].split('.')[0].replace('mode', ''))

def read_traj_and_freq(traj_filename, vib_xyz_filename):
    # Read trajectories from the .traj file
    trajs = read(traj_filename, index=':10')       #Trajectory(traj_filename)

    # Extract mode number from the traj_filename
    mode_numbers = [extract_mode_number(traj_filename)]

    # Read frequencies from the vib.xyz file for specified mode numbers
    freq_dict = {}
    with open(vib_xyz_filename, 'r') as file:
        lines = file.readlines()
        for mode_number in mode_numbers:
            frequencies = []
            for line in lines:
                if f'Mode #{mode_number}' in line:
                    freq = float(line.split('=')[1].split()[0])
                    frequencies.append(freq)
            freq_dict[mode_number] = np.array(frequencies)

    return trajs #, freq_dict


def farthest_point_sampling(similarities, n_points=None):
    N = similarities.shape[0]
    farthest_indices = np.zeros(N, int)
    ds = similarities[0, :]
    for i in range(1, N):
        idx = np.argmin(ds)
        farthest_indices[i] = idx
        ds = np.maximum(ds, similarities[idx, :])
    #print(distances, farthest_indices)
    if n_points is None:
        return farthest_indices
    else:
        return farthest_indices[:n_points]
    

def one_sampling(base_path, saving_path):
    # create the output directory if it doesn't exist
    os.makedirs(saving_path, exist_ok=True)
    opt_atoms = read(f'{base_path}/geometry.xyz')
    n_atoms = len(opt_atoms)
    trajs = [opt_atoms]
    indices = range(6, 3*n_atoms) # Since the indexing starts from 6
    n_points = 1 + 3*n_atoms
    plot_indices = [0]

    for i, idx in enumerate(indices):
        traj_file_path = f'{base_path}/modified_traj_mode{idx}.traj'

        # Check if the trajectory file exists
        if os.path.exists(traj_file_path):
            # Read the trajectory and frequencies
            trajs += read_traj_and_freq(traj_file_path, f'{base_path}/vib.xyz')
            # Update plot_indices
            plot_indices += [i + 1] * 10 # Since the number of structures which are read from traj files are 10
        else:
            print(f"File {traj_file_path} does not exist, skipping index {idx}")

    
    print(f'number of structures: {len(trajs)}')
    print(f'number of atoms: {n_atoms}')
    np.save(f'{saving_path}/mode_plot_indices.npy', plot_indices)

    # First we will have to create the features for atomic environments. Lets
    # use SOAP.
    local_desc = SOAP(species=["H", "C", "N", "O"], r_cut=4.0, n_max=4, l_max=5, sigma=0.2)
    global_desc = SOAP(species=["H", "C", "N", "O"], r_cut=4.0, n_max=4, l_max=5, sigma=0.2, average='inner')

    print(f'calculating {local_desc.get_number_of_features()} local features ')
    all_local_features = local_desc.create(trajs)
    print(f'calculating {global_desc.get_number_of_features()} global features ')
    all_global_features = global_desc.create(trajs)

    # Before passing the features we normalize them. Depending on the metric, the
    # REMatch kernel can become numerically unstable if some kind of normalization
    # is not done.
    print('normalizing local features')
    all_normalized_features = [normalize(features) for features in all_local_features]
    

    # Calculates the similarity with the REMatch kernel and a linear metric. The
    # result will be a full similarity matrix.
    print('REMatchKernel')
    re = REMatchKernel(metric="linear", alpha=1e-2, threshold=1e-7, normalize_kernel=True)
    re_kernel = re.create(all_normalized_features)
    #distances = StandardScaler().fit_transform(re_kernel)
    all_indices = farthest_point_sampling(re_kernel)
    print(all_indices)
    unique_all_indices = np.unique(all_indices)
    print(unique_all_indices)

    selected_indices = all_indices[:n_points]
    selected_trajs = [trajs[i] for i in selected_indices]
    print(len(all_indices), len(selected_indices))
    np.save(f'{saving_path}/farthest_all_indices.npy', all_indices)
    np.save(f'{saving_path}/farthest_selected_indices.npy', selected_indices)

    all_trajs = [trajs[i] for i in unique_all_indices] 
    write(f'{saving_path}/farthest_selected_all_indices.traj', all_trajs)
    write(f'{saving_path}/farthest_selected_non_opt_with_velocities.traj', selected_trajs)

    write(f'{saving_path}/farthest_selected_all_indices.extxyz', all_trajs, format='extxyz')
    write(f'{saving_path}/farthest_selected_indices.extxyz', selected_trajs,format='extxyz')
    
    # PCA analysis
    print('PCA analysis')
    model = PCA(n_components=2, random_state=1000).fit(all_global_features)
    reduced_global_features = model.transform(all_global_features)
    np.save(f'{saving_path}/global_features_PCA.npy', reduced_global_features)

one_sampling("./mode_analysis/ase_traj/" , './result/')


## code for launching the dft calculation

def copy_file_and_submit_job(filename):
    source_dir = './result'
    destination_dir = './DFT_cal/all_indices'
    job_script = 'run_new_modules.sh'

    # Construct full file paths
    source_path = os.path.join(source_dir, filename)
    destination_path = os.path.join(destination_dir, filename)

    try:
        # Copy the file
        shutil.copy(source_path, destination_path)
        print(f"Copied {filename} to {destination_path}")

        # Change directory
        os.chdir(destination_dir)
        print(f"Changed directory to {os.getcwd()}")

        # Submit the job
        result = subprocess.run(['sbatch', job_script])

    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
#copy_file_and_submit_job('farthest_selected_all_indices.extxyz')
