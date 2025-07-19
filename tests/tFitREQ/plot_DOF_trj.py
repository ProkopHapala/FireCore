#!/usr/bin/python3

import enum
import numpy as np
import matplotlib.pyplot as plt

def extract_dof_names(output_file):
    """
    Extract DOF names from the output file.
    Returns a list of DOF names in order.
    """
    dof_names = []
    reading_dofs = False
    with open(output_file, 'r') as f:
        for line in f:
            if 'printDOFsToTypes()' in line:
                dof_names = []
                reading_dofs = True
                continue
            if reading_dofs:
                if '->' in line:
                    # Extract the DOF name after the '=' sign
                    name = line.split('|')[1][1:]
                    #print (name, "    line: ", line, )
                    dof_names.append( name.strip() )
                elif line.strip() == '':
                    reading_dofs = False
                    break
    #print(dof_names)
    for i,name in enumerate(dof_names):
        print(i,name)
    return dof_names

def read_DOF_trj(output_file):
    # Initialize lists to store data
    steps      = []
    dofs_data  = []
    fdofs_data = []
    # Read the output file
    with open(output_file, 'r') as f:
        for line in f:
            if '_REG:  step:' not in line: continue
            parts = line.strip().split()
            step = int(parts[parts.index('step:') + 1])
            if ' DOFs:' in line:  # Note the space before DOFs to avoid matching fDOFs
                # Extract DOFs values
                dof_idx = parts.index('DOFs:') + 1
                dofs = [float(x) for x in parts[dof_idx:]]
                if len(steps) == 0 or steps[-1] != step:
                    steps.append(step)
                    dofs_data.append(dofs)
            elif 'fDOFs:' in line:
                # Extract fDOFs values
                fdof_idx = parts.index('fDOFs:') + 1
                fdofs = [float(x) for x in parts[fdof_idx:]]
                fdofs_data.append(fdofs)
    # Convert to numpy arrays for easier plotting
    if( len(steps) == 0 ):
        print("ERROR in read_DOF_trj() len(steps) == 0")
        exit()
    steps      = np.array(steps)
    dofs_data  = np.array(dofs_data)
    fdofs_data = np.array(fdofs_data)
    return steps, dofs_data, fdofs_data


def unselect_dofs(dof_names, unselect_strings : set ):
    '''
    if name in dof_selection is containing sub-string from unselect_strings then it will be unselected, remaining will be copied to dof_selection
    '''
    dof_selection = []
    for name in dof_names:
        if not any(substring in name for substring in unselect_strings):
            dof_selection.append(name)
    return dof_selection

def select_only_dofs(dof_names, select_strings : set ):
    dof_selection = []
    for name in dof_names:
        if any(substring in name for substring in select_strings):
            dof_selection.append(name)
    return dof_selection
    
def plot_dofs_fdofs(output_file="OUT", figsize=(12, 8), dof_selection=None, dof_names=None, unselect_strings=None, select_strings=None, bPrint=True):
    """
    Plot DOFs and fDOFs from the output file with better control and nicer appearance.
    
    Args:
        output_file (str): Path to the output file containing DOFs and fDOFs data
        figsize (tuple): Figure size in inches (width, height)
    """

    steps, dofs_data, fdofs_data = read_DOF_trj(output_file)

    # Get DOF names
    if dof_names is None:
        dof_names = extract_dof_names(output_file)
    # if dof_names contains .H should be full line '-' if not should doted ':'
    dof_ls = [ '--' if 'E'==name[0] else '-' for name in dof_names ]

    if select_strings is not None:
        dof_selection = set(select_only_dofs(dof_names, select_strings))
    else:
        dof_selection = set(dof_names)

    if unselect_strings is not None:
        if dof_selection is None:
            dof_selection = set(unselect_dofs(dof_names,     unselect_strings))
        else:
            dof_selection = set(unselect_dofs(dof_selection, unselect_strings))
        
    # Filter DOFs if dof_selection is provided
    if dof_selection is not None:
        selected_indices = [i for i, name in enumerate(dof_names) if name in dof_selection]
        dofs_data  = dofs_data[:, selected_indices]
        fdofs_data = fdofs_data[:, selected_indices]
        dof_names  = [f"{i}:{dof_names[i]}" for i in selected_indices]
        dof_ls = [dof_ls[i] for i in selected_indices]
    else:
        dof_names = [f"{i}:{name}" for i, name in enumerate(dof_names)]

    if bPrint:
        for i,name in enumerate(dof_names):
            print(i,name)
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
    
    # Plot DOFs
    for i in range(dofs_data.shape[1]):
        ax1.plot(steps, dofs_data[:, i], ls=dof_ls[i], label=dof_names[i])
    
    ax1.set_xlabel('Step')
    ax1.set_ylabel('DOFs')
    ax1.set_title('Degrees of Freedom (DOFs) vs Step')
    ax1.grid(True)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Plot fDOFs
    for i in range(fdofs_data.shape[1]):
        ax2.plot(steps, fdofs_data[:, i], ls=dof_ls[i], label=dof_names[i])
    
    ax2.set_xlabel('Step')
    ax2.set_ylabel('fDOFs')
    ax2.set_title('Variational Forces (fDOFs) vs Step')
    ax2.grid(True)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    '''
    steps, dofs_data, fdofs_data = read_DOF_trj(output_file)

    # Get DOF names
    dof_names = extract_dof_names(output_file)
    # if dof_names contains .H should be full line '-' if not should doted ':'
    dof_ls = [ '--' if 'E'==name[0] else '-' for name in dof_names ]
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
    
    # Plot DOFs
    for i in range(dofs_data.shape[1]):
        label = dof_names[i] if i < len(dof_names) else f'DOF {i+1}'
        ax1.plot(steps, dofs_data[:, i], ls=dof_ls[i], label=label)
    
    ax1.set_xlabel('Step')
    ax1.set_ylabel('DOFs')
    ax1.set_title('Degrees of Freedom (DOFs) vs Step')
    ax1.grid(True)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Plot fDOFs
    for i in range(fdofs_data.shape[1]):
        label = dof_names[i] if i < len(dof_names) else f'fDOF {i+1}'
        ax2.plot(steps, fdofs_data[:, i], ls=dof_ls[i], label=label)
    
    ax2.set_xlabel('Step')
    ax2.set_ylabel('fDOFs')
    ax2.set_title('Variational Forces (fDOFs) vs Step')
    ax2.grid(True)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ''' 
    # Adjust layout to prevent overlapping
    plt.tight_layout()
    return fig, (ax1, ax2)

if __name__ == "__main__":
    plot_dofs_fdofs(output_file="OUT", figsize=(12, 8), select_strings={'E_H'} ) 
    plot_dofs_fdofs(output_file="OUT", figsize=(12, 8), select_strings={'E_O'} )
    plot_dofs_fdofs(output_file="OUT", figsize=(12, 8), select_strings={'E_N'} )
    plot_dofs_fdofs(output_file="OUT", figsize=(12, 8), unselect_strings={'E_H','E_O','E_N'} )
    plt.show()