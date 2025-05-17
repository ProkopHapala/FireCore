import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

sys.path.append("../../")
from pyBall import eFF as eff

from pyBall import eFF_terms as pyeff
const_bohr = 0.5291772105638411
eff.setVerbosity(0) # If 1: it will  write more stuff to console; If 0: it will wirte less stuff to console

def save_fgo(filename, atoms, electrons, closed_shell=False):
    """
    Save molecular system in .fgo format
    
    Parameters:
    -----------
    filename : str
        Path to save the .fgo file
    atoms : list of dicts
        Each dict contains:
            'pos': [x,y,z] - position
            'Z': int - atomic number (nuclear charge)
            'sQ': float - charge distribution width
            'sP': float - Pauli repulsion radius
            'cP': float - number of core electrons
    electrons : list of dicts
        Each dict contains:
            'pos': [x,y,z] - position
            'size': float - electron radius
            'coef': float - expansion coefficient (usually 1.0)
            'spin': int - spin value (1 or -1)
    closed_shell : bool
        Whether system is closed shell
    """
    
    # Create directory if it doesn't exist
    os.makedirs('data', exist_ok=True)
    
    with open(f'data/{filename}', 'w') as f:
        # Write header
        n_atoms = len(atoms)
        n_electrons = len(electrons)
        per_orb = 1
        f.write(f"{n_atoms} {n_electrons} {per_orb} {1 if closed_shell else 0}\n")
        
        # Write atoms
        for atom in atoms:
            pos = atom['pos']
            f.write(f"{pos[0]:8.6f} {pos[1]:8.6f} {pos[2]:8.6f}   ")
            f.write(f"{-atom['Z']:8.6f} {atom['sQ']:8.6f} {atom['sP']:8.6f} {atom['cP']:8.6f}\n")
        
        # Write electrons
        for e in electrons:
            pos = e['pos']
            f.write(f"{pos[0]:8.6f} {pos[1]:8.6f} {pos[2]:8.6f} ")
            f.write(f"{e['size']:8.6f} {e['coef']:8.6f} {e['spin']:3d}\n")

# Example usage - creating H2 molecule
def createH2(r): # r is half the distance between the two H atoms
    # Define H2 molecule
    h2_atoms = [
        {'pos': [r, 0.0, 0.0], 'Z': 1, 'sQ': 0.2, 'sP': 1.0, 'cP': 0.0},
        {'pos': [-r, 0.0, 0.0], 'Z': 1, 'sQ': 0.2, 'sP': 1.0, 'cP': 0.0}
    ]

    h2_electrons = [
        {'pos': [1.1, -0.5, 0.0], 'size': 0.5, 'coef': 1.0, 'spin': -1},
        {'pos': [1.0, 0.5, 0.0], 'size': 0.5, 'coef': 1.0, 'spin': 1},
        {'pos': [-1.1, -0.5, 0.0], 'size': 0.5, 'coef': 1.0, 'spin': -1},
        {'pos': [-1.0, 0.5, 0.0], 'size': 0.5, 'coef': 1.0, 'spin': 1}
    ]

    save_fgo('H2.fgo', h2_atoms, h2_electrons)

# def createCH4(r): # r is half the distance between the two H atoms
#     # Define H2 molecule
#     h2_atoms = [
#         {'pos': [r, 0.0, 0.0], 'Z': 1, 'sQ': 0.2, 'sP': 1.0, 'cP': 0.0},
#         {'pos': [-r, 0.0, 0.0], 'Z': 1, 'sQ': 0.2, 'sP': 1.0, 'cP': 0.0}
#     ]

#     h2_electrons = [
#         {'pos': [1.1, -0.5, 0.0], 'size': 0.5, 'coef': 1.0, 'spin': -1},
#         {'pos': [1.0, 0.5, 0.0], 'size': 0.5, 'coef': 1.0, 'spin': 1},
#         {'pos': [-1.1, -0.5, 0.0], 'size': 0.5, 'coef': 1.0, 'spin': -1},
#         {'pos': [-1.0, 0.5, 0.0], 'size': 0.5, 'coef': 1.0, 'spin': 1}
#     ]


# outE=eff.relax_mol("H2", dt=0.03, outE=True)  

radiusMin = 0.1
radiusMax = 3.5
radiusSteps = 100
radiusSpace = np.linspace(radiusMin, radiusMax, radiusSteps)
outESpace = []

i = 0
for radius in tqdm(radiusSpace, desc="Computing Progress of energy to radius function"):
    input()
    createH2(radius)
    print("CreatedH2")
    outE=eff.relax_mol("H2", dt=0.03, outE=True)  
    print("relaxed")
    outESpace.append(outE)
    i += 1
    print(f"Progress: {i}/{radiusSteps}")

Etot_value = [arr[0] for arr in outESpace]
Ek_value = [arr[1] for arr in outESpace]
Eee_value = [arr[2] for arr in outESpace]
EeePaul_value = [arr[3] for arr in outESpace]
Eae_value = [arr[4] for arr in outESpace]
EaePaul_value = [arr[5] for arr in outESpace]
Eaa_value = [arr[6] for arr in outESpace]
eV_value = [arr[7] for arr in outESpace]

plt.plot(radiusSpace, Etot_value,"o--", label='Etot', color='blue')
# plt.plot(radiusSpace, Ek_value, label='Ek', color='red')
# plt.plot(radiusSpace, Eee_value, label='Eee', color='green')
# plt.plot(radiusSpace, EeePaul_value, label='EeePaul', color='orange')
# plt.plot(radiusSpace, Eae_value, label='Eae', color='purple')
# plt.plot(radiusSpace, EaePaul_value, label='EaePaul', color='brown')
# plt.plot(radiusSpace, Eaa_value, label='Eaa', color='pink')
# plt.plot(radiusSpace, eV_value, label='eV', color='gray')

plt.xlabel('Radius [A]')
plt.ylabel('Energy [eV]')
plt.legend()
plt.grid()
plt.show()