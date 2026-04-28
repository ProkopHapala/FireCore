import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import atomicUtils as au

type_map = {
    '1': 2,   # Hydrogen
    '8': 1,   # Oxygen
    '92':  -1,  # Spin-down electron
    '109': -3,  # Spin-up electron
}

mass_map = {
    1: 15.99903,  # Oxygen
    2: 1.000794,   # Hydrogen
    3: 1.0         # Electron
}


def write_lammps_data(filename, elements, positions, charges, radii, comment=None, box_size=500, element_to_atom_type=None, elec_type=3  ):
    """Write LAMMPS data file directly from eFF data"""
    atoms = []
    atom_types = set()
    
    # Create atom records with proper electron atom style format
    for (element, pos, q, radius) in zip(elements, positions, charges, radii):
        atom_type = type_map[element]
        element = int(element)
        if atom_type < 0:
            element = 0
            spin = atom_type + 2
            atom_type = elec_type
        else:
            spin = 0
        atom_types.add(atom_type)
        atoms.append((atom_type, element, spin, radius, pos, q))  # spin=0 by default
    
    # Convert to sorted list of atom types
    atom_types = sorted(atom_types)
    
    with open(filename, 'w') as f:
        # Header
        f.write("LAMMPS data file generated from eFF output\n\n")
        f.write(f"{len(elements)} atoms\n")
        f.write(f"{len(atom_types)} atom types\n\n")
        
        # Box dimensions
        f.write(f"{-box_size} {box_size} xlo xhi\n")
        f.write(f"{-box_size} {box_size} ylo yhi\n")
        f.write(f"{-box_size} {box_size} zlo zhi\n\n")
        
        # Masses
        f.write("Masses\n\n")
        for atype in atom_types:
            f.write(f"{atype} {mass_map[atype]}\n")
        f.write("\n")
        
        # Atoms section - electron style format
        f.write("Atoms # electron\n\n")
        for atom_id, (atype, element, spin, radius, pos, q) in enumerate(atoms, 1):
            x, y, z = pos
            if np.isnan(radius):
                radius = 0.0
            # Format for atom_style electron: id type q spin eradius x y z
            #f.write(f"{atom_id} {atype} {q} {spin} {radius} {x} {y} {z}\n")
            f.write( f"{atom_id:5d} {atype:5d} {element:5d} {spin:5d} {radius:10.6f} {x:10.6f} {y:10.6f} {z:10.6f} \n" )
        
        # # Velocities (optional)
        # f.write("\nVelocities\n\n")
        # for atom_id in range(1, len(elements)+1):
        #     f.write(f"{atom_id} 0 0 0\n")


def run_lmp_script(lmp, script, bPrint=False):
    for line in script.strip().split('\n'):
            line = line.strip()
            if line and not line.startswith('#'):
                if bPrint:
                    print("init cmd:", line)
                lmp.command(line)

def run_lammps_simulation(lmp, data_file, init_commands=None, post_commands=None):
    """Run LAMMPS simulation using initialization strings"""
    if init_commands:
        run_lmp_script(lmp, init_commands)
    else:
        lmp.command(f"clear")
    
    # Read data file
    print(f"read_data {data_file}")
    lmp.command(f"read_data {data_file}")
    if post_commands:
        run_lmp_script(lmp, post_commands)
    energies = {
        'total':  lmp.extract_variable("etotal", "all", 1),
        'eke':    lmp.extract_variable("v_eke", "all", 1),
        'epauli': lmp.extract_variable("v_epauli", "all", 1),
        'ecoul':  lmp.extract_variable("v_ecoul", "all", 1),
        'erres':  lmp.extract_variable("v_erres", "all", 1)
    }
    # energies = lmp.numpy.extract_compute("energies", lmp.style.compute, 1)
    # return {
    #     'total': lmp.get_thermo("etotal"),
    #     'eke': energies[0],
    #     'epauli': energies[1],
    #     'ecoul': energies[2],
    #     'erres': energies[3] if len(energies) > 3 else None
    # }

    print(energies)
    return energies

# ---- Body

imgs = au.load_xyz_movie("H2O_noe2_relax.xyz")


# for i,img in enumerate(imgs):
#     print("---------------- ", i)
#     print(img)
#     write_lammps_data(f"H2O_noe2_relax_{i}.data", img[0], img[1], img[2], img[3])

#img = imgs[0]
#write_lammps_data( "H2O_noe2_relax_0.lammps", img[0], img[1], img[2], img[3] )

init='''
units           electron
newton          on
boundary        f f f
atom_style      electron
'''

after='''
pair_style      eff/cut 100.0 ecp 1 O
pair_coeff      * *

compute         energies all pair eff/cut
variable        eke equal c_energies[1]
variable        epauli equal c_energies[2]
variable        ecoul equal c_energies[3]
variable        erres equal c_energies[4]

thermo          1
thermo_style    custom step etotal v_eke v_epauli v_ecoul v_erres
thermo_modify   format float %23.15g

group           atoms type 1 2
fix             freeze atoms setforce 0.0 0.0 0.0

min_style       cg
compute         1 all property/atom spin eradius erforce
dump            2 all custom 1 mini.lammpstrj id type q c_1[1] c_1[2] x y z fx fy fz c_1[3]
minimize        0 1e-6 2000 4000
'''

from lammps import lammps
lmp = lammps()
for i, img in enumerate(imgs[:1]):
    #data_file = f"H2O_noe2_relax_{i}.data"
    data_file = "init-.data"
    print("========== ", i ) 
    write_lammps_data("init.data", img[0], img[1], img[2], img[3], box_size=500)
    # Run LAMMPS simulation
    lmp.command(f"clear")
    energies = run_lammps_simulation(lmp, data_file, init_commands=init, post_commands=after)
    #init = None
    #exit()
    print(f"Frame {i}: Total energy = {energies['total']}")
    print(f"  E_ke: {energies['eke']}, E_pauli: {energies['epauli']}, E_coul: {energies['ecoul']}")