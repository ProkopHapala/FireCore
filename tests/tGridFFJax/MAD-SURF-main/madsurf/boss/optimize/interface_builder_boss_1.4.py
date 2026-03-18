#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 11:39:24 2022
@author: jestilj1
"""
import argparse
from ase import Atoms
from ase.build import fcc111, fcc110, fcc100, fcc211, bcc111, bcc110, bcc100
from ase.build import hcp0001, hcp10m10, diamond100, diamond111, mx2, graphene
from ase.build import add_adsorbate
from ase.constraints import FixAtoms
from ase.io import write, read, formats
from ase.io.aims import read_aims, write_aims
from ase.visualize import view
#from aimstools.preparation.aims_setup import FHIAimsSetup
## Parsing in command line arguments as parameters for the construction of the interface ##
parser = argparse.ArgumentParser()
parser.add_argument('--input_adsorbate', '-ia', type=open, required=True,
                    help="Input file specifying the adsorbate. aims .in files are the easiest to use"
                    " at the moment.")

parser.add_argument('--input_substrate', '-is', type=open, default=None,
                    help="Optional input file specifying the substrate. If not provided, provide --substrate_atoms and --crystal_structure." 
                    "input as aims in-file. NOTE: Substrate size must be provided as it is used for constraining atoms")                                 

parser.add_argument('--adsorption_height', '-z', type=float, default=3.0,
                    help="Adsorption height (in A) measured from the lowest adsorbate"
                    " atom to top of substrate layer. DEFAULT: 3.0 A")

# parser.add_argument('--adjacent_distance', '-d', type=float, default=5.0,
#                     help="Distance (in A) from adsorbate to unit cell borders. For x-,"
#                     " -y, and z-directions. DEFAULT: 5.0 A")

parser.add_argument('--substrate_atoms', '-s', default='Cu',
                    help="Substrate atoms as a string, examples: 'Cu', 'Ag', 'H2O', "
                    "'AlCO2', 'SiC'. These atoms are used to construct the substrate."
                    " DEFAULT: 'Cu'")

parser.add_argument('--crystal_structure', '-cs', default='fcc111',
                    help="This argument describes the crystal structure of the"
                    " substrate in terms of the Bravais lattice and miller indices "
                    "(for the most common low-index surfaces). Available options: 'fcc111', 'fcc110',"
                    " 'fcc100', 'fcc211', 'bcc111', 'bcc110', 'bcc100', 'hcp0001', 'hcp10m10',"
                    " 'diamond100', 'diamond111', 'mx2' and 'graphene'. DEFAULT: 'fcc111'")

parser.add_argument('--substrate_size', '-ss', type=int, nargs='+', default=[12, 12, 4],
                    help="Specifies the size of the substrate supercell by number of substrate atoms"
                    " in each direction: size = [X, Y, Z]. Input as space-separated integers:"
                    " --substrate_size X Y Z. DEFAULT: [12, 12, 4]")

parser.add_argument('--lattice_constant', '-lc', type=float,
                    help="Defines the lattice constant, a, for the unit cell."
                    " If omitted, the script will use ASE default values.")

parser.add_argument('--molindex', '-mi', type=int, default=0,
                    help="Defines the atom (by index) in the adsorbate to be positioned above the location specified by the"
                    " xy-translation argument. DEFAULT: 0"
                    " If omitted, the script will use ASE default values.")

parser.add_argument('--vacuum_length', '-v', type=float, default=30,
                    help="Defines the length of the vacuum region of the interface."
                    " Due to how ASE adds vacuum (on both sides of the slab), the total length of"
                    " the vacuum region is twice the value of this parameter (plus the substrate)."
                    " DEFAULT: 30 A, i.e. total vacuum region is (60 A - Length of substrate in z-direction).")

parser.add_argument('--orthogonal_cell', '-o', default='True',
                    help="Specifies whether the unit cell is orthogonal or not (True or False)."
                    " Caution: some common crystal structures do not support the default choice. DEFAULT: True")

parser.add_argument('--periodic_cell', '-p', default='True',
                    help="Specifies whether the unit cell is periodic or not (True or False)."
                    " Caution: some common crystal structures do not support the default choice. DEFAULT: True")

parser.add_argument('--xy_translation', '-t', type=float, nargs=2, default=[0, 0],
                    help="If the user wants to translate the adsorbate on the substrate surface,"
                    " this argument enables that. Specifies the translation by x- and y- coordinates in Angstrom. DEFAULT: (0, 0)")

parser.add_argument('--xyz_rotation', '-r', type=float, nargs=3, default=[0, 0, 0],
                    help="If the user wants to rotate the adsorbate in the unit cell, this argument enables that."
                    " Specifies the angle of rotation around the center of positions. Mind the order of rotations,"
                    " which is x -> y -> z. DEFAULT: (0, 0, 0)")

parser.add_argument('--file_format', '-ff', type=str, default='xyz',
                    help="Specifies the file format in which to save the interface. Support for most typical "
                    " computational chemistry/physics program packages, such as (FHI-)aims, Gaussian, NWChem, VASP, LAMMPS,"
                    " both input- and output-files. Confer with ase.io.formats for exact usage. DEFAULT: xyz)")

parser.add_argument('--constrained_atoms', '-ca', type=int, default=0,
                    help="Specifies the number of constrained atoms (or stricly speaking, number of layers) for"
                    " relaxation jobs in FHI-aims. Constrains the N first layers of the substrate (in z-direction)."
                    " For example: -ca 2 means that the two lowest layers of substrate atoms are fixed"
                    " during relaxation runs. Produces corresponding geometry.in files. DEFAULT = 0 ")

parser.add_argument('--magnetic_moment_atoms', '-mma', type=str, nargs='+', default=None,
                    help="Specifies which atoms to give an initial magnetic moment. Input as atom str together with the"
                    " corresponding to initial --initial_magnetic_moments parameter as explained next. Default=None")

parser.add_argument('--initial_magnetic_moments', '-im', type=int, nargs='+', default=0,
                    help="Specifies initial magnetic moments for the atoms given by --magnetic_moment_atoms. Input as integers (spin up - spin down) "
                    "corresponding to magnetic_moment_atoms: For example, to specify an initial moment of 1 for all Na atoms,"
                    " type -mma Na -im 1. Default=0")

args = parser.parse_args()


## Reading the input file and converting it into an Atoms object ##
adsorbate = read_aims(args.input_adsorbate)
#view(adsorbate)

#adsorbate = read(args.input_adsorbate) # use this if another file type is desired

## Voxels for the adsorbate defined by distance between the furthest away pair of atoms on each axis ##
x_vox = max(adsorbate.positions[:, 0]) - \
  min(adsorbate.positions[:, 0])
y_vox = max(adsorbate.positions[:, 1]) - \
    min(adsorbate.positions[:, 1])
z_vox = max(adsorbate.positions[:, 2]) - \
    min(adsorbate.positions[:, 2]) 


## Rotating the adsorbate so that most of its surface is (more or less) parallel to the substrate surface ##
if x_vox < z_vox:
    adsorbate.rotate(90, 'y', rotate_cell=True, center='COP')
elif y_vox < z_vox:
    adsorbate.rotate(90, 'x', rotate_cell=True, center='COP')


## Defining the substrate and adding the adsorbate to produce the interface ##
if args.input_substrate == None:
        substrate = eval(args.crystal_structure)(args.substrate_atoms, size=args.substrate_size,
                                                 a=args.lattice_constant, orthogonal=eval(
                                                     args.orthogonal_cell),
                                                 vacuum=args.vacuum_length, periodic=eval(args.periodic_cell))
else:
        substrate = read_aims(args.input_substrate) 

 
## Rotations taken as command line arguments ##
adsorbate.rotate(args.xyz_rotation[0], 'x', center='COP')
adsorbate.rotate(args.xyz_rotation[1], 'y', center='COP')
adsorbate.rotate(args.xyz_rotation[2], 'z', center='COP')


## Moving the substrate to the bottom of the unit cell and ensuring the adsorbate never moves into the surface ##
substrate.positions += (0, 0, -substrate.positions[:, 2].min())
subMax = substrate.positions[:, 2].max() 


## Constructing the interface by adding the adsorbate to the substrate. ##
## Point of rotation = point of translation = geometric center of the adsorbate ## 
# Defining the geometric center #
avg_pos_x = sum(adsorbate.positions[:, 0])/len(adsorbate.positions) 
avg_pos_y = sum(adsorbate.positions[:, 1])/len(adsorbate.positions)
avg_pos_z = sum(adsorbate.positions[:, 2])/len(adsorbate.positions)


# finding the distance (on each axis) between geometric center and indexed atom (mol_index) applied as an offset for add_adsorbate #
d_avg_index_corr_x = (adsorbate[args.molindex].position[0] - avg_pos_x)
d_avg_index_corr_y = (adsorbate[args.molindex].position[1] - avg_pos_y) 
d_avg_index_corr_z = (adsorbate[args.molindex].position[2] - avg_pos_z)   

add_adsorbate(substrate, adsorbate, d_avg_index_corr_z + args.adsorption_height, 
              position=(args.xy_translation[0] + d_avg_index_corr_x, args.xy_translation[1] + d_avg_index_corr_y), mol_index=(args.molindex))


## Ensuring the adsorbate does not rotate into the surface ##
AdSurf = adsorbate.positions[:, 2].min() - subMax 
aimsFailThreshold = 0.50
if AdSurf < 1.0 and AdSurf >= aimsFailThreshold:
    cost = 100
elif AdSurf <= aimsFailThreshold:
    cost = 1000
    E = 10000
    with open('energy.out', 'w') as f:
        f.write('{:15.7f}'.format(E)+"\n")
        f.write("Clash!")
        f.close()
else:
    cost = 0


## Constrain atoms ##
c_index = list(range(0,(args.substrate_size[0]*args.substrate_size[1]*args.constrained_atoms)))
c = FixAtoms(indices=(c_index))
substrate.set_constraint(c)


## Set magnetic moments (spin) ##
if args.magnetic_moment_atoms != None:
    magatoms = args.magnetic_moment_atoms
    defmom = args.initial_magnetic_moments
    magmom_list = [0 for x in range(len(substrate))]
    for i in range(0, len(substrate)):
        for j in range(0, len(magatoms)):
            if substrate[i].symbol in magatoms[j]:
                magmom_list[i] = defmom[j]
    substrate.set_initial_magnetic_moments(magmom_list)

## Set atomic charges ##
# Not implemented yet. Might not be feasible for periodic systems


## Inspect the interface visually ##
#view(adsorbate)
#view(substrate)


## Saving the file if the interface is satisfactory ##
substrate.write("geometry.in",  str(args.file_format)) # Change file extension/name as needed
#print("File saved successfully with", args.file_format, "formatting")