#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall.OCL.NonBondFitting import run_energy_imshow, setup_driver

np.set_printoptions(linewidth=200)

ROOT = Path(__file__).parent.parent.parent

p = argparse.ArgumentParser(description="Evaluate energy-only OCL kernel and plot imshow vs reference")
p.add_argument('--model',  default='ENERGY_MorseQ_PAIR', help='Energy model macro to inject (e.g., ENERGY_MorseQ_PAIR, ENERGY_LJQH2_PAIR)')
p.add_argument('--xyz',    default=str(ROOT / 'tests/tFitREQ/HHalogens/porcessed/HBr-D1_HBr-A1.xyz'),  help='Path to packed XYZ scan file with Etot headers')
p.add_argument('--atypes', default=str(ROOT / 'cpp/common_resources/AtomTypes.dat'), help='Path to AtomTypes.dat with base parameters')
#p.add_argument('--kcal',           action='store_true', help='Plot energies in kcal/mol (default: eV)')
#p.add_argument('--sym',            action='store_true', help='Use symmetric color scale around 0 for ref/model panels')
p.add_argument('--emin',   default=2.0, type=float,     help='Magnitude for symmetric color scale if --sym (default 2.0)')
p.add_argument('--cmap',   default='bwr',               help='Matplotlib colormap for imshow panels (default: bwr)')
#p.add_argument('--lines',          action='store_true', help='Show additional line profiles figure comparing ref vs model')
p.add_argument('--rmax',   default=8.0, type=float,     help='Max radius shown on line-plot x-axis (default 8.0)')
p.add_argument('--no-colorbar',    action='store_true', help='Disable colorbars')
p.add_argument('--save',   default=None,                help='Path to save the figure (PNG). If not set, do not save')
p.add_argument('--no-show',       action='store_true',  help='Do not display the figure')
p.add_argument('-v', '--verbose',default=0,  action='count',  help='Increase verbosity (repeat for more)')

args = p.parse_args()


os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
os.environ['PYOPENCL_CTX'] = '0'


drv = setup_driver(model_name=args.model, atom_types_file=args.atypes, verbose=args.verbose)

data_path=ROOT / 'tests/tFitREQ/HHalogens/porcessed/'


#xyz_path=data_path/"HBr-D1_HBr-A1.xyz"; drv.load_data(xyz_path); drv.type_set=[('H' ,'H',0.85),('Br','H',-0.85)]; drv.alphaMorse = 1.8; q=0.25;drv.set_charges_for_sample_atoms([ (0,-q),(1,q), (2,-q),(3,q) ])
#xyz_path=data_path/"HCl-D1_HCl-A1.xyz"; drv.load_data(xyz_path); drv.type_set=[('H' ,'H',0.90),('Cl','H',-0.90)]; drv.alphaMorse = 1.8; q=0.2;drv.set_charges_for_sample_atoms([ (0,-q),(1,q), (2,-q),(3,q) ])
#xyz_path=data_path/"HF-D1_HF-A1.xyz";   drv.load_data(xyz_path); drv.type_set=[('H' ,'H',0.85),('F' ,'H',-0.80)]; drv.alphaMorse = 1.8; # q=0.2;drv.set_charges_for_sample_atoms([ (0,-q),(1,q), (2,-q),(3,q) ])

#xyz_path=data_path/"HBr-A1_HBr-D1.xyz";  drv.load_data(xyz_path); drv.type_set= [('H' ,'H',0.85),('Br','H',-0.85)]; drv.alphaMorse = 1.8; q=0.25; drv.set_charges_for_sample_atoms([ (0,-q),(1,q), (2,-q),(3,q) ])
#xyz_path=data_path/"HCl-A1_HCl-D1.xyz"; drv.load_data(xyz_path); drv.type_set= [('H' ,'H',0.90),('Cl','H',-0.90)]; drv.alphaMorse = 1.8; q=0.2;drv.set_charges_for_sample_atoms([ (0,-q),(1,q), (2,-q),(3,q) ])
xyz_path=data_path/"HF-A1_HF-D1.xyz";   drv.load_data(xyz_path); drv.type_set= [('H' ,'H',0.85),('F' ,'H',-0.80)]; drv.alphaMorse = 1.8; # q=0.2;drv.set_charges_for_sample_atoms([ (0,-q),(1,q), (2,-q),(3,q) ])

run_energy_imshow(
    drv=drv,
    #model_name=args.model,
    xyz_file=xyz_path,
    #atom_types_file=args.atypes,
    #kcal=args.kcal,
    #sym=args.sym,
    #Emin=args.emin,
    bColorbar=not args.no_colorbar,
    verbose=args.verbose,
    show=not args.no_show,
    save=args.save,
    cmap=args.cmap,
    #lines=args.lines,
    rmax=args.rmax,
)