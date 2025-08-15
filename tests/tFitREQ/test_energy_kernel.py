#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys

sys.path.append("../../")
from pyBall.OCL.NonBondFitting import run_energy_imshow

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

run_energy_imshow(
    model_name=args.model,
    xyz_file=args.xyz,
    atom_types_file=args.atypes,
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