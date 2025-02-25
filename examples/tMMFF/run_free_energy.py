import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

import argparse

# Create parser for command-line arguments
parser = argparse.ArgumentParser(description='Run simulation with specified parameters')
parser.add_argument('--nMDsteps', type=int, help='Number of MD steps')
parser.add_argument('--nEQsteps', type=int, help='Number of equilibration steps')
parser.add_argument('--t_damp', type=float, help='Damping time')
parser.add_argument('--T', type=float, help='Temperature')
parser.add_argument('--dt', type=float, help='Time step')
parser.add_argument('--nbSteps', type=int, help='Number of steps')
parser.add_argument('--N', type=int, help='System size')
parser.add_argument('--K', type=float, help='Spring constant')
parser.add_argument('--lamda1', type=float, help='Lambda 1')
parser.add_argument('--lamda2', type=float, help='Lambda 2')

# Parse arguments
args = parser.parse_args()

# Assign values from command line arguments
MY_N = args.N
MY_K = args.K 
MY_nbStep = args.nbSteps
MY_nMDsteps = args.nMDsteps
MY_nEQsteps = args.nEQsteps
MY_t_damp = args.t_damp
MY_T = args.T
MY_dt = args.dt
MY_lamda1 = args.lamda1
MY_lamda2 = args.lamda2

xyz_name = "data/nucleobasis_AT"
constr_name = "data/nucleobasis_AT.cons"

mmff.setVerbosity( verbosity=5, idebug=1 )

mmff.init( xyz_name=xyz_name, constr_name=constr_name ,bMMFF=True)
colectiveVariable = np.array([0,1,2], dtype=np.int32)
E = mmff.compute_Free_energy(MY_lamda1, MY_lamda2, colectiveVariable, nbStep=MY_nbStep, nMDsteps=MY_nMDsteps, nEQsteps=MY_nEQsteps, t_damp=MY_t_damp, T=MY_T, dt=MY_dt)
print("E=", E)