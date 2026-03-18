
from ase.io import write, read
from ase.constraints import FixAtoms
from ase.io.trajectory import Trajectory
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation,
)
from ase.md import MDLogger
from mace.calculators import MACECalculator
from ase import units
from sys import argv
import os

# =======================
# PARAMETERS
# =======================
timestep = 1.0 * units.fs

interval = 1 # printing interval
nsteps_anneal = 10000 # number of steps taken for the annealing
nsteps_quench = 200 # number of steps per temperature increment (temperature_anneal-temperature_quench/|dT|)

temperature_anneal = 400  # K
temperature_quench = 10   # K
dT = -10                  # temperature step

friction = 0.01 / units.fs  # weak coupling

z_threshold = 9.0 # z-threshold for fixing atoms
fixed_atom_indices = [] # if any other atoms need fixing, define their indices

# =======================
# INPUT
# =======================
atoms = read(argv[1])
save_file = os.path.splitext(argv[1])[0]

# =======================
# CONSTRAINTS (by z or idx)
# =======================
fixed_indices = [
    i for i, atom in enumerate(atoms)
    if atom.position[2] < z_threshold
]

atoms.set_constraint([
    FixAtoms(indices=fixed_indices),
    FixAtoms(indices=fixed_atom_indices),
])

# =======================
# CALCULATOR (with model)
# =======================
atoms.calc = MACECalculator(
    model_paths="./MAD-SURF.model",
    device="cuda",
    default_dtype="float64",
)

# =======================
# INITIAL VELOCITIES
# =======================
MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_anneal)
Stationary(atoms)
ZeroRotation(atoms)

# =======================
# LANGEVIN DYNAMICS (NVT)
# =======================
dyn = Langevin(
    atoms,
    timestep=timestep,
    temperature_K=temperature_anneal,
    friction=friction,
)

def write_frame():
    # === Logging and saving ===
    dyn.atoms.write(f"{save_file}_{temperature_anneal}K_anneal_OH_network.xyz", append=True, format="extxyz")

dyn.attach(MDLogger(dyn, atoms, f"{save_file}_{temperature_anneal}K_anneal_OH_network.log", header=True, stress=False,
           peratom=True, mode="a"), interval=interval)
dyn.attach(write_frame, interval=interval)

# =======================
# ANNEAL
# =======================
print("Starting annealing...")
dyn.run(nsteps_anneal)

# =======================
# QUENCH
# =======================
print("Starting quench...")
for T in range(int(atoms.get_temperature()), temperature_quench - 1, dT):
    print(f"Quenching to {T} K")
    dyn.set_temperature(temperature_K=T)
    dyn.run(nsteps_quench)

print("Done.")
