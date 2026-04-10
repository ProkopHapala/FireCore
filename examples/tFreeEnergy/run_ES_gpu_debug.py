"""
run_ES_gpu_debug.py  –  Debug test for entropic-spring thermodynamic integration on GPU.

Runs entropic_spring_TI_gpu_debug() over a lambda grid and compares the result
against the analytical Gaussian-chain free energy.  Produces one data file
and one PNG plot.

Usage (from examples/tFreeEnergy/):
    python run_ES_gpu_debug.py [--nSys N] [--nbStep N] [--nMDsteps N] [--nEQsteps N] [-T K]
"""

import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../../")
from pyBall import MMFF_multi as mmff


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_constraint_indices(path: Path) -> np.ndarray:
    """Read bond-constraint indices from a simple text file (one index per line, # comments)."""
    idxs = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if (not line) or line.startswith("#"):
                continue
            idxs.append(int(line.split()[0]))
    if not idxs:
        raise ValueError(f"No bond-constraint indices found in {path}")
    return np.array(idxs, dtype=np.int32)


def infer_n_segments(xyz_path: Path) -> int:
    """Infer number of chain segments from filename, e.g. entropic_spring_20.xyz -> 19."""
    m = re.search(r"entropic_spring_(\d+)", xyz_path.name)
    if not m:
        raise ValueError(f"Cannot infer chain size from filename: {xyz_path.name}")
    return int(m.group(1)) - 1


def gaussian_chain_reference(lam: np.ndarray, lam0: float, T: float,
                              n_segments: int, bond_length: float = 1.198,
                              bond_k: float = 80.0) -> np.ndarray:
    """
    Analytical free energy of a Gaussian chain with harmonic bonds.

    Uses the effective spring constant
        k_eff = 3 k_B T / (N * <l^2>)
    where <l^2> = l0^2 + k_B T / k_bond (thermal correction from bond stiffness).

    Returns  F(lam) - F(lam0)  in eV.
    """
    k_b    = 8.617333262145e-5          # eV / K
    k_eff  = 3.0 * k_b * T / (n_segments * bond_length**2)
    return 0.5 * k_eff * (lam**2 - lam0**2)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="GPU debug TI for entropic spring vs analytical Gaussian-chain reference."
    )
    parser.add_argument("--nSys",        type=int,   default=1, help="Number of parallel GPU replicas")
    parser.add_argument("--xyz_name",    default="data/entropic_spring_20.xyz", help="Path to XYZ file")
    parser.add_argument("--constr_name", default="data/entropic_spring_20.cons", help="Path to .cons file")
    parser.add_argument("--constraints", default="constraints_ES.txt", help="Text file listing bond-constraint indices")
    parser.add_argument("--lamda1",      type=float, default=1.0, help="Start of lambda range [Å]")
    parser.add_argument("--lamda2",      type=float, default=11.0, help="End of lambda range [Å]")
    parser.add_argument("--nbStep",      type=int,   default=25, help="Number of lambda windows")
    parser.add_argument("--nMDsteps",    type=int,   default=10000, help="Production MD steps per lambda window")
    parser.add_argument("--nEQsteps",    type=int,   default=5000, help="Equilibration steps per lambda window")
    parser.add_argument("--t_damp",      type=float, default=100.0, help="Langevin damping time [steps]")
    parser.add_argument("-T", "--temperature", type=float, default=300.0, help="Temperature [K]")
    parser.add_argument("--dt",          type=float, default=0.05, help="MD time step [fs]")
    args = parser.parse_args()

    xyz_name    = Path(args.xyz_name)
    constr_name = Path(args.constr_name)
    dc          = load_constraint_indices(Path(args.constraints))
    n_segments  = infer_n_segments(xyz_name)
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)

    # Shared stem for output files, e.g. "entropic_spring_20_TI_gpu_debug_T300K"
    out_stem = f"{xyz_name.stem}_TI_gpu_debug_T{args.temperature:.0f}K"

    # --- Init GPU system ---
    mmff.setVerbosity(verbosity=0, idebug=0)
    mmff.init(
        nSys_=args.nSys,
        xyz_name=str(xyz_name),
        constr_name=str(constr_name),
        sElementTypes="../../cpp/common_resources/ElementTypes.dat",
        sAtomTypes="../../cpp/common_resources/AtomTypes.dat",
        sBondTypes="../../cpp/common_resources/BondTypes.dat",
        sAngleTypes="../../cpp/common_resources/AngleTypes.dat",
        bMMFF=True, bEpairs=False,
        nPBC=(0, 0, 0), gridnPBC=(0, 0, 0), GridFF=0,
        T=args.temperature, gamma=1.0 / (args.t_damp * args.dt),
        nExplore=0, nRelax=0,
    )

    # --- Run GPU debug TI ---
    mmff.entropic_spring_TI_gpu_debug(
        args.lamda1, args.lamda2, dc=dc,
        nbStep=args.nbStep, nMDsteps=args.nMDsteps, nEQsteps=args.nEQsteps,
        tdamp=args.t_damp, T=args.temperature, dt=args.dt,
    )

    # --- Load raw TI output written by C++ entropic_spring_TI_debug() ---
    raw_data = np.loadtxt(results_dir / "TI_plot_ES_GPU_DEBUG.dat", skiprows=1)
    if raw_data.ndim == 1:
        raw_data = raw_data[None, :]
    lam   = raw_data[:, 0]
    ti    = raw_data[:, 1]
    sigma = raw_data[:, 2] if raw_data.shape[1] > 2 else np.zeros_like(lam)

    # --- Analytical Gaussian-chain reference ---
    ref = gaussian_chain_reference(lam, args.lamda1, args.temperature, n_segments)

    # --- Save one data file ---
    out_dat = results_dir / f"{out_stem}.dat"
    header  = "# lambda[A]  TI_F[eV]  TI_sigma[eV]  Ref_F[eV]  TI-Ref[eV]"
    np.savetxt(out_dat,
               np.column_stack([lam, ti, sigma, ref, ti - ref]),
               header=header, comments="", fmt="%.8f")

    # --- Save one PNG ---
    out_png = results_dir / f"{out_stem}.png"
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(lam, ti,  label="GPU TI (debug)")
    ax.fill_between(lam, ti - sigma, ti + sigma, alpha=0.25, label="±1σ")
    ax.plot(lam, ref, "--", label="Gaussian chain (analytical)")
    ax.set_xlabel("λ [Å]")
    ax.set_ylabel("ΔF [eV]")
    ax.set_title(f"Entropic spring TI debug — {xyz_name.stem},  T = {args.temperature:.0f} K")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

    print(f"ΔF (GPU TI, last window) : {ti[-1]:.6f} eV")
    print(f"ΔF (analytical reference): {ref[-1]:.6f} eV")
    print(f"Difference               : {ti[-1] - ref[-1]:.6f} eV")
    print(f"Data : {out_dat}")
    print(f"Plot : {out_png}")


if __name__ == "__main__":
    main()

