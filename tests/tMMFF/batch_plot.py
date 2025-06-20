#!/usr/bin/env python3
"""Batch driver for visualize_top_layer_xy.py

Searches for trajectory files in directory trees matching
    relax_*/dir_*/cons_*/PTCDA_*_total_trajectory.xyz
and runs the visualiser for each.

Usage (from tests/tMMFF directory):
    python batch_plot.py   # generates png next to each trajectory

Modify constants at the top as needed.
"""
import subprocess
import pathlib
import re
import argparse

VIS = "/home/indranil/git/FireCore/doc/py/visualize_top_layer_xy.py"  # path to visualiser
ROOT = pathlib.Path(__file__).resolve().parent  # tMMFF directory

def find_trajectories():
    """Yield Path objects pointing to trajectory xyz files."""
    pattern = "relax_*/*/cons_*/*_total_trajectory.xyz"
    yield from ROOT.glob(pattern)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry", action="store_true", help="Just list commands, do not run")
    parser.add_argument("--projections", default="xy")
    args = parser.parse_args()

    for traj in find_trajectories():
        m = re.search(r"cons_(\d+)", traj.as_posix())
        if not m:
            continue
        fixed = int(m.group(1))
        # opposite = fixed + 3  # adjust if needed
        if fixed == 24:
            opposite = 25
        elif fixed == 26:
            opposite = 29
        else:
            print(f"Unknown fixed atom {fixed} in {traj}")
            continue

        out_png = traj.with_name(traj.stem + "_plot_xyz.png")

        cmd = [
            "python",VIS,
            "--traj",str(traj),
            "--fixed",str(fixed),
            "--opposite",str(opposite),
            "--projections",args.projections,
            "--out",str(out_png),
            "--no-show",
        ]

        if args.dry:
            print("DRY:", " ".join(cmd))
        else:
            print("Running:", traj)
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()

