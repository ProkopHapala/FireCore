"""
Fireball Fdata plotting CLI

Examples:
- Plot one wavefunction (H, shell 1) only:
    python plot_fdata_wf_2c.py --fdata Fdata --nz 1 --shell 1 --jobs wf --wf_out h_s1.png
- Plot all wavefunctions for C only:
    python plot_fdata_wf_2c.py --fdata Fdata --nz 6 --plot_all_wf --jobs wf
- Plot 2c overlap for C-C and H-C only:
    python plot_fdata_wf_2c.py --fdata Fdata --jobs c2 --c2_list 6:6,1:6 --channels 4
- Compare s/p across elements (by L panels):
    python plot_fdata_wf_2c.py --fdata Fdata --jobs wf_comp_l --elements 1,6,7,8 --deg_min 4 --deg_max 8
- Stacked per-element s+p panels:
    python plot_fdata_wf_2c.py --fdata Fdata --jobs wf_comp_elem --elements 1,6,7,8 --deg_min 4 --deg_max 8
- Fit-only summary (no plots):
    python plot_fdata_wf_2c.py --fdata Fdata --jobs wf_fit --elements 1,6,7,8 --deg_min 4 --deg_max 8
- Everything (wf + 2c + vnl + pp + wf comparisons):
    python plot_fdata_wf_2c.py --fdata Fdata --jobs wf,c2,vnl,pp,wf_comp_l,wf_comp_elem
"""

import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Utility to visualize Fdata contents:
# - Wavefunctions: one-body radial (basis/*.wf*)
# - 2-center integrals: overlap/kinetic/vna/vnl/etc. (columns are mu-nu matrix elements)
# - vnl: non-local pseudopotential projectors; each column is one projector-channel pair
#   (channels argument controls how many columns to plot for readability)
# - PP cutoff: rc_PP read from info.dat (no channels; just a radius marker)

# Add pyBall to path
sys.path.append(os.path.abspath('../../'))
from pyBall.FireballOCL.FdataParser import FdataParser
from pyBall.FireballOCL.FireballPlots import (
    plot_wavefunctions,
    plot_wavefunctions_all,
    plot_2c,
    plot_2c_list,
    plot_2c_compare,
    plot_vnl_pp,
    plot_pp_radii,
    plot_pp,
    plot_na,
    plot_wf_compare_by_L,
    plot_wf_compare_by_elem,
    wf_fit_summary,
    plot_wf_fits,
)

def main():
    ap = argparse.ArgumentParser(description="Plot wavefunctions and 2c integrals from Fdata")
    # General
    ap.add_argument("--fdata", type=str, default="Fdata", help="Fdata directory")
    ap.add_argument("--jobs", type=str, default="wf_comp_l", help="comma list of jobs: wf,c2,vnl,pp,wf_comp_l,wf_comp_elem,wf_fit,wf_fit_plot")

    # Wavefunction options
    ap.add_argument("--nz", type=int, default=1, help="nuclear charge for wavefunction")
    ap.add_argument("--shell", type=int, default=1, help="which wavefunction shell (1-based)")
    ap.add_argument("--wf_out", default="wf.png", help="output image for wavefunction")
    ap.add_argument("--plot_all_wf", action="store_true", help="plot all wavefunction shells for nz")

    # 2c integrals
    ap.add_argument("--root", default="overlap", help="2c root (overlap, kinetic, vna, vnl, ...)")
    ap.add_argument("--nz1", type=int, default=6, help="nuclear charge 1 for 2c")
    ap.add_argument("--nz2", type=int, default=6, help="nuclear charge 2 for 2c")
    ap.add_argument("--channels", type=int, default=6, help="how many 2c channels to plot")
    ap.add_argument("--c2_out", default="c2.png", help="output image for 2c integrals")
    ap.add_argument("--c2_list", default='1:1,6:6,1:6', help="comma-separated list of nz1:nz2 pairs (e.g., '1:1,6:6') to batch plot 2c")

    # vnl / PP options
    ap.add_argument("--vnl_channels", type=int, default=6, help="how many vnl channels to plot")
    ap.add_argument("--vnl_out", default="vnl.png", help="output image for vnl integrals")
    ap.add_argument("--pp_cutoff_out", default="pp_cutoff.png", help="output image for PP cutoff")
    ap.add_argument("--pp_out", default="pp.png", help="output image for one-body PP")
    ap.add_argument("--na_out", default="na.png", help="output image for neutral-atom potentials")

    # Wavefunction comparison / fitting
    ap.add_argument("--elements", type=str, default="1,6,7,8", help="comma-separated nuclear charges for comparison, e.g., '1,6,7,8'")
    ap.add_argument("--deg_min", type=int, default=4, help="min polynomial degree for wf fit")
    ap.add_argument("--deg_max", type=int, default=8, help="max polynomial degree for wf fit")

    args = ap.parse_args()

    parser = FdataParser(args.fdata)
    # this populates species_info and keeps loud if missing
    parser.parse_info()

    jobs = {j.strip() for j in args.jobs.split(",") if j.strip()}

    # Wavefunctions
    if "wf" in jobs:
        if args.plot_all_wf:
            plot_wavefunctions_all(parser, args.nz, os.path.splitext(args.wf_out)[0])
        else:
            plot_wavefunctions(parser, args.nz, args.shell, args.wf_out)

    # 2c integrals
    if "c2" in jobs:
        if args.c2_list:
            specs = []
            for pair in args.c2_list.split(","):
                if ":" in pair:
                    a, b = pair.split(":")
                    try:
                        specs.append((int(a), int(b)))
                    except ValueError:
                        print(f"[DEBUG] skipping malformed pair '{pair}'")
            if specs:
                plot_2c_list(parser, args.root, specs, args.channels, os.path.splitext(args.c2_out)[0])
        else:
            plot_2c(parser, args.root, args.nz1, args.nz2, args.channels, args.c2_out)

    # vnl / PP
    if "vnl" in jobs:
        plot_vnl_pp(parser, args.nz1, args.nz2, args.vnl_channels, args.vnl_out)
    if "pp" in jobs:
        plot_pp_radii(parser, args.nz1, args.pp_cutoff_out)
        plot_pp(parser, args.nz1, args.pp_out)
        plot_na(parser, args.nz1, args.na_out)

    elems = [int(x) for x in args.elements.split(',') if x.strip() != ""]
    degs = range(args.deg_min, args.deg_max + 1)
    # Multi-element wf comparison with fits
    if "wf_comp_l" in jobs:
        plot_wf_compare_by_L(parser, elements=elems, degs=degs)
    if "wf_comp_elem" in jobs:
        plot_wf_compare_by_elem(parser, elements=elems, degs=degs)
    if "wf_fit" in jobs:
        wf_fit_summary(parser, elements=elems, degs=degs)
    if "wf_fit_plot" in jobs:
        plot_wf_fits(parser, elements=elems, degs=degs)
    plt.show()


if __name__ == "__main__":
    main()
