#!/usr/bin/env python3
"""Test linear hybrid potential fitting (B-spline + compact radial atomic basis) on NaCl 1x1 L3.

Fits against erf-damped Ewald2D Coulomb + Morse reference with Boltzmann-weighted
least squares. Generates diagonal XZ cuts and vertical z-scan plots.
See: doc/Topics/FastCollisionSplitNonbond/hybrid_fitting_design_gpt56sol.md
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from pyBall.OCL.Surface_utils import test_hybrid_potential_nacl

if __name__ == '__main__':
    results = test_hybrid_potential_nacl(
        save_dir='results_hybrid',
        Ng=80,
        R_cut_hc=4.0,
        spline_step=1.0,
        sample_step=0.5,
        force_weight=0.25,
        T_eff=300.,
        sigma=2.0,
    )
    print("\nDone. Results saved to results_hybrid/")
