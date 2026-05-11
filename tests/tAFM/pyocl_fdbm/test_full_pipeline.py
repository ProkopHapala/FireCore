#!/usr/bin/env python3
"""
Test full AFM pipeline from .xyz to AFM images using AFM_utils.py
"""

import sys, os
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.OCL import AFM_utils as afm_utils
from pyBall.DFTB.DFTBplusParser import parse_basis_hsd_ang
from pyBall import dftb_utils as du
import numpy as np

SK = 'mio-1-1'
basis = parse_basis_hsd_ang(du.WFC_HSD_PATHS[SK])
slako_prefix = du.SK_PATHS[SK]

results = afm_utils.run_afm_from_xyz(
    xyz_file='pentacene.xyz',
    output_dir='test_pipeline_output',
    basis=basis,
    slako_prefix=slako_prefix,
    co_tip_dir='old_2/debug_dftb_mio_1_1/co_tip',
    step=0.1, margin=4.0, z_extra=6.0,
    scan_range=3.0, scan_step=0.1,
    height_range=(2.8, 3.6), height_step=0.1,
    pauli_params={'A': 965.0, 'beta': 0.871},
    vdw_params={'C6_CO': 30.0},
    relax_params={'K_LAT': 0.5},
    plot_steps=True
)

print(f"\ndf shape: {results['df'].shape}")
print(f"df range: [{results['df'].min():.4f}, {results['df'].max():.4f}]")

# Plot diagnostic panel with field components
inter = results['intermediates']
heights = np.arange(2.8, 3.6, 0.1)
afm_utils.plot_diagnostic_panel(
    inter['E_pauli_field'], inter['E_ES_field'], inter['E_vdw'],
    inter['E_pauli_field'] + inter['E_ES_field'] + inter['E_vdw'],
    results['grid_spec']['origin'], 0.1, heights,
    'test_pipeline_output'
)
