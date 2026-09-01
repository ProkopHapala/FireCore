"""XRD simulation utilities: pair-distance histogram + Debye transform on GPU."""
from .form_factors import get_form_factor_table, cromer_mann
from .debye_histogram import (
    XRDDebye, load_xyz_atoms,
    compute_sigma_from_sparse_blocks, compute_sigma_exact,
)
