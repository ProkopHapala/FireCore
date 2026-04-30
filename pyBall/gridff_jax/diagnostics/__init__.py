"""MLIP diagnostic and benchmarking tools for GridFF.

Provides rigorous assessment of any MLIP (MAD-SURF, MACE, etc.) and
comparison with DFT or GridFF, including:

- Energy surface smoothness analysis
- Force-energy derivative consistency checks
- Lateral corrugation analysis
- Multi-method comparison framework
- Automated report generation

All tools are substrate- and molecule-agnostic. JAX is used where
beneficial for batch finite-difference computation.
"""

from .comparison import MethodResult, compare_methods, parity_statistics
from .corrugation import corrugation_analysis
from .force_consistency import force_energy_consistency
from .report import generate_diagnostic_report
from .smoothness import smoothness_analysis

__all__ = [
    "MethodResult",
    "compare_methods",
    "parity_statistics",
    "corrugation_analysis",
    "force_energy_consistency",
    "generate_diagnostic_report",
    "smoothness_analysis",
]
