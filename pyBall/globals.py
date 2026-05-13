"""
Global debugging/verbosity controls used across FireCore.

This module is intentionally tiny and dependency-free so it can be imported
from performance-sensitive code without side effects.
"""

from __future__ import annotations

import os
from typing import Any


def _get_int_env(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return default


# Print verbosity: debug_print(level, msg) will print when DEBUG_PRINT_LEVEL >= level.
# Default is 0 to avoid changing existing console output unless code explicitly
# uses debug_print().
DEBUG_PRINT_LEVEL: int = _get_int_env("AFM_DEBUG_PRINT_LEVEL", 1)

# Controls how much auxiliary data gets written to disk during AFM runs.
# Default is 2 to preserve current behavior in this codebase, where most debug
# outputs are always saved today.
DEBUG_SAVE_LEVEL: int = _get_int_env("AFM_DEBUG_SAVE_LEVEL", 2)

# Controls how many auxiliary figures (.png) are generated.
# Default is 2 to preserve current behavior.
DEBUG_PLOT_LEVEL: int = _get_int_env("AFM_DEBUG_PLOT_LEVEL", 2)


def debug_print(level: int, message: str) -> None:
    """Conditional print helper (used by new debug-only code)."""
    if DEBUG_PRINT_LEVEL >= int(level):
        print(message)


def debug_save_enabled(level: int) -> bool:
    """Return True if auxiliary array saving should happen."""
    return DEBUG_SAVE_LEVEL >= int(level)


def debug_plot_enabled(level: int) -> bool:
    """Return True if auxiliary plotting should happen."""
    return DEBUG_PLOT_LEVEL >= int(level)


def debug_summarize_array(x: Any, max_len: int = 120) -> str:
    """Small helper for debug logs; never used in hot paths."""
    try:
        import numpy as np  # local import to keep this module light

        if isinstance(x, np.ndarray):
            return f"array(shape={x.shape}, dtype={x.dtype}, min={x.min()}, max={x.max()})"
        return f"{type(x).__name__}"
    except Exception:
        return type(x).__name__

