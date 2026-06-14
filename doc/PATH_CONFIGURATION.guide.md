# FireCore Path Configuration Guide

This document explains how to configure paths for FireCore on different systems. FireCore uses a flexible configuration system that supports both `config.json` files and environment variables.

## Problem

FireCore previously had hardcoded paths in the code (e.g., `/home/prokop/SIMULATIONS/...`) which made it difficult to use on different systems. The code was looking for basis files at specific locations like `pyBall/DFTB/data/wfc.mio-1-1.hsd` which didn't exist.

## Solution

FireCore now uses a hierarchical configuration system:

1. **Environment variables** (highest priority)
2. **config.json file** (medium priority)
3. **Hardcoded defaults** (lowest priority, for backward compatibility)

## Configuration Methods

### Method 1: Using config.json (Recommended)

Create or edit `firecore_config.json` in the FireCore root directory:

```json
{
  "paths": {
    "dftb_sk_path": "/path/to/your/dftbplus/slakos",
    "dftb_exe": "/path/to/dftb+/dftb+",
    "dftb_basis_path": "/path/to/FireCore/pyBall/DFTB/data",
    "fdata_dir": "/path/to/FireCore/tests/pyFireball/Fdata",
    "firecore_root": "/path/to/FireCore"
  },
  "dftb": {
    "default_sk_set": "mio-1-1",
    "available_sk_sets": ["mio-1-1", "3ob-3-1"],
    "basis_sets": {
      "mio-1-1": {
        "sk_path": "/path/to/slakos/mio/mio-1-1",
        "wfc_path": "/path/to/FireCore/pyBall/DFTB/data/wfc.mio-1-1.hsd"
      },
      "3ob-3-1": {
        "sk_path": "/path/to/slakos/3ob-3-1",
        "wfc_path": "/path/to/FireCore/pyBall/DFTB/data/wfc.3ob-3-1.hsd"
      }
    }
  }
}
```

### Method 2: Using Environment Variables

Set environment variables in your shell (e.g., in `~/.bashrc` or `~/.zshrc`):

```bash
# DFTB+ Slater-Koster files directory
export DFTB_SK_PATH=/home/yourname/SIMULATIONS/dftbplus/slakos

# DFTB+ executable (optional, only needed for DFTB backend)
export DFTB_EXE=/path/to/dftb+/dftb+

# DFTB basis files directory (wfc.*.hsd files)
export DFTB_BASIS_PATH=/home/yourname/git/FireCore/pyBall/DFTB/data

# FireCore root directory (optional, auto-detected)
export FIRECORE_ROOT=/home/yourname/git/FireCore

# Fdata directory (optional, for AFM simulations)
export FDATA_DIR=/home/yourname/git/FireCore/tests/pyFireball/Fdata
```

**Note:** Environment variables override config.json settings.

## Required Data Files

### 1. Slater-Koster (SK) Files

Download DFTB+ parameter sets from https://dftb.org/parameters/download.html

Example for mio-1-1:
```bash
cd ~/SIMULATIONS/dftbplus/slakos
wget https://github.com/dftbparams/mio/releases/download/v1.1/mio-1-1.tar.xz
tar xf mio-1-1.tar.xz
```

Example for 3ob-3-1:
```bash
cd ~/SIMULATIONS/dftbplus/slakos
wget https://github.com/dftbparams/3ob/releases/download/v3.1.0/3ob-3-1.tar.xz
tar xf 3ob-3-1.tar.xz
```

### 2. Wavefunction Coefficient (WFC) Files

The WFC files (`wfc.*.hsd`) contain Slater-type orbital basis parameters. They are typically found in the DFTB+ source distribution or in the `recipes/slakos/wfc/` directory.

**For your system:**
```bash
# Copy WFC files to FireCore data directory
cp ~/SIMULATIONS/dftbplus/recipes/slakos/wfc/wfc.mio-1-1.hsd ~/git/FireCore/pyBall/DFTB/data/
```

If you don't have the WFC files, you can:
1. Download DFTB+ source from https://github.com/dftbplus/dftbplus
2. Extract and find them in `recipes/slakos/wfc/`
3. Or use the waveplot_in.hsd files from DFTB+ test directories

### 3. CO Tip Geometry File (for AFM simulations)

The CO tip geometry file is required for AFM/STM simulations. This file defines the CO molecule used as the AFM tip.

**Location:** `tests/tAFM/pyocl_fdbm/CO.xyz`

**If missing, copy from:**
```bash
cp ~/git/FireCore/cpp/common_resources/xyz/CO.xyz ~/git/FireCore/tests/tAFM/pyocl_fdbm/CO.xyz
```

Alternative locations in the repository:
- `cpp/common_resources/xyz/CO.xyz`
- `tests/pyFireball/relaxed_mols/CO.xyz`
- `tests/tEFF/export/molecules/CO.xyz`

## Current Configuration for carbsisYoga

Based on your system, the current `firecore_config.json` is configured as:

```json
{
  "paths": {
    "dftb_sk_path": "/home/prokophapala/SIMULATIONS/dftbplus/slakos",
    "dftb_basis_path": "/home/prokophapala/git/FireCore/pyBall/DFTB/data",
    "fdata_dir": "/home/prokophapala/git/FireCore/tests/pyFireball/Fdata",
    "firecore_root": "/home/prokophapala/git/FireCore"
  },
  "dftb": {
    "basis_sets": {
      "mio-1-1": {
        "sk_path": "/home/prokophapala/SIMULATIONS/dftbplus/slakos/mio/mio-1-1",
        "wfc_path": "/home/prokophapala/git/FireCore/pyBall/DFTB/data/wfc.mio-1-1.hsd"
      },
      "3ob-3-1": {
        "sk_path": "/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1",
        "wfc_path": "/home/prokophapala/git/FireCore/pyBall/DFTB/data/wfc.3ob-3-1.hsd"
      }
    }
  }
}
```

## Testing Your Configuration

You can test if the configuration is working by running:

```python
from pyBall.config_utils import print_config
print_config()
```

Or check specific paths:

```python
from pyBall.config_utils import get_path, get_dftb_basis_path

# Check SK path
sk_path = get_path('dftb_sk_path')
print(f"SK path: {sk_path}")

# Check basis file for mio-1-1
basis_path = get_dftb_basis_path('mio-1-1')
print(f"mio-1-1 basis: {basis_path}")
```

## Troubleshooting

### Error: "Basis file not found: pyBall/DFTB/data/wfc.mio-1-1.hsd"

**Cause:** The WFC file is missing or the path is incorrect.

**Solution:**
1. Check if the file exists: `ls ~/git/FireCore/pyBall/DFTB/data/`
2. Copy the WFC file from your DFTB+ installation
3. Update `firecore_config.json` to point to the correct location
4. Or set `DFTB_BASIS_PATH` environment variable

### Error: "SK library not found"

**Cause:** The Slater-Koster files directory is not found.

**Solution:**
1. Check if `DFTB_SK_PATH` is set: `echo $DFTB_SK_PATH`
2. Verify the directory exists and contains `.skf` files
3. Update `firecore_config.json` with the correct path
4. Download the parameter sets if missing

### Error: "DFTB_EXE not set"

**Cause:** DFTB+ executable path is not configured (only needed for DFTB backend).

**Solution:**
1. Install DFTB+ from https://github.com/dftbplus/dftbplus
2. Set `DFTB_EXE` in config.json or environment variable
3. Or use pySCF backend instead (doesn't require DFTB+)

## Porting to Other Systems

To set up FireCore on a different system:

1. **Copy the FireCore repository** to the new system
2. **Download required data files** (SK files, WFC files)
3. **Create/Edit `firecore_config.json`** with paths for the new system
4. **Optionally set environment variables** if you prefer that method
5. **Test configuration** using the Python commands above

Example for a new system at `/local/data/user/FireCore`:

```json
{
  "paths": {
    "dftb_sk_path": "/local/data/params/dftbplus/slakos",
    "dftb_basis_path": "/local/data/user/FireCore/pyBall/DFTB/data",
    "fdata_dir": "/local/data/user/FireCore/tests/pyFireball/Fdata",
    "firecore_root": "/local/data/user/FireCore"
  }
}
```

## Code Changes Made

The following files were modified to support the new configuration system:

1. **`pyBall/config_utils.py`** (new file)
   - Configuration loading utility
   - Supports config.json and environment variables
   - Helper functions for getting specific paths

2. **`pyBall/dftb_utils.py`**
   - Updated `_check_sk_path()` to use config system
   - Replaced hardcoded `SK_PATHS` and `WFC_HSD_PATHS` with config-based loading
   - Maintains backward compatibility with hardcoded paths

3. **`pyBall/OCL/ModularPipeline.py`**
   - Added import for config utilities
   - Updated basis file path resolution to use config system
   - Falls back to hardcoded paths if config fails

## Summary

- **Primary method:** Edit `firecore_config.json` in FireCore root
- **Alternative method:** Set environment variables
- **Environment variables override config.json**
- **Backward compatible:** Falls back to old hardcoded paths if config fails
- **Your system is now configured** with paths from `/home/prokophapala/SIMULATIONS/dftbplus/slakos`
