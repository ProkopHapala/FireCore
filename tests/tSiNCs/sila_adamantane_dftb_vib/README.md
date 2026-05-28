# DEPRECATED - Cache Files Only

This directory contains only DFTB force/energy cache files from finite difference vibrational calculations.

**These are generated artifacts, not source code.**

## Contents
- `cache.*.json` - Force/energy cache for finite difference displacements
- `cache.eq.json` - Equilibrium geometry forces

## Status
This directory is deprecated. The cache files are kept for reference but should not be used in production. Use the main `sila_adamantane_dftb_*` files in the parent directory instead.

## Recommendation
Consider adding this directory to `.gitignore` to avoid tracking cache files.
