---
description: Consolidation and refactoring review after coding to remove duplication and enforce separation of concerns
---

## Purpose

After implementing a feature, review the codebase for duplication, ad-hoc utilities, and backend/frontend entanglement.

## Workflow

1. **Find duplication**
   - Search for functions/classes with similar names or logic across languages (Fortran ↔ C++ ↔ Python).
   - Identify legacy code branches that are superseded by your new implementation.
   - Mark old implementations as deprecated in `doc/topical_audit/`.

2. **Separate backend from utilities**
   - Move plotting, diagnostics, and visualization out of core compute modules into shared utilities (`pyBall/plot_utils.py` or similar).
   - Move GUI/CLI-specific logic out of backend libraries.
   - Ensure test scripts do not reimplement debug/plotting functions that belong in shared modules.

3. **Consolidate test scripts**
   - Check if your new test duplicates an existing test harness.
   - Factor common test logic into a unified CLI with routing parameters rather than isolated scripts.
   - Ensure tests use shared utilities for diagnostics and plotting.

4. **Update documentation**
   - Add/update the topic in `doc/topical_audit/` with new implementation location and status.
   - Update `CODEMAP.md` if module boundaries changed.

## Stop Conditions

- Do NOT consolidate if it breaks backward compatibility without explicit approval.
- Do NOT move compute logic into plotting modules or vice versa.
