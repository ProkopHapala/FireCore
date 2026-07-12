---
description: Inventory existing code before writing new code to avoid duplication and find reusable components
---

## Purpose

Prevent reimplementation. Before writing any new function, class, or module, confirm the topic does not already exist somewhere in the codebase.

## Workflow

1. **Query the topic**
   - Search `doc/topical_audit/` for the scientific topic (e.g., "gridff", "afm", "kpoints").
   - Search `CODEMAP.md` for file locations and module relationships.
   - Search the codebase (grep) for function/class names matching your planned implementation.

2. **Inventory existing implementations**
   - List all locations where similar logic exists (Fortran, C++, Python, OpenCL, etc.).
   - Note status: active, experimental, deprecated, unfinished.
   - Check if an existing function can be generalized to fit your needs.

3. **Assess reuse potential**
   - If an existing implementation fits: generalize it instead of duplicating.
   - If generalization requires risky architectural changes that threaten backward compatibility: **stop and report for approval**.
   - If no suitable implementation exists: proceed, but document the new topic in `doc/topical_audit/`.

4. **Verify separation of concerns**
   - Ensure compute logic is not mixed with plotting/diagnostics.
   - Ensure backend modules are not entangled with GUI/CLI/test scripts.
   - If you need plotting/debugging: check `pyBall/` for shared utilities before writing ad-hoc code.
