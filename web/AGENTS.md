# Web & Visualization

## Purpose

WebGL/WebGPU-based molecular visualization, physics simulation, and interactive GUI components. Browser-based frontends for molecular dynamics, force-field editing, and STM/AFM visualization.

## Ownership

- `molgui_web/` — WebGL molecular GUI (legacy)
- `molgui_webgpu/` — WebGPU molecular GUI and physics simulation
- `NBody2D_WebGPU/` — 2D N-body WebGPU simulation
- `common_js/` — Shared JavaScript utilities
- `ppstm_web/` — Web-based STM visualization

## Local Contracts

- **Run from this directory** — web assets use relative paths.
- **JavaScript fail loudly:** Throw errors instead of providing silent defaults or suppressing issues. Never add geometry-based or default-value fallbacks that hide missing parameter data.
- **Debuggability > UX:** If required data (bond/angle/atom types, l0/k, Ass/Kss, etc.) is missing/invalid, print a descriptive diagnostic and throw an Error.

## Work Guidance

### WebGPU Molecular GUI
- `molgui_webgpu/` — Main WebGPU frontend for molecular visualization and physics
- Uses XPBD/MMFFL force fields in compute shaders
- Supports real-time MD, editing, and parameter inspection

### WebGL (Legacy)
- `molgui_web/` — Original WebGL implementation
- Still used for some visualization tasks

### N-Body 2D
- `NBody2D_WebGPU/` — 2D gravitational/coulomb N-body simulation

## Verification

- Serve via `run_web_server.sh` from repo root and test in browser
- Check browser console for loud errors (should not silently fail)

## Child DOX Index

- `molgui_webgpu/` — WebGPU molecular physics frontend
- `molgui_web/` — WebGL molecular visualization frontend
- `NBody2D_WebGPU/` — 2D N-body WebGPU simulation
- `ppstm_web/` — Web-based STM visualization
