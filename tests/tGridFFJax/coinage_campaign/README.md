# Coinage Campaign Orchestration

This directory contains the current DFT database orchestration for the coinage-metal GridFF project. The HPC defaults are aligned to the verified Ag reference under `ORR_HER_Ag_Colab/results`.

Main entry points:

- `setup_coinage_campaign.py`
  - creates the external campaign tree
  - writes bulk, Ag representative 3x3x4 reference jobs, gas references, and when `--selected-slab 3x3x4` is provided also writes:
    - `02_clean_slab_final`
    - `04_ads_seed_library`
    - `05_ads_relax`
- `run_campaign_queue.py`
  - walks the generated tree sequentially, preferring `run_pipeline_local.sh` at each job root
  - skips stages that already look complete
  - supports `--pilot-only`, `--phase`, `--dry-run`, and `--wait-if-busy`
- `setup_scans_from_minima.py`
  - derives `rigid`, `relaxed`, `relaxed_slab`, `both`, or `all` scan families from one converged minimum or an entire minima root
  - writes a top-level `phase_manifest.json` for queue-aware follow-up
- `setup_interaction_references.py`
  - derives `slab_only` and `molecule_only` reference jobs from completed scan geometries
  - writes a top-level `phase_manifest.json`

Current practical definition of the three scan families:

- `rigid`
  - slab fixed
  - all adsorbate atoms fixed
  - static single-point job at each z
- `relaxed`
  - slab fixed
  - one chemically meaningful anchor atom fixed at the site and z
  - remaining adsorbate atoms free
  - relaxation job plus final static at each z
- `relaxed_slab`
  - bottom slab layers fixed
  - top slab layers free
  - one chemically meaningful anchor atom fixed at the site and z
  - remaining adsorbate atoms free
  - relaxation job plus final static at each z

The `relaxed` family is the VASP-compatible implementation of the plan’s “rigid slab, relaxed molecule” scan family. The `relaxed_slab` family is the higher-cost benchmark tier for “top slab layers + molecule relaxed”.

Current reference-aligned defaults:

- slab: `3x3x4`
- HPC launcher:
  - slab, adsorbate, and scan jobs: `1` node, `96` CPUs
  - gas-reference jobs: `1` node, `4` CPUs
  - direct `mpiexec` launch from `run.pbs`, aligned to the ORR results-folder pattern
- INCAR parallel settings:
  - bulk / slab / adsorbate / scan: `NCORE = 16`
  - gas references: `NCORE = 1`
  - atomic `H` gas reference: `ISPIN = 2`, `NUPDOWN = 1`, `MAGMOM = 1`
- HPC binaries:
  - slab / adsorbate / scan / bulk: `vasp_std`
  - gamma-only gas references: `vasp_gam`
- bulk: no surface-screened dispersion tags
- clean slab / adsorbates:
  - `relax_stage1_nodipole`
  - `relax_stage2_cycle*_dipole`
  - `final_static`
  - `workfunction`
  - `bader`
- canonical campaign molecule name: `HCONH2`
