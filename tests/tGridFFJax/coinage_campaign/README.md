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
  - derives `rigid`, `relaxed`, or `both` scan families from one converged minimum or an entire minima root
  - writes a top-level `phase_manifest.json` for queue-aware follow-up
- `setup_interaction_references.py`
  - derives `slab_only` and `molecule_only` reference jobs from completed scan geometries
  - writes a top-level `phase_manifest.json`

Current practical definition of the two scan families:

- `rigid`
  - slab fixed
  - all adsorbate atoms fixed
  - static single-point job at each z
- `relaxed`
  - slab fixed
  - one chemically meaningful anchor atom fixed at the site and z
  - remaining adsorbate atoms free
  - relaxation job plus final static at each z

This anchor-fixed relaxed scan is the VASP-compatible implementation of the plan’s “rigid slab, relaxed molecule” scan family.

Current reference-aligned defaults:

- slab: `3x3x4`
- HPC launcher:
  - `1` node
  - `64` CPUs
  - direct `mpiexec` launch from `run.pbs`, aligned to the ORR results-folder pattern
- bulk: no surface-screened dispersion tags
- clean slab / adsorbates:
  - `relax_stage1_nodipole`
  - `relax_stage2_cycle*_dipole`
  - `final_static`
  - `workfunction`
  - `bader`
- canonical campaign molecule name: `HCONH2`
