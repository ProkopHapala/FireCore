# Local-First Coinage-Metal DFT Database Campaign for GridFF, With Rigid and Relaxed Scan Families

## Summary
- Build this as a **database-creation campaign** with two execution tiers:
  - **local workstation** for protocol validation, Ag pilot, and universal-slab selection
  - **HPC** for the full production matrix across `Ag`, `Cu`, and `Au`
- Use **`PBE + vdWsurf` as the primary production protocol** and keep **`RPBE-D3(BJ)` as a smaller audit subset** after the main workflow is stable.
- Use **`3x3x4` as the fixed production slab** for the current HPC campaign, following the verified Ag reference workflow in `ORR_HER_Ag_Colab/results`.
- The database must contain **two scan families** for every retained adsorption minimum:
  - **rigid scan**: rigid substrate, rigid molecule
  - **relaxed scan**: rigid substrate, molecule relaxed at each scan point
- For GridFF, the **rigid scan family is the primary fitting target**. The relaxed scan family is a required benchmark/sensitivity layer, not the baseline GridFF target.

## Scientific Scope
- Metals:
  - `Ag(111)`, `Cu(111)`, `Au(111)`
- Adsorbates:
  - `H`, `CO`, `H2O`, `NH3`, `methanol`, `HCONH2`, `pyridine`
- Representative adsorbates for universal-slab selection:
  - `H`
  - `CO`
  - `H2O`
  - `HCONH2`
  - `pyridine`
- Production slab:
  - `3x3x4`
- Reference source:
  - use the verified Ag workflow under `ORR_HER_Ag_Colab/results/Ag_ORR_HER`
  - bulk templates follow the verified `Ag_bulk`
  - slab/adsorbate phase settings follow the verified `phase1 -> phase2 -> final_scf_12x12x1 -> workfunc_12x12x1 -> bader_12x12x1` pattern

## Local vs HPC Boundary
- Local workstation is used only for:
  - bulk template validation
  - Ag bulk convergence
  - Ag universal-slab selection
  - Ag clean-slab final `CHGCAR/LOCPOT`
  - a narrow Ag representative adsorption pilot
  - one rigid and one relaxed scan pilot per representative adsorbate
- HPC is required for:
  - full Ag adsorption library
  - all Cu and Au production jobs
  - full rigid and relaxed scan databases
  - volumetric subsets for all retained minima
- Keep one manifest schema and one directory layout for both local and HPC; only the launcher changes.

## Database Structure and Reuse
- External campaign root:
  - `/home/niel/git/coinage_gridff_dft`
- Repo-side orchestration only:
  - `tests/tGridFFJax/coinage_campaign/`
- Campaign phases:
  - `00_bulk`
  - `01_universal_slab_selection`
  - `02_clean_slab_final`
  - `03_gas_refs`
  - `04_ads_seed_library`
  - `05_ads_relax`
  - `06_scan_rigid`
  - `07_scan_relaxed`
  - `08_volumetrics`
  - `09_analysis`
- Reuse from `/home/niel/git/ORR_HER_Ag_Colab`:
  - bulk-to-slab construction logic
  - multi-stage metallic slab relaxation logic
  - dipole-corrected slab/static workflow
  - convergence-study structure
- Treat existing Ag templates as references only, not production defaults.

## Fixed Slab Reference
- The current production campaign uses the verified `Ag(111) 3x3x4` slab pattern from `ORR_HER_Ag_Colab/results`.
- Clean slab stages:
  - `relax_stage1_nodipole`
  - `relax_stage2_cycle*_dipole`
  - `final_static`
  - `workfunction`
  - `bader`
- Bulk is generated separately and used to supply the lattice constant to slab construction.
- Bulk uses no surface-screened dispersion tags in the generated INCAR.

## Adsorption Search and Retained Minima
- Search sites:
  - `top`, `bridge`, `fcc`, `hcp`
- Use exhaustive seed-library generation at the setup stage.
- Use MAD-SURF only for:
  - seed prescreening
  - duplicate pruning
  - exact-geometry comparison after DFT data exists
- Run DFT relaxation on all survivors within the MAD-SURF energy window.
- Keep all unique relaxed minima within the production retention window.

## Scan Families
- Both scan families use a **rigid slab**.
- **Rigid scan**
  - molecule geometry and orientation fixed
  - scan only along the surface normal from the retained adsorption minimum
  - this is the primary GridFF fitting dataset
- **Relaxed scan**
  - slab remains fixed
  - at each scan point, constrain the molecule’s scan coordinate and lateral site anchor
  - allow the molecule to relax internally and reorient
  - this is a benchmark/sensitivity dataset for how much rigid GridFF misses due to molecular relaxation
- Common scan range:
  - `2.0–12.0 Å`
- Common scan spacing:
  - `2.0–3.5 Å` in `0.10 Å`
  - `3.75–6.0 Å` in `0.25 Å`
  - `6.5–12.0 Å` in `0.50 Å`

## GridFF-Relevant Targets
- For each scan geometry, store:
  - DFT total energy
  - DFT adsorbate forces
  - MAD-SURF energy and forces on the exact same geometry
- Primary GridFF fitting targets:
  - rigid-scan interaction energy
  - rigid-scan interaction forces
- Derived interaction quantities:
  - `E_int = E_slab+mol - E_slab - E_mol`
  - `F_int = F_slab+mol(on adsorbate atoms) - F_mol`
- Relaxed scans are not used as the primary strict-GridFF fit target in phase 1.
- Relaxed scans are used to measure:
  - where rigid GridFF remains sufficient
  - where molecule relaxation contributes a systematic residual
  - which systems may later need a relaxation-aware correction layer

## Data Products
- Clean slab:
  - relaxed geometry
  - final static energy
  - `CHGCAR`
  - `LOCPOT` from `workfunction`
  - `AECCAR0/AECCAR2` from `bader`
  - work function
- Adsorption minima:
  - relaxed geometry
  - total energy
  - adsorbate forces
  - selected volumetric outputs
- Rigid scan:
  - geometry series
  - total and interaction energies
  - adsorbate and interaction forces
  - MAD-SURF comparison on identical poses
- Relaxed scan:
  - geometry series
  - relaxed molecular geometries
  - total and interaction energies
  - adsorbate and interaction forces
  - MAD-SURF comparison on identical constrained geometries
- Analysis outputs:
  - DFT vs MAD-SURF parity
  - rigid vs relaxed scan differences
  - metal-by-metal and adsorbate-by-adsorbate residual tables
  - GridFF fit inputs for rigid scans

## Assumptions
- This project’s goal is a **reusable DFT database for GridFF**, not only a single Ag test.
- One universal slab is required for the production campaign.
- The workstation is a **protocol-validation platform**, not the final production engine.
- The substrate remains rigid in both scan families because that matches the current GridFF framework.
- The molecule is:
  - rigid in the rigid scan family
  - relaxed in the relaxed scan family
- The rigid scan family is the **phase-1 GridFF fitting target**; the relaxed scan family is the required benchmark layer for assessing what rigid GridFF cannot capture.

## Implemented Orchestration
- package:
  - `tests/tGridFFJax/coinage_campaign/`
- entry points:
  - `setup_coinage_campaign.py`
  - `setup_scans_from_minima.py`
  - `setup_interaction_references.py`
- current HPC production generation:
  - `setup_coinage_campaign.py --selected-slab 3x3x4`
  - materializes `02_clean_slab_final`, `04_ads_seed_library`, and `05_ads_relax`
- current relaxed-scan implementation:
  - rigid slab
  - one chemically meaningful adsorbate anchor atom fixed at the selected site and scan z
  - remaining adsorbate atoms relaxed
  - this is the VASP-compatible realization of the “rigid substrate, relaxed molecule” scan family

## Example Commands
```bash
python3 tests/tGridFFJax/coinage_campaign/setup_coinage_campaign.py \
  --campaign-root /tmp/coinage_gridff_dft_dryrun

python3 tests/tGridFFJax/coinage_campaign/setup_coinage_campaign.py \
  --campaign-root /tmp/coinage_gridff_dft_dryrun \
  --selected-slab 3x3x4

python3 tests/tGridFFJax/coinage_campaign/setup_scans_from_minima.py \
  --minimum-dir /tmp/coinage_gridff_dft_dryrun/01_universal_slab_selection/Ag/3x3x4/representative_adsorbates/CO/C_down_ontop_180_180/final_static \
  --family both \
  --out-root /tmp/coinage_gridff_scans

python3 tests/tGridFFJax/coinage_campaign/setup_interaction_references.py \
  --scan-root /tmp/coinage_gridff_scans_rigid \
  --out-root /tmp/coinage_gridff_scan_refs
```
