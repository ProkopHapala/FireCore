# MAD-SURF

MAD-SURF is a machine learning interatomic potential specifically tailored for molecular adsorption on coinage metal
surfaces (Cu, Ag, Au). It is trained on a broad dataset covering diverse organic molecules, adsorption motifs,
surfaces, molecular dynamics trajectories and non‑covalent aggregates. The resulting potential reproduces DFT
reference energies, forces, and adsorption geometries at a fraction of the computational cost, enabling large‑scale,
data‑driven simulations of complex interfaces such as organic monolayers, polycyclic aggregates, flexible biomolecules,
and reconstructed metal surfaces.

This repository contains:
- The trained MAD-SURF models.
- The training and evaluation dataset.
- Example scripts and notebooks for using the potential, including geometry optimization of adsorbate–surface systems.

If you use MAD-SURF in your work, please cite the associated paper (MAD-SURF: a general ML interatomic potential for
molecular adsorption on coinage metal surfaces; full citation details to be added here).

---

## Repository layout

- `madsurf/` – core code for the MAD-SURF paper
- `dataset/` – training and test data in `extxyz` format.
- `models/` – pre-trained MAD-SURF models (downloaded via `zenodo_download.sh`).
- `examples/` – example notebooks and scripts showcasing use cases for MAD-SURF:
  - `run_geom_op.ipynb` – end-to-end example of setting up and relaxing an adsorbate on a metal surface.
- `zenodo_download.sh` – helper script to download dataset and model files from Zenodo.

---

## Installation

Create a Python environment (e.g. with `conda` or `venv`) and install the required packages. At minimum you will need:

- Python ≥ 3.9
- PyTorch (CPU or GPU build)
- [`mace`](https://github.com/ACEsuit/mace) / `mace-torch`
- `ase`
- `aalto-boss`
- `numpy`, `scipy`, `matplotlib`

A minimal example with `pip` (adapt to your setup):

```bash
pip install torch  # choose the build appropriate for your hardware
pip install mace-torch ase aalto-boss numpy scipy matplotlib
```

In addition, we highly suggest you install the libraries for CUDA acceleration (if building for GPU) for faster training and inference of MACE:

```bash
pip install cuequivariance cuequivariance-torch cuequivariance-ops-torch-cu12 
```
where 'cu12' corresponds to your CUDA-version, e.g. 12 in this case.

---

## Downloading data and models from Zenodo

The script `zenodo_download.sh` automates downloading the dataset and pre‑trained MAD-SURF models from the Zenodo
record **18312238**.

From the repository root (`CODE/MAD-SURF`):

```bash
chmod +x zenodo_download.sh
./zenodo_download.sh
```

Default behaviour (no flags):
- If `./dataset` already exists: the dataset download/extraction is skipped.
- If `./models` already exists: the model download/extraction is skipped.
- Otherwise:
  - Downloads `dataset.zip` and extracts it into `./dataset`.
  - Downloads the two core model files into `./models/`:
    - `MAD-SURF.model`
    - `MAD-SURF_fewshot.model`

To download all available models instead of just the two core ones:

```bash
./zenodo_download.sh --all_models
```

This will download `models.zip` and extract it into `./models` (unless that directory already exists).


## License

Apache 2.0 License.

