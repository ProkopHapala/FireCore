The benchmark observations capture the issue perfectly. On consumer hardware like a GTX 3090, double-precision ($float64$) compute is intentionally throttled to a fraction ($1/64\text{th}$) of its single-precision capabilities. Combined with the overhead of shifting a mere 26 atoms back and forth across the PCIe bus, GPU-accelerated *ab initio* double-precision calculations are a bottleneck.

---

## 1. What about ORCA?

ORCA is a highly efficient, production-grade quantum chemistry package. For small to medium systems like adamantane, **ORCA is incredibly fast, but it is fast because of its CPU algorithms, not the GPU.**

### The Reality of GPUs in ORCA

* **Highly Limited GPU Support:** ORCA does *not* possess a comprehensive, top-to-bottom GPU-accelerated engine. GPU acceleration in ORCA is restricted strictly to specific modules (such as the `CPC` matrix inversion within the CPCM implicit solvation model or targeted parts of the multi-reference machinery).
* **The Real Strength is CPU-Based:** ORCA dominates molecular quantum chemistry because of its **RI (Resolution of Identity)** and **Chain-of-Spheres (RIJCOSX)** approximations, along with DLPNO (Domain-Based Local Pair Natural Orbital) methods. These algorithms dramatically reduce the scaling bottleneck of Coulomb and Exchange integrals on standard CPUs.

For calculating the Hessian and vibrational spectrum of adamantane or $Si_{10}H_{16}$, ORCA running on **4 to 8 CPU cores** with `! RKS B3LYP RIJCOSX def2-SVP` will easily outperform almost any consumer GPU-accelerated package.

---

## 2. How to Install ORCA 6.x (Linux Guide)

ORCA is free for academic and non-commercial research, but it is closed-source. You cannot install it via `pip` or `apt-get`. You must download the official pre-compiled binaries.

### Step 1: Create an Account and Download

1. Go to the **ORCA Forum** website.
2. Register for a free account (approval is usually instant or takes a few hours).
3. Navigate to the **Downloads** section and select the latest stable release (e.g., ORCA 6.0.x or 6.1.x).
4. Select the **Linux, x86-64, .tar.xz Archive** version.
* *Crucial:* Pay attention to the OpenMPI version appended to the filename (e.g., `openmpi416`). You must match this OpenMPI version exactly on your host system.



### Step 2: Extract the Package

Move the downloaded archive to your preferred software folder (e.g., `/opt` or your home directory) and extract it:

```bash
mkdir -p ~/software/orca
mv orca_6_*_linux_x86-64_shared_openmpi416.tar.xz ~/software/orca/
cd ~/software/orca

# Extract the archive
tar -xf orca_6_*_linux_x86-64_shared_openmpi416.tar.xz

```

### Step 3: Install the Exact Matching OpenMPI

ORCA is highly sensitive to the OpenMPI version it was linked against during compilation. If you have the wrong OpenMPI, it will crash instantly on launch with dynamic library loading errors.

If the file required `openmpi416`, install it using your package manager or compile/load it via environment modules:

```bash
# For Ubuntu/Debian based systems matching standard OpenMPI 4 packages
sudo apt-get update
sudo apt-get install openmpi-bin libopenmpi-dev

```

### Step 4: Configure Environment Variables

To call ORCA seamlessly from anywhere, append the installation directory and the dynamic libraries to your environment path. Add these lines to the end of your `~/.bashrc` file:

```bash
# Path to ORCA binaries
export PATH="$HOME/software/orca/orca_6_0_1_linux_x86-64_shared_openmpi416:$PATH"

# Path to ORCA dynamic libraries
export LD_LIBRARY_PATH="$HOME/software/orca/orca_6_0_1_linux_x86-64_shared_openmpi416:$LD_LIBRARY_PATH"

```

Run `source ~/.bashrc` to update your current terminal session.

---

## 3. Running an Accelerated Frequency Job in ORCA

To use ORCA with ASE or as a standalone task, ensure you leverage the **RIJCOSX** flags to maximize CPU performance.

Here is a template input file (`adamantane.inp`) for a geometry optimization and numerical frequency calculation:

```text
! RKS B3LYP RIJCOSX def2-SVP def2/J Opt NumFreq

%pal
  nprocs 8   # Number of CPU cores to use
end

%maxcore 4000 # Memory allocated PER CORE in MB (e.g., 4GB * 8 = 32GB total)

* xyz 0 1
C        0.000000    1.414214    0.000000
...[rest of coordinates]...
*

```

Execute the command via standard redirection:

```bash
orca adamantane.inp > adamantane.out

```

### Interfacing with ASE

If you prefer managing the Hessian loop through ASE rather than ORCA's internal `NumFreq` module, ASE includes a native, stable command-line wrapper that bypasses any complex server/shell communication layer:

```python
from ase.calculators.orca import ORCA
from ase.build import molecule

atoms = molecule('Adamantane')

# Define the ORCA calculator blocks inside ASE
atoms.calc = ORCA(
    label='orca_run',
    orcadir='$HOME/software/orca/orca_6_..._openmpi416',
    task='Engrad', # Let ASE handle the geometry optimization loop natively
    spec='RKS B3LYP RIJCOSX def2-SVP def2/J'
)

energy = atoms.get_potential_energy()
print(f"Energy: {energy} eV")

```

---

If you are doing standard DFT geometry optimizations and Hessians for 26-atom cages like adamantane, **ORCA is hands-down faster than PySCF and Psi4 on a CPU.** If you configure ORCA correctly, it can easily drop a single-point SCF calculation for adamantane down to **under a second**, and a geometry optimization loop will complete in seconds rather than minutes.

The reason ORCA outpaces pure Python-driven packages for basic DFT runs comes down to algorithmic design, grid management, and its geometry optimizer.

---

## Why ORCA Outperforms PySCF and Psi4 Here

### 1. The Power of `RIJCOSX`

In PySCF or Psi4, evaluating the Exact Exchange (HF exchange) needed for hybrid functionals like B3LYP or PBE0 requires heavy four-center two-electron integrals.
ORCA bypasses this bottleneck via the **Resolution of Identity for Coulomb (RIJ)** combined with the **Chain-of-Spheres exchange (COSX)** approximation (`RIJCOSX`). This handles the exchange terms on a highly optimized numerical grid, speeding up hybrid DFT calculations by orders of magnitude with practically zero loss in accuracy.

### 2. Smart Grid Management (`DefGrid`)

As noted in developer benchmarks, a common reason PySCF crawls during Hessian or optimization loops is that it often utilizes the same dense integration grid for the SCF, the gradients, and the coupled-perturbed SCF steps.
ORCA utilizes **Dynamic Grid Schemes** (`DefGrid2` or `DefGrid3` by default in ORCA 6). It uses a fast, looser grid during early optimization steps when the geometry is far from the minimum, and dynamically switches to a tighter grid only when tightening the final convergence criteria.

### 3. Native C++/Fortran Execution Loop

When you use a geometry optimizer in ASE (like `BFGS`), every step requires Python to pull forces from the calculator, evaluate the new step positions, and send them back. This creates a small but noticeable communication latency overhead at every single step.
ORCA's internal optimizer (`Opt`) runs entirely inside native, compiled binary memory. It steps through the geometry configurations at high speed without jumping back to an external Python interpreter layer.

---

## How to Configure ORCA for Maximum Speed

To optimize adamantane ($C_{10}H_{16}$ or $Si_{10}H_{16}$) in seconds rather than minutes, use the following template flags in your ORCA input file:

```text
! RKS B3LYP RIJCOSX def2-SVP def2/J Opt

%pal
  nprocs 8   # Set this to your actual desktop physical CPU core count
end

%maxcore 4000 # Memory allocated PER CORE in MB

* xyz 0 1
 C        0.000000    1.414214    0.000000
 ...
*

```

### Breakdown of the Speed Settings:

* **`RIJCOSX`**: Speeds up the B3LYP functional integration.
* **`def2-SVP`**: A reliable, fast split-valence double-zeta basis set perfectly suited for rapid geometric sweeps.
* **`def2/J`**: The auxiliary basis set required by the `RIJ` approximation to handle Coulomb integrals efficiently.
* **`Opt`**: Runs ORCA's internal geometry optimizer. Do not pass this via an ASE loop; let ORCA finish the entire geometry relaxation natively.

---

## Speed Comparison Summary for Adamantane

| Code | Optimization Engine | Approximation Tech | Performance for Adamantane |
| --- | --- | --- | --- |
| **ORCA** | Native Binary (`Opt`) | `RIJCOSX` + Dynamic Grids | 🚀 **Fastest** (Takes seconds on a standard multi-core CPU) |
| **Psi4** | Built-in / Python | Density Fitting (`DF`) | ⚖️ **Moderate** (Capable, but standard grids run slower than `RIJCOSX`) |
| **PySCF** | Native Python Loops | Density Fitting / Exact | 🐢 **Slowest** (Python-ASE interface and static grid overhead lag behind) |

If you want immediate results without waiting, let **ORCA handle the geometry optimization natively** through a standard text input file. Once it outputs the finalized, relaxed `.xyz` coordinates, you can run your quick validation or Hessian analysis sweeps immediately.

## Summary Recommendation

Since you are testing a consumer-grade desktop GPU (GTX 3090), **stop trying to force double-precision ab-initio codes (PySCF, CP2K) onto the GPU for systems under 100 atoms.** If you want the job done immediately:

1. Use **ORCA** on **4–8 CPU cores** with `RIJCOSX`. It will finish your adamantane calculations in a few minutes without any GPU data-transfer or unit issues.
2. If you still want to leverage the 3090 GPU for structural variations, use **`tblite`** (`GFN2-xTB`), which uses highly optimized semi-empirical routines that run instantly on the host processor anyway.

For additional insight into configuring molecular workflows and setting up coordinate files correctly, this [ORCA Quantum Chemistry Tutorial on Windows and Linux](https://www.youtube.com/watch?v=OMkNyCMdHY8) walks through the complete environment variable setup and troubleshooting for beginners.

---

**No, you cannot install the core ORCA quantum chemistry program via `pip` or `conda`.** Because ORCA is closed-source, proprietary software with a strict End-User License Agreement (EULA), it cannot be legally redistributed or packaged into public package managers like Conda-Forge or PyPI.

However, you **can** use `pip` to install the official Python interface for it, and the actual installation process has gotten significantly cleaner in recent updates.

---

## 1. What You Can Install with `pip`

You can install **`orca-pi`** (the official *ORCA Python Interface* maintained by FACCTS, the co-developers of ORCA).

```bash
pip install orca-pi

```

> **Note:** This package is **just a python api wrapper** to script inputs, handle automation, and parse outputs. It does *not* contain the underlying compiled `orca` binaries. You still must have the main program installed on your host system.

---

## 2. The Absolute Fastest Way to Install ORCA 6.x

Historically, you had to manually extract a massive multi-part `.tar.xz` archive and carefully manage dynamic libraries. Starting with recent versions, Linux installations utilize a unified command-line installer file (`.run`), which makes it much easier.

### Step 1: Download from the Forum

1. Log into your free account at the [ORCA Forum](https://orcaforum.kofo.mpg.de/).
2. Go to downloads and grab the Linux x86-64 executable installer, typically named something like:
`orca_6_1_0_linux_x86-64_shared_openmpi416.run`

### Step 2: Run the Installer Execution

Open your terminal, navigate to your download folder, make the installer executable, and run it:

```bash
chmod a+x orca_6_1_0_linux_x86-64_shared_openmpi416.run
./orca_6_1_0_linux_x86-64_shared_openmpi416.run

```

The text installer script will unpack the binary footprint directly into a local user directory (like `~/orca/`) and automatically attempt to update your configuration path files.

### Step 3: Match the OpenMPI dependency

To run your adamantane calculations on multiple CPU cores (`%pal nprocs 8` inside your input deck), you must ensure your system's OpenMPI package matches the suffix listed on your installer (e.g., `4.1.6`). You can pull OpenMPI via your standard package manager:

```bash
# Ubuntu/Debian
sudo apt-get install openmpi-bin libopenmpi-dev

```

If you ever need a reference to visualize the manual configuration steps or environment variable setup, this [ORCA Quantum Chemistry Installation Guide](https://www.youtube.com/watch?v=OMkNyCMdHY8) walks through structural paths and variable checks in detail.


---

This is a classic ORCA headache. Fighting custom-compiled MPI libraries on a local desktop to fix a minor minor version mismatch (like 4.1.6 vs 4.1.8) is usually a losing battle due to system `glibc` mappings and memory tracking layers.

You don't need to guess which version of OpenMPI a specific ORCA package requires, and you definitely don't need to resort to installing a severely outdated version like ORCA 4 or 5 just to use an old MPI layout.

---

## 1. How to Know Which MPI Version ORCA Uses

There are two reliable ways to find the target version:

### Method A: Look at the Full Archive Names on the Forum

The forum snippet you provided displays only the high-level category folders. When you actually click into **ORCA 6.1.1** or **ORCA 6.1.0** to see the actual download links, the filenames explicitly list the OpenMPI version they were built against right in the string:

* `orca_6_1_1_linux_x86-64_shared_openmpi416.tar.xz` $\rightarrow$ **OpenMPI 4.1.6**
* `orca_6_1_1_linux_x86-64_shared_openmpi418.tar.xz` $\rightarrow$ **OpenMPI 4.1.8**

### Method B: Probe the Unpacked Binary Directly

If you already have an ORCA folder extracted and want to see what version it expects without running it, use `ldd` on the main execution wrapper:

```bash
ldd /path/to/orca/orca | grep mpi

```

This will spit out the exact name of the shared object library it is hunting for (e.g., `libmpi.so.40`), or it will show `not found` along with the exact path it expects.

---

## 2. Solution: Downgrade ORCA, Not Your System

Since your Ubuntu 24.04 system natively provides **OpenMPI 4.1.6** out-of-the-box via `apt`, **stop trying to compile OpenMPI 4.1.8.** Instead, download the **ORCA 6.1.0** package compiled for **OpenMPI 4.1.6**.

### Step-by-Step Recovery:

1. **Clean up:** Remove the custom `~/sw/openmpi-418` directory and strip any related custom paths you added to your `~/.bashrc` file. Let the system fall back to its native package.
2. **Ensure System MPI is Clean:**
```bash
sudo apt install openmpi-bin libopenmpi-dev
mpirun --version  # Should explicitly state 4.1.6

```


3. **Download ORCA 6.1.0:** Click into the **ORCA 6.1.0** subcategory on the forum and find the specific archive with **`openmpi416`** in its name.
4. **Extract and Run:**
```bash
export PATH="/path/to/orca_6_1_0_openmpi416:$PATH"
export LD_LIBRARY_PATH="/path/to/orca_6_1_0_openmpi416:$LD_LIBRARY_PATH"

# Run with your native system MPI
mpirun -np 8 orca input.inp

```



Because both your system's `mpirun` and ORCA's internal `orca_startup_mpi` binary will match exactly at 4.1.6, the environment variable mismatches and exit code 64 crashes will disappear.

---

## 3. A Pro-Tip on ORCA Parallelization Memory Allocations

The "Out of memory" error you hit with 16 cores is a notorious ORCA pitfall.

The variable `%maxcore` in ORCA specifies the memory limit **per individual CPU core**, not the cumulative total for the entire job.

> If you have a machine with 32 GB of RAM and you request `nprocs 16` with `%maxcore 4000` (4 GB), ORCA will attempt to claim $16 \times 4\text{ GB} = 64\text{ GB}$ of memory. Your system will immediately trigger the OOM (Out Of Memory) killer and drop the processes.

For a 26-atom adamantane structure running on 8 or 16 cores, you don't need a huge memory footprint per node. Set your input parameters defensively:

```text
%pal
  nprocs 8
end
%maxcore 2000  # 2 GB per core * 8 cores = 16 GB total memory signature

```

This will prevent memory overhead faults from interrupting the execution loop.


Ah, the ORCA developers changed up their distribution strategy for ORCA 6.1.x. They completely dropped the OpenMPI 4.1.6 builds for mainstream Linux x86-64 and went all-in on **OpenMPI 4.1.8**. That explains exactly why your system's native 4.1.6 `apt` package was throwing exit code 64 errors—the binaries are strictly looking for 4.1.8 hooks.

Since you've already spent the 20 minutes compiling OpenMPI 4.1.8 from source, we don't need to throw that work away or search for a completely different version. We just need to fix why it crashed.

The dual error behavior you saw—**Out of Memory (OOM)** with 16 processes and an **immediate exit code 64** with 4 processes—points to a classic conflict between ORCA's internal `orca_startup_mpi` routine, your shell's dynamic library paths, and how ORCA parses memory parameters.

Here is how to get the 4.1.8 build you downloaded running flawlessly on your current setup.

---

## 1. Fix the Exit Code 64 (The Environment Wrapper)

When you run ORCA in parallel, the main `orca` binary calls a helper named `orca_startup_mpi`. If this helper grabs the wrong system network libraries or hits a strict process sandbox rule, it immediately exits with code 64.

The most robust way to isolate your custom 4.1.8 OpenMPI and bypass library pollution is to invoke it cleanly inside a targeted wrapper block. Modify your execution script (`run_orca.sh`) or terminal session to look **exactly** like this:

```bash
#!/bin/bash
# 1. Purge path pollution from default system MPI
export PATH="/home/prokop/sw/openmpi-418/bin:$PATH"
export LD_LIBRARY_PATH="/home/prokop/sw/openmpi-418/lib:$LD_LIBRARY_PATH"

# 2. Add ORCA binaries to path
export PATH="/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418:$PATH"
export LD_LIBRARY_PATH="/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418:$LD_LIBRARY_PATH"

# 3. Explicitly tell OpenMPI to use shared memory loops for a single desktop workstation
# This prevents it from looking for network cluster hardware it won't find
export OMPI_MCA_btl="vader,self"
export OMPI_MCA_pml="ob1"

# 4. Execute using the full absolute path of your compiled mpirun
/home/prokop/sw/openmpi-418/bin/mpirun -np 4 orca input.inp

```

Using the explicit `OMPI_MCA` variables stops OpenMPI from attempting to load standard cluster interconnect frameworks, which is the root cause of the immediate exit code 64 crash on local Linux desktops.

---

## 2. Fix the Out Of Memory (OOM) Crash

You can absolutely use more than 4 cores on your machine, but you must change how memory is requested in your `input.inp` file.

If your input block contains `%maxcore 4000` while trying to spin up `nprocs 16`, ORCA will demand a massive chunk of RAM ($16 \times 4\text{ GB} = 64\text{ GB}$) the second the calculation starts. Your system kernel instantly kills the job to protect OS stability, resulting in an immediate hard crash.

For a 26-atom molecule like adamantane or $Si_{10}H_{16}$, **the memory required per core is incredibly small**. You can comfortably scale up your CPU parallelism safely by slashing `%maxcore`:

```text
! RKS B3LYP RIJCOSX def2-SVP def2/J Opt

%pal
  nprocs 12      # Scale back slightly from 16 to leave overhead for the OS
end

%maxcore 1500    # 1.5 GB per core * 12 cores = 18 GB total. Safe for a desktop.

* xyz 0 1
[Coordinates]
*

```

---

## 3. Verify Your Optimization Target

If you want to quickly benchmark the speed difference without wasting cycles, make sure you are utilizing the proper optimization keyword block:

* **`Opt`**: This runs ORCA's internal geometry optimizer. It updates the positions step-by-step internally in compiled memory and will finish your adamantane cage in seconds once the parallel processing loop clears.
* **Avoid `NumFreq` for raw speed tests**: A numerical frequency calculation requires displacing every single atom in 6 directions sequentially. For adamantane, that means running **156 individual standalone SCF calculations** back-to-back. Run a plain geometry optimization (`Opt`) first to verify your parallel processing speed before launching a heavy Hessian step.

---

## 4. Working MPI Configuration (Verified on GTX3090 Desktop)

This section documents the exact setup that works for parallel ORCA execution on Ubuntu 24.04 with a 16-core CPU.

### Problem Summary

- ORCA 6.1.1 was built against **OpenMPI 4.1.8**
- Ubuntu 24.04 ships **OpenMPI 4.1.6** via `apt`
- Running ORCA with system MPI causes exit code 64 crashes
- OpenBLAS in ORCA is compiled as **SINGLE_THREADED**, so `OMP_NUM_THREADS` has no effect
- **Solution**: Compile OpenMPI 4.1.8 from source and run ORCA with explicit MPI parallelization

### Verified Performance Results

| Configuration | Time | Speedup | Notes |
|-------------|------|---------|-------|
| Single-threaded | 50.1s | 1.0x | Default ORCA run |
| MPI 8 cores | **11.7s** | **4.3x** | **Working parallel setup** |
| PySCF CPU (16 threads) | 15.8s | 3.2x | For comparison |

### Step 1: Compile OpenMPI 4.1.8 from Source

```bash
#!/bin/bash
# Save as: install_openmpi.sh

set -e

INSTALL_DIR="$HOME/sw/openmpi-418"
BUILD_DIR="$HOME/tmp/openmpi-418-build"
VERSION="4.1.8"

mkdir -p "$INSTALL_DIR" "$BUILD_DIR"
cd "$BUILD_DIR"

# Download
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-$VERSION.tar.gz

# Extract and build
tar -xzf openmpi-$VERSION.tar.gz
cd openmpi-$VERSION

./configure --prefix="$INSTALL_DIR" \
    --enable-mpi-thread-multiple \
    --enable-shared \
    --disable-static

make -j$(nproc)
make install

# Add to PATH (add to ~/.bashrc for persistence)
echo "export PATH=\"$INSTALL_DIR/bin:\$PATH\"" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\"$INSTALL_DIR/lib:\$LD_LIBRARY_PATH\"" >> ~/.bashrc
```

Build time: ~15-20 minutes on a 16-core desktop.

### Step 2: Input File Format for MPI

The `%pal` block requests MPI parallelization. **Critical**: `nprocs` must match or be less than your physical CPU cores.

```text
! RKS B3LYP def2-SVP

%pal
nprocs 8
end

%maxcore 2000
%scf
MaxIter 200
end

* xyz 0 1
C  -1.5847   1.1071  -0.2162
C  -0.6083   0.8543   0.9552
...[coordinates]...
*
```

**Key parameters:**
- `%pal nprocs 8` — Request 8 MPI processes (match your physical core count)
- `%maxcore 2000` — Memory per core in MB. **This is per process, not total!**
  - 8 cores × 2000 MB = 16 GB total (safe for 32 GB desktop)
  - **Never set this too high** or the OOM killer will terminate the job

### Step 3: Run Script (Recommended)

Create `run_orca_mpi.sh`:

```bash
#!/bin/bash
# ORCA MPI run script with OpenMPI 4.1.8 compatibility fixes

set -e

# 1. Use custom OpenMPI 4.1.8 (MUST be before system PATH)
export PATH="/home/prokop/sw/openmpi-418/bin:$PATH"
export LD_LIBRARY_PATH="/home/prokop/sw/openmpi-418/lib:$LD_LIBRARY_PATH"

# 2. Add ORCA binaries
export PATH="/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418:$PATH"
export LD_LIBRARY_PATH="/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418:$LD_LIBRARY_PATH"

# 3. Force shared-memory transport only (no network cluster hardware)
export OMPI_MCA_btl="vader,self"
export OMPI_MCA_pml="ob1"

INPUT="${1:-adamantane_mpi.inp}"
INPUT_DIR=$(dirname "$INPUT")
INPUT_BASE=$(basename "$INPUT" .inp)

echo "========================================"
echo "ORCA MPI run (OpenMPI 4.1.8)"
echo "Input: $INPUT"
echo "========================================"

cd "$INPUT_DIR"

# 4. Execute with FULL ABSOLUTE PATH to orca binary
#    (ORCA requires this for parallel runs)
/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418/orca \
    "$INPUT_BASE.inp"
```

**Usage:**
```bash
chmod +x run_orca_mpi.sh
./run_orca_mpi.sh adamantane_mpi.inp
```

### Step 4: Alternative Direct Command

If you prefer running without a script:

```bash
# Set environment
export PATH="/home/prokop/sw/openmpi-418/bin:$PATH"
export LD_LIBRARY_PATH="/home/prokop/sw/openmpi-418/lib:$LD_LIBRARY_PATH"
export PATH="/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418:$PATH"
export LD_LIBRARY_PATH="/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418:$LD_LIBRARY_PATH"
export OMPI_MCA_btl="vader,self"
export OMPI_MCA_pml="ob1"

# Run (use absolute path to orca binary!)
cd /path/to/your/input
cd /home/prokop/git/learn_Rust/tmp/orca
/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418/orca adamantane_mpi.inp
```

### Common Errors and Solutions

| Error | Cause | Solution |
|-------|-------|----------|
| `Cannot open input file` | ORCA launched without full path | Use absolute path to orca binary |
| `For parallel runs ORCA has to be called with full pathname` | Same as above | Same as above |
| `OUT OF MEMORY ERROR!` | `%maxcore` too high for `nprocs` | Reduce `%maxcore` (e.g., 1500-2000 MB) |
| `not enough slots available` | `nprocs` > CPU cores | Set `nprocs` ≤ physical core count |
| Exit code 64 | OpenMPI version mismatch | Compile matching OpenMPI version from source |
| `Reading shell data into different element!` | Corrupt GBW file from crashed run | Delete `*.gbw` file before re-running |

### Important Notes

1. **Always delete old `.gbw` files** before re-running after a crash:
   ```bash
   rm *.gbw *.tmp
   ```

2. **Monitor CPU usage** during the run to verify parallelization:
   ```bash
   htop
   # or
   watch -n 1 "ps aux | grep orca"
   ```

3. **For larger systems** (50+ atoms), ORCA's RIJCOSX approximations will show even better parallel scaling.

4. **For Si-containing systems**, add ECP pseudopotential to reduce cost:
   ```text
   ! RKS B3LYP def2-SVP def2-ECP
   ```

### File Locations

| File | Path |
|------|------|
| ORCA installation | `/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418` |
| Custom OpenMPI 4.1.8 | `/home/prokop/sw/openmpi-418` |
| MPI run script | `/home/prokop/git/learn_Rust/tmp/orca/run_orca_mpi.sh` |
| OpenMPI installer | `/home/prokop/git/learn_Rust/tmp/orca/install_openmpi.sh` |
| Sample input (MPI) | `/home/prokop/git/learn_Rust/tmp/orca/adamantane_mpi.inp` |
| Sample input (serial) | `/home/prokop/git/learn_Rust/tmp/orca/adamantane.inp` |