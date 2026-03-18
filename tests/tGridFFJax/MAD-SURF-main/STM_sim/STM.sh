#!/bin/bash
#SBATCH --job-name=FDBM_stm
#SBATCH --chdir=.
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=64
#SBATCH -N 1
#SBATCH -p low_priority
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=3900
#SBATCH --get-user-env
#NOTE: 
# - '#SBATCH' is NOT a commented line
# - '##SBATCH' this one YES
# - '# SBATCH' this one also.

# [DO NOT CHANGE]
# Use '/scratch' FS - file system
# -------------------------------
source scratchfs
# -------------------------------
# SYNC DATA: uwd SlurmJobName.out
# -------------------------------

module load intel/oneAPI/conda/ vasp/6.3.0
export VASP_PP_PATH='/apps/exported/installed/software/vasp/PP/'
export DBSPM_PATH='/home/mgonzalez/DBSPM/'

# --- 1) DBSPM sample + VASP ---

TOTAL_CORES=$((SLURM_NTASKS * SLURM_CPUS_PER_TASK))
echo "Total cores available: $TOTAL_CORES"

mpirun -n 1 "$DBSPM_PATH/dbspm" sample -i input.in \
    --dft-cmd="mpirun -np $TOTAL_CORES vasp_std" > prepdft_sample.out

# --- 2) Prepare STM folder from LDIPOL outputs ---

mkdir -p sample/STM

for f in POSCAR INCAR WAVECAR CONTCAR KPOINTS POTCAR; do
    cp "sample/LDIPOL/$f" "sample/STM/$f"
done

cd sample/STM

# --- 3) STM calculations (FDBM STM) ---

declare -A voltage_ranges
voltage_ranges=(
    ["-0.5_0.0"]="-0.5 0.0"
    ["-0.3_0.0"]="-0.3 0.0"
    ["-0.1_0.0"]="-0.1 0.0"
    ["0.0_0.1"]="0.0 0.1"
    ["0.0_0.3"]="0.0 0.3"
    ["0.0_0.5"]="0.0 0.5"
)

mkdir -p filled_states empty_states

# keep original INCAR from LDIPOL
cp INCAR INCAR_orig

for key in "${!voltage_ranges[@]}"; do
    voltages=${voltage_ranges[$key]}

    if [[ $voltages == -* ]]; then
        state_dir="filled_states/${key}"
    else
        state_dir="empty_states/${key}"
    fi

    mkdir -p "$state_dir"

    cp INCAR_orig INCAR
    {
        echo "# STM"
        echo "LPARD = .TRUE."
        echo "LPARDH5 = .TRUE."
        echo "NBMOD = -3"
        echo "EINT = $voltages"
        echo "LSEPB = .FALSE."
        echo "LSEPK = .FALSE."
    } >> INCAR

    srun vasp_std >> vasp.out

    mv vaspout.h5 "$state_dir/"
    mv PARCHG "$state_dir/"
    mv vasprun.xml "$state_dir/"

    rm DOSCAR

    echo "Finished computation for voltage range: $voltages"
done
