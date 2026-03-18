#!/bin/bash
#SBATCH --account=XXX
#SBATCH --partition=YYY #small #gputest, gpusmall
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:15:00
#SBATCH --gres=gpu:a100:1

export PATH="PATH_TO_MACE_ENV:$PATH"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo $OMP_NUM_THREADS

mkdir -p mlip_eval/evaluations/all_substrates/test_mols && cd mlip_eval/evaluations/all_substrates/test_mols

echo "Running MACE-evaluation calcs..."
srun python mlip_eval/scripts/run_MACE_comparison_and_analysis_per_atom_final_static_comp.py mlip_eval/ref_data/all_substrates/test_mols/collected_opt_data_all_subs.extxyz mlip_eval/ref_data/all_substrates/test_mols/collected_opt_data_all_subs.extxyz  # input mace, input reference
echo "Done with calcs!"

echo "plotting results"
cd mlip_eval/evaluations/all_substrates/test_mols/relaxed_MACE_results_final_comparison_static
srun python mlip_eval/scripts/plot_model_analysis_updated_per_atom_pd_summary_static.py
echo "done plotting"


