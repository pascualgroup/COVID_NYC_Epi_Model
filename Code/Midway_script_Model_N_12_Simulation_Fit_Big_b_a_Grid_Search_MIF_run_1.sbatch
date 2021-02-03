#!/bin/bash
#SBATCH --job-name=Sim_Fit_big_b_a_Grid_Search_MIF_run_1_N_12
#SBATCH --output=Sim_Fit_big_b_a_Grid_Search_MIF_run_1_N_12_%A_%a.out
#SBATCH --error=error_Sim_Fit_big_b_a_Grid_Search_MIF_run_1_N_12_%A_%a.err
#SBATCH --array=1-500
#SBATCH --partition=broadwl
#SBATCH --account=covid-19
#SBATCH --qos=covid-19
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=2000
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000

echo $SLURM_ARRAY_TASK_ID

module load gcc
module load R/3.5.1
R CMD BATCH --vanilla '--args  N_12' MIF_run_Model_N_12_sim_data_big_b_a.R   OUT_Sim_Fit_Big_b_a_Grid_Search_MIF_run_1/out.$SLURM_ARRAY_TASK_ID 
