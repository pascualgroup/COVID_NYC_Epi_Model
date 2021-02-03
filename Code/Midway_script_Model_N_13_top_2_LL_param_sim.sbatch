#!/bin/bash
#SBATCH --job-name=top_2_LL_sim__N_13
#SBATCH --output=top_2_LL_sim_N_13_%A_%a.out
#SBATCH --error=error_top_2_LL_sim__N_13_%A_%a.err
#SBATCH --array=1-500
#SBATCH --partition=broadwl
#SBATCH --account=covid-19
#SBATCH --qos=covid-19
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000

echo $SLURM_ARRAY_TASK_ID

module load gcc
module load R/3.5.1
R CMD BATCH --vanilla '--args  N_13' Sim_Top_2_LL_Model_N_12_Midway_test.R    OUT_top_2_LL_sim/out.$SLURM_ARRAY_TASK_ID 
