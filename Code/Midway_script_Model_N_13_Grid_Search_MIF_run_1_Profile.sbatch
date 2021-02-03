#!/bin/bash
#SBATCH --job-name=Grid_Search_MIF_run_1_N_13
#SBATCH --output=Grid_Search_MIF_run_1_N_13_%A_%a.out
#SBATCH --error=error_Grid_Search_MIF_run_1_N_13_%A_%a.err
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
R CMD BATCH --vanilla '--args  N_13' MIF_run_Model_N_13.R   OUT_Grid_Search_MIF_run_1_N_13/out.$SLURM_ARRAY_TASK_ID 
