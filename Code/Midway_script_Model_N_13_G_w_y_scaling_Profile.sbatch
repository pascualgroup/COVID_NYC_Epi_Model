#!/bin/bash
#SBATCH --job-name=G_w_y_scaling_N_13
#SBATCH --output=G_w_y_scaling_N_13_%A_%a.out
#SBATCH --error=error_G_w_y_scaling_N_13_%A_%a.err
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
R CMD BATCH --vanilla '--args  G_w_y_scaling N_13' MIF_run_Profile_Model_N_13.R    OUT_G_w_y_scaling/out.$SLURM_ARRAY_TASK_ID 
