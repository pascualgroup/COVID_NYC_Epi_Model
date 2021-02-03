#!/bin/bash
#SBATCH --job-name=b_a_N_12
#SBATCH --output=b_a_N_12_%A_%a.out
#SBATCH --error=error_b_a_N_12_%A_%a.err
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
R CMD BATCH --vanilla '--args  b_a N_12' MIF_run_Profile_Model_N_12.R    OUT_b_a/out.$SLURM_ARRAY_TASK_ID 
