#!/bin/bash

#SBATCH --job-name=BLPAC
#SBATCH --output=/data/user/anaskhan/PRJ_Stroop/cluster_logs/logfile_%A_%a.out
#SBATCH --error=/data/user/anaskhan/PRJ_Stroop/cluster_errors/results_%A_%a.err
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --array=0-16

# Load MATLAB module (if required)
module load rc/matlab

# Export additional variables (if needed)
export PROC_ID="main_ft"

# Read subject ID from the file using $SLURM_ARRAY_TASK_ID
SUBJECT=$(awk "NR==${SLURM_ARRAY_TASK_ID}+1" subject_ids.txt)

# Run BaselineAnalysis
matlab -nodisplay -nosplash -r "BaselinePACAnalysis('$SUBJECT','$PROC_ID'); exit"