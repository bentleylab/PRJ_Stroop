#!/bin/bash

#SBATCH --job-name=Stroop_Baseline_Analysis
#SBATCH --output=/data/user/anaskhan/PRJ_Stroop/cluster_logs/logfile_%A_%a.out
#SBATCH --error=/data/user/anaskhan/PRJ_Stroop/cluster_errors/results_%A_%a.err
#SBATCH --partition=short
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200G
#SBATCH --array=0-11

# Load MATLAB module (if required)
module load rc/matlab

# Export additional variables (if needed)
export PROC_ID="main_ft"
export AN_ID="ISPC_S_wvlt_f2to40Log_zbt5to2_trl5to2"
export COND='CI'

# Read subject ID from the file using $SLURM_ARRAY_TASK_ID
SUBJECT=$(awk "NR==${SLURM_ARRAY_TASK_ID}+1" new_sbjs.txt)

# Run BaselineAnalysis
matlab -nodisplay -nosplash -r "SequenceGrangerAnalysis('$SUBJECT','$PROC_ID', '$AN_ID','$COND'); exit"