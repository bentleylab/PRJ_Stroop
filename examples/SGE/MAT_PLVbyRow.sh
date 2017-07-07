#!/bin/sh
# This script analyzes subject data in ~/MATLAB/DATA/Avgusta for each task
# The SGE_TASK variable specifies the data sets, as in $HOME/test.qsub

cd /home/despo/hoycw/KnightLab/PRJ_DecisionMaking/Scripts/

# define where the datasets are located
PRCSDATADIR="/home/despo/hoycw/KnightLab/PRJ_DecisionMaking/PrcsData/"

# don't change this variable
# used by the submit script to define which data sets to analyze
DATASET="${SGE_TASK}"

# define function
FUNCTION='func_PLV_confusion_row'

# set up matlab function call
func_call="${FUNCTION}('${PRCSDATADIR}', '${DATASET}', '${TIMING}', '${PERIOD}', '${LOFREQ}', '${HIFREQ}')"

# define commands to execute via SGE
echo ${DATASET}
echo ${func_call}
echo ${func_call} > MAT_PLVbyRow_${DATASET}.m
matlab -nodesktop -nosplash -nojvm -nodisplay < MAT_PLVbyRow_${DATASET}.m
rm MAT_PLVbyRow_${DATASET}.m
