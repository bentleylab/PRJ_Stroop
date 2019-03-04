#!/bin/sh
# This script analyzes subject data in ~/MATLAB/DATA/Avgusta for each task
# The SGE_TASK variable specifies the data sets, as in $HOME/test.qsub

cd /home/knight/hoycw/PRJ_Stroop/scripts/

# define where the datasets are located
#PRCSDATADIR="/home/knight/hoycw/PRJ_Stroop/data/"

# don't change this variable
# used by the submit script to define which data sets to analyze
SBJ="${SGE_TASK}"

# define function
FUNCTION='check_neuralynx_validsamples'

# set up matlab function call
func_call="${FUNCTION}('${SBJ}')"

# define commands to execute via SGE
echo ${DATASET}
echo ${func_call}
echo $$
echo ${func_call} > NotBackedUp/tmpSGE/check_neuralynx_validsamples_${SBJ}.m
time matlab -nodesktop -nosplash -nodisplay < NotBackedUp/tmpSGE/check_neuralynx_validsamples_${SBJ}.m
rm NotBackedUp/tmpSGE/check_neuralynx_validsamples_${SBJ}.m
