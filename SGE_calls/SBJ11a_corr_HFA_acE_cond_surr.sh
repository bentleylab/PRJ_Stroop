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
FUNCTION='SBJ11a_corr_HFA_acE_cond_surr'

# set up matlab function call
func_call="${FUNCTION}('${SBJ}', '${pipeline_id}', '${stat_id}', '${an_id}', '${plt_id}')"

# define commands to execute via SGE
echo ${SBJ}
echo ${func_call}
echo $$
echo ${func_call} > NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}_${pipeline_id}_${stat_id}_${an_id}_${plt_id}.m
time matlab -nodesktop -nosplash -nodisplay < NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}_${pipeline_id}_${stat_id}_${an_id}_${plt_id}.m
rm NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}_${pipeline_id}_${stat_id}_${an_id}_${plt_id}.m
