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
FUNCTION='SBJ08d_HFA_ROL'

# set up matlab function call
func_call="${FUNCTION}('${SBJ}', '${an_id}', '${rol_id}', '${plot_qa_trl}', '${plot_qa_summary}', '${plot_stack}', 'off', '${fig_ftype}')"

# define commands to execute via SGE
echo ${SBJ}
echo ${func_call}
echo $$
echo ${func_call} > NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}_${an_id}_${rol_id}.m
time matlab -nodesktop -nosplash -nodisplay < NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}_${an_id}_${rol_id}.m
rm NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}_${an_id}_${rol_id}.m
