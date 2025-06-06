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
FUNCTION='SBJ08b_HFA_plot_SR_actv'

# set up matlab function call
func_call="${FUNCTION}('${SBJ}', '${an_id_s}', '${an_id_r}', '${actv_win}', '${plt_id}', '${save_fig}', 'off')"

# define commands to execute via SGE
echo ${SBJ}
echo ${func_call}
echo $$
echo ${func_call} > NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}_${an_id_s}_${an_id_r}.m
time matlab -nodesktop -nosplash -nodisplay < NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}_${an_id_s}_${an_id_r}.m
rm NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}_${an_id_s}_${an_id_r}.m
