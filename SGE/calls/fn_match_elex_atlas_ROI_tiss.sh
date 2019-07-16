#!/bin/sh
# The SGE_TASK variable specifies the data sets, as in $HOME/test.qsub

cd /home/knight/hoycw/PRJ_Stroop/scripts/utils/

# don't change this variable, used by the submit script to define which data sets to analyze
SBJ="${SGE_TASK}"

# define function
FUNCTION='fn_match_elec_atlas_ROI_tiss'

# set up matlab function call
func_call="${FUNCTION}('${SBJ}','${proc_id}','${view_space}','${reg_type}','${atlas_id}','${reref}')"

# define commands to execute via SGE
echo ${DATASET}
echo ${func_call}
echo $$
echo ${func_call} > ../NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}.m
time matlab -nodesktop -nosplash -nodisplay < ../NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}.m
rm ../NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}.m
