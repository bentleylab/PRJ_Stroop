#!/bin/sh
# The SGE_TASK variable specifies the data sets, as in $HOME/test.qsub

cd /home/knight/hoycw/PRJ_Stroop/scripts/utils/

# don't change this variable, used by the submit script to define which data sets to analyze
SBJ="${SGE_TASK}"

# define function
FUNCTION='fn_save_elec_atlas'

# set up matlab function call
func_call="${FUNCTION}('${SBJ}','${pipeline_id}','${view_space}','${reg_type}','${atlas_id}','${reref}')"

# define commands to execute via SGE
echo ${DATASET}
echo ${func_call}
echo $$
echo ${func_call} > ../NotBackedUp/tmpSGE/fn_save_elec_atlas_${SBJ}.m
time matlab -nodesktop -nosplash -nodisplay < ../NotBackedUp/tmpSGE/fn_save_elec_atlas_${SBJ}.m
rm ../NotBackedUp/tmpSGE/fn_save_elec_atlas_${SBJ}.m
