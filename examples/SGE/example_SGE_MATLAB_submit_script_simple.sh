#!/bin/sh
# Edit this to be the name of your MATLAB script
FUNCTION='script_name'

# This assumes your script is in your home directory.
#   Otherwise, change this to the directory with the script.
cd ${HOME}

# Don't edit anything else!
DATASET="${SGE_TASK}"
func_call="${FUNCTION}('${DATASET}')"
echo ${func_call} > ${HOME}/${FUNCTION}_${DATASET}.m
matlab -nodesktop -nosplash -nodisplay < ${HOME}/${FUNCTION}_${DATASET}.m
rm ${HOME}/${FUNCTION}_${DATASET}.m

