#!/bin/sh
# ^-- this "shebang" indicates this is a bash (shell) script
# Also, ${var} is how shell scripts call variables.

# For general documentation on all SGE functions on the HWNI cluster:
#   https://www.neuro.berkeley.edu/resources/software/sge/scripting.html
# For documentation on commands to check the status of your jobs:
#   https://www.neuro.berkeley.edu/resources/software/sge/commands.html
# For documentation on reserving memory:
#   https://www.neuro.berkeley.edu/resources/software/sge/memory.html

# This script runs one MATLAB script for a list of inputs (patients) 
# The SGE_TASK variable specifies the string inputs
# Here, the inputs will be a list of patient IDs

# Edit this line to move inside the job to wherever your .m script is
cd /home/knight/hoycw/PRJ_Stroop/scripts/

# Don't change this variable! (SGE_TASK is the specific name
# used by the -f flag of the submit command on the cluster)
#   (Actually, you could change the "DATASET" name of the variable,
#   the ${SGE_TASK} part is the one that SGE uses.)
DATASET="${SGE_TASK}"

# Edit this to be the name of your MATLAB script
# The .m script should be written as a function,
# so that the func_call string could be run in the
# command window of MATLAB once the strings are substituted
FUNCTION='script_name'

# This string will be run in MATLAB
func_call="${FUNCTION}('${DATASET}')"
# If you had other inputs to your function call, you could
# edit the line like the example below to use the options file
# to specify the other parameters.
#   func_call="${FUNCTION}('${DATASET}', '${option_var1}', '${option_var2}')"

# These echo commands print text to the .o output file,
# so it's nice to print what input and the exact script
# you're running in any particular job
echo ${DATASET}
echo ${func_call}
# This one prints the job ID, which is helpful for checking progress
echo $$

# This creates a temporary .m script with only the single line containing func_call
# I have a directory (tmpSGE/) where these temporary scripts live while the job runs,
# but they could be in your home or scripts directory if you wanted.
#   Side note: Julie's book keeping scripts that create daily/weekly/monthly backups
#   will skip any files within/underneath any directories called "NotBackedUp/" 
#   It's not necessary to put your temporary SGE scripts or the SGE text/error output files in a 
#   "NotBackedUp" folder, but I figured it'd be worth mentioning as a cool feature of the cluster.
echo ${func_call} > NotBackedUp/tmpSGE/${FUNCTION}_${DATASET}.m
# WARNING: Be careful not to create multiple temporary files with the same name!
#   This could happen if you used 2 values for ${option_var1} and ran them for the same
#   patient at the same time, so then maybe name the temporary script ${FUNCTION}_${DATASET}_${option_var1}.m

# This will start a MATLAB session without any graphics, just a command line, and
# will run the ${func_call} line in our temporary script on the MATLAB command line.
#   Side note: the "time" part at the beginning is a shell command that will spit out
#   how long it took for the following command to run, so again, not necessary for 
#   running SGE jobs, but a nice trick I thought I'd share.
time matlab -nodesktop -nosplash -nodisplay < NotBackedUp/tmpSGE/${FUNCTION}_${DATASET}.m

# This will delete that temporary script we made after the job is done.
# Note that if this shell script gets killed in the job before this line is run,
# that temporary script might still exist and could help you see what happened.
rm NotBackedUp/tmpSGE/${FUNCTION}_${DATASET}.m

