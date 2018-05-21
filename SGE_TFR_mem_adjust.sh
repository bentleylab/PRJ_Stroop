#!/bin/bash
#
# Adjust the memory allocation of SGE jobs for TFR comutations

# appropritate amoun ti product of time and n channels

# get output of qstat

# grab jobs with TFRstat names

# for each job:
# run qstat -j ${jobid}
# grab SBJ_ID
# qalt to change memory allocation of that job to that SBJ_ID
