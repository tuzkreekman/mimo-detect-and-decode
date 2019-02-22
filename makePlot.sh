#!/bin/bash

#
# set the name of the job; this will appear in the job listing
#$ -N PlotBERBaseline
#

#
# set the maximum memory usage (per slot)
#$ -l mem_free=2G
#
# on other clusters this memory resource may have a different name

#
# set the number of slots, replace '1' with a larger number if needed
#$ -pe shm 1
#
# on other clusters this pe may have a different name

#
# set the maximum run time, hh:mm:ss, default is 48hrs on FarmShare
#$ -l h_rt=12:00:00
#

#
# send mail when job ends or aborts
#$ -m ea
#

#
# specify an email address
#$ -M $USER@stanford.edu
#

# check for errors in the job submission options
#$ -w e
#

##We strongly discourage users from exporting their environment onto the compute node. 
##Doing this pretty much means the job is non-reproductible, 
##because all the required settings are not captured in the job script.
##
## pass the current environment variables
##$ -V
##

# join the stdout and stderr streams into one file
#$ -j y
#

echo "Start test"

matlab -nodesktop -r "run('PlotBaseline.m'); exit(0);"

echo "End test"

