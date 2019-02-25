#!/bin/bash -l

# set the name of the job; this will appear in the job listing
#SBATCH --job-name=PlotBERBaseline

# set the maximum memory usage (per slot)
#SBATCH --mem=2G

##SBATCH --partition=gpu
##SBATCH --gres=gpu:1

#SBATCH -o test."%j".out
#SBATCH -e test."%j".err
# Default in slurm
#SBATCH --mail-user $USER@stanford.edu
#SBATCH --mail-type=ALL
# Request 5 hours run time
#SBATCH -t 5:0:0
#SBATCH -p normal


echo "Start test"
module load matlab 

matlab -nodesktop -r "run('PlotPerfectBaseline.m'); run('PlotEstimateBaseline.m'); exit(0);"

echo "End test"

