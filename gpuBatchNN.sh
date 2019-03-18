#!/bin/bash -l

# set the name of the job; this will appear in the job listing
#SBATCH --job-name=Polar

# set the maximum memory usage (per slot)
#SBATCH --mem=8G

#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu

#SBATCH -o polar."%j".out
#SBATCH -e polar."%j".err
# Default in slurm
#SBATCH --mail-user=$USER@stanford.edu
#SBATCH --mail-type=ALL
# Request 5 hours run time
#SBATCH -t 48:0:0

echo "Start test"
module load matlab 

matlab -nodesktop -r "run('NeuralPolarDecode.m'); exit(0);"

echo "End test"

