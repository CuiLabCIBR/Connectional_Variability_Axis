#!/bin/bash
#SBATCH -p q_cn
#SBATCH --ntasks=1 # Run a single serial task
#SBATCH --cpus-per-task=1
#SBATCH -e job.%j.log # Standard error
#SBATCH --job-name=calc_FC

for subj in {1..275}
do
   sbatch step_01_calculate_FC.sh $subj
done