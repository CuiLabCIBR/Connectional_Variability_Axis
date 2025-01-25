#!/bin/bash
#SBATCH -p q_fat_c
#SBATCH -q high_c
#SBATCH --ntasks=1 # Run a single serial task
#SBATCH --cpus-per-task=3
#SBATCH --job-name=calc_FC

module load MATLAB/R2019a
cd /ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/hcp/
subj=$1
echo $subj
matlab -singleCompThread -nodisplay -nosplash -r "step_01_calculate_fc_hcp('$subj')"
