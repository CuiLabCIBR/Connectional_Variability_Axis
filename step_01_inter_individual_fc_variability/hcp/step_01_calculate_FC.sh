#!/bin/bash
#SBATCH -p q_cn
#SBATCH --ntasks=1 # Run a single serial task
#SBATCH --cpus-per-task=1
#SBATCH -e job.%j.log # Standard error
#SBATCH --job-name=calc_FC

module load MATLAB/R2019a
cd /ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/hcp/
subi=$1 
matlab -singleCompThread -nodisplay -nosplash -r "step_01_calculate_fc_schaefer400($subi)"