#!/bin/bash
#SBATCH -p q_fat_c
#SBATCH -q high_c
#SBATCH --ntasks=1 # Run a single serial task
#SBATCH --cpus-per-task=3
#SBATCH -e job.%j.log # Standard error
#SBATCH --job-name=calc_FC

module load MATLAB/R2019a
cd /ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/hcpd/
subi=$1
matlab -singleCompThread -nodisplay -nosplash -r "step_01_calculate_fc_hcpd('$subi')" 
