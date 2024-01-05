% ---------------------------------------------------------------------------------------------------------------
% This script is used to get the FC variability matrix.
% The intra-class correlation (ICC) was used as 
% the normalized inter-subject FC variability in this study.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Gradient/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_01_inter_individual_fc_variability/hcpd/'];
var_dir = [root_dir 'data/fc_variability/schaefer400/'];

%% schaefer400
load([working_dir, 'lme_hcpd_schaefer400.mat'])

icc = lme_hcpd_schaefer400.ICC_c;
icc = squareform(icc);
fc_variability_hcpd.schaefer400 = icc;

%%
save([var_dir 'fc_variability_hcpd.mat'],'fc_variability_hcpd')
