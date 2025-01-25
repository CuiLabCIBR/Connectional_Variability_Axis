% ---------------------------------------------------------------------------------------------------------------
% This script was used to generate figures and results of Figure 1.
% <Individual variability in edge-level FC declines along a connectional axis>
% The fc variability matrix was plotted at the nodal-level and network-level for HCP-D and HCP-YA.
% The similarity between the two datasets were measured by spearman's rank correlation.
% The within-network and between-network fc variability values were compared using permutation test.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_02_connectional_axis_of_fc_variability/sensitivity_analyses/'];
var_dir = [root_dir 'data/fc_variability/schaefer400/'];

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.
half_flag = 0; % if plot the lower triangle of the matrix, set it to 1.

%% plot the fc variability matrix for hcp-d
fc_var_yen_schaefer400 = load([var_dir 'fc_variability_yen.mat']);
fc_var_yen_schaefer400 = fc_var_yen_schaefer400.fc_variability;

% nodal level
half_flag = 0;
plot_matrix(fc_var_yen_schaefer400,net_label,net_order,half_flag)
clim([0,0.7]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_yen_schaefer400_full.png'])
close all

half_flag = 1;
plot_matrix(fc_var_yen_schaefer400,net_label,net_order,half_flag)
clim([0,0.7]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_yen_schaefer400_half.png'])
close all

% network level
half_flag = 0;
yen_schaefer400_net = plot_matrix_mean(fc_var_yen_schaefer400,net_label,net_order,half_flag);
clim([0.1,0.5]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_yen_schaefer400_net_full.png'])
close all

half_flag = 1;
plot_matrix_mean(fc_var_yen_schaefer400,net_label,net_order,half_flag);
clim([0.1,0.5]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_yen_schaefer400_net_half.png'])
close all

% used in step_04_barplot_connectional_variability_axis_single_net.R
save([var_dir '/yen_schaefer400_net.mat'],'yen_schaefer400_net')

connectional_variability_axis_yen = get_connectional_axis(yen_schaefer400_net);
% used in step_05_barplot_connectional_axis_all_net.R
writecell(connectional_variability_axis_yen,[working_dir 'connectional_variability_axis_yen.xlsx'])

%% permutation test between A-A, S-S and S-A fc variability
outpath = [var_dir, 'yen_'];
% used in step_06_boxplot_S_A_axis.R
[S_S,A_A,S_A] = get_sensorimotor_association_fc_var(fc_var_yen_schaefer400,net_label,net_order,outpath);

M = 10000;
[p,diff_true] = permutation_test(A_A,S_S,M)
[p,diff_true] = permutation_test(A_A,S_A,M)
[p,diff_true] = permutation_test(S_S,S_A,M)

%%
var_temp = yen_schaefer400_net(idx_tril);
var_temp = sort(var_temp,'descend');
mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
slope_yen = abs(table2array(mdl.Coefficients(2,1)));