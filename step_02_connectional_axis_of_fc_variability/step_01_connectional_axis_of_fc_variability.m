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

working_dir = [root_dir 'step_02_connectional_axis_of_fc_variability/'];
var_dir = [root_dir 'data/fc_variability/schaefer400/'];

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.
half_flag = 0; % if plot the lower triangle of the matrix, set it to 1.

%% plot the fc variability matrix for hcp-d
fc_var_hcpd_schaefer400 = load([var_dir 'fc_variability_hcpd.mat']);
fc_var_hcpd_schaefer400 = fc_var_hcpd_schaefer400.fc_variability;

% nodal level
half_flag = 0;
plot_matrix(fc_var_hcpd_schaefer400,net_label,net_order,half_flag)
clim([0,0.5]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_hcpd_schaefer400_full.png'])
close all

half_flag = 1;
plot_matrix(fc_var_hcpd_schaefer400,net_label,net_order,half_flag)
clim([0,0.5]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_hcpd_schaefer400_half.png'])
close all

% network level
half_flag = 0;
hcpd_schaefer400_net = plot_matrix_mean(fc_var_hcpd_schaefer400,net_label,net_order,half_flag);
clim([0.05,0.35]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_hcpd_schaefer400_net_full.png'])
close all

half_flag = 1;
plot_matrix_mean(fc_var_hcpd_schaefer400,net_label,net_order,half_flag);
clim([0.05,0.35]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_hcpd_schaefer400_net_half.png'])
close all

% used in step_04_barplot_connectional_variability_axis_single_net.R
save([var_dir '/hcpd_schaefer400_net.mat'],'hcpd_schaefer400_net')

connectional_variability_axis_hcpd = get_connectional_axis(hcpd_schaefer400_net);
% used in step_05_barplot_connectional_axis_all_net.R
writecell(connectional_variability_axis_hcpd,[working_dir 'connectional_variability_axis_hcpd.xlsx'])

%% plot the fc variability matrix for hcp-ya
fc_var_hcp_schaefer400 = load([var_dir 'fc_variability_hcp.mat']);
fc_var_hcp_schaefer400 = fc_var_hcp_schaefer400.fc_variability;

% nodal level
half_flag = 0;
plot_matrix(fc_var_hcp_schaefer400,net_label,net_order,half_flag)
clim([0,0.6]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_hcp_schaefer400_full.png'])
close all

half_flag = 1;
plot_matrix(fc_var_hcp_schaefer400,net_label,net_order,half_flag)
clim([0,0.6]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_hcp_schaefer400_half.png'])
close all

% network level
half_flag = 0;
hcp_schaefer400_net = plot_matrix_mean(fc_var_hcp_schaefer400,net_label,net_order,half_flag);
clim([0.10,0.40]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_hcp_schaefer400_net_full.png'])
close all

half_flag = 1;
plot_matrix_mean(fc_var_hcp_schaefer400,net_label,net_order,half_flag);
clim([0.10,0.40]);
print(gcf,'-dpng','-r300',[working_dir '/matrix_plot_hcp_schaefer400_net_half.png'])
close all

% used in step_04_barplot_connectional_variability_axis_single_net.R
save([var_dir '/hcp_schaefer400_net.mat'],'hcp_schaefer400_net')

connectional_variability_axis_hcp = get_connectional_axis(hcp_schaefer400_net);
% used in step_05_barplot_connectional_variability_axis_all_net.R
writecell(connectional_variability_axis_hcp,[working_dir 'connectional_variability_axis_hcp.xlsx'])

%% calculate the similarity between hcp-d and hcp-ya
corr_method = 'spearman';
load([root_dir 'data/parcellation_files/perm_id_schaefer400.mat'])

% nodal-level similarity
[r_node,p_node,hcpd_hcp_schaefer400,r_spin,p_spin] = corr_matrix(fc_var_hcpd_schaefer400,fc_var_hcp_schaefer400,net_label,net_order,corr_method,perm_id);

writetable(hcpd_hcp_schaefer400,[working_dir '/hcpd_hcp_fc_variability_schaefer400.csv'])

% network-level similarity
idx_tril = find(tril(ones(6,6)) > 0);
[r_net,p_net] = corr(hcpd_schaefer400_net(idx_tril),hcp_schaefer400_net(idx_tril),'type',corr_method);

hcpd_hcp_schaefer400_net.mat_a = hcpd_schaefer400_net(idx_tril);
hcpd_hcp_schaefer400_net.mat_b = hcp_schaefer400_net(idx_tril);
hcpd_hcp_schaefer400_net = struct2table(hcpd_hcp_schaefer400_net);

writetable(hcpd_hcp_schaefer400_net,[working_dir 'hcpd_hcp_fc_variability_schaefer400_net.csv'])

%% permutation test between A-A, S-S and S-A fc variability
outpath = [var_dir, 'hcpd_'];
% used in step_06_boxplot_S_A_axis.R
[S_S,A_A,S_A] = get_sensorimotor_association_fc_var(fc_var_hcpd_schaefer400,net_label,net_order,outpath);

M = 10000;
[p,diff_true] = permutation_test(A_A,S_S,M)
[p,diff_true] = permutation_test(A_A,S_A,M)
[p,diff_true] = permutation_test(S_S,S_A,M)

outpath = [var_dir, 'hcp_'];
% used in step_06_boxplot_S_A_axis.R
[S_S,A_A,S_A] = get_sensorimotor_association_fc_var(fc_var_hcp_schaefer400,net_label,net_order,outpath);

[p,diff_true] = permutation_test(A_A,S_S,M)
[p,diff_true] = permutation_test(A_A,S_A,M)
[p,diff_true] = permutation_test(S_S,S_A,M)

%% connectional axis slope
var_temp = hcpd_schaefer400_net(idx_tril);
var_temp = sort(var_temp,'descend');
mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
slope_hcpd = abs(table2array(mdl.Coefficients(2,1)));

var_temp = hcp_schaefer400_net(idx_tril);
var_temp = sort(var_temp,'descend');
mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
slope_hcp = abs(table2array(mdl.Coefficients(2,1)));
