 % ---------------------------------------------------------------------------------------------------------------
% This script was used to generate figures and results of Figure 2.
% <Individual variability in structural connectivity communicability is associated with
% the variability axis in functional connectivity across connectome edges>
% Communicability (CMY) was used to measure the indirect structural connections.
% Mean absolute deviation (MAD) was used to measure the CMY variability.
% The similarity between fc variability matrix and CMY variability matrix was measured by spearman's rank correlation.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_03_white_matter_structural_connectome/'];
data_dir = [root_dir 'data/fc_variability/schaefer400/'];
mat_dir = [root_dir 'data/connectome_matrix/schaefer400/'];

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.
half_flag = 1; % if plot the lower triangle of the matrix, set it to 1.

corr_method = 'spearman';

load([data_dir 'fc_variability_hcpd.mat']);
load([data_dir 'fc_variability_hcp.mat']);

%% HCP-D
load([mat_dir 'sc_hcpd.mat'])
[~,~,n_hcpd] = size(sc_all_hcpd);
W_thr = threshold_consistency(sc_all_hcpd, 0.75);
sc_mask_hcpd = double(W_thr > 0);

% plot the sc matrix at nodal level
for i = 1:n_hcpd
    sc_all_hcpd_mask = sc_all_hcpd(:,:,i) .* sc_mask_hcpd;
end
sc_hcpd_group_atlas_mask = mean(sc_all_hcpd_mask,3);
sc_hcpd = log10(sc_hcpd_group_atlas_mask);

half_flag = 0;
plot_matrix(sc_hcpd,net_label,net_order,half_flag)
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_sc_network_hcpd_full.png'])
close all

% calculate communicability for each individual
for i = 1:n_hcpd
    sc_temp = sc_all_hcpd(:,:,i) .* sc_mask_hcpd;
    G_hcpd(:,:,i) = get_communicability(sc_temp,1);
    G_vec_hcpd(:,i) = mat2vec(G_hcpd(:,:,i));
end

% mean communicability (cmy)
G_mean_vec_hcpd = mean(G_vec_hcpd,2);
G_mean_mat_hcpd = squareform(G_mean_vec_hcpd);
sc_communicability_hcpd = G_mean_mat_hcpd;
save([mat_dir 'sc_communicability_hcpd.mat'],'sc_communicability_hcpd')

% plot mean communicability full
plot_matrix(log10(G_mean_mat_hcpd),net_label,net_order,half_flag)
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_cmy_hcpd_full.png'])
close all

% mad communicability
G_mad_vec_hcpd = mad(G_vec_hcpd')';
G_mad_mat_hcpd = squareform(G_mad_vec_hcpd);
sc_variability_hcpd = G_mad_mat_hcpd;
save([mat_dir 'sc_variability_hcpd.mat'],'sc_variability_hcpd')

% plot mad communicability half
plot_matrix(log10(G_mad_mat_hcpd),net_label,net_order,1)
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_cmy_mad_hcpd_half.png'])
close all

% plot network-level mad communicability half
plot_sc_mean(G_mad_mat_hcpd,net_label,net_order,1);
caxis([-4.5,-2.5])
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_cmy_mad_net_hcpd_half.png'])
close all

% used in step_02_density_plot.py
load([mat_dir 'sc_variability_hcpd.mat'],'sc_variability_hcpd')

[r_hcpd,p_hcpd,hcpd_fc_variability_sc_variability] = corr_matrix(fc_variability_hcpd.schaefer400,log10(sc_variability_hcpd),net_label,net_order,corr_method);
writetable(hcpd_fc_variability_sc_variability,[working_dir '/hcpd_fc_variability_sc_variability.csv']) 

%% HCP-YA
load([mat_dir 'sc_hcp.mat'])
[~,~,n_hcp] = size(sc_all_hcp);
W_thr = threshold_consistency(sc_all_hcp, 0.75);
sc_mask_hcp = double(W_thr > 0);

% calculate communicability for each individual
for i = 1:n_hcp
    sc_temp = sc_all_hcp(:,:,i) .* sc_mask_hcp;
    G_hcp(:,:,i) = get_communicability(sc_temp,1);
    G_vec_hcp(:,i) = mat2vec(G_hcp(:,:,i));
end

% mean communicability (cmy)
G_mean_vec_hcp = mean(G_vec_hcp,2);
G_mean_mat_hcp = squareform(G_mean_vec_hcp);
sc_communicability_hcp = G_mean_mat_hcp;
save([mat_dir 'sc_communicability_hcp.mat'],'sc_communicability_hcp')

% mad communicability
G_mad_vec_hcp = mad(G_vec_hcp')';
G_mad_mat_hcp = squareform(G_mad_vec_hcp);
sc_variability_hcp = G_mad_mat_hcp;
save([mat_dir 'sc_variability_hcp.mat'],'sc_variability_hcp')

% plot mad communicability half
plot_matrix(log10(G_mad_mat_hcp),net_label,net_order,1)
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_cmy_mad_hcp_half.png'])
close all

% plot network-level mad communicability half
plot_sc_mean(G_mad_mat_hcp,net_label,net_order,1);
caxis([-4.5,-2.5])
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_cmy_mad_net_hcp_half.png'])
close all

% used in step_02_density_plot.py
load([mat_dir 'sc_variability_hcp.mat'],'sc_variability_hcp')

[r_hcp,p_hcp,hcp_fc_variability_sc_variability] = corr_matrix(fc_variability_hcp.schaefer400,log10(sc_variability_hcp),net_label,net_order,corr_method);
writetable(hcp_fc_variability_sc_variability,[working_dir '/hcp_fc_variability_sc_variability.csv']) 
