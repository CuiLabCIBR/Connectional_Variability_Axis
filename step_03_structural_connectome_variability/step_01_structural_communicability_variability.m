% ---------------------------------------------------------------------------------------------------------------
% This script was used to generate figures and results of Figure 2.
% <Individual variability in structural connectivity communicability is associated with 
% the connectional axis pattern in FC variability across connectome edges>
% Communicability (CMY) was used to measure the indirect structural connections.
% Mean absolute deviation (MAD) was used to measure the CMY variability.
% The similarity between fc variability matrix and CMY variability matrix was measured by spearman's rank correlation.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_03_structural_connectome_variability/'];
var_dir = [root_dir 'data/fc_variability/schaefer400/'];
mat_dir = [root_dir 'data/connectome_matrix/schaefer400/'];

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.
idx_nolimbic = find(net_label ~= 5);

half_flag = 1; % if plot the lower triangle of the matrix, set it to 1.

corr_method = 'spearman';
load([root_dir 'data/parcellation_files/perm_id_schaefer400.mat'])

consistency_thr = 0.75;

%% HCP-D
load([mat_dir 'sc_hcpd.mat'])
[~,~,n_hcpd] = size(sc_all_hcpd);
sc_all_hcpd = sc_all_hcpd(idx_nolimbic, idx_nolimbic, :);

W_thr = threshold_consistency(sc_all_hcpd, consistency_thr);
sc_mask_hcpd = double(W_thr > 0);

% plot the sc matrix at nodal level
for i = 1:n_hcpd
    sc_all_hcpd_mask(:,:,i) = sc_all_hcpd(:,:,i) .* sc_mask_hcpd;
end
sc_hcpd_group_atlas_mask = mean(sc_all_hcpd_mask,3);
sc_hcpd = zeros(400,400);
sc_hcpd(idx_nolimbic, idx_nolimbic) = log10(sc_hcpd_group_atlas_mask);

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
sc_communicability_hcpd = zeros(400, 400);
sc_communicability_hcpd(idx_nolimbic, idx_nolimbic) = G_mean_mat_hcpd;
save([mat_dir 'sc_communicability_hcpd.mat'],'sc_communicability_hcpd')

% plot mean communicability full
plot_matrix(log10(sc_communicability_hcpd),net_label,net_order,half_flag)
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_cmy_hcpd_full.png'])
close all

% mad communicability
G_mad_vec_hcpd = mad(G_vec_hcpd')';
G_mad_mat_hcpd = squareform(G_mad_vec_hcpd);
sc_variability_hcpd = zeros(400, 400);
sc_variability_hcpd(idx_nolimbic, idx_nolimbic) = G_mad_mat_hcpd;
save([mat_dir 'sc_variability_hcpd.mat'],'sc_variability_hcpd')

% plot mad communicability half
plot_matrix(log10(sc_variability_hcpd),net_label,net_order,1)
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_sc_variability_hcpd_half.png'])
close all

% plot network-level mad communicability half
plot_sc_mean(sc_variability_hcpd,net_label,net_order,1);
caxis([-4.5,-2.5])
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_sc_variability_net_hcpd_half.png'])
close all

% used in step_02_density_plot.py
load([mat_dir 'sc_variability_hcpd.mat'],'sc_variability_hcpd')

fc_variability_hcpd = load([var_dir 'fc_variability_hcpd.mat']);
fc_variability_hcpd = fc_variability_hcpd.fc_variability;

[r_hcpd,p_hcpd,hcpd_fc_variability_sc_variability,r_hcpd_spin,p_hcpd_spin] = corr_matrix(fc_variability_hcpd,log10(sc_variability_hcpd),net_label,net_order,corr_method,perm_id);
writetable(hcpd_fc_variability_sc_variability,[working_dir '/hcpd_fc_variability_sc_variability.csv']) 

%% HCP-YA
load([mat_dir 'sc_hcp.mat'])
[~,~,n_hcp] = size(sc_all_hcp);
sc_all_hcp = sc_all_hcp(idx_nolimbic, idx_nolimbic, :);
W_thr = threshold_consistency(sc_all_hcp, consistency_thr);
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
sc_communicability_hcp = zeros(400, 400);
sc_communicability_hcp(idx_nolimbic, idx_nolimbic) = G_mean_mat_hcp;
save([mat_dir 'sc_communicability_hcp.mat'],'sc_communicability_hcp')

% mad communicability
G_mad_vec_hcp = mad(G_vec_hcp')';
G_mad_mat_hcp = squareform(G_mad_vec_hcp);
sc_variability_hcp = zeros(400, 400);
sc_variability_hcp(idx_nolimbic, idx_nolimbic) = G_mad_mat_hcp;
save([mat_dir 'sc_variability_hcp.mat'],'sc_variability_hcp')

% plot mad communicability half
plot_matrix(log10(sc_variability_hcp),net_label,net_order,1)
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_sc_variability_hcp_half.png'])
close all

% plot network-level mad communicability half
plot_sc_mean(sc_variability_hcp,net_label,net_order,1);
caxis([-4.5,-2.5])
print(gcf,'-dpng','-r300',[working_dir 'matrix_plot_sc_variability_net_hcp_half.png'])
close all

% used in step_02_density_plot.py
load([mat_dir 'sc_variability_hcp.mat'],'sc_variability_hcp')

fc_variability_hcp = load([var_dir 'fc_variability_hcp.mat']);
fc_variability_hcp = fc_variability_hcp.fc_variability;

[r_hcp,p_hcp,hcp_fc_variability_sc_variability,r_hcp_spin,p_hcp_spin] = corr_matrix(fc_variability_hcp,log10(sc_variability_hcp),net_label,net_order,corr_method,perm_id);
writetable(hcp_fc_variability_sc_variability,[working_dir '/hcp_fc_variability_sc_variability.csv']) 

%% save results
corr_hcpd = cell(2,4);
corr_hcpd(1,:) = {'hcp-d','r-value','p-value','p-spin'};
corr_hcpd(2,1) = {'schaefer400'};
corr_hcpd(2,2:end) = num2cell([r_hcpd, p_hcpd, p_hcpd_spin]);

corr_hcp = corr_hcpd;
corr_hcp(1,1) = {'hcp-ya'};
corr_hcp(2,2:end) = num2cell([r_hcp, p_hcp, p_hcp_spin]);

corr_results = [corr_hcpd;corr_hcp];
writecell(corr_results,[working_dir 'fc_var_sc_var_corr_schaefer400.csv'])

