clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_03_structural_connectome_variability/sensitivity_analyses/'];
var_dir = [root_dir 'data/fc_variability/schaefer400/'];
mat_dir = [root_dir 'data/connectome_matrix/schaefer400/'];

load('7net_label_schaefer400.mat')
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.
idx_nolimbic = find(net_label ~= 5); % We excluded the limbic network in this study.

mask = zeros(400,400);
mask(idx_nolimbic, idx_nolimbic) = 1;

corr_method = 'spearman';
load([root_dir 'data/parcellation_files/perm_id_schaefer400.mat'])

%% Initialization inputs and outputs
load([mat_dir 'sc_variability_hcpd.mat'],'sc_variability_hcpd')
sc_var_hcpd = log10(sc_variability_hcpd);

load([mat_dir 'sc_variability_hcp.mat'],'sc_variability_hcp')
sc_var_hcp = log10(sc_variability_hcp);

%% noGSR
% HCP-D
fc_var_hcpd_noGSR = load([var_dir 'fc_variability_hcpd_noGSR.mat']);
fc_var_hcpd_noGSR = fc_var_hcpd_noGSR.fc_variability;

[r_hcpd,p_hcpd,hcpd_fc_var_sc_var,~,p_hcpd_spin] = corr_matrix(fc_var_hcpd_noGSR,sc_var_hcpd,net_label,net_order,corr_method,perm_id);
writetable(hcpd_fc_var_sc_var,[working_dir '/hcpd_fc_var_sc_var_noGSR.csv'])

% HCP-YA
fc_var_hcp_noGSR = load([var_dir 'fc_variability_hcp_noGSR.mat']);
fc_var_hcp_noGSR = fc_var_hcp_noGSR.fc_variability;

[r_hcp,p_hcp,hcp_fc_var_sc_var,~,p_hcp_spin] = corr_matrix(fc_var_hcp_noGSR,sc_var_hcp,net_label,net_order,corr_method,perm_id);
writetable(hcp_fc_var_sc_var,[working_dir '/hcp_fc_var_sc_var_noGSR.csv'])

%% save results
corr_hcpd = cell(2,4);
corr_hcpd(1,:) = {'hcp-d','r-value','p-value','p-spin'};
corr_hcpd(2,1) = {'noGSR'};
corr_hcpd(2,2:end) = num2cell([r_hcpd, p_hcpd, p_hcpd_spin]);

corr_hcp = corr_hcpd;
corr_hcp(1,1) = {'hcp-ya'};
corr_hcp(2,2:end) = num2cell([r_hcp, p_hcp, p_hcp_spin]);

corr_results = [corr_hcpd;corr_hcp];
writecell(corr_results,[working_dir 'fc_var_sc_var_corr_noGSR.csv'])
