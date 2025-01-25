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
fc_var_hcpd = load([var_dir 'fc_variability_hcpd.mat']);
fc_var_hcpd = fc_var_hcpd.fc_variability;

fc_var_hcp = load([var_dir 'fc_variability_hcp.mat']);
fc_var_hcp = fc_var_hcp.fc_variability;

load([mat_dir 'sc_variability_hcpd.mat'],'sc_variability_hcpd')
sc_var_hcpd = log10(sc_variability_hcpd);

load([mat_dir 'sc_variability_hcp.mat'],'sc_variability_hcp')
sc_var_hcp = log10(sc_variability_hcp);

corr_hcpd = cell(3,4);
corr_hcpd(1,:) = {'hcp-d','r-value','p-value','p-spin'};
corr_hcpd(2:end,1) = {'reg_dist','reg_strength'};

corr_hcp = corr_hcpd;
corr_hcp(1,1) = {'hcp-ya'};

r_mat = zeros(2,2);
p_mat = zeros(2,4);

%% regress distance
dist_dir = [root_dir 'data/parcellation_files/'];
load([dist_dir 'dist_schaefer400.mat'],'dist_schaefer400')

% HCP-D
fc_var_hcpd_reg_dist = reg_cov_from_mat(fc_var_hcpd,dist_schaefer400,mask);
sc_var_hcpd_reg_dist = reg_cov_from_mat(sc_var_hcpd,dist_schaefer400,mask);

[r_mat(1,1),p_mat(1,1),hcpd_fc_var_sc_var,~,p_mat(1,2)] = corr_matrix(fc_var_hcpd_reg_dist,sc_var_hcpd_reg_dist,net_label,net_order,corr_method,perm_id);
writetable(hcpd_fc_var_sc_var,[working_dir '/hcpd_fc_var_sc_var_reg_dist.csv'])

% HCP-YA
fc_var_hcp_reg_dist = reg_cov_from_mat(fc_var_hcp,dist_schaefer400,mask);
sc_var_hcp_reg_dist = reg_cov_from_mat(sc_var_hcp,dist_schaefer400,mask);

[r_mat(1,2),p_mat(1,3),hcp_fc_var_sc_var,~,p_mat(1,4)] = corr_matrix(fc_var_hcp_reg_dist,sc_var_hcp_reg_dist,net_label,net_order,corr_method,perm_id);
writetable(hcp_fc_var_sc_var,[working_dir '/hcp_fc_var_sc_var_reg_dist.csv'])

%% regress strength
% HCP-D
fc_hcpd = load([mat_dir 'fc_hcpd.mat']);
fc_hcpd = fc_hcpd.fc_mat;

load([mat_dir 'sc_communicability_hcpd.mat'])
cov_hcpd(:,:,1) = fc_hcpd;
cov_hcpd(:,:,2) = log10(sc_communicability_hcpd);

sc_var_min = min(mat2vec(sc_var_hcpd(idx_nolimbic, idx_nolimbic)));
sc_var_max = max(mat2vec(sc_var_hcpd(idx_nolimbic, idx_nolimbic)));

fc_var_hcpd_reg_strength = reg_cov_from_mat(fc_var_hcpd,cov_hcpd,mask);
sc_var_hcpd_reg_strength = reg_cov_from_mat(sc_var_hcpd,cov_hcpd,mask);

[r_mat(2,1),p_mat(2,1),hcpd_fc_var_sc_var,~,p_mat(2,2)] = corr_matrix(fc_var_hcpd_reg_strength,sc_var_hcpd_reg_strength,net_label,net_order,corr_method,perm_id);

hcpd_fc_var_sc_var.mat_b = mapminmax(hcpd_fc_var_sc_var.mat_b',sc_var_min,sc_var_max)';
writetable(hcpd_fc_var_sc_var,[working_dir '/hcpd_fc_var_sc_var_reg_strength.csv'])

% HCP-YA
fc_hcp = load([mat_dir 'fc_hcp.mat']);
fc_hcp = fc_hcp.fc_mat;

load([mat_dir 'sc_communicability_hcp.mat'])
cov_hcp(:,:,1) = fc_hcp;
cov_hcp(:,:,2) = log10(sc_communicability_hcp);

sc_var_min = min(mat2vec(sc_var_hcp(idx_nolimbic, idx_nolimbic)));
sc_var_max = max(mat2vec(sc_var_hcp(idx_nolimbic, idx_nolimbic)));

fc_var_hcp_reg_strength = reg_cov_from_mat(fc_var_hcp,cov_hcp,mask);
sc_var_hcp_reg_strength = reg_cov_from_mat(sc_var_hcp,cov_hcp,mask);

[r_mat(2,2),p_mat(2,3),hcp_fc_var_sc_var,~,p_mat(2,4)] = corr_matrix(fc_var_hcp_reg_strength,sc_var_hcp_reg_strength,net_label,net_order,corr_method,perm_id);

hcp_fc_var_sc_var.mat_b = mapminmax(hcp_fc_var_sc_var.mat_b',sc_var_min,sc_var_max)';
writetable(hcp_fc_var_sc_var,[working_dir '/hcp_fc_var_sc_var_reg_strength.csv'])

%% save results
corr_hcpd(2:end,2:end) = [num2cell(r_mat(:,1)),num2cell(p_mat(:,1:2))];
corr_hcp(2:end,2:end) = [num2cell(r_mat(:,2)),num2cell(p_mat(:,3:4))];
corr_results = [corr_hcpd;corr_hcp];
writecell(corr_results,[working_dir 'fc_var_sc_var_corr_reg_dist_strength.csv'])
