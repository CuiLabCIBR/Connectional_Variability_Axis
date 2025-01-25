clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_03_structural_connectome_variability/sensitivity_analyses/'];
var_dir = [root_dir 'data/fc_variability/schaefer400/'];
mat_dir = [root_dir 'data/connectome_matrix/schaefer400/'];

load('7net_label_schaefer400.mat')
load([root_dir 'data/parcellation_files/perm_id_schaefer400.mat'])

% Initialization Parameters
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.
idx_nolimbic = find(net_label ~= 5); % We excluded the limbic network in this study.

corr_method = 'spearman';
thr_range = 0.05:0.05:0.25;

% Load SC
load([mat_dir 'sc_hcpd.mat'])
[~,~,n_hcpd] = size(sc_all_hcpd);
sc_all_hcpd = sc_all_hcpd(idx_nolimbic, idx_nolimbic, :);

load([mat_dir 'sc_hcp.mat'])
[~,~,n_hcp] = size(sc_all_hcp);
sc_all_hcp = sc_all_hcp(idx_nolimbic, idx_nolimbic, :);

for thr_i = 1:length(thr_range)

    prop_thr = thr_range(thr_i)

    %% HCP-D
    % proportional threshold for each individual
    for i = 1:n_hcpd
        sc_all_hcpd_bi(:,:,i) = threshold_proportional(sc_all_hcpd(:,:,i), prop_thr);
    end

    % calculate communicability for each individual
    for i = 1:n_hcpd
        sc_temp = double(sc_all_hcpd_bi(:,:,i)>0);
        G_hcpd(:,:,i) = get_communicability(sc_temp,1);
        G_vec_hcpd(:,i) = mat2vec(G_hcpd(:,:,i));
    end

    % mean communicability (cmy)
    G_mean_vec_hcpd = mean(G_vec_hcpd,2);
    G_mean_mat_hcpd = squareform(G_mean_vec_hcpd);
    sc_communicability_hcpd = zeros(400, 400);
    sc_communicability_hcpd(idx_nolimbic, idx_nolimbic) = G_mean_mat_hcpd;

    % mad communicability
    G_mad_vec_hcpd = mad(G_vec_hcpd')';
    G_mad_mat_hcpd = squareform(G_mad_vec_hcpd);
    sc_variability_hcpd = zeros(400, 400);
    sc_variability_hcpd(idx_nolimbic, idx_nolimbic) = G_mad_mat_hcpd;

    % association with FC variability
    fc_variability_hcpd = load([var_dir 'fc_variability_hcpd.mat']);
    fc_variability_hcpd = fc_variability_hcpd.fc_variability;

    [r_hcpd(thr_i,1),p_hcpd(thr_i,1),hcpd_fc_var_sc_var,~,p_hcpd_spin(thr_i,1)] = corr_matrix(fc_variability_hcpd,log10(sc_variability_hcpd),net_label,net_order,corr_method,perm_id);
    % writetable(hcpd_fc_var_sc_var,[working_dir '/hcpd_fc_var_sc_var_thr' num2str(100*prop_thr) '.csv'])

    %% HCP-YA
    for i = 1:n_hcp
        sc_all_hcp_bi(:,:,i) = threshold_proportional(sc_all_hcp(:,:,i), prop_thr);
    end

    % calculate communicability for each individual
    for i = 1:n_hcp
        sc_temp = double(sc_all_hcp_bi(:,:,i)>0);
        G_hcp(:,:,i) = get_communicability(sc_temp,1);
        G_vec_hcp(:,i) = mat2vec(G_hcp(:,:,i));
    end

    % mean communicability (cmy)
    G_mean_vec_hcp = mean(G_vec_hcp,2);
    G_mean_mat_hcp = squareform(G_mean_vec_hcp);
    sc_communicability_hcp = zeros(400, 400);
    sc_communicability_hcp(idx_nolimbic, idx_nolimbic) = G_mean_mat_hcp;
    
    % mad communicability
    G_mad_vec_hcp = mad(G_vec_hcp')';
    G_mad_mat_hcp = squareform(G_mad_vec_hcp);
    sc_variability_hcp = zeros(400, 400);
    sc_variability_hcp(idx_nolimbic, idx_nolimbic) = G_mad_mat_hcp;
  
    % association with FC variability
    fc_variability_hcp = load([var_dir 'fc_variability_hcp.mat']);
    fc_variability_hcp = fc_variability_hcp.fc_variability;

    [r_hcp(thr_i,1),p_hcp(thr_i,1),hcp_fc_var_sc_var,~,p_hcp_spin(thr_i,1)] = corr_matrix(fc_variability_hcp,log10(sc_variability_hcp),net_label,net_order,corr_method,perm_id);
    % writetable(hcp_fc_var_sc_var,[working_dir '/hcp_fc_var_sc_var_thr' num2str(100*prop_thr) '.csv'])

end

%% save results
corr_hcpd = cell(6,4);
corr_hcpd(1,:) = {'hcp-d','r-value','p-value','p-spin'};
corr_hcpd(2:end,1) = num2cell(thr_range);
corr_hcpd(2:end,2:end) = num2cell([r_hcpd, p_hcpd, p_hcpd_spin]);

corr_hcp = corr_hcpd;
corr_hcp(1,1) = {'hcp-ya'};
corr_hcp(2:end,2:end) = num2cell([r_hcp, p_hcp, p_hcp_spin]);

corr_results = [corr_hcpd;corr_hcp];
writecell(corr_results,[working_dir 'fc_var_sc_var_corr_binary.csv'])
