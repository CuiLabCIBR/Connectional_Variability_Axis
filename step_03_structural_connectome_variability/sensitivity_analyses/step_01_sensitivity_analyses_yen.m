clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_03_structural_connectome_variability/sensitivity_analyses/'];
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
load([mat_dir 'sc_yen.mat'])
[~,~,n_yen] = size(sc_all_yen);
sc_all_yen = sc_all_yen(idx_nolimbic, idx_nolimbic, :);

W_thr = threshold_consistency(sc_all_yen, consistency_thr);
sc_mask_yen = double(W_thr > 0);

% calculate communicability for each individual
for i = 1:n_yen
    sc_temp = sc_all_yen(:,:,i) .* sc_mask_yen;
    G_yen(:,:,i) = get_communicability(sc_temp,1);
    G_vec_yen(:,i) = mat2vec(G_yen(:,:,i));
end

% mean communicability (cmy)
G_mean_vec_yen = mean(G_vec_yen,2);
G_mean_mat_yen = squareform(G_mean_vec_yen);
sc_communicability_yen = zeros(400, 400);
sc_communicability_yen(idx_nolimbic, idx_nolimbic) = G_mean_mat_yen;

% mad communicability
G_mad_vec_yen = mad(G_vec_yen')';
G_mad_mat_yen = squareform(G_mad_vec_yen);
sc_variability_yen = zeros(400, 400);
sc_variability_yen(idx_nolimbic, idx_nolimbic) = G_mad_mat_yen;

save([mat_dir 'sc_variability_yen.mat'],'sc_variability_yen')

% association with fc variability
% load([mat_dir 'sc_variability_yen.mat'],'sc_variability_yen')

fc_variability_yen = load([var_dir 'fc_variability_yen.mat']);
fc_variability_yen = fc_variability_yen.fc_variability;

[r_yen,p_yen,yen_fc_var_sc_var,~,p_yen_spin] = corr_matrix(fc_variability_yen,log10(sc_variability_yen),net_label,net_order,corr_method,perm_id);
writetable(yen_fc_var_sc_var,[working_dir '/yen_fc_var_sc_var_schaefer400.csv']) 

%% save results
corr_yen = cell(2,4);
corr_yen(1,:) = {'yen','r-value','p-value','p-spin'};
corr_yen(2,1) = {'schaefer400'};
corr_yen(2,2:end) = num2cell([r_yen, p_yen, p_yen_spin]);

corr_results = corr_yen;
writecell(corr_results,[working_dir 'fc_var_sc_var_corr_yen.csv'])