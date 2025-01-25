clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_04_developmental_effects/sensitivity_analyses/'];
cd(working_dir)

dist_dir = [root_dir 'data/parcellation_files/'];
load([dist_dir 'dist_schaefer400.mat'],'dist_schaefer400')

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.W

%% for debug
fc_dir = [root_dir 'data/fc/schaefer400/'];
load([fc_dir 'FC_schaefer400_8run_hcpd.mat'],'all_session');
hcpd_fc = all_session;

info_dir = [root_dir 'step_04_developmental_effects/'];
load([info_dir 'hcpd_development_info.mat'])
group_num = length(Age);

%% calculate the fc variability and axis slope for each group
idx_tril = find(tril(ones(6,6)));

slope_label = {'axis_slope_reg_dist','axis_slope_reg_fc'};
slope_mat_reg_cov = zeros(group_num,2);

for group_i = 1:group_num
    group_i
    group_idx = age_sort_idx(age_group_idx{group_i});
    fc = hcpd_fc(group_idx,:,:);
    fc_mat = squareform(squeeze(mean(mean(fc))));

    [~,~,edge_num] = size(hcpd_fc);
    icc = zeros(edge_num,1);
    
    for edge_i = 1:edge_num
        data_now = fc(:,:,edge_i);
        [icc(edge_i,1),~,~] = ICC(3,'single',data_now);
    end
    icc(icc < 0) = 0;
    fc_variability_mat = squareform(icc);

    %% regress distance
    fc_variability_mat_reg_dist = reg_cov_from_mat(fc_variability_mat,dist_schaefer400);
    fc_variability_mat_mean_reg_dist = get_matrix_mean(fc_variability_mat_reg_dist,net_label);
    fc_variability_mean_reg_dist = fc_variability_mat_mean_reg_dist(idx_tril);

    % axis slope
    var_temp = sort(fc_variability_mean_reg_dist,'descend');
    mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
    slope_mat_reg_cov(group_i,1) = abs(table2array(mdl.Coefficients(2,1)));

    %% regress FC
    fc_variability_mat_reg_fc = reg_cov_from_mat(fc_variability_mat,fc_mat);
    fc_variability_mat_mean_reg_fc = get_matrix_mean(fc_variability_mat_reg_fc,net_label);
    fc_variability_mean_reg_fc = fc_variability_mat_mean_reg_fc(idx_tril);

    % axis slope
    var_temp = sort(fc_variability_mean_reg_fc,'descend');
    mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
    slope_mat_reg_cov(group_i,2) = abs(table2array(mdl.Coefficients(2,1)));
end

%%
clc

for i = 1:2
    Slope = slope_mat_reg_cov(:,i);
    tbl = table(Age,Sex,HeadMotion,Slope,'VariableNames',{'Age','Sex','HeadMotion','Slope'});
    writetable(tbl,[working_dir 'gam_age_' slope_label{i} '.csv'])
end
