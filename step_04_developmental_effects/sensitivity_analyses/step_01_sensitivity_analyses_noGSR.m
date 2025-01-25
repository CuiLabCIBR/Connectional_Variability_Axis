clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_04_developmental_effects/sensitivity_analyses/'];
cd(working_dir)

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN

%% 
fc_dir = [root_dir 'data/fc/schaefer400/'];
load([fc_dir 'FC_schaefer400_8run_hcpd_noGSR.mat'],'all_session');
hcpd_fc = all_session;

info_dir = [root_dir 'step_04_developmental_effects/'];
load([info_dir 'hcpd_development_info.mat'])
group_num = length(Age);

%% calculate the fc variability and axis slope for each group
idx_tril = find(tril(ones(6,6)));

slope_mat = zeros(group_num,1);

for group_i = 1:group_num
    group_i
    group_idx = age_sort_idx(age_group_idx{group_i});
    fc = hcpd_fc(group_idx,:,:);

    [~,~,edge_num] = size(hcpd_fc);
    icc = zeros(edge_num,1);

    for edge_i = 1:edge_num
        data_now = fc(:,:,edge_i);
        [icc(edge_i,1),~,~] = ICC(3,'single',data_now);
    end
    icc(icc < 0) = 0;

    fc_variability_mat = squareform(icc);
    fc_variability_mat_mean = get_matrix_mean(fc_variability_mat,net_label,net_order);
    fc_variability_mean = fc_variability_mat_mean(idx_tril);

    % axis slope
    var_temp = sort(fc_variability_mean,'descend');
    mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
    slope_mat(group_i,1) = abs(table2array(mdl.Coefficients(2,1)));
end

%%
clc

Slope = slope_mat;
tbl = table(Age,Sex,HeadMotion,Slope,'VariableNames',{'Age','Sex','HeadMotion','Slope'});
writetable(tbl,[working_dir '/gam_age_axis_slope_noGSR.csv'])
