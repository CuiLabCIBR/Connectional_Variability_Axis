function fc_var_p_factor_sliding_window(sub_fc,data_year,window_length,step_length,atlas,term_label,reg_flag)

% clear
% clc

if ~exist('reg_flag')
    reg_flag = 0;
end

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

%%
load(['7net_label_' atlas '.mat'])
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

net_name = {'VS','SM','DA','VA','FP','DM'};
net_name_mat = cell(6,6);

for i = 1:6
    for j = 1:6
        net_name_mat{i,j} = strcat(net_name{j},'_',net_name{i});
    end
end

idx_tril = find(tril(ones(6,6)));
net_name_vec = net_name_mat(idx_tril);

%%
sub_info = readtable([root_dir 'step_06_psychopathology_effects/behav/sub_p_factor_' data_year '.csv']);

outpath = [root_dir 'step_06_psychopathology_effects/' term_label '/' atlas '/no_reg_' data_year '/'];
if ~exist(outpath)
    mkdir(outpath)
end

if reg_flag == 1
    dist_dir = [root_dir 'data/parcellation_files/'];
    load([dist_dir 'dist_schaefer400.mat'],'dist_schaefer400')
    outpath_reg_fc = [root_dir 'step_06_psychopathology_effects/' term_label '/' atlas '/reg_fc_' data_year '/'];
    if ~exist(outpath_reg_fc)
        mkdir(outpath_reg_fc)
    end

    outpath_reg_dist = [root_dir 'step_06_psychopathology_effects/' term_label '/' atlas '/reg_dist_' data_year '/'];
    if ~exist(outpath_reg_dist)
        mkdir(outpath_reg_dist)
    end
end

%% sort subjects based on P_factor score
% 5.general 6.external 7.internal
p_factor = sub_info.(term_label);
[pfactor_sort,pfactor_sort_idx] = sort(p_factor,'ascend');
sub_info_sort = sub_info(pfactor_sort_idx,:);

sub_fc = sub_fc(pfactor_sort_idx,:,:);
[N,~,~] = size(sub_info_sort);
pfactor_group_idx = get_sliding_windows(N,window_length,step_length);
group_num = length(pfactor_group_idx);

for i = 1:group_num
    idx = pfactor_group_idx{i};
    Age(i,1) = mean(sub_info_sort.age(idx));
    Gender(i,1) = sum(sub_info_sort.sex(idx))./length(idx);
    Pfactor(i,1) = mean(sub_info_sort.(term_label)(idx));
    HeadMotion(i,1) = mean(sub_info_sort.mean_fd(idx));
end

%%
idx_tril = find(tril(ones(6,6)));

slope_label = {'slope_net'};
slope_mat = zeros(group_num,1);

variability_label = net_name_vec;
variability_mat = zeros(group_num,length(idx_tril));

if reg_flag == 1
    slope_mat_reg_fc = zeros(group_num,1);
    slope_mat_reg_dist = zeros(group_num,1);
    variability_mat_reg_fc = zeros(group_num,length(idx_tril));
    variability_mat_reg_dist = zeros(group_num,length(idx_tril));
end

%%

for group_i = 1:group_num
    group_i
    group_idx = pfactor_group_idx{group_i};
    fc = sub_fc(group_idx,:,:);
    fc_mat = squareform(squeeze(mean(mean(fc))));

    [~,~,edge_num] = size(sub_fc);

    for edge_i = 1:edge_num
        data_now = fc(:,:,edge_i);
        [icc(edge_i,1),sigma2_b(edge_i,1),sigma2_w(edge_i,1)] = ICC(3,'single',data_now);
    end
    icc(icc < 0) = 0;

    %%
    fc_variability_mat = squareform(icc);
    fc_variability_mat_mean = get_matrix_mean(fc_variability_mat,net_label);
    fc_variability_mean = fc_variability_mat_mean(idx_tril);

    variability_mat(group_i,:) = fc_variability_mean;

    % Axis slope
    var_temp = fc_variability_mean;
    var_temp = sort(var_temp,'descend');
    mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
    slope_mat(group_i,1) = abs(table2array(mdl.Coefficients(2,1)));

    if reg_flag == 1
        %% regress FC
        fc_variability_mat = squareform(icc);
        fc_variability_mat = reg_cov_from_mat(fc_variability_mat,fc_mat);
        fc_variability_mat_mean = get_matrix_mean(fc_variability_mat,net_label);
        fc_variability_mean = fc_variability_mat_mean(idx_tril);

        variability_mat_reg_fc(group_i,:) = fc_variability_mean;

        % Axis slope
        var_temp = fc_variability_mean;
        var_temp = sort(var_temp,'descend');
        mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
        slope_mat_reg_fc(group_i,1) = abs(table2array(mdl.Coefficients(2,1)));

        %% regress distance
        fc_variability_mat = squareform(icc);
        fc_variability_mat = reg_cov_from_mat(fc_variability_mat,dist_schaefer400);
        fc_variability_mat_mean = get_matrix_mean(fc_variability_mat,net_label);
        fc_variability_mean = fc_variability_mat_mean(idx_tril);

        variability_mat_reg_dist(group_i,:) = fc_variability_mean;

        % Axis slope
        var_temp = fc_variability_mean;
        var_temp = sort(var_temp,'descend');
        mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
        slope_mat_reg_dist(group_i,1) = abs(table2array(mdl.Coefficients(2,1)));
    end
end

%%
clc

Slope = slope_mat;

tbl = table(Age,Gender,HeadMotion,Pfactor,Slope,'VariableNames',{'Age','Gender','HeadMotion','Pfactor','Slope'});
writetable(tbl,[outpath 'gam_pfactor_axis_slope_' num2str(window_length) '_' num2str(step_length) '.csv'])

if reg_flag == 1
    Slope = slope_mat_reg_fc;
    tbl = table(Age,Gender,HeadMotion,Pfactor,Slope,'VariableNames',{'Age','Gender','HeadMotion','Pfactor','Slope'});
    writetable(tbl,[outpath_reg_fc 'gam_pfactor_axis_slope_' num2str(window_length) '_' num2str(step_length) '.csv'])

    Slope = slope_mat_reg_dist;
    tbl = table(Age,Gender,HeadMotion,Pfactor,Slope,'VariableNames',{'Age','Gender','HeadMotion','Pfactor','Slope'});
    writetable(tbl,[outpath_reg_dist 'gam_pfactor_axis_slope_' num2str(window_length) '_' num2str(step_length) '.csv'])
end

for i = 1:length(variability_label)
    Variability = variability_mat(:,i);
    tbl = table(Age,Gender,HeadMotion,Pfactor,Variability,'VariableNames',{'Age','Gender','HeadMotion','Pfactor','Variability'});
    writetable(tbl,[outpath 'gam_pfactor_edge' num2str(i,'%.2d') '_' variability_label{i} '_' num2str(window_length) '_' num2str(step_length) '.csv'])

    if reg_flag == 1
        Variability = variability_mat_reg_fc(:,i);
        tbl = table(Age,Gender,HeadMotion,Pfactor,Variability,'VariableNames',{'Age','Gender','HeadMotion','Pfactor','Variability'});
        writetable(tbl,[outpath_reg_fc 'gam_pfactor_edge'  num2str(i,'%.2d') '_' variability_label{i} '_' num2str(window_length) '_' num2str(step_length) '.csv'])

        Variability = variability_mat_reg_dist(:,i);
        tbl = table(Age,Gender,HeadMotion,Pfactor,Variability,'VariableNames',{'Age','Gender','HeadMotion','Pfactor','Variability'});
        writetable(tbl,[outpath_reg_dist 'gam_pfactor_edge'  num2str(i,'%.2d') '_' variability_label{i} '_' num2str(window_length) '_' num2str(step_length) '.csv'])
    end
end

