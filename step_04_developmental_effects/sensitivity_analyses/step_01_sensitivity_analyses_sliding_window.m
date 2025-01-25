clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))
working_dir = [root_dir 'step_04_developmental_effects/sensitivity_analyses/'];
cd(working_dir)

atlas = 'schaefer400';
window_list = 40:10:60;
step_list = [5,10];

load(['7net_label_' atlas '.mat'])
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

%% sort subjects based on age
load([root_dir 'data/sub_info/hcpd_sublist.mat'])
N = length(hcpd_sublist_id);

hcpd_info = readtable([root_dir 'data/sub_info/hcpd_subinfo.csv']);
hcpd_sex = hcpd_info.sex;
hcpd_sex_num = ones(length(hcpd_sex),1);
idx = find(ismember(hcpd_sex,'F'));
hcpd_sex_num(idx) = 0;

% get the averaged head motion and gender proportion in each age group
hcpd_HM = hcpd_info.mean_FD;

hcpd_age_raw = hcpd_info.interview_age / 12;
hcpd_info = [hcpd_age_raw,hcpd_sex_num,hcpd_HM];

[age_sort,age_sort_idx] = sort(hcpd_age_raw,'ascend');
hcpd_info = hcpd_info(age_sort_idx,:);

%% HCP-D
fc_dir = [root_dir 'data/fc/' atlas '/'];
load([fc_dir '/FC_' atlas '_8run_hcpd.mat'],'all_session');
hcpd_fc = all_session;

for window_i = 1:length(window_list)
    window_length = window_list(window_i);

    for step_i = 1:length(step_list)
        step_length = step_list(step_i);
        progressbar(window_i/length(window_list), step_i/length(step_list))

        %% divide the dataset using sliding windows
        age_group_idx = get_sliding_windows(N,window_length,step_length);
        group_num = length(age_group_idx);

        Age = zeros(group_num,1);
        Sex = zeros(group_num,1);
        HeadMotion = zeros(group_num,1);

        for i = 1:group_num
            idx = age_group_idx{i};
            Age(i,1) = mean(hcpd_info(idx,1));
            Sex(i,1) = sum(hcpd_info(idx,2))./length(idx);
            HeadMotion(i,1) = mean(hcpd_info(idx,3));
        end

        %% calculate the fc variability and axis slope for each group
        idx_tril = find(tril(ones(6,6)));
        group_num = length(Age);

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

        clc

        Slope = slope_mat;
        tbl = table(Age,Sex,HeadMotion,Slope,'VariableNames',{'Age','Sex','HeadMotion','Slope'});
        writetable(tbl,[working_dir '/gam_age_axis_slope_window_' num2str(window_length) '_step_' num2str(step_length) '.csv'])
    end
end
