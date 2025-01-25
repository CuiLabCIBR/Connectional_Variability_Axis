clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_04_developmental_effects/sensitivity_analyses/'];
cd(working_dir)

subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'yen_sublist.mat'])

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN

%% sort subjects based on age
yen_info = readtable([root_dir 'data/sub_info/yen_subinfo.csv']);
yen_sex = yen_info.sex;
yen_sex_num = ones(length(yen_sex),1);
idx = find(ismember(yen_sex,'F'));
yen_sex_num(idx) = 0;

% get the averaged head motion and gender proportion in each age group
yen_HM = yen_info.mean_FD;

yen_age_raw = yen_info.age;
yen_info = [yen_age_raw,yen_sex_num,yen_HM];

[age_sort,age_sort_idx] = sort(yen_age_raw,'ascend');
yen_info = yen_info(age_sort_idx,:);

%% divide the dataset using sliding windows
window_length = 50;
N = length(yen_sublist_id);
step_length = 5;

age_group_idx = get_sliding_windows(N,window_length,step_length);
group_num = length(age_group_idx);

Age = zeros(group_num,1);
Sex = zeros(group_num,1);
HeadMotion = zeros(group_num,1);

for i = 1:group_num
    idx = age_group_idx{i};
    Age(i,1) = mean(yen_info(idx,1));
    Sex(i,1) = sum(yen_info(idx,2))./length(idx);
    HeadMotion(i,1) = mean(yen_info(idx,3));
end

%% calculate the fc variability and axis slope for each group
group_num = length(Age);

fc_dir = [root_dir 'data/fc/schaefer400/'];
load([fc_dir 'FC_schaefer400_3run_yen.mat'],'all_session');
yen_fc = all_session;

idx_tril = find(tril(ones(6,6)));
slope_mat = zeros(group_num,1);

for group_i = 1:group_num
    group_i
    group_idx = age_sort_idx(age_group_idx{group_i});
    fc = yen_fc(group_idx,:,:);

    [~,~,edge_num] = size(yen_fc);
    icc = zeros(edge_num,1);

    for edge_i = 1:edge_num
        data_now = fc(:,:,edge_i);
        [icc(edge_i,1),~,~] = ICC(3,'single',data_now);
    end
    icc(icc < 0) = 0;

    fc_variability_mat = squareform(icc);
    fc_variability_mat_mean = get_matrix_mean(fc_variability_mat,net_label,net_order);
    fc_variability_mean(group_i,:) = fc_variability_mat_mean(idx_tril);

    % axis slope
    var_temp = fc_variability_mean(group_i,:);
    var_temp = sort(var_temp,'descend');
    mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
    slope_mat(group_i,1) = abs(table2array(mdl.Coefficients(2,1)));
end

%%
clc

Slope = slope_mat;
tbl = table(Age,Sex,HeadMotion,Slope,'VariableNames',{'Age','Sex','HeadMotion','Slope'});
writetable(tbl,[working_dir '/gam_age_axis_slope_yen.csv'])
