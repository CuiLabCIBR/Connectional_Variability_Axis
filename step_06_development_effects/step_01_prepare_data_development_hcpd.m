% ---------------------------------------------------------------------------------------------------------------
% This script was used to prepare data of Figure 5
% <Development of the connectional variability axis during youth>
% Participants from HCP-D were divided using sliding windows, length = 50, step = 5.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))
mat_dir = [root_dir 'data/connectome_matrix/schaefer400/'];

working_dir = [root_dir 'step_06_development_effects/'];
data_dir = [root_dir 'data/fc_variability/schaefer400/age_effects/'];
subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'hcpd_sublist.mat'])

%% sort subjects based on age
fc_dir = [root_dir 'data/fc/schaefer400/'];
load([fc_dir '/FC_schaefer400_8run_hcpd.mat'],'all_session');
hcpd_fc = permute(all_session,[2 1 3]);

load([root_dir 'data/sub_info/hcpd_subject_info.mat'])
hcpd_gender = hcpd_info(:,3);
hcpd_gender_num = ones(length(hcpd_gender),1);
idx = find(ismember(hcpd_gender,'F'));
hcpd_gender_num(idx) = 0;

hcpd_age_raw = cell2mat(hcpd_info(:,2));
[age_sort,age_sort_idx] = sort(hcpd_age_raw,'ascend');

%% get the averaged head motion and gender proportion in each age group
[~,~,HM_rest] = xlsread([subinfo_dir 'hcpd_rest_head_motion.xlsx']);
HM_rest = HM_rest(2:end,[1,3]);
[~,~,HM_task] = xlsread([subinfo_dir 'hcpd_task_head_motion.xlsx']);
HM_task = HM_task(2:end,[1,4]);

HM_all = [HM_rest;HM_task];

for sub_i = 1:length(hcpd_sublist_id)
    sub_name = hcpd_sublist_id{sub_i};
    idx = find(ismember(HM_all(:,1),sub_name));
    hcpd_HM(sub_i,1) = mean(cell2mat(HM_all(idx,2)));
end

hcpd_info = [hcpd_age_raw,hcpd_gender_num,hcpd_HM];
hcpd_info = hcpd_info(age_sort_idx,:);
hcpd_info(:,1) = hcpd_info(:,1) ./ 12;

%% divide the dataset using sliding windows
window_length = 50;
N = 415;
step_length = 5;

age_group_idx = get_sliding_windows(N,window_length,step_length);
group_num = length(age_group_idx);

for i = 1:group_num
    idx = age_group_idx{i};
    Age(i,1) = mean(hcpd_info(idx,1));
    Gender(i,1) = sum(hcpd_info(idx,2))./length(idx);
    HeadMotion(i,1) = mean(hcpd_info(idx,3));
end

save([working_dir 'hcpd_development_info.mat'],'Age','Gender','HeadMotion')

%% save fc for each age group
for i = 1:group_num
    i
    idx = age_sort_idx(age_group_idx{i});
    fc = hcpd_fc(:,idx,:);
    fc = reshape(fc,[8*length(idx),79800]);
    save([fc_dir 'age_effects/fc_' num2str(i) '.mat'],'fc')
end

%% subject/session information
for i = 1:50
    subID_str{i,1} = ['subID-' num2str(i)];
end
subID = repmat(subID_str,[1,8])';
subID = subID(:);
save([fc_dir 'age_effects/subID_hcpd.mat'],'subID')

for i = 1:8
    sess_str{i,1} = ['sess-' num2str(i)];
end
session = repmat(sess_str,[50,1]);
save([fc_dir 'age_effects/session_hcpd.mat'],'session')