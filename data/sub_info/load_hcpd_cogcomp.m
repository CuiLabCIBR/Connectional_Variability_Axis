clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'hcpd_sublist.mat'])
load([subinfo_dir 'hcpd_subject_info.mat'])

hcpd_cog = {'sub_id','age','gender','fluid_cogcomp','fd_mean'};

%% Intersection
[~,~,cogcomp] = xlsread([subinfo_dir 'hcpd_cogcomp.csv']);
cogcomp_sub = cogcomp(2:end,1);
sub_cogcomp = intersect(cogcomp_sub,hcpd_sublist_id);

%% EF
idx = find(ismember(cogcomp_sub,sub_cogcomp));
hcpd_cog(2:length(idx)+1,1:4) = cogcomp(idx+1,[1:3,6]);

cogcomp_score = cell2mat(hcpd_cog(2:end,4));
idx_nan = find(cogcomp_score == 999); % exclude invalid values

hcpd_cog(idx_nan+1,:) = [];
cogcomp_score(idx_nan) = [];
sub_cogcomp(idx_nan) = [];

%% Age and Gender
hcpd_cog(2:end,2) = num2cell(cell2mat(hcpd_cog(2:end,2))./12); % Age, month2year

hcpd_gender = hcpd_cog(2:end,3);
hcpd_gender_num = ones(length(hcpd_gender),1);
idx = find(ismember(hcpd_gender,'F'));
hcpd_gender_num(idx) = 0;
hcpd_cog(2:end,3) = num2cell(hcpd_gender_num);

%% HeadMotion
[~,~,HM_rest] = xlsread([subinfo_dir 'hcpd_rest_head_motion.xlsx']);
HM_rest = HM_rest(2:end,[1,3]);
[~,~,HM_task] = xlsread([subinfo_dir 'hcpd_task_head_motion.xlsx']);
HM_task = HM_task(2:end,[1,4]);

HM_all = [HM_rest;HM_task];

for sub_i = 1:length(sub_cogcomp)
    sub_name = sub_cogcomp{sub_i};
    idx = find(ismember(HM_all(:,1),sub_name));
    hcpd_cog{sub_i+1,5} = mean(cell2mat(HM_all(idx,2)));
end
    
%%
for sub_i = 1:length(sub_cogcomp)
    sub_name = sub_cogcomp{sub_i};
    idx_hcpd_cog(sub_i,1) = find(ismember(hcpd_sublist_id,sub_name));
end

save([subinfo_dir 'hcpd_cog.mat'],'hcpd_cog','idx_hcpd_cog')
