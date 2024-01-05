clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'hcp_sublist.mat'])
load([subinfo_dir 'hcp_subject_info.mat'])

hcp_cog = {'sub_id','age','gender','fluid_cogcomp','fd_mean'};

%% Intersection
[~,~,cogcomp] = xlsread([subinfo_dir 'hcp_cogcomp.csv']);
cogcomp_sub = cogcomp(2:end,1);
for i = 1:length(cogcomp_sub)
    cogcomp_sub{i} = num2str(cogcomp_sub{i});
end

sub_cog = intersect(cogcomp_sub,hcp_sublist_id);
N = length(sub_cog);

%% EF
idx = find(ismember(hcp_sublist_id,sub_cog));
hcp_cog(2:N+1,1:3) = hcp_info(idx,1:3);

idx = find(ismember(cogcomp_sub,sub_cog));
hcp_cog(2:N+1,4) = cogcomp(idx+1,4);

cog_score = cell2mat(hcp_cog(2:end,4));
idx_nan = find(isnan(cog_score));

hcp_cog(idx_nan+1,:) = [];
cog_score(idx_nan) = [];
sub_cog(idx_nan) = [];

%% Gender
hcp_gender = hcp_cog(2:end,3);
hcp_gender_num = ones(length(hcp_gender),1);
idx = find(ismember(hcp_gender,'F'));
hcp_gender_num(idx) = 0;
hcp_cog(2:end,3) = num2cell(hcp_gender_num);

%% HeadMotion
[~,~,HM_rest] = xlsread([subinfo_dir 'hcp_rest_head_motion.xlsx']);
HM_rest = HM_rest(2:end,[1,4]);
for i = 1:length(HM_rest)
    HM_rest{i,1} = num2str(HM_rest{i,1});
end

HM_all = HM_rest;

for sub_i = 1:length(sub_cog)
    sub_name = sub_cog{sub_i};
    idx = find(ismember(HM_all(:,1),sub_name));
    hcp_cog{sub_i+1,5} = mean(cell2mat(HM_all(idx,2)));
end

%%
for sub_i = 1:length(sub_cog)
    sub_name = sub_cog{sub_i};
    idx_hcp_cog(sub_i,1) = find(ismember(hcp_sublist_id,sub_name));
end

save([subinfo_dir 'hcp_cog.mat'],'hcp_cog','idx_hcp_cog')
