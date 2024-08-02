clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'hcp_sublist.mat'])
[~,~,hcp_info] = xlsread([root_dir 'data/sub_info/hcp_subinfo.csv']);

hcp_cog = {'sub_id','age','gender','flanker','lswmt','dccs','total_ef','fluid_cogcomp','fd_mean'};

%% Intersection
% 3.flanker 4.lswmt 5.dccs 6.total_ef 7.fluid_cogcomp
[~,~,cog] = xlsread([subinfo_dir 'behav/hcp_cog_info.csv']);
cog_score = cell2mat(cog(2:end,3:7));
idx_nan = find(isnan(sum(cog_score,2)));
cog_score(idx_nan,:)
cog(idx_nan+1,:) = []; %remove subjects without scores

cog_sub = cog(2:end,1);
for i = 1:length(cog_sub)
    cog_sub{i} = num2str(cog_sub{i});
end

sub_cog = intersect(cog_sub,hcp_sublist_id);
N = length(sub_cog);

%% cog
idx = find(ismember(hcp_sublist_id,sub_cog));
hcp_cog(2:N+1,1:3) = hcp_info(idx+1,1:3);

idx = find(ismember(cog_sub,sub_cog));
% 3.flanker 4.lswmt 5.dccs 6.total_ef 7.fluid_cogcomp
hcp_cog(2:N+1,4:8) = cog(idx+1,3:7);

%% Gender
hcp_gender = hcp_cog(2:end,3);
hcp_gender_num = ones(length(hcp_gender),1);
idx = find(ismember(hcp_gender,'F'));
hcp_gender_num(idx) = 0;
hcp_cog(2:end,3) = num2cell(hcp_gender_num);

%% HeadMotion
idx = find(ismember(hcp_sublist_id,sub_cog));
fd_mean = hcp_info(2:end,end);
hcp_cog(2:end,end) = fd_mean(idx);

%%
for sub_i = 1:length(sub_cog)
    sub_name = sub_cog{sub_i};
    idx_hcp_cog(sub_i,1) = find(ismember(hcp_sublist_id,sub_name));
end

save([subinfo_dir 'hcp_cog.mat'],'hcp_cog','idx_hcp_cog')
