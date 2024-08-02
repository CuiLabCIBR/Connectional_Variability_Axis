clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'hcpd_sublist.mat'])

[~,~,hcpd_info] = xlsread([root_dir 'data/sub_info/hcpd_subinfo.csv']);

hcpd_cog = {'sub_id','age','gender','flanker','lswmt','dccs','total_ef','fluid_cogcomp','fd_mean'};

%% Intersection
[~,~,flanker] = xlsread([subinfo_dir 'behav/flanker01.xlsx']);
flanker_sub = flanker(2:end,1);

[~,~,lswmt] = xlsread([subinfo_dir 'behav/lswmt01.xlsx']);
lswmt_sub = lswmt(2:end,1);

[~,~,dccs] = xlsread([subinfo_dir 'behav/dccs01.xlsx']);
dccs_sub = dccs(2:end,1);

[~,~,cogcomp] = xlsread([subinfo_dir 'behav/cogcomp01.xlsx']);
cogcomp_sub = cogcomp(2:end,1);

sub_cog = intersect(flanker_sub,dccs_sub);
sub_cog = intersect(sub_cog,lswmt_sub);
sub_cog = intersect(sub_cog,cogcomp_sub);
sub_cog = intersect(sub_cog,hcpd_sublist_id);

N = length(sub_cog);

%%
idx = find(ismember(hcpd_sublist_id,sub_cog));
hcpd_cog(2:N+1,1:3) = hcpd_info(idx+1,1:3);

% 4.flanker 5.lswmt 6.dccs 7.total_ef 8.fluid_cogcomp
flanker_idx = find(ismember(flanker_sub,sub_cog));
hcpd_cog(2:N+1,4) = flanker(flanker_idx+1,4);

lswmt_idx = find(ismember(lswmt_sub,sub_cog));
hcpd_cog(2:N+1,5) = lswmt(lswmt_idx+1,4);

dccs_idx = find(ismember(dccs_sub,sub_cog));
hcpd_cog(2:N+1,6) = dccs(dccs_idx+1,4);

hcpd_cog(2:N+1,7) = num2cell(sum(cell2mat(hcpd_cog(2:N+1,4:6)),2));

cogcomp_idx = find(ismember(cogcomp_sub,sub_cog));
hcpd_cog(2:N+1,8) = cogcomp(cogcomp_idx+1,6);

%% EF
total_ef = cell2mat(hcpd_cog(2:end,7));
idx_nan = find(isnan(total_ef)); % exclude invalid values

cogcomp_score = cell2mat(hcpd_cog(2:end,8));
idx_outlier = find(cogcomp_score == 999); % exclude invalid values

idx_rm = union(idx_nan,idx_outlier);

hcpd_cog(idx_rm+1,:) = [];
sub_cog(idx_rm) = [];

%% Age and Gender
hcpd_cog(2:end,2) = num2cell(cell2mat(hcpd_cog(2:end,2))./12); % Age, month2year

hcpd_gender = hcpd_cog(2:end,3);
hcpd_gender_num = ones(length(hcpd_gender),1);
idx = find(ismember(hcpd_gender,'F'));
hcpd_gender_num(idx) = 0;
hcpd_cog(2:end,3) = num2cell(hcpd_gender_num);

%% HeadMotion
idx = find(ismember(hcpd_sublist_id,sub_cog));
fd_mean = hcpd_info(2:end,end);
hcpd_cog(2:end,end) = fd_mean(idx);
    
%%
for sub_i = 1:length(sub_cog)
    sub_name = sub_cog{sub_i};
    idx_hcpd_cog(sub_i,1) = find(ismember(hcpd_sublist_id,sub_name));
end

save([subinfo_dir 'hcpd_cog.mat'],'hcpd_cog','idx_hcpd_cog')
