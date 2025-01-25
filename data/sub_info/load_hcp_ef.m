clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'hcp_sublist.mat'])
sub_demo = hcp_sublist_id;

hcp_info = readtable([root_dir 'data/sub_info/hcp_subinfo.csv']);

%% Intersection
% We used fluid_cogcomp score in this study
hcp_cog_info = readtable([subinfo_dir 'behav/hcp_cog_info.csv']);
fluid_cogcomp = hcp_cog_info.CogFluidComp_Unadj;
idx_nan = find(isnan(fluid_cogcomp));
hcp_cog_info(idx_nan,:) = []; %remove subjects without scores

sub_cog = num2cell(hcp_cog_info.Subject);
for i = 1:length(sub_cog)
    sub_cog{i} = num2str(sub_cog{i});
end

sub_demo_cog = intersect(sub_cog,sub_demo);
N = length(sub_demo_cog);

%% demo
[~, idx_demo] = ismember(sub_demo_cog,sub_demo);
sub_id = hcp_info.Subject(idx_demo);
age = hcp_info.Age_in_Yrs(idx_demo);
sex = hcp_info.Sex(idx_demo);

%% cog
[~, idx_cog] = ismember(sub_demo_cog,sub_cog);
fluid_cogcomp = hcp_cog_info.CogFluidComp_Unadj(idx_cog);

%% motion
fd_mean = hcp_info.mean_FD(idx_demo);

%% save results
idx_hcp_cog = idx_demo;
hcp_coginfo = table(sub_id, age, sex, fluid_cogcomp, fd_mean, idx_hcp_cog, ...
    'VariableNames', {'sub_id', 'age', 'sex', 'fluid_cogcomp', 'fd_mean', 'idx_hcp_cog'});
writetable(hcp_coginfo, [subinfo_dir 'hcp_coginfo.csv'])
save([subinfo_dir 'hcp_coginfo.mat'],'hcp_coginfo','idx_hcp_cog')