clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'hcpd_sublist.mat'])
sub_demo = hcpd_sublist_id;

hcpd_info = readtable([root_dir 'data/sub_info/hcpd_subinfo.csv']);

%% Intersection
% We used fluid_cogcomp score in this study
hcpd_cog_info = readtable([subinfo_dir 'behav/hcpd_cogcomp01.xlsx']);
fluid_cogcomp = hcpd_cog_info.nih_fluidcogcomp_unadjusted;
idx_nan = find(isnan(fluid_cogcomp));
idx_999 = find(fluid_cogcomp == 999);

hcpd_cog_info(idx_nan,:) = []; % remove subjects without scores
hcpd_cog_info(idx_999,:) = []; % exclude invalid values

sub_cog = hcpd_cog_info.src_subject_id;
sub_demo_cog = intersect(sub_cog,sub_demo);
N = length(sub_demo_cog);

%% demo
[~, idx_demo] = ismember(sub_demo_cog,sub_demo);
sub_id = hcpd_info.src_subject_id(idx_demo);
age = hcpd_info.interview_age(idx_demo);
sex = hcpd_info.sex(idx_demo);

%% cog
[~, idx_cog] = ismember(sub_demo_cog,sub_cog);
fluid_cogcomp = hcpd_cog_info.nih_fluidcogcomp_unadjusted(idx_cog);

%% motion
fd_mean = hcpd_info.mean_FD(idx_demo);

%% save results
idx_hcpd_cog = idx_demo;
hcpd_coginfo = table(sub_id, age, sex, fluid_cogcomp, fd_mean, idx_hcpd_cog, ...
    'VariableNames', {'sub_id', 'age', 'sex', 'fluid_cogcomp', 'fd_mean', 'idx_hcpd_cog'});
writetable(hcpd_coginfo, [subinfo_dir 'hcpd_coginfo.csv'])
save([subinfo_dir 'hcpd_coginfo.mat'],'hcpd_coginfo','idx_hcpd_cog')
