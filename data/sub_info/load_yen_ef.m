clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'yen_sublist.mat'])
sub_demo = yen_sublist;

yen_info = readtable([root_dir 'data/sub_info/yen_subinfo.csv']);

%% Intersection
% We used fluid_cogcomp score in this study
yen_cog_info = readtable([subinfo_dir 'behav/yen_cog.xlsx']);
fluid_cogcomp = yen_cog_info.standardized_score;

sub_cog = yen_cog_info.MRIdata_ID;
sub_demo_cog = intersect(sub_cog,sub_demo);
N = length(sub_demo_cog);

%% demo
[~, idx_demo] = ismember(sub_demo_cog,sub_demo);
sub_id = yen_info.sub_id(idx_demo);
age = yen_info.age(idx_demo);
sex = yen_info.sex(idx_demo);

%% cog
[~, idx_cog] = ismember(sub_demo_cog,sub_cog);
fluid_cogcomp = yen_cog_info.standardized_score(idx_cog);

%% motion
fd_mean = yen_info.mean_FD(idx_demo);

%% save results
idx_yen_cog = idx_demo;
yen_coginfo = table(sub_id, age, sex, fluid_cogcomp, fd_mean, idx_yen_cog, ...
    'VariableNames', {'sub_id', 'age', 'sex', 'fluid_cogcomp', 'fd_mean', 'idx_yen_cog'});
writetable(yen_coginfo, [subinfo_dir 'yen_coginfo.csv'])
save([subinfo_dir 'yen_coginfo.mat'],'yen_coginfo','idx_yen_cog')
