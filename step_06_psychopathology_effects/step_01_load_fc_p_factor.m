clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))
working_dir = [root_dir 'step_06_psychopathology_effects/'];

sub_func_info = readtable([working_dir 'behav/sub_func_info.csv']);
sublist_func = sub_func_info.sub_id;

%% schaefer400
load([working_dir 'fc/fc_covbat_schaefer400.mat'])

% Baseline 0-year
sub_p_factor = readtable([working_dir 'behav/sub_p_factor_0y.csv']);
sublist_p = sub_p_factor.sub_id;

idx = cellfun(@(x) find(ismember(sublist_func,x)),sublist_p);
fc_p_factor = fc_covbat(idx,:,:);
save([working_dir 'fc/fc_p_factor_0y_schaefer400.mat'],'fc_p_factor','-v7.3')

% Follow-up 2-year
sub_p_factor = readtable([working_dir 'behav/sub_p_factor_2y.csv']);
sublist_p = sub_p_factor.sub_id;

idx = cellfun(@(x) find(ismember(sublist_func,x)),sublist_p);
fc_p_factor = fc_covbat(idx,:,:);
save([working_dir 'fc/fc_p_factor_2y_schaefer400.mat'],'fc_p_factor','-v7.3')

%% glasser360
load([working_dir 'fc/fc_covbat_glasser360.mat'])

% Baseline
sub_p_factor = readtable([working_dir 'behav/sub_p_factor_0y.csv']);
sublist_p = sub_p_factor.sub_id;

idx = cellfun(@(x) find(ismember(sublist_func,x)),sublist_p);
fc_p_factor = fc_covbat(idx,:,:);
save([working_dir 'fc/fc_p_factor_0y_glasser360.mat'],'fc_p_factor','-v7.3')

% Follow-up 2-year
sub_p_factor = readtable([working_dir 'behav/sub_p_factor_2y.csv']);
sublist_p = sub_p_factor.sub_id;

idx = cellfun(@(x) find(ismember(sublist_func,x)),sublist_p);
fc_p_factor = fc_covbat(idx,:,:);
save([working_dir 'fc/fc_p_factor_2y_glasser360.mat'],'fc_p_factor','-v7.3')