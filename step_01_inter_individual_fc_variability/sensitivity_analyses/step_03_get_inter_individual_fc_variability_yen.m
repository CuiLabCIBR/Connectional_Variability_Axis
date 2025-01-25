% ---------------------------------------------------------------------------------------------------------------
% This script is used to get the FC variability matrix.
% The intra-class correlation (ICC) was used as
% the normalized inter-subject FC variability in this study.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

atlas_list = {'schaefer400'};
dataset_list = {'yen'};

%%
for dataset_i = 1:length(dataset_list)
    dataset = dataset_list{dataset_i};

    for atlas_i = 1:length(atlas_list)
        atlas = atlas_list{atlas_i};
        var_dir = [root_dir 'data/fc_variability/' atlas '/'];

        dataset_atlas = [dataset '_' atlas];

        load([var_dir, 'lme_' dataset_atlas '.mat'])
        icc = lme_results.ICC_c;
        icc = squareform(icc);
        fc_variability = icc;

        save([var_dir 'fc_variability_' dataset '.mat'],'fc_variability')
    end
end

