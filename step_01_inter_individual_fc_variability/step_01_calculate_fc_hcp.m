function step_01_calculate_fc_hcp(sub_name)
% ---------------------------------------------------------------------------------------------------------------
% This script is used to calculate the functional connectivity.
% Schaefer400 atlas was used, generating a 400*400 matrix.
% The ROI_ID was ordered by Yeo7networks.
% ---------------------------------------------------------------------------------------------------------------

display(['Processing: ' sub_name])

addpath(genpath('/ibmgpfs/cuizaixu_lab/yanghang/code/cifti-matlab-master/'))
addpath(genpath('/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/functions/'))

root_dir = '/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/';
working_dir = [root_dir 'hcp/'];
atlas_path = [root_dir 'parcellation_files/atlas-4S456Parcels/atlas-4S456Parcels_space-fsLR_den-91k_dseg.dlabel.nii'];
data_path = '/ibmgpfs/cuizaixu_lab/yanghang/data/xcpd_0.2.2/xcp_d/params_36P/hcp/xcp_d/';

cd(working_dir)

sess_list = {'REST1_acq-LR','REST1_acq-RL','REST2_acq-LR','REST2_acq-RL'};

data_all = [];
for sess_i = 1:length(sess_list)
    sess_name = sess_list{sess_i};
    nii_name = [sub_name '_task-' sess_name  '_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
    file_path = [data_path sub_name '/func/' nii_name];
    
    data = parcellate_cifti(file_path, atlas_path);
    data = data(1:400, :); % exclude subcortical regions
    data_all = [data_all,data];
end

[ROI_num,TP_num] = size(data_all);

%% 12 run
outpath_12run = [root_dir '/results/fc_matrix_xcpd_0.2.2/hcp/schaefer400/twelve_run/'];
if ~exist(outpath_12run)
    mkdir(outpath_12run)
end

TP_12run = round(TP_num/12);

for sess_i = 1:12
    sess_now = data_all(:, (sess_i-1)*TP_12run+1 : sess_i*TP_12run);
    FC = atanh(corrcoef(sess_now'));
    FC(1:ROI_num+1:end) = 0;
    FC_12run(:, :, sess_i) = FC;
end

save([outpath_12run filesep sub_name '.mat'], 'FC_12run');
