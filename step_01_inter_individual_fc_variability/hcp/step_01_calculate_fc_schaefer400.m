function step_01_calculate_fc_schaefer400(sub_i)

% ---------------------------------------------------------------------------------------------------------------
% This script is used to calculate the functional connectivity.
% Schaefer400 atlas was used, generating a 400*400 matrix.
% The ROI_ID was reordered from Yeo17networks to Yeo7networks.
% ---------------------------------------------------------------------------------------------------------------

addpath(genpath('/ibmgpfs/cuizaixu_lab/yanghang/code/cifti-matlab-master/'))
addpath(genpath('/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/functions/'))

root_dir = '/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/';
working_dir = '/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/hcp/';
data_path = '/ibmgpfs/cuizaixu_lab/yanghang/data/xcpd_0.1.0/no_censor/hcp/xcp_d/';

cd(working_dir)
Schaefer400_17to7 = load('Schaefer400_17to7.mat');
Schaefer400_17to7 = Schaefer400_17to7.Schaefer400_17to7;

hcp_sublist = load([working_dir 'hcp_sublist.mat']);
hcp_sublist = hcp_sublist.hcp_sublist;

if isstr(sub_i)
    sub_i = str2num(sub_i);
end
sub_id = hcp_sublist{sub_i};
sub_name = ['sub-' sub_id];

sess_list = {'REST1_acq-LR','REST1_acq-RL','REST2_acq-LR','REST2_acq-RL'};

data_all = [];
for sess_i = 1:length(sess_list)
    sess_name = sess_list{sess_i};
    nii_name = [sub_name '_task-' sess_name '_space-fsLR_atlas-Schaefer417_den-91k_bold.ptseries.nii'];
    file_path = [data_path sub_name '/func/' nii_name];
    
    data = cifti_read(file_path);
    data = data.cdata;
    data = data(Schaefer400_17to7,:);
    data_all = [data_all,data];
end

[ROI_num,TP_num] = size(data_all);

%% 12 run
outpath_12run = [root_dir '/results/fc_matrix/hcp/schaefer400/twelve_run/'];
if ~exist(outpath_12run)
    mkdir(outpath_12run)
end

TP_12run = round(TP_num/12);

for sess_i = 1:12
    sess_now = data_all(:,(sess_i-1)*TP_12run+1:sess_i*TP_12run);
    FC = atanh(corrcoef(sess_now'));
    FC(1:ROI_num+1:end) = 0;
    FC_12run(:,:,sess_i) = FC;
end

save([outpath_12run filesep sub_name '.mat'],'FC_12run');

