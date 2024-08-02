function step_01_calculate_fc_schaefer400(sub_i)

% ---------------------------------------------------------------------------------------------------------------
% This script is used to calculate the functional connectivity.
% Schaefer400 atlas was used, generating a 400*400 matrix.
% The ROI_ID was reordered from Yeo17networks to Yeo7networks.
% ---------------------------------------------------------------------------------------------------------------

addpath(genpath('/ibmgpfs/cuizaixu_lab/yanghang/code/cifti-matlab-master/'))
addpath(genpath('/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/functions/'))

root_dir = '/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/';
working_dir = '/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/hcpd/';
data_path = '/ibmgpfs/cuizaixu_lab/yanghang/data/xcpd_0.1.0/no_censor/hcpd/xcp_d/';

cd(working_dir)
Schaefer400_17to7 = load('Schaefer400_17to7.mat');
Schaefer400_17to7 = Schaefer400_17to7.Schaefer17to7;

hcpd_sublist = load([working_dir 'hcpd_sublist.mat']);
hcpd_sublist = hcpd_sublist.hcpd_sublist;

if isstr(sub_i)
    sub_i = str2num(sub_i);
end

sub_name = hcpd_sublist{sub_i};
sess_list = {'REST1_acq-AP','REST1_acq-PA','REST2_acq-AP','REST2_acq-PA',...
    'CARIT_acq-AP','CARIT_acq-PA','EMOTION_acq-PA','GUESSING_acq-AP','GUESSING_acq-PA'};

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

%% 8 run
outpath_8run = [root_dir '/results/fc_matrix/hcpd/schaefer400/eight_run/'];

for sess_i = 1:8
    sess_now = data_all(:,(sess_i-1)*400+1:sess_i*400);
    FC = atanh(corrcoef(sess_now'));
    FC(1:401:end) = 0;
    FC_8run(:,:,sess_i) = FC;
end

save([outpath_8run filesep sub_name '.mat'],'FC_8run');

