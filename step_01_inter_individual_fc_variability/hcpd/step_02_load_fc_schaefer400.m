% ---------------------------------------------------------------------------------------------------------------
% This script is used to load the FC from all participants.
% The data is save as a 3D matrix, dim: subji * sess * edges.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

addpath(genpath('/ibmgpfs/cuizaixu_lab/yanghang/code/cifti-matlab-master/'))
addpath(genpath('/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/functions/'))

root_dir = '/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/';
working_dir = [root_dir,'hcpd/'];
data_path = '/ibmgpfs/cuizaixu_lab/yanghang/data/xcpd_0.1.0/no_censor/hcpd/xcp_d/';
out_path = [root_dir,'results/fc_matrix/hcpd/schaefer400/'];

cd(working_dir)

hcpd_sublist = load([working_dir 'hcpd_sublist.mat']);
hcpd_sublist = hcpd_sublist.hcpd_sublist;
sub_num = length(hcpd_sublist);

if isempty(gcp('nocreate'))
    % If not, start a parallel pool with 24 workers
    parpool(24);
end

%% FC estimation based on 4 runs
FC_4run = [root_dir '/results/fc_matrix/hcpd/schaefer400/four_run/'];

clear all_session
all_session = zeros(sub_num,4,79800);

parfor sub_i = 1:sub_num  % for each subject
    sub_i
    sub_name = hcpd_sublist{sub_i};
    FC_path = [FC_4run sub_name '.mat'];
    sub_data = load(FC_path);
    sub_data = sub_data.FC_4run;
       
    for run_i = 1:4     
        all_session(sub_i,run_i,:) = mat2vec(sub_data(:,:,run_i));
    end
       
end

save([out_path,'/FC_schaefer400_4run_hcpd.mat'],'all_session','-v7.3');

%% FC estimation based on 8 runs
FC_8run = [root_dir '/results/fc_matrix/hcpd/schaefer400/eight_run/'];

clear all_session
all_session = zeros(sub_num,8,79800);

parfor sub_i = 1:sub_num  % for each subject
    sub_i
    sub_name = hcpd_sublist{sub_i};
    FC_path = [FC_8run sub_name '.mat'];
    sub_data = load(FC_path);
    sub_data = sub_data.FC_8run;
       
    for run_i = 1:8     
        all_session(sub_i,run_i,:) = mat2vec(sub_data(:,:,run_i));
    end
       
end

save([out_path,'/FC_schaefer400_8run_hcpd.mat'],'all_session','-v7.3');
