% ---------------------------------------------------------------------------------------------------------------
% This script is used to load the FC from all participants.
% The data is save as a 3D matrix, dim: subji * sess * edges.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

addpath(genpath('/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/functions/'))

root_dir = '/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/';
working_dir = [root_dir,'hcp/'];
out_path = [root_dir,'results/fc_matrix_xcpd_0.2.2/hcp/schaefer400/'];

cd(working_dir)

hcp_sublist = load([working_dir 'hcp_sublist.mat']);
hcp_sublist = hcp_sublist.hcp_sublist;
sub_num = length(hcp_sublist);

if isempty(gcp('nocreate'))
    % If not, start a parallel pool with 24 workers
    parpool(24);
end

%% FC estimation based on 12 runs
FC_12run = [root_dir '/results/fc_matrix_xcpd_0.2.2/hcp/schaefer400/twelve_run/'];

clear all_session
all_session = zeros(sub_num,12,79800);

parfor sub_i = 1:sub_num  % for each subject
    sub_i
    sub_name = hcp_sublist{sub_i};
    FC_path = [FC_12run sub_name '.mat'];
    sub_data = load(FC_path);
    sub_data = sub_data.FC_12run;
       
    for run_i = 1:12     
        all_session(sub_i,run_i,:) = mat2vec(sub_data(:,:,run_i));
    end
       
end

save([out_path,'/FC_schaefer400_12run_hcp.mat'],'all_session','-v7.3');
