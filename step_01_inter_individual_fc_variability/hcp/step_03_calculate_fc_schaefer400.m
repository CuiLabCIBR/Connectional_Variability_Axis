% ---------------------------------------------------------------------------------------------------------------
% This script is used to calculate the functional connectivity.
% Schaefer400 atlas was used, generating a 400*400 matrix.
% Timepoints with FD > 0.3 were excluded from the correlation analysis.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

addpath(genpath('/GPFS/cuizaixu_lab_permanent/yanghang/projects/p01_connectional_variability_gradient/functions/'));

sub_dir = '/GPFS/cuizaixu_lab_permanent/yanghang/data/xcpd_0.1.0/hcp_d/xcp_d/';
working_dir = '/GPFS/cuizaixu_lab_permanent/yanghang/projects/p01_connectional_variability_gradient/hcp';
bold_dir = [working_dir '/bold_signal/'];

load('/GPFS/cuizaixu_lab_permanent/yanghang/projects/p01_connectional_variability_gradient/functions/hcp_sublist.mat')
sub_num = length(hcp_sublist);

%% FC estimation based on 12 runs
clear all_session_mat
clear all_session

fc_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Gradient/data/fc/schaefer400/';

for sub_i = 1:sub_num  % for each subject
    sub_i
    sub_name = hcp_sublist{sub_i};
    parcel_name = [bold_dir '/schaefer400_12run/' sub_name filesep 'data_parcel.mat'];    
    sub_data = load(parcel_name);sub_data = sub_data.data_parcel_12run;
    
    load([sub_dir sub_name '/func/' sub_name '_censor.mat'],'censor_12run')
       
    for run_i = 1:12
        censor_flag = censor_12run(:,run_i);
        schaefer400 = squeeze(sub_data(2:end,censor_flag,run_i));
        
        z_corr_data = atanh(corr(schaefer400'));
        all_session_mat(sub_i,run_i,:,:) = z_corr_data;
        all_session(sub_i,run_i,:) = mat2vec(z_corr_data);
    end
       
end

save([fc_dir '/FC_schaefer400_12run_hcp.mat'],'all_session','-v7.3');

%%
% load([fc_dir '/FC_schaefer400_12run_hcp.mat'],'all_session');
fc_mat = squareform(squeeze(mean(mean(all_session))));
save('F:/Cui_Lab/Projects/Connectional_Variability_Gradient/data/connectome_matrix/schaefer400/fc_hcp.mat','fc_mat');