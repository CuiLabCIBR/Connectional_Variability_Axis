% ---------------------------------------------------------------------------------------------------------------
% This script is used to prepare data for fc variability calculation.
% The data is save as a 2D matrix, dim: (sess * subji) * edges.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Gradient/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_01_inter_individual_fc_variability/hcpd/'];
fc_dir = [root_dir 'data/fc/schaefer400/'];

%% schaefer400
load([fc_dir '/FC_schaefer400_8run_hcpd.mat'],'all_session');
hcpd_fc = permute(all_session,[2 1 3]);
[n1,n2,n3] = size(hcpd_fc);
% n1 = 8 sessions, n2 = 415 subjects, n3 = 79800 edges
hcpd_fc = reshape(hcpd_fc,[n1*n2,n3]);
save([fc_dir 'hcpd_fc_schaefer400.mat'],'hcpd_fc')

%% subject/session information
load([root_dir 'data/sub_info/hcpd_sublist.mat'],'hcpd_sublist_id')
subID = repmat(hcpd_sublist_id,[1,n1])';
subID = subID(:);
save([fc_dir 'subID_hcpd.mat'],'subID')

for i = 1:n1
    sess_str{i,1} = ['sess-' num2str(i)];
end
session = repmat(sess_str,[n2,1]);
save([fc_dir 'session_hcpd.mat'],'session')

