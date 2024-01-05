% ---------------------------------------------------------------------------------------------------------------
% This script is used to prepare data for fc variability calculation.
% The data is save as a 2D matrix, dim: (sess * subji) * edges.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Gradient/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_01_inter_individual_fc_variability/hcp/'];
fc_dir = [root_dir 'data/fc/schaefer400/'];

%% schaefer400
load([fc_dir '/FC_schaefer400_12run_hcp.mat'],'all_session');
hcp_fc = permute(all_session,[2 1 3]);
[n1,n2,n3] = size(hcp_fc);
% n1 = 12 sessions, n2 = 245 subjects, n3 = 79800 edges
hcp_fc = reshape(hcp_fc,[n1*n2,n3]);
save([fc_dir 'hcp_fc_schaefer400.mat'],'hcp_fc')

%% subject/session information
load([root_dir 'data/sub_info/hcp_sublist.mat'],'hcp_sublist_id')
subID = repmat(hcp_sublist_id,[1,n1])';
subID = subID(:);
save([fc_dir 'subID_hcp.mat'],'subID')

for i = 1:n1
    sess_str{i,1} = ['sess-' num2str(i)];
end
session = repmat(sess_str,[n2,1]);
save([fc_dir 'session_hcp.mat'],'session')

