% ---------------------------------------------------------------------------------------------------------------
% This script was used to prepare HCP-YA data of Figure 6
% <The connectional variability gradient is associated with the higher-order cognitions>
% Participants from HCP-YA were divided using sliding windows, length = 50, step = 5.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Gradient/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_07_cognitive_effects/'];
data_dir = [root_dir 'data/fc_variability/schaefer400/cognitive_effects/hcp/'];
subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'hcp_sublist.mat'])
load([subinfo_dir 'hcp_cog.mat'])

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

%% separate the subjects using sliding windows
fc_dir = [root_dir 'data/fc/schaefer400/'];
load([fc_dir '/FC_schaefer400_12run_hcp.mat'],'all_session');
hcp_fc = all_session;
hcp_fc = hcp_fc(idx_hcp_cog,:,:);
hcp_fc = permute(hcp_fc,[2 1 3]);
[n1,n2,n3] = size(hcp_fc);

% 4.fluid_cogcomp
term_i = 4;

cog = cell2mat(hcp_cog(2:end,term_i));
[cog_sort,cog_sort_idx] = sort(cog,'ascend');
hcp_cog_sort = cell2mat(hcp_cog(2:end,2:end));
hcp_cog_sort = hcp_cog_sort(cog_sort_idx,:);

window_length = 50;
N = length(idx_hcp_cog);
step_length = 5;

cog_group_idx = get_sliding_windows(N,window_length,step_length);
group_num = length(cog_group_idx);

for i = 1:group_num
    idx = cog_group_idx{i};
    Age(i,1) = mean(hcp_cog_sort(idx,1));
    Gender(i,1) = sum(hcp_cog_sort(idx,2))./length(idx);
    Cognition(i,1) = mean(hcp_cog_sort(idx,term_i-1));
    HeadMotion(i,1) = mean(hcp_cog_sort(idx,end));
end

save([working_dir 'hcp_cognitive_info.mat'],'Age','Gender','HeadMotion','Cognition')

%% save fc for each group
mkdir([fc_dir 'cognitive_effects/hcp/'])
for i = 1:group_num
    i
    idx = cog_sort_idx(cog_group_idx{i});
    fc = hcp_fc(:,idx,:);
    fc = reshape(fc,[n1*length(idx),n3]);
    save([fc_dir 'cognitive_effects/hcp/fc_' num2str(i) '.mat'],'fc')
end

%% subject/session information
for i = 1:window_length
    subID_str{i,1} = ['subID-' num2str(i)];
end
subID = repmat(subID_str,[1,n1])';
subID = subID(:);
save([fc_dir 'cognitive_effects/hcp/subID_hcp.mat'],'subID')

for i = 1:n1
    sess_str{i,1} = ['sess-' num2str(i)];
end
session = repmat(sess_str,[window_length,1]);
save([fc_dir 'cognitive_effects/hcp/session_hcp.mat'],'session')
