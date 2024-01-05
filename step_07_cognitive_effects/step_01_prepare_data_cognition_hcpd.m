% ---------------------------------------------------------------------------------------------------------------
% This script was used to prepare HCP-D data of Figure 6
% <The connectional variability gradient is associated with the higher-order cognitions>
% Participants from HCP-D were divided using sliding windows, length = 50, step = 5.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Gradient/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_07_cognitive_effects/'];
data_dir = [root_dir 'data/fc_variability/schaefer400/cognitive_effects/hcpd/'];
subinfo_dir = [root_dir 'data/sub_info/'];
load([subinfo_dir 'hcpd_sublist.mat'])
load([subinfo_dir 'hcpd_cog.mat'])

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

%% separate the subjects using sliding windows
fc_dir = [root_dir 'data/fc/schaefer400/'];
load([fc_dir '/FC_schaefer400_8run_hcpd.mat'],'all_session');
hcpd_fc = all_session;
hcpd_fc = hcpd_fc(idx_hcpd_cog,:,:);
hcpd_fc = permute(hcpd_fc,[2 1 3]);
[n1,n2,n3] = size(hcpd_fc);

% 4.fluid_cogcomp
term_i = 4;

cog = cell2mat(hcpd_cog(2:end,term_i));
[cog_sort,cog_sort_idx] = sort(cog,'ascend');
hcpd_cog_sort = cell2mat(hcpd_cog(2:end,2:end));
hcpd_cog_sort = hcpd_cog_sort(cog_sort_idx,:);

window_length = 50;
N = length(idx_hcpd_cog);
step_length = 5;

cog_group_idx = get_sliding_windows(N,window_length,step_length);
group_num = length(cog_group_idx);

for i = 1:group_num
    idx = cog_group_idx{i};
    Age(i,1) = mean(hcpd_cog_sort(idx,1));
    Gender(i,1) = sum(hcpd_cog_sort(idx,2))./length(idx);
    Cognition(i,1) = mean(hcpd_cog_sort(idx,term_i-1));
    HeadMotion(i,1) = mean(hcpd_cog_sort(idx,end));
end

save([working_dir 'hcpd_cognitive_info.mat'],'Age','Gender','HeadMotion','Cognition')

%% save fc for each group
mkdir([fc_dir 'cognitive_effects/hcpd/'])
for i = 1:group_num
    i
    idx = cog_sort_idx(cog_group_idx{i});
    fc = hcpd_fc(:,idx,:);
    fc = reshape(fc,[n1*length(idx),n3]);
    save([fc_dir 'cognitive_effects/hcpd/fc_' num2str(i) '.mat'],'fc')
end

%% subject/session information
for i = 1:window_length
    subID_str{i,1} = ['subID-' num2str(i)];
end
subID = repmat(subID_str,[1,n1])';
subID = subID(:);
save([fc_dir 'cognitive_effects/hcpd/subID_hcpd.mat'],'subID')

for i = 1:n1
    sess_str{i,1} = ['sess-' num2str(i)];
end
session = repmat(sess_str,[window_length,1]);
save([fc_dir 'cognitive_effects/hcpd/session_hcpd.mat'],'session')
