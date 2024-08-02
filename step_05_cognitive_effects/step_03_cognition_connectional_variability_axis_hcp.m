% ---------------------------------------------------------------------------------------------------------------
% This script was used to generate HCP-YA results of Figure 6
% <The connectional variability axis is associated with the higher-order cognitions>
% Participants from HCP-YA were divided using sliding windows, length = 50, step = 5.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_07_cognitive_effects/'];
data_dir = [root_dir 'data/fc_variability/schaefer400/cognitive_effects/hcp/'];

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

%% calculate the fc variability and inter-edge Axis for each group
load([root_dir,'step_07_cognitive_effects/hcp_cognitive_info.mat'])
group_num = length(Cognition);

%%
idx = find(tril(ones(6,6)));

axis_label = {'slope_net'};
axis_mat = zeros(group_num,1);

for group_i = 1:group_num
    load([data_dir 'lme_' num2str(group_i) '.mat'])
    icc = lme.ICC_c;
    fc_variability_mat = squareform(icc);
    fc_variability_mat_mean = get_matrix_mean(fc_variability_mat,net_label);
    fc_variability_mean(group_i,:) = fc_variability_mat_mean(idx);

    % Axis slope
    var_temp = fc_variability_mean(group_i,:);
    var_temp = sort(var_temp,'descend');
    mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
    axis_mat(group_i,1) = abs(table2array(mdl.Coefficients(2,1)));
end

%%
clc

for i = 1
    Axis = axis_mat(:,i);
    tbl = table(Age,Gender,HeadMotion,Cognition,Axis,'VariableNames',{'Age','Gender','HeadMotion','Cognition','Axis'});

    % step_04_plot_gam_cognition.R
    writetable(tbl,[working_dir 'gam_cog_hcp_' axis_label{i} '.csv'])
end