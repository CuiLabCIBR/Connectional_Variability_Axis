% ---------------------------------------------------------------------------------------------------------------
% This script was used to prepare yen-YA data of Figure 4
% <The connectional variability axis pattern is associated with the individual differences in higher-order cognitive functions>
% Participants from yen-YA were divided using sliding windows, length = 50, step = 5.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_05_cognitive_effects/sensitivity_analyses/'];
subinfo_dir = [root_dir 'data/sub_info/'];

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

%% separate the subjects using sliding windows
yen_cog = readtable([subinfo_dir 'yen_coginfo.csv']);
yen_cog.sex = categorical(yen_cog.sex);
yen_cog.sex = double(yen_cog.sex == 'M');

cog = yen_cog.fluid_cogcomp;
[cog_sort,cog_sort_idx] = sort(cog,'ascend');
yen_cog_sort = yen_cog(cog_sort_idx,:);

window_length = 50;
N = length(cog);
step_length = 5;

cog_group_idx = get_sliding_windows(N,window_length,step_length);
group_num = length(cog_group_idx);

Age = zeros(group_num,1);
Sex = zeros(group_num,1);
Cognition = zeros(group_num,1);
HeadMotion = zeros(group_num,1);

for i = 1:group_num
    idx = cog_group_idx{i};
    Age(i,1) = mean(yen_cog_sort.age(idx));
    Sex(i,1) = sum(yen_cog_sort.sex(idx))./length(idx);
    Cognition(i,1) = mean(yen_cog_sort.fluid_cogcomp(idx));
    HeadMotion(i,1) = mean(yen_cog_sort.fd_mean(idx));
end

% save([working_dir 'yen_cognitive_info.mat'],'Age','Sex','HeadMotion','Cognition','cog_group_idx','cog_sort_idx')

%% calculate the fc variability and inter-edge Axis for each group
% load([working_dir, 'yen_cognitive_info.mat'])
idx_yen_cog = yen_cog.idx_yen_cog;

fc_dir = [root_dir 'data/fc/schaefer400/'];
load([fc_dir '/FC_schaefer400_3run_yen.mat'],'all_session');

yen_fc = all_session;
yen_fc = yen_fc(idx_yen_cog,:,:); %get subjects with cog

group_num = length(Cognition);
idx_tril = find(tril(ones(6,6)));

slope_mat = zeros(group_num,1);

%%
for group_i = 1:group_num
    group_i
    group_idx = cog_sort_idx(cog_group_idx{group_i});
    fc = yen_fc(group_idx,:,:);

    [~,~,edge_num] = size(yen_fc);
    icc = zeros(edge_num,1);

    for edge_i = 1:edge_num
        data_now = fc(:,:,edge_i);
        [icc(edge_i,1),~,~] = ICC(3,'single',data_now);
    end
    icc(icc < 0) = 0;

    fc_variability_mat = squareform(icc);
    fc_variability_mat_mean = get_matrix_mean(fc_variability_mat,net_label);
    fc_variability_mean(group_i,:) = fc_variability_mat_mean(idx_tril);

    % axis slope
    var_temp = fc_variability_mean(group_i,:);
    var_temp = sort(var_temp,'descend');
    mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
    slope_mat(group_i,1) = abs(table2array(mdl.Coefficients(2,1)));
end

%%
clc

Slope = slope_mat;
tbl = table(Age,Sex,HeadMotion,Cognition,Slope,'VariableNames',{'Age','Sex','HeadMotion','Cognition','Slope'});
writetable(tbl,[working_dir 'gam_cog_yen_axis_slope.csv'])
