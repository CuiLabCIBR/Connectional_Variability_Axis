clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_05_cognitive_effects/sensitivity_analyses/'];
subinfo_dir = [root_dir 'data/sub_info/'];
info_dir = [root_dir 'step_05_cognitive_effects/'];

load('7net_label_glasser360.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; % 1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

dataset_list = {'hcpd','hcp'};

for dataset_i = 1:length(dataset_list)
    
    dataset = dataset_list{dataset_i};

    data_cog = readtable([subinfo_dir, dataset, '_coginfo.csv']);
    idx_data_cog = table2array(data_cog(:,end));

    load([info_dir, dataset, '_cognitive_info.mat'])
    
    fc_dir = [root_dir 'data/fc/glasser360/'];

    if strcmp(dataset, 'hcpd')
        load([fc_dir '/FC_glasser360_8run_hcpd.mat'],'all_session');
    elseif strcmp(dataset, 'hcp')
        load([fc_dir '/FC_glasser360_12run_hcp.mat'],'all_session');
    end

    data_fc = all_session;
    data_fc = data_fc(idx_data_cog,:,:); % get subjects with cog

    group_num = length(Cognition);
    idx_tril = find(tril(ones(6,6)));

    slope_mat = zeros(group_num,1);

    %% calculate the fc variability and axis slope for each group
    for group_i = 1:group_num
        group_i
        group_idx = cog_sort_idx(cog_group_idx{group_i});
        fc = data_fc(group_idx,:,:);

        [~,~,edge_num] = size(data_fc);
        icc = zeros(edge_num,1);
        
        for edge_i = 1:edge_num
            data_now = fc(:,:,edge_i);
            [icc(edge_i,1),~,~] = ICC(3,'single',data_now);
        end
        icc(icc < 0) = 0;

        fc_variability_mat = squareform(icc);
        fc_variability_mat_mean = get_matrix_mean(fc_variability_mat,net_label);
        fc_variability_mean = fc_variability_mat_mean(idx_tril);

        % axis slope
        var_temp = sort(fc_variability_mean,'descend');
        mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
        slope_mat(group_i,1) = abs(table2array(mdl.Coefficients(2,1)));
    end

    %%
    clc

    Slope = slope_mat;
    tbl = table(Age,Sex,HeadMotion,Cognition,Slope,'VariableNames',{'Age','Sex','HeadMotion','Cognition','Slope'});
    writetable(tbl,[working_dir 'gam_cog_', dataset, '_axis_slope_glasser360.csv'])

end