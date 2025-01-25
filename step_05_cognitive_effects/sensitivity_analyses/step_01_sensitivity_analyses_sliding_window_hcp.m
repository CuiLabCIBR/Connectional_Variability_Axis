clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))
working_dir = [root_dir 'step_05_cognitive_effects/sensitivity_analyses/'];
cd(working_dir)

subinfo_dir = [root_dir 'data/sub_info/'];
info_dir = [root_dir 'step_05_cognitive_effects/'];

window_list = 40:10:60;
step_list = [5,10];

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

%% sort subjects based on cognition
hcp_cog = readtable([subinfo_dir 'hcp_coginfo.csv']);
hcp_cog.sex = categorical(hcp_cog.sex);
hcp_cog.sex = double(hcp_cog.sex == 'M');

cog = hcp_cog.fluid_cogcomp;
N = length(cog);

[cog_sort,cog_sort_idx] = sort(cog,'ascend');
hcp_cog_sort = hcp_cog(cog_sort_idx,:);

idx_hcp_cog = hcp_cog.idx_hcp_cog;

fc_dir = [root_dir 'data/fc/schaefer400/'];
load([fc_dir '/FC_schaefer400_12run_hcp.mat'],'all_session');
hcp_fc = all_session;
hcp_fc = hcp_fc(idx_hcp_cog,:,:); % get subjects with cog

%%
for window_i = 1:length(window_list)
    window_length = window_list(window_i);

    for step_i = 1:length(step_list)
        step_length = step_list(step_i);
        progressbar(window_i/length(window_list), step_i/length(step_list))

        %% divide the dataset using sliding windows
        cog_group_idx = get_sliding_windows(N,window_length,step_length);
        group_num = length(cog_group_idx);

        Age = zeros(group_num,1);
        Sex = zeros(group_num,1);
        Cognition = zeros(group_num,1);
        HeadMotion = zeros(group_num,1);

        for i = 1:group_num
            idx = cog_group_idx{i};
            Age(i,1) = mean(hcp_cog_sort.age(idx));
            Sex(i,1) = sum(hcp_cog_sort.sex(idx))./length(idx);
            Cognition(i,1) = mean(hcp_cog_sort.fluid_cogcomp(idx));
            HeadMotion(i,1) = mean(hcp_cog_sort.fd_mean(idx));
        end

        %%
        idx_tril = find(tril(ones(6,6)));
        group_num = length(Age);

        slope_mat = zeros(group_num,1);

        for group_i = 1:group_num
            group_i
            group_idx = cog_sort_idx(cog_group_idx{group_i});
            fc = hcp_fc(group_idx,:,:);

            [~,~,edge_num] = size(hcp_fc);
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

        clc

        Slope = slope_mat;
        tbl = table(Age,Sex,HeadMotion,Cognition,Slope,'VariableNames',{'Age','Sex','HeadMotion','Cognition','Slope'});

        writetable(tbl,[working_dir '/gam_cog_hcp_axis_slope_window_' num2str(window_length) '_step_' num2str(step_length) '.csv'])
    end
end
