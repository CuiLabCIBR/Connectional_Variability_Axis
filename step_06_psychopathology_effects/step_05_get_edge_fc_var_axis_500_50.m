clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))
working_dir = [root_dir 'step_06_psychopathology_effects/'];
cd(working_dir)

net_label = {'VS','SM','DA','VA','FP','DM'};
edge_label = cell(6,6);
for i = 1:6
    for j = 1:6
        edge_label{i,j} = [net_label{i} '_' net_label{j}];
    end
end

idx_tril = find(tril(ones(6,6)));
edge_label = edge_label(idx_tril);

load(['7net_label_schaefer400.mat'])
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

%%
atlas_list = {'schaefer400'};
year_list = {'no_reg_0y','no_reg_2y'};

for atlas_i = 1:length(atlas_list)
    atlas = atlas_list{atlas_i};

    for year_i = 1:length(year_list)
        data_year = year_list{year_i};

        [edge_var,~,raw] = xlsread([working_dir '/results/axis/gam_edge_general_' atlas '_' data_year '.csv']);
        edge_var_t = edge_var(:,3);
        edge_var_p = edge_var(:,4);

        [pthr,pcor,padj] = fdr_adjust(edge_var_p);
        edge_label(find(edge_var_p <= 0.05))
        edge_label(find(padj <= 0.05))

        t_mat = zeros(6,6);
        t_mat(idx_tril) = edge_var_t;
        t_mat = (t_mat + t_mat');
        t_mat(1:7:end) = t_mat(1:7:end)/2;

        t_7mat = [t_mat(1:4,:);zeros(1,6);t_mat(5:6,:)];
        t_7mat = [t_7mat(:,1:4),zeros(7,1),t_7mat(:,5:6)];

        plot_net_mat(t_7mat, net_label, net_order)
        thr = round(max(abs(edge_var_t)));
        caxis([-thr,thr])
        % colorbar

        print(gcf,'-dpng','-r300',[working_dir '/results/axis/edge_fc_var_net_' atlas '_' data_year '.png'])
        close all

        slope_effects_axis = get_connectional_axis(t_mat,'ascend');
        writecell(slope_effects_axis,[working_dir '/results/axis/edge_fc_var_axis_' atlas '_' data_year '.xlsx'])

    end

end
