clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))
working_dir = [root_dir 'step_05_cognitive_effects/'];
cd(working_dir)

%% hcp
cog_effcts = readtable('gam_cognition_edge_stats_hcp.csv');
cog_effcts_t = cog_effcts.t_gam;
cog_effcts_p = cog_effcts.p_gam;

[pthr,pcor,padj] = fdr_adjust(cog_effcts_p);
cog_effcts.matrices(find(padj <= 0.05))

t_mat = zeros(6,6);
idx_tril = find(tril(ones(6,6)));

t_mat(idx_tril) = cog_effcts_t;
t_mat = (t_mat + t_mat');
t_mat(1:7:end) = t_mat(1:7:end)/2;

t_7mat = [t_mat(1:4,:);zeros(1,6);t_mat(5:6,:)];
t_7mat = [t_7mat(:,1:4),zeros(7,1),t_7mat(:,5:6)];

load(['7net_label_schaefer400.mat'])
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

plot_net_mat(t_7mat, net_label, net_order)
caxis([-6.5,6.5])

print(gcf,'-dpng','-r300',[working_dir 'edge_fc_var_cog_effects_hcp.png'])
close all

cog_effects_axis = get_connectional_axis(t_mat);
writecell(cog_effects_axis,[working_dir 'cog_effects_axis_hcp.xlsx'])

%% hcpd
cog_effcts = readtable('gam_cognition_edge_stats_hcpd.csv');
cog_effcts_t = cog_effcts.t_gam;
cog_effcts_p = cog_effcts.p_gam;

[pthr,pcor,padj] = fdr_adjust(cog_effcts_p);
cog_effcts.matrices(find(padj <= 0.05))

t_mat = zeros(6,6);
idx_tril = find(tril(ones(6,6)));

t_mat(idx_tril) = cog_effcts_t;
t_mat = (t_mat + t_mat');
t_mat(1:7:end) = t_mat(1:7:end)/2;

t_7mat = [t_mat(1:4,:);zeros(1,6);t_mat(5:6,:)];
t_7mat = [t_7mat(:,1:4),zeros(7,1),t_7mat(:,5:6)];

load(['7net_label_schaefer400.mat'])
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

plot_net_mat(t_7mat, net_label, net_order)
caxis([-6,6])

print(gcf,'-dpng','-r300',[working_dir 'edge_fc_var_cog_effects_hcpd.png'])
close all

cog_effects_axis = get_connectional_axis(t_mat);
writecell(cog_effects_axis,[working_dir 'cog_effects_axis_hcpd.xlsx'])
