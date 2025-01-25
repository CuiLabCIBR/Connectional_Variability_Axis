clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))
working_dir = [root_dir 'step_04_developmental_effects/'];
cd(working_dir)

%% 
age_effcts = readtable('gam_development_edge_stats.csv');
age_effcts_R2 = age_effcts.partial_R2;
age_effcts_p = age_effcts.p_anova;

[pthr,pcor,padj] = fdr_adjust(age_effcts_p);
age_effcts.matrices(find(padj <= 0.05))

R2_mat = zeros(6,6);
idx_tril = find(tril(ones(6,6)));

R2_mat(idx_tril) = age_effcts_R2;
R2_mat = (R2_mat + R2_mat');
R2_mat(1:7:end) = R2_mat(1:7:end)/2;

R2_7mat = [R2_mat(1:4,:);zeros(1,6);R2_mat(5:6,:)];
R2_7mat = [R2_7mat(:,1:4),zeros(7,1),R2_7mat(:,5:6)];

load(['7net_label_schaefer400.mat'])
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

plot_net_mat(R2_7mat, net_label, net_order)
clim([-0.3,0.3])

print(gcf,'-dpng','-r300',[working_dir 'edge_fc_var_age_effects.png'])
close all

age_effects_axis = get_connectional_axis(R2_mat,'ascend');
writecell(age_effects_axis,[working_dir 'age_effects_axis.xlsx'])
