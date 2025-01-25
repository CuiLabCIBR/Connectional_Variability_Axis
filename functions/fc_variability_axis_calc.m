function [fc_variability_mat,fc_variability_mean,axis_slope] = fc_variability_axis_calc(fc)

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

%%
[~,~,edge_num] = size(fc);

for edge_i = 1:edge_num
    data_now = fc(:,:,edge_i);
    [icc(edge_i,1),sigma2_b(edge_i,1),sigma2_w(edge_i,1)] = ICC(3,'single',data_now);
end
icc(icc < 0) = 0;

fc_variability_mat = squareform(icc);
fc_variability_mat_mean = get_matrix_mean(fc_variability_mat,net_label);
fc_variability_mean = fc_variability_mat_mean(tril(ones(6,6))==1);

% gradient slope
var_temp = fc_variability_mean;
var_temp = sort(var_temp,'descend');
mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
axis_slope = abs(table2array(mdl.Coefficients(2,1)));
