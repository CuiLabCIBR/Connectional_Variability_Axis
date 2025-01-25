function plot_net_vec(net_vec)

var_mat = zeros(6,6);
idx_tril = find(tril(ones(6,6)));

var_mat(idx_tril) = net_vec;
var_mat = (var_mat + var_mat');
var_mat(1:7:end) = var_mat(1:7:end)/2;

var_7mat = [var_mat(1:4,:);zeros(1,6);var_mat(5:6,:)];
var_7mat = [var_7mat(:,1:4),zeros(7,1),var_7mat(:,5:6)];

load(['7net_label_schaefer400.mat'])
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN.

plot_net_mat(var_7mat, net_label, net_order)