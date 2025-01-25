function [data_mat_mean,net_edges] = get_matrix_mean(data_mat,net_label,net_order,data_mask)
% This function is used to get the network-avearged values of a nodal-level matrix.
% input:
%       data_mat:  N*N data matrix, N is the node number.
%       net_label: N*1 vector assigns each node to a network.
%       net_order: M*1 vector indicates the network order, M is the network number.
%       data_mask: N*N matrix, analyze the data within the mask.
% output£º
%       data_mat_mean: M*M network-avearged matrix.
%       net_edges: M*M cell matrix which stores edges within and between networks.

N = length(data_mat);
t = [];
start = 1;
lines = 1;

if exist('net_order','var')
    order = net_order;
else
    order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN
end

for i = 1:length(order)
    add = find(net_label==order(i));
    t = [t;add];
    start = start + length(add);
    lines(i+1) = start;
end

data_reorder = data_mat(t,t);

if ~exist('data_mask','var')
    data_mask = ones(N,N);
end

data_mask = data_mask(t,t);
%%
idx_begin = lines(1:end-1);
idx_end = lines(2:end)-1;

for i = 1:length(lines)-1
   for j = 1:length(lines)-1
       data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
       mask_temp = data_mask(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
       if i == j
           data_temp = mat2vec(data_temp);
           mask_temp = mat2vec(mask_temp);
       end
       data_temp = data_temp(mask_temp > 0);
       data_mat_mean(i,j) = mean(data_temp);
       net_edges{i,j} = data_temp;
   end
end