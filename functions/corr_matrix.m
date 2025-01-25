function [r,p,mat_a_b,r_spin,p_spin] = corr_matrix(mat_a,mat_b,net_label,net_order,corr_method,perm_id)
%-------------------------------------------------------------------------------------%
% This function is used to calculate the correlation between two matrices.
% input:
%       mat_a: N*N data matrix a, N is the node number.
%       mat_b: N*N data matrix b.
%       perm_id: N*K, K is the repetition number.
%       net_label: N*1 vector assigns each node to a network.
%       net_order: M*1 vector indicates the network order, M is the network number.
%       corr_method: the correlation method, default: 'spearman'.
%       perm_id: N*K, K is the repetition number.
%                shuffle the label of the matrix, and do permutation test.
% output：
%       r: the correlation coefficient between mat_a and mat_b.
%       p: the significance level of the correlation.
%       mat_a_b: a table variabile includes the unique elements from mat_a and mat_b.
%       r_spin: r values of the permutation test.
%       p_spin: p value of the permutation test.
%-------------------------------------------------------------------------------------%

mat_a_raw = mat_a;
net_label_raw = net_label;

if ~exist('corr_method')
    corr_method = 'spearman';
end

%% reorder the data matrix
remove_net = setdiff(unique(net_label),net_order);

t = 1:length(net_label);
if ~isempty(remove_net)
    for i = 1:length(remove_net)
        idx_remove = find(net_label_raw == remove_net(i));
        t = setdiff(t,idx_remove);
    end

    mat_a = mat_a(t,t);
    mat_b = mat_b(t,t);
end
mat_a_vec = mat2vec(mat_a);
mat_b_vec = mat2vec(mat_b);

[r,p] = corr(mat_a_vec,mat_b_vec,'type',corr_method);

%% save the results
mat_a_b.mat_a = mat_a_vec;
mat_a_b.mat_b = mat_b_vec;
mat_a_b = struct2table(mat_a_b);

%%
if exist('perm_id')
    [~,rand_num] = size(perm_id);
    for i = 1:rand_num
        idx = perm_id(:,i);
        net_label_rand = net_label(idx);
        mat_a_rand = mat_a_raw(idx,idx);
        if ~isempty(remove_net)
            t = find(net_label_rand ~= remove_net);
            mat_a_rand = mat_a_rand(t,t);
        end
        mat_a_rand_vec = mat2vec(mat_a_rand);
        r_spin(i,1) = corr(mat_a_rand_vec,mat_b_vec,'type',corr_method);
    end

    if r > 0
        p_spin = (1+length(find(r_spin >= r)))/(rand_num+1);
    else
        p_spin = (1+length(find(r_spin <= r)))/(rand_num+1);
    end

end
