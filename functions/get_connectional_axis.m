function connectional_axis = get_connectional_axis(data_mat,direction)

if ~exist('direction')
    direction = 'descend';
end

net_label = {'VS','SM','DA','VA','FP','DM'};

idx_mat = reshape([1:36],[6,6]);
tril_idx = find(tril(ones(6,6)));
linear_idx = idx_mat(tril_idx);
[idx_i,idx_j] = ind2sub([6,6],linear_idx);
sub_idx = [idx_i,idx_j];

%%
data = data_mat(tril_idx);
[data_sort,data_idx] = sort(data,direction);
data_sub_idx = sub_idx(data_idx,:);

for i = 1:21
        data_NetA{i,1} = net_label{data_sub_idx(i,1)};
        data_NetB{i,1} = net_label{data_sub_idx(i,2)};
        data_Net{2*i-1,1} = strcat(data_NetA{i,1}, '&', data_NetB{i,1});
        data_Net{2*i-1,2} = data_NetA{i,1};
        data_Net{2*i-1,3} = data_sort(i)./2;
        data_Net{2*i,1} = strcat(data_NetA{i,1}, '&', data_NetB{i,1});
        data_Net{2*i,2} = data_NetB{i,1};
        data_Net{2*i,3} = data_sort(i)./2;    
end

connectional_axis = cell(43,3);
connectional_axis(1,:) = {'Net_All','Net_Single','Value'};
connectional_axis(2:end,:) = data_Net;