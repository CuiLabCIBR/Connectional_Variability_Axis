% ---------------------------------------------------------------------------------------------------------------
% This script is used to prepare data for fc variability calculation.
% The data is save as a 2D matrix, dim: (sess * subji) * edges.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

atlas_list = {'schaefer400'};

%% yen
for i = 1:length(atlas_list)
    atlas = atlas_list{i};
    fc_dir = [root_dir 'data/fc/' atlas '/'];

    load([fc_dir '/FC_' atlas '_3run_yen.mat'],'all_session');
    fc_data = permute(all_session,[2 1 3]);
    [n1,n2,n3] = size(fc_data);
    % n1 = session num, n2 = subject num, n3 = edge num
    fc_data = reshape(fc_data,[n1*n2,n3]);
    save([fc_dir 'fc_yen_' atlas '.mat'],'fc_data','-v7.3')

    % subject/session information
    load([root_dir 'data/sub_info/yen_sublist.mat'],'yen_sublist_id')
    subID = repmat(yen_sublist_id,[1,n1])';
    subID = subID(:);
    save([fc_dir 'subID_yen.mat'],'subID')

    for i = 1:n1
        sess_str{i,1} = ['sess-' num2str(i)];
    end
    session = repmat(sess_str,[n2,1]);
    save([fc_dir 'session_yen.mat'],'session')
end
