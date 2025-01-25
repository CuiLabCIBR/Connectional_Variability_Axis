% ---------------------------------------------------------------------------------------------------------------
% This script is used to prepare data for fc variability calculation.
% The data is save as a 2D matrix, dim: (sess * subji) * edges.
% ---------------------------------------------------------------------------------------------------------------

clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

atlas_list = {'schaefer400', 'glasser360'};

%% hcp
for i = 1:length(atlas_list)
    atlas = atlas_list{i};
    fc_dir = [root_dir 'data/fc/' atlas '/'];

    load([fc_dir '/FC_' atlas '_12run_hcp.mat'],'all_session');
    fc_data = permute(all_session,[2 1 3]);
    [n1,n2,n3] = size(fc_data);
    % n1 = session num, n2 = subject num, n3 = edge num
    fc_data = reshape(fc_data,[n1*n2,n3]);
    save([fc_dir 'fc_hcp_' atlas '.mat'],'fc_data','-v7.3')

    % save group averaged FC
    out_dir = [root_dir 'data/connectome_matrix/' atlas '/'];
    fc_mat = squareform(squeeze(mean(mean(all_session))));
    save([out_dir 'fc_hcp.mat'],'fc_mat')

    % for noGSR data
    if strcmp(atlas, 'schaefer400')
        load([fc_dir '/FC_' atlas '_12run_hcp_noGSR.mat'],'all_session');
        fc_data = permute(all_session,[2 1 3]);
        [n1,n2,n3] = size(fc_data);
        % n1 = session num, n2 = subject num, n3 = edge num
        fc_data = reshape(fc_data,[n1*n2,n3]);
        save([fc_dir 'fc_hcp_' atlas '_noGSR.mat'],'fc_data','-v7.3')
    end

    % subject/session information
    load([root_dir 'data/sub_info/hcp_sublist.mat'],'hcp_sublist_id')
    subID = repmat(hcp_sublist_id,[1,n1])';
    subID = subID(:);
    save([fc_dir 'subID_hcp.mat'],'subID')

    for i = 1:n1
        sess_str{i,1} = ['sess-' num2str(i)];
    end
    session = repmat(sess_str,[n2,1]);
    save([fc_dir 'session_hcp.mat'],'session')
end

%% hcpd
clear sess_str

for i = 1:length(atlas_list)
    atlas = atlas_list{i};
    fc_dir = [root_dir 'data/fc/' atlas '/'];

    load([fc_dir '/FC_' atlas '_8run_hcpd.mat'],'all_session');
    fc_data = permute(all_session,[2 1 3]);
    [n1,n2,n3] = size(fc_data);
    % n1 = session num, n2 = subject num, n3 = edge num
    fc_data = reshape(fc_data,[n1*n2,n3]);
    save([fc_dir 'fc_hcpd_' atlas '.mat'],'fc_data','-v7.3')

    % save group averaged FC
    out_dir = [root_dir 'data/connectome_matrix/' atlas '/'];
    fc_mat = squareform(squeeze(mean(mean(all_session))));
    save([out_dir 'fc_hcpd.mat'],'fc_mat')

    % for noGSR data
    if strcmp(atlas, 'schaefer400')
        load([fc_dir '/FC_' atlas '_8run_hcpd_noGSR.mat'],'all_session');
        fc_data = permute(all_session,[2 1 3]);
        [n1,n2,n3] = size(fc_data);
        % n1 = session num, n2 = subject num, n3 = edge num
        fc_data = reshape(fc_data,[n1*n2,n3]);
        save([fc_dir 'fc_hcpd_' atlas '_noGSR.mat'],'fc_data','-v7.3')
    end

    % subject/session information
    load([root_dir 'data/sub_info/hcpd_sublist.mat'],'hcpd_sublist_id')
    subID = repmat(hcpd_sublist_id,[1,n1])';
    subID = subID(:);
    save([fc_dir 'subID_hcpd.mat'],'subID')

    for i = 1:n1
        sess_str{i,1} = ['sess-' num2str(i)];
    end
    session = repmat(sess_str,[n2,1]);
    save([fc_dir 'session_hcpd.mat'],'session')
end
