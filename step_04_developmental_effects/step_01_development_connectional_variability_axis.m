clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))

working_dir = [root_dir 'step_04_developmental_effects/'];
cd(working_dir)

subinfo_dir = [root_dir 'data/sub_info/'];

load('7net_label_schaefer400.mat')
% We excluded the limbic network in this study.
net_order = [1 2 3 4 6 7]; %1 VIS 2 SMN 3 DAN 4 VAN 5 LIM 6 FPN 7 DMN

%% sort subjects based on age
hcpd_info = readtable([root_dir 'data/sub_info/hcpd_subinfo.csv']);
hcpd_sex = hcpd_info.sex;
hcpd_sex_num = ones(length(hcpd_sex),1);
idx = find(ismember(hcpd_sex,'F'));
hcpd_sex_num(idx) = 0;

% get the averaged head motion and gender proportion in each age group
hcpd_HM = hcpd_info.mean_FD;

hcpd_age_raw = hcpd_info.interview_age / 12;
hcpd_info = [hcpd_age_raw,hcpd_sex_num,hcpd_HM];

[age_sort,age_sort_idx] = sort(hcpd_age_raw,'ascend');
hcpd_info = hcpd_info(age_sort_idx,:);

%% divide the dataset using sliding windows
window_length = 50;
N = length(hcpd_info);
step_length = 5;

age_group_idx = get_sliding_windows(N,window_length,step_length);
group_num = length(age_group_idx);

Age = zeros(group_num,1);
Sex = zeros(group_num,1);
HeadMotion = zeros(group_num,1);

for i = 1:group_num
    idx = age_group_idx{i};
    Age(i,1) = mean(hcpd_info(idx,1));
    Sex(i,1) = sum(hcpd_info(idx,2))./length(idx);
    HeadMotion(i,1) = mean(hcpd_info(idx,3));
end

save([working_dir 'hcpd_development_info.mat'],'Age','Sex','HeadMotion','age_group_idx','age_sort_idx')

%% calculate the fc variability and axis slope for each group
% load([working_dir 'hcpd_development_info.mat'])
group_num = length(Age);

fc_dir = [root_dir 'data/fc/schaefer400/'];
load([fc_dir 'FC_schaefer400_8run_hcpd.mat'],'all_session');
hcpd_fc = all_session;

idx_tril = find(tril(ones(6,6)));
slope_mat = zeros(group_num,1);

for group_i = 1:group_num
    group_i
    group_idx = age_sort_idx(age_group_idx{group_i});
    fc = hcpd_fc(group_idx,:,:);

    [~,~,edge_num] = size(hcpd_fc);
    icc = zeros(edge_num,1);

    for edge_i = 1:edge_num
        data_now = fc(:,:,edge_i);
        [icc(edge_i,1),~,~] = ICC(3,'single',data_now);
    end
    icc(icc < 0) = 0;

    fc_variability_mat = squareform(icc);
    fc_variability_mat_mean = get_matrix_mean(fc_variability_mat,net_label,net_order);
    fc_variability_mean(group_i,:) = fc_variability_mat_mean(idx_tril);

    % axis slope
    var_temp = fc_variability_mean(group_i,:);
    var_temp = sort(var_temp,'descend');
    mdl = fitlm((1:length(var_temp))./length(var_temp),var_temp);
    slope_mat(group_i,1) = abs(table2array(mdl.Coefficients(2,1)));
end

%%
clc

Slope = slope_mat;
tbl = table(Age,Sex,HeadMotion,Slope,'VariableNames',{'Age','Sex','HeadMotion','Slope'});
writetable(tbl,[working_dir '/gam_age_axis_slope.csv'])

%% for single edge
net_name = {'VS','SM','DA','VA','FP','DM'};
net_name_mat = cell(6,6);

for i = 1:6
    for j = 1:6
        net_name_mat{i,j} = strcat(net_name{i},'_',net_name{j});
    end
end

net_name_vec = net_name_mat(idx_tril);

for i = 1:length(net_name_vec)
    Variability = fc_variability_mean(:,i);
    tbl = table(Age,Sex,HeadMotion,Variability,'VariableNames',{'Age','Sex','HeadMotion','Variability'});
    writetable(tbl,[working_dir '/gam_age_edge'  num2str(i,'%.2d') '_'  net_name_vec{i} '.csv'])
end

%% connectional axis similarity between age groups
corr_axis = corr(fc_variability_mean',fc_variability_mean',"type","Spearman");
imagesc(corr_axis)
colorbar
clim([0.9,1])

load('cmap.mat','cmap')
colormap(cmap); 

set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
axis square;
set(gcf, 'units', 'inches', 'position', [0, 0, 5, 5], 'PaperUnits', 'inches', 'PaperSize', [5, 5])

print(gcf,'-dpng','-r300',[working_dir '/corr_age_groups_connectional_axis.png'])
close all
