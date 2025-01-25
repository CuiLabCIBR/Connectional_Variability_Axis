clear
clc

working_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/data/parcellation_files/';
cd(working_dir)

glasser = cifti_read('glasser_space-fsLR_den-32k_desc-atlas.dlabel.nii');
glasser_left = cifti_struct_dense_extract_surface_data(glasser,'CORTEX_LEFT');
glasser_right = cifti_struct_dense_extract_surface_data(glasser,'CORTEX_RIGHT');
glasser_cdata = [glasser_left;glasser_right];

yeo_7networks = cifti_read('Yeo2011_7Networks_N1000.dlabel.nii');
yeo_7networks_cdata = yeo_7networks.cdata;

for i = 1:360
    idx = find(glasser_cdata == i);
    tbl = tabulate(yeo_7networks_cdata(idx));

    if any(tbl(:,1) == 0)
        tbl(tbl(:,1) == 0,:)=[];
    end

    [~,max_idx] = max(tbl(:,2));
    tbl_max(i,:) = tbl(max_idx,:);
    glasser_yeo(i) = tbl(max_idx,1);
end

net_label = glasser_yeo';
save([working_dir '7net_label_glasser360.mat'],'net_label')