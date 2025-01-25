function ptseries = parcellate_cifti(dtseries, atlas)

% Inputs:
% dtseries: the dtseries.nii cifti.cdata or filepath
% atlas: the dlabel.nii cifti.cdata or filepath
% Outputs:
% ptseries: the mean timeseries parcellated by the atlas

% addpath(genpath('/ibmgpfs/cuizaixu_lab/yanghang/code/cifti-matlab-master/')) % change to your cifti path

% if the dtseries is a filepath but not exist
if (isa(dtseries, 'char') || isa(dtseries, 'string')) && exist(dtseries,'file') ~= 2 
    msg = ['The file: ' dtseries ' does not exist!']; 
    error(msg)
end

% if the dtseries is a filepath
if (isa(dtseries, 'char') || isa(dtseries, 'string')) && exist(dtseries,'file') == 2 
    dtseries_cifti = cifti_read(dtseries);
    dtseries_cdata = dtseries_cifti.cdata;
% if the dtseries is a cifti.cdata
else                           
    dtseries_cdata = dtseries;
end

% if the atlas is a filepath but not exist
if (isa(atlas, 'char') || isa(atlas, 'string')) && exist(atlas,'file') ~= 2
    msg = ['The file: ' atlas ' does not exist!'];  
    error(msg)
end

% if the atlas is a filepath
if (isa(atlas, 'char') || isa(atlas, 'string')) && exist(atlas,'file') == 2 
    atlas_cifti = cifti_read(atlas);
    atlas_cdata = atlas_cifti.cdata;
% if the atlas is a cifti.cdata
else                           
    atlas_cdata = atlas;
end

%%
if size(dtseries_cdata, 1) ~= size(atlas_cdata, 1)
    
    dtseries_left = cifti_struct_dense_extract_surface_data(dtseries_cifti, 'CORTEX_LEFT');
    dtseries_right = cifti_struct_dense_extract_surface_data(dtseries_cifti, 'CORTEX_RIGHT');
    dtseries_cdata = [dtseries_left; dtseries_right];
    
    atlas_left = cifti_struct_dense_extract_surface_data(atlas_cifti, 'CORTEX_LEFT');
    atlas_right = cifti_struct_dense_extract_surface_data(atlas_cifti, 'CORTEX_RIGHT');
    atlas_cdata = [atlas_left; atlas_right];
    
    if size(dtseries_cdata, 1) ~= size(atlas_cdata, 1)
        msg = 'The vertex number does match between the data and atlas!';
        error(msg)
    end
    
end

N = length(unique(atlas_cdata)) - 1; % Number of ROIs, exclude id = 0;
T = size(dtseries_cdata, 2);         % Number of timepoints
ptseries = nan(N, T);

for roi_i = 1:N
    idx = find(atlas_cdata == roi_i);
    ptseries(roi_i, :) = mean(dtseries_cdata(idx,:));
end
