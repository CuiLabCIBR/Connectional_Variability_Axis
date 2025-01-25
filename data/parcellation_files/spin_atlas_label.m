clear
clc
addpath(genpath('.../ENIGMA/matlab'))

%% schaefer_400
lsphere = SurfStatReadSurf1('fsa5_sphere_lh');
rsphere = SurfStatReadSurf1('fsa5_sphere_rh');

lh_centroid = centroid_extraction_sphere(lsphere.coord.', 'fsa5_lh_schaefer_400.annot');
rh_centroid = centroid_extraction_sphere(rsphere.coord.', 'fsa5_rh_schaefer_400.annot');

perm_id = rotate_parcellation(lh_centroid, rh_centroid, 10000);

save('F:/Cui_Lab/Projects/Connectional_Variability_Axis/data/parcellation_files/perm_id_schaefer400.mat','perm_id')

%% glasser_360
lsphere = SurfStatReadSurf1('fsa5_sphere_lh');
rsphere = SurfStatReadSurf1('fsa5_sphere_rh');

lh_centroid = centroid_extraction_sphere(lsphere.coord.', 'fsa5_lh_glasser_360.annot');
rh_centroid = centroid_extraction_sphere(rsphere.coord.', 'fsa5_rh_glasser_360.annot');

perm_id = rotate_parcellation(lh_centroid, rh_centroid, 10000);

save('F:/Cui_Lab/Projects/Connectional_Variability_Axis/data/parcellation_files/perm_id_glasser360.mat','perm_id')
