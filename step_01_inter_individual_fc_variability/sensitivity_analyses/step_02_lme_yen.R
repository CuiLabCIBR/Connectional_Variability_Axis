# The Rex package was used to estimate the inter- and intra-individual variability.
# lme_ICC_2wayM	ICC (ICC3) calculation using 2-way Mixed model.
# References:
# Xu, T., Kiar, G., Cho, J. W., Bridgeford, E. W., Nikolaidis, A., Vogelstein, J. T., & Milham, M. P. (2023). 
# ReX: an integrative tool for quantifying and optimizing measurement reliability for the study of individual differences. 
# Nature Methods, 20(7), 1025-1028.

library(ReX)
library(R.matlab)
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("rhdf5")
library(rhdf5)

rm(list = ls())

root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'

dataset_list <- c('yen')
atlas_list <- c('schaefer400')

# atlas <- 'schaefer400'
# dataset <- 'yen'

for (dataset in dataset_list) {
  
  for (atlas in atlas_list) {
    
    data_dir <- paste0(root_dir, 'data/fc/', atlas, '/')
    out_dir <- paste0(root_dir, 'data/fc_variability/', atlas, '/')
    
    subID <- as.matrix(unlist(readMat(paste0(data_dir, 'subID_', dataset, '.mat'))))
    session <- as.matrix(unlist(readMat(paste0(data_dir, 'session_', dataset, '.mat'))))
    
    dataset_atlas <- paste0(dataset, '_', atlas)
    
    fc_path <- paste0(data_dir, 'fc_', dataset_atlas, '.mat')
    h5ls(fc_path)
    data <- h5read(fc_path, '/fc_data')
    
    lme_results <- data.frame(lme_ICC_2wayM(data, subID, session))
    writeMat(paste0(out_dir, 'lme_', dataset_atlas, '.mat'), lme_results = lme_results)
    
  }
  
}
