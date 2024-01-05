# This script is used to estimate the inter- and intra-individual variability

library(ReX)
library(R.matlab)

rm(list = ls())

root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Gradient/'
data_dir <- paste0(root_dir,'data/fc/schaefer400/')
out_dir <- paste0(root_dir,'step_01_inter_individual_fc_variability/hcp/')

subID <-  as.matrix(unlist(readMat(paste0(data_dir,'subID_hcp.mat'))))
session <- as.matrix(unlist(readMat(paste0(data_dir,'session_hcp.mat'))))

data <- readMat(paste0(data_dir,'hcp_fc_schaefer400.mat'))
data <- data$hcp.fc

lme_hcp_schaefer400 <- data.frame(lme_ICC_2wayM(data, subID, session))
writeMat(paste0(out_dir,'lme_hcp_schaefer400.mat'),lme_hcp_schaefer400 = lme_hcp_schaefer400)