# This script is used to estimate the inter- and intra-individual variability
# install.packages("devtools")
# devtools::install_github("TingsterX/ReX")
library(ReX)
library(R.matlab)

rm(list = ls())

root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Gradient/'
data_dir <- paste0(root_dir,'data/fc/schaefer400/')
out_dir <- paste0(root_dir,'step_01_inter_individual_fc_variability/hcpd/')

subID <-  as.matrix(unlist(readMat(paste0(data_dir,'subID_hcpd.mat'))))
session <- as.matrix(unlist(readMat(paste0(data_dir,'session_hcpd.mat'))))

data <- readMat(paste0(data_dir,'hcpd_fc_schaefer400.mat'))
data <- data$hcpd.fc

lme_hcpd_schaefer400 <- data.frame(lme_ICC_2wayM(data, subID, session))
writeMat(paste0(out_dir,'lme_hcpd_schaefer400.mat'),lme_hcpd_schaefer400 = lme_hcpd_schaefer400)