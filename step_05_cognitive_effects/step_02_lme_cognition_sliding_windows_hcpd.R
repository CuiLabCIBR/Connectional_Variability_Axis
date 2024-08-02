# 52 windows

get_fc_variability <- function(i) {
  
  library(ReX)
  library(R.matlab)
  
  root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
  working_dir <- paste0(root_dir,'data/fc/schaefer400/cognitive_effects/hcpd/')
  out_dir <- paste0(root_dir,'data/fc_variability/schaefer400/cognitive_effects/hcpd/')
  
  subID <-  as.matrix(unlist(readMat(paste0(working_dir,'subID_hcpd.mat'))))
  session <- as.matrix(unlist(readMat(paste0(working_dir,'session_hcpd.mat'))))
  
  data <- readMat(paste0(working_dir,'fc_', as.character(i), '.mat'))
  data <- data$fc
  
  rm(df_lme_all)
  df_lme_all <- data.frame(lme_ICC_2wayM(data, subID, session))
  
  writeMat(paste0(out_dir,'lme_',as.character(i),'.mat'),lme = df_lme_all)  
}


#################### run the lme in parallel ####################
library(foreach)
library(doParallel)
num_cores <- 6
cluster <- makeCluster(num_cores)
registerDoParallel(cluster)

foreach(i = 1:52) %dopar% {
  get_fc_variability(i)
}

#################### run the lme in serial ####################
# library(ReX)
# library(R.matlab)
# 
# root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
# working_dir <- paste0(root_dir,'data/fc/schaefer400/cognitive_effects/hcpd/')
# out_dir <- paste0(root_dir,'data/fc_variability/schaefer400/cognitive_effects/hcpd/')
# 
# subID <-  as.matrix(unlist(readMat(paste0(working_dir,'subID_hcpd.mat'))))
# session <- as.matrix(unlist(readMat(paste0(working_dir,'session_hcpd.mat'))))
# 
# for (i in 1:57)
# {
#   data <- readMat(paste0(working_dir,'fc_', as.character(i), '.mat'))
#   data <- data$fc
# 
#   rm(df_lme_all)
#   df_lme_all <- data.frame(lme_ICC_2wayM(data, subID, session))
# 
#   writeMat(paste0(out_dir,'lme_',as.character(i),'.mat'),lme = df_lme_all)
# }
