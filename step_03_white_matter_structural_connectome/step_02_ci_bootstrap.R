library(bootcorci)
library(tools)

rm(list = ls())
working_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Gradient/step_03_white_matter_structural_connectome/'

file_list <- dir(path = working_dir, pattern = "^hcp.*\\.csv$")
corr_results <- data.frame()

# reproducability of results
set.seed(716)

for (file_name in file_list)
{
  data <- read.csv(paste0(working_dir,file_name));
  
  res <- corci(data$mat_a, data$mat_b, method = "spearman", nboot = 1000, 
               alpha = 0.05, alternative = "two.sided")
  
  corr_results <- rbind(corr_results,c(file_path_sans_ext(file_name),res$estimate,res$conf.int))
}
colnames(corr_results) <- c('file','rho','ci_low','ci_up')
# write.table(corr_results,paste0(working_dir,'corr_results.csv'),sep = ",",col.names = TRUE,row.names = FALSE)