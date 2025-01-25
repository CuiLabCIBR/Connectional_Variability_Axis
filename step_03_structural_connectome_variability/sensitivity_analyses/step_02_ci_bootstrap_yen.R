library(bootcorci)
library(tools)

rm(list = ls())
working_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/step_03_structural_connectome_variability/sensitivity_analyses/'

file_list <- dir(path = working_dir, pattern = "^yen.*\\.csv$")
corr_results <- data.frame()

# reproducibility of results
set.seed(410)

for (file_name in file_list)
{
  data <- read.csv(paste0(working_dir, file_name));
  
  res <- corci(data$mat_a, data$mat_b, method = "spearman", nboot = 1000, 
               alpha = 0.05, alternative = "two.sided")
  
  corr_results <- rbind(corr_results, c(file_path_sans_ext(file_name), res$estimate, res$conf.int))
}
colnames(corr_results) <- c('file', 'rho', 'ci_low', 'ci_up')
write.table(corr_results, paste0(working_dir, 'bootstrap_corr_results_yen.csv'), sep = ",", col.names = TRUE, row.names = FALSE)