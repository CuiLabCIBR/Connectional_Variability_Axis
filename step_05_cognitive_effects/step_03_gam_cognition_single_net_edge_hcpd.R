library(visreg)
library(mgcv)
library(ggplot2)
library(tools)

rm(list = ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir, 'step_05_cognitive_effects/')

csv_dir <- working_dir
file_list <- dir(path = csv_dir, pattern = '^gam_cog_hcpd_edge.*.csv$')

N = length(file_list)
p_anova <- matrix(0, N, 1)
partial_R2 <- matrix(0, N, 1)
t_gam <- matrix(0, N, 1)
p_gam <- matrix(0, N, 1)

i = 1
for (file_name in file_list)
{
  # file_name = file_list[1]
  gam.data <- read.csv(paste0(csv_dir, '/', file_name))
  gam.model <- gam(Variability ~ Cognition + s(Age, k = 3) +  Sex + HeadMotion, method = "REML", data = gam.data)
  gam.model.results <- summary(gam.model)
  
  t_gam[i] <- gam.model.results$p.table[2, 3]
  p_gam[i] <- gam.model.results$p.table[2, 4]
  
  # reduced model without cognition
  gam.nullmodel <- gam(Variability ~ s(Age, k = 3) +  Sex + HeadMotion, method = "REML", data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  ##Full versus reduced model anova p-value
  p_anova[i] <- anova.gam(gam.nullmodel, gam.model, test = 'Chisq')$`Pr(>Chi)`[2]
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partial_R2[i] <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  if (t_gam[i] < 0) {
    partial_R2[i] <- -partial_R2[i]
  }
  
  i = i+1
}

results_table <- data.frame(
  matrices = file_list, 
  t_gam = t_gam, 
  p_gam = p_gam, 
  p_gam_sig = as.double(p_gam < 0.05), 
  partial_R2 = partial_R2, 
  p_anova = p_anova, 
  p_anova_sig = as.double(p_anova < 0.05)
)

# results_table <- t(results_table)
csv_outpath <- paste0(working_dir, 'gam_cognition_edge_stats_hcpd.csv')
write.csv(results_table, csv_outpath)
