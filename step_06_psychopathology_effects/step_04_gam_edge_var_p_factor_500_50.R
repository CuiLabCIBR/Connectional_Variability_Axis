library('visreg');
library('mgcv')
library('ggplot2')
library('tools')

rm(list = ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir,'step_06_psychopathology_effects/')

pfactor_term <- c('general')
atlas_term <- c('schaefer400')
reg_list <- c('no_reg_0y','no_reg_2y')

# pfactor_term = 'general'
# atlas_term = 'schaefer400'
# reg_term = 'no_reg_0y'

for (reg_term in reg_list)
{
  csv_dir <- file.path(working_dir,pfactor_term,atlas_term,reg_term)
  file_list <- dir(path = csv_dir, pattern = '^gam_pfactor_edge.*_500_50.csv$')
  
  N = length(file_list)
  p_anova <- matrix(0, N, 1);
  partial_R2 <- matrix(0, N, 1);
  t_gam <- matrix(0, N, 1);
  p_gam <- matrix(0, N, 1);
  
  i = 1
  for (file_name in file_list)
  {
    # file_name = file_list[1]
    gam.data <- read.csv(paste0(csv_dir,'/',file_name));
    gam.model <- gam(Variability ~ Pfactor + s(Age, k=3) + Gender + HeadMotion, method = "REML", data = gam.data);
    gam.model.results <- summary(gam.model)
    
    t_gam[i] <- gam.model.results$p.table[2,3]
    p_gam[i] <- gam.model.results$p.table[2,4]
    
    # reduced model without Pfactor
    gam.nullmodel <- gam(Variability ~ s(Age, k=3) + Gender + HeadMotion, method = "REML", data = gam.data);
    gam.nullmodel.results <- summary(gam.nullmodel)
    
    ##Full versus reduced model anova p-value
    p_anova[i] <- anova.gam(gam.nullmodel,gam.model,test='Chisq')$`Pr(>Chi)`[2]
    
    ##Full versus reduced model direction-dependent partial R squared
    ### effect size
    sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
    sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
    partial_R2[i] <- (sse.nullmodel - sse.model)/sse.nullmodel
    
    t <- gam.model.results$p.table[2,3];
    
    if (t < 0) {
      partial_R2[i] <- -partial_R2[i];
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
  csv_outpath <- paste0(working_dir,'results/axis/gam_edge_',pfactor_term,'_',atlas_term,'_',reg_term,'.csv')
  write.csv(results_table,csv_outpath)

}

