library(visreg)
library(mgcv)
library(ggplot2)
library(tools)

rm(list = ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir, 'step_05_cognitive_effects/sensitivity_analyses/')

csv_dir <- working_dir
file_list <- dir(path = csv_dir, pattern = '^gam_cog_.*.csv$')

N = length(file_list)
p_anova <- matrix(0, N, 1)
partial_R2 <- matrix(0, N, 1)
t_gam <- matrix(0, N, 1)
p_gam <- matrix(0, N, 1)

i = 1
for (file_name in file_list)
{
  # file_name = file_list[1]
  print(file_name)
  gam.data <- read.csv(paste0(csv_dir, '/', file_name))
  gam.model <- gam(Slope ~ Cognition + s(Age, k = 3) +  Sex + HeadMotion, method = "REML", data = gam.data)
  gam.model.results <- summary(gam.model)
  
  t_gam[i] <- gam.model.results$p.table[2, 3]
  p_gam[i] <- gam.model.results$p.table[2, 4]
  
  # reduced model without cognition
  gam.nullmodel <- gam(Slope ~ s(Age, k = 3) +  Sex + HeadMotion, method = "REML", data = gam.data)
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
  
  #############################
  # do not plot for sliding window results
  if (!grepl("window", file_name, ignore.case = TRUE))
  {
    v <- visreg(gam.model, "Cognition", plot = FALSE)
    color_min = min(v$res$visregRes)
    color_max = max(v$res$visregRes)
    color_mid = median(v$res$visregRes)
    
    Fig <- plot(v, xlab = "Fluid cognition composite score", ylab = "Connectional axis slope", 
                line.par = list(col = '#7499C2'), fill = list(fill = '#D9E2EC'), gg = TRUE,  rug = FALSE)
    Fig <- Fig + geom_point(data = v$res, aes(x = Cognition, y = visregRes, color = visregRes), size = 3, alpha = 1, shape = 19)+
      scale_color_gradient2(low = "#2473B5", high = "#CF1A1D", mid = "#F6FBFF", 
                            midpoint = color_mid, limit = c(color_min, color_max))+
      theme_classic() +
      theme(axis.text = element_text(size = 18, color = 'black'), axis.title = element_text(size = 18), aspect.ratio = 0.8) +
      theme(legend.position = "none") +
      scale_y_continuous(labels = number_format(accuracy = 0.01)) 
    
    if (grepl("hcpd", file_name, ignore.case = TRUE)) {
      Fig <- Fig + scale_x_continuous(expand = c(0.01, 0), limits = c(88, 132), breaks = seq(90, 130, by = 10))
    } else if (grepl("hcp", file_name, ignore.case = TRUE)) {
      Fig <- Fig + scale_x_continuous(expand = c(0.01, 0), limits = c(100, 135), breaks = seq(100, 130, by = 10))
    } 
    else {
      Fig <- Fig + scale_x_continuous(expand = c(0.01, 0), limits = c(80, 125), breaks = seq(80, 120, by = 10))
    }
    
    
    Fig
    
    ggsave(paste0(csv_dir, '/', file_path_sans_ext(file_name), '.png'), plot = Fig, width = 12, height = 10, units = "cm", dpi = 1200)
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
csv_outpath <- paste0(working_dir, 'gam_cognition_stats.csv')
write.csv(results_table, csv_outpath)

GAM_results <- cbind(round(partial_R2, digits = 2), p_anova)
GAM_results
