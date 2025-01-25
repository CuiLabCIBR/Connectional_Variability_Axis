library(visreg)
library(mgcv)
library(ggplot2)
library(tools)

rm(list = ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir,'step_04_developmental_effects/sensitivity_analyses/')

file_list <- dir(path = working_dir, pattern = '^gam_age_axis_slope.*.csv$')
N <- length(file_list)

p_value <- matrix( 0, N, 1)
partial_R2 <- matrix( 0, N, 1)

i = 1

for (file_name in file_list)
{
  gam.data <- read.csv(paste0(working_dir, file_name))
  gam.model <- gam(Slope ~ s(Age, k = 3) + Sex + HeadMotion, method = "REML", data = gam.data)
  gam.model.results <- summary(gam.model)
  
  gam.nullmodel <- gam(Slope ~ Sex + HeadMotion, method = "REML", data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  ##Full versus reduced model anova p-value
  p_value[i] <- anova.gam(gam.nullmodel, gam.model, test = 'Chisq')$`Pr(>Chi)`[2]
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partial_R2[i] <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  lm.model <- lm(Slope ~ Age + HeadMotion, data = gam.data);
  t <- summary(lm.model)$coefficients[2, 3]
  
  if (t < 0) {
    partial_R2[i] <- -partial_R2[i]
  }
  
  ###################
  v <- visreg(gam.model, "Age", plot = FALSE)
  
  color_min = min(v$res$visregRes)
  color_max = max(v$res$visregRes)
  color_mid = mean(v$res$visregRes)
  
  Fig <- plot(v, xlab = "Age (years)", ylab = "Connectional axis slope", 
              line.par = list(col = '#7499C2'), fill = list(fill = '#D9E2EC'), gg = TRUE,  rug = FALSE) + 
    theme_classic() +
    geom_point(data = v$res, aes(x = Age, y = visregRes, color = visregRes), size = 1.5, alpha = 1, shape = 19)+
    scale_color_gradient2(low = "#2473B5", high = "#CF1A1D", mid = "#F6FBFF", 
                          midpoint = color_mid, limit = c(color_min, color_max))+
    theme(legend.position = "none") +
    theme(axis.text = element_text(size = 18, color = 'black'), axis.title = element_text(size = 18), aspect.ratio = 0.8) +
    scale_x_continuous(expand = c(0.01, 0), limits = c(8.5, 22), breaks = seq(9, 21, by = 2))
  
  Fig
  
  if (file_name == "gam_age_axis_slope_brainproject.csv") 
  {
    Fig <- Fig + scale_x_continuous(expand = c(0.01, 0), limits = c(7.5, 20.5), breaks = seq(8, 20, by = 2)) +
      scale_y_continuous(expand = c(0.01, 0), limits = c(0.235, 0.345), breaks = seq(0.24, 0.34, by = 0.02))
  }
  
  ggsave(paste0(working_dir, file_path_sans_ext(file_name), '.png'), plot = Fig, width = 12, height = 10, units = "cm", dpi = 1200)
  
  i = i+1
}

results_table <- data.frame(
  matrices = file_list, 
  partial_R2 = partial_R2, 
  p_value = p_value
)

write.table(results_table, paste0(working_dir, 'gam_development_stats.csv'), sep = ",", row.names = FALSE)

GAM_results <- cbind(round(partial_R2, digits = 2), p_value)
GAM_results
