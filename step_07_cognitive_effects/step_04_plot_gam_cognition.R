library('visreg');
library('mgcv')
library('ggplot2')
library('effectsize')

rm(list = ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir,'step_07_cognitive_effects/')

t_gam <- matrix(0, 2, 1);
p_gam <- matrix(0, 2, 1);
r_gam <- matrix(0, 2, 1);
p_anova <- matrix(0, 2, 1);
partial_R2 <- matrix(0, 2, 1);

file_list <- c('gam_cog_hcpd_slope_net','gam_cog_hcp_slope_net')
y_name <- c('Variability axis slope','Variability axis slope')
i = 1

for (file_name in file_list)
{
  gam.data <- read.csv(paste0(working_dir,file_name,'.csv'));
  gam.model <- gam(Axis ~ Cognition + s(Age, k=3) + Gender + HeadMotion, method = "REML", data = gam.data);
  gam.model.results <- summary(gam.model)
  
  t_gam[i] <- gam.model.results$p.table[2,3]
  p_gam[i] <- gam.model.results$p.table[2,4]
  
  # convert T to partial r
  n <- nrow(gam.data)
  df_error <- n - 3 - summary(gam.model)$edf
  r_gam[i] <- t_to_r(t_gam[i],df_error)$r
  
  # reduced model without cognition
  gam.nullmodel <- gam(Axis ~ s(Age, k=3) + Gender + HeadMotion, method = "REML", data = gam.data);
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
  
  #############################
  
  v <- visreg(gam.model, "Cognition", plot=FALSE)
  color_min = min(v$res$visregRes)
  color_max = max(v$res$visregRes)
  color_mid = median(v$res$visregRes)
  
  Fig <- plot(v, xlab = "Fluid cognition composite score", ylab = y_name[i],
              line.par = list(col = '#7499C2'), fill = list(fill = '#D9E2EC'), gg=TRUE,  rug=FALSE)
  Fig <- Fig + geom_point(data=v$res, aes(x=Cognition, y=visregRes, color=visregRes), size=3, alpha = 1, shape=19)+
    scale_color_gradient2(low = "#2473B5", high = "#CF1A1D", mid = "#F6FBFF",
                          midpoint = color_mid, limit = c(color_min,color_max))+
    theme_classic() +
    theme(axis.text=element_text(size=16, color='black'), axis.title=element_text(size=16), aspect.ratio = 0.8) +
    theme(legend.position="none")
  
  Fig <- switch(
    i,
    "1" = Fig + scale_y_continuous(expand = c(0, 0),limits = c(0.255, 0.325), breaks = seq(0.26, 0.32, by = 0.02))+
      scale_x_continuous(expand = c(0.01, 0),limits = c(88, 132), breaks = seq(90, 130, by = 10)),
    "2" = Fig + scale_x_continuous(expand = c(0.01, 0),limits = c(100, 135), breaks = seq(100, 130, by = 10)) +
      scale_y_continuous(expand = c(0, 0),limits = c(0.205, 0.255), breaks = seq(0.21, 0.25, by = 0.01))
  )
  
  Fig
  
  ggsave(paste0(working_dir,file_name, '.png'),plot=Fig,width = 13,height = 11,units = "cm",dpi = 1200)
  
  i = i+1
}

########
results_table <- data.frame(
  matrices = y_name,
  t_gam = t_gam,
  p_gam = p_gam,
  r_gam = r_gam,
  partial_R2 = partial_R2,
  p_anova = p_anova
)

results_table <- t(results_table)
write.table(results_table,paste0(working_dir,'gam_cognition.csv'),sep = ",",col.names = FALSE)

GAM_results <- cbind(round(partial_R2,digits = 2),p_anova)
GAM_results
