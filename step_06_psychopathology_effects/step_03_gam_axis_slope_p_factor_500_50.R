library('visreg');
library('mgcv')
library('ggplot2')
library('tools')

rm(list = ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir,'step_06_psychopathology_effects/')

pfactor_term <- c('general')

atlas_list <- c(rep('schaefer400', 6), rep('glasser360', 2))

reg_list <- c('no_reg_0y','no_reg_2y','reg_fc_0y','reg_fc_2y',
              'reg_dist_0y','reg_dist_2y','no_reg_0y','no_reg_2y')

x_lable <- c(rep(c('General psychopathology\n(baseline)','General psychopathology\n(2-year follow-up)'),4))

N = length(atlas_list)
p_anova <- matrix(0, N, 1)
partial_R2 <- matrix(0, N, 1)
t_gam <- matrix(0, N, 1)
p_gam <- matrix(0, N, 1)

for (i in 1:N)
{
  csv_dir <- file.path(working_dir,pfactor_term,atlas_list[i],reg_list[i])
  file_name = 'gam_pfactor_axis_slope_500_50.csv'
  gam.data <- read.csv(file.path(csv_dir,file_name))
  gam.model <- gam(Slope ~ Pfactor + s(Age, k=3) + Gender + HeadMotion, method = "REML", data = gam.data)
  gam.model.results <- summary(gam.model)
  
  t_gam[i] <- gam.model.results$p.table[2,3]
  p_gam[i] <- gam.model.results$p.table[2,4]
  
  # reduced model without Pfactor
  gam.nullmodel <- gam(Slope ~ s(Age, k=3) + Gender + HeadMotion, method = "REML", data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  ##Full versus reduced model anova p-value
  p_anova[i] <- anova.gam(gam.nullmodel,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partial_R2[i] <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  t <- gam.model.results$p.table[2,3]
  
  if (t < 0) {
    partial_R2[i] <- -partial_R2[i]
  }
  
  #############################
  v <- visreg(gam.model, "Pfactor", plot=FALSE)
  color_min = min(v$res$visregRes)
  color_max = max(v$res$visregRes)
  color_mid = median(v$res$visregRes)
  
  Fig <- plot(v, xlab = x_lable[i], ylab = 'Connectional axis slope',
              line.par = list(col = '#7499C2'), fill = list(fill = '#D9E2EC'), gg=TRUE,  rug=FALSE)
  
  Fig <- Fig + geom_point(data=v$res, aes(x=Pfactor, y=visregRes, color=visregRes), size=3, alpha = 1, shape=19)+
    scale_color_gradient2(low = "#2473B5", high = "#CF1A1D", mid = "#F6FBFF",
                          midpoint = color_mid, limit = c(color_min,color_max))+
    theme_classic() +
    theme(axis.text=element_text(size=16, color='black'), axis.title=element_text(size=16), aspect.ratio = 0.7) +
    theme(legend.position="none")

  switch(
    i,
    "1" = Fig <- Fig + scale_y_continuous(expand = c(0, 0),limits = c(0.257, 0.275), breaks = seq(0.260, 0.272, by = 0.004)),
    "2" = Fig <- Fig + scale_y_continuous(expand = c(0, 0),limits = c(0.258, 0.271), breaks = seq(0.260, 0.269, by = 0.003)),
    "3" = Fig <- Fig + scale_y_continuous(expand = c(0, 0),limits = c(0.254, 0.276), breaks = seq(0.258, 0.273, by = 0.005)),
    {Fig <- Fig}
  )
   
  Fig <- Fig + scale_x_continuous(expand = c(0.01, 0),limits = c(-1.3, 1.52), breaks = seq(-1, 1.5, by = 0.5))
  
  Fig
  
  Fig_path <- paste0(working_dir,'results/slope/',file_path_sans_ext(file_name),'_',atlas_list[i],'_',reg_list[i], '.png')
  ggsave(Fig_path,plot=Fig,width = 13,height = 11,units = "cm",dpi = 1200)
}

results_table <- data.frame(
  matrices = paste(atlas_list, reg_list, sep = "_"),
  t_gam = t_gam,
  p_gam = p_gam,
  p_gam_sig = as.double(p_gam < 0.05),
  partial_R2 = partial_R2,
  p_anova = p_anova,
  p_anova_sig = as.double(p_anova < 0.05)
)

csv_outpath <- paste0(working_dir,'results/slope/gam_results_pfactor_axis_slope_500_50.csv')
write.csv(results_table,csv_outpath)
