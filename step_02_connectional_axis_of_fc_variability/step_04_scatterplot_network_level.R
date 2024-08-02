library(ggplot2)

rm(list=ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir,'step_02_connectional_axis_of_fc_variability/')
FC_Var_net <- import(paste0(working_dir,'hcpd_hcp_fc_variability_schaefer400_net.csv'))

r <- cor(FC_Var_net$mat_b,FC_Var_net$mat_a,method = "spearman")

p<-ggplot(data = FC_Var_net, aes(x = mat_b, y = mat_a))+geom_point(size=2.5,color="#7CAED8",alpha=0.8)+
  geom_smooth(method = lm, se=FALSE, color = "#4682B4", size=2) +
  theme_classic()+labs(y = "HCP-D \n FC variability", x = "FC variability \n HCP-YA") +
  theme(axis.text=element_text(size=20, color='black'), axis.title=element_text(size=20), aspect.ratio = 1.3)+
  theme(legend.position="none") + scale_y_continuous(breaks = seq(0.1, 0.5, by = 0.1), limits = c(0.15, 0.51)) +
  scale_x_continuous(breaks = seq(0.1, 0.5, by = 0.1), limits = c(0.15, 0.51))
p
ggsave(paste0(working_dir,'scatter_plot_hcpd_hcp_schaefer400_net.png'),plot=p,width = 12,height = 12,units = "cm",dpi = 600)
