library(R.matlab)
library(ggplot2)
library(dplyr)

rm(list = ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir, 'step_02_connectional_axis_of_fc_variability/')
data_dir <- paste0(root_dir, 'data/fc_variability/schaefer400/')

#####################################################################################
# plot the connectional axis from within-network to between-network
plot_connectional_axis <- function(FC_Var, plot_path){
  
  myPalette <- c("#974DA1FF", "#89B4D8FF", "#398E43FF", "#DE89FFFF", "#FFC87DFF", "#F68F9BFF")
  netOrder <- c("VS", "SM", "DA", "VA", "FP", "DM")
  netWidth <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
  
  for (i in 1:6) 
  {
    idx_within = i;
    variability_all = FC_Var[i, ]
    variability_within = FC_Var[i, i]
    idx_between = setdiff(1:6, idx_within)
    variability_between = FC_Var[i, idx_between]
    
    idx_between_sort = order(variability_between, decreasing = TRUE)
    netOrder_idx = c(idx_within, idx_between[idx_between_sort])
    netOrder_i = netOrder[netOrder_idx]
    variability_i = variability_all[netOrder_idx]
    df <- data.frame(network = netOrder_i, variability = variability_i)
    
    edgeColor = myPalette[netOrder_idx]
    fillColor = c(myPalette[i], "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF")
    
    p <- ggplot(data = df, aes(x = network, y = variability)) +
      geom_bar(stat = "identity", color = edgeColor, fill = fillColor, width = netWidth, lwd = 1.5) + 
      scale_x_discrete(limits = netOrder_i) +
      theme_classic()+labs(y = "", x = "") +
      theme(axis.text = element_text(size = 12, color = 'black'), axis.title = element_text(size = 12))+
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      theme(axis.line = element_line(linewidth = 0.6)) +
      scale_y_continuous(expand = c(0, 0))+
      coord_cartesian(ylim = c(0, 0.5))
    if (i != 6) 
    {
      p <- p + theme(axis.line.y = element_blank())
    }
    p
    
    outPath = paste0(plot_path, netOrder[i], '.png')
    ggsave(outPath, plot = p, width = 6, height = 6, units = "cm", dpi = 600)
  }
  
}

# HCP-D
plot_path <- paste0(working_dir, 'HCPD/')
if (!file.exists(plot_path)){
  dir.create(plot_path)
}

fc_data <- readMat(paste0(data_dir, 'hcpd_schaefer400_net.mat'))
fc_var <- fc_data$hcpd.schaefer400.net
plot_connectional_axis(fc_var, plot_path)

# HCP-YA
plot_path <- paste0(working_dir, 'HCP/')
if (!file.exists(plot_path)){
  dir.create(plot_path)
}

fc_data <- readMat(paste0(data_dir, 'hcp_schaefer400_net.mat'))
fc_var <- fc_data$hcp.schaefer400.net
plot_connectional_axis(fc_var, plot_path)
