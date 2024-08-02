library(ggplot2)
library(dplyr)
library(openxlsx)
library(tools)

rm(list = ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir,'step_06_psychopathology_effects/results/axis/')

netOrder <- c("DA","DM","FP","SM","VA","VS");
myPalette <- c("#398E43FF", "#F68F9BFF", "#FFC87DFF", "#89B4D8FF", "#DE89FFFF", "#974DA1FF");

file_list <- c('edge_fc_var_axis_schaefer400_no_reg_0y.xlsx','edge_fc_var_axis_schaefer400_no_reg_2y.xlsx')

i = 1

for (file_name in file_list) {
  NetInfo <- read.xlsx(paste0(working_dir,file_name))
  axis <- NetInfo$Value;
  Net_All = NetInfo$Net_All;
  Net_Single = NetInfo$Net_Single;
  df <- data.frame(Net_All=Net_All, axis=axis, Net_Single=Net_Single);
  
  p <- ggplot(data=df, aes(x=Net_All, y=axis, fill=Net_Single)) + 
    geom_col(width=1.5, lwd=1) + scale_x_discrete(expand=c(0.03, 0),limits=Net_All) +
    theme_classic() + labs(y = "",x = "") + scale_fill_manual(values=myPalette) +
    theme(axis.text=element_text(size=9, color='black'), axis.title=element_text(size=12)) +
    theme(axis.line=element_line(size=0.3), aspect.ratio = 0.6) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    # scale_y_continuous(expand = c(0, 0), limits = c(-0.031, 0.01)) +
    # theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme(legend.position="none")
  p
  
  switch(
    i,
    "1" = p <- p + scale_y_continuous(expand = c(0, 0),limits = c(-7, 4), breaks = seq(-6, 3, by = 3)),
    "2" = p <- p + scale_y_continuous(expand = c(0, 0),limits = c(-5, 4.5), breaks = seq(-4, 4, by = 2))
  )
  
  ggsave(paste0(working_dir,file_path_sans_ext(file_name),'.png'),plot=p,width = 13,height = 6,units = "cm",dpi = 1200)
  
  i = i+1
}

