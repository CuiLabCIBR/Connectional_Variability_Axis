library(R.matlab)
library(ggplot2)
library(dplyr)
library(openxlsx)

rm(list = ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir, 'step_02_connectional_axis_of_fc_variability/sensitivity_analyses/')

netOrder <- c("DA", "DM", "FP", "SM", "VA", "VS");
myPalette <- c("#398E43FF", "#F68F9BFF", "#FFC87DFF", "#89B4D8FF", "#DE89FFFF", "#974DA1FF");

NetInfo <- read.xlsx(paste0(working_dir, 'connectional_variability_axis_yen.xlsx'))
axis <- NetInfo$Value;
Net_All = NetInfo$Net_All;
Net_Single = NetInfo$Net_Single;
df <- data.frame(Net_All = Net_All, axis = axis, Net_Single = Net_Single);

p <- ggplot(data = df, aes(x = Net_All, y = axis, fill = Net_Single)) + 
  geom_col(width = 1.5, lwd = 1) + scale_x_discrete(expand = c(0.03, 0), limits = Net_All) +
  theme_classic() + labs(y = "", x = "") + scale_fill_manual(values = myPalette) +
  theme(axis.text = element_text(size = 12, color = 'black'), axis.title = element_text(size = 12)) +
  theme(axis.line = element_line(size = 0.3)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.6)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(legend.position = "none")
p

ggsave(paste0(working_dir, 'connectional_variability_axis_yen.png'), plot = p, width = 13, height = 6, units = "cm", dpi = 1200)
