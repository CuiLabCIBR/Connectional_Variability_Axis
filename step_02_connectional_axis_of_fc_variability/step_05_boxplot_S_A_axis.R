library(ggplot2)
library(dplyr)
library(bruceR)
library(rmatio)
library(psych)
library(effectsize)

rm(list = ls())
root_dir <- 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/'
working_dir <- paste0(root_dir, 'step_02_connectional_axis_of_fc_variability/')
data_dir <- paste0(root_dir, 'data/fc_variability/schaefer400/')

S_S <- as.matrix(read.mat(paste0(data_dir, 'hcpd_fc_variability_S_S.mat'))$S_S)
A_A <- as.matrix(read.mat(paste0(data_dir, 'hcpd_fc_variability_A_A.mat'))$A_A)
S_A <- as.matrix(read.mat(paste0(data_dir, 'hcpd_fc_variability_S_A.mat'))$S_A)

hcpd_cohens_d <- matrix(0, 3, 1)
hcpd_cohens_d[1] <- cohens_d(A_A, S_S)$Cohens_d
hcpd_cohens_d[2] <- cohens_d(A_A, S_A)$Cohens_d
hcpd_cohens_d[3] <- cohens_d(S_S, S_A)$Cohens_d

all_tabel <- data.frame(type = c(rep("S_S", nrow(S_S)), 
                                 rep("A_A", nrow(A_A)), 
                                 rep("S_A", nrow(S_A))), 
                        variability = c(S_S, A_A, S_A))

myPalette <- c("#F68F9B", "#B8C9E8", "#E6C6C8");
netOrder <- c("A_A", "S_A", "S_S");

# Change color by groups
p1 <- ggplot(all_tabel, aes(x = type, y = variability, fill = type)) + 
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = "#E3E3E3") + 
  theme_classic() + labs(y = "", x = "") + scale_x_discrete(limits = netOrder) + 
  coord_cartesian(ylim = c(0, 0.7)) + theme(legend.position = "none") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p1

ggsave(paste0(working_dir, 'boxplot_S_A_hcpd.png'), plot = p1, width = 3, height = 4, units = "cm", dpi = 600)

############HCP-YA############
S_S <- as.matrix(read.mat(paste0(data_dir, 'hcp_fc_variability_S_S.mat'))$S_S)
A_A <- as.matrix(read.mat(paste0(data_dir, 'hcp_fc_variability_A_A.mat'))$A_A)
S_A <- as.matrix(read.mat(paste0(data_dir, 'hcp_fc_variability_S_A.mat'))$S_A)

hcp_cohens_d <- matrix(0, 3, 1)
hcp_cohens_d[1] <- cohens_d(A_A, S_S)$Cohens_d
hcp_cohens_d[2] <- cohens_d(A_A, S_A)$Cohens_d
hcp_cohens_d[3] <- cohens_d(S_S, S_A)$Cohens_d

all_tabel <- data.frame(type = c(rep("S_S", nrow(S_S)), 
                                 rep("A_A", nrow(A_A)), 
                                 rep("S_A", nrow(S_A))), 
                        variability = c(S_S, A_A, S_A))

myPalette <- c("#F68F9B", "#B8C9E8", "#E6C6C8");
netOrder <- c("A_A", "S_A", "S_S");

# Change color by groups
p2 <- ggplot(all_tabel, aes(x = type, y = variability, fill = type)) + 
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = "#E3E3E3") + 
  theme_classic() + labs(y = "", x = "") + scale_x_discrete(limits = netOrder) + 
  coord_cartesian(ylim = c(0, 0.7)) + theme(legend.position = "none") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p2

ggsave(paste0(working_dir, 'boxplot_S_A_hcp.png'), plot = p2, width = 3, height = 4, units = "cm", dpi = 600)


###
cohens_d <- data.frame(hcpd = hcpd_cohens_d, hcp = hcp_cohens_d)
write.csv(cohens_d, paste0(working_dir, 'S_A_cohens_d.csv'), row.names = FALSE)
