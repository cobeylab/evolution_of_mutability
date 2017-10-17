# Produces figure for the manuscript showing the results of the simulations on MCC trees
# Uses the following MCMC chains: CH103_con_run1a, CH103L_con_run1a, VRC26int_con_run1a, VRC26L_con_run1a, VRC01_01_log_run1a, VRC01_13_log_run1a, VRC01_19_log_run1a,

library('ggplot2')
library('cowplot')
#library('reshape')
#library('gridExtra')
library('grid')
#library('coda')
library('lattice')
#library('gridBase')

results_directory <- '../results/relative_mutability/observed_lineages/'

# ======= GET DATAFRAMES WITH OBSERVED MCC MUTABILITY RESULTS FOR ALL LINEAGES

CH103_observed_dataframe <- read.table(paste(results_directory, 'CH103_constant/CH103_con_run1a_observed_mutability_MCC.csv', sep = ''), header = T, sep = ',')
CH103L_observed_dataframe <- read.table(paste(results_directory, 'CH103L_constant/CH103L_con_run1a_observed_mutability_MCC.csv', sep = ''), header = T, sep = ',')
VRC26int_observed_dataframe <- read.table(paste(results_directory, 'VRC26int_constant/VRC26int_con_run1a_observed_mutability_MCC.csv', sep = ''), header = T, sep = ',')
VRC26L_observed_dataframe <- read.table(paste(results_directory, 'VRC26L_constant/VRC26L_con_run1a_observed_mutability_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_01_observed_dataframe <- read.table(paste(results_directory, 'VRC01_01_logistic/VRC01_01_log_run1a_observed_mutability_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_13_observed_dataframe <- read.table(paste(results_directory, 'VRC01_13_logistic/VRC01_13_log_run1a_observed_mutability_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_19_observed_dataframe <- read.table(paste(results_directory, 'VRC01_19_logistic/VRC01_19_log_run1a_observed_mutability_MCC.csv', sep = ''), header = T, sep = ',')

observed_dataframe_list <- list('CH103' = CH103_observed_dataframe, 'CH103L' = CH103L_observed_dataframe,
                                'VRC26int' = VRC26int_observed_dataframe, 'VRC26L' = VRC26L_observed_dataframe,
                                'VRC01_01' = VRC01_01_observed_dataframe, 'VRC01_13' = VRC01_13_observed_dataframe,
                                'VRC01_19' = VRC01_19_observed_dataframe)


# ====== COMBINED DATAFRAME WITH RESULTS FOR ALL LINEAGES =======
lineage_vector <- c()
logS5F_CDR <- c()
logS5F_FR <- c()
time_from_root <- c()
distance_from_root <- c()
for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  clone_dataframe_obs <- observed_dataframe_list[[clone]] 
  
  logS5F_FR <- c(logS5F_FR,clone_dataframe_obs$observed_logS5F_FR)
  logS5F_CDR <- c(logS5F_CDR,clone_dataframe_obs$observed_logS5F_CDR)
  
  time_from_root <- c(time_from_root, clone_dataframe_obs$time_from_root)
  distance_from_root <- c(distance_from_root, clone_dataframe_obs$distance_to_root)
  lineage_vector <- c(lineage_vector, nrow(clone_dataframe_obs))
}

lineage_vector <- factor(lineage_vector, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))

combined_dataframe <- data.frame(lineage=lineage_vector,time_from_root,distance_from_root,
                                 logS5F_FR, logS5F_CDR, CDR_FR_logS5F_difference = logS5F_CDR - logS5F_FR)

# ================ IMPORTING GGPLOT2 PARAMETERS ===============
source('ggplot_parameters.R')

# =====================================================================================================================
# ================================================== PLOTS ============================================================
# =====================================================================================================================

# ======== PLOT OF CDR/FR S5F difference VS distance.
pl_ratio_vs_distance <- ggplot(combined_dataframe, 
                               aes(x = distance_from_root, y = CDR_FR_logS5F_difference)) +
  geom_hline(yintercept=0,linetype = 2) +
  theme_bw() +
  #theme_classic() + 
  ylab('Relative difference between CDR and FR mutability') + 
  xlab('Genetic distance from root') +
  theme(axis.title.y = element_text(size = axis_title_size,
                                    margin = margin(0,ylab_distance,0,0)),
        axis.title.x = element_text(size = axis_title_size,
                                    margin = margin(xlab_distance,0,0,0)),
        axis.text.x = element_text(size = x_axis_text_size, angle = 0, vjust = 0.6),
        axis.text.y = element_text(size = y_axis_text_size),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        plot.margin = unit(c(1, 1, 1, 1),"cm"),
        title = element_text(size = title_size),
        strip.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = 'top',
        legend.text=element_text(size=12)
  ) +
  
  geom_point(alpha = 0.8) +
  geom_smooth(method = lm, se = FALSE, col= 'red') +
  facet_grid(lineage~., scales = 'free_y')
#geom_point(aes(colour = lineage),alpha = 0.8,size = 2.5) +

#scale_color_brewer(palette = 'Set2')


pdf('CDR_FR_S5F_relative_difference.pdf', height=12, width=6)
plot(pl_ratio_vs_distance)
dev.off()

