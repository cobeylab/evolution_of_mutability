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

results_directory <- '../../results/simulations_MCC/observed_lineages/'

# ======= GET DATAFRAMES WITH OBSERVED MCC MUTABILITY RESULTS FOR ALL LINEAGES

CH103_dataframe_obs <- read.table(paste(results_directory, 'CH103_constant/CH103_con_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
CH103L_dataframe_obs <- read.table(paste(results_directory, 'CH103L_constant/CH103L_con_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
VRC26_dataframe_obs <- read.table(paste(results_directory, 'VRC26int_constant/VRC26int_con_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
VRC26L_dataframe_obs <- read.table(paste(results_directory, 'VRC26L_constant/VRC26L_con_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_01_dataframe_obs <- read.table(paste(results_directory, 'VRC01_01_logistic/VRC01_01_log_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_13_dataframe_obs <- read.table(paste(results_directory, 'VRC01_13_logistic/VRC01_13_log_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_19_dataframe_obs <- read.table(paste(results_directory, 'VRC01_19_logistic/VRC01_19_log_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')

dataframe_list_obs <- list('CH103' = CH103_dataframe_obs, 'CH103L' = CH103L_dataframe_obs,
                           'VRC26' = VRC26_dataframe_obs, 'VRC26L' = VRC26L_dataframe_obs,
                           'VRC01_01' = VRC01_01_dataframe_obs, 'VRC01_13' = VRC01_13_dataframe_obs,
                           'VRC01_19' = VRC01_19_dataframe_obs)                        


# ====== COMBINED DATAFRAME WITH RESULTS FOR ALL LINEAGES =======
lineage_vector <- c()
S5F_CDR <- c()
S5F_FR <- c()
time_from_root <- c()
distance_from_root <- c()
for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  clone_dataframe_obs <- dataframe_list_obs[[clone]] 
  
  # Get points for internal nodes:
  for(internal_node in unique(clone_dataframe_obs$parent)){
    node_dataframe <- clone_dataframe_obs[clone_dataframe_obs$parent == internal_node,]
    
    S5F_FR <- c(S5F_FR,node_dataframe$S5F_FR_parent[1])
    S5F_CDR <- c(S5F_CDR,node_dataframe$S5F_CDR_parent[1])
    
    if(grepl('VRC01',clone)){
      time_from_root <- c(time_from_root, 4*node_dataframe$parent_time_to_root[1])
    }else{
      time_from_root <- c(time_from_root, node_dataframe$parent_time_to_root[1])
    }
    
    distance_from_root <- c(distance_from_root, node_dataframe$parent_distance_to_root[1])
    lineage_vector <- c(lineage_vector, clone)
  }
  # Get points for terminal nodes:
  for(terminal_node in unique(clone_dataframe_obs$child)){
    node_dataframe <- clone_dataframe_obs[clone_dataframe_obs$child == terminal_node,]
    
    S5F_FR <- c(S5F_FR,node_dataframe$S5F_FR_child)
    S5F_CDR <- c(S5F_CDR,node_dataframe$S5F_CDR_child)
    
    if(grepl('VRC01',clone)){
      time_from_root <- c(time_from_root, 4*node_dataframe$parent_time_to_root[1])
    }else{
      time_from_root <- c(time_from_root, node_dataframe$parent_time_to_root[1])
    }
    distance_from_root <- c(distance_from_root, node_dataframe$parent_distance_to_root + node_dataframe$branch_exp_subs)
    lineage_vector <- c(lineage_vector, clone)
  }
}

lineage_vector <- factor(lineage_vector, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))

combined_dataframe <- data.frame(lineage=lineage_vector,plot_row,time_from_root,distance_from_root,
                                 S5F_FR, S5F_CDR, CDR_FR_S5F_ratio = S5F_CDR/S5F_FR, CDR_FR_S5F_difference = S5F_CDR - S5F_FR)

# ================ DEFINING GGPLOT2 PARAMETERS ===============

{
  title_size <- 20
  axis_title_size <- 18
  y_axis_text_size <- 12
  x_axis_text_size <- 12
  
  axis_title_size_subplot <- 20
  y_axis_text_size_subplot <- 15
  x_axis_text_size_subplot <- 16
  
  ylab_distance <- 16
  xlab_distance <- 14
  plot_margins <- c(0.3, 1, 0.3, 2)

  lineage_names <- c(
    'CH103' = "CH103 (H)",
    'CH103L' = "CH103 (L)",
    'VRC26' = "VRC26 (H)",
    'VRC26L' = "VRC26 (L)",
    'VRC01_13' = 'VRC01-13 (H)',
    'VRC01_01' = 'VRC01-01 (H)',
    'VRC01_19' = 'VRC01-19 (H)'
  )
  
}

# =====================================================================================================================
# ================================================== PLOTS ============================================================
# =====================================================================================================================

# ======== PLOT OF CDR/FR S5F difference VS distance.
pl_ratio_vs_difference <- ggplot(combined_dataframe, 
                               aes(x = distance_from_root, y = CDR_FR_S5F_difference/S5F_FR)) +
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
plot(pl_ratio_vs_difference)
dev.off()

