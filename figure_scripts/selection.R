# Produces figure for the manuscript showing the results of the selection analysis using Baseline
# Uses the following MCMC chains: CH103_con_run1a, CH103L_con_run1a, VRC26int_con_run1a, VRC26L_con_run1a, VRC01_01_log_run1a, VRC01_13_log_run1a, VRC01_19_log_run1a,

library('ggplot2')
library('cowplot')
library('grid')
library('lattice')

results_directory <- '../results/selection/'

source('ggplot_parameters.R')

# ======= GET DATAFRAMES WITH SELECTION RESULTS FOR ALL LINEAGES

CH103_dataframe <- read.table(paste(results_directory, 'CH103_con_run1a_baseline.tsv', sep = ''), header = T, sep = '\t')
CH103L_dataframe <- read.table(paste(results_directory, 'CH103L_con_run1a_baseline.tsv', sep = ''), header = T, sep = '\t')
VRC26_dataframe <- read.table(paste(results_directory, 'VRC26int_con_run1a_baseline.tsv', sep = ''), header = T, sep = '\t')
VRC26L_dataframe <- read.table(paste(results_directory, 'VRC26L_con_run1a_baseline.tsv', sep = ''), header = T, sep = '\t')
VRC01_01_dataframe <- read.table(paste(results_directory, 'VRC01_01_log_run1a_baseline.tsv', sep = ''), header = T, sep = '\t')
VRC01_13_dataframe <- read.table(paste(results_directory, 'VRC01_13_log_run1a_baseline.tsv', sep = ''), header = T, sep = '\t')
VRC01_19_dataframe <- read.table(paste(results_directory, 'VRC01_19_log_run1a_baseline.tsv', sep = ''), header = T, sep = '\t')

dataframe_list <- list('CH103' = CH103_dataframe, 'CH103L' = CH103L_dataframe,
                           'VRC26' = VRC26_dataframe, 'VRC26L' = VRC26L_dataframe,
                           'VRC01_01' = VRC01_01_dataframe, 'VRC01_13' = VRC01_13_dataframe,
                           'VRC01_19' = VRC01_19_dataframe)    

# ======= MAKE COMBINED DATAFRAME WITH FR AND CDR RESULTS FOR EACH LINEAGE =========
lineage <- c()
region_vector <- c()
sigma <- c()
sigma_CI_lower <- c()
sigma_CI_upper <- c()

for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  clone_dataframe <- dataframe_list[[clone]] 
  
  for(region in c('FWR','CDR')){
    lineage <- c(lineage, clone)
    region_vector <- c(region_vector, region)
    sigma <- c(sigma, clone_dataframe[1, paste('Focused_Sigma_',region,sep='')])
    sigma_CI_lower <- c(sigma_CI_lower, clone_dataframe[1, paste('Focused_CIlower_',region,sep='')])
    sigma_CI_upper <- c(sigma_CI_upper, clone_dataframe[1, paste('Focused_CIupper_',region,sep='')])
    
  }
}
  
lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
region_vector <- factor(region_vector, levels = c('FWR','CDR'))

combined_dataframe <- data.frame(lineage, region = region_vector, sigma, sigma_CI_lower, sigma_CI_upper)


# =====================================================================================================================
# ================================================== PLOTS ============================================================
# =====================================================================================================================

pl <- ggplot(combined_dataframe, 
                 aes(x = lineage, y = sigma)) +
  geom_hline(yintercept=0,linetype = 2) +
  theme_bw() +
  #theme_classic() + 
  ylab('Selection strength') + 
  xlab('Lineage') +
  ggplot_theme + 
  theme(legend.position = 'top',
        legend.text=element_text(size=12)) +
  
  theme(axis.text.x = element_text(angle=90, size = 8)) +
  
  geom_linerange(aes(ymin = sigma_CI_lower, ymax = sigma_CI_upper,
                     colour = region)) +
  
  geom_point(aes(fill = region), size = 2.5, shape = 21,alpha = 0.7) + 
  
  scale_x_discrete(labels = c('CH103' = 'CH103 (H)', 'CH103L' = 'CH103 (L)',
                              'VRC26' = 'VRC26 (H)', 'VRC26L' = 'VRC26 (L)',
                              'VRC01_13' = 'VRC01-13 (H)','VRC01_01' = 'VRC01-01 (H)',
                              'VRC01_19' = 'VRC01-19 (H)'
  )) +
  scale_fill_manual(values = c('#f1a340', '#998ec3'),
                    labels = c('FRs','CDRs')) +
  
  scale_colour_manual(values = c('#f1a340', '#998ec3'),
                    labels = c('FRs','CDRs')) +
  
  
  guides(fill=guide_legend(title=NULL), colour = guide_legend(title=NULL))
  

pdf('selection.pdf', width = 3.43, height = 4)
plot(pl)
dev.off()