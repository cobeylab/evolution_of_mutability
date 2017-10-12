# Produces figure for the manuscript showing the distribution of the fraction of negative changes in S5F, for all branches (Fig.2) and separately for different branch positions (Fig 5)
# Uses the following MCMC chains: CH103_con_run1a, CH103L_con_run1a, VRC26int_con_run1a, VRC26L_con_run1a, VRC01_01_log_run1a, VRC01_13_log_run1a, VRC01_19_log_run1a,

library('ggplot2')
library('grid')
library('coda')
#library('lattice')

# ================ LOADING GGPLOT2 PARAMETERS ===============
source('ggplot_parameters.R')

results_directory <- '../results/contrasts/observed_lineages/'

# ======= GET DATAFRAMES WITH FRACTION NEGATIVE LOSSES FROM ALL LINEAGES

CH103_dataframe <- read.table(paste(results_directory, 'CH103_constant/CH103_con_run1a_mutability_contrasts.csv', sep = ''), header = T, sep = ',')
CH103L_dataframe <- read.table(paste(results_directory, 'CH103L_constant/CH103L_con_run1a_mutability_contrasts.csv', sep = ''), header = T, sep = ',')
VRC26int_dataframe <- read.table(paste(results_directory, 'VRC26int_constant/VRC26int_con_run1a_mutability_contrasts.csv', sep = ''), header = T, sep = ',')
VRC26L_dataframe <- read.table(paste(results_directory, 'VRC26L_constant/VRC26L_con_run1a_mutability_contrasts.csv', sep = ''), header = T, sep = ',')
VRC01_01_dataframe <- read.table(paste(results_directory, 'VRC01_01_logistic/VRC01_01_log_run1a_mutability_contrasts.csv', sep = ''), header = T, sep = ',')
VRC01_13_dataframe <- read.table(paste(results_directory, 'VRC01_13_logistic/VRC01_13_log_run1a_mutability_contrasts.csv', sep = ''), header = T, sep = ',')
VRC01_19_dataframe <- read.table(paste(results_directory, 'VRC01_19_logistic/VRC01_19_log_run1a_mutability_contrasts.csv', sep = ''), header = T, sep = ',')

dataframe_list <- list('CH103' = CH103_dataframe, 'CH103L' = CH103L_dataframe,
                       'VRC26int' = VRC26int_dataframe, 'VRC26L' = VRC26L_dataframe,
                       'VRC01_01' = VRC01_01_dataframe, 'VRC01_13' = VRC01_13_dataframe,
                       'VRC01_19' = VRC01_19_dataframe)

# ==== MAKE DATAFRAME WITH FRACTION OF NEGATIVE CHANGES FOR EACH TREE, LINEAGE, REGION AND METRIC
fraction_negative_all <- c()
fraction_negative_terminal <- c()
fraction_negative_internal <- c()
fraction_negative_trunk <- c()
fraction_negative_nontrunk <- c()

lineage <- c()
region_column <- c()
metric_column <- c()
#branch_is_terminal <- c()
  
for(clone in c('CH103','CH103L','VRC26int','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
#for(clone in c('CH103','CH103L','VRC26L','VRC01_01','VRC01_13','VRC01_19')){  
  clone_dataframe <- dataframe_list[[clone]]
  print(clone)
  #For each tree from the posterior distribution for the clone:
  for(tree in unique(dataframe_list[[clone]]$tree)[1:1000]){
    
    # For each metric
    for(metric in c('S5F','HS','OHS','X7M','log_S5F')){
    
      # For each region
      for(region in c('WS','FR','CDR')){
      
        #print(paste(clone,metric,region,tree))
        # Find contrasts for the specified metric and region
        contrasts_all <- clone_dataframe[clone_dataframe$tree == tree,paste(metric,'_',region,'_contrast',sep='')]
        contrasts_terminal <- clone_dataframe[clone_dataframe$tree == tree & clone_dataframe$branch_is_terminal == T,
                                              paste(metric,'_',region,'_contrast',sep='')]
        
        # Contrast for all internal branches (both trunk and non-trunk)
        contrasts_internal <- clone_dataframe[clone_dataframe$tree == tree & clone_dataframe$branch_is_terminal == F,
                                              paste(metric,'_',region,'_contrast',sep='')]
        contrasts_trunk <- clone_dataframe[clone_dataframe$tree == tree & clone_dataframe$branch_in_trunk == T,
                                           paste(metric,'_',region,'_contrast',sep='')]
        
        # Contrasts for branches that are internal but not in the trunk
        contrasts_nontrunk <- clone_dataframe[clone_dataframe$tree == tree & clone_dataframe$branch_in_trunk == F &
                                                clone_dataframe$branch_is_terminal == F,
                                              paste(metric,'_',region,'_contrast',sep='')]
    
        f_negative_all <- sum(contrasts_all < 0) / sum(contrasts_all != 0)
        f_negative_terminal <- sum(contrasts_terminal < 0) / sum(contrasts_terminal != 0)
        f_negative_internal <- sum(contrasts_internal < 0) / sum(contrasts_internal != 0)
        f_negative_trunk <- sum(contrasts_trunk < 0) / sum(contrasts_trunk != 0)
        f_negative_nontrunk <- sum(contrasts_nontrunk < 0) / sum(contrasts_nontrunk != 0)
    
        fraction_negative_all <- c(fraction_negative_all, f_negative_all)
        fraction_negative_terminal <- c(fraction_negative_terminal, f_negative_terminal)
        fraction_negative_internal <- c(fraction_negative_internal, f_negative_internal)
        fraction_negative_trunk <- c(fraction_negative_trunk, f_negative_trunk)
        fraction_negative_nontrunk <- c(fraction_negative_nontrunk, f_negative_nontrunk)
        
        lineage <- c(lineage, clone)
        metric_column <- c(metric_column, metric)
        region_column <- c(region_column, region)
        }
    }
  }
}

lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26int','VRC26L','VRC01_13',
                                           'VRC01_01','VRC01_19'))
metric_column <- factor(metric_column, levels = c('S5F','HS','OHS','X7M','log_S5F'))
region_column <- factor(region_column, levels = c('WS','FR','CDR'))

  
global_dataframe <- data.frame(fraction_negative_all = fraction_negative_all, fraction_negative_terminal = fraction_negative_terminal,
                               fraction_negative_internal = fraction_negative_internal, fraction_negative_trunk = fraction_negative_trunk,
                               fraction_negative_nontrunk = fraction_negative_nontrunk, lineage = lineage,
                               region = region_column, metric = metric_column)

# ==== MAKE DATAFRAME WITH MEAN FRACTION NEG. AND HPD LIMITS FOR EACH LINEAGE, METRIC AND REGION

mean_fraction_negative_all <- c()
mean_fraction_negative_terminal <- c()
mean_fraction_negative_internal <- c()
mean_fraction_negative_trunk <- c()
mean_fraction_negative_nontrunk <- c()

lineage <- c()
metric_column <- c()
region_column <- c()
HPD_llim_all <- c()
HPD_ulim_all <- c()
HPD_llim_terminal <- c()
HPD_ulim_terminal <- c()
HPD_llim_internal <- c()
HPD_ulim_internal <- c()
HPD_ulim_trunk <- c()
HPD_llim_trunk <- c()
HPD_ulim_nontrunk <- c()
HPD_llim_nontrunk <- c()

for(clone in c('CH103','CH103L','VRC26int','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
#for(clone in c('CH103','CH103L','VRC26L','VRC01_01','VRC01_13','VRC01_19')){  
  for(metric in c('S5F','HS','OHS','X7M','log_S5F')){
    for(region in c('WS','FR','CDR')){
      values_all <- global_dataframe$fraction_negative_all[global_dataframe$lineage == clone &
                                                     global_dataframe$region == region &
                                                     global_dataframe$metric == metric]
      
      values_terminal <- global_dataframe$fraction_negative_terminal[global_dataframe$lineage == clone &
                                                             global_dataframe$region == region &
                                                             global_dataframe$metric == metric]
      
      values_internal <- global_dataframe$fraction_negative_internal[global_dataframe$lineage == clone &
                                                                       global_dataframe$region == region &
                                                                       global_dataframe$metric == metric]
      
      values_trunk <- global_dataframe$fraction_negative_trunk[global_dataframe$lineage == clone &
                                                                    global_dataframe$region == region &
                                                                    global_dataframe$metric == metric]
      
      values_nontrunk <- global_dataframe$fraction_negative_nontrunk[global_dataframe$lineage == clone &
                                                                 global_dataframe$region == region &
                                                                 global_dataframe$metric == metric]
      
      # Mean fraction
      mean_fraction_negative_all <- c(mean_fraction_negative_all, mean(values_all,na.rm=T))
      mean_fraction_negative_terminal <- c(mean_fraction_negative_terminal, mean(values_terminal,na.rm=T))
      mean_fraction_negative_internal <- c(mean_fraction_negative_internal, mean(values_internal,na.rm=T))
      mean_fraction_negative_trunk <- c(mean_fraction_negative_trunk, mean(values_trunk,na.rm=T))
      mean_fraction_negative_nontrunk <- c(mean_fraction_negative_nontrunk, mean(values_nontrunk,na.rm=T))
      
      
      # Temporary sub-dataframes for computing HPD interval
      sub_dat_all <- data.frame(values_all)
      fraction_HPD_limits_all <- HPDinterval(as.mcmc(sub_dat_all$values_all), 0.95)
      HPD_llim_all <- c(HPD_llim_all, fraction_HPD_limits_all[1])
      HPD_ulim_all <- c(HPD_ulim_all, fraction_HPD_limits_all[2])
      
      sub_dat_terminal <- data.frame(values_terminal)
      fraction_HPD_limits_terminal <- HPDinterval(as.mcmc(sub_dat_terminal$values_terminal), 0.95)
      HPD_llim_terminal <- c(HPD_llim_terminal, fraction_HPD_limits_terminal[1])
      HPD_ulim_terminal <- c(HPD_ulim_terminal, fraction_HPD_limits_terminal[2])
      
      sub_dat_internal <- data.frame(values_internal)
      fraction_HPD_limits_internal <- HPDinterval(as.mcmc(sub_dat_internal$values_internal), 0.95)
      HPD_llim_internal <- c(HPD_llim_internal, fraction_HPD_limits_internal[1])
      HPD_ulim_internal <- c(HPD_ulim_internal, fraction_HPD_limits_internal[2])
      
      sub_dat_trunk <- data.frame(values_trunk)
      fraction_HPD_limits_trunk <- HPDinterval(as.mcmc(sub_dat_trunk$values_trunk), 0.95)
      HPD_llim_trunk <- c(HPD_llim_trunk, fraction_HPD_limits_trunk[1])
      HPD_ulim_trunk <- c(HPD_ulim_trunk, fraction_HPD_limits_trunk[2])
      
      sub_dat_nontrunk <- data.frame(values_nontrunk)
      fraction_HPD_limits_nontrunk <- HPDinterval(as.mcmc(sub_dat_nontrunk$values_nontrunk), 0.95)
      HPD_llim_nontrunk <- c(HPD_llim_nontrunk, fraction_HPD_limits_nontrunk[1])
      HPD_ulim_nontrunk <- c(HPD_ulim_nontrunk, fraction_HPD_limits_nontrunk[2])
      
      metric_column <- c(metric_column, metric)
      region_column <- c(region_column, region)
      lineage <- c(lineage, clone)
    }
  }
}

lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26int','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
metric_column <- factor(metric_column, levels = c('S5F','HS','OHS','X7M','log_S5F'))
region_column <- factor(region_column, levels = c('WS','FR','CDR'))  

  
summary_dataframe <- data.frame(lineage,metric=metric_column,region=region_column,
                                mean_fraction_negative_all,HPD_llim_all,HPD_ulim_all,
                                mean_fraction_negative_terminal,HPD_llim_terminal,HPD_ulim_terminal,
                                mean_fraction_negative_internal,HPD_llim_internal,HPD_ulim_internal,
                                mean_fraction_negative_trunk,HPD_llim_trunk,HPD_ulim_trunk,
                                mean_fraction_negative_nontrunk,HPD_llim_nontrunk,HPD_ulim_nontrunk)


# ================================= PLOTS =====================================
metric_names <- c(
  'S5F' = 'S5F',
  'HS' = ' WRCH/DGYW hotspots',
  'OHS' = 'Overlapping hotspots',
  'X7M' = '7-mer mutability',
  'log_S5F' = 'log-S5F mutability'
)
region_names <- c(
  'WS'="Whole sequence",
  'FR'="FRs",
  'CDR'="CDRs"
)

#================= Plot with S5F contrasts for all branches combined (main text):  

base_plot <- ggplot(data = summary_dataframe, aes(x=lineage,y=mean_fraction_negative_all)) + 
    scale_y_continuous(expand = c(0,0), limits = c(-0.05,1.05)) + 
    theme_bw() +
    #theme_classic() + 
    ylab('Frequency of losses') + 
    xlab('Lineage') +
    ggplot_theme +
    
    # Add horizontal line at 0.5
    geom_hline(yintercept=0.5,linetype=5, size = 0.5)

    # Change x-axis labels
    x_labels <- c()
    for(level in levels(global_dataframe$lineage)){
      x_labels <- c(x_labels, 
                    switch(level,
                           CH103 = 'CH103 (H)',
                           CH103L = 'CH103 (L)',
                           VRC26int = 'VRC26 (H)',
                           VRC26 = 'VRC26 (H)',
                           VRC26L = 'VRC26 (L)',
                           VRC01_13 = 'VRC01-13(H)',
                           VRC01_01 = 'VRC01-01(H)',
                           VRC01_19 = 'VRC01-19(H)'
                    )
      )
    }
    base_plot <- base_plot + scale_x_discrete(labels=x_labels) + theme(axis.text.x = element_text(angle=90, size = 8))

pl_S5F <- base_plot + 
  
  geom_linerange(data=subset(summary_dataframe, metric == 'S5F'),
                 aes(ymin=HPD_llim_all, 
                     ymax = HPD_ulim_all), colour = 'red4') +
    
  geom_point(data=subset(summary_dataframe, metric == 'S5F'),aes(y=mean_fraction_negative_all), 
             shape = 21, fill = 'firebrick1', alpha = 0.7,size = 2.5)  +
  
  facet_grid(.~region, labeller = labeller(region=region_names)) + 
  theme(strip.text.x = element_text(size = 12))
  
pdf('contrasts_observed.pdf', width = 7.01, height = 3.43)
plot(pl_S5F)
dev.off()

#======= Plot with mean-log-S5F contrasts
pl_logS5F <- base_plot + 
  
  geom_linerange(data=subset(summary_dataframe, metric == 'log_S5F'),
                 aes(ymin=HPD_llim_all, 
                     ymax = HPD_ulim_all), colour = 'red4') +
  
  geom_point(data=subset(summary_dataframe, metric == 'log_S5F'),aes(y=mean_fraction_negative_all), 
             shape = 21, fill = 'firebrick1', alpha = 0.7,size = 2.5)  +
  
  facet_grid(.~region, labeller = labeller(region=region_names)) + 
  theme(strip.text.x = element_text(size = 12))

pdf('geom_and_log_S5F_plots/contrasts_observed_logS5F.pdf', width = 7.01, height = 3.43)
plot(pl_logS5F)
dev.off()


#======= Plot with results for other metrics (supplement)


 pl_other_metrics <- ggplot(data=subset(summary_dataframe, metric %in% c('HS','OHS','X7M')), 
                            aes(x=lineage,y=mean_fraction_negative_all)) + 
   scale_y_continuous(expand = c(0,0), limits = c(-0.05,1.05)) + 
   theme_bw() +
   #theme_classic() + 
   ylab('Frequency of losses') + 
   xlab('Lineage') +
   ggplot_theme +
   
   # Add horizontal line at 0.5
   geom_hline(yintercept=0.5,linetype=5, size = 0.5)
 
   # Change x-axis labels
   x_labels <- c()
   for(level in levels(global_dataframe$lineage)){
     x_labels <- c(x_labels, 
                   switch(level,
                          CH103 = 'CH103 (H)',
                          CH103L = 'CH103 (L)',
                          VRC26int = 'VRC26 (H)',
                          VRC26 = 'VRC26 (H)',
                          VRC26L = 'VRC26 (L)',
                          VRC01_13 = 'VRC01-13(H)',
                          VRC01_01 = 'VRC01-01(H)',
                          VRC01_19 = 'VRC01-19(H)'
                   )
     )
   }
   
pl_other_metrics <- pl_other_metrics + scale_x_discrete(labels=x_labels) + 
  theme(axis.text.x = element_text(angle=90, size = 8)) +
  geom_point(data=subset(summary_dataframe, metric %in% c('HS','OHS','X7M')),
             aes(y=mean_fraction_negative_all),
             shape = 21, fill = 'firebrick1', alpha = 0.7,size = 3)  +

  geom_linerange(data=subset(summary_dataframe, metric %in% c('HS','OHS','X7M')),
                 aes(ymin=HPD_llim_all,
                     ymax = HPD_ulim_all), colour = 'red4') +
  facet_grid(metric~region, labeller = labeller(metric = metric_names, region = region_names)) +
  theme(strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 9))

pdf('contrasts_observed_supp.pdf', width = 7.01, height = 7.01)
plot(pl_other_metrics)
dev.off()

