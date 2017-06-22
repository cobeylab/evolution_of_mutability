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

source('ggplot_parameters.R')

results_directory <- '../results/S_NS_mutability_changes/observed_lineages/'

# ======= GET DATAFRAMES WITH OBSERVED MCC MUTABILITY RESULTS FOR ALL LINEAGES

CH103_dataframe_obs <- read.table(paste(results_directory, 'CH103_constant/CH103_con_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
CH103L_dataframe_obs <- read.table(paste(results_directory, 'CH103L_constant/CH103L_con_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
VRC26_dataframe_obs <- read.table(paste(results_directory, 'VRC26int_constant/VRC26int_con_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
VRC26L_dataframe_obs <- read.table(paste(results_directory, 'VRC26L_constant/VRC26L_con_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_01_dataframe_obs <- read.table(paste(results_directory, 'VRC01_01_logistic/VRC01_01_log_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_13_dataframe_obs <- read.table(paste(results_directory, 'VRC01_13_logistic/VRC01_13_log_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_19_dataframe_obs <- read.table(paste(results_directory, 'VRC01_19_logistic/VRC01_19_log_run1a_observed_MCC.csv', sep = ''), header = T, sep = ',')


# ======= GET DATAFRAMES WITH SIMULATED MCC MUTABILITY RESULTS FOR ALL LINEAGES
CH103_dataframe_sim <- read.table(paste(results_directory, 'CH103_constant/CH103_con_run1a_simulated_MCC.csv', sep = ''), header = T, sep = ',')
CH103L_dataframe_sim <- read.table(paste(results_directory, 'CH103L_constant/CH103L_con_run1a_simulated_MCC.csv', sep = ''), header = T, sep = ',')
VRC26_dataframe_sim <- read.table(paste(results_directory, 'VRC26int_constant/VRC26int_con_run1a_simulated_MCC.csv', sep = ''), header = T, sep = ',')
VRC26L_dataframe_sim <- read.table(paste(results_directory, 'VRC26L_constant/VRC26L_con_run1a_simulated_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_01_dataframe_sim <- read.table(paste(results_directory, 'VRC01_01_logistic/VRC01_01_log_run1a_simulated_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_13_dataframe_sim <- read.table(paste(results_directory, 'VRC01_13_logistic/VRC01_13_log_run1a_simulated_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_19_dataframe_sim <- read.table(paste(results_directory, 'VRC01_19_logistic/VRC01_19_log_run1a_simulated_MCC.csv', sep = ''), header = T, sep = ',')


dataframe_list_obs <- list('CH103' = CH103_dataframe_obs, 'CH103L' = CH103L_dataframe_obs,
                           'VRC26' = VRC26_dataframe_obs, 'VRC26L' = VRC26L_dataframe_obs,
                           'VRC01_01' = VRC01_01_dataframe_obs, 'VRC01_13' = VRC01_13_dataframe_obs,
                           'VRC01_19' = VRC01_19_dataframe_obs)                        

dataframe_list_sim <- list('CH103' = CH103_dataframe_sim, 'CH103L' = CH103L_dataframe_sim,
                        'VRC26' = VRC26_dataframe_sim, 'VRC26L' = VRC26L_dataframe_sim,
                        'VRC01_01' = VRC01_01_dataframe_sim, 'VRC01_13' = VRC01_13_dataframe_sim,
                        'VRC01_19' = VRC01_19_dataframe_sim)


# ====== COMBINED DATAFRAME WITH RESULTS FOR ALL LINEAGES =======
lineage_vector <- c()
metric_vector <- c()
substitution_class <- c()

mean_contrast_true <- c()

mean_contrast_simulation_S5F_S5F <- c()
contrast_simulation_S5F_S5F_llim <- c()
contrast_simulation_S5F_S5F_ulim <- c()
mean_correlation_S5F_S5F <- c()
correlation_S5F_S5F_llim <- c()
correlation_S5F_S5F_ulim <- c()

mean_contrast_simulation_uniform_S5F <- c()
contrast_simulation_uniform_S5F_llim <- c()
contrast_simulation_uniform_S5F_ulim <- c()
mean_correlation_uniform_S5F <- c()
correlation_uniform_S5F_llim <- c()
correlation_uniform_S5F_ulim <- c()

mean_contrast_simulation_CP_S5F <- c()
contrast_simulation_CP_S5F_llim <- c()
contrast_simulation_CP_S5F_ulim <- c()
mean_correlation_CP_S5F <- c()
correlation_CP_S5F_llim <- c()
correlation_CP_S5F_ulim <- c()

region <- 'WS'
for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  clone_dataframe_obs <- dataframe_list_obs[[clone]] 
  clone_dataframe_sim <- dataframe_list_sim[[clone]]
  
  #for(metric in c('S5F','HS','OHS')){
  for(metric in c('S5F')){
    for(sub_type in c('total','syn_only','nonsyn_only')){
      
      parent_mutability <- clone_dataframe_obs[,paste(metric, '_parent', sep ='')]
      
      child_mutability <- clone_dataframe_obs[,paste(metric, '_child', sep ='')]

      contrast_syn <- clone_dataframe_obs$S5F_change_syn
      contrast_nonsyn <- clone_dataframe_obs$S5F_change_nonsyn
      
      if(sub_type == 'total'){
        contrast_true <- contrast_syn + contrast_nonsyn
      }
      if(sub_type == 'syn_only'){
        contrast_true <- contrast_syn
      }
      if(sub_type == 'nonsyn_only'){
        contrast_true <- contrast_nonsyn
      }
      
      mean_contrast_true <- c(mean_contrast_true, mean(contrast_true))
      
      simulated_means_S5F_S5F <- c()
      simulated_means_uniform_S5F <- c()
      simulated_means_CP_S5F <- c()
      
      spearman_coefficient_S5F_S5F <- c()
      spearman_coefficient_uniform_S5F <- c()
      spearman_coefficient_CP_S5F <- c()
      
      for(replicate in 1:max(clone_dataframe_sim$replicate)){
        
        if(sub_type == 'total'){
          contrasts_S5F_S5F <- clone_dataframe_sim[clone_dataframe_sim$replicate == replicate,
                                                   paste(metric,'_',region, '_change_S5FMut_S5FTrans_total', sep ='')]
          contrasts_uniform_S5F <- clone_dataframe_sim[clone_dataframe_sim$replicate == replicate,
                                                   paste(metric,'_',region, '_change_uniformMut_S5FTrans_total', sep ='')]
          contrasts_CP_S5F <- clone_dataframe_sim[clone_dataframe_sim$replicate == replicate,
                                                       paste(metric,'_',region, '_change_CPMut_S5FTrans_total', sep ='')]
        }
        if(sub_type == 'syn_only'){
          contrasts_S5F_S5F <- clone_dataframe_sim[clone_dataframe_sim$replicate == replicate,
                                                   paste(metric,'_',region, '_change_S5FMut_S5FTrans_syn', sep ='')]
          contrasts_uniform_S5F <- clone_dataframe_sim[clone_dataframe_sim$replicate == replicate,
                                                       paste(metric,'_',region, '_change_uniformMut_S5FTrans_syn', sep ='')]
          contrasts_CP_S5F <- clone_dataframe_sim[clone_dataframe_sim$replicate == replicate,
                                                  paste(metric,'_',region, '_change_CPMut_S5FTrans_syn', sep ='')]
        }
        
        if(sub_type == 'nonsyn_only'){
          contrasts_S5F_S5F <- clone_dataframe_sim[clone_dataframe_sim$replicate == replicate,
                                                   paste(metric,'_',region, '_change_S5FMut_S5FTrans_nonsyn', sep ='')]
          contrasts_uniform_S5F <- clone_dataframe_sim[clone_dataframe_sim$replicate == replicate,
                                                       paste(metric,'_',region, '_change_uniformMut_S5FTrans_nonsyn', sep ='')]
          contrasts_CP_S5F <- clone_dataframe_sim[clone_dataframe_sim$replicate == replicate,
                                                  paste(metric,'_',region, '_change_CPMut_S5FTrans_nonsyn', sep ='')]
        }
        
        simulated_means_S5F_S5F <- c(simulated_means_S5F_S5F, mean(contrasts_S5F_S5F))
        simulated_means_uniform_S5F <- c(simulated_means_uniform_S5F, mean(contrasts_uniform_S5F))
        simulated_means_CP_S5F <- c(simulated_means_CP_S5F, mean(contrasts_CP_S5F))
        
        spearman_corr_S5F_S5F <- cor.test(contrasts_S5F_S5F,contrast_true,method = 'spearman')$estimate
        spearman_corr_uniform_S5F <- cor.test(contrasts_uniform_S5F,contrast_true,method = 'spearman')$estimate
        spearman_corr_CP_S5F <- cor.test(contrasts_CP_S5F,contrast_true,method = 'spearman')$estimate
        
        spearman_coefficient_S5F_S5F <- c(spearman_coefficient_S5F_S5F,  spearman_corr_S5F_S5F)
        spearman_coefficient_uniform_S5F <- c(spearman_coefficient_uniform_S5F,  spearman_corr_uniform_S5F)
        spearman_coefficient_CP_S5F <- c(spearman_coefficient_CP_S5F,  spearman_corr_CP_S5F)
        
      }
 
      ##### COMPUTE AVERAGE OF MEAN CONTRASTS AND 95% INTERVAL
      mean_contrast_simulation_S5F_S5F <- c(mean_contrast_simulation_S5F_S5F, mean(simulated_means_S5F_S5F))
      mean_contrast_simulation_uniform_S5F <- c(mean_contrast_simulation_uniform_S5F, mean(simulated_means_uniform_S5F))
      mean_contrast_simulation_CP_S5F <- c(mean_contrast_simulation_CP_S5F, mean(simulated_means_CP_S5F))
      
      # 95% INTERVAL
      contrast_simulation_S5F_S5F_llim <- c(contrast_simulation_S5F_S5F_llim, quantile(simulated_means_S5F_S5F,0.025))
      contrast_simulation_S5F_S5F_ulim <- c(contrast_simulation_S5F_S5F_ulim,quantile(simulated_means_S5F_S5F,0.975))
      
      contrast_simulation_uniform_S5F_llim <- c(contrast_simulation_uniform_S5F_llim, quantile(simulated_means_uniform_S5F,0.025))
      contrast_simulation_uniform_S5F_ulim <- c(contrast_simulation_uniform_S5F_ulim,quantile(simulated_means_uniform_S5F,0.975))
      
      contrast_simulation_CP_S5F_llim <- c(contrast_simulation_CP_S5F_llim, quantile(simulated_means_CP_S5F,0.025))
      contrast_simulation_CP_S5F_ulim <- c(contrast_simulation_CP_S5F_ulim,quantile(simulated_means_CP_S5F,0.975))
      
      substitution_class <- c(substitution_class, sub_type)
      metric_vector <- c(metric_vector, metric)
      lineage_vector <- c(lineage_vector, clone)
      
      ##### COMPUTE MEAN CORRELATION WITH TRUE VALUES
      mean_correlation_S5F_S5F <- c(mean_correlation_S5F_S5F, mean(spearman_coefficient_S5F_S5F))
      mean_correlation_uniform_S5F <- c(mean_correlation_uniform_S5F, mean(spearman_coefficient_uniform_S5F))
      mean_correlation_CP_S5F <- c(mean_correlation_CP_S5F, mean(spearman_coefficient_CP_S5F))
      
      correlation_S5F_S5F_llim <- c(correlation_S5F_S5F_llim, quantile(spearman_coefficient_S5F_S5F,0.025))
      correlation_uniform_S5F_llim <- c(correlation_uniform_S5F_llim, quantile(spearman_coefficient_uniform_S5F,0.025))
      correlation_CP_S5F_llim <- c(correlation_CP_S5F_llim, quantile(spearman_coefficient_CP_S5F,0.025))
      
      correlation_S5F_S5F_ulim <- c(correlation_S5F_S5F_ulim, quantile(spearman_coefficient_S5F_S5F,0.975))
      correlation_uniform_S5F_ulim <- c(correlation_uniform_S5F_ulim, quantile(spearman_coefficient_uniform_S5F,0.975))
      correlation_CP_S5F_ulim <- c(correlation_CP_S5F_ulim, quantile(spearman_coefficient_CP_S5F,0.975))
      
    }
      
  }
}

lineage_vector <- factor(lineage_vector, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
metric_column <- factor(metric_column, levels = c('S5F','HS','OHS'))
region_column <- factor(region_column, levels = c('WS','FR','CDR'))

combined_dataframe <- data.frame(lineage=lineage_vector, metric = metric_vector, substitution_class, mean_contrast_true, 
                                 mean_contrast_simulation_S5F_S5F, contrast_simulation_S5F_S5F_llim, contrast_simulation_S5F_S5F_ulim,
                                 mean_contrast_simulation_uniform_S5F, contrast_simulation_uniform_S5F_llim, contrast_simulation_uniform_S5F_ulim,
                                 mean_contrast_simulation_CP_S5F, contrast_simulation_CP_S5F_llim, contrast_simulation_CP_S5F_ulim,
                                 mean_correlation_S5F_S5F, correlation_S5F_S5F_llim, correlation_S5F_S5F_ulim,
                                 mean_correlation_uniform_S5F, correlation_uniform_S5F_llim, correlation_uniform_S5F_ulim,
                                 mean_correlation_CP_S5F, correlation_CP_S5F_llim, correlation_CP_S5F_ulim)
  

# DATAFRAME WITH SUMMED CHANGES IN MUTABILITY ACROSS THE TREE
lineage_vector <- c()
region_vector <- c()
cumulative_change <- c()
substitution_class <- c()
group <- c()

for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){

  clone_dataframe_obs <- dataframe_list_obs[[clone]] 
  
  for(region in c('WS','FR','CDR')){
    region_id <- ifelse(region == 'WS','', paste('_',region,sep=''))
    for(class in c('syn_only','nonsyn_only')){
  
      if(class == 'syn_only'){
        contrasts <- clone_dataframe_obs[,paste('S5F_change_syn',region_id, sep='')]
      }else{
        contrasts <- clone_dataframe_obs[,paste('S5F_change_nonsyn',region_id, sep='')]
      }
      
      lineage_vector <- c(lineage_vector, clone)
      region_vector <- c(region_vector, region)
      substitution_class <- c(substitution_class, class)
      cumulative_change <- c(cumulative_change, sum(contrasts))
      
      group <- c(group, paste(region, class, sep ='_'))
    }
  }
} 
lineage_vector <- factor(lineage_vector, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                  'VRC01_01','VRC01_19'))
group <- factor(group, levels = c('WS_syn_only','WS_nonsyn_only','FR_syn_only','FR_nonsyn_only',
                                  'CDR_syn_only','CDR_nonsyn_only'))
region_vector <- factor(region_vector, levels = c('WS','FR','CDR'))
substitution_class <- factor(substitution_class, levels = c('syn_only','nonsyn_only'))
total_change_dataframe <- data.frame(lineage=lineage_vector, region = region_vector, 
                                     substitution_class,cumulative_change, group)

# DATAFRAME TO PLOT CONTRASTS VS. TIME FOR ALL LINEAGES TOGETHER:
lineage <- c()
change_total <- c()
change_syn <- c()
change_nonsyn <- c()
parent_time_to_root <- c()
parent_distance_to_root <- c()

for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  clone_dataframe <- dataframe_list_obs[[clone]]

  change_syn <- c(change_syn, clone_dataframe$S5F_WS_change_syn)
  change_nonsyn <- c(change_nonsyn, clone_dataframe$S5F_WS_change_nonsyn)
  change_total <- change_syn + change_nonsyn
    
  lineage <- c(lineage, rep(clone, nrow(clone_dataframe)))
  
  time <- clone_dataframe$parent_time_to_root
  if(grepl('VRC01',clone)){
    time <- 4 * time
  }
  parent_time_to_root <- c(parent_time_to_root, time)
  parent_distance_to_root <- c(parent_distance_to_root, clone_dataframe$parent_distance_to_root)
  
}

lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                                    'VRC01_01','VRC01_19'))
contrasts_vs_time_dataframe <- data.frame(lineage, parent_time_to_root, parent_distance_to_root,
                                         change_total, change_syn, change_nonsyn,
                                         cumulative_change_up_to_time, cumulative_change_up_to_distance)

# ===== DATAFRAMES WITH CUMULATIVE MUTABILITY CHANGE VS TIME AND VS DISTANCE :
lineage_time <- c()
lineage_distance <- c()
time_from_root <- c()
distance_to_root <- c()
cumulative_change_up_to_time <- c()
cumulative_change_up_to_distance <- c()

for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  clone_dataframe <- dataframe_list_obs[[clone]]
  
  change_syn <- clone_dataframe$S5F_WS_change_syn
  change_nonsyn <- clone_dataframe$S5F_WS_change_nonsyn
  change_total <- change_syn + change_nonsyn
  
  distance <- clone_dataframe$parent_distance_to_root
  time <- clone_dataframe$parent_time_to_root
  if(grepl('VRC01',clone)){
    time <- 4 * time
  }
  
  time_from_root <- c(time_from_root, sort(unique(time)))
  distance_to_root <- c(distance_to_root,sort(unique(distance)))
  
  lineage_time <- c(lineage_time, rep(clone, length(unique(time))))
  lineage_distance <- c(lineage_distance, rep(clone, length(unique(distance))))
  
  for(t in sort(unique(time))){
    cumulative_change_up_to_time <- c(cumulative_change_up_to_time,
                                      sum(change_total[time < t]))
  }
  for(d in sort(unique(distance))){
    cumulative_change_up_to_distance <- c(cumulative_change_up_to_distance,
                                      sum(change_total[distance < d]))
  }
  
}

lineage_time <- factor(lineage_time, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
lineage_distance <- factor(lineage_distance, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                                'VRC01_01','VRC01_19'))


cumulative_change_time_dataframe <- data.frame(lineage_time, time_from_root,
                                                  cumulative_change_up_to_time)
cumulative_change_distance_dataframe <- data.frame(lineage_distance, distance_to_root,
                                               cumulative_change_up_to_distance)


# COMPUTING AVERAGE FRACTION EXPLAINED BY NON-SYN.
syn <- c()
nonsyn <- c()
syn_fraction <- c()
nonsyn_fraction <- c()
total <- c()

for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_13','VRC01_01','VRC01_19')){
  syn_change <- total_change_dataframe[total_change_dataframe$lineage == clone & substitution_class == 'syn_only',
                                       'cumulative_change']
  
  nonsyn_change <- total_change_dataframe[total_change_dataframe$lineage == clone & substitution_class == 'nonsyn_only',
                                          'cumulative_change']
  
  syn <- c(syn, syn_change)
  nonsyn <- c(nonsyn, nonsyn_change)
  
  total_change <- syn_change + nonsyn_change
  
  total <- c(total, total_change)
  
  syn_fraction <- c(syn_fraction, syn_change / total_change)
  nonsyn_fraction <- c(nonsyn_fraction, nonsyn_change / total_change)
  
}



# =====================================================================================================================
# ================================================== PLOTS ============================================================
# =====================================================================================================================

# ====================================== CUMULATIVE CHANGES IN MUTABILITY =============================================
total_changes_pl <- ggplot(subset(total_change_dataframe, region != 'WS'),
                           aes(y=cumulative_change,x=region)) +

  geom_col(aes(fill = substitution_class), 
           width = 0.7, position=position_stack()) + 
  
  facet_grid(~lineage) +
  geom_hline(yintercept=0,linetype = 2) +
  theme_bw() +
  #theme_classic() + 
  ylab('Cumulative change in average mutability') + 
  xlab('Region') +
  theme(legend.position = 'top',
        legend.text=element_text(size=11)
  ) +
  #scale_y_continuous(limits = c(-10,2)) + 
  scale_fill_manual(values = c('gray80','gray40'),
                    labels = c('Synonymous  ','Non-synonynous' )) +
  
  scale_x_discrete(labels = c('CH103' = 'CH103 (H)', 'CH103L' = 'CH103 (L)',
                              'VRC26' = 'VRC26 (H)', 'VRC26L' = 'VRC26 (L)',
                              'VRC01_13' = 'VRC01-13 (H)','VRC01_01' = 'VRC01-01 (H)',
                              'VRC01_19' = 'VRC01-19 (H)')
                  ) +
  guides(fill=guide_legend(title=NULL))
  
  

# =======================  BRANCH-LEVEL CHANGES IN MUTABILITY DUE TO SYN. SUBSTITUTIONS - DATA VS. MODELS  =========================
pl_syn <- ggplot(subset(combined_dataframe, substitution_class == 'syn_only'), 
                 aes(x = lineage, y = mean_contrast_true)) +
  geom_hline(yintercept=0,linetype = 2) +
  theme_bw() +
  #theme_classic() + 
  ylab('Synonymous change in mutability') + 
  xlab('Lineage') +
  theme(legend.position = 'top',
        legend.text=element_text(size=11)
  ) +
  
  scale_y_continuous(limits = c(-0.01,0.01)) +
  
  # Observed
  geom_point(aes(y=mean_contrast_true), shape = 21, fill = 'gray80', alpha = 0.8,
             size = 3) + 
  
  # Simulated under S5F model
  geom_linerange(aes(ymin=contrast_simulation_S5F_S5F_llim, 
                     ymax = contrast_simulation_S5F_S5F_ulim), colour = 'red4',
                 position = position_nudge(x=-0.1,y=0)) +
  geom_point(aes(y=mean_contrast_simulation_S5F_S5F), shape = 21, fill = 'firebrick1', alpha = 0.8,
             size = 3,position = position_nudge(x=-0.1,y=0)) +
  
  # Simulated under uniform model
  geom_linerange(aes(ymin=contrast_simulation_uniform_S5F_llim, 
                     ymax = contrast_simulation_uniform_S5F_ulim), colour = 'blue4',
                 position = position_nudge(x=0.1,y=0)) +
  geom_point(aes(y=mean_contrast_simulation_uniform_S5F), shape = 21, fill = 'royalblue3', alpha = 0.8,
             size = 3, position = position_nudge(x=0.1,y=0)) +
  
  # Simulated under codon-position model
  geom_linerange(aes(ymin=contrast_simulation_CP_S5F_llim, 
                     ymax = contrast_simulation_CP_S5F_ulim), colour = 'skyblue',
                 position = position_nudge(x=0.2,y=0)) +
  geom_point(aes(y=mean_contrast_simulation_CP_S5F), shape = 21, fill = 'skyblue', alpha = 0.8,
             size = 3, position = position_nudge(x=0.2,y=0)) +
  

  scale_x_discrete(labels = c('CH103' = 'CH103 (H)', 'CH103L' = 'CH103 (L)',
                              'VRC26' = 'VRC26 (H)', 'VRC26L' = 'VRC26 (L)',
                              'VRC01_13' = 'VRC01-13 (H)','VRC01_01' = 'VRC01-01 (H)',
                              'VRC01_19' = 'VRC01-19 (H)'
  )) + 

  # Dummy points (plotted outside of plotting region, only purpose is to quickly produce ggplot legend)
  geom_point(data = data.frame(lineage = c('CH103','CH103','CH103','CH103','CH103L','CH103L','CH103L','CH103L'),
                               change = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),
                               type = factor(c('observed','S5F model','uniform model','codon-position model',
                                               'observed','S5F model','uniform model','codon-position model'),
                                             levels = c('observed','S5F model','uniform model','codon-position model'))
                               ),
             aes(x = lineage, y = change, fill = type), shape = 21, size = 3
             
             ) +
  scale_fill_manual(values = c('gray80', 'firebrick1', 'royalblue3','skyblue'),
                    labels = c('observed','S5F\nmodel','Uniform\nmodel','Codon-position\nmodel')) +
  
  guides(fill=guide_legend(title=NULL))
  

#pdf('temp.pdf', height=6, width=6)
#plot(pl_syn)
#dev.off()

# ================================================ COMBINING PLOTS ====================================================
combined_plot <- plot_grid(total_changes_pl, pl_syn, labels = c("a)", "b)"),
                           rel_widths = c(1.5,1), label_size = 20, nrow =2)

save_plot("S_NS_mutability_changes.pdf", combined_plot,
          base_height = 10, base_width = 7.5
)