# Produces figure for the manuscript showing the results of the simulations on MCC trees
# Uses the following MCMC chains: CH103_con_run1a, CH103L_con_run1a, VRC26int_con_run1a, VRC26L_con_run1a, VRC01_01_log_run1a, VRC01_13_log_run1a, VRC01_19_log_run1a,

library('ggplot2')
library('cowplot')
library('grid')
library('lattice')
source('ggplot_parameters.R')

results_directory <- '../results/S_NS_mutability_changes/observed_lineages/'

observed_results_files <- c(CH103 = 'CH103_constant/CH103_con_run1a_observed_MCC.csv',
                  CH103L = 'CH103L_constant/CH103L_con_run1a_observed_MCC.csv',
                  VRC26 = 'VRC26int_constant/VRC26int_con_run1a_observed_MCC.csv',
                  VRC26L = 'VRC26L_constant/VRC26L_con_run1a_observed_MCC.csv',
                  VRC01_13 = 'VRC01_13_logistic/VRC01_13_log_run1a_observed_MCC.csv',
                  VRC01_01 = 'VRC01_01_logistic/VRC01_01_log_run1a_observed_MCC.csv',
                  VRC01_19 = 'VRC01_19_logistic/VRC01_19_log_run1a_observed_MCC.csv'
)

for(i in 1:length(observed_results_files)){
  observed_results_files[i] <- paste(results_directory, observed_results_files[i], sep = '')
}

constrained_simulations_files <- gsub('observed_MCC','simulated_MCC_constrained', observed_results_files)
unconstrained_simulations_files <- gsub('observed_MCC','simulated_MCC_unconstrained', observed_results_files)

dataframe_list_obs <- lapply(observed_results_files, FUN = read.csv, header = T)
dataframe_list_sim_constrained <- lapply(constrained_simulations_files, FUN = read.csv, header = T)
dataframe_list_sim_unconstrained <- lapply(unconstrained_simulations_files, FUN = read.csv, header = T)


# ====== COMBINED DATAFRAME WITH RESULTS FOR ALL LINEAGES =======
lineage_vector <- c()
metric_vector <- c()
region_vector <- c()
substitution_class <- c()

mean_contrast_true <- c()

# Vectors storing mean and CIs for constrained simulation results
mean_contrast_simulation_constrained_S5F_S5F <- c()
contrast_simulation_constrained_S5F_S5F_llim <- c()
contrast_simulation_constrained_S5F_S5F_ulim <- c()

mean_contrast_simulation_constrained_uniform_S5F <- c()
contrast_simulation_constrained_uniform_S5F_llim <- c()
contrast_simulation_constrained_uniform_S5F_ulim <- c()

mean_contrast_simulation_constrained_CP_S5F <- c()
contrast_simulation_constrained_CP_S5F_llim <- c()
contrast_simulation_constrained_CP_S5F_ulim <- c()

# Vectors storing mean and CIs for unconstrained simulation results
mean_contrast_simulation_unconstrained_S5F_S5F <- c()
contrast_simulation_unconstrained_S5F_S5F_llim <- c()
contrast_simulation_unconstrained_S5F_S5F_ulim <- c()

mean_contrast_simulation_unconstrained_uniform_S5F <- c()
contrast_simulation_unconstrained_uniform_S5F_llim <- c()
contrast_simulation_unconstrained_uniform_S5F_ulim <- c()

mean_contrast_simulation_unconstrained_CP_S5F <- c()
contrast_simulation_unconstrained_CP_S5F_llim <- c()
contrast_simulation_unconstrained_CP_S5F_ulim <- c()

#region <- 'WS'
for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  for(region in c('WS','FR','CDR')){  
    clone_dataframe_obs <- dataframe_list_obs[[clone]]
    
    #for(metric in c('S5F','HS','OHS')){
    for(metric in c('S5F', 'logS5F')){
      for(sub_type in c('total','syn_only','nonsyn_only')){
          
        contrast_syn <- clone_dataframe_obs[,paste(metric, '_change_syn',ifelse(region=='WS','',
                                                         paste('_', region, sep ='')), sep ='')]
        
        contrast_nonsyn <-clone_dataframe_obs[,paste(metric, '_change_nonsyn',ifelse(region=='WS','',
                                                                           paste('_', region, sep ='')), sep ='')]
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
        
        substitution_class <- c(substitution_class, sub_type)
        metric_vector <- c(metric_vector, metric)
        lineage_vector <- c(lineage_vector, clone)
        region_vector <- c(region_vector, region)
        
        # Get simulation results
        for(simulation_type in c('constrained','unconstrained')){
          
          if(simulation_type == 'constrained'){
            clone_dataframe_sim <- dataframe_list_sim_constrained[[clone]]
          }else{
            clone_dataframe_sim <- dataframe_list_sim_unconstrained[[clone]]
          }
          
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
            
          }
          ##### COMPUTE AVERAGE OF MEAN CONTRASTS AND 95% INTERVAL
          if(simulation_type == 'constrained'){
            mean_contrast_simulation_constrained_S5F_S5F <- c(mean_contrast_simulation_constrained_S5F_S5F, mean(simulated_means_S5F_S5F))
            mean_contrast_simulation_constrained_uniform_S5F <- c(mean_contrast_simulation_constrained_uniform_S5F, mean(simulated_means_uniform_S5F))
            mean_contrast_simulation_constrained_CP_S5F <- c(mean_contrast_simulation_constrained_CP_S5F, mean(simulated_means_CP_S5F))
            
            # 95% INTERVAL
            contrast_simulation_constrained_S5F_S5F_llim <- c(contrast_simulation_constrained_S5F_S5F_llim, quantile(simulated_means_S5F_S5F,0.025))
            contrast_simulation_constrained_S5F_S5F_ulim <- c(contrast_simulation_constrained_S5F_S5F_ulim,quantile(simulated_means_S5F_S5F,0.975))
            
            contrast_simulation_constrained_uniform_S5F_llim <- c(contrast_simulation_constrained_uniform_S5F_llim, quantile(simulated_means_uniform_S5F,0.025))
            contrast_simulation_constrained_uniform_S5F_ulim <- c(contrast_simulation_constrained_uniform_S5F_ulim,quantile(simulated_means_uniform_S5F,0.975))
            
            contrast_simulation_constrained_CP_S5F_llim <- c(contrast_simulation_constrained_CP_S5F_llim, quantile(simulated_means_CP_S5F,0.025))
            contrast_simulation_constrained_CP_S5F_ulim <- c(contrast_simulation_constrained_CP_S5F_ulim,quantile(simulated_means_CP_S5F,0.975))
          }else{
            mean_contrast_simulation_unconstrained_S5F_S5F <- c(mean_contrast_simulation_unconstrained_S5F_S5F, mean(simulated_means_S5F_S5F))
            mean_contrast_simulation_unconstrained_uniform_S5F <- c(mean_contrast_simulation_unconstrained_uniform_S5F, mean(simulated_means_uniform_S5F))
            mean_contrast_simulation_unconstrained_CP_S5F <- c(mean_contrast_simulation_unconstrained_CP_S5F, mean(simulated_means_CP_S5F))
            
            # 95% INTERVAL
            contrast_simulation_unconstrained_S5F_S5F_llim <- c(contrast_simulation_unconstrained_S5F_S5F_llim, quantile(simulated_means_S5F_S5F,0.025))
            contrast_simulation_unconstrained_S5F_S5F_ulim <- c(contrast_simulation_unconstrained_S5F_S5F_ulim,quantile(simulated_means_S5F_S5F,0.975))
            
            contrast_simulation_unconstrained_uniform_S5F_llim <- c(contrast_simulation_unconstrained_uniform_S5F_llim, quantile(simulated_means_uniform_S5F,0.025))
            contrast_simulation_unconstrained_uniform_S5F_ulim <- c(contrast_simulation_unconstrained_uniform_S5F_ulim,quantile(simulated_means_uniform_S5F,0.975))
            
            contrast_simulation_unconstrained_CP_S5F_llim <- c(contrast_simulation_unconstrained_CP_S5F_llim, quantile(simulated_means_CP_S5F,0.025))
            contrast_simulation_unconstrained_CP_S5F_ulim <- c(contrast_simulation_unconstrained_CP_S5F_ulim,quantile(simulated_means_CP_S5F,0.975))
          }
        }
      }
      
    }
  }
}

lineage_vector <- factor(lineage_vector, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
metric_vector <- factor(metric_vector, levels = c('S5F','logS5F'))
region_vector <- factor(region_vector, levels = c('WS','FR','CDR'))

combined_dataframe <- data.frame(lineage=lineage_vector, metric = metric_vector, region = region_vector, substitution_class, mean_contrast_true, 
                                 mean_contrast_simulation_constrained_S5F_S5F, contrast_simulation_constrained_S5F_S5F_llim, contrast_simulation_constrained_S5F_S5F_ulim,
                                 mean_contrast_simulation_constrained_uniform_S5F, contrast_simulation_constrained_uniform_S5F_llim, contrast_simulation_constrained_uniform_S5F_ulim,
                                 mean_contrast_simulation_constrained_CP_S5F, contrast_simulation_constrained_CP_S5F_llim, contrast_simulation_constrained_CP_S5F_ulim,
                                 mean_contrast_simulation_unconstrained_S5F_S5F, contrast_simulation_unconstrained_S5F_S5F_llim, contrast_simulation_unconstrained_S5F_S5F_ulim,
                                 mean_contrast_simulation_unconstrained_uniform_S5F, contrast_simulation_unconstrained_uniform_S5F_llim, contrast_simulation_unconstrained_uniform_S5F_ulim,
                                 mean_contrast_simulation_unconstrained_CP_S5F, contrast_simulation_unconstrained_CP_S5F_llim, contrast_simulation_unconstrained_CP_S5F_ulim
                                 )
                                 
# DATAFRAME WITH SUMMED CHANGES IN MUTABILITY ACROSS THE TREE
lineage_vector <- c()
metric_vector <- c()
region_vector <- c()
cumulative_change <- c()
substitution_class <- c()
group <- c()

for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){

  clone_dataframe_obs <- dataframe_list_obs[[clone]] 
  
  for(region in c('WS','FR','CDR')){
    for(metric in c('S5F','logS5F')){
      
      region_id <- ifelse(region == 'WS','', paste('_',region,sep=''))
      for(class in c('syn_only','nonsyn_only')){
    
        if(class == 'syn_only'){
          contrasts <- clone_dataframe_obs[,paste(metric, '_change_syn',region_id, sep='')]
        }else{
          contrasts <- clone_dataframe_obs[,paste(metric, '_change_nonsyn',region_id, sep='')]
        }
        
        lineage_vector <- c(lineage_vector, clone)
        region_vector <- c(region_vector, region)
        metric_vector <- c(metric_vector, metric)
        substitution_class <- c(substitution_class, class)
        cumulative_change <- c(cumulative_change, sum(contrasts))
        
        group <- c(group, paste(region, class, sep ='_'))
      }
    }
  }
} 
lineage_vector <- factor(lineage_vector, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                  'VRC01_01','VRC01_19'))
group <- factor(group, levels = c('WS_syn_only','WS_nonsyn_only','FR_syn_only','FR_nonsyn_only',
                                  'CDR_syn_only','CDR_nonsyn_only'))
region_vector <- factor(region_vector, levels = c('WS','FR','CDR'))
metric_vector <- factor(metric_vector, levels = c('S5F','logS5F'))
substitution_class <- factor(substitution_class, levels = c('syn_only','nonsyn_only'))
total_change_dataframe <- data.frame(lineage=lineage_vector, region = region_vector, metric = metric_vector, 
                                     substitution_class,cumulative_change, group)


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
total_changes_plot <- function(metric){
  subset_dataframe <- total_change_dataframe[total_change_dataframe$region != 'WS' &
                                              total_change_dataframe$metric == metric, ]
  
  ylabel <- switch(metric,
                   S5F = 'Cumulative change in mean S5F mutability',
                   logS5F = 'Cumulative change in mean log-S5F mutability')
  
  pl <- ggplot(subset_dataframe,aes(y=cumulative_change,x=region)) +

  geom_col(aes(fill = substitution_class), 
           width = 0.7, position=position_stack()) + 
  
  facet_grid(~lineage) +
  geom_hline(yintercept=0,linetype = 2) +
  theme_bw() +
  #theme_classic() + 
  ylab(ylabel) + 
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
  
  return(pl)
}

total_changes_plot_S5F <- total_changes_plot('S5F')
total_changes_plot_logS5F <- total_changes_plot('logS5F')

pdf('geom_and_log_S5F_plots/S_NS_total_changes_logS5F.pdf', height = 5, width = 8)
  plot(total_changes_plot_logS5F)
dev.off()

# =======================  BRANCH-LEVEL CHANGES IN MUTABILITY DUE - DATA VS. MODELS  =========================
# Basic plotting function
base_plot <- function(dataframe, substitution_class, region, simulation_type, ylims = c(-0.03,0.03)){
  
  dataframe <- dataframe[dataframe$substitution_class == substitution_class & dataframe$region == region,]
  
  ylabel_subst_class <- switch(substitution_class,
                               syn_only = 'Synonymous change in mutability',
                               nonsyn_only = 'Non-synonymous change in mutability',
                               total = 'Total change in mutability')
  
  # Generate dataframe to pass to ggplot
  ggplot_dataframe <- data.frame(lineage = dataframe$lineage, mean_contrast_true = dataframe$mean_contrast_true,
                                 mean_contrast_simulation_S5F_S5F = dataframe[,paste('mean_contrast_simulation_',simulation_type,'_S5F_S5F', sep = '')],
                                 contrast_simulation_S5F_S5F_llim = dataframe[,paste('contrast_simulation_',simulation_type,'_S5F_S5F_llim', sep = '')],
                                 contrast_simulation_S5F_S5F_ulim = dataframe[,paste('contrast_simulation_',simulation_type,'_S5F_S5F_ulim', sep = '')],
                                 mean_contrast_simulation_uniform_S5F = dataframe[,paste('mean_contrast_simulation_',simulation_type,'_uniform_S5F', sep = '')],
                                 contrast_simulation_uniform_S5F_llim = dataframe[,paste('contrast_simulation_',simulation_type,'_uniform_S5F_llim', sep = '')],
                                 contrast_simulation_uniform_S5F_ulim = dataframe[,paste('contrast_simulation_',simulation_type,'_uniform_S5F_ulim', sep = '')],
                                 mean_contrast_simulation_CP_S5F = dataframe[,paste('mean_contrast_simulation_',simulation_type,'_CP_S5F', sep = '')],
                                 contrast_simulation_CP_S5F_llim = dataframe[,paste('contrast_simulation_',simulation_type,'_CP_S5F_llim', sep = '')],
                                 contrast_simulation_CP_S5F_ulim = dataframe[,paste('contrast_simulation_',simulation_type,'_CP_S5F_ulim', sep = '')]
                                 )
  
  pl <- ggplot(ggplot_dataframe, aes(x = lineage, y = mean_contrast_true)) +
    geom_hline(yintercept=0,linetype = 2) +
    theme_bw() +
    ylab(paste(ylabel_subst_class, ' (', region,')', sep = '')) + 
    xlab('Lineage') +
    theme(legend.position = 'top',
          legend.text=element_text(size=11)
    ) +
    
    scale_y_continuous(limits = ylims) +
    
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
                                 change = c(10,10,10,10,10,10,10,10),
                                 type = factor(c('observed','S5F model','uniform model','codon-position model',
                                                 'observed','S5F model','uniform model','codon-position model'),
                                               levels = c('observed','S5F model','uniform model','codon-position model'))
    ),
    aes(x = lineage, y = change, fill = type), shape = 21, size = 3
    
    ) +
    scale_fill_manual(values = c('gray80', 'firebrick1', 'royalblue3','skyblue'),
                      labels = c('observed','S5F\nmodel','Uniform\nmodel','Codon-position\nmodel')) +
    
    guides(fill=guide_legend(title=NULL))
  
    
  return(pl)
}

# Changes compared to aa-unconstrained simulations
pl_syn_WS_unconstrained <- base_plot(combined_dataframe, substitution_class = 'syn_only', region = 'WS',
                                     simulation_type = 'unconstrained', ylims = c(-0.03,0.02))
pl_nonsyn_WS_unconstrained <- base_plot(combined_dataframe, substitution_class = 'nonsyn_only', region = 'WS',
                                        simulation_type = 'unconstrained', ylims = c(-0.03,0.02))
pl_syn_FR_unconstrained <- base_plot(combined_dataframe, substitution_class = 'syn_only', region = 'FR',
                                     simulation_type = 'unconstrained', ylims = c(-0.03,0.02))
pl_nonsyn_FR_unconstrained <- base_plot(combined_dataframe, substitution_class = 'nonsyn_only', region = 'FR',
                                        simulation_type = 'unconstrained', ylims = c(-0.03,0.02))
pl_syn_CDR_unconstrained <- base_plot(combined_dataframe, substitution_class = 'syn_only', region = 'CDR',
                                     simulation_type = 'unconstrained', ylims = c(-0.035,0.01))
pl_nonsyn_CDR_unconstrained <- base_plot(combined_dataframe, substitution_class = 'nonsyn_only', region = 'CDR',
                                        simulation_type = 'unconstrained', ylims = c(-0.035,0.01))

unconstrained_plot <- plot_grid(pl_syn_WS_unconstrained, pl_nonsyn_WS_unconstrained,
                                pl_syn_FR_unconstrained, pl_nonsyn_FR_unconstrained, 
                                 pl_syn_CDR_unconstrained, pl_nonsyn_CDR_unconstrained, 
                                labels = c("Whole sequence", "","FRs","","CDRs",""), 
                                label_size = 15, nrow =3, hjust = 0)

save_plot('randomization_figures/unconstrained.pdf',
          unconstrained_plot,
          base_height = 18, base_width = 14
          )


# ================================================ COMBINING PLOTS ====================================================
#combined_plot <- plot_grid(total_changes_pl, pl_syn, pl_nonsyn_FR, pl_nonsyn_CDR, 
#                           labels = c("a)", "b)","c)","d)"),
#                           label_size = 20, nrow =2)
print('WARNING: PLOT NOT SAVED')
#save_plot("S_NS_mutability_changes.pdf", combined_plot,
#          base_height = 10, base_width = 14
#)