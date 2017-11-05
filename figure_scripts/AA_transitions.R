# Produces figure for the manuscript showing the results of the analysis of S. and NS mutability changes in MCC trees and Liao ancestors
# Uses the following MCMC chains: CH103_con_run1a, CH103L_con_run1a, VRC26int_con_run1a, VRC26L_con_run1a, VRC01_01_log_run1a, VRC01_13_log_run1a, VRC01_19_log_run1a,

library('ggplot2')
library('cowplot')
library('grid')
library('lattice')
source('ggplot_parameters.R')

results_directory <- '../results/S_NS_mutability_changes/observed_lineages/'

observed_results_files <- c(CH103 = 'CH103_constant/CH103_con_run1a_aa_transitions_obs.csv',
                  CH103L = 'CH103L_constant/CH103L_con_run1a_aa_transitions_obs.csv',
                  VRC26 = 'VRC26int_constant/VRC26int_con_run1a_aa_transitions_obs.csv',
                  VRC26L = 'VRC26L_constant/VRC26L_con_run1a_aa_transitions_obs.csv',
                  VRC01_13 = 'VRC01_13_logistic/VRC01_13_log_run1a_aa_transitions_obs.csv',
                  VRC01_01 = 'VRC01_01_logistic/VRC01_01_log_run1a_aa_transitions_obs.csv',
                  VRC01_19 = 'VRC01_19_logistic/VRC01_19_log_run1a_aa_transitions_obs.csv'
)
# observed_results_files <- c(CH103 = 'CH103L_constant/CH103L_con_run1a_aa_transitions_obs.csv',
#                   CH103L = 'CH103L_constant/CH103L_con_run1a_aa_transitions_obs.csv',
#                   VRC26 = 'CH103L_constant/CH103L_con_run1a_aa_transitions_obs.csv',
#                   VRC26L = 'CH103L_constant/CH103L_con_run1a_aa_transitions_obs.csv',
#                   VRC01_13 = 'CH103L_constant/CH103L_con_run1a_aa_transitions_obs.csv',
#                   VRC01_01 = 'CH103L_constant/CH103L_con_run1a_aa_transitions_obs.csv',
#                   VRC01_19 = 'CH103L_constant/CH103L_con_run1a_aa_transitions_obs.csv'
# )

for(i in 1:length(observed_results_files)){
  observed_results_files[i] <- paste(results_directory, observed_results_files[i], sep = '')
}

unconstrained_simulations_files <- gsub('_obs','_unconstrained', observed_results_files)

dataframe_list_obs <- lapply(observed_results_files, FUN = read.csv, header = T)
dataframe_list_sim_unconstrained <- lapply(unconstrained_simulations_files, FUN = read.csv, header = T,row.names = NULL)


lineage_vector <- c()
value_type_vector <- c()
total_change_in_mean_log_S5F <- c()
total_change_in_mean_log_S5F_llim <- c()
total_change_in_mean_log_S5F_ulim <- c()
ancestor_aa_vector <- c()
descendant_aa_vector <- c()
ntrans_vector <- c()
ntrans_llim_vector <- c()
ntrans_ulim_vector <- c()

for(lineage in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  print(lineage)
  for(value_type in c('observed', 'S5F', 'CP','uniform')){
    
    if(value_type == 'observed'){
      dataframe <- dataframe_list_obs[[lineage]]
    }else{
      dataframe <- dataframe_list_sim_unconstrained[[lineage]] 
      dataframe <- dataframe[dataframe$mutability_model == value_type,]
    }
    
    # Aggregate data across all nodes:
    for(ancestor_aa in unique(dataframe$ancestor_aa)){
      for(descendant_aa in unique(dataframe$descendant_aa)){
        pair_data <- dataframe[dataframe$ancestor_aa == ancestor_aa & dataframe$descendant_aa == descendant_aa,]
        
        # If this transition was observed at least once:
        if(nrow(pair_data)>0){
        
          # If value_type is a simulation model, average over replicates, summing across nodes within replicates:
          if(value_type != 'observed'){
            replicate_ntrans <- c()
            replicate_meanLogS5F_changes <- c()
            for(replicate in unique(pair_data$replicate)){
              replicate_ntrans <- c(replicate_ntrans, sum(pair_data[pair_data$replicate == replicate,
                                                                'n_trans']))
              replicate_meanLogS5F_changes <- c(replicate_meanLogS5F_changes, 
                                                sum(pair_data[pair_data$replicate == replicate,
                                                              'total_meanlogS5Fchange']))
            }
            ntrans <- mean(replicate_ntrans)
            ntrans_llim <- quantile(replicate_ntrans,0.025)
            ntrans_ulim <-quantile(replicate_ntrans,0.975)
            
            meanlogS5F_change <- mean(replicate_meanLogS5F_changes)
            meanlogS5F_change_llim <- quantile(replicate_meanLogS5F_changes, 0.025)
            meanlogS5F_change_ulim <- quantile(replicate_meanLogS5F_changes, 0.975)
            
          }else{
            # If values are from the observed MCC, just sum across nodes.
            ntrans <- sum(pair_data[,'n_trans'])
            meanlogS5F_change <- sum(pair_data[,'total_meanlogS5Fchange'])
            ntrans_llim <- NA
            ntrans_ulim <- NA
            meanlogS5F_change_llim <- NA
            meanlogS5F_change_ulim <- NA
          }
          
          lineage_vector <- c(lineage_vector, lineage)
          value_type_vector <- c(value_type_vector, value_type)
          ancestor_aa_vector <- c(ancestor_aa_vector, ancestor_aa)
          descendant_aa_vector <- c(descendant_aa_vector, descendant_aa)
          ntrans_vector <- c(ntrans_vector, ntrans)
          ntrans_llim_vector <- c(ntrans_llim_vector, ntrans_llim)
          ntrans_ulim_vector <- c(ntrans_ulim_vector, ntrans_ulim)
          
          total_change_in_mean_log_S5F <- c(total_change_in_mean_log_S5F, meanlogS5F_change)
          total_change_in_mean_log_S5F_llim <- c(total_change_in_mean_log_S5F_llim, meanlogS5F_change_llim)
          total_change_in_mean_log_S5F_ulim <- c(total_change_in_mean_log_S5F_ulim, meanlogS5F_change_ulim)
          
        }
      }
    }
  }
}


lineage_vector <- factor(lineage_vector, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                  'VRC01_01','VRC01_19'))

combined_dataframe <- data.frame(lineage = lineage_vector, model = value_type_vector,
                                 ancestor_aa = ancestor_aa_vector, descendant_aa = descendant_aa_vector,
                                 n_trans = ntrans_vector, n_trans_llim = ntrans_llim_vector, n_trans_ulim = ntrans_ulim_vector,
                                 total_change_in_mean_log_S5F,total_change_in_mean_log_S5F_llim, total_change_in_mean_log_S5F_ulim
                                 )

  # transition_freqs_pl <- ggplot(subset(combined_dataframe, model %in% c('observed','S5F') & 
  #                                        descendant_aa != '*' & ancestor_aa != '*'), aes(y = ancestor_aa, x = descendant_aa)) +
  #   #geom_tile(aes(fill = n_trans)) +
  #   geom_tile(aes(fill = total_change_in_mean_log_S5F)) +
  #   facet_grid(model~lineage) +
  #   scale_fill_gradient(low = 'red', high = 'blue')
  # 
  # pdf('AA_transition_figs/AA_transition_freqs.pdf', height = 7, width = 21)
  # plot(transition_freqs_pl)
  # dev.off()

# DATAFRAME FOR SCATTERPLOTS
lineage_vector <- c()
model_vector <- c()
ancestor_aa_vector <- c()
descendant_aa_vector <- c()
sim_ntrans_vector <- c()
obs_ntrans_vector <- c()
sim_ntrans_llim_vector <- c()
sim_ntrans_ulim_vector <- c()
obs_total_change_in_mean_log_S5F <- c()
sim_total_change_in_mean_log_S5F <- c()
sim_total_change_in_mean_log_S5F_llim <- c()
sim_total_change_in_mean_log_S5F_ulim <- c()
for(lineage in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  for(model in c('S5F', 'CP','uniform')){
    sim_dataframe <- combined_dataframe[combined_dataframe$model == model &
                                          combined_dataframe$lineage == lineage,]
   
    for(ancestor_aa in unique(sim_dataframe$ancestor_aa)){
      for(descendant_aa in unique(sim_dataframe$descendant_aa)){
        
        obs_ntrans <- combined_dataframe$n_trans[combined_dataframe$lineage == lineage &
                                                   combined_dataframe$model == 'observed'&
                                                   combined_dataframe$ancestor_aa == ancestor_aa &
                                                   combined_dataframe$descendant_aa == descendant_aa]
        obs_ntrans <- ifelse(length(obs_ntrans) > 0, obs_ntrans, 0)
        
        obs_change <- combined_dataframe$total_change_in_mean_log_S5F[combined_dataframe$lineage == lineage&
                                                                        combined_dataframe$model == 'observed'&
                                                                        combined_dataframe$ancestor_aa == ancestor_aa &
                                                                        combined_dataframe$descendant_aa == descendant_aa]
        obs_change <- ifelse(length(obs_change) > 0, obs_change, 0)
        
        sim_ntrans <- sim_dataframe$n_trans[sim_dataframe$ancestor_aa == ancestor_aa &
                                              sim_dataframe$descendant_aa == descendant_aa]
        sim_ntrans <- ifelse(length(sim_ntrans) > 0, sim_ntrans, 0)
        
        sim_change <- sim_dataframe$total_change_in_mean_log_S5F[sim_dataframe$ancestor_aa == ancestor_aa &
                                                            sim_dataframe$descendant_aa == descendant_aa]
        sim_change <- ifelse(length(sim_change) > 0, sim_change, 0)
        
  
        sim_ntrans_llim <- sim_dataframe$n_trans_llim[sim_dataframe$ancestor_aa == ancestor_aa &
                                              sim_dataframe$descendant_aa == descendant_aa]
        sim_ntrans_llim <- ifelse(length(sim_ntrans_llim) > 0, sim_ntrans_llim, 0)
  
        sim_ntrans_ulim <- sim_dataframe$n_trans_ulim[sim_dataframe$ancestor_aa == ancestor_aa &
                                                        sim_dataframe$descendant_aa == descendant_aa]
        sim_ntrans_ulim <- ifelse(length(sim_ntrans_ulim) > 0, sim_ntrans_ulim, 0)
        
        sim_change_llim <- sim_dataframe$total_change_in_mean_log_S5F_llim[sim_dataframe$ancestor_aa == ancestor_aa &
                                                                   sim_dataframe$descendant_aa == descendant_aa]
        sim_change_llim <- ifelse(length(sim_change_llim) > 0, sim_change_llim, 0)
        
        sim_change_ulim <- sim_dataframe$total_change_in_mean_log_S5F_ulim[sim_dataframe$ancestor_aa == ancestor_aa &
                                                                             sim_dataframe$descendant_aa == descendant_aa]
        sim_change_ulim <- ifelse(length(sim_change_ulim) > 0, sim_change_ulim, 0)
        
        lineage_vector <- c(lineage_vector, lineage)
        model_vector <- c(model_vector, model)
        sim_ntrans_vector <- c(sim_ntrans_vector, sim_ntrans)
        obs_ntrans_vector <- c(obs_ntrans_vector, obs_ntrans)
        ancestor_aa_vector <- c(ancestor_aa_vector, ancestor_aa)
        descendant_aa_vector <- c(descendant_aa_vector, descendant_aa)
        sim_ntrans_llim_vector <- c(sim_ntrans_llim_vector, sim_ntrans_llim)
        sim_ntrans_ulim_vector <- c(sim_ntrans_ulim_vector, sim_ntrans_ulim)
        obs_total_change_in_mean_log_S5F <- c(obs_total_change_in_mean_log_S5F, obs_change)
        sim_total_change_in_mean_log_S5F <- c(sim_total_change_in_mean_log_S5F, sim_change)
        sim_total_change_in_mean_log_S5F_llim <- c(sim_total_change_in_mean_log_S5F_llim, sim_change_llim)
        sim_total_change_in_mean_log_S5F_ulim <- c(sim_total_change_in_mean_log_S5F_ulim, sim_change_ulim)
      
      }
    }
  }
}
lineage_vector <- factor(lineage_vector, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                                    'VRC01_01','VRC01_19'))

scatterplot_dataframe <- data.frame(lineage = lineage_vector, 
                                    transition = paste(ancestor_aa_vector, descendant_aa_vector, sep = '-'),
                                    model = model_vector, obs_ntrans = obs_ntrans_vector,
                                    obs_total_change_in_mean_log_S5F, sim_ntrans = sim_ntrans_vector,
                                    sim_ntrans_llim = sim_ntrans_llim_vector, sim_ntrans_ulim = sim_ntrans_ulim_vector,
                                    sim_total_change_in_mean_log_S5F, sim_total_change_in_mean_log_S5F_llim, sim_total_change_in_mean_log_S5F_ulim
                                    )

scatterplot_dataframe <- data.frame(scatterplot_dataframe, 
                                    sim_average_change_in_mean_log_S5F = 
                                      scatterplot_dataframe$sim_total_change_in_mean_log_S5F / scatterplot_dataframe$sim_ntrans,
                                    obs_average_change_in_mean_log_S5F = 
                                      scatterplot_dataframe$obs_total_change_in_mean_log_S5F / scatterplot_dataframe$obs_ntrans
)

obs_vs_sim_ntrans <- function(lineage,scat_dataframe = scatterplot_dataframe){
  lineage_dataframe = scat_dataframe[scat_dataframe$lineage == lineage &
                                     scat_dataframe$model == 'S5F',]
  
  pl <- ggplot(data = lineage_dataframe, aes(x = obs_ntrans, y = sim_ntrans)) + 
    #geom_point() +
    ggplot_theme +
    theme(plot.margin = margin(20,10,10,1,'pt')) +
    theme(legend.position = 'top') +
    theme(legend.text = element_text(size = 8)) +
    theme(legend.title = element_text(size = 10)) +
    xlab('Observed number of transitions') +
    ylab('Mean number of transitions under S5F model') +
    geom_label(aes(label = transition, fill = sim_average_change_in_mean_log_S5F),
    size = 2, fontface = "bold", alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    scale_fill_gradient2(low = 'red', mid = 'gray',high = 'blue',
                        name = 'Average change in mean log-S5F') +
    guides(fill = guide_colorbar(title.position = 'top', title.hjust = -1,
                                barwidth = 7, barheight = 0.8, nbin = 4))
  return(pl)
}

obs_vs_sim_ntrans_plots <- sapply(c('CH103','CH103L','VRC26','VRC26L','VRC01_13','VRC01_01','VRC01_19'),
                                  FUN = obs_vs_sim_ntrans, simplify = FALSE)

lineage_labels <- c()
for(name in names(obs_vs_sim_ntrans_plots)){
  lineage_labels <- c(lineage_labels, switch(name,
       CH103 = 'CH103 (H)',
       CH103L = 'CH103 (L)',
       VRC26 = 'VRC26 (H)',
       VRC26L = 'VRC26 (L)',
       VRC01_13 = 'VRC01-13 (H)',
       VRC01_01 = 'VRC01-01 (H)',
       VRC01_19 = 'VRC01-19 (H)'
       ))
}
sim_vs_obs_ntrans_figure <- plot_grid(plotlist = obs_vs_sim_ntrans_plots, nrow = 3, ncol = 3, 
                          labels = lineage_labels, label_size = 16)

pdf('sim_vs_obs_ntrans.pdf', height = 16, width = 15)
plot(sim_vs_obs_ntrans_figure)
dev.off()



# Summary dataframe
lineage_vector <- c()
mean_change_of_top_sim_transitions_vector <- c()
mean_change_of_top_obs_transitions_vector <- c()
top_sim_transitions_fraction_in_sim_vector <- c()
top_sim_transitions_fraction_in_obs_vector <- c()

for(lineage in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  
  transition <- scatterplot_dataframe[scatterplot_dataframe$lineage == lineage & scatterplot_dataframe$model == 'S5F',
                                      'transition']
  sim_ntrans <- scatterplot_dataframe[scatterplot_dataframe$lineage == lineage & scatterplot_dataframe$model == 'S5F',
                                       'sim_ntrans']
  obs_ntrans <- scatterplot_dataframe[scatterplot_dataframe$lineage == lineage & scatterplot_dataframe$model == 'S5F',
                                       'obs_ntrans']
  sim_changes <- scatterplot_dataframe[scatterplot_dataframe$lineage == lineage & scatterplot_dataframe$model == 'S5F',
                                       'sim_average_change_in_mean_log_S5F']
  obs_changes <- scatterplot_dataframe[scatterplot_dataframe$lineage == 'CH103' & scatterplot_dataframe$model == 'S5F',
                                       'obs_average_change_in_mean_log_S5F']
  
  top_sim_transitions <- transition[sim_ntrans %in% sort(sim_ntrans,decreasing=T)[1:10]]
  top_obs_transitions <- transition[obs_ntrans %in% sort(obs_ntrans,decreasing=T)[1:10]]
  
  mean_change_of_top_sim_transitions <- mean(sim_changes[transition %in% top_sim_transitions],na.rm =T)
  mean_change_of_top_obs_transitions <- mean(obs_changes[transition %in% top_obs_transitions], na.rm = T)
  
  top_sim_transitions_fraction_in_sim <- 0
  top_sim_transitions_fraction_in_obs <- 0
  
  for(trans in top_sim_transitions){
    obs_value <- obs_ntrans[transition == trans]
    sim_value <- sim_ntrans[transition == trans]
    
    top_sim_transitions_fraction_in_obs <- top_sim_transitions_fraction_in_obs + obs_value / sum(obs_ntrans)
    top_sim_transitions_fraction_in_sim <- top_sim_transitions_fraction_in_sim + sim_value / sum(sim_ntrans)
    
    obs_rank <- rank(obs_ntrans)[transition == trans]
    obs_rank_of_top_sim_transitions <- c(obs_rank_of_top_sim_transitions, obs_rank)
    
    sim_rank <- rank(sim_ntrans)[transition == trans]
    sim_rank_of_top_sim_transitions <- c(sim_rank_of_top_sim_transitions, sim_rank)
    
  }
  lineage_vector <- c(lineage_vector, lineage)
  mean_change_of_top_sim_transitions_vector <- c(mean_change_of_top_sim_transitions_vector, mean_change_of_top_sim_transitions)
  mean_change_of_top_obs_transitions_vector <- c(mean_change_of_top_obs_transitions_vector, mean_change_of_top_obs_transitions)
  top_sim_transitions_fraction_in_sim_vector <- c(top_sim_transitions_fraction_in_sim_vector, top_sim_transitions_fraction_in_sim)
  top_sim_transitions_fraction_in_obs_vector <- c(top_sim_transitions_fraction_in_obs_vector, top_sim_transitions_fraction_in_obs)
  
}
summary_dataframe <- data.frame(lineage = lineage_vector,
                                mean_change_of_top_sim_transitions = mean_change_of_top_sim_transitions_vector,
                                mean_change_of_top_obs_transitions = mean_change_of_top_obs_transitions_vector,
                                top_sim_transitions_fraction_in_sim = top_sim_transitions_fraction_in_sim_vector,
                                top_sim_transitions_fraction_in_obs = top_sim_transitions_fraction_in_obs_vector)


mean(summary_dataframe[c(1,2,4,5,7),]$mean_change_of_top_sim_transitions / summary_dataframe[c(1,2,4,5,7),]$mean_change_of_top_obs_transitions)
# sim_vs_obs_change <- function(lineage,scat_dataframe = scatterplot_dataframe){
#   lineage_dataframe = scat_dataframe[scat_dataframe$lineage == lineage &
#                                        scat_dataframe$model == 'S5F',]
#   
#   pl <- ggplot(data = lineage_dataframe, aes(x = obs_total_change_in_mean_log_S5F, 
#                                              y = sim_total_change_in_mean_log_S5F)) + 
#     #geom_point() +
#     ggplot_theme +
#     theme(plot.margin = margin(30,10,1,1,'pt'))+
#     geom_abline(intercept = 0, slope = 1, linetype = 2) +
#     xlab('Total mutability change in observed tree') +
#     ylab('Total mutability change simulated under S5F model') +
#     geom_text(aes(label = transition), size = 2)
#   
#   return(pl)
# }
# obs_vs_sim_change_plots <- sapply(c('CH103','CH103L','VRC26','VRC26L','VRC01_13','VRC01_01','VRC01_19'),
#                                   FUN = sim_vs_obs_change, simplify = FALSE)
# 
# sim_vs_obs_change_figure <- plot_grid(plotlist = obs_vs_sim_change_plots, nrow = 3, ncol = 3, 
#                                       labels = names(obs_vs_sim_change_plots))
# 
# pdf('AA_transition_figs/sim_vs_obs_change.pdf', height = 12, width = 12)
# plot(sim_vs_obs_change_figure)
# dev.off()



# delta_change_vs_obs_ntrans <- function(lineage,scat_dataframe = scatterplot_dataframe){
#   lineage_dataframe = scat_dataframe[scat_dataframe$lineage == lineage &
#                                        scat_dataframe$model == 'S5F',]
#   
#   pl <- ggplot(data = lineage_dataframe, aes(x = obs_ntrans, 
#                                              y = sim_total_change_in_mean_log_S5F - obs_total_change_in_mean_log_S5F)) + 
#     #geom_point() +
#     ggplot_theme +
#     theme(plot.margin = margin(30,10,1,1,'pt'))+
#     xlab('Observed number of transitions') +
#     ylab('Simulated - observed mutability change') +
#     geom_text(aes(label = transition), size = 2)
#   return(pl)
# }
# 
# delta_change_vs_obs_ntrans_plots <- sapply(c('CH103','CH103L','VRC26','VRC26L','VRC01_13','VRC01_01','VRC01_19'),
#                                   FUN = delta_change_vs_obs_ntrans, simplify = FALSE)
# 
# delta_change_vs_obs_ntrans_figure <- plot_grid(plotlist = delta_change_vs_obs_ntrans_plots, nrow = 3, ncol = 3, 
#                                       labels = names(delta_change_vs_obs_ntrans_plots))
# 
# pdf('AA_transition_figs/delta_change_vs_obs_ntrans.pdf', height = 12, width = 12)
# plot(delta_change_vs_obs_ntrans_figure)
# dev.off()
# 
# average_change_vs_delta_ntrans <- function(lineage,scat_dataframe = scatterplot_dataframe){
#   lineage_dataframe = scat_dataframe[scat_dataframe$lineage == lineage &
#                                        scat_dataframe$model == 'S5F',]
#   
#   pl <- ggplot(data = lineage_dataframe, aes(x = sim_ntrans - obs_ntrans, 
#                                              y = sim_total_change_in_mean_log_S5F/sim_ntrans)) + 
#     #geom_point() +
#     geom_hline(yintercept = 0, linetype = 2) +
#     ggplot_theme +
#     theme(plot.margin = margin(30,10,1,1,'pt'))+
#     xlab('Simulated - observed number of transitions') +
#     ylab('Average simulated change in mean-log S5F') +
#     geom_text(aes(label = transition), size = 2) +
#     ylim(-0.03,0.03)
#   return(pl)
# }
# 
# average_change_vs_delta_ntrans_plots <- sapply(c('CH103','CH103L','VRC26','VRC26L','VRC01_13','VRC01_01','VRC01_19'),
#                                            FUN = average_change_vs_delta_ntrans, simplify = FALSE)
# 
# average_change_vs_delta_ntrans_figure <- plot_grid(plotlist = average_change_vs_delta_ntrans_plots, nrow = 3, ncol = 3, 
#                                                labels = names(average_change_vs_delta_ntrans_plots))
# 
# pdf('AA_transition_figs/average_change_vs_delta_ntrans.pdf', height = 12, width = 12)
# plot(average_change_vs_delta_ntrans_figure)
# dev.off()





