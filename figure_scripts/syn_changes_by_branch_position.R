# Produces figure for the manuscript showing the distribution of the fraction of negative changes in S5F, HS and OHS in observed lineages in the whole sequence
# Main text : S5F, HS and OHS
# Uses the following MCMC chains: CH103_con_run1a, CH103L_con_run1a, VRC26int_con_run1a, VRC26L_con_run1a, VRC01_01_log_run1a, VRC01_13_log_run1a, VRC01_19_log_run1a,

library('ggplot2')
library('grid')
library('coda')
library('lattice')

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
# USING ONLY S5F CHANGES DUE TO SYN. SUBSTITUTIONS.
fraction_negative_all <- c()
fraction_negative_terminal <- c()
fraction_negative_internal <- c()
fraction_negative_trunk <- c()
fraction_negative_nontrunk <- c()

lineage <- c()
region_column <- c()
#branch_is_terminal <- c()
  
for(clone in c('CH103','CH103L','VRC26int','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
#for(clone in c('CH103','CH103L','VRC26L','VRC01_01','VRC01_13','VRC01_19')){  
  clone_dataframe <- dataframe_list[[clone]]
  print(clone)
  #For each tree from the posterior distribution for the clone:
  for(tree in unique(dataframe_list[[clone]]$tree)[1:1000]){
      # For each region
      for(region in c('WS','FR','CDR')){
      
        #print(paste(clone,metric,region,tree))
        # Find contrasts for the specified metric and region
        contrasts_all <- clone_dataframe[clone_dataframe$tree == tree,paste('S5F_',region,'_contrast_syn',sep='')]
        contrasts_terminal <- clone_dataframe[clone_dataframe$tree == tree & clone_dataframe$branch_is_terminal == T,
                                              paste('S5F_',region,'_contrast_syn',sep='')]
        
        # Contrast for all internal branches (both trunk and non-trunk)
        contrasts_internal <- clone_dataframe[clone_dataframe$tree == tree & clone_dataframe$branch_is_terminal == F,
                                              paste('S5F_',region,'_contrast_syn',sep='')]
        contrasts_trunk <- clone_dataframe[clone_dataframe$tree == tree & clone_dataframe$branch_in_trunk == T,
                                           paste('S5F_',region,'_contrast_syn',sep='')]
        
        # Contrasts for branches that are internal but not in the trunk
        contrasts_nontrunk <- clone_dataframe[clone_dataframe$tree == tree & clone_dataframe$branch_in_trunk == F &
                                                clone_dataframe$branch_is_terminal == F,
                                              paste('S5F_',region,'_contrast_syn',sep='')]
    
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
        region_column <- c(region_column, region)
        }
    }
}

lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26int','VRC26L','VRC01_13',
                                           'VRC01_01','VRC01_19'))
region_column <- factor(region_column, levels = c('WS','FR','CDR'))

  
global_dataframe <- data.frame(fraction_negative_all = fraction_negative_all, fraction_negative_terminal = fraction_negative_terminal,
                               fraction_negative_internal = fraction_negative_internal, fraction_negative_trunk = fraction_negative_trunk,
                               fraction_negative_nontrunk = fraction_negative_nontrunk, lineage = lineage,
                               region = region_column)

# ==== MAKE DATAFRAME WITH MEAN FRACTION NEG. AND HPD LIMITS FOR EACH LINEAGE, METRIC AND REGION

mean_fraction_negative_all <- c()
mean_fraction_negative_terminal <- c()
mean_fraction_negative_internal <- c()
mean_fraction_negative_trunk <- c()
mean_fraction_negative_nontrunk <- c()

lineage <- c()
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
    for(region in c('WS','FR','CDR')){
      values_all <- global_dataframe$fraction_negative_all[global_dataframe$lineage == clone &
                                                     global_dataframe$region == region]
      
      values_terminal <- global_dataframe$fraction_negative_terminal[global_dataframe$lineage == clone &
                                                             global_dataframe$region == region]
      
      values_internal <- global_dataframe$fraction_negative_internal[global_dataframe$lineage == clone &
                                                                       global_dataframe$region == region]
      
      values_trunk <- global_dataframe$fraction_negative_trunk[global_dataframe$lineage == clone &
                                                                    global_dataframe$region == region]
      
      values_nontrunk <- global_dataframe$fraction_negative_nontrunk[global_dataframe$lineage == clone &
                                                                 global_dataframe$region == region]
      
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
      
      region_column <- c(region_column, region)
      lineage <- c(lineage, clone)
    }
}

lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26int','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
region_column <- factor(region_column, levels = c('WS','FR','CDR'))  

  
summary_dataframe <- data.frame(lineage,region=region_column,
                                mean_fraction_negative_all,HPD_llim_all,HPD_ulim_all,
                                mean_fraction_negative_terminal,HPD_llim_terminal,HPD_ulim_terminal,
                                mean_fraction_negative_internal,HPD_llim_internal,HPD_ulim_internal,
                                mean_fraction_negative_trunk,HPD_llim_trunk,HPD_ulim_trunk,
                                mean_fraction_negative_nontrunk,HPD_llim_nontrunk,HPD_ulim_nontrunk)


# =============== PLOTS =================
region_names <- c(
  'WS'="Whole sequence",
  'FR'="FRs",
  'CDR'="CDRs"
)

  # Base Plot
  pl <- ggplot(data = subset(summary_dataframe, region != 'WS'), 
               aes(x=lineage,y=mean_fraction_negative_internal)) + 
    scale_y_continuous(expand = c(0,0), limits = c(-0.05,1.05)) + 
    theme_bw() +
    #theme_classic() + 
    ylab('Frequency of mutability losses') + 
    xlab('Lineage') +
    ggplot_theme +
    theme(axis.text.x = element_text(angle = 90,vjust=0.6, size = 9),
          legend.position = 'top',
          legend.text=element_text(size=11),
          axis.title.x = element_text(size = 13),
          strip.text.x = element_text(size = 13)
    ) +
    
    # Add violin plots
    #geom_violin(size = 0.6, colour = 'gray40') +
    #geom_boxplot(colour = 'gray40') +

    # Add horizontal line at 0.5
    geom_hline(yintercept=0.5,linetype=5, size = 0.5) + 

    facet_grid(.~region, labeller = labeller(region = region_names)) + 
    theme(panel.margin.y = unit(1.3, 'lines'))
  
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
  
# ======= Plot with separate results for terminal and internal branches
    
    pl_terminal_vs_internal <- pl + 
      geom_linerange(data=subset(summary_dataframe, region != 'WS'),
                     aes(ymin=HPD_llim_terminal, ymax = HPD_ulim_terminal), 
                     colour = 'blue') +
      
      geom_point(data=subset(summary_dataframe, region != 'WS'),
                 aes(y=mean_fraction_negative_terminal),
                 shape = 21, fill = 'royalblue1', alpha = 0.7,size = 2.5)  +

      geom_linerange(data=subset(summary_dataframe, region != 'WS'),
                     aes(ymin=HPD_llim_internal, ymax = HPD_ulim_internal), 
                     colour = 'grey40', position = position_nudge(x=-0.3,y=0)) +
      
      geom_point(data=subset(summary_dataframe, region != 'WS'),
                 aes(y=mean_fraction_negative_internal),
                 shape = 21, fill = 'grey25', alpha = 0.7,size = 2.5,
                 position = position_nudge(x=-0.3,y=0)) +
      
    # Artificial points plotted outside of plotting region, only purpose is to quickly produce ggplot legend
    geom_pointrange(data = data.frame(lineage = c('CH103','CH103','CH103L','CH103L'),
                                 change = c(-0.1,-0.1,-0.1,-0.1),
                                 type = factor(c('Internal branches','Terminal branches',
                                                 'Internal branches','Terminal branches'),
                                               levels = c('Internal branches','Terminal branches'))
    ),
    aes(x = lineage, y = change, ymin = change, ymax = change, colour = type, fill = type), 
        shape = 21, size = 1, alpha = 0.7) +
      scale_colour_manual(values = c('grey25','royalblue1'),
                        labels = c('Internal\nbranches','Terminal\nbranches')) +
      scale_fill_manual(values = c('grey25','royalblue1'),
                        labels = c('Internal\nbranches','Terminal\nbranches')) +
      
      guides(colour=guide_legend(title=NULL)) + 
      guides(fill=guide_legend(title=NULL)) + 
      
      scale_x_discrete(labels=x_labels)
  
####
pdf('syn_S5F_changes_terminal_vs_internal.pdf',width = 5.9, height = 4.5)
plot(pl_terminal_vs_internal)
dev.off()
####
    
    
    
    
    # Plot results with trunk vs tips vs non-trunk internal
    
    pl_trunk_vs_nontrunk <- pl + 
      
      geom_linerange(data=subset(summary_dataframe, region != 'WS'),
                     aes(ymin=HPD_llim_trunk, ymax = HPD_ulim_trunk), 
                     colour = 'darkorange2', position = position_nudge(x=0.2,y=0)) +
      
      geom_point(data=subset(summary_dataframe, region != 'WS'),
                 aes(y=mean_fraction_negative_trunk),
                 shape = 21, fill = 'darkorange3', alpha = 0.7,size = 2.5, position = position_nudge(x=0.2,y=0)) +
      
      
      geom_linerange(data=subset(summary_dataframe, region != 'WS'),
                     aes(ymin=HPD_llim_terminal, ymax = HPD_ulim_terminal), 
                     colour = 'blue') +
      
      geom_point(data=subset(summary_dataframe, region != 'WS'),
                 aes(y=mean_fraction_negative_terminal),
                 shape = 21, fill = 'royalblue1', alpha = 0.7,size = 2.5)  +
      
      geom_linerange(data=subset(summary_dataframe, region != 'WS'),
                     aes(ymin=HPD_llim_nontrunk, ymax = HPD_ulim_nontrunk), 
                     colour = 'grey40', position = position_nudge(x=-0.2,y=0)) +
      
      geom_point(data=subset(summary_dataframe, region != 'WS'),
                 aes(y=mean_fraction_negative_nontrunk),
                 shape = 21, fill = 'grey25', alpha = 0.7,size = 2.5,
                 position = position_nudge(x=-0.2,y=0)) +

      # Artificial points plotted outside of plotting region, only purpose is to quickly produce ggplot legend
      geom_pointrange(data = data.frame(lineage = c('CH103','CH103','CH103','CH103L','CH103L','CH103L'),
                                        change = c(-0.1,-0.1,-0.1,-0.1,-0.1,-0.1),
                                        type = factor(c('Internal branches','trunk','Terminal branches',
                                                        'Internal branches','trunk','Terminal branches'),
                                                      levels = c('Internal branches','trunk','Terminal branches'))
      ),
      aes(x = lineage, y = change, ymin = change, ymax = change, colour = type, fill = type), 
      shape = 21, size = 1, alpha = 0.7) +
      scale_colour_manual(values = c('royalblue1','darkorange2','grey25'),
                          labels = c('Terminal\nbranches','Trunk','Non-trunk\ninternal branches')) +
      scale_fill_manual(values = c('royalblue1','darkorange2','grey25'),
                        labels =  c('Terminal\nbranches','Trunk','Non-trunk\ninternal branches')) +
      guides(colour=guide_legend(title=NULL)) + 
      guides(fill=guide_legend(title=NULL)) +
      
      scale_x_discrete(labels=x_labels)
    
  

  ####
  pdf('syn_S5F_changes_trunk_vs_nontrunk.pdf', width = 5.9, height = 4.5)
  plot(pl_trunk_vs_nontrunk)
  dev.off()
  ####
  

