# Produces figure for the manuscript showing the results of the selection analysis using Baseline
# Uses the following MCMC chains: CH103_con_run1a, CH103L_con_run1a, VRC26int_con_run1a, VRC26L_con_run1a, VRC01_01_log_run1a, VRC01_13_log_run1a, VRC01_19_log_run1a,

library('ggplot2')
library('cowplot')
library('grid')
library('lattice')

results_dataframe <- read.csv('../results/Liao_ancestors_analysis/Liao_ancestors_S_NS_changes.csv',
                              header = T)

# Reordering factors to plot pairs roughly in order from root to tips
results_dataframe$pair <- factor(results_dataframe$pair, 
                                 levels = c('UCA_VH-I8_VH', 'I8_VH-I7_VH',
                                            'I7_VH-1AZCETI5_VH', 'I7_VH-1A102RI6_VH',
                                            'I8_VH-I4_VH', 'I4_VH-I3_VH',
                                            'I4_VH-1AH92U_VH', 'I3_VH-I2_VH',
                                            'I3_VH-CH105_VH', 'I2_VH-I1_VH',
                                            'I2_VH-CH104_VH', 'I1_VH-CH103_VH',
                                            'I1_VH-CH106_VH'
                                            ))
results_dataframe$substitution_class <- factor(results_dataframe$substitution_class,
                                               levels = c('total','syn','nonsyn'))

source('ggplot_parameters.R')



# =====================================================================================================================
# ================================================== PLOT ============================================================
# =====================================================================================================================

pl_CDR <- ggplot(subset(results_dataframe, region == 'CDR' & substitution_class != 
                      'total'), 
                 aes(x = pair, y = S5F_mutability_change)) +
  ggplot_theme +
  xlab('Ancestor-descendant pair') +
  ylab('Change in mean CDR mutability') +
  ylim(-0.17,0.05) +
  geom_hline(yintercept = 0, linetype= 2) +
  geom_col(aes(fill = factor(substitution_class)), position = position_stack(),
           width = 0.7) +

  geom_text(data = subset(results_dataframe, region == 'CDR' & substitution_class == 
                            'syn'), aes(x = pair, label=n_codon_changes, y = -0.16), colour = 'gray80') +
  
  geom_text(data = subset(results_dataframe, region == 'CDR' & substitution_class == 
                     'nonsyn'), aes(x = pair, label=n_codon_changes, y = -0.17), colour = 'gray40') +
  
  geom_text(data = data.frame(x=0.5,y=-0.15, text='Number of differences'), 
            aes(x = x, y = y,label = text), hjust = -0.1, size = 3.5) +
  theme(axis.text.x = element_text(size = 9),
        legend.position = 'top') +
  scale_x_discrete(labels = gsub('-','-\n',gsub('_VH','',levels(results_dataframe$pair))),
                   expand = c(0,0)) +
  scale_fill_manual(name = '', values = c('gray80','gray40'),
                    labels = c('Synonymous  ','Non-synonynous' )) +
  scale_colour_manual(name = '', values = c('gray80','gray40'),
                    labels = c('Synonymous  ','Non-synonynous' ))
  
pl_FR <- ggplot(subset(results_dataframe, region == 'FR' & substitution_class != 
                          'total'), 
                 aes(x = pair, y = S5F_mutability_change)) +
  ggplot_theme +
  xlab('Ancestor-descendant pair') +
  ylab('Change in mean FR mutability') +
  ylim(-0.09,0.025) +
  geom_hline(yintercept = 0, linetype= 2) +
  geom_col(aes(fill = factor(substitution_class)), position = position_stack(),
           width = 0.7) +
  
  geom_text(data = subset(results_dataframe, region == 'FR' & substitution_class == 
                            'syn'), aes(x = pair, label=n_codon_changes, y = -0.085), colour = 'gray80') +
  
  geom_text(data = subset(results_dataframe, region == 'FR' & substitution_class == 
                            'nonsyn'), aes(x = pair, label=n_codon_changes, y = -0.09), colour = 'gray40') +
  
  geom_text(data = data.frame(x=0.5,y=-0.08, text='Number of differences'), 
            aes(x = x, y = y,label = text), hjust = -0.1, size = 3.5) +
  theme(axis.text.x = element_text(size = 9),
        legend.position = 'top') +
  scale_x_discrete(labels = gsub('-','-\n',gsub('_VH','',levels(results_dataframe$pair)))) +
  scale_fill_manual(name = '', values = c('gray80','gray40'),
                    labels = c('Synonymous  ','Non-synonynous' )) +
  scale_colour_manual(name = '', values = c('gray80','gray40'),
                      labels = c('Synonymous  ','Non-synonynous' ))

 