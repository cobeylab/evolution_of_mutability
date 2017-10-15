# Produces figure for the manuscript showing the mutability of the inferred ancestral sequence for chain "con_run1a" of each lineage, together with the distribution of mutability values we get by randomizing that sequence while keeping the amino acid sequence constant".
library('ggplot2')
library('reshape')
library('gridExtra')
library('grid')
library('coda')
library('lattice')
library('cowplot')

results_directory <- '../results/relative_mutability/observed_lineages/'

# ======= GET DATAFRAMES WITH OBSERVED MUTABILITY FOR ALL LINEAGES

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
  
  
  #  ======= GET DATAFRAMES WITH RANDOMIZED MUTABILITY FOR ALL LINEAGES
CH103_randomized_dataframe <- read.table(paste(results_directory, 'CH103_constant/CH103_con_run1a_randomized_mutability_MCC.csv', sep = ''), header = T, sep = ',')
CH103L_randomized_dataframe <- read.table(paste(results_directory, 'CH103L_constant/CH103L_con_run1a_randomized_mutability_MCC.csv', sep = ''), header = T, sep = ',')
VRC26int_randomized_dataframe <- read.table(paste(results_directory, 'VRC26int_constant/VRC26int_con_run1a_randomized_mutability_MCC.csv', sep = ''), header = T, sep = ',')
VRC26L_randomized_dataframe <- read.table(paste(results_directory, 'VRC26L_constant/VRC26L_con_run1a_randomized_mutability_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_01_randomized_dataframe <- read.table(paste(results_directory, 'VRC01_01_logistic/VRC01_01_log_run1a_randomized_mutability_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_13_randomized_dataframe <- read.table(paste(results_directory, 'VRC01_13_logistic/VRC01_13_log_run1a_randomized_mutability_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_19_randomized_dataframe <- read.table(paste(results_directory, 'VRC01_19_logistic/VRC01_19_log_run1a_randomized_mutability_MCC.csv', sep = ''), header = T, sep = ',')
  
  
randomized_dataframe_list <- list('CH103' = CH103_randomized_dataframe, 'CH103L' = CH103L_randomized_dataframe,
                                 'VRC26int' = VRC26int_randomized_dataframe, 'VRC26L' = VRC26L_randomized_dataframe,
                                 'VRC01_01' = VRC01_01_randomized_dataframe, 'VRC01_13' = VRC01_13_randomized_dataframe,
                                 'VRC01_19' = VRC01_19_randomized_dataframe)
  

# ================ IMPORTING GGPLOT2 PARAMETERS ===============
source('ggplot_parameters.R')


# ==== MAKE GLOBAL DATAFRAME WITH RANDOMIZED ANCESTRAL VALUES FOR EACH LINEAGE, REPLICATE, REGION AND METRIC (FOR VIOLIN PLOTS)
lineage <- c()
metric_vector <- c()
region_vector <- c()
mutability <- c()

for(clone in c('CH103','CH103L','VRC26int','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  clone_dataframe <- randomized_dataframe_list[[clone]]
  print(clone)
  for(metric in c('S5F','HS','OHS','geomS5F')){
    for(region in c('WS','FR','CDR')){
        lineage <- c(lineage, rep(clone, sum(clone_dataframe$sequence_id == 'Node_0')))
        metric_vector <- c(metric_vector, rep(metric, sum(clone_dataframe$sequence_id == 'Node_0')))
        region_vector <- c(region_vector, rep(region, sum(clone_dataframe$sequence_id == 'Node_0')))

        mutability <- c(mutability,
                        clone_dataframe[clone_dataframe$sequence_id == 'Node_0',
                                        paste('randomized_',metric,'_',region,'_allsites',sep='')]) 
      }
   }
}
lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26int','VRC26L','VRC01_13',
                                           'VRC01_01','VRC01_19'))
metric_vector <- factor(metric_vector, levels = c('S5F','HS','OHS','geomS5F'))
region_vector <- factor(region_vector, levels = c('WS','FR','CDR'))  
  
global_randomized_dataframe <- data.frame(lineage=lineage, metric=metric_vector, region=region_vector, mutability=mutability)

# ==== MAKE GLOBAL DATAFRAME WITH OBSERVED VALUE FOR EACH LINEAGE, REPLICATE, REGION AND METRIC (FOR VIOLIN PLOTS)
lineage <- c()
metric_vector <- c()
region_vector <- c()
mutability <- c()
percentile <- c()

for(clone in c('CH103','CH103L','VRC26int','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  clone_dataframe <- observed_dataframe_list[[clone]]
    for(metric in c('S5F','HS','OHS','geomS5F')){
      for(region in c('WS','FR','CDR')){
        
        lineage <- c(lineage, clone)
        metric_vector <- c(metric_vector, metric)
        region_vector <- c(region_vector, region)
        
        obs_mutability <-  clone_dataframe[clone_dataframe$sequence_id == 'Node_0',
                                           paste('observed_',metric,'_',region,sep='')]
        
        mutability <- c(mutability, obs_mutability)
        
        randomized_values <- global_randomized_dataframe$mutability[global_randomized_dataframe$lineage == clone &
                                                                    global_randomized_dataframe$metric == metric &
                                                                    global_randomized_dataframe$region == region]
        
        percentile <- c(percentile, sum(obs_mutability >= randomized_values)/length(randomized_values))
        
        
      }
  }
}

lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26int','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
metric_vector <- factor(metric_vector, levels = c('S5F','HS','OHS','geomS5F'))
region_vector <- factor(region_vector, levels = c('WS','FR','CDR'))  

global_observed_dataframe <- data.frame(lineage=lineage, metric=metric_vector, region=region_vector, 
                                        mutability=mutability, percentile)


# Plot
  pl <- ggplot(data = subset(global_randomized_dataframe, metric == 'geomS5F'), 
               aes(x=lineage,y=mutability)) + 
    #scale_y_continuous(expand = c(0,0), limits = c(0,1)) + 
    theme_bw() +
    #theme_classic() + 
    ylab('Geometric mean of S5F mutability') + 
    xlab('Lineage') +
    theme(axis.title.y = element_text(size = axis_title_size,
                                      margin = margin(0,ylab_distance,0,0)),
          axis.title.x = element_text(size = axis_title_size,
                                      margin = margin(xlab_distance,0,0,0)),
          axis.text.x = element_text(size = x_axis_text_size,angle = 90,vjust=0.6),
          axis.text.y = element_text(size = y_axis_text_size),
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(0.5, "cm"),
          plot.margin = unit(c(1, 1, 0.2, 0.2),"cm"),
          title = element_text(size = title_size),
          strip.text = element_text(size = 20)
          #axis.line.x = element_line(colour="black", size = 1.5),
          #axis.line.y = element_line(colour="black", size = 1.5)
    ) +
    
    # Add violin plots
    geom_violin(size = 0.6, colour = 'gray40') +
    #geom_boxplot(colour = 'gray40') +

    facet_grid(.~region,scales='free_y') 
    #theme(panel.margin.y = unit(1.3, 'lines')) +
    

    # Add observed value for each lineage, metric and region
    
    pl <- pl + geom_point(data=subset(global_observed_dataframe, metric == 'geomS5F'),
                               aes(y=mutability),
                               shape = 22,
                               fill = 'firebrick1',
                               alpha =0.8,
                               size = 3)
    
    
  # Change x-axis labels
  x_labels <- c()
  for(level in levels(global_randomized_dataframe$lineage)){
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
  pl <- pl + scale_x_discrete(labels=x_labels)
  
  ####
  pdf('relative_ancestral_mutability.pdf',width=10,height=5)
  plot(pl)
  dev.off()
  ####

  # CDR FR difference in mean log S5F
  mean_logS5F_CDR <- log(subset(global_observed_dataframe, metric == 'geomS5F' & region == 'CDR')$mutability)
  mean_logS5F_FR <- log(subset(global_observed_dataframe, metric == 'geomS5F' & region == 'FR')$mutability)
  
  mean(mean_logS5F_CDR - mean_logS5F_FR)
  
  
  
  # CDR / FR ratios:
  #ratios <- subset(global_observed_dataframe, metric == 'S5F' & region == 'CDR')$mutability / subset(global_observed_dataframe, metric == 'S5F' & region == 'FR')$mutability -1
  
  #mean(ratios)
  #max(ratios)
  #min(ratios)
  
  percentiles_CDR <- subset(global_observed_dataframe,  metric == 'S5F' & region == 'CDR')$percentile
  percentiles_FR <- subset(global_observed_dataframe,  metric == 'S5F' & region == 'FR')$percentile
