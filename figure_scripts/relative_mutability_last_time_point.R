library('ggplot2')
library('reshape')
library('gridExtra')
library('grid')
library('coda')
library('lattice')
library('gridBase')

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
  

# ==== MAKE GLOBAL DATAFRAME WITH FR AND CDR geomS5F-MUTABILITY PERCENTILES FOR SEQUENCES AT THE LAST TIME POINT, FOR ALL LINEAGES ====
lineage <- c()
region_vector <- c()
max_time <- c()
percentile <- c()

  
for(clone in c('CH103','CH103L','VRC26int','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  clone_dataframe_obs <- observed_dataframe_list[[clone]]
  clone_dataframe_rdm <- randomized_dataframe_list[[clone]]
  
  mtime <- max(clone_dataframe_obs$time_from_root)
  
  clone_dataframe_obs <- subset(clone_dataframe_obs, time_from_root == mtime)
  
  #percents_FR <- c()
  #percents_CDR <- c()
  
  for(node in clone_dataframe_obs$sequence_id){
    for(region in c('FR','CDR')){
      
      observed <- clone_dataframe_obs[clone_dataframe_obs$sequence_id == node, paste('observed_geomS5F_',region,sep='')]
      randomized <- clone_dataframe_rdm[clone_dataframe_rdm$sequence_id == node, paste('randomized_geomS5F_',region, '_allsites',sep='')]
      
      
      lineage <- c(lineage, clone)
      region_vector <- c(region_vector, region)
      
      if(grepl('VRC01', clone)){
        max_time <- c(max_time, 4*mtime)
      }else{
        max_time <- c(max_time, mtime)
      }
      
      percentile <- c(percentile, sum(observed >= randomized) / length(randomized))
     
    }
  }
  
}
  
lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26int','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
region_vector <- factor(region_vector, levels = c('FR','CDR'))

combined_dataframe <- data.frame(lineage, max_time, region = region_vector, percentile)

# ==== MAKE GLOBAL DATAFRAME WITH FR AND CDR S5F-MUTABILITY PERCENTILES FOR THE ANCESTRAL SEQUENCE, FOR ALL LINEAGES ====
lineage <- c()
region_vector <- c()
ancestral_percentile <- c()


for(clone in c('CH103','CH103L','VRC26int','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  clone_dataframe_obs <- observed_dataframe_list[[clone]]
  clone_dataframe_rdm <- randomized_dataframe_list[[clone]]
  
  for(node in c('Node_0')){
    for(region in c('FR','CDR')){
      
      observed <- clone_dataframe_obs[clone_dataframe_obs$sequence_id == node, paste('observed_geomS5F_',region,sep='')]
      randomized <- clone_dataframe_rdm[clone_dataframe_rdm$sequence_id == node, paste('randomized_geomS5F_',region, '_allsites',sep='')]
      
      
      lineage <- c(lineage, clone)
      region_vector <- c(region_vector, region)
      
      if(grepl('VRC01', clone)){
        max_time <- c(max_time, 4*mtime)
      }else{
        max_time <- c(max_time, mtime)
      }
      
      ancestral_percentile <- c(ancestral_percentile, sum(observed >= randomized) / length(randomized))
      
    }
  }
  
}

lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26int','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
region_vector <- factor(region_vector, levels = c('FR','CDR'))

combined_dataframe_ancestor <- data.frame(lineage, region = region_vector, ancestral_percentile)




# PLOT
pl <- ggplot(data = combined_dataframe, 
               aes(x=lineage,y=percentile)) + 
    #scale_y_continuous(expand = c(0,0), limits = c(0,1)) + 
    theme_bw() +
    #theme_classic() + 
    ylab('Mutability percentile') + 
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
    
    # Add boxplot with randomized percentiles
    geom_boxplot(colour = 'gray40', alpha = 0.7) +

    # Add ancestral percentiles:
    geom_point(data = combined_dataframe_ancestor, 
               aes(x= lineage, y = ancestral_percentile), size = 2, colour = 'red') +
    
    facet_grid(.~region,scales='free_y') 
    #theme(panel.margin.y = unit(1.3, 'lines')) +
    
  # Change x-axis labels
  x_labels <- c()
  for(level in levels(combined_dataframe$lineage)){
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
pdf('relative_mutability_last_time_point.pdf',width=10,height=7)
plot(pl)
dev.off()
####

FR_percentiles <- c()
CDR_percentiles <- c()
for(clone in c('CH103','CH103L','VRC26int','VRC26L','VRC01_13','VRC01_01','VRC01_19')){
  FR_percent <- subset(combined_dataframe, region == 'FR' & lineage == clone)$percentile
  FR_percentiles <- c(FR_percentiles, mean(FR_percent))
    
  CDR_percent <- subset(combined_dataframe, region == 'CDR' & lineage == clone)$percentile
  CDR_percentiles <- c(CDR_percentiles, mean(CDR_percent))
}
  
mean(FR_percentiles)
max(FR_percentiles)
min(FR_percentiles)
  
mean(CDR_percentiles)
  max(CDR_percentiles)
  min(CDR_percentiles)
  