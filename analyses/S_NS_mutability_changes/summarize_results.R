# Summarizes results of the analyses and simulations on MCC trees
# Uses the following MCMC chains: CH103_con_run1a, CH103L_con_run1a, VRC26int_con_run1a, VRC26L_con_run1a, VRC01_01_log_run1a, VRC01_13_log_run1a, VRC01_19_log_run1a,

library('ggplot2')
library('reshape')
library('gridExtra')
library('grid')
library('coda')
library('lattice')
library('gridBase')
options(expressions=10000)

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

# ======= GET DATAFRAMES WITH SIMULATION RESULTS FOR ALL LINEAGES
CH103_dataframe_sim <- read.table(paste(results_directory, 'CH103_constant/CH103_con_run1a_simulations_MCC.csv', sep = ''), header = T, sep = ',')
CH103L_dataframe_sim <- read.table(paste(results_directory, 'CH103L_constant/CH103L_con_run1a_simulations_MCC.csv', sep = ''), header = T, sep = ',')
VRC26int_dataframe_sim <- read.table(paste(results_directory, 'VRC26int_constant/VRC26int_con_run1a_simulations_MCC.csv', sep = ''), header = T, sep = ',')
VRC26L_dataframe_sim <- read.table(paste(results_directory, 'VRC26L_constant/VRC26L_con_run1a_simulations_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_01_dataframe_sim <- read.table(paste(results_directory, 'VRC01_01_logistic/VRC01_01_log_run1a_simulations_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_13_dataframe_sim <- read.table(paste(results_directory, 'VRC01_13_logistic/VRC01_13_log_run1a_simulations_MCC.csv', sep = ''), header = T, sep = ',')
VRC01_19_dataframe_sim <- read.table(paste(results_directory, 'VRC01_19_logistic/VRC01_19_log_run1a_simulations_MCC.csv', sep = ''), header = T, sep = ',')

dataframe_sim_list <- list('CH103' = CH103_dataframe_sim, 'CH103L' = CH103L_dataframe_sim,
                       'VRC26' = VRC26int_dataframe_sim, 'VRC26L' = VRC26L_dataframe_sim,
                       'VRC01_01' = VRC01_01_dataframe_sim, 'VRC01_13' = VRC01_13_dataframe_sim,
                       'VRC01_19' = VRC01_19_dataframe_sim)

# ======= PARTITION POINTS DELIMITING FRS and CDRS FOR EACH LINEAGE:
# c(First position FR1, 1st position CDR1, 1st position FR2, 1st position CDR2, 1st position FR3, 1st position CDR3, LAST POSITION CDR3)
partition_points <- list('CH103' = c(1, 34, 58, 109, 130, 244, 288),
                         'CH103L' = c(1, 76, 94, 145, 154, 262, 285),
                         'VRC01_01' = c(1, 76, 100, 151, 175, 289, 363),
                         'VRC01_13' = c(1, 76, 100, 151, 175, 289, 342),
                         'VRC01_19' = c(1, 76, 100, 151, 175, 310,429),
                         'VRC26' = c(1, 82, 106, 157, 181, 295, 453),
                         'VRC26L' = c(1, 76, 100, 151, 160, 268, 294)
                          )


# ====== COMBINED DATAFRAME WITH SIMULATION RESULTS FOR ALL LINEAGES =======

time_from_UCA <- c()
lineage <- c()
metric_vector <- c()
region_vector <- c()
contrast_true <- c()
contrast_simulation_S5F_S5F <- c()
contrast_simulation_uniform_S5F <- c()
contrast_simulation_CP_S5F <- c()
n_NT_diffs <- c()
n_syn_diffs <- c()
n_nonsyn_diffs <- c()
branch_exp_subs <- c()

for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
#for(clone in c('VRC01_13')){
  print(clone)
  clone_dataframe <- dataframe_list[[clone]] 

  for(metric in c('S5F','HS','OHS')){
    #for(metric in c('S5F')){
    for(region in c('WS','FR','CDR')){
      for(row in 1:nrow(clone_dataframe)){
        
        metric_vector <- c(metric_vector, metric)
        region_vector <- c(region_vector, region)
        lineage <- c(lineage, clone)
        
        time_from_UCA <- c(time_from_UCA, clone_dataframe[row, 'parent_time_to_root'])
        
        n_NT_diffs <- c(n_NT_diffs, clone_dataframe[row, 'n_NT_diffs'])
        n_syn_diffs <- c(n_syn_diffs, clone_dataframe[row, 'n_Syn_diffs'])
        n_nonsyn_diffs <- c(n_nonsyn_diffs, clone_dataframe[row, 'n_NonSyn_diffs'])
        branch_exp_subs <- c(branch_exp_subs, clone_dataframe[row, 'branch_exp_subs'])
        
        parent_mutability <- clone_dataframe[row,paste(metric,'_',region, '_parent', sep ='')]
        
        child_mutability_true <- clone_dataframe[row,paste(metric,'_',region, '_child_observed', sep ='')]
        child_mutability_simulation_S5F_S5F <-  clone_dataframe[row,paste(metric,'_',region, '_child_sim_S5FMut_S5FTrans', sep ='')]
        child_mutability_simulation_uniform_S5F <-  clone_dataframe[row,paste(metric,'_',region, '_child_sim_uniformMut_S5FTrans', sep ='')]
        child_mutability_simulation_CP_S5F <- clone_dataframe[row,paste(metric,'_',region, '_child_sim_CPMut_S5FTrans', sep ='')]
        
        contrast_true <- c(contrast_true, child_mutability_true - parent_mutability)
        contrast_simulation_S5F_S5F <- c(contrast_simulation_S5F_S5F, child_mutability_simulation_S5F_S5F - parent_mutability)
        contrast_simulation_uniform_S5F <- c(contrast_simulation_uniform_S5F, child_mutability_simulation_uniform_S5F - parent_mutability)
        contrast_simulation_CP_S5F <- c(contrast_simulation_CP_S5F, child_mutability_simulation_CP_S5F - parent_mutability)
        
      }
    }
  }
}

lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
metric_vector <- factor(metric_vector, levels = c('S5F','HS','OHS'))
region_vector <- factor(region_vector, levels = c('WS','FR','CDR'))  


combined_dataframe <- data.frame(lineage, metric = metric_vector, region = region_vector, time_from_UCA,
                                 contrast_true, contrast_simulation_S5F_S5F,
                                 contrast_simulation_uniform_S5F,contrast_simulation_CP_S5F,
                                 n_NT_diffs,n_syn_diffs,n_nonsyn_diffs,branch_exp_subs)

# ===== DATAFRAME WITH T-TEST RESULTS COMPARING EACH SIMULATION MODEL TO OBSERVED CONTRASTS =======

lineage <- c()
metric_vector <- c()
region_vector <- c()
type_vector <- c()

t_statistic <- c()
mean_absolute_error <- c()
p_value <-
fraction_of_change_explained <- c()
deg_freedom <- c()

mean_contrast <- c()

for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
  for(metric in c('S5F')){
    for(region in c('WS')){
      for(type in c('true','simulation_S5F_S5F','simulation_uniform_S5F','simulation_CP_S5F')){
        metric_vector <- c(metric_vector, metric)
        region_vector <- c(region_vector, region)
        type_vector <- c(type_vector, type)
        lineage <- c(lineage, clone)
        
        contrasts <- combined_dataframe[combined_dataframe$metric == metric & combined_dataframe$region == region &
                                          combined_dataframe$lineage == clone, paste('contrast_',type,sep='')]
        
        # To compute correlation...
        true_contrasts <- combined_dataframe[combined_dataframe$metric == metric & combined_dataframe$region == region &
                                               combined_dataframe$lineage == clone, 'contrast_true']
        
        mean_contrast <- c(mean_contrast, mean(contrasts))
        
        if(type == 'true'){
          t_statistic <- c(t_statistic, NA)
          p_value <- c(p_value, NA)
          fraction_of_change_explained <- c(fraction_of_change_explained, NA)
          mean_absolute_error <- c(mean_absolute_error, NA)
          deg_freedom <- c(deg_freedom, NA)
          
        }else{
          
          t_test <- t.test(contrasts - true_contrasts)
          
          t_statistic <- c(t_statistic, t_test$statistic)
          p_value <- c(p_value, t_test$p.value)
          # For some reason df is stored under 'parameter'
          deg_freedom <- c(deg_freedom, t_test$parameter)
          
          mean_absolute_error <- c(mean_absolute_error, 
                                   mean(abs(contrasts-true_contrasts)))
          
          # If mean contrast under model is in the same direction as observed mean contrast:
          fraction_of_change_explained <- c(fraction_of_change_explained, mean(contrasts)/mean(true_contrasts))

        }
      }
    }
  }
}

lineage <- factor(lineage, levels = c('CH103','CH103L','VRC26','VRC26L','VRC01_13',
                                      'VRC01_01','VRC01_19'))
metric_vector <- factor(metric_vector, levels = c('S5F','HS','OHS'))
region_vector <- factor(region_vector, levels = c('WS','FR','CDR'))  
type_vector <- factor(type_vector, levels = c('true','simulation_S5F_S5F','simulation_uniform_S5F','simulation_CP_S5F'))

ttest_dataframe <- data.frame(lineage, metric = metric_vector, region = region_vector, type = type_vector,
                              mean_contrast, mean_absolute_error, t_statistic, p_value, deg_freedom, fraction_of_change_explained)

names(ttest_dataframe)[names(ttest_dataframe) == 'type'] <- 'model'

write.table(ttest_dataframe, '../../results/simulations_MCC/simulation_MCC_summary.csv', 
            sep = ',', row.names = F)

# ===== DATAFRAME WITH CUMULATIVE CHANGES BY REGION AND SUBSTITUTION TYPE (SYN. / NON-SYN.) =======
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


# SUMMARY TABLE ON CHANGES ATTRIBUTABLE TO SELECTION
lineage_vector <- c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')
CDR_change_total <- c()
CDR_change_syn <- c()
CDR_change_nonsyn <- c()
for(clone in lineage_vector){
  p_points <- partition_points[[clone]]

  CDR_change_total = c(CDR_change_total, sum(subset(total_change_dataframe, lineage == clone & region == 'CDR')$cumulative_change))
  CDR_change_syn = c(CDR_change_syn, subset(total_change_dataframe, lineage == clone & region == 'CDR' & substitution_class == 'syn_only')$cumulative_change)
  CDR_change_nonsyn = c(CDR_change_nonsyn, subset(total_change_dataframe, lineage == clone & region == 'CDR' & substitution_class == 'nonsyn_only')$cumulative_change)

}

selection_dataframe <- data.frame(lineage = lineage_vector, CDR_change_total, CDR_change_nonsyn, CDR_change_syn)
write.table(selection_dataframe, '../../results/simulations_MCC/CDR_changes_summary.csv', 
            sep = ',', row.names = F)



# SUMMARY TABLE WITH SYN. AND NON-SYN CONTRIBUTIONS BY REGION:
# lineage_vector <- c()
# fraction_FR <- c()
# fraction_CDR <- c()
# WS_change_total <- c()
# WS_change_syn <- c()
# WS_change_nonsyn <- c()
# FR_change_total <- c()
# FR_change_syn <- c()
# FR_change_nonsyn <- c()
# CDR_change_total <- c()
# CDR_change_syn <- c()
# CDR_change_nonsyn <- c()
# 
# for(clone in c('CH103','CH103L','VRC26','VRC26L','VRC01_01','VRC01_13','VRC01_19')){
#   p_points <- partition_points[[clone]]
#   
#   length_FR <- length(c(p_points[1]:(p_points[2] -1), p_points[3]:(p_points[4] -1), p_points[5]:(p_points[6] -1)))
#   length_CDR <- length(c(p_points[2]:(p_points[3] -1), p_points[4]:(p_points[5] -1), p_points[6]:p_points[7]))
#   
#   fraction_FR <- c(fraction_FR, length_FR / (length_FR + length_CDR))
#   fraction_CDR <- c(fraction_CDR, length_CDR / (length_FR + length_CDR))
#   
#   WS_change_total = c(WS_change_total, sum(subset(total_change_dataframe, lineage == clone & region == 'WS')$cumulative_change))
#   WS_change_syn = c(WS_change_syn, subset(total_change_dataframe, lineage == clone & region == 'WS' & substitution_class == 'syn_only')$cumulative_change)
#   WS_change_nonsyn = c(WS_change_nonsyn, subset(total_change_dataframe, lineage == clone & region == 'WS' & substitution_class == 'nonsyn_only')$cumulative_change)
#   
#   FR_change_total = c(FR_change_total, sum(subset(total_change_dataframe, lineage == clone & region == 'FR')$cumulative_change))
#   FR_change_syn = c(FR_change_syn, subset(total_change_dataframe, lineage == clone & region == 'FR' & substitution_class == 'syn_only')$cumulative_change)
#   FR_change_nonsyn = c(FR_change_nonsyn, subset(total_change_dataframe, lineage == clone & region == 'FR' & substitution_class == 'nonsyn_only')$cumulative_change)
#   
#   CDR_change_total = c(CDR_change_total, sum(subset(total_change_dataframe, lineage == clone & region == 'CDR')$cumulative_change))
#   CDR_change_syn = c(CDR_change_syn, subset(total_change_dataframe, lineage == clone & region == 'CDR' & substitution_class == 'syn_only')$cumulative_change)
#   CDR_change_nonsyn = c(CDR_change_nonsyn, subset(total_change_dataframe, lineage == clone & region == 'CDR' & substitution_class == 'nonsyn_only')$cumulative_change)
#   
# }
