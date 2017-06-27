library('ggplot2')
library('reshape')
library('gridExtra')
library('grid')
library('coda')
library('scales')
library('cowplot')

source('ggplot_parameters.R')

# Function for making ggplot show axis labels in exponent notation
# From http://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

# ===== DEFINING BASE PLOTTING FUNCTION FOR RATES VS TIME (POINTS) =====

base_plot_time <- function(results_dataframe, rate, time_units){
  
  # Sample of 100 trees
  sample_size <- 100
  tree_sample <- sample(unique(results_dataframe$tree), sample_size, replace = F)
  
  ylabel <- switch(rate,
                   total_rate = "Total substitution rate",
                   total_rate_RC = "Total substitution rate (RC)",
                   S_rate = "Synonymous substitution rate",
                   N_rate = "Non-synonymous rate")
  
  # Plot:
  
  pl <- ggplot() + 
    # Internal: black (grey25), Terminal: 'royalblue1'
    scale_colour_manual(values=c("grey25",'royalblue1')) +
    xlab(paste("Time since MRCA (", time_units, ')', sep = '')) + 
    ylab(ylabel) +
    ggplot_theme +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10),
          plot.margin = unit(c(1, 1, 1, 1),"cm")) +
    scale_y_continuous(label=scientific_10, expand = c(0.006,0.001)) +
    scale_x_continuous(expand = c(0.006,0.006) )
  
  for(tree in tree_sample){
    dataframe_subset <- results_dataframe[results_dataframe$tree == tree,]
    x <- dataframe_subset[,'parent_time']
    y <- dataframe_subset[,paste(rate,'_to_child',sep='')]
    factor <- factor(dataframe_subset[,'branch_is_terminal'])
    
    pl <- pl + geom_point(data=data.frame(x,y,factor), 
                          aes(x=x,y=y, colour = factor),
                          size = 4, alpha = 0.8)
    
    
  }
  
  return(pl)
}

plot_list_total_rate_RC <- c()
plot_list_S_rate <- c()
plot_list_N_rate <- c()


results_directory <- '../results/rates_vs_time/observed_lineages/'

clone_labels <- c('CH103 (H)', 'CH103 (L)', 'VRC26 (H)', 'VRC26 (L)', 
                  'VRC01-13 (H)', 'VRC01-01 (H)', 'VRC01-19 (H)')

time_units <- c('weeks','weeks','weeks','weeks','months','months','months')

paths <- c('CH103_constant/CH103_con_run1a', 'CH103L_constant/CH103L_con_run1a',
           'VRC26int_constant/VRC26int_con_run1a','VRC26L_constant/VRC26L_con_run1a',
           'VRC01_13_logistic/VRC01_13_log_run1a','VRC01_01_logistic/VRC01_01_log_run1a',
           'VRC01_19_logistic/VRC01_19_log_run1a'
           )
paths <- paste(results_directory,paths,'_rates_vs_time_points.csv', sep = '')


for(i in 1:length(clone_labels)){
  
  print(clone_labels[i])
  
  results_dataframe <- read.table(paths[i], header = T, sep = ',')
  
  plot_list_total_rate_RC[[clone_labels[i]]] <- base_plot_time(results_dataframe, 'total_rate_RC', time_units = time_units[i])
  plot_list_S_rate[[clone_labels[i]]] <- base_plot_time(results_dataframe, 'S_rate', time_units = time_units[i])
  plot_list_N_rate[[clone_labels[i]]] <- base_plot_time(results_dataframe, 'N_rate', time_units = time_units[i])
  
}

# Plot with results from uniform mutability simulations, w/ 30 (2c), 40 (2d) and 50% (2e) decrease in total rate


uniform_mutability_plot <- plot_grid(plot_list_total_rate_RC[['2c']],plot_list_S_rate[['2c']], plot_list_N_rate[['2c']], 
          plot_list_total_rate_RC[['2d']],plot_list_S_rate[['2d']], plot_list_N_rate[['2d']],
          plot_list_total_rate_RC[['2e']],plot_list_S_rate[['2e']], plot_list_N_rate[['2e']],
          labels = c('', '30% simulated decrease','','', '40% simulated decrease','','', '50% simulated decrease','')
          )


save_plot("rates_vs_time_uniform.pdf", uniform_mutability_plot,
          base_height = 10, base_width = 11
)


nonuniform_mutability_plot <- plot_grid(plot_list_total_rate_RC[['3a']],plot_list_S_rate[['3a']], plot_list_N_rate[['3a']], 
                                     plot_list_total_rate_RC[['3b']],plot_list_S_rate[['3b']], plot_list_N_rate[['3b']],
                                     plot_list_total_rate_RC[['3c']],plot_list_S_rate[['3c']], plot_list_N_rate[['3c']],
                                     labels = c('', 'S5F-based model','','', 'Hotspot-based model (3x)','','', 'Hotspot-based model (30x)','')
)

save_plot("rates_vs_time_nonuniform.pdf", nonuniform_mutability_plot,
          base_height = 10, base_width = 12
)





dataf <- data.frame(x,y,z)

plot_list <- list()

for(var in c('y','z')){
plot_list[[var]] <- ggplot(dataf, aes_string(x=x,y=var)) + geom_point() 
}




args <- commandArgs(trailingOnly = TRUE)

chain_id <- as.character(args[1])

options(expressions=15000)



# ========================= RETRIEVE DIRECTORIES FROM CHAIN ID =========================
# If chain is not from a simulated scenario
if(grepl("scenario", chain_id) == FALSE){
  time_units <- 'weeks'
  if(grepl('CH103', chain_id)){
    if(grepl('CH103L', chain_id)){
      clone <- 'CH103L'
    }else{
      clone <- 'CH103'
    }
  }else{
    if(grepl('VRC26', chain_id)){
      if(grepl('VRC26L', chain_id)){
        clone <- 'VRC26L'
      }else{
        clone <- 'VRC26'
      }
    }else{
    if(grepl('VRC01', chain_id)){
      if(grepl('VRC01_01', chain_id)){
        clone <- 'VRC01_01'}
      if(grepl('VRC01_13', chain_id)){
        clone <- 'VRC01_13'}
      if(grepl('VRC01_19', chain_id)){
        clone <- 'VRC01_19'}
      if(grepl('VRC01_H0306', chain_id)){
        clone <- 'VRC01_H0306'}
      if(grepl('VRC01_L0306', chain_id)){
        clone <- 'VRC01_L0306'}
      if(grepl('VRC01_H08', chain_id)){
        clone <- 'VRC01_H08'}
      if(grepl('VRC01_L08', chain_id)){
        clone <- 'VRC01_L08'}
      }
    }
  }
  # (TEMPORARY) to accomodate interrupted ('int') chains:
  if(grepl('int', chain_id)){
    clone <- paste(clone,'int', sep ='')
  }
  
  if(grepl('_con_', chain_id)){
    prior <- 'constant'
  }
  if(grepl('_log_', chain_id)){
    prior <- 'logistic'
  }
  if(grepl('_exp_', chain_id)){
    prior <- 'exponential'
  }
  results_file_path <- paste('../../results/rate_correlations/observed_lineages/', clone, '_', prior, '/',
                             chain_id,'_rate_correlations.csv',sep = '')
  plot_directory <- paste('../../figures/rate_correlations/observed_lineages/', clone, '_', prior, '/', sep = '')
  
  
}else{
  # If chain is from a simulated scenario
  time_units <- 'generations'
  scenario <- regmatches(chain_id,regexpr("scenario[0-9a-z]*", chain_id, perl=TRUE))
  results_file_path <- paste('../../results/rate_correlations/simulated_alignments/', chain_id,'/', 
                             chain_id, '_rate_correlations.csv',sep = '')
  plot_directory <- paste('../../figures/rate_correlations/simulated_alignments/', chain_id, '/', sep = '')
}

results_dataframe <- read.table(results_file_path, header = T, sep = ',')


print("PLOTTING 800 LINES PER PLOT, GGPLOT2 CAN'T SEEM TO HANDLE MANY MORE")

# ================ DEFINING GGPLOT2 PARAMETERS ===============


{
  title_size <- 20
  axis_title_size <- 30
  y_axis_text_size <- 25
  x_axis_text_size <- 25
  
  ylab_distance <- 20
  xlab_distance <- 30
  plot_margins <- c(0.3, 1, 0.3, 2)
}



# Master plotting function:
master_plot_time <- function(pars){
  do.call(base_plot_time, pars)
}

# ===== DEFINING BASE PLOTTING FUNCTION FOR RATES VS MOL-CLOCK DISTANCE LINES =====

base_plot_distance <- function(rate){
  
  ylabel <- switch(rate,
                   total_rate = "Total substitution rate",
                   total_rate_RC = "Total substitution rate (RC)",
                   S_rate = "Synonymous substitution rate",
                   N_rate = "Non-synonymous rate")
  
  intercept_column <- which(names(results_dataframe) == paste("intercept_",rate,"_vs_parent_distance", sep = ''))
  slope_column <- which(names(results_dataframe) == paste("slope_",rate,"_vs_parent_distance", sep = ''))
  r_column <- which(names(results_dataframe) == paste("r_",rate,"_vs_parent_distance", sep = ''))
  
  
  # Find mean and HPD for slope, intercept and R2
  mean_slope <- mean(results_dataframe[, slope_column])
  mean_intercept <- mean(results_dataframe[, intercept_column])
  mean_r2 <- mean((results_dataframe[, r_column])^2)
  
  slope_HPD_limits <- HPDinterval(as.mcmc(results_dataframe[, slope_column]), 0.95)
  intercept_HPD_limits <- HPDinterval(as.mcmc(results_dataframe[, intercept_column]), 0.95)
  r2_HPD_limits <- HPDinterval(as.mcmc((results_dataframe[, r_column])^2), 0.95)
  
  # Plot:
  
  # Subsampling: ggplot can't seem to plot a lot more than 800 lines:
  subsample <- sample(1:nrow(results_dataframe), min(800,nrow(results_dataframe)), replace = F)
  
  ymax <- mean(results_dataframe[subsample,paste('max_',rate,sep='')])
  if(mean_intercept >= 0){
    ymax <- min(c(ymax,1.5*mean_intercept))
  }
  ymin <- max(mean(results_dataframe[subsample,paste('min_',rate,sep='')]),
              min(results_dataframe[subsample, intercept_column]))
  if(mean_intercept>=0){
    ymin <- max(c(ymin,0.5*mean_intercept))
  }
  
  xmax <- mean(results_dataframe[subsample,'max_parent_distance'])
  xmin <- mean(results_dataframe[subsample,'min_parent_distance'])
  
  pl <- ggplot() + 
    scale_x_continuous(expand = c(0,0), limits = c(xmin,xmax)) +
    scale_y_continuous(expand = c(0,0), limits = c(ymin,ymax),label=scientific_10) + 
    theme_classic() + 
    xlab("Genetic distance from MRCA") + 
    ylab(ylabel) +
    theme(axis.title.y = element_text(size = axis_title_size,
                                      margin = margin(0,ylab_distance,0,0)),
          axis.title.x = element_text(size = axis_title_size,
                                      margin = margin(xlab_distance,0,0,0)),
          axis.text.x = element_text(size = x_axis_text_size),
          axis.text.y = element_text(size = y_axis_text_size),
          axis.ticks = element_line(size = 1.5),
          axis.ticks.length = unit(0.5, "cm"),
          plot.margin = unit(c(1, 1, 0.2, 0.2),"cm"),
          title = element_text(size = title_size),
          axis.line.x = element_line(colour="black", size = 1.5),
          axis.line.y = element_line(colour="black", size = 1.5)
          
    )
  
  for(i in subsample){
    slope = results_dataframe[i, slope_column]
    intercept <- results_dataframe[i, intercept_column]
    
    pl <- pl + geom_abline(intercept = intercept, slope = slope, alpha = 0.1)
    
  }
  
  pl <- pl +  geom_abline(intercept = mean_intercept, slope = mean_slope, 
                          colour = 'red', size = 2)
  
  pl <- pl + ggtitle(paste('slope = ', signif(mean_slope,2), 
                           " [", signif(slope_HPD_limits[1],3), 
                           ' , ', signif(slope_HPD_limits[2],3),
                           ']', sep = ''))
  
  
  #annotate("text", x = 75, y = 0.00175,
  #         label = paste('slope = ', signif(mean_slope,2), 
  #                       " [", signif(slope_HPD_limits[1],3), ' , ', signif(slope_HPD_limits[2],3), ']'
  #                       , sep = '')) + 
  
  #annotate("text", x = 75, y = 0.0017,
  #         label = paste('intercept = ', signif(mean_intercept,2), 
  #                       " [", signif(intercept_HPD_limits[1],3), ' , ', signif(intercept_HPD_limits[2],3), ']'
  #                       , sep = '')) + 
  
  #annotate("text", x = 75, y = 0.00165,
  #         label = paste('R-squared = ', signif(mean_r2,2), 
  #                       " [", signif(r2_HPD_limits[1],3), ' , ', signif(r2_HPD_limits[2],3), ']'
  #                       , sep = ''))
  
  return(pl)
}

master_plot_distance <- function(pars){
  do.call(base_plot_distance, pars)
}

# ===== DEFINING BASE PLOTTING FUNCTION FOR RATES VS ROBUST COUNTING DISTANCE LINES =====
base_plot_distance_RC <- function(rate){
  
  ylabel <- switch(rate,
                   total_rate = "Total substitution rate",
                   total_rate_RC = "Total substitution rate (RC)",
                   S_rate = "Synonymous substitution rate",
                   N_rate = "Non-synonymous rate")
  
  intercept_column <- which(names(results_dataframe) == paste("intercept_",rate,"_vs_parent_distance_RC", sep = ''))
  slope_column <- which(names(results_dataframe) == paste("slope_",rate,"_vs_parent_distance_RC", sep = ''))
  r_column <- which(names(results_dataframe) == paste("r_",rate,"_vs_parent_distance_RC", sep = ''))
  
  
  # Find mean and HPD for slope, intercept and R2
  mean_slope <- mean(results_dataframe[, slope_column])
  mean_intercept <- mean(results_dataframe[, intercept_column])
  mean_r2 <- mean((results_dataframe[, r_column])^2)
  
  slope_HPD_limits <- HPDinterval(as.mcmc(results_dataframe[, slope_column]), 0.95)
  intercept_HPD_limits <- HPDinterval(as.mcmc(results_dataframe[, intercept_column]), 0.95)
  r2_HPD_limits <- HPDinterval(as.mcmc((results_dataframe[, r_column])^2), 0.95)
  
  # Plot:
  
  # Subsampling: ggplot can't seem to plot a lot more than 800 lines:
  subsample <- sample(1:nrow(results_dataframe), min(800,nrow(results_dataframe)), replace = F)
  
  ymax <- mean(results_dataframe[subsample,paste('max_',rate,sep='')])
  if(mean_intercept >= 0){
    ymax <- min(c(ymax,1.5*mean_intercept))
  }
  ymin <- max(mean(results_dataframe[subsample,paste('min_',rate,sep='')]),
              min(results_dataframe[subsample, intercept_column]))
  if(mean_intercept>=0){
    ymin <- max(c(ymin,0.5*mean_intercept))
  }
  
  xmax <- mean(results_dataframe[subsample,'max_parent_distance_RC'])
  xmin <- mean(results_dataframe[subsample,'min_parent_distance_RC'])
  
  pl <- ggplot() + 
    scale_x_continuous(expand = c(0,0), limits = c(xmin,xmax)) +
    scale_y_continuous(expand = c(0,0), limits = c(ymin, ymax),label=scientific_10) + 
    theme_classic() + 
    xlab("Genetic distance (RC) from MRCA") + 
    ylab(ylabel) +
    theme(axis.title.y = element_text(size = axis_title_size,
                                      margin = margin(0,ylab_distance,0,0)),
          axis.title.x = element_text(size = axis_title_size,
                                      margin = margin(xlab_distance,0,0,0)),
          axis.text.x = element_text(size = x_axis_text_size),
          axis.text.y = element_text(size = y_axis_text_size),
          axis.ticks = element_line(size = 1.5),
          axis.ticks.length = unit(0.5, "cm"),
          plot.margin = unit(c(1, 1, 0.2, 0.2),"cm"),
          title = element_text(size = title_size),
          axis.line.x = element_line(colour="black", size = 1.5),
          axis.line.y = element_line(colour="black", size = 1.5)
          
    )
  
  for(i in subsample){
    slope = results_dataframe[i, slope_column]
    intercept <- results_dataframe[i, intercept_column]
    
    pl <- pl + geom_abline(intercept = intercept, slope = slope, alpha = 0.1)
    
  }
  
  pl <- pl +  geom_abline(intercept = mean_intercept, slope = mean_slope, 
                          colour = 'red', size = 2)
  pl <- pl + ggtitle(paste('slope = ', signif(mean_slope,2), 
                           " [", signif(slope_HPD_limits[1],3), 
                           ' , ', signif(slope_HPD_limits[2],3),
                           ']', sep = ''))
  
  #annotate("text", x = 75, y = 0.00175,
  #         label = paste('slope = ', signif(mean_slope,2), 
  #                       " [", signif(slope_HPD_limits[1],3), ' , ', signif(slope_HPD_limits[2],3), ']'
  #                       , sep = '')) + 
  
  #annotate("text", x = 75, y = 0.0017,
  #         label = paste('intercept = ', signif(mean_intercept,2), 
  #                       " [", signif(intercept_HPD_limits[1],3), ' , ', signif(intercept_HPD_limits[2],3), ']'
  #                       , sep = '')) + 
  
  #annotate("text", x = 75, y = 0.00165,
  #         label = paste('R-squared = ', signif(mean_r2,2), 
  #                       " [", signif(r2_HPD_limits[1],3), ' , ', signif(r2_HPD_limits[2],3), ']'
  #                       , sep = ''))
  
  return(pl)
}

master_plot_distance_RC <- function(pars){
  do.call(base_plot_distance_RC, pars)
}

# ===== DEFINING BASE PLOTTING FUNCTION FOR RATES VS MUTABILITY LINES

base_plot_mutability <- function(rate, metric,region){
   ylabel <- switch(rate,
                   total_rate = "Total substitution rate",
                   total_rate_RC = "Total substitution rate (RC)",
                   S_rate = "Synonymous substitution rate",
                   N_rate = "Non-synonymous rate")
   xlabel <- switch(metric,
                         S5F = paste('Mean S5F mutability in ', region, sep =''),
                         '7M' = paste('Mean 7-mer mutability in ', region, sep =''),
                         HS = paste('Number of hotspots in ', region, sep =''),
                         CS = paste('Number of coldspots in ', region, sep =''),
                         OHS = paste('Number of overlapping hotspots in ', region, sep ='')
                         )
   
   
   
  
    subsample <- sample(1:nrow(results_dataframe), min(800,nrow(results_dataframe)), replace = F)
   
    intercept_column <- which(names(results_dataframe) == paste('intercept_',rate,'_vs_',metric,'_',region, sep =''))
    slope_column <- which(names(results_dataframe) == paste('slope_',rate,'_vs_',metric,'_',region, sep =''))
    r_column <- which(names(results_dataframe) == paste('r_',rate,'_vs_',metric,'_',region, sep =''))
    
    # Find mean and HPD for slope, intercept and R2
    mean_slope <- mean(results_dataframe[, slope_column])
    mean_intercept <- mean(results_dataframe[, intercept_column])
    mean_r2 <- mean((results_dataframe[, r_column])^2)
    
    slope_HPD_limits <- HPDinterval(as.mcmc(results_dataframe[, slope_column]), 0.95)
    intercept_HPD_limits <- HPDinterval(as.mcmc(results_dataframe[, intercept_column]), 0.95)
    r2_HPD_limits <- HPDinterval(as.mcmc((results_dataframe[, r_column])^2), 0.95)
    
    ymax <- mean(results_dataframe[subsample,paste('max_',rate,sep='')])
    #ymax <- 1.1*(mean_intercept)
    ymin <- mean(results_dataframe[subsample,paste('min_',rate,sep='')])
    ymin <- ifelse(ymin<ymax,ymin,0)
    
    xmax <- mean(results_dataframe[subsample,paste('max_',metric,'_',region,sep='')])
    xmin <- mean(results_dataframe[subsample,paste('min_',metric,'_',region,sep='')])
    
    # Generate base plot
    pl <- ggplot() +
      xlab(xlabel) + 
      ylab(ylabel) +
      scale_x_continuous(expand = c(0,0), limits = c(xmin,xmax)) +
      scale_y_continuous(expand = c(0,0), limits = c(ymin,ymax),label=scientific_10) + 
      theme_classic() +
      theme(axis.title.y = element_text(size = axis_title_size,
                                        margin = margin(0,ylab_distance,0,0)),
            axis.title.x = element_text(size = axis_title_size,
                                        margin = margin(xlab_distance,0,0,0)),
            axis.text.x = element_text(size = x_axis_text_size),
            axis.text.y = element_text(size = y_axis_text_size),
            axis.ticks = element_line(size = 1.5),
            axis.ticks.length = unit(0.5, "cm"),
            plot.margin = unit(c(1, 1, 0.2, 0.2),"cm"),
            title = element_text(size = title_size),
            axis.line.x = element_line(colour="black", size = 1.5),
            axis.line.y = element_line(colour="black", size = 1.5)
            
      )
    
    
    for(i in subsample){
      slope = results_dataframe[i, slope_column]
      intercept <- results_dataframe[i, intercept_column]
      
      pl <- pl + geom_abline(intercept = intercept, slope = slope, alpha = 0.1)
      
    }
    
    pl <- pl +  geom_abline(intercept = mean_intercept, slope = mean_slope, 
                            colour = 'red', size = 2)
    pl <- pl + ggtitle(paste('slope = ', signif(mean_slope,2), 
                             " [", signif(slope_HPD_limits[1],3), 
                             ' , ', signif(slope_HPD_limits[2],3),
                             ']', sep = ''))
    

  return(pl)
}

master_plot_mutability <- function(pars){
  do.call(base_plot_mutability, pars)
}



# ============== PLOTTING RATES VS TIME AND DISTANCE LINES ===============
for(rate in c('total_rate','total_rate_RC','S_rate','N_rate')){
  pdf_path_time <- paste(plot_directory, chain_id,'_', rate, '_vs_time_lines.pdf', sep = '')
  pdf_path_distance <- paste(plot_directory, chain_id,'_',rate, '_vs_distance_lines.pdf', sep = '')
  pdf_path_distance_RC <- paste(plot_directory, chain_id,'_',rate, '_vs_distance_RC_lines.pdf', sep = '')
  
  pdf(pdf_path_time, width = 10, height = 10)
    plot(base_plot_time(rate))
  dev.off()
#   
#   pdf(pdf_path_distance, width = 10, height = 10)
#     plot(base_plot_distance(rate))
#   dev.off()
#   
#   pdf(pdf_path_distance_RC, width = 10, height = 10)
#     plot(base_plot_distance_RC(rate))
#   dev.off()
}




# # ============== PLOTTING RATES VS MUTABILITY LINES ==============
# 
# for(metric in c('S5F','7M','HS','CS','OHS')){
#   # If chain is from a simulated scenario, skip plots for FRs and CDRs
#   if(grepl("scenario", chain_id) == TRUE){
# 
#     pars_list <- list(list('total_rate',metric,'WS'),
#                              list('total_rate_RC',metric,'WS'),
#                              list('S_rate',metric,'WS'),
#                              list('N_rate',metric,'WS')
#                             )
#     
#     mutability_plots <- lapply(pars_list, master_plot_mutability)
#     
#     args_list_mutability <- c(mutability_plots, 2, 2)
#     names(args_list_mutability) <- c(letters[1:4], "nrow", "ncol")
#     
#     pdf_path_mutability <- paste(plot_directory, chain_id, '_rates_vs_', metric,'_WS',sep='')
#     pdf(pdf_path_mutability, width = 16, height = 16)
#       plot(mutability_plots)
#     dev.off()
#   }else{
#     for(region in c('WS','FR','CDR')){
#       pars_list <- list(list('total_rate',metric,region),
#                          list('total_rate_RC',metric,region),
#                          list('S_rate',metric,region),
#                          list('N_rate',metric,region)
#                          )
#  
#       mutability_plots <- lapply(pars_list, master_plot_mutability)
#       args_list_mutability <- c(mutability_plots, 2, 2)
#       names(args_list_mutability) <- c(letters[1:4], "nrow", "ncol")
#       
#       pdf_path_mutability <- paste(plot_directory, chain_id, '_rates_vs_', metric ,'_',region,'_lines.pdf', sep = '')
#       pdf(pdf_path_mutability, width = 16, height = 16)
#         do.call(grid.arrange, args_list_mutability)
#       dev.off()
#     }
#   }
# }



