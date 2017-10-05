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

# ===== DEFINING BASE PLOTTING FUNCTION FOR RATES VS TIME LINES =====

base_plot_time <- function(results_dataframe, rate){
  
  time_units <- 'generations'
  
  ylabel <- switch(rate,
                   total_rate = "Total substitution rate",
                   total_rate_RC = "Total substitution rate (RC)",
                   S_rate = "Synonymous substitution rate",
                   N_rate = "Non-synonymous rate")
  
  intercept_column <- which(names(results_dataframe) == paste("intercept_",rate,"_vs_parent_time", sep = ''))
  slope_column <- which(names(results_dataframe) == paste("slope_",rate,"_vs_parent_time", sep = ''))
  r_column <- which(names(results_dataframe) == paste("r_",rate,"_vs_parent_time", sep = ''))
  
  
  # Find mean and HPD for slope, intercept and R2
  mean_slope <- mean(results_dataframe[, slope_column])
  mean_intercept <- mean(results_dataframe[, intercept_column])
  mean_r2 <- mean((results_dataframe[, r_column])^2)
  
  slope_HPD_limits <- HPDinterval(as.mcmc(results_dataframe[, slope_column]), 0.95)
  intercept_HPD_limits <- HPDinterval(as.mcmc(results_dataframe[, intercept_column]), 0.95)
  r2_HPD_limits <- HPDinterval(as.mcmc((results_dataframe[, r_column])^2), 0.95)
  
  # Plot:
  
  # Subsampling: ggplot can't seem to plot a lot more than 500 lines:
  subsample <- sample(1:nrow(results_dataframe), min(500,nrow(results_dataframe)), replace = F)
  
  ymax <- mean(results_dataframe[subsample,paste('max_',rate,sep='')])
  if(mean_intercept >= 0){
    ymax <- min(c(ymax,1.5*mean_intercept))
  }
  ymin <- max(mean(results_dataframe[subsample,paste('min_',rate,sep='')]),
              min(results_dataframe[subsample, intercept_column]))
  if(mean_intercept>=0){
    ymin <- max(c(ymin,0.5*mean_intercept))
  }
  
  xmax <- mean(results_dataframe[subsample,'max_parent_time'])
  xmin <- mean(results_dataframe[subsample,'min_parent_time'])
  
  
  
  pl <- ggplot() + 
    scale_x_continuous(expand = c(0,0), limits = c(xmin,xmax)) +
    scale_y_continuous(expand = c(0,0), limits = c(ymin, ymax),label=scientific_10) + 
    xlab(paste("Time since MRCA (", time_units, ')', sep = '')) + 
    ylab(ylabel) +
    ggplot_theme +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10),
          plot.margin = unit(c(1, 1, 1, 1),"cm"))
    
  for(i in subsample){
    slope = results_dataframe[i, slope_column]
    intercept <- results_dataframe[i, intercept_column]
    
    pl <- pl + geom_abline(intercept = intercept, slope = slope, alpha = 0.1)
    
  }
  
  pl <- pl +  geom_abline(intercept = mean_intercept, slope = mean_slope, 
                          colour = 'red', size = 2)
  
  pl <- pl + ggtitle(paste('slope = ', signif(mean_slope,2), 
                           " [", signif(slope_HPD_limits[1],3), 
                           ', ', signif(slope_HPD_limits[2],3),
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

plot_list_total_rate_RC <- c()
plot_list_S_rate <- c()
plot_list_N_rate <- c()




for(scenario in c('2c','2d','2e','3a','3b','3c','4a','4b','4c','4d')){
  file_path <- paste('../results/rates_vs_time/simulated_alignments/scenario',
                     scenario, '_rep2/scenario', scenario, '_rep2_rates_vs_time_correlations.csv', sep = '')
  
  results_dataframe <- read.table(file_path, header = T, sep = ',')
  
  plot_list_total_rate_RC[[scenario]] <- base_plot_time(results_dataframe, 'total_rate_RC')
  plot_list_S_rate[[scenario]] <- base_plot_time(results_dataframe, 'S_rate')
  plot_list_N_rate[[scenario]] <- base_plot_time(results_dataframe, 'N_rate')
  
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

# Plot for uniform mutability simulations obtained under different r, K, s:
uniform_mutability_plot_newpars <- plot_grid(plot_list_total_rate_RC[['4c']],plot_list_S_rate[['4c']], plot_list_N_rate[['4c']], 
                                     plot_list_total_rate_RC[['4d']],plot_list_S_rate[['4d']], plot_list_N_rate[['4d']],
                                     plot_list_total_rate_RC[['4e']],plot_list_S_rate[['4e']], plot_list_N_rate[['4e']],
                                     labels = c('', '30% simulated decrease','','', '40% simulated decrease','','', '50% simulated decrease','')
)
save_plot("rates_vs_time_uniform_newpars.pdf", uniform_mutability_plot,
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