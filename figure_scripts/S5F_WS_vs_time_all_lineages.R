library('ggplot2')
library('coda')
# Load ggplot parameters
source('ggplot_parameters.R')

points_file_path <- '../results/mutability_vs_time/observed_lineages/CH103_constant/CH103_con_run1a_mutability_vs_time_points.csv'



regressions_file_path <- '../results/mutability_vs_time/observed_lineages/CH103_constant/CH103_con_run1a_mutability_vs_time_correlations.csv'

points_dataframe <- read.csv(points_file_path)
regressions_dataframe <- read.csv(regressions_file_path)

# Function for generating the plot for each dataset
base_plot <- function(points_dataframe, regressions_dataframe, metric, 
                      region, theme_specs, time_units, tree_sample_size = 100){
  ylabel <- switch(metric,
                   S5F = 'Mean S5F mutability',
                   X7M = 'Mean 7-mer mutability',
                   HS = 'Number of hotspots',
                   CS = 'Number of coldspots',
                   OHS = 'Number of overlapping hotspots')
  ylabel <- paste(ylabel, ' (', region, ')', sep = '')
  
  tree_sample <- sample(unique(points_dataframe$tree), tree_sample_size, replace = F)
  
  # Finding mean slope and mean intercept from the posterior distribution
  # (with a little workaround for the naming issue with 7M)
  mean_intercept_all <- mean(regressions_dataframe[, paste('intercept_',
                                                       ifelse(metric=='X7M',substr(metric,2,3),metric),
                                                       '_',region,'_vs_node_time',sep='')])
  
  #Mean intercept for regressions including observed nodes only
  mean_intercept_obs <- mean(regressions_dataframe[, paste('intercept_',
                                                       ifelse(metric=='X7M',substr(metric,2,3),metric),
                                                       '_',region,'_vs_node_time_obs_only',sep='')])
  
  # Slope values for regressions including all nodes
  slope_values_all <- regressions_dataframe[, paste('slope_',ifelse(metric=='X7M',substr(metric,2,3),metric)
                                                ,'_',region,'_vs_node_time',sep='')]

  # Slope values for regressions including observed nodes only
  slope_values_obs <- regressions_dataframe[, paste('slope_',ifelse(metric=='X7M',substr(metric,2,3),metric)
                                                         ,'_',region,'_vs_node_time_obs_only',sep='')]
  
  # Slope values for obs. nodes only should all be the same. Stop if they are not
  stopifnot(sd(slope_values_obs) == 0)
  
  mean_slope_all <- mean(slope_values_all)
  mean_slope_obs <- mean(slope_values_obs)
  
  # Find 95% HPD interval for slope (all nodes only):
  slope_all_HPD_limits <- HPDinterval(as.mcmc(slope_values_all), 
                                  0.95)
  
  # Use a solid line for the mean regression line including all nodes
  line_type_all <- 1
  
  # Use a dashed line for the regression line including obs nodes only
  line_type_obs <- 5
  
  # Plot:
  pl <- ggplot() + 
    xlab(paste("Time at node (",time_units,')',sep='')) + 
    ylab(ylabel) +
    ggplot_theme +
    # Internal: black (grey25), Terminal: blue ('royalblue1')
    scale_colour_manual(values=c("grey25",'royalblue1'))

  for(tree in tree_sample){
    dataframe_subset <- points_dataframe[points_dataframe$tree == tree,]
    x <- dataframe_subset[,'node_time']
    y <- dataframe_subset[,paste(metric,'_',region, sep='')]
    factor <- factor(dataframe_subset[,'node_is_tip'])
    
    pl <- pl + geom_point(data=data.frame(x,y,factor), 
                          aes(x=x,y=y, colour = factor), size = 1,
                          alpha = 0.8)
    
    
  }
  
  # Add "mean line" (line with mean slope and mean intercept) for regressions including all nodes
  pl <- pl + geom_abline(slope=mean_slope_all, 
                         intercept=mean_intercept_all, 
                         colour = 'red', linetype = line_type_all) 
  
  # Add mean line for regressions including observed nodes only
  
  pl <- pl + geom_abline(slope=mean_slope_obs, 
                         intercept=mean_intercept_obs, 
                         colour = 'red', linetype = line_type_obs) 
  
  subpl <-   subpl <- ggplot(data.frame(slope=slope_values), aes(x=slope)) +
    ylim <- c(0,)
    geom_segment(aes(x = slope_all_HPD_limits[1], ))
  
  
  # Subplot with posterior distribution of slope
  subpl <- ggplot(data.frame(slope=slope_values), aes(x=slope)) + 
    geom_density() + 
    theme(axis.text=element_text(size=x_axis_text_size_subplot),
          axis.title=element_text(size=axis_title_size_subplot, margin = margin(0,0,0,0)),
          axis.line.x = element_line(colour="black", size = 0.3),
          axis.line.y = element_line(colour="black",  size = 0.3),
          plot.title = element_text(size=2)
    ) + 
    ylab("Density") +
    xlab("Slope") +
    expand_limits(x = 0) +
    scale_y_continuous(expand = c(0,0)) +
    geom_vline(xintercept = 0, linetype = 2)

  
  return(pl)
}

# Master plotting function:
master_plot_time <- function(pars){
  do.call(base_plot_time, pars)
}
