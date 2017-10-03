library('ggplot2')
library('coda')
library('cowplot')
# Load ggplot parameters
source('ggplot_parameters.R')

results_directory <- '../results/mutability_vs_time/observed_lineages/'

points_files <- c(CH103 = 'CH103_constant/CH103_con_run1a_mutability_vs_time_points.csv',
                  CH103L = 'CH103L_constant/CH103L_con_run1a_mutability_vs_time_points.csv',
                  VRC26 = 'VRC26int_constant/VRC26int_con_run1a_mutability_vs_time_points.csv',
                  VRC26L = 'VRC26L_constant/VRC26L_con_run1a_mutability_vs_time_points.csv',
                  VRC01_13 = 'VRC01_13_logistic/VRC01_13_log_run1a_mutability_vs_time_points.csv',
                  VRC01_01 = 'VRC01_01_logistic/VRC01_01_log_run1a_mutability_vs_time_points.csv',
                  VRC01_19 = 'VRC01_19_logistic/VRC01_19_log_run1a_mutability_vs_time_points.csv'
)

for(i in 1:length(points_files)){
  points_files[i] <- paste(results_directory, points_files[i], sep = '')
}

regressions_files <- gsub('points','correlations', points_files) 

points_dataframes_list <- lapply(points_files, FUN = read.csv, header = T)
regressions_dataframes_list <- lapply(regressions_files, FUN = read.csv, header = T)


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
  
  mean_slope_all <- mean(slope_values_all)
  mean_slope_obs <- mean(slope_values_obs)
  
  # Find 95% HPD interval for slope (all nodes only):
  slope_all_HPD_limits <- HPDinterval(as.mcmc(slope_values_all), 
                                  0.95)
  slope_all_llim <- slope_all_HPD_limits[1]
  slope_all_ulim <- slope_all_HPD_limits[2]
  
  # Use a solid line for the mean regression line including all nodes
  line_type_all <- 1
  
  # Use a dashed line for the regression line including obs nodes only
  line_type_obs <- 5
  
  # Maximum time
  max_time <-  max(points_dataframe[points_dataframe$tree %in% tree_sample,
                                    'node_time'])
  
  # Plot:
  pl <- ggplot() + 
    xlab(paste("Time at node (",time_units,')',sep='')) + 
    ylab(ylabel) +
    ggplot_theme +
    theme(plot.title = element_text(size = 10),
          plot.margin = margin(24,12,6,6,'pt')) +
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
  
  # Add mean slope and 95% Cred.I (for regressions involving all points)
  pl <- pl + ggtitle(parse(text = paste('Slope == {',
                                        scientific_10(mean_slope_all),
                                        '} (',
                                        scientific_10(slope_all_llim),
                                        ',',
                                        scientific_10(slope_all_ulim),
                                        ')',
                                        sep = ''
                                        )
                           )
                     )
  
  # Subplot with posterior distribution of slope
  
  # # Expand axis to fit inset plots
  # pl <- pl + scale_x_continuous(limits = c(0, 1.25*max_time))
  # 
  # # Data-frame for shading region under density curve in inset plot
  # shading_dataframe <- data.frame(slope_values_all)
  # shading_dataframe <- with(density(shading_dataframe$slope_values_all,na.rm=T), data.frame(x, y))
  # 
  # subpl <- ggplot(data.frame(slope=slope_values_all), aes(x=slope)) + 
  #   
  #   # Added shaded region under slope posterior density curve
  #   geom_area(data = shading_dataframe, 
  #                     mapping = aes(y = ifelse(x>slope_all_llim & x< slope_all_ulim, y, 0),x), 
  #                     fill = "snow2") +
  #   # Add curve:
  #   geom_density() + 
  #   
  #   # Add line indicating mean slope for all points:
  #   geom_vline(xintercept = mean_slope_all, linetype = line_type_all, colour = 'red') +
  #   
  #   # Add line indicating slope for obs points only:
  #   geom_vline(xintercept = mean_slope_obs, linetype = line_type_obs, colour = 'red') +
  # 
  #   theme(axis.text=element_text(size=x_axis_text_size_subplot),
  #         axis.title=element_text(size=axis_title_size_subplot, margin = margin(0,0,0,0)),
  #         axis.line.x = element_line(colour="black", size = 0.3),
  #         axis.line.y = element_line(colour="black",  size = 0.3),
  #         plot.title = element_text(size=2)
  #   ) + 
  #   ylab("Density") +
  #   xlab("Slope") +
  #   #expand_limits(x = 0) +
  #   scale_x_continuous(labels = function(x) format(x, scientific=TRUE)) +
  #   #scale_x_continuous(limits = c(-0.001,0.001)) +
  #   scale_y_continuous(expand = c(0,0)) +
  #   geom_vline(xintercept = 0, linetype = 2)
  # 
  # 
  # return(list('main_plot' = pl, 'subplot' = subpl))
  return(pl)
}

# Time units list (for mapply):
time_units_list <- list(CH103 = 'weeks', CH103L = 'weeks', VRC26 = 'weeks', VRC26L = 'weeks',
                        VRC01_13 ='months', VRC01_19 = 'months', VRC01_01 = 'months')

# Labels for each clone
plot_grid_labels <- c() 
for(clone_name in names(points_dataframes_list)){
  plot_grid_labels <- c(plot_grid_labels, switch(clone_name, 
                                                 'CH103' = 'CH103 (H)',
                                                 'CH103L' = 'CH103 (L)',
                                                 'VRC26' = 'VRC26 (H)',
                                                 'VRC26L' = 'VRC26 (L)',
                                                 'VRC01_13' = 'VRC01-13 (H)',
                                                 'VRC01_01' = 'VRC01-01 (H)',
                                                 'VRC01_19' = 'VRC01-19 (H)'
  ))
}


# For all metrics and regions:
for(metric in c('S5F','HS','OHS','CS','X7M')){
  for(region in c('WS','FR','CDR')){
    plots <- mapply(FUN = base_plot, points_dataframe = points_dataframes_list, 
                    regressions_dataframe = regressions_dataframes_list,
                    time_units = time_units_list,
                    MoreArgs = list(metric = metric, region = region, theme_specs = ggplot_theme),
                    SIMPLIFY = FALSE
                    )
    
    png(paste(ifelse(metric == 'X7M','7M',metric), '_', region, '_vs_time_all_lineages.png', sep  =''), 
        height = 3500, width = 3500, res = 300)
    plot(plot_grid(plotlist = plots, labels = plot_grid_labels))
    dev.off()
  }
}


