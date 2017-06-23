library('ggplot2')


# Load ggplot parameters
source('ggplot_parameters.R')

points_file_path <- '../results/mutability_vs_time/observed_lineages/CH103_constant/CH103_con_run1a_mutability_vs_time_points.csv'
regressions_file_path <- '../results/mutability_vs_time/observed_lineages/CH103_constant/CH103_con_run1a_mutability_vs_time_correlations.csv'
  
# ====== READ DATAFRAMES WITH POINTS and REGRESSION RESULTS =====
points_dataframe <- read.table(points_file_path, header = T, sep = ',')

# Sample of 100 trees
sample_size <- 100
tree_sample <- sample(unique(points_dataframe$tree), sample_size, replace = F)

# ====== READ DATAFRAME WITH REGRESSION RESULTS ======
regressions_dataframe <- read.table(regressions_file_path, header = T, sep = ',')


mean_intercept <- mean(regressions_dataframe[, 'intercept_S5F_WS_vs_node_time'])

mean_slope <- mean(regressions_dataframe[, 'slope_S5F_WS_vs_node_time'])


pl <- ggplot(points_dataframe, aes(x = node_time, y = S5F_WS)) +
  xlab("Time since unmutated ancestor (weeks)") +
  ylab("Mean S5F mutability") + 
  ggplot_theme

for(tree in tree_sample){
  dataframe_subset <- points_dataframe[points_dataframe$tree == tree,]
  x <- dataframe_subset[,'node_time']
  y <- dataframe_subset[,'S5F_WS']
  factor <- factor(dataframe_subset[,'node_is_tip'])
  
  pl <- pl + geom_point(data=data.frame(x,y,factor), 
                        aes(x=x,y=y, colour = factor), size = 0.4,
                        alpha = 0.8)
  
  
}

pl <- pl + scale_colour_manual(values=c("grey25",'royalblue1'))

# Add "mean line" (line with mean slope and mean intercept)
pl <- pl + geom_abline(slope=mean_slope, 
                       intercept=mean_intercept, 
                       colour = 'red', linetype = 1, size = 0.5) 

png('CH103_overview_S5F_vs_time.png', width = 3.43, height = 3, units = 'in', res = 400)
plot(pl)
dev.off()

pdf('CH103_overview_S5F_vs_time.pdf', width = 3.43, height = 3)
plot(pl)
dev.off()

