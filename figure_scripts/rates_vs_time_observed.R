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
          plot.margin = unit(c(1, 0.1, 0.1, 0.1),"cm")
          )
    #scale_y_continuous(label=scientific_10)
  
  for(tree in tree_sample){
    dataframe_subset <- results_dataframe[results_dataframe$tree == tree,]
    x <- dataframe_subset[,'parent_time']
    y <- dataframe_subset[,paste(rate,'_to_child',sep='')]
    factor <- factor(dataframe_subset[,'branch_is_terminal'])
    
    pl <- pl + geom_point(data=data.frame(x,y,factor), 
                          aes(x=x,y=y, colour = factor),
                          size = 2, alpha = 0.8)
    
    
  }
  
  return(pl)
}

plot_list_total_rate <- c()
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
  
  plot_list_total_rate[[clone_labels[i]]] <- base_plot_time(results_dataframe, 'total_rate', time_units = time_units[i])
  plot_list_total_rate_RC[[clone_labels[i]]] <- base_plot_time(results_dataframe, 'total_rate_RC', time_units = time_units[i])
  plot_list_S_rate[[clone_labels[i]]] <- base_plot_time(results_dataframe, 'S_rate', time_units = time_units[i])
  plot_list_N_rate[[clone_labels[i]]] <- base_plot_time(results_dataframe, 'N_rate', time_units = time_units[i])
  
}

# Plot with results for total substitution rates from molecular clock
total_rate_row_1 <- plot_grid(plot_list_total_rate[['CH103 (H)']], plot_list_total_rate[['CH103 (L)']],  plot_list_total_rate[['VRC26 (H)']],
          nrow = 1, labels = clone_labels[1:3])
total_rate_row_2 <- plot_grid(plot_list_total_rate[['VRC26 (L)']], plot_list_total_rate[['VRC01-13 (H)']],  plot_list_total_rate[['VRC01-01 (H)']],
                              nrow = 1, labels = clone_labels[4:6])
total_rate_row_3 <- plot_grid(plot_list_total_rate[['VRC01-19 (H)']], ggplot() + geom_blank(), ggplot() + geom_blank(), labels = c(clone_labels[7],'',''), nrow = 1)

total_rate_plot <- plot_grid(total_rate_row_1, total_rate_row_2, total_rate_row_3, nrow = 3)
         
png('total_rate_vs_time_observed.png', width = 10, height =10, units = 'in', res = 400)
plot(total_rate_plot)
dev.off()

# Plot with results for total substitution rates from robust counting
total_rate_RC_row_1 <- plot_grid(plot_list_total_rate_RC[['CH103 (H)']], plot_list_total_rate_RC[['CH103 (L)']],  plot_list_total_rate_RC[['VRC26 (H)']],
                              nrow = 1, labels = clone_labels[1:3])
total_rate_RC_row_2 <- plot_grid(plot_list_total_rate_RC[['VRC26 (L)']], plot_list_total_rate_RC[['VRC01-13 (H)']],  plot_list_total_rate_RC[['VRC01-01 (H)']],
                              nrow = 1, labels = clone_labels[4:6])
total_rate_RC_row_3 <- plot_grid(plot_list_total_rate_RC[['VRC01-19 (H)']], ggplot() + geom_blank(), ggplot() + geom_blank(), labels = c(clone_labels[7],'',''), nrow = 1)

total_rate_RC_plot <- plot_grid(total_rate_RC_row_1, total_rate_RC_row_2, total_rate_RC_row_3, nrow = 3)

png('total_rate_RC_vs_time_observed.png', width = 10, height = 10, units = 'in', res = 400)
plot(total_rate_RC_plot)
dev.off()

# Plot with results for synonymous substitution rates from robust counting
S_rate_row_1 <- plot_grid(plot_list_S_rate[['CH103 (H)']], plot_list_S_rate[['CH103 (L)']],  plot_list_S_rate[['VRC26 (H)']],
                              nrow = 1, labels = clone_labels[1:3])
S_rate_row_2 <- plot_grid(plot_list_S_rate[['VRC26 (L)']], plot_list_S_rate[['VRC01-13 (H)']],  plot_list_S_rate[['VRC01-01 (H)']],
                              nrow = 1, labels = clone_labels[4:6])
S_rate_row_3 <- plot_grid(plot_list_S_rate[['VRC01-19 (H)']], ggplot() + geom_blank(), ggplot() + geom_blank(), labels = c(clone_labels[7],'',''), nrow = 1)

S_rate_plot <- plot_grid(S_rate_row_1, S_rate_row_2, S_rate_row_3, nrow = 3)

png('S_rate_vs_time_observed.png', width = 10, height = 10, units = 'in', res = 400)
plot(S_rate_plot)
dev.off()

# Plot with results for non-synonymous substitution rates from robust counting
N_rate_row_1 <- plot_grid(plot_list_N_rate[['CH103 (H)']], plot_list_N_rate[['CH103 (L)']],  plot_list_N_rate[['VRC26 (H)']],
                          nrow = 1, labels = clone_labels[1:3])
N_rate_row_2 <- plot_grid(plot_list_N_rate[['VRC26 (L)']], plot_list_N_rate[['VRC01-13 (H)']],  plot_list_N_rate[['VRC01-01 (H)']],
                          nrow = 1, labels = clone_labels[4:6])
N_rate_row_3 <- plot_grid(plot_list_N_rate[['VRC01-19 (H)']], ggplot() + geom_blank(), ggplot() + geom_blank(), labels = c(clone_labels[7],'',''), nrow = 1)

N_rate_plot <- plot_grid(N_rate_row_1, N_rate_row_2, N_rate_row_3, nrow = 3)

png('N_rate_vs_time_observed.png', width = 10, height = 10, units = 'in', res = 400)
plot(N_rate_plot)
dev.off()
