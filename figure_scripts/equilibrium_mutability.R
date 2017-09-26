library(cowplot)
source('ggplot_parameters.R')

results <- read.csv('../results/equilibrium_mutability/equilibrium_mutability.csv', header = T)
results$replicate_trajectory <- factor(results$replicate_trajectory)

# Computing average mutability at each generation for each initial mutability:
average_mutability <- c()
initial_mutability <- c()
for(init_mutability in c("low","high")){
  for(gen in unique(results$generation)){
    mutability_values <- results$S5F_mutability[results$generation == gen & results$initial_mutability == init_mutability]
    average_mutability <- c(average_mutability, mean(mutability_values))
    initial_mutability <- c(initial_mutability, init_mutability)
  }
}
means_dataframe <- data.frame(average_mutability, initial_mutability, generation = unique(results$generation))
rm(initial_mutability)
rm(average_mutability)



base_pl <- ggplot(results, aes(x = generation, y = S5F_mutability,
                               colour = initial_mutability)) + xlim(0,5000)
for(rep in unique(results$replicate_trajectory)){
  
  results_subset <- results[results$replicate_trajectory == rep,]
  
  base_pl <- base_pl + geom_line(data=results_subset, alpha = 0.1)
}

base_pl <- base_pl + scale_color_manual(name = 'Initial mutability',
                                        values = c('firebrick1','royalblue1'))

base_pl <- base_pl + geom_point(data = data.frame(x=0, 
                                                  y = results$S5F_mutability[results$generation==0 & results$initial_mutability == 'low'][1]),
                                colour = 'royalblue1', aes(x=x,y=y), size = 2)

base_pl <- base_pl + geom_point(data = data.frame(x=0, 
                                                  y = results$S5F_mutability[results$generation==0 & results$initial_mutability == 'high'][1]),
                                colour = 'firebrick1', aes(x=x,y=y), size = 2)

base_pl <- base_pl + xlab('Generation') + ylab('S5F mutability')

base_pl <- base_pl + theme(legend.position = 'top')

# Add mean trajectories
base_pl <- base_pl + geom_line(data = means_dataframe,
                               aes(x = generation, y = average_mutability, colour = initial_mutability))


pdf('equilibrium_mutability.pdf',height = 7, width = 7)
plot(base_pl)
dev.off()
