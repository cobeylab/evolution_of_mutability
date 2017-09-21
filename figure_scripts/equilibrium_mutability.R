library(cowplot)
source('ggplot_parameters.R')

results <- read.csv('../results/equilibrium_mutability/equilibrium_mutability.csv', header = T)
results$replicate_trajectory <- factor(results$replicate_trajectory)

base_pl <- ggplot(results, aes(x = generation, y = S5F_mutability,
                               colour = initial_mutability)) + xlim(0,500)
for(rep in unique(results$replicate_trajectory)[5:10]){
  
  results_subset <- results[results$replicate_trajectory == rep,]
  
  base_pl <- base_pl + geom_line(data=results_subset)
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



pdf('mutability_equilibrium.pdf',height = 7, width = 7)
plot(base_pl)
dev.off()
