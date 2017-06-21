library('ggplot2')
library('coda')
library('grid')

# Load ggplot parameters
source('ggplot_parameters.R')

contrasts_file_path <- '../results/contrasts/observed_lineages/CH103_constant/CH103_con_run1a_mutability_contrasts.csv'

# ====== READ DATAFRAMES WITH MUTABILITY CONTRASTS FOR CH103 HEAVY CHAIN =====
contrasts_dataframe <- read.table(contrasts_file_path, header = T, sep = ',')

# Sample of 100 trees
sample_size <- 100
tree_sample <- sample(unique(contrasts_dataframe$tree), sample_size, replace = F)

pl <- ggplot(contrasts_dataframe, aes(x = S5F_WS_contrast)) +
  xlab('Change in S5F mutability') +
  ylab('Density') +
  ggplot_theme + 
  scale_y_continuous(expand = c(0,0), limits = c(0,80))

for(tree in tree_sample){

  x <- contrasts_dataframe[contrasts_dataframe$tree == tree,
                            'S5F_WS_contrast']
  
  temp_dataframe <- data.frame(x)
  
  #pl <- pl + geom_histogram(data=temp_dataframe,aes(x=x), fill = 'white',
                                      #colour = 'black',alpha = 0.3)
  pl <- pl + geom_density(data=temp_dataframe,aes(x=x), fill = 'white',
                                      colour = 'black',alpha = 0.1, size = 0.1)

}

pl <- pl + geom_vline(xintercept = 0, linetype = 2, size = 0.2)
  
  
# plot means on top of histograms as vertical red lines
# Get distribuiton of the fraction of negative changes
fraction_negative <- c()

  contrast_means <- c()
for(tree in unique(contrasts_dataframe$tree)[1:1000]){
    
  x <- contrasts_dataframe[contrasts_dataframe$tree == tree,
                                'S5F_WS_contrast']

  contrast_means <- c(contrast_means, mean(x))

  if(tree %in% tree_sample){
    mean_value <- mean(x)
    pl <- pl + geom_vline(xintercept = mean_value, colour = 'red',
                                    alpha = 0.3, size = 0.1)
  }
  
  
  # Get overall fraction of negative changes (disregarding 0s)
  if(length(x[x!=0]) > 0){
    fraction_negative <- c(fraction_negative,sum(x[x!=0]<0)/length(x[x!=0]))
  }else{
    fraction_negative <- c(fraction_negative, NA)
  }

}  
  # Get average mean contrast (and HPD)
  means_dat <- data.frame(contrast_means)
  means_HPD_limits <- HPDinterval(as.mcmc(means_dat$contrast_means), 0.95)
  means_llim <- means_HPD_limits[1]
  means_ulim <- means_HPD_limits[2]
  
  temp_dataframe <- data.frame(fraction_negative)
  fraction_HPD_limits <- HPDinterval(as.mcmc(temp_dataframe[, 1]), 0.95)
  fraction_llim <- fraction_HPD_limits[1]
  fraction_ulim <- fraction_HPD_limits[2]
  
  # Add inset plot showing distribuiton of changes that are negative (disregarding 0s):
  dat <- with(density(temp_dataframe$fraction_negative,na.rm=T), data.frame(x, y))
    

  subpl <- ggplot(data=dat, aes(x = x)) + theme_classic() + xlab('Fraction of mutability losses')
  
  limit_factor <- max(abs(0.5 - max(dat$x)),
                      abs(0.5 - min(dat$x)))
  
  subpl <- subpl +
    geom_area(data = dat, mapping = aes(x = ifelse(x>fraction_llim & x< fraction_ulim , x, 0),y), fill = "grey",alpha = 0.5) +
    geom_line(data = dat, mapping = aes(x = x, y = y),size=0.3) +
    scale_y_continuous(expand = c(0,0.01), limits = c(0,1.1*max(dat$y))) +
    geom_vline(xintercept=mean(fraction_negative,na.rm=T), size = 0.3,colour='red') + 
    xlim(c(0.3,0.7)) +
    geom_vline(xintercept = 0.5, linetype=2, size = 0.3) + 
    ylab('Density') + 
    theme(axis.text=element_text(size=x_axis_text_size_subplot),
          axis.title=element_text(size=axis_title_size_subplot, margin = margin(0,0,0,0)),
          axis.line.x = element_line(colour="black", size = 0.3),
          axis.line.y = element_line(colour="black",  size = 0.3),
          plot.title = element_text(size=2)
    )
  
  
# Plotting...
  
vp <- viewport(width = 0.38, height = 0.38, x = 0.185,
               y = unit(8.5, "lines"), just = c("left",
                                               "bottom"))
  
pdf('CH103_overview_S5F_contrasts.pdf', width=3.43, height = 3)
  plot(pl)
  plot(subpl, vp = vp)
dev.off()
  