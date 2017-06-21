# ================ DEFINING GGPLOT2 PARAMETERS ===============
title_size <- 7
axis_title_size <- 11
y_axis_text_size <- 11
x_axis_text_size <- 11
  
axis_title_size_subplot <- 5
y_axis_text_size_subplot <- 5
x_axis_text_size_subplot <- 5

ylab_distance <- 7
xlab_distance <- 7
plot_margins <- c(0.1, 0.1, 0.1, 0.1)

lineage_names <- c(
  'CH103' = "CH103 (H)",
  'CH103L' = "CH103 (L)",
  'VRC26' = "VRC26 (H)",
  'VRC26L' = "VRC26 (L)",
  'VRC01_13' = 'VRC01-13 (H)',
  'VRC01_01' = 'VRC01-01 (H)',
  'VRC01_19' = 'VRC01-19 (H)'
)
  
ggplot_theme <- theme_bw() +
  theme(axis.title.y = element_text(size = axis_title_size,
                                    margin = margin(0,ylab_distance,0,0)),
        axis.title.x = element_text(size = axis_title_size,
                                    margin = margin(xlab_distance,0,0,0)),
        axis.text.x = element_text(size = x_axis_text_size, vjust = 0.6),
        axis.text.y = element_text(size = y_axis_text_size),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.y = element_line(size = 0.3),
        axis.ticks.length = unit(0.2, "cm"),
        #plot.margin = unit(c(1, 1, 1, 1),"cm"),
        title = element_text(size = title_size),
        strip.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = 'None',
        legend.text=element_text(size=12),
        panel.grid.minor = element_line(size = 0.1),
        panel.grid.major = element_line(size = 0.2)
  )

