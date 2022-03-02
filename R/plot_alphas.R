
plot_alphas <- function(one_locus_data, measures, colours,
                       expected = T, errorbars = T, insert = F) {
  #select data 
  data <- one_locus_data$data_summary
  
  #select variables to print at the end
  variables <- one_locus_data$variables
  
  #calculate n for standard error 
  n.se <- variables$n.loci
  
  #create base plot
  p <- ggplot(data) +
    ylab("Alpha Diversity") +
    xlab("Distance (0 to 1)") +
    xlim(0, 1) +
    theme_classic()
  
   #Function for adding error bars to plots
  gg_errorbars <- function(measure, colour){
    
    # Calculate standard errors
    data <- data %>%
      mutate(measure_se = !!sym(paste0(measure, "_sd"))/sqrt(n.se),
             measure_mean = !!sym(paste0(measure, "_mean")))
    
    p <- p +
      geom_errorbar( data = data,
                     aes(
                       x = d,
                       y = measure_mean,
                       ymin = ifelse(measure_mean - measure_se < 0, 0, 
                                     measure_mean - measure_se),
                       ymax = measure_mean + measure_se
                     ),
                     colour = colour,
                     width = .05,
                     position = position_dodge(10))
    
    return(p)
  }
  
  
  for (i in seq(1, length(measures))) {
    
    p <- p + 
      geom_point(aes(d, !!sym((paste0(measures[i], "_mean")))), colour = colours[i]) +
      geom_line(aes(d, !!sym((paste0(measures[i], "_mean")))), colour = colours[i])
    
    if(errorbars == T)  p <- gg_errorbars(measures[i], colours[i])
    
  }
  
  #Toggle of an iset plot of the underlying allele proportions
  if(insert == T){
    
    inset.plot <- ggplot(data) +
      geom_point(aes(d, p_mean))+
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    
    p <- ggdraw(p) +
      draw_plot(inset.plot, x = 0.7, y = .7, width = .3, height = .3)
    
  }
  
  return(alpha_plot = p)
}

