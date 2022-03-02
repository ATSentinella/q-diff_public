##' Take one_locus_data and return a plot of beta diversities
##'
##' .. content for \details{} ..
##'
##' @title

plot_betas <- function(one_locus_data, measures, colours,
                               expected = T, errorbars = T, insert = T) {
   #select data 
   data <- one_locus_data$data_summary
   
   #select variables to print at the end
   variables <- one_locus_data$variables
   
   #calculate n for standard error 
   n.se <- variables$n.loci
   
   #create base plot
   p <- ggplot(data) +
     ylab("Adjacent Beta Diversity") +
     xlab("Distance (0 to 1)") +
     xlim(0, 1) +
     theme_classic()
   
   
   if("all.AvLast" %in% measures) {measures <- c(measures, 
                                       "H0b.Jac.AvLast", "H0b.Sor.AvLast",
                                       "H1b.MI.AvLast", "H1b.ShD.AvLast",
                                       "H2b.JOST.AvLast", "H2b.GST.AvLast",
                                       "D0b.A.AvLast", "D0b.B.AvLast",
                                       "D1b.A.AvLast", "D1b.B.AvLast",
                                       "D2b.A.AvLast", "D2b.B.AvLast", "BC.AvLast")}
   
   if("all.H.AvLast" %in% measures) {measures <- c(measures,
                                        "H0b.Jac.AvLast", "H0b.Sor.AvLast",
                                        "H1b.MI.AvLast", "H1b.ShD.AvLast",
                                        "H2b.JOST.AvLast", "H2b.GST.AvLast")}
   
   if("all.D.AvLast" %in% measures) {measures <- c(measures,
                                        "D0b.A.AvLast", "D0b.B.AvLast",
                                        "D1b.A.AvLast", "D1b.B.AvLast",
                                        "D2b.A.AvLast", "D2b.B.AvLast")}
   
   if("all.AvFirst" %in% measures) {measures <- c(measures, 
                                                "H0b.Jac.AvFirst", "H0b.Sor.AvFirst",
                                                "H1b.MI.AvFirst", "H1b.ShD.AvFirst",
                                                "H2b.JOST.AvFirst", "H2b.GST.AvFirst",
                                                "D0b.A.AvFirst", "D0b.B.AvFirst",
                                                "D1b.A.AvFirst", "D1b.B.AvFirst",
                                                "D2b.A.AvFirst", "D2b.B.AvFirst")}
   
   if("all.H.AvFirst" %in% measures) {measures <- c(measures,
                                                  "H0b.Jac.AvFirst", "H0b.Sor.AvFirst",
                                                  "H1b.MI.AvFirst", "H1b.ShD.AvFirst",
                                                  "H2b.JOST.AvFirst", "H2b.GST.AvFirst")}
   
   if("all.D.AvFirst" %in% measures) {measures <- c(measures,
                                                  "D0b.A.AvFirst", "D0b.B.AvFirst",
                                                  "D1b.A.AvFirst", "D1b.B.AvFirst",
                                                  "D2b.A.AvFirst", "D2b.B.AvFirst")}
   
   #Function for adding error bars to plots
   gg_errorbars <- function(measure, colour){
     
      D_adjust <- ifelse(startsWith(measures[i], "D") == T, -1, 0)
      
      
     # Calculate standard errors
     data <- data %>%
       mutate(measure_se = !!sym(paste0(measure, "_sd"))/sqrt(n.se),
              measure_mean = !!sym(paste0(measure, "_mean")))
     
     p <- p +
       geom_errorbar( data = data,
                      aes(
                        x = d + i_mean / 2,
                        y = measure_mean + D_adjust,
                        ymin = ifelse(measure_mean - measure_se + D_adjust < 0, 0, 
                                      measure_mean - measure_se + D_adjust),
                        ymax = measure_mean + measure_se +D_adjust
                      ),
                      colour = colour,
                      width = .05,
                      position = position_dodge(10))
     
     return(p)
   }
   
   
   
   
   for (i in seq(1, length(measures))) {
      
      D_adjust <- ifelse(startsWith(measures[i], "D") == T, -1, 0)

      p <- p + 
         geom_point(aes(d + i_mean / 2, !!sym((paste0(measures[i], "_mean"))) + !!D_adjust), colour = colours[i]) +
         geom_line(aes(d + i_mean / 2, !!sym((paste0(measures[i], "_mean"))) + !!D_adjust), colour = colours[i])
      
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

  return(list(beta_plot = p, variables = variables))
  }
  
  