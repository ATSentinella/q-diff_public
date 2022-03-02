##' For Supplement
##' Creates five figures, one for each of the allele proportion treatments
##' Each compares the best 6 measures overlayed with each other

create_comparison_figures <- function() {
  
  measures = c("D2b.A.AvFirst", "H1b.MI.AvFirst", "H1b.MI.AvLast", 
               "H2b.GST.AvFirst", "H2b.GST.AvLast", "BC.AvLast")
  
  colours = c("purple", "lightblue",  "blue",
              "darkgreen", "green","orange")
  
  fig_0_1 <- create_sensitivities_figure(
    measures = measures,
    colours = colours,
    p.start = 0,
    p.end = 1
  )
  
  fig_0.1_0.9 <- create_sensitivities_figure(
    measures = measures,
    colours = colours,
    p.start = 0.1,
    p.end = 0.9
  )
  
  fig_0_0.5 <- create_sensitivities_figure(
    measures = measures,
    colours = colours,
    p.start = 0,
    p.end = 0.5
  )
  
  
  fig_0_0.2 <- create_sensitivities_figure(
    measures = measures,
    colours = colours,
    p.start = 0,
    p.end = 0.2
  )
  
  fig_0.3_0.5 <- create_sensitivities_figure(
    measures = measures,
    colours = colours,
    p.start = 0.3,
    p.end = 0.5
  )
  
  comparison_figures <- list(fig_0_1 = fig_0_1, 
                            fig_0.1_0.9 = fig_0.1_0.9, 
                            fig_0_0.5 = fig_0_0.5, 
                            fig_0_0.2 = fig_0_0.2,
                            fig_0.3_0.5 = fig_0.3_0.5)

  
  return(comparison_figures)

}
