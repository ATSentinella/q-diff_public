##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

generate_sensitivity_figures <- function(step_sensitivities_samples,
                                         step_sensitivities_loci,
                                         step_sensitivities_pops) {
  
  measures <- c("H0b.Jac.AvLast", "H0b.Sor.AvLast",
                  "D0b.A.AvLast", "D0b.B.AvLast",
                  "H0b.Jac.AvFirst", "H0b.Sor.AvFirst",
                  "D0b.A.AvFirst", "D0b.B.AvFirst",
                  "H1b.MI.AvLast", "H1b.ShD.AvLast", 
                  "D1b.A.AvLast", "D1b.B.AvLast",	
                  "H1b.MI.AvFirst", "H1b.ShD.AvFirst", 
                  "D1b.A.AvFirst", "D1b.B.AvFirst",
                  "H2b.GST.AvLast", "H2b.JOST.AvLast", 
                  'D2b.A.AvLast', 'D2b.B.AvLast',
                  "H2b.GST.AvFirst", "H2b.JOST.AvFirst", 
                  'D2b.A.AvFirst', 'D2b.B.AvFirst',
                  "BC.AvLast") #skip relative measures
  
  #Create pdf of results for each measure, for each allele treatment
  figs1 <- map(measures, create_sensitivities_figure, p.start = 0, p.end = 1)
  
  figs2 <- map(measures, create_sensitivities_figure, p.start = 0.1, p.end = 0.9)
  
  figs3 <- map(measures, create_sensitivities_figure, p.start = 0, p.end = 0.5)
  
  figs4 <- map(measures, create_sensitivities_figure, p.start = 0, p.end = 0.2)
  
  figs5 <- map(measures, create_sensitivities_figure, p.start = 0.3, p.end = 0.5)
  
  #After saving the initial files, put them all together in one pdf (130 pages)
  
  #add them all to a list
  l <- c(figs1, figs2, figs3, figs4, figs5)
  
  pdf("Outputs/all_supp_figs.pdf")
  invisible(lapply(l, print))
  dev.off()
  
}
