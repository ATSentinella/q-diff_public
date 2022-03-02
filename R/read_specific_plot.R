##' Get a plot of beta diversities from one_locus_data
##' 
##'
##'
##'
##' @title

read_specific_plot <- function(step, n.samples, p.start, p.end, n.loci, 
                               n.pops, ...) {


plot_index <- detect_index(one_locus_data[c(T,F,F)], ~ 
             ((.x$step == step) &           #options: 0, 1, 5, 10, 50
             (.x$n.samples == n.samples) &  #options: 5, 10, 25
             (.x$p.start == p.start) &      #options: 0, 0.1, 0.3
             (.x$p.end == p.end) &          #options: 0.4, 0.5, 0.9, 1
             (.x$n.loci == n.loci) &        #options: 10, 50, 1000
             (.x$n.pops == n.pops)))        #options: 4, 5, 6, 7, 10, 11


beta_plot <- get.one.locus.plot(readd(one_locus_data, subtargets = plot_index), ...)$beta_plot

return(beta_plot)
}


#read_specific_plot(step = 0, n.samples = 10, p.start = 0.3, 
 #                p.end = 0.4, n.loci = 1000, n.pops = 11, 
  #             div_type = "D", expected = T, errorbars = T, ratios = F)

