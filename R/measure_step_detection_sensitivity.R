##' Takes inputs of: step, n.samples, p.start, p.end, 
##'                  n.samples, n.loci, n.pops, n_reps
##'                  
##' And tallies number of correct and incorrect steps detected for each 
##' diversity measure
##' 
##' If any input is a vector, runs a set of simulations for each combination
##' 
##' and checks for a step 1000 times for each sample size
##' returns data and plot
##' 
##' To allow rerunning of later functions, this is skipped if 
##' "merged_results" already exists
##' 
##' @param step Intensity of step (default = 0). Can be vector of values.
##' @param p.start Starting allele proportion (default = 0). Can be vector of values.
##' @param p.end End allele proportion (default = 1). Can be vector of values.
##' @param n.samples Number of genomes sampled (default = 0). Can be vector of values.
##' @param n.loci Number of loci (default = 1000). Can be vector of values.
##' @param n.pops Number of localities (default = 10). Can be vector of values.
##' @param n_reps Number of times to repeat simulation (default = 10). Can be vector of values.

measure_step_detection_sensitivity <- function(step = 0, p.start = 0, p.end = 1, 
                                               n.samples = 20, n.loci = 1000, 
                                               n.pops = 10, n_reps = 10) {

  #Create a data frame of every possible combination of input variables
  input_combinations <- purrr::cross_df(list(step = step, 
                                             n.samples = n.samples, 
                                             p.start = p.start, 
                                             p.end = p.end, 
                                             n.loci = n.loci,
                                             n.pops = n.pops,
                                             n_reps = n_reps))
  
  
  # Count number of steps detected for each combinations (run n_reps times)
  # Create a table with the model inputs, number of model replicates
  # steps detected for each measure and correct steps detected for each measure
  replication_table_step <- future_pmap_dfr(input_combinations,
                                            count_step_detections) %>%
                            add_column(input_combinations, .before = 1)
  
  return(replication_table_step)

}
