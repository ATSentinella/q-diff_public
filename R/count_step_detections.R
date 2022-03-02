##' Takes specific inputs of:
##' step, n.samples, p.start, p.end, n.loci, n.pops
##' and counts how many times each of the measures detects 
##' correct/incorrect steps
##' repeats this "n_reps" times
##' 
##'
##' Simulates new data with these variables each time 
##' and returns the sum of correct/incorrect steps across all replicates
##' 
##' functions it calls:
##' simulate_data - main simulation and calculation script

count_step_detections <- function(step, n.samples, p.start, p.end,
                                  n.loci, n.pops, n_reps = 10) {
  
  # Simple function for running a simulation, then extracting the no. of steps
  extract_steps <- function(x) {
    
    sim_results <- simulate_data(step = step, #strength of step(s)
                                 n.samples = n.samples, #number of samples from each site
                                 p.start = p.start, #allele proportion at start
                                 p.end = p.end, #allele proportion at end
                                 n.loci = n.loci, #number of loci
                                 n.pops = n.pops, #number of populations
                                 detect_step = T)
    
    step_results <- sim_results$step_results
    
    total_steps <- as.numeric(step_results[[2]])
    
    correct_steps <- ifelse(step_results[[3]] <= 0.5 & step_results[[4]] >= 0.5, 1, 0)
    
    correct_steps <- ifelse(is.na(correct_steps), 0, correct_steps)
    
    names(total_steps) <- step_results[[1]]
    
    names(correct_steps) <- paste0(step_results[[1]], "_correct")
    
    output <- c(total_steps, correct_steps)
    
    return(output)
  }

  # Add up total number of steps
  number_of_steps <- summarise(map_dfr(1:n_reps, extract_steps), 
                               across(1:76, sum))

  return(number_of_steps)
}





