##' Simulate genetic data and check for steps 
##' 
##' 
##' 
##' @title simulate_data
##' @param step strength of step: 0 (linear) to x (steep step e.g. 50)
##' @param p.start starting allele proportion, 0 to 1
##' @param p.end end allele proportion, 0 to 1
##' @param n.samples number of genomes sampled from a locality 
##' @param n.loci number of loci per genome
##' @param n.pops number of localities sampled along gradient 
##' @param detect_step T/F for detecting steps
##' 
##' Calculated variables
##' d = distance from 0 to 1, increments of increment.size - 1/(n.pops - 1)
##' 
##' 
##' @return output$variables table of inputted variables
##' @return output$data_summary table of simulated data
##' @return output$step_results table of step detection results
##' 
##' @author Alex Sentinella

simulate_data <- function(step = 0, p.start = 0, p.end = 1, n.samples = 20, 
                          n.loci = 1000, n.pops = 10, detect_step = T) {
  
  #Make a table of the variables used for this run
  data_variables <- tibble(step, n.samples, p.start, p.end, n.loci, n.pops)
  
  #Add variables to output
  output <- list(variables = data_variables)
  
  #Make a variable which is the distance between each value of d (distance)
  increment.size = 1/(n.pops - 1)

  #Create a table with input variables and calculate p along distance (d)
  data <- tibble(
    step = step,
    n.samples = n.samples,
    i = increment.size,
    n.pops = n.pops,
    p.start = p.start,
    p.end = p.end,
    d = seq(0, 1, increment.size)) %>% #Location of population along distance (0 to 1)
    mutate(p = qbeta(d, 1 / (1 + step), 1 / (1 + step)) * #true population allele frequency
        (p.end - p.start) + p.start)  #offsets allele proportion from 0 and 1
  
  ## Take random samples from p to get variable allele frequencies
  # Take n.samples from a binomial distribution around allele proportion p
  # Divide by n.samples to get 'measured' allele frequency
  # e.g. n.samples = 5, allele proportion p = 0.1
  # rbinom(1, 5, 0.1) / 5
  # Can only be 0, 0.2, 0.4, 0.6, 0.8, 1
  # But would more likely be 0/0.2  
  # repeat over n.loci
  data <- data %>% 
  .[rep(1:nrow(.), times = n.loci),] %>% #replicate over multiple loci with same p
  rowwise() %>% #allows for mutate to work row by row (rather than as a vector)
  mutate(p.binom = rbinom(1, n.samples, p)/n.samples) %>% #p.binom
  ungroup() #stops rowwise operations
  
  
  ## Alpha diversities

  #Calculate alpha diversities of each p.binom (i.e. sampled allele frequency)
  #Calculated at the locus level
  data <- data %>%
      mutate(
        H0a = get.Hq.alpha(p.binom, 0),
        H1a = get.Hq.alpha(p.binom, 1),
        H2a = get.Hq.alpha(p.binom, 2),
        D0a = H.to.D.alpha(H0a, 0),
        D1a = H.to.D.alpha(H1a, 1),
        D2a = H.to.D.alpha(H2a, 2)
      )

  ## Beta diversities - per locus - AvLast variants
  
  # Calculate beta diversities of each p.binom (i.e. sampled allele frequencies)
  # with the lead p.binom (distance + increment.size)
  #
  # Calculated at the locus level
  #
  # Don't calculate betas comparing distance 1 and 0 (d = 1), return NA instead
  data <- data %>%
    mutate(p.lead = lead(p.binom)) %>% #Add a column containing next p.binom
    mutate(
      H0b.Jac.AvLast  = if_else(d == 1, NA_real_, get.Hq.beta(p.binom, p.lead, 0, q0measure = "Jaccard")),
      H0b.Sor.AvLast  = if_else(d == 1, NA_real_, get.Hq.beta(p.binom, p.lead, 0, q0measure = "Sorenson")),
      H1b.MI.AvLast   = if_else(d == 1, NA_real_, get.Hq.beta(p.binom, p.lead, 1, q1measure = "Mutual Information")),
      H1b.ShD.AvLast  = if_else(d == 1, NA_real_, get.Hq.beta(p.binom, p.lead, 1, q1measure = "Shannon Differentiation")),
      H2b.JOST.AvLast = if_else(d == 1, NA_real_, get.Hq.beta(p.binom, p.lead, 2, q2measure = "Jost-D")),
      H2b.GST.AvLast  = if_else(d == 1, NA_real_, get.Hq.beta(p.binom, p.lead, 2, q2measure = "GST")),
      D0b.A.AvLast    = if_else(d == 1, NA_real_, get.Dq.beta(p.binom, p.lead, 0)),
      D0b.B.AvLast    = H.to.D.alpha(H0b.Jac.AvLast, 0),
      D1b.A.AvLast    = if_else(d == 1, NA_real_, get.Dq.beta(p.binom, p.lead, 1)),
      D1b.B.AvLast    = H.to.D.alpha(H1b.MI.AvLast, 1),
      D2b.A.AvLast    = if_else(d == 1, NA_real_, get.Dq.beta(p.binom, p.lead, 2)),
      D2b.B.AvLast    = H.to.D.alpha(H2b.JOST.AvLast, 2),
      BC.AvLast      = if_else(d == 1, NA_real_, get.BC(p.binom, p.lead)),
      RBC.AvLast      = if_else(d == 1, NA_real_, get.RBC(p.binom, p.lead)),
      H0b.Jac.rel.AvLast  = if_else(d == 1, NA_real_, get.Hq.relative.beta(p.binom, p.lead, 0, q0measure = "Jaccard")),
      H0b.Sor.rel.AvLast  = if_else(d == 1, NA_real_, get.Hq.relative.beta(p.binom, p.lead, 0, q0measure = "Sorenson")),
      H1b.MI.rel.AvLast   = if_else(d == 1, NA_real_, get.Hq.relative.beta(p.binom, p.lead, 1, q1measure = "Mutual Information")),
      H1b.ShD.rel.AvLast  = if_else(d == 1, NA_real_, get.Hq.relative.beta(p.binom, p.lead, 1, q1measure = "Shannon Differentiation")),
      H2b.JOST.rel.AvLast = if_else(d == 1, NA_real_, get.Hq.relative.beta(p.binom, p.lead, 2, q2measure = "Jost-D")),
      H2b.GST.rel.AvLast  = if_else(d == 1, NA_real_, get.Hq.relative.beta(p.binom, p.lead, 2, q2measure = "GST")),
      D0b.A.rel.AvLast    = if_else(d == 1, NA_real_, get.Dq.relative.beta(p.binom, p.lead, 0)),
      D0b.B.rel.AvLast    = H.to.D.alpha(H0b.Jac.rel.AvLast, 0),
      D1b.A.rel.AvLast    = if_else(d == 1, NA_real_, get.Dq.relative.beta(p.binom, p.lead, 1)),
      D1b.B.rel.AvLast    = H.to.D.alpha(H1b.MI.rel.AvLast, 1),
      D2b.A.rel.AvLast    = if_else(d == 1, NA_real_, get.Dq.relative.beta(p.binom, p.lead, 2)),
      D2b.B.rel.AvLast    = H.to.D.alpha(H2b.JOST.rel.AvLast, 2),
    )
    
 ## Beta diversities - AvFirst (average gamma, alphas before calculating beta)
  
 # Calculate beta diversities of each p.binom (i.e. sampled allele frequencies)
 # with the lead p.binom (distance + increment.size)
 #
 # Calculated overall per group (each site)
 #
 #Don't calculate betas comparing distance 1 and 0 (d =1), return NA instead  
 AvFirst_data <- data %>%
   group_by(d) %>%
   filter(d != 1) %>% #Avoid calculating at d = 1
   mutate(
     H0b.Jac.AvFirst_mean =  get.Hq.beta(p.binom, p.lead, 0, per.locus = F, q0measure = "Jaccard"),
     H0b.Sor.AvFirst_mean =  get.Hq.beta(p.binom, p.lead, 0, per.locus = F, q0measure = "Sorenson"),
     H1b.MI.AvFirst_mean =   get.Hq.beta(p.binom, p.lead, 1, per.locus = F, q1measure = "Mutual Information"),
     H1b.ShD.AvFirst_mean =  get.Hq.beta(p.binom, p.lead, 1, per.locus = F, q1measure = "Shannon Differentiation"),
     H2b.JOST.AvFirst_mean = get.Hq.beta(p.binom, p.lead, 2, per.locus = F, q2measure = "Jost-D"),
     H2b.GST.AvFirst_mean =  get.Hq.beta(p.binom, p.lead, 2, per.locus = F, q2measure = "GST"),
     D0b.A.AvFirst_mean =    get.Dq.beta(p.binom, p.lead, 0, per.locus = F),
     D0b.B.AvFirst_mean = H.to.D.alpha(H0b.Jac.AvFirst_mean, 0),
     D1b.A.AvFirst_mean =    get.Dq.beta(p.binom, p.lead, 1, per.locus = F),
     D1b.B.AvFirst_mean = H.to.D.alpha(H1b.MI.AvFirst_mean, 1),
     D2b.A.AvFirst_mean =    get.Dq.beta(p.binom, p.lead, 2, per.locus = F),
     D2b.B.AvFirst_mean = H.to.D.alpha(H2b.JOST.AvFirst_mean, 2)
   ) %>%
   mutate(
     H0b.Jac.AvFirst_sd  = get.Hq.beta.sd(p.binom, p.lead, 0, q0measure = "Jaccard"),
     H0b.Sor.AvFirst_sd  = get.Hq.beta.sd(p.binom, p.lead, 0, q0measure = "Sorenson"),
     H1b.MI.AvFirst_sd   = get.Hq.beta.sd(p.binom, p.lead, 1, q1measure = "Mutual Information"),
     H1b.ShD.AvFirst_sd  = get.Hq.beta.sd(p.binom, p.lead, 1, q1measure = "Shannon Differentiation"),
     H2b.JOST.AvFirst_sd = get.Hq.beta.sd(p.binom, p.lead, 2, q2measure = "Jost-D"),
     H2b.GST.AvFirst_sd  = get.Hq.beta.sd(p.binom, p.lead, 2, q2measure = "GST"),
     D0b.A.AvFirst_sd    = get.Dq.beta.sd(p.binom, p.lead, 0),
     D0b.B.AvFirst_sd    = H0b.Jac.AvFirst_sd, #same sd as H
     D1b.A.AvFirst_sd    = get.Dq.beta.sd(p.binom, p.lead, 1),
     D1b.B.AvFirst_sd    = H1b.MI.AvFirst_sd, #same sd as H
     D2b.A.AvFirst_sd    = get.Dq.beta.sd(p.binom, p.lead, 2),
     D2b.B.AvFirst_sd    = H2b.JOST.AvFirst_sd #same sd as H
   ) %>%
   summarise(across( everything(), mean)) %>%
   select("d", ends_with("_mean"), ends_with("_sd"))
  

  #Create data summary table, 
  #what are the mean and sd of each measure at each distance (d)
  data_summary <- data %>%
    group_by(d) %>%
    summarise_each(list(mean = mean, sd = sd, var = var)) %>%
    left_join(AvFirst_data, by = "d")
  
  #Add data summary to output
  output <- c(output, list(data_summary = data_summary))

  #Check for presence and location of step (optional)
  if(detect_step == T){
 
    #Names of each beta measure
    beta_measure_names <- c("H0b.Jac", "H0b.Sor", "H1b.MI", "H1b.ShD", 
                            "H2b.JOST", "H2b.GST", "D0b.A", "D0b.B", 
                            "D1b.A", "D1b.B", "D2b.A", "D2b.B")

    #Names of each beta measure including their by AvLast and AvFirst variant
    beta_measures <- c(paste0(beta_measure_names, ".AvLast"),
                       paste0(beta_measure_names, ".rel.AvLast"), 
                       paste0(beta_measure_names, ".AvFirst"), 
                       "BC.AvLast", "RBC.AvLast")

    #Detect for a step for each measure, also calculate coefficient of variation
    step_results <- map_dfr(beta_measures, step_check, data = data_summary, n.loci)

    #Add step results to output
    output <- c(output, list(step_results = step_results))
  }

return(output)
}
