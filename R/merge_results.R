##' Merges the results from
##' step_sensitivities_samples, step_sensitivities_loci, step_sensitivities_pops
##' into a single table - merged_results
##'
##' @title merge_results
##' @param step_sensitivities_samples
##' @param step_sensitivities_loci
##' @param step_sensitivities_pops
##' @return merged_results


merge_results <- function(step_sensitivities_samples,
                          step_sensitivities_loci,
                          step_sensitivities_pops) {

  #Convert to long format
  data_samples_long <- step_sensitivities_samples %>%
    #Take columns with "correct" and make new columns: Measure, Correct.Steps
    pivot_longer(contains("_correct"), 
                 names_to = "Measure.Correct", values_to = "Correct.Steps") %>%
    #Take columns with "coefvar" and make new columns: Measure, coefvar
    pivot_longer(contains("_coefvar"), 
                 names_to = "Measure.coefvar", values_to = "coefvar") %>%
    #Take columns with "correct" and make new columns: Measure, Correct.Steps
    pivot_longer(-c(step:n_reps, Measure.Correct, Correct.Steps, Measure.coefvar, coefvar), 
                 names_to = "Measure.Total", values_to = "Total.Steps") %>%
    #Remove "_correct" from  Measure.Correct column
    mutate(Measure.Correct = sub("_correct", "", Measure.Correct)) %>%
    #Remove "_coefvar" from  Measure.coefvar column
    mutate(Measure.coefvar = sub("_coefvar", "", Measure.coefvar)) %>%
    #just keep rows with matching Measures
    filter(Measure.Correct == Measure.Total) %>%
    #just keep rows with matching Measures
    filter(Measure.Correct == Measure.coefvar) %>%
    #Rename measure column
    rename(Measure = Measure.Correct) %>%
    #Select specific columns for output
    select(c(Measure, step, p.start, p.end, n.samples, 
             n.loci, n.pops, Total.Steps, Correct.Steps, coefvar))
  
  
  data_loci_long <- step_sensitivities_loci %>%
    #Take columns with "correct" and make new columns: Measure, Correct.Steps
    pivot_longer(contains("_correct"), 
                 names_to = "Measure.Correct", values_to = "Correct.Steps") %>%
    #Take columns with "coefvar" and make new columns: Measure, coefvar
    pivot_longer(contains("_coefvar"), 
                 names_to = "Measure.coefvar", values_to = "coefvar") %>%
    #Take columns with "correct" and make new columns: Measure, Correct.Steps
    pivot_longer(-c(step:n_reps, Measure.Correct, Correct.Steps, Measure.coefvar, coefvar), 
                 names_to = "Measure.Total", values_to = "Total.Steps") %>%
    #Remove "_correct" from  Measure.Correct column
    mutate(Measure.Correct = sub("_correct", "", Measure.Correct)) %>%
    #Remove "_coefvar" from  Measure.coefvar column
    mutate(Measure.coefvar = sub("_coefvar", "", Measure.coefvar)) %>%
    #just keep rows with matching Measures
    filter(Measure.Correct == Measure.Total) %>%
    #just keep rows with matching Measures
    filter(Measure.Correct == Measure.coefvar) %>%
    #Rename measure column
    rename(Measure = Measure.Correct) %>%
    #Select specific columns for output
    select(c(Measure, step, p.start, p.end, n.samples, 
             n.loci, n.pops, Total.Steps, Correct.Steps, coefvar))
  
  data_pops_long <- step_sensitivities_pops %>%
    #Take columns with "correct" and make new columns: Measure, Correct.Steps
    pivot_longer(contains("_correct"), 
                 names_to = "Measure.Correct", values_to = "Correct.Steps") %>%
    #Take columns with "coefvar" and make new columns: Measure, coefvar
    pivot_longer(contains("_coefvar"), 
                 names_to = "Measure.coefvar", values_to = "coefvar") %>%
    #Take columns with "correct" and make new columns: Measure, Correct.Steps
    pivot_longer(-c(step:n_reps, Measure.Correct, Correct.Steps, Measure.coefvar, coefvar), 
                 names_to = "Measure.Total", values_to = "Total.Steps") %>%
    #Remove "_correct" from  Measure.Correct column
    mutate(Measure.Correct = sub("_correct", "", Measure.Correct)) %>%
    #Remove "_coefvar" from  Measure.coefvar column
    mutate(Measure.coefvar = sub("_coefvar", "", Measure.coefvar)) %>%
    #just keep rows with matching Measures
    filter(Measure.Correct == Measure.Total) %>%
    #just keep rows with matching Measures
    filter(Measure.Correct == Measure.coefvar) %>%
    #Rename measure column
    rename(Measure = Measure.Correct) %>%
    #Select specific columns for output
    select(c(Measure, step, p.start, p.end, n.samples, 
             n.loci, n.pops, Total.Steps, Correct.Steps, coefvar))

  
  #Merge all data together
  data <- data_samples_long %>%
    full_join(data_loci_long) %>%
    full_join(data_pops_long) %>%
    distinct(across(Measure:n.pops), .keep_all = T)
  
  #Add in incorrect steps column
  merged_results <- data %>%
    mutate(Incorrect.Steps = Total.Steps - Correct.Steps)
  
  saveRDS(merged_results, "Outputs/merged_results")
  
  return(merged_results)
}
