##' Checks for step in data by finding 
##' beta values at point:
##' d being significantly lower than d + 1
##' AND
##' d + 1 being significantly higher than d + 2
##' 
##' OR
##' 
##' d being significantly lower than d + 1
##' AND
##' d + 1 being not different to d + 2
##' AND
##' d + 2 being significantly higher than d + 3
##' 
##' 
##' @param measure Diversity measure to check for step
##' @param data one_locus_data
##' @param n.loci Number of loci
##'
##' @return A tibble containing:
##' measure = name of measure
##' step_present = Is the step present (T/F)
##' step_location_start = location of start of step range
##' step_location_end = location of end of step range
##'
##' @title Step check

step_check <- function(measure, data, n.loci) {
  
  #Column containing mean beta values for the measure being tested
  mean_vec <- paste0(measure, "_mean")
  
  #Column containing sd beta values for the measure being tested
  sd_vec <- paste0(measure, "_sd")

  #Run higher_lower_same function (t.test) over each beta diversity value
  beta_diff <- pmap_chr(list(m1 = data[[mean_vec]], m2 = lead(data[[mean_vec]]), 
                             s1 = data[[sd_vec]], s2 = lead(data[[sd_vec]])),
                        higher_lower_same,
                             n1 = n.loci, n2 = n.loci)

  x <- data[c("d", mean_vec, sd_vec)] %>% #select just d, mean, sd of measure
       cbind(beta_diff) %>% #add higher_lower_same results
       mutate(step_location_start = ifelse((lag(beta_diff) == "less" &  #Is there a step at this d
                                   beta_diff == "greater"), 
                                d, #Step between d and d + 1
                         ifelse(lag(beta_diff) == "less" & #Is there a step between this d
                                  beta_diff == "notdifferent" & # and the next d
                                lead(beta_diff) == "greater", 
                                d, #Step between d and d +2
                                NA)), #else NA
              step_location_end = ifelse((lag(beta_diff) == "less" &  #Is there a step at this d
                                              beta_diff == "greater"), 
                                           lead(d), #Step between d and d + 1
                                           ifelse(lag(beta_diff) == "less" & #Is there a step between this d
                                                    beta_diff == "notdifferent" & # and the next d
                                                    lead(beta_diff) == "greater", 
                                                  lead(d, 2), #Step between d and d +2
                                                  NA)))
  
  #Replace any Inf/ -Inf values with NA
  step_location_start <- ifelse(!is.null((na.omit(x$step_location_start))),
                                 (na.omit(x$step_location_start)), NA)
  
  #Replace any Inf/ -Inf values with NA
  step_location_end <- ifelse(!is.null((na.omit(x$step_location_end))),
                                (na.omit(x$step_location_end)), NA)
  
  #T/F column indictating if a step is present
  step_present <- ifelse(!is.na(step_location_start), T, F)
  
  return(tibble(measure = measure, 
                step_present = step_present, 
                step_location_start = step_location_start,
                step_location_end = step_location_end
                ))
}
