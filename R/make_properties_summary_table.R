##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param merged_results
make_properties_summary_table <- function(merged_results) {

## Rows - Measures
  
## Columns - Properties of measures
 # True Positive Detection - Sum all true positives
 # False Positive Detection - Sum all false positives
 # True Negative Detection - Sum all true negatives
 # Effect of Allele Proportion Position 
 # Effect of Allele proportion difference
 

  #Get full data set of all simulations
  data <- merged_results
  
  #Vector stating the order you want the measures to be displayed
  measure_order <- c("H1b.MI.AvLast", "H1b.MI.AvFirst",
                     "H2b.GST.AvLast", "H2b.GST.AvFirst",
                     "D2b.A.AvFirst","BC.AvLast")
  
  #Filter down to the six best measures
  data <- data %>%
    filter(Measure %in% measure_order)
  
  #Make a column for false positives
  #When step is zero, all steps detected are false
  #Otherwise only those not at 0.5 are false
  data <- data %>%
    mutate(false_pos = ifelse(step == 0, Total.Steps, Incorrect.Steps),
           true_pos = ifelse(step == 0, NA, Correct.Steps),
           true_neg = ifelse(step == 0, 100 - Total.Steps, NA))

  
  
  data_filt <- data %>%
    filter(
      (p.start == 0 & p.end == 1) |
        (p.start == 0.1 & p.end == 0.9) |
        (p.start == 0 & p.end == 0.5) |
        (p.start == 0 & p.end == 0.2) |
        (p.start == 0.3 & p.end == 0.5)
    ) %>%
    group_by(Measure) %>%
    summarise(true_pos_mean = round(mean(true_pos, na.rm = T), digits = 1), 
              false_pos_mean = round(mean(false_pos, na.rm = T), digits = 1), 
              true_neg_mean = round(mean(true_neg, na.rm = T), digits = 1))
  
  
  data_filt <- data_filt %>%
    mutate(Measure = factor(Measure, levels = measure_order)) %>%
    arrange(Measure)
  

  #create basic table
  results_table <- data_filt %>%
    kbl(col.names = c("Measures", "True Positive Detection",
                      "False Positive Detection", 
                      "True Negative Detection")) %>%
    kable_classic(full_width = F, html_font = "Cambria") 
  
  save_kable(results_table, file = "./Outputs/properties_summary_table.html")
  
    return(results_table)
  }
