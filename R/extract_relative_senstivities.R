##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param merged_results
extract_relative_senstivities <- function(merged_results) {

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
           true_pos = ifelse(step == 0, NA, Correct.Steps))
  
  #Split into individual allele proportion treatments
  data_0_1 <- data %>%
    filter(p.start == 0) %>%
    filter(p.end == 1) %>%
    group_by(Measure, step) %>%
    summarise(true_pos = round(sum(true_pos)/(n()), digits = 1), 
              false_pos = round(sum(false_pos)/(n()), digits = 1)) %>%
    pivot_longer(names_to = "step_type", cols = c(true_pos, false_pos))  %>%
    pivot_wider(names_from = step, values_from = value) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste(.x, "0_1", sep = "_"))
  
  data_0.1_0.9 <- data %>%
    filter(p.start == 0.1) %>%
    filter(p.end == 0.9) %>%
    group_by(Measure, step) %>%
    summarise(true_pos = round(sum(true_pos)/(n()), digits = 1), 
              false_pos = round(sum(false_pos)/(n()), digits = 1)) %>%
    pivot_longer(names_to = "step_type", cols = c(true_pos, false_pos))  %>%
    pivot_wider(names_from = step, values_from = value) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste(.x, "0.1_0.9", sep = "_"))
  
  
  data_0_0.5 <- data %>%
    filter(p.start == 0) %>%
    filter(p.end == 0.5)%>%
    group_by(Measure, step) %>%
    summarise(true_pos = round(sum(true_pos)/(n()), digits = 1), 
              false_pos = round(sum(false_pos)/(n()), digits = 1)) %>%
    pivot_longer(names_to = "step_type", cols = c(true_pos, false_pos))  %>%
    pivot_wider(names_from = step, values_from = value) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste(.x, "0_0.5", sep = "_"))
  
  
  data_0_0.2 <- data %>%
    filter(p.start == 0) %>%
    filter(p.end == 0.2) %>%
    group_by(Measure, step) %>%
    summarise(true_pos = round(sum(true_pos)/(n()), digits = 1), 
              false_pos = round(sum(false_pos)/(n()), digits = 1)) %>%
    pivot_longer(names_to = "step_type", cols = c(true_pos, false_pos))  %>%
    pivot_wider(names_from = step, values_from = value) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste(.x, "0_0.2", sep = "_"))
  
  
  data_0.3_0.5 <- data %>%
    filter(p.start == 0.3) %>%
    filter(p.end == 0.5) %>%
    group_by(Measure, step) %>%
    summarise(true_pos = round(sum(true_pos)/(n()), digits = 1), 
              false_pos = round(sum(false_pos)/(n()), digits = 1)) %>%
    pivot_longer(names_to = "step_type", cols = c(true_pos, false_pos))  %>%
    pivot_wider(names_from = step, values_from = value) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste(.x, "0.3_0.5", sep = "_"))
  
  
  #Put all the data back together
  data_merge <- data_0_1 %>%
  left_join(data_0.1_0.9, by = c("Measure", "step_type"))%>%
  left_join(data_0_0.5, by = c("Measure", "step_type"))%>%
  left_join(data_0_0.2, by = c('Measure', "step_type"))%>%
  left_join(data_0.3_0.5, by = c("Measure", "step_type")) %>%
  rename_with(.cols = -c(Measure, step_type), ~ paste("step", .x, sep = "_"))

  
  
  
  data_merge_format <- data_merge
  
  
  data_merge_format <- data_merge_format %>%
    mutate(Measure = factor(Measure, levels = measure_order)) %>%
    # mutate(Measure = fct_relevel(Measure, measure_order)) %>%
    arrange(Measure)
  
  
  data_merge_format$Measure <- c("MI \n (AvLast)", " ", "MI \n (AvFirst)", " ", 
                                 "GST \n (AvLast)", " ", "GST \n (AvFirst)", " ", 
                                  "D2b.A (AvFirst)", " ", "BC (AvLast)", " ")
  
  data_merge_format <- data_merge_format %>%
  mutate_all(~replace(., is.na(.), "-"))
  
  
  data_merge_format$step_type[data_merge_format$step_type == "true_pos"] <- "True Positives"
  
  
  data_merge_format$step_type[data_merge_format$step_type == "false_pos"] <- "False Positives"
  
  
  
  #correct steps formatting
  data_merge_format[3:22] <- lapply(data_merge_format[3:22], function(x) {
    v <- 1:14
    v[c(1,3,5,7,9,11,13,14)] <- color_tile("transparent", "green")(c(x[c(1,3,5,7,9,11)], 0, 100))
    v[c(2,4,6,8,10,12,13,14)] <- color_tile("transparent", "red")(c(x[c(2,4,6,8,10,12)], 0, 100))
    
    return(v[1:12])                                             
  })
  
  
  
  options(knitr.kable.NA = '-')
  
  
  #create basic table
  results_table <- data_merge_format %>%
    kbl(col.names = c("Measure", "Type of Step", "0", "1", "5", "50", "0", "1", "5", "50",
                 "0", "1", "5", "50", "0", "1", "5", "50", "0", "1", "5", "50"),
        escape = FALSE) %>%
    kable_classic(full_width = F, html_font = "Cambria")  %>%
    #Add heading
    
    #Remove every second Measure name
    
    #Replace Measure name with better formatting
    #Add stub head for each allele treatment
    add_header_above(c(" " = 2, 
                       "Maximal  \n range: \n p = 0 - 1" = 4, 
                       "Maximal range \n without fixation: \n p = 0.1 - 0.9" = 4, 
                       "Halfmaximal  \n range: \n p = 0 - 0.5" = 4,
                       "Narrow range \n near fixation: \n p = 0 - 0.2" = 4,
                       "Narrow range \n far from fixation: \n p = 0.3 - 0.5" = 4)
                     )
    #Make Step names shorter
    #Add conditional formatting
  
  results_table
  
  save_kable(results_table, file = "./Outputs/relative_sensitivities_table.html")
  
  
  return(results_table)
}
