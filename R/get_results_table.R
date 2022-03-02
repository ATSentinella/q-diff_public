##' Takes merged_results
##' and creates Table 2
##' results under standard conditions

get_results_table <- function(merged_results){
  
  #Get full data set of all simulations
  data <- merged_results
  
  #Vector stating the order you want the measures to be displayed
  measure_order <-  c("H0b.Jac.AvLast", "H0b.Sor.AvLast",
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
                      "BC.AvLast")
  
  #Filter data down to default variables
  data <- data %>%
    filter(Measure %in% measure_order,
           n.samples == 20,
           n.loci == 1000,
           n.pops == 10)
  
  #Make a column for false positives
  #When step is zero, all steps detected are false
  #Otherwise only those not at 0.5 are false
  data <- data %>%
    mutate(false_pos = ifelse(step == 0, Total.Steps, Incorrect.Steps),
           Correct.Steps = ifelse(step == 0, NA, Correct.Steps))
  
  #Split into individual allele proportion treatments
  data_0_1 <- data %>%
    filter(p.start == 0) %>%
    filter(p.end == 1) %>%
    group_by(Measure, step) %>%
    summarise(correct = round(sum(Correct.Steps)/(n()), digits = 1), 
              false_pos = round(sum(false_pos)/(n()), digits = 1)) %>%
    pivot_longer(names_to = "step_type", cols = c(correct, false_pos))  %>%
    pivot_wider(names_from = step, values_from = value) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste(.x, "0_1", sep = "_"))
  
  data_0.1_0.9 <- data %>%
    filter(p.start == 0.1) %>%
    filter(p.end == 0.9) %>%
    group_by(Measure, step) %>%
    summarise(correct = round(sum(Correct.Steps)/(n()), digits = 1), 
              false_pos = round(sum(false_pos)/(n()), digits = 1)) %>%
    pivot_longer(names_to = "step_type", cols = c(correct, false_pos))  %>%
    pivot_wider(names_from = step, values_from = value) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste(.x, "0.1_0.9", sep = "_"))
  
  
  data_0_0.5 <- data %>%
    filter(p.start == 0) %>%
    filter(p.end == 0.5)%>%
    group_by(Measure, step) %>%
    summarise(correct = round(sum(Correct.Steps)/(n()), digits = 1), 
              false_pos = round(sum(false_pos)/(n()), digits = 1)) %>%
    pivot_longer(names_to = "step_type", cols = c(correct, false_pos))  %>%
    pivot_wider(names_from = step, values_from = value) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste(.x, "0_0.5", sep = "_"))
  
  
  data_0_0.2 <- data %>%
    filter(p.start == 0) %>%
    filter(p.end == 0.2) %>%
    group_by(Measure, step) %>%
    summarise(correct = round(sum(Correct.Steps)/(n()), digits = 1), 
              false_pos = round(sum(false_pos)/(n()), digits = 1)) %>%
    pivot_longer(names_to = "step_type", cols = c(correct, false_pos))  %>%
    pivot_wider(names_from = step, values_from = value) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste(.x, "0_0.2", sep = "_"))
  
  
  data_0.3_0.5 <- data %>%
    filter(p.start == 0.3) %>%
    filter(p.end == 0.5) %>%
    group_by(Measure, step) %>%
    summarise(correct = round(sum(Correct.Steps)/(n()), digits = 1), 
              false_pos = round(sum(false_pos)/(n()), digits = 1)) %>%
    pivot_longer(names_to = "step_type", cols = c(correct, false_pos))  %>%
    pivot_wider(names_from = step, values_from = value) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste(.x, "0.3_0.5", sep = "_"))
  
  
  #Put all the data back together
  data_merge <- data_0_1 %>%
    left_join(data_0.1_0.9, by = c("Measure", "step_type"))%>%
    left_join(data_0_0.5, by = c("Measure", "step_type"))%>%
    left_join(data_0_0.2, by = c('Measure', "step_type"))%>%
    left_join(data_0.3_0.5, by = c("Measure", "step_type")) %>%
    rename_with(.cols = -c(Measure, step_type), ~ paste("step", .x, sep = "_"))
  
  #Make table in correct order
  data_merge <- data_merge %>%
    mutate(Measure = factor(Measure, levels = measure_order)) %>%
    arrange(Measure) 
  
  #Format NAs as "-"
  data_merge <- data_merge %>%
    mutate_all(~replace(., is.na(.), "-"))
  
  # Change Labels for T/F positives
  data_merge$step_type[data_merge$step_type == "correct"] <- "True Positives"
  data_merge$step_type[data_merge$step_type == "false_pos"] <- "False Positives"
  
  
  #T/F positive steps formatting
  data_merge[3:22] <- lapply(data_merge[3:22], function(x) {
    
    v <- 1:(length(x) + 4)
    
    v[c(T, F)] <- color_tile("transparent", "lightblue")(c(x[c(T, F)], 0, 100))
    
    v[c(F, T)] <- color_tile("transparent", "tomato")(c(x[c(F, T)], 0, 100))
    
    return(v[1:length(x)])                                             
  })
  
  #Make blank measure rows for easier reading
  data_merge$Measure <-c("H0b Jac. (AvLast)", " ", "H0b Sor. (AvLast)"," ", 
                                "D0b.A (AvLast)"," ",  "D0b.B (AvLast)"," ", 
                                "H0b Jac. (AvFirst)"," ",  "H0b Sor. (AvFirst)"," ", 
                                "D0b.A (AvFirst)"," ",  "D0b.B (AvFirst)"," ", 
                                "H1b MI (AvLast)"," ",  "H1b ShD (AvLast)"," ",  
                                "D1b.A (AvLast)"," ",  "D1b.B (AvLast)"," ", 	
                                "H1b MI (AvFirst)"," ",  "H1b ShD (AvFirst)"," ",  
                                "D1b.A (AvFirst)"," ",  "D1b.B (AvFirst)"," ", 
                                "H2b GST (AvLast)"," ",  "H2b Jost-D (AvLast)"," ",  
                                'D2b.A (AvLast)'," ",  'D2b.B (AvLast)'," ", 
                                "H2b GST (AvFirst)"," ",  "H2b Jost-D (AvFirst)"," ",  
                                'D2b.A (AvFirst)'," ",  'D2b.B (AvFirst)'," ", 
                                "BC (AvLast)", " ")
  
  
  #create one big table 
  results_table <- data_merge %>%
    ungroup %>%
    kbl(col.names = c("Measure", "Type of Step", "0", "1", "5", "50", "0", "1", "5", "50",
                      "0", "1", "5", "50", "0", "1", "5", "50", "0", "1", "5", "50"),
        escape = FALSE) %>%
    kable_classic(full_width = F, html_font = "Cambria")  %>%
    #Add stub head for each allele treatment
    add_header_above(c(" " = 2, 
                       "Maximal  \n range: \n p = 0 - 1" = 4, 
                       "Maximal range \n without fixation: \n p = 0.1 - 0.9" = 4, 
                       "Halfmaximal  \n range: \n p = 0 - 0.5" = 4,
                       "Narrow range \n near fixation: \n p = 0 - 0.2" = 4,
                       "Narrow range \n far from fixation: \n p = 0.3 - 0.5" = 4)
    ) %>%
    #Add conditional formatting
    pack_rows("q = 0 Measures", 1, 16)%>%
    pack_rows("q = 1 Measures", 17, 32) %>%
    pack_rows("q = 2 Measures", 33, 48)%>%
    pack_rows("Bray-Curtis", 49, 50)
  
  save_kable(results_table, file = "./Outputs/standard_treatment_table.html")
  
  return(results_table)
}
