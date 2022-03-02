#variable - one of: n.samples, n.loci, n.pops

visualise_step_sensitivities <- function(measure, variable, colour) {
  
  if (variable == "n.samples") data <- readd(step_sensitivities_samples)
  
  if (variable == "n.loci") data <- readd(step_sensitivities_loci)
  
  if (variable == "n.pops") data <- readd(step_sensitivities_pops)

  #Creates 4 plots for each value of step (0, 1, 5, 50)
  # as a horizontal grid  
  #Takes inputs:
  # measure - genetic diversity measure results to extract
  # colour - colour of lines to be plotted
  vis.beta.step.grid <- function(variable, p.start, p.end){
    
    #Base properties of each plot
    p <-  ggplot()+
      ylab("Step detections (out of 100)") +
      ylim(0,100) +
      theme_classic() 
    
    #step == 0
    p1 <- p +
      geom_line(data  = filter(data, step == 0 & p.start == !!p.start  & p.end == !!p.end),
                aes(x= !!sym(variable), y = !!sym(measure)), colour = colour,  linetype = "dashed")
    
    #step == 1
    p2 <- p +
      geom_line(data  = filter(data, step == 1 & p.start == !!p.start  & p.end == !!p.end),
                aes(x= !!sym(variable), y = !!sym(measure)), colour = colour,  linetype = "dashed") +
      geom_line(data  = filter(data, step == 1 & p.start == !!p.start  & p.end == !!p.end),
                aes(x = !!sym(variable), y = !!sym(paste0(measure, "_correct"))), colour = colour)
    
    #step == 5
    p3 <- p +
      geom_line(data  = filter(data, step == 5 & p.start == !!p.start  & p.end == !!p.end),
                aes(x= !!sym(variable), y = !!sym(measure)), colour = colour,  linetype = "dashed") +
      geom_line(data  = filter(data, step == 5 & p.start == !!p.start  & p.end == !!p.end),
                aes(x = !!sym(variable), y = !!sym(paste0(measure, "_correct"))), colour = colour)
    
    #step == 50
    p4 <- p +
      geom_line(data  = filter(data, step == 50 & p.start == !!p.start  & p.end == !!p.end),
                aes(x= !!sym(variable), y = !!sym(measure)), colour = colour,  linetype = "dashed") +
      geom_line(data  = filter(data, step == 50 & p.start == !!p.start  & p.end == !!p.end),
                aes(x = !!sym(variable), y = !!sym(paste0(measure, "_correct"))), colour = colour)

    #Merge the four plots together
    p_grid <- plot_grid(p1, p2, p3, p4, ncol = 4,
                        labels = c('Step = 0', 'Step = 1', 'Step = 5', 'Step = 50'), 
                        label_size = 12, hjust = -1, vjust= 5)
    
    
    return(p_grid)
  }
  
  
  
  final_plot <- plot_grid(
    vis.beta.step.grid(variable,  0, 1),
    vis.beta.step.grid(variable,  0.1, 0.9),
    vis.beta.step.grid(variable,  0, 0.5),
    vis.beta.step.grid(variable,  0, 0.2),
    vis.beta.step.grid(variable,  0.3, 0.5),
    ncol = 1, 
    labels = c("0 to 1", "0.1 to 0.9", "0 to 0.5", "0 to 0.2", "0.3 to 0.5"))
  
  
  ggsave(final_plot, 
         filename = paste0("./Outputs/", measure, "_", variable, ".pdf"),
         height = 297,  width = 210, unit = "mm")

}

create_measure_pdfs <- function() {
  
  beta_measure_names <- c("H0b.Jac", "H0b.Sor", "H1b.MI", "H1b.ShD", 
                          "H2b.JOST", "H2b.GST", "D0b.A", "D0b.B", 
                          "D1b.A", "D1b.B", "D2b.A", "D2b.B")
  
  #Names of each beta measure including their by locus and global variant
  beta_measures <- c(paste0(beta_measure_names, ".locus"),
                     paste0(beta_measure_names, ".rel.locus"), 
                     paste0(beta_measure_names, ".global"), 
                     "BC.locus", "RBC.locus")
  
  test_variables <- c("n.samples", "n.loci", "n.pops")
  
  sens_vars <- cross_df(list(beta_measures = beta_measures, 
                             test_variables = test_variables))
  
  map2(sens_vars$beta_measures, sens_vars$test_variables, 
       visualise_step_sensitivities, "black")
  
}