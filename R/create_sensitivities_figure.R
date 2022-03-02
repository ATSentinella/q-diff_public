##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param variable
##' 
##' @param variable
create_sensitivities_figure <- function(measures = c("H1b.MI.AvLast"), 
                                        colours = c("black"),
                                        p.start = 0,
                                        p.end = 1,
                                        n.samples = 20,
                                        n.loci = 1000,
                                        n.pops = 10) {
    
    #Plot the intensity of a step as Distance (d) vs. Allele Proportion (p)
  vis.steps <- function(step) {
      
      step_plot <- ggplot(data =  simulate_data(step = step,
                                                p.start = p.start, 
                                                p.end = p.end,
                                                n.samples = n.samples,
                                                n.loci = n.loci,
                                                n.pops = n.pops)[[2]],
                          aes(d, p_mean)) +
        geom_point() +
        geom_line() +
        xlab("Distance (d)") +
        ylab("Allele Proportion (p)") +
        theme_classic()
      
      return(step_plot)
    }
    
    #Plot the step sensitivities while varying an input (n.samples, n.loci, n.pops)
  vis.sensitivities <- function(step, variable) {
      
      #Which tested variable is being visualised?
      if (variable == "n.samples") {
        data <- readd(step_sensitivities_samples)
        
        label_for_x <- "Number of samples (n)"
      }
      
      if (variable == "n.loci") {
        data <- readd(step_sensitivities_loci)
        
        label_for_x <- "Number of loci (L)"
      }
      
      if (variable == "n.pops") {
        data <- readd(step_sensitivities_pops)
        
        label_for_x <- "Number of populations (K)"
      }
      
      p <-  ggplot()+
        ylab("Step detections (out of 100)") +
        geom_vline(xintercept = get(variable), 
                   linetype = "dotted") +
        ylim(0,100) +
        theme_classic()
      
      for (i in seq(1, length(measures))) {
        p <- p +
          geom_line(data  = filter(data, step == !!step & p.start == !!p.start  & p.end == !!p.end),
                    aes(x= !!sym(variable), y = !!sym(measures[[i]])), colour = colours[[i]],  linetype = "dashed") +
          xlab(label = label_for_x)
        
        if (step != 0) { #this makes sure no "correct" steps are plotted for step = 0
          p <- p +
            geom_line(data  = filter(data, step == !!step & p.start == !!p.start  & p.end == !!p.end),
                      aes(x = !!sym(variable), y = !!sym(paste0(measures[[i]], "_correct"))), colour = colours[[i]]) +
            xlab(label = label_for_x)
          
        }
      }
      
      return(p)
    }
    
    #Plot an example single simulation with the specified inputs
  vis.betas  <- function(step){
      
      beta_plot <- plot_betas(simulate_data(step = step,
                                            p.start = p.start, 
                                            p.end = p.end,
                                            n.samples = n.samples,
                                            n.loci = n.loci,
                                            n.pops = n.pops),
                              measures = measures,
                              colours = colours,
                              insert = F)[[1]]
      
      return(beta_plot)
    }
    
    #Make a 4 by 3 grid of plots
  step_row <- plot_grid(vis.steps(0), 
                          vis.steps(1), 
                          vis.steps(5), 
                          vis.steps(50),
                          nrow = 1)
    
  samples_row <- plot_grid(
      vis.sensitivities(0, "n.samples"),
      vis.sensitivities(1, "n.samples"),
      vis.sensitivities(5, "n.samples"),
      vis.sensitivities(50, "n.samples"),
      nrow = 1
    )
    
  loci_row <- plot_grid(
      vis.sensitivities(0, "n.loci"),
      vis.sensitivities(1, "n.loci"),
      vis.sensitivities(5, "n.loci"),
      vis.sensitivities(50, "n.loci"),
      nrow = 1
    )
    
  pops_row <- plot_grid(
      vis.sensitivities(0, "n.pops"),
      vis.sensitivities(1, "n.pops"),
      vis.sensitivities(5, "n.pops"),
      vis.sensitivities(50, "n.pops"),
      nrow = 1
    )
    
    
    beta_row <- plot_grid(vis.betas(0),
                          vis.betas(1),
                          vis.betas(5),
                          vis.betas(50),
                          nrow = 1)
    
    title_row <- ggdraw() +
      draw_label(paste0("Measure = ", str_flatten(measures, "_"), 
                        ", p.start = ", p.start,
                        ", p.end = ", p.end), 
                 fontface = 'bold', x = 0, hjust = 0) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
    
   
  final_figure <- plot_grid(title_row, step_row, beta_row,
                            samples_row, loci_row, pops_row,
                            nrow = 6, rel_heights = c(0.3, 1, 1, 1, 1, 1))
  
  folder_name <- ifelse(length(measures) > 1, 
                        "Comparison", 
                        str_flatten(measures, "_"))
  
  dir.create(paste0("./Outputs/", folder_name))
  
  ggsave(final_figure, 
         filename = paste0("./Outputs/", folder_name,
                           "/", folder_name, 
                           "_", p.start, "_", p.end, ".pdf"),
         height = 297,  width = 210, unit = "mm")
  
  
  ggsave(final_figure, 
         filename = paste0("./Outputs/", folder_name,
                           "/", folder_name, 
                           "_", p.start, "_", p.end, ".png"),
         height = 410 *0.8,  width = 250, unit = "mm")
  
  return(final_figure)
}


