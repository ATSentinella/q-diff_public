##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

vis_alpha_peaks <- function() {

  alpha_peaks <- plot_grid(plot_alphas(simulate_data(step =0,),"H0a", "red"),
                               plot_alphas(simulate_data(step =50),"H0a", "red"),
                               plot_alphas(simulate_data(step =0),"H1a", "blue"),
                               plot_alphas(simulate_data(step =50),"H1a", "blue"),
                               plot_alphas(simulate_data(step =0),"H2a", "darkgreen"),
                               plot_alphas(simulate_data(step =50),"H2a", "darkgreen"),
                               nrow = 3,
                               labels = c("A", "B", "C", "D", "E", "F"))
  
  return(alpha_peaks)

  }
