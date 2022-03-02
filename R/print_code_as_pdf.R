##' To Print all of the code if needed
##'
##' @title

print_code_as_pdf <- function() {
  
  lapply(list.files("./R", full.names = TRUE), knitr::stitch)
 
}
