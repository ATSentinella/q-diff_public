## These are all the packages used for this project

library(conflicted) #to avoid package conflicts
library(dotenv)
library(drake) #for project management and caching of results
library(tidyverse) #for data wrangling, plotting, iteration (purrr)
library(cowplot) #for insets in ggplot, and merging figures
library(future) #for parallel computing
library(furrr) #for using parallel computing in purrr
library(future.batchtools) #for parallel computing
library(gt)
library(kableExtra)
library(formattable)

#Which packages to prefer when there is a conflict 
conflict_prefer("filter", "dplyr")
conflict_prefer("expand", "tidyr")
conflict_prefer("lag", "dplyr")

