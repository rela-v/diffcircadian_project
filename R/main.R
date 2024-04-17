# main.R

# Load the required libraries
library(tidyverse)
library(devtools)
dir.create("./site-library", showWarnings = FALSE)
withr::with_libpaths(new = "./site-library/", install_github('diffCircadian/diffCircadian'))

# Load the package
withr::with_libpaths(new = "./site-library/", source("./R/main_analysis.R"))


