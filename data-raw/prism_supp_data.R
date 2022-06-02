## code to prepare `prism_supp_data` dataset goes here
library(magrittr)
source("R/fct_reading_processing.R")
prism_supp_data <- process_prism_json_suppl( rjson::fromJSON(file = "https://raw.githubusercontent.com/pavlohrab/BGCViz/raw_datasets/datasets/sco_prism.json"))[[2]]
usethis::use_data(prism_supp_data, overwrite = TRUE)
