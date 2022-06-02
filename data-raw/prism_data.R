## code to prepare `prism_data` dataset goes here
# Function ro parse prism_json
library(magrittr)
source("R/fct_reading_processing.R")
prism_data <- process_prism_json_suppl( rjson::fromJSON(file = "https://raw.githubusercontent.com/pavlohrab/BGCViz/raw_datasets/datasets/sco_prism.json"))[[1]]
usethis::use_data(prism_data, overwrite = TRUE)
