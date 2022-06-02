## code to prepare `arts_data` dataset goes here
arts_data <- utils::read.csv("https://raw.githubusercontent.com/pavlohrab/BGCViz/raw_datasets/datasets/sco_arts.csv")
usethis::use_data(arts_data, overwrite = TRUE)
