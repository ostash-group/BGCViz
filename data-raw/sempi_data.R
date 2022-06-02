## code to prepare `sempi_data` dataset goes here
sempi_data <- utils::read.csv("https://raw.githubusercontent.com/pavlohrab/BGCViz/raw_datasets/datasets/sco_sempi.csv")
sempi_data['Type2'] <- stringr::str_trim(tolower(sempi_data$Type))
usethis::use_data(sempi_data, overwrite = TRUE)
