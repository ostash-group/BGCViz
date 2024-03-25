## code to prepare `sempi_data` dataset goes here
sempi_data <- utils::read.csv("https://raw.githubusercontent.com/ostash-group/BGCViz-datasets/main/example_data/sco_sempi.csv")
sempi_data["Type2"] <- stringr::str_trim(tolower(sempi_data$Type))
# usethis::use_data(sempi_data, overwrite = TRUE)
## Look at use_data_internally.R file!

