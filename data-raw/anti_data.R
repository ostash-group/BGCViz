## code to prepare `anti_data` dataset goes here
anti_data <- utils::read.csv("https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_antismash.csv")
# Add chromosome column
anti_data$chromosome <- rep("A", length(anti_data$Cluster))
# Type magic
anti_data$Type <- stringr::str_trim(tolower(anti_data$Type))
anti_data["Type2"] <- stringr::str_trim(tolower(anti_data$Type))
usethis::use_data(anti_data, overwrite = TRUE)
