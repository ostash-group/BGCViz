## code to prepare `ripp_data` dataset goes here
ripp_data <- utils::read.table("https://raw.githubusercontent.com/2061Tsarin/BGCViz-datasets/main/example_data/sco_ripp.txt")

colnames(ripp_data) <-c("Cluster", "Type", "Start", "Stop")

#ADDING CHROMOSOME COLUMN
ripp_data$chromosome <- rep("GF", length(ripp_data$Cluster))
#Type magic
ripp_data$Type <- stringr::str_trim(tolower(ripp_data$Type))
ripp_data["Type2"] <- stringr::str_trim(tolower(ripp_data$Type))
#Mutate NAs
ripp_data <- dplyr::mutate(ripp_data, Cluster = 1:length(ripp_data$Type))


# usethis::use_data(ripp_data, overwrite = TRUE)
## Look at use_data_internally.R file!



