library("dplyr")
## code to prepare RippMining-Genome goes here
ripp_data <- utills::read.table("/home/lev/Uni/genome-finder/RiPP_masterfile.txt", header = FALSE)

colnames(ripp_data) <-c("Cluster", "Type", "Start", "Stop")


#ADDING CHROMOSOME COLUMN
ripp_data$Chromosome <- rep("GF", length(ripp_data$Cluster))
#Type magic
ripp_data$Type <- stringr::str_trim(tolower(ripp_data$Type))
ripp_data["Type2"] <- stringr::str_trim(tolower(ripp_data$Type))
#Mutate NAs
ripp_data <- mutate(ripp_data, Cluster = 1:length(ripp_data$Type))
usethis::use_data(ripp_data, overwrite = TRUE)