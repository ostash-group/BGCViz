## code to prepare `rre_data` dataset goes here
library(magrittr)
Gene.name <- Coordinates <- NULL # Silence R CMD error
rre_data <- utils::read.delim("https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_rre.txt")
# Clean RRE data. Extract coordinates and Locus tag with double underscore delimiter (__)
rre_data <- rre_data %>%
    tidyr::separate(Gene.name, c("Sequence", "Coordinates", "Locus_tag"), sep = "__") %>%
    tidyr::separate(Coordinates, c("Start", "Stop"), sep = "-")
# Add chromosome info column
rre_data$chromosome <- rep("RRE", length(rre_data$Sequence))
# Add ID column
rre_data$ID <- seq(1:length(rre_data$Sequence))
rre_data$Cluster <- rre_data$ID
rre_data <- data.frame(rre_data)
rre_data["Type"] <- "ripp"
rre_data["Type2"] <- "ripp"
rre_data$Start <- as.numeric(rre_data$Start)
rre_data$Stop <- as.numeric(rre_data$Stop)
# Store rre data into local variable
rre_data <- data.frame(rre_data)
usethis::use_data(rre_data, overwrite = TRUE)
