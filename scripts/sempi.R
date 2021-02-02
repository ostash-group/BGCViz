library(RSQLite)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

conn <- dbConnect(RSQLite::SQLite(), args[1])

data <- dbGetQuery(conn, "SELECT * FROM tbl_segments")


data <- data %>%
  filter(trackid==6)

types <- sapply(data$name, function(x){
  tmp <- str_trim(x)
  tmp <- gsub(", ", "", tmp)
  gsub(" ", "__", tmp)
})

sempi_data <- data.frame(cbind(seq(1:length(data$trackid)),data$start, data$end,as.character(types)))
colnames(sempi_data) <- c("Cluster", "Start", "Stop", "Type")
sempi_data$Cluster <- as.numeric(sempi_data$Cluster)
sempi_data$Start <- as.numeric(sempi_data$Start)
sempi_data$Stop <- as.numeric(sempi_data$Stop)
sempi_data$Type <- str_trim(tolower(sempi_data$Type))
write.csv(sempi_data, "SEMPI.csv", row.names = F)