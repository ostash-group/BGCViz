library(dplyr)
library(tidyr)
library(stringr)
library(rjson)

args = commandArgs(trailingOnly=TRUE)

data <- fromJSON(file = args[1])
types <- sapply(data$records, function(y){
  lapply(y$features, function(x){
    if (unlist(x$type == 'region')){
      tolower(x$qualifiers$product)
    }
  })
})

types <-  Filter(Negate(is.null), types)

types <- sapply(types, function(x){
  if (length(unlist(x))>1){
    tmp <- str_trim(paste0(unlist(x), collapse = '', sep = " "))
    gsub(" ", "__", tmp)
  }else{
    x
  }
})

location <- sapply(data$records, function(y){
  unlist(sapply(y$features, function(x){
    if (unlist(x$type == 'region')){
      unlist(x$location)
    }
  })
  )
})


location <- gsub("\\[", "", location)
location <- gsub("\\]", "", location)
location <- data.frame(location)
colnames(location) <- "split"
anti_data <- location %>%
  separate(split, c("Start", "Stop")) %>%
  transmute(ID = rownames(location), Start, Stop)

anti_data <- cbind(anti_data, types)
colnames(anti_data) <- c("Cluster", "Start", "Stop", "Type")
anti_data$Cluster <- as.numeric(anti_data$Cluster)
anti_data$Start <- as.numeric(anti_data$Start)
anti_data$Stop <- as.numeric(anti_data$Stop)
write.csv(anti_data, "antismash.csv", row.names = F)