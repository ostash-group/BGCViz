library(stringr)
library(dplyr)
library(tidyr)
library(rjson)

args = commandArgs(trailingOnly=TRUE)

data <- fromJSON(file =args[1])


types <- sapply(data$prism_results$clusters, function(x){
  tolower(x$type)
})

types <- sapply(types, function(x){
  if (length(unlist(x))>1){
    tmp <- str_trim(paste0(unlist(x), collapse = '', sep = " "))
    gsub(" ", "__", tmp)
  }else{
    x
  }
})

start <- sapply(data$prism_results$clusters, function(x){
  x$start
  
})
end <- sapply(data$prism_results$clusters, function(x){
  x$end
  
})


prism_data <- data.frame(cbind(start, end, types))
prism_data <- prism_data %>%
  transmute(Cluster=as.numeric(rownames(prism_data)), Start=as.numeric(start), Stop = as.numeric(end), Type = types)

write.csv(prism_data, "prism.csv", row.names = F)