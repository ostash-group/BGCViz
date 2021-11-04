#' sempi_to_df 
#'
#' @description Function, which transforms Track.db file into dataframe, which could be then written to csv
#' 
#' @param file - path to a json file,
#' @param write_to - path where to write generated csv file
#'
#' @return csv file in specified location
#'
#' @export
sempi_to_csv <- function(file, write_to = getwd()){
conn <- RSQLite::dbConnect(RSQLite::SQLite(), file)

data <- RSQLite::dbGetQuery(conn, "SELECT * FROM tbl_segments")
RSQLite::dbDisconnect(conn)

data <- data %>%
  dplyr::filter(trackid==6)

types <- sapply(data$name, function(x){
  tmp <- stringr::str_trim(x)
  tmp <- gsub(", ", "", tmp)
  gsub(" ", "__", tmp)
})

sempi_data <- data.frame(cbind(seq(1:length(data$trackid)),data$start, data$end,as.character(types)))
colnames(sempi_data) <- c("Cluster", "Start", "Stop", "Type")
sempi_data$Cluster <- as.numeric(sempi_data$Cluster)
sempi_data$Start <- as.numeric(sempi_data$Start)
sempi_data$Stop <- as.numeric(sempi_data$Stop)
sempi_data$Type <- stringr::str_trim(tolower(sempi_data$Type))
write.csv(sempi_data, paste0(write_to,"/sempi.csv"), row.names = FALSE)
}

#' prism_to_df 
#'
#' @description Function, that transforms prism json object into dataframe, which could be written to the csv file
#' 
#' @param file - path to a json file,
#' @param write_to - path where to write generated csv file
#'
#' @return csv file in specified location
#'
#' @export
prism_to_csv <- function(file, write_to = getwd()){
data <- rjson::fromJSON(file =file)


types <- sapply(data$prism_results$clusters, function(x){
  tolower(x$type)
})

types <- sapply(types, function(x){
  if (length(unlist(x))>1){
    tmp <- stringr::str_trim(paste0(unlist(x), collapse = '', sep = " "))
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
  dplyr::transmute(Cluster=as.numeric(rownames(prism_data)), Start=as.numeric(start), Stop = as.numeric(end), Type = types)
write.csv(prism_data, paste0(write_to,"/prism.csv"), row.names = FALSE)

}

#' antismash_to_df 
#'
#' @description Function, that returns dataframe, out of supplied antismash json file
#' 
#' @param file - path to a json file,
#' @param write_to - path where to write generated csv file
#'
#' @return csv file in specified location
#'
#' @export
antismash_to_csv <- function(file, write_to = getwd()){
  data <- rjson::fromJSON(file = file)
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
      tmp <- stringr::str_trim(paste0(unlist(x), collapse = '', sep = " "))
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
    tidyr::separate(split, c("Start", "Stop")) %>%
    dplyr::transmute(ID = rownames(location), Start, Stop)
  
  anti_data <- cbind(anti_data, types)
  colnames(anti_data) <- c("Cluster", "Start", "Stop", "Type")
  anti_data$Cluster <- as.numeric(anti_data$Cluster)
  anti_data$Start <- as.numeric(anti_data$Start)
  anti_data$Stop <- as.numeric(anti_data$Stop)
  write.csv(anti_data, paste0(write_to,"/antismash.csv"), row.names = FALSE)
}