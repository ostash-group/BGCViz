#' sempi_to_csv
#'
#' @description Function, which transforms Track.db file into dataframe, which could be then written to csv.
#' Download project folder from SEMPI and supply as `project_archive` argument to a function
#'
#' @param project_archive - path to project.zip file, downloaded from SEMPI
#' @param write_to - path where to write generated csv file
#'
#' @return csv file in specified location
#' @examples
#' \dontrun{
#' sempi_to_csv(<zip-file>)
#' }
#' @export
sempi_to_csv <- function(project_archive, write_to = getwd()) {
    trackid <- NULL # Silence R CMD note
    utils::unzip(project_archive, files = "genome_browser/main/Tracks.db", exdir = paste0(write_to, "/SEMPI_TracksDB"), junkpaths = T)
    fl <- paste0(stringr::str_extract(write_to, ".*/"), "/SEMPI_TracksDB/Tracks.db")
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)
  
    data <- RSQLite::dbGetQuery(conn, "SELECT * FROM tbl_segments")
    RSQLite::dbDisconnect(conn)
    unlink(paste0(stringr::str_extract(write_to, ".*/"), "/SEMPI_TracksDB"), recursive = T)
    data <- data %>%
        dplyr::filter(trackid == 6)
  
    types <- sapply(data$name, function(x) {
        tmp <- stringr::str_trim(x)
        tmp <- gsub(", ", "", tmp)
        gsub(" ", "__", tmp)
    })
  
    sempi_data <- data.frame(cbind(seq(1:length(data$trackid)), data$start, data$end, as.character(types)))
    colnames(sempi_data) <- c("Cluster", "Start", "Stop", "Type")
    sempi_data$Cluster <- as.numeric(sempi_data$Cluster)
    sempi_data$Start <- as.numeric(sempi_data$Start)
    sempi_data$Stop <- as.numeric(sempi_data$Stop)
    sempi_data$Type <- stringr::str_trim(tolower(sempi_data$Type))
    utils::write.csv(sempi_data, paste0(write_to, "/sempi.csv"), row.names = FALSE)
}

#' prism_to_csv
#'
#' @description Function, that transforms prism json object into dataframe, which could be written to the csv file
#'
#' @param file - path to a json file,
#' @param write_to - path where to write generated csv file
#'
#' @return csv file in specified location
#' @examples
#' \dontrun{
#' prism_to_csv(<json-file>)
#' }
#' @export
prism_to_csv <- function(file, write_to = getwd()) {
    data <- rjson::fromJSON(file = file)
  
  
    types <- sapply(data$prism_results$clusters, function(x) {
       tolower(x$type)
    })
  
    types <- sapply(types, function(x) {
        if (length(unlist(x)) > 1) {
            tmp <- stringr::str_trim(paste0(unlist(x), collapse = "", sep = " "))
            gsub(" ", "__", tmp)
        } else {
           x
        }
    })
  
    start <- sapply(data$prism_results$clusters, function(x) {
       x$start
    })
    end <- sapply(data$prism_results$clusters, function(x) {
        x$end
    })
  
    prism_data <- data.frame(cbind(start, end, types))
    prism_data <- prism_data %>%
       dplyr::transmute(Cluster = as.numeric(rownames(prism_data)), Start = as.numeric(start), Stop = as.numeric(end), Type = types)
    utils::write.csv(prism_data, paste0(write_to, "/prism.csv"), row.names = FALSE)
}

#' antismash_to_csv
#'
#' @description Function, that returns dataframe, out of supplied antismash json file
#'
#' @param file - path to a json file,
#' @param write_to - path where to write generated csv file
#'
#' @return csv file in specified location
#' @examples
#' \dontrun{
#' antismash_to_csv(<json-file>)
#' }
#' @export
antismash_to_csv <- function(file, write_to = getwd()) {
    Start <- Stop <- NULL # To silence R CMD notes
    data <- rjson::fromJSON(file = file)
    types <- sapply(data$records, function(y) {
        lapply(y$features, function(x) {
            if (unlist(x$type == "region")) {
                tolower(x$qualifiers$product)
            }
        })
    })
  
    types <- Filter(Negate(is.null), types)
  
    types <- sapply(types, function(x) {
      if (length(unlist(x)) > 1) {
          tmp <- stringr::str_trim(paste0(unlist(x), collapse = "", sep = " "))
          gsub(" ", "__", tmp)
      } else {
         x
      }
    })
  
    location <- sapply(data$records, function(y) {
      unlist(sapply(y$features, function(x) {
          if (unlist(x$type == "region")) {
             unlist(x$location)
          }
      }))
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
    utils::write.csv(anti_data, paste0(write_to, "/antismash.csv"), row.names = FALSE)
}

#' arts_to_csv
#'
#' @description Function, which extracts tables from arts result zip archive and transforms them into BGCViz input
#'
#' @param project_archive - path to zip file, downloaded from ARTS
#' @param write_to - path where to write generated csv file
#'
#' @return csv file in specified location
#' @examples
#' \dontrun{
#' arts_to_csv(<zip-file>)
#' }
#'
#' @export
arts_to_csv <- function(project_archive, write_to = getwd()) {
    Start <- NULL # Silence R CMD note
    utils::unzip(project_archive, files = c("tables/duptable.tsv", "tables/knownhits.tsv"), exdir = paste0(write_to, "/ARTS_tables"), junkpaths = T)
    known_hits <- utils::read.delim(paste0(stringr::str_extract(write_to, ".*/"), "/ARTS_tables/knownhits.tsv"))
    dupl_table <- utils::read.delim(paste0(stringr::str_extract(write_to, ".*/"), "/ARTS_tables/duptable.tsv"))
    locations <- sapply(known_hits$Sequence.description, function(x) {
        utils::tail(stringr::str_split(x, "\\|")[[1]], 1)
    })
  
    start <- sapply(locations, function(x) {
        stringr::str_split(x, "_")[[1]][1]
    })
    stop <- sapply(locations, function(x) {
        stringr::str_split(x, "_")[[1]][2]
    })
    # Parse known_hits data
    known_table <- data.frame(cbind(start, stop))
    colnames(known_table) <- c("Start", "Stop")
    rownames(known_table) <- seq(1:dim(known_table)[1])
    known_table$Start <- as.numeric(known_table$Start)
    known_table$Stop <- as.numeric(known_table$Stop)
    known_table$Description <- known_hits$Description
    known_table$Model <- known_hits$X.Model
    known_table$Evalue <- known_hits$evalue
    known_table$Bitscore <- known_hits$bitscore
    known_table$ID <- seq(1:dim(known_table)[1])
    known_table$Cluster <- known_table$ID
    known_table$Type <- "resistance"
    known_table$Type2 <- known_table$Type
    known_table$Hit <- NA
    known_table$Core <- "Not_core"
    known_table$Count <- 1
    # Parse duplication data
    get_location_duptable <- function(x, y) {
        test <- stringr::str_split(x, ";")
        test2 <- sub(".*loc\\|", "", test[[1]])
        test3 <- stringr::str_split(test2, " ")
        res <- list()
        for (i in seq(1:length(test3))) {
            id <- paste("hit", as.character(i), sep = "_")
            start <- test3[[i]][1]
            stop <- test3[[i]][2]
            res_1 <- list(id, start, stop)
            res <- append(res, list(res_1))
        }
        return(res)
    }
  
    dup_table <- data.frame()
    for (i in seq(1:dim(dupl_table)[1])) {
        lst <- get_location_duptable(dupl_table$X.Hits_listed.[i])
        fin_data <- data.frame(do.call("rbind", lst))
        fin_data$Core_gene <- dupl_table$X.Core_gene[i]
        fin_data$Description <- dupl_table$Description[i]
        fin_data$Count <- dupl_table$Count[i]
        colnames(fin_data) <- c("Hit", "Start", "Stop", "Core", "Description", "Count")
        dup_table <- rbind(dup_table, fin_data)
    }
    dup_table$Hit <- unlist(dup_table$Hit)
    dup_table$Start <- unlist(dup_table$Start)
    dup_table$Stop <- unlist(dup_table$Stop)
    dup_table$Start <- as.numeric(dup_table$Start)
    dup_table$Stop <- as.numeric(dup_table$Stop)
    dup_table$ID <- seq(1:dim(dup_table)[1])
    dup_table$Cluster <- dup_table$ID
    dup_table$Type <- "core"
    dup_table$Type2 <- dup_table$Type
    dup_table$Evalue <- NA
    dup_table$Bitscore <- NA
    dup_table$Model <- "Core"
    arts_data <- rbind(dup_table, known_table)
    arts_data <- arts_data %>%
      dplyr::arrange(Start)
    arts_data$ID <- seq(1:dim(arts_data)[1])
    arts_data$Cluster <- arts_data$ID
    utils::write.csv(arts_data, paste0(write_to, "/arts.csv"), row.names = FALSE)
}
