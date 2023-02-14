
#' #' @description A function, that reads RiPPMiner-Genome file txt
#' #'
#' #' @return csv file
#' #'
#' #' @noRd

read_ripp <- function(data) {
    ripp_data <- data
    # CHANGING COLNAMES -- temporary
    colnames(ripp_data) <-c("Cluster", "Type", "Start", "Stop")
    #Validation of input
    res_validation <- validate_basic_input(ripp_data)
    if (!(res_validation[[1]])) {
      ripp_data <- NULL
      return(NULL)
    } else {
      ripp_data <- res_validation[[2]]
    }
    #ADDING CHROMOSOME COLUMN
    ripp_data$Chromosome <- rep("GF", length(ripp_data$Cluster))
    #Type magic
    ripp_data$Type <- stringr::str_trim(tolower(ripp_data$Type))
    ripp_data["Type2"] <- stringr::str_trim(tolower(ripp_data$Type))
    #Mutate NAs
    ripp_data <- mutate(ripp_data, Cluster = 1:length(ripp_data$Type))

  return(ripp_data)
  
}



#' read_anti
#'
#' @description A function, that reads RRE-finder file
#'
#' @return csv file
#'
#' @noRd
read_anti <- function(data) {
    anti_data <- data
    res_validation <- validate_basic_input(anti_data)
    if (!(res_validation[[1]])) {
        anti_data <- NULL
        return(NULL)
    } else {
        anti_data <- res_validation[[2]]
    }
    # Add chromosome column
    anti_data$chromosome <- rep("A", length(anti_data$Cluster))
    # Type magic
    anti_data$Type <- stringr::str_trim(tolower(anti_data$Type))
    anti_data["Type2"] <- stringr::str_trim(tolower(anti_data$Type))
    return(anti_data)
}
#' read_anti
#'
#' @description A function, that reads antismash file
#'
#' @return csv file
#'
#' @noRd
read_gecco <- function(data) {
    # Silence R CMD note
    polyketide_probability <- other_probability <-
        nrp_probability <- alkaloid_probability <-
        terpene_probability <- saccharide_probability <-
        ripp_probability <- NULL
    # Add chromosome column
    gecco_data <- data

    gecco_data$chromosome <- rep("G", length(gecco_data$type))
    # Type magic
    gecco_data$Cluster <- seq(1:length(gecco_data$chromosome))
    gecco_data$ID <- gecco_data$Cluster
    gecco_data$Type <- stringr::str_trim(tolower(gecco_data$type))
    gecco_data$Type <- gsub("polyketide", "pks", gecco_data$Type)
    gecco_data$Type <- gsub("nrp", "nrps", gecco_data$Type)
    gecco_data$Type <- gsub("unknown", "under_threshold", gecco_data$Type)
    gecco_data["Type2"] <- stringr::str_trim(tolower(gecco_data$Type))
    drop_cols <- c(
        "alkaloid_probability", "polyketide_probability", "ripp_probability", "saccharide_probability",
        "terpene_probability", "nrp_probability"
    )
    # Read data
    gecco_data <- gecco_data %>%
        dplyr::mutate(
            pks = polyketide_probability,  nrps = nrp_probability, alkaloid = alkaloid_probability,
            terpene = terpene_probability, saccharide = saccharide_probability, ripp = ripp_probability
        ) %>%
        dplyr::select(-dplyr::one_of(drop_cols))
    gecco_data$num_prot <- sapply(stringr::str_split(as.character(gecco_data$proteins), ";"), length)
    gecco_data$num_domains <- sapply(stringr::str_split(as.character(gecco_data$domains), ";"), length)
    names(gecco_data)[names(gecco_data) == "start"] <- "Start"
    names(gecco_data)[names(gecco_data) == "end"] <- "Stop"
    return(gecco_data)
}
read_prism <- function(data, json = TRUE) {
    if (json == TRUE) {
        processed_data <- process_prism_json_suppl(data)
        prism_data <- processed_data[[1]]
        prism_supp_data <- processed_data[[2]]
    } else {
        prism_data <- data
        prism_supp_data <- NULL
    }
    res_validation <- validate_basic_input(prism_data)
    if (!(res_validation[[1]])) {
        prism_data <- NULL
        return(NULL)
    } else {
        prism_data <- res_validation[[2]]
    }
    prism_data$Type <- stringr::str_trim(tolower(prism_data$Type))
    prism_data["Type2"] <- stringr::str_trim(tolower(prism_data$Type))
    return(list(prism_data, prism_supp_data))
}
read_sempi <- function(data, zip = TRUE) {
    # Silence R CMD note
    trackid <- NULL
    if (zip == TRUE) {
        utils::unzip(data, files = "genome_browser/main/Tracks.db", exdir = "./SEMPI_TracksDB", junkpaths = TRUE)
        fl <- "./SEMPI_TracksDB/Tracks.db"
        conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

        data <- RSQLite::dbGetQuery(conn, "SELECT * FROM tbl_segments")
        RSQLite::dbDisconnect(conn)
        unlink("./SEMPI_TracksDB", recursive = TRUE)
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
    } else {
        sempi_data <- data
    }
    res_validation <- validate_basic_input(sempi_data)
    if (!(res_validation[[1]])) {
        sempi_data <- NULL
        return(NULL)
    } else {
        sempi_data <- res_validation[[2]]
    }
    sempi_data["Type2"] <- stringr::str_trim(tolower(sempi_data$Type))
    return(sempi_data)
}
read_arts_archive <- function(archive, zip = TRUE) {
    # Silence R CMD note
    Start <- Core <- NULL
    if (zip == TRUE) {
        utils::unzip(archive, files = c("tables/duptable.tsv", "tables/knownhits.tsv"), exdir = "./ARTS_tables", junkpaths = TRUE)
        known_hits <- utils::read.delim("./ARTS_tables/knownhits.tsv")
        dupl_table <- utils::read.delim("./ARTS_tables/duptable.tsv")
        unlink("./ARTS_tables", recursive = TRUE)
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
    } else {
        arts_data <- archive
    }
    return(arts_data)
}
read_deep <- function(data) {
    polyketide <- nrp <- NULL # Silence R CMD error
    # Fix colnames in deepbgc data
    colnames(data) <- stringr::str_to_lower(colnames(data))
    res_validation <- validate_deep_input(data)
    if (!(res_validation[[1]])) {
        deep_data <- NULL
        return(NULL)
    } else {
        deep_data <- res_validation[[2]]
    }
    drop_cols <- c("nrp", "polyketide")
    # Read data
    deep_data <- deep_data %>%
        dplyr::mutate(pks = polyketide, nrps = nrp) %>%
        dplyr::select(-dplyr::one_of(drop_cols))
    return(deep_data)
}
read_rre <- function(data) {
    Gene.name <- Coordinates <- NULL # Silence R CMD error
    res_validation <- validate_rre_input(data)
    if (!(res_validation[[1]])) {
        data <- NULL
        return(NULL)
    } else {
        data <- res_validation[[2]]
    }
    # Clean RRE data. Extract coordinates and Locus tag with double underscore delimiter (__)
    rre_data <- data %>%
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
}
