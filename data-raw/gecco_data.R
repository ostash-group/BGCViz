## code to prepare `gecco_data` dataset goes here
gecco_data <- utils::read.delim("https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_gecco.tsv")
# Silence R CMD note
polyketide_probability <- other_probability <-
    nrp_probability <- alkaloid_probability <-
    terpene_probability <- saccharide_probability <-
    ripp_probability <- NULL
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
    "terpene_probability", "nrp_probability", "other_probability"
)
# Read data
gecco_data <- gecco_data %>%
    dplyr::mutate(
        pks = polyketide_probability, other = other_probability, nrps = nrp_probability, alkaloid = alkaloid_probability,
        terpene = terpene_probability, saccharide = saccharide_probability, ripp = ripp_probability
    ) %>%
    dplyr::select(-dplyr::one_of(drop_cols))
gecco_data$num_prot <- sapply(stringr::str_split(as.character(gecco_data$proteins), ";"), length)
gecco_data$num_domains <- sapply(stringr::str_split(as.character(gecco_data$domains), ";"), length)
names(gecco_data)[names(gecco_data) == "start"] <- "Start"
names(gecco_data)[names(gecco_data) == "end"] <- "Stop"
usethis::use_data(gecco_data, overwrite = TRUE)
