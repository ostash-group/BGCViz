## code to prepare `deep_data` dataset goes here
deep_data <- utils::read.delim("https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_deep.tsv")
polyketide <- nrp <- NULL # Silence R CMD error
drop_cols <- c("nrp", "polyketide")
colnames(deep_data) <- stringr::str_to_lower(colnames(deep_data))
# Read data
deep_data <- deep_data %>%
  dplyr::mutate(pks = polyketide, nrps = nrp) %>%
  dplyr::select(-dplyr::one_of(drop_cols))
usethis::use_data(deep_data, overwrite = TRUE)
