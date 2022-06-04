#' filter_deepbgc
#'
#' @description Function, that given the deepbgc dataframe and input values as filter, will return filtered dataframe
#'
#' @return Filtered deepbgc dataframe
#'
#' @noRd
filter_deepbgc <- function(deep_data, cluster_type, score_a_input, score_c_input, score_d_input, domains_filter, biodomain_filter, gene_filter) {
    # Silence R CMD note
    alkaloid <- nrps <- other <-
        pks <- ripp <- saccharide <-
        terpene <- score <- Cluster_type <-
        num_domains <- num_bio_domains <-
        num_proteins <- NULL
    score_a <- apply(deep_data %>% dplyr::select(c("antibacterial", "cytotoxic", "inhibitor", "antifungal")), 1, function(x) max(x))
    score_d <- apply(deep_data %>% dplyr::select(c("deepbgc_score")), 1, function(x) max(x))
    score_c <- apply(deep_data %>% dplyr::select(c("alkaloid", "nrps", "other", "pks", "ripp", "saccharide", "terpene")), 1, function(x) max(x))
    if (is.null(cluster_type)) {
        deep_data_chromo <- deep_data %>%
            dplyr::mutate(score = apply(deep_data %>%
                dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene), 1, function(x) max(x)))
        # Cluster_type column. Here extract colnames, and assign max value to a new column
        deep_data_chromo$Cluster_type <- colnames(deep_data_chromo %>% dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene))[apply(deep_data_chromo %>% dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene), 1, which.max)]
        # If max score is under_threshold, print "under_threshold"
        deep_data_chromo <- deep_data_chromo %>%
            dplyr::mutate(Cluster_type = ifelse(score > 50 / 100, Cluster_type, "under_threshold"))
        # Finally store deepbgc data in plotting variable. Do final scores processing
        biocircos_deep <- deep_data_chromo %>%
            dplyr::mutate(product_class = Cluster_type, score_a = score_a, score_d = score_d, score_c = score_c) %>%
            dplyr::filter(
                score_a >= 50 / 100, score_c >= 50 / 100,
                score_d >= 50 / 100, num_domains >= 5,
                num_bio_domains >= 1, num_proteins >= 1
            )
    } else {
        deep_data_chromo <- deep_data %>%
            dplyr::mutate(score = apply(deep_data %>%
                dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene), 1, function(x) max(x)))
        # Cluster_type column. Here extract colnames, and assign max value to a new column
        deep_data_chromo$Cluster_type <- colnames(deep_data_chromo %>% dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene))[apply(deep_data_chromo %>% dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene), 1, which.max)]
        # If max score is under_threshold, print "under_threshold"
        deep_data_chromo <- deep_data_chromo %>%
            dplyr::mutate(Cluster_type = ifelse(score > as.numeric(cluster_type) / 100, Cluster_type, "under_threshold"))
        # Finally store deepbgc data in plotting variable. Do final scores processing
        biocircos_deep <- deep_data_chromo %>%
            dplyr::mutate(product_class = Cluster_type, score_a = score_a, score_d = score_d, score_c = score_c) %>%
            dplyr::filter(
                score_a >= as.numeric(score_a_input) / 100, score_c >= as.numeric(score_c_input) / 100,
                score_d >= as.numeric(score_d_input) / 100, num_domains >= domains_filter,
                num_bio_domains >= biodomain_filter, num_proteins >= gene_filter
            )
    }

    biocircos_deep["Start"] <- biocircos_deep$nucl_start
    biocircos_deep["Stop"] <- biocircos_deep$nucl_end
    biocircos_deep["Type"] <- biocircos_deep$product_class
    biocircos_deep["Type2"] <- biocircos_deep$product_class
    biocircos_deep["Cluster"] <- biocircos_deep$ID
    return(biocircos_deep)
}
#' filter_gecco
#'
#' @description Function, that given the gecco dataframe and input values as filter, will return filtered dataframe
#'
#' @return Filtered gecco dataframe
#'
#' @noRd
filter_gecco <- function(gecco_data, score_cluster_gecco, score_average_gecco, domains_filter_gecco, prot_filter_gecco) {
    # Silence R CMD note
    alkaloid <- nrps <- other <-
        pks <- ripp <- saccharide <-
        terpene <- score <- Type2 <-
        Cluster_type <- score_a <- score_c <-
        num_domains <- num_prot <- NULL
    score_a_gecco <- apply(gecco_data %>% dplyr::select(c("average_p")), 1, function(x) max(x))
    score_c_gecco <- apply(gecco_data %>% dplyr::select(c("alkaloid", "nrps", "other", "pks", "ripp", "saccharide", "terpene")), 1, function(x) max(x))
    if (is.null(score_cluster_gecco)) {
        gecco_data <- gecco_data %>%
            dplyr::mutate(score = apply(gecco_data %>%
                dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene), 1, function(x) max(x))) %>%
            dplyr::mutate(Cluster_type = ifelse(score > 50 / 100, Type2, "under_threshold")) %>%
            dplyr::mutate(Type2 = Cluster_type, score_a = score_a_gecco, score_c = score_c_gecco) %>%
            dplyr::filter(
                score_a >= 50 / 100, score_c >= 50 / 100,
                num_domains >= 1, num_prot >= 1
            )
    } else {
        gecco_data <- gecco_data %>%
            dplyr::mutate(score = apply(gecco_data %>%
                dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene), 1, function(x) max(x))) %>%
            dplyr::mutate(Cluster_type = ifelse(score > as.numeric(score_cluster_gecco) / 100, Type2, "under_threshold")) %>%
            dplyr::mutate(Type2 = Cluster_type, score_a = score_a_gecco, score_c = score_c_gecco) %>%
            dplyr::filter(
                score_a >= as.numeric(score_average_gecco) / 100, score_c >= as.numeric(score_cluster_gecco) / 100,
                num_domains >= domains_filter_gecco, num_prot >= prot_filter_gecco
            )
    }
    return(gecco_data)
}
