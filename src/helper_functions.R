is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}
# Fix the duplicates in PRISM-Supp data. Therefore 1 row for 1 orf
fix_duplicates <-  function(test_score, order_vec, regul_genes_orfs, test_name){
  dupl_names <- regul_genes_orfs[duplicated(regul_genes_orfs)]
  duplicated_values <- which(duplicated(regul_genes_orfs[order_vec]))
  test_score <- test_score[order_vec]
  to_add <- test_score[(which(duplicated(regul_genes_orfs[order_vec])))]
  test_score <- test_score[-(which(duplicated(regul_genes_orfs[order_vec])))]
  iterate_one_more_time <- c()
  should_iterate = F
  for (i in seq(1:length(test_name))){
    if (length(dupl_names) == 0){
      should_iterate = F
      break
    }
    if (test_name[i]==dupl_names[1]){
      dupl_names = dupl_names[-1]
      test_score[i] = paste0(test_score[i], "/" ,to_add[1])
      to_add = to_add[-1]
      iterate_one_more_time <- c(iterate_one_more_time, i)
    }
  }
  if ((length(iterate_one_more_time)>1) && (length(dupl_names) != 0)){
    should_iterate = T
  }
  while (should_iterate == T) {
    for (i in iterate_one_more_time){
      if (test_name[i] == dupl_names[1]){
        dupl_names = dupl_names[-1]
        test_score[i] = paste0(test_score[i], "/" ,to_add[1])
        to_add = to_add[-1]
      } 
      if (length(dupl_names)==0){
        should_iterate = F
        break
      }
    }
  }
  return(test_score)
}
# PRISM JSON data processing
process_prism_json_suppl <- function(data){
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
  
  
  prism_data <-  data.frame(Cluster=as.numeric(seq(1:length(start))), Start=as.numeric(start), Stop = as.numeric(end), Type = types)
  
  regul_genes_orfs <- sapply(data$prism_results$regulatory_genes, function(x){
    x$orf
  })
  
  names <- sapply(data$prism_results$orfs[[1]]$orfs, function(y){
    y$name
  })
  coordinates <- sapply(data$prism_results$orfs[[1]]$orfs, function(y){
    y$coordinates
  })
  
  test_coords <- as.matrix(coordinates[, names %in% regul_genes_orfs ])
  
  
  reg_genes <-data.frame(t(test_coords))
  colnames(reg_genes) <- c("Start", "Stop")
  reg_genes$Type <- 'regulatory'
  reg_genes$Type2 <- reg_genes$Type
  
  test_name <- names[names %in% regul_genes_orfs]
  ref_names <- test_name
  order_vec <- order(match(regul_genes_orfs, test_name))
  
  
  test_score <- sapply(data$prism_results$regulatory_genes, function(x){
    x$score
  })
  reg_genes$Score  <- fix_duplicates(test_score , order_vec, regul_genes_orfs, ref_names)
  test_name <- sapply(data$prism_results$regulatory_genes, function(x){
    x$name
  })
  reg_genes$Name  <- fix_duplicates(test_name , order_vec, regul_genes_orfs, ref_names)
  test_full_name<- sapply(data$prism_results$regulatory_genes, function(x){
    x$full_name
  })
  reg_genes$Full_name  <- fix_duplicates(test_full_name , order_vec, regul_genes_orfs, ref_names)
  resist_genes_orfs <- sapply(data$prism_results$resistance_genes, function(x){
    x$orf
  })
  
  test_coords_res <- as.matrix(coordinates[, names %in% resist_genes_orfs ])
  
  
  res_genes <-data.frame(t(test_coords_res))
  
  colnames(res_genes) <- c("Start", "Stop")
  res_genes$Type <- 'resistance'
  res_genes$Type2 <- res_genes$Type
  test_name <- names[names %in% resist_genes_orfs]
  order_vec <- order(match(resist_genes_orfs, test_name))
  
  
  test_score <- sapply(data$prism_results$resistance_genes, function(x){
    x$score
  })
  res_genes$Score  <- fix_duplicates(test_score , order_vec, resist_genes_orfs, ref_names)
  test_name <- sapply(data$prism_results$resistance_genes, function(x){
    x$name
  })
  res_genes$Name  <- fix_duplicates(test_name , order_vec, resist_genes_orfs, ref_names)
  test_full_name<- sapply(data$prism_results$resistance_genes, function(x){
    x$full_name
  })
  res_genes$Full_name  <- fix_duplicates(test_full_name , order_vec, resist_genes_orfs, ref_names)
  
  final_reg <- rbind(res_genes, reg_genes) %>% dplyr::arrange(Start)
  final_reg$ID <- seq(1:dim(final_reg)[1])
  final_reg$Cluster <- final_reg$ID
  rownames(final_reg) <- as.numeric(seq(1:dim(final_reg)[1]))
  return(list(prism_data, final_reg))
}
# Filtering the DeepBGC
filter_deepbgc <- function(deep_data,cluster_type,score_a_input,score_c_input,score_d_input,domains_filter,biodomain_filter,gene_filter){
  score_a <- apply(deep_data %>% dplyr::select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
  score_d <- apply(deep_data %>% dplyr::select(c("deepbgc_score")),1, function(x) max(x))
  score_c <- apply(deep_data %>% dplyr::select(c("alkaloid", "nrps","other","pks","ripp","saccharide","terpene")),1, function(x) max(x))
  if (is.null(cluster_type)){
    deep_data_chromo <- deep_data %>%
      dplyr::mutate(score = apply(deep_data %>%
                                    dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, function(x) max(x))) 
    # Cluster_type column. Here extract colnames, and assign max value to a new column
    deep_data_chromo$Cluster_type <- colnames(deep_data_chromo %>% dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene))[apply(deep_data_chromo%>%dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, which.max) ]
    # If max score is under_threshold, print "under_threshold"
    deep_data_chromo <- deep_data_chromo%>%
      dplyr::mutate(Cluster_type = ifelse(score>50/100, Cluster_type, "under_threshold"))
    #Finally store deepbgc data in plotting variable. Do final scores processing 
    biocircos_deep <- deep_data_chromo%>%
      dplyr::mutate( product_class = Cluster_type, score_a = score_a, score_d = score_d, score_c = score_c) %>%
      dplyr::filter(score_a >= 50/ 100, score_c >=50/100 , 
                    score_d >= 50/100,  num_domains >= 5,
                    num_bio_domains>=1, num_proteins>=1)
  } else {
    deep_data_chromo <- deep_data %>%
      dplyr::mutate(score = apply(deep_data %>%
                                    dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, function(x) max(x))) 
    # Cluster_type column. Here extract colnames, and assign max value to a new column
    deep_data_chromo$Cluster_type <- colnames(deep_data_chromo %>% dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene))[apply(deep_data_chromo%>%dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, which.max) ]
    # If max score is under_threshold, print "under_threshold"
    deep_data_chromo <- deep_data_chromo%>%
      dplyr::mutate(Cluster_type = ifelse(score>as.numeric(cluster_type)/100, Cluster_type, "under_threshold"))
    #Finally store deepbgc data in plotting variable. Do final scores processing 
    biocircos_deep <- deep_data_chromo%>%
      dplyr::mutate( product_class = Cluster_type, score_a = score_a, score_d = score_d, score_c = score_c) %>%
      dplyr::filter(score_a >= as.numeric(score_a_input )/ 100, score_c >=as.numeric(score_c_input)/100 , 
                    score_d >= as.numeric(score_d_input)/100,  num_domains >= domains_filter,
                    num_bio_domains>=biodomain_filter, num_proteins>=gene_filter)
  }
  
  biocircos_deep['Start'] <- biocircos_deep$nucl_start
  biocircos_deep['Stop'] <- biocircos_deep$nucl_end
  biocircos_deep['Type'] <- biocircos_deep$product_class
  biocircos_deep['Type2'] <- biocircos_deep$product_class
  biocircos_deep['Cluster'] <- biocircos_deep$ID
  return(biocircos_deep)
}
# Filtering GECCO
filter_gecco <- function(gecco_data,score_cluster_gecco,score_average_gecco,domains_filter_gecco,prot_filter_gecco){
  score_a_gecco <- apply(gecco_data %>% dplyr::select(c("average_p")),1, function(x) max(x))
  score_c_gecco <- apply(gecco_data %>% dplyr::select(c("alkaloid", "nrps","other","pks","ripp","saccharide","terpene")),1, function(x) max(x))
  if (is.null(score_cluster_gecco)){
    gecco_data <- gecco_data %>%
      dplyr::mutate(score = apply(gecco_data %>%
                                    dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, function(x) max(x))) %>%
      dplyr::mutate(Cluster_type = ifelse(score>50/100, Type2, "under_threshold")) %>%
      dplyr::mutate( Type2 = Cluster_type, score_a = score_a_gecco, score_c = score_c_gecco) %>%
      dplyr::filter(score_a >= 50/ 100, score_c >=50/100 ,
                    num_domains >= 1, num_prot>=1)
  } else{
    gecco_data <- gecco_data %>%
      dplyr::mutate(score = apply(gecco_data %>%
                                    dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, function(x) max(x))) %>%
      dplyr::mutate(Cluster_type = ifelse(score>as.numeric(score_cluster_gecco)/100, Type2, "under_threshold")) %>%
      dplyr::mutate( Type2 = Cluster_type, score_a = score_a_gecco, score_c = score_c_gecco) %>%
      dplyr::filter(score_a >= as.numeric(score_average_gecco )/ 100, score_c >=as.numeric(score_cluster_gecco)/100 ,
                    num_domains >= domains_filter_gecco, num_prot>=prot_filter_gecco)
  }
  return(gecco_data)
}
# Renaming the vector for inut$rename event
rename_vector <- function(data, renamed_dataframe, renaming_notification){
  type <- stringr::str_split(data$Type2, "__")
  type_2 <- sapply(type, function(x){
    sapply(x, function(y){
      if (y %in% renamed_dataframe$Code){
        renamed <- as.character(renamed_dataframe$Group[renamed_dataframe$Code == y])
        if( (length(renamed) >1) & (!( as.character(y) %in% names(renaming_notification)))){
          shiny::showNotification(paste("The ", as.character(y), " type have multiple renaming options: ", paste(renamed, collapse = ", ")), 
                                  type = "warning", duration = NULL)
          shiny::showNotification(paste("The  ", renamed[[1]], " was chosen."), type = "warning", duration=NULL)
          renaming_notification[[as.character(y)]] <- renamed[[1]] 
        }
        renamed[[1]]
      } else {
        y
      }
      
    })
    
  })
  type_3 <- sapply(type_2, function(x){
    dupl <-  x[!duplicated(x)]
    paste(dupl, collapse = "__")
  })
  type_4 <- sapply(type_3, function(y){
    if (y %in% as.character(renamed_dataframe$Code)){
      as.character(renamed_dataframe$Group[renamed_dataframe$Code == y])[[1]]
    } else {
      y
    }
  })
  return(list(as.character(type_4), renaming_notification))
}
# Adding the thickness to the visualization for SEMPI, 
# ARTS, RRE, PRISN-Supp data
correct_width <- function(data, label, sempi_width,prism_supp_data_input_width,arts_width,rre_width){
  if ((label == 'SEMPI')&(sempi_width == T)){
    data$Stop <- data$Stop + 30000
  } else if ((label == 'PRISM-Supp')&(prism_supp_data_input_width == T)){
    data$Stop <- data$Stop + 20000
  } else if ((label == 'ARTS')&(arts_width == T)){
    data$Stop <- data$Stop + 30000
  } else if ((label == 'RRE-Finder')&(rre_width == T)){
    data$Stop <- data$Stop + 50000
  }
  return(data)
}

hybrid_col <- function(data){
  data_split <- stringr::str_split(data$Type2, "__")
  types <- sapply(data_split, function(x){
    if (length(unlist(x))>1){
      "hybrid"
    } else{
      x
    }
  })
  return(types)
}