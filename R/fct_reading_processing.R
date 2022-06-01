#' fix_duplicates 
#'
#' @description Function, that fix duplicates in Prism_supp data. (Is certain orf have 2 hits to a model in a data, then 
#' they will be combined under one orf, but not the 2 separate records)
#'
#' @return Vector without duplicates
#'
#' @noRd
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
#' process_prism_json_suppl 
#'
#' @description Function, that processes the json prism file.
#'
#' @return list of prism data and prism_supplement data
#'
#' @noRd
process_prism_json_suppl <- function(data){
  Start <- NULL # Silence R CMD note
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