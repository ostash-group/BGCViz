# Basic check if column exists for undelying functions
check_if_column_exists <- function(data_names,column_name){
  if (column_name %in% stringr::str_to_lower(data_names)){
    return(TRUE)
  } else {
    shiny::showNotification( paste0(column_name, " column does not exist. Data was not integrated into analysis. Please recheck your data and try one more time"),type = "warning")
    return(FALSE)
  }
}
# The functions below are used to validate if all the columns exist,
# that are used for underlying analysis.
# Also basic check of data is numeric is done for respective colunms

# Validation of basic input like PRISM, Antismash and SEMPI
validate_basic_input <- function(data){
  data_names <- names(data)
  if (!(check_if_column_exists(data_names, 'cluster'))){
    shiny::showNotification( paste0("Cluster columns was created on the fly."),type = "message")
    data$Cluster <- seq(1:dim(data)[1])
  }
  if (!(check_if_column_exists(data_names, 'start'))){
    return(FALSE)
  }
  if (!(check_if_column_exists(data_names, 'stop'))){
    return(FALSE)
  }
  if (!(check_if_column_exists(data_names, 'type'))){
    return(FALSE)
  }
  if (length(unique(data$Cluster)) != length(data$Cluster)){
    shiny::showNotification( paste0("Cluster columns contains non unique values. It was regenerated"),type = "message")
    data$Cluster <- seq(1:dim(data)[1])
  }
  if (( T %in% is.na(data$Start)) | (T %in% is.na(data$Stop))){
    shiny::showNotification( paste0(" Start or Stop columns contain missing values. Please fix this and redownload dataframe"),type = "error")
    return(FALSE)
  }
  if ((T %in% is.na(data$Type)) | ("" %in% data$Type)){
    shiny::showNotification( paste0("Type column contain empty data. It was populated with 'unknown' "),type = "warning")
    data$Type[is.na(data$Type)] <- 'unknown'
    data$Type["" %in% data$Type] <- 'unknown'
  }
  if (!(is.numeric(data$Cluster))){
    data$Cluster <- as.numeric(data$Cluster)
  }
  if (!(is.numeric(data$Start))){
    data$Start <- as.numeric(data$Start)
  }
  if (!(is.numeric(data$Stop))){
    data$Stop <- as.numeric(data$Stop)
  }
  if (!(is.character(data$Type))){
    data$Type <- as.character(data$Type)
  }
  return(list(TRUE, data))
}
#Validation of RRE-Finder input
validate_rre_input <- function(data){
  data_names <- names(data)
  if (!(check_if_column_exists(data_names, 'gene.name'))){
    return(FALSE)
  }
  if (F %in% grepl("__", data$Gene.name)){
    return(FALSE)
  }
  if (!(check_if_column_exists(data_names, 'e.value'))){
    return(FALSE)
  } else {
    data$E.value <- as.numeric(data$E.value)
  }
  if (!is.null(data$Probability)){
    if (!(check_if_column_exists(data_names, 'score'))){
      return(FALSE)
    } else{
      data$Score <- as.numeric(data$Score)
    }
    if (!(check_if_column_exists(data_names, 'p.value'))){
      return(FALSE)
    } else {
      data$P.value <- as.numeric(data$P.value)
    }
    if (!(check_if_column_exists(data_names, 'rre.start'))){
      return(FALSE)
    } else {
      data$RRE.start <- as.numeric(data$RRE.start)
    }
    if (!(check_if_column_exists(data_names, 'rre.end'))){
      return(FALSE)
    } else {
      data$RRE.end <- as.numeric(data$RRE.end)
    }
    if (!(check_if_column_exists(data_names, 'probability'))){
      return(FALSE)
    } else {
      data$Probability <- as.numeric(data$Probability)
    }
  }
  return(list(TRUE, data))
}
#Validation of DeepBGC input
validate_deep_input <- function(data){
  data_names <- names(data)
  col_names <- c("nucl_start", "nucl_end","num_proteins", "num_domains", "num_bio_domains","deepbgc_score","antibacterial",
                 "cytotoxic","inhibitor","antifungal","alkaloid","nrp","other","polyketide","ripp","saccharide","terpene",
                 "bgc_candidate_id", "sequence_id")
  num_columns <- c("nucl_start", "nucl_end","num_proteins", "num_domains", "num_bio_domains","deepbgc_score","antibacterial",
                   "cytotoxic","inhibitor","antifungal","alkaloid","nrp","other","polyketide","ripp","saccharide","terpene")
  if (!('cluster' %in% stringr::str_to_lower(data_names))){
    data$Cluster <- seq(1:dim(data)[1])
  }
  for (column_name in col_names){
    if (!(check_if_column_exists(data_names, column_name))){
      return(FALSE)
    }
    if ( T %in% is.na(data[[column_name]])){
      return(FALSE)
    }
    if ( "" %in% data[[column_name]]){
      return(FALSE)
    }
    if (column_name %in% num_columns){
      names(data)[stringi::stri_trans_tolower(names(data)) == column_name] <- column_name
      data[[column_name]] <- as.numeric(data[[column_name]])
    }
  }
  return(list(TRUE, data))
}
#Validation of GECCO input
validate_gecco_input <- function(data){
  data_names <- names(data)
  col_names <- c("start", "end","average_p", "max_p", "type","alkaloid_probability","polyketide_probability",
                 "ripp_probability","saccharide_probability","terpene_probability","nrp_probability","other_probability",
                 "proteins","domains")
  num_columns <- c("start", "end","average_p", "max_p", "alkaloid_probability","polyketide_probability",
                   "ripp_probability","saccharide_probability","terpene_probability","nrp_probability","other_probability")
  if (!('cluster' %in% stringr::str_to_lower(data_names))){
    data$Cluster <- seq(1:dim(data)[1])
  }
  for (column_name in col_names){
    if (!(check_if_column_exists(data_names, column_name))){
      return(FALSE)
    }
    if ( T %in% is.na(data[[column_name]])){
      return(FALSE)
    }
    if ( "" %in% data[[column_name]]){
      return(FALSE)
    }
    if (column_name %in% num_columns){
      names(data)[stringi::stri_trans_tolower(names(data)) == column_name] <- column_name
      data[[column_name]] <- as.numeric(data[[column_name]])
    }
  }
  return(list(TRUE, data))
}