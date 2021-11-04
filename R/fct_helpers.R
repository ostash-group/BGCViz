#' rename_vector 
#'
#' @description Function, that given the dataframe, and renaming dataframe, returns renamed vector.
#'
#' @return list of renamed values in a vector + notification
#'
#' @noRd
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
#' correct_width 
#'
#' @description Function, that checks if the width of RRE-Finder, ARTS< SEMPI and PRISM-Supp data should be increases
#' to improve visualization
#'
#' @return data with the corrected width
#'
#' @noRd
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
#' hybrid_col 
#'
#' @description Function, that substitute type with "hybrid", if it contains "__" (therefore cluster have multiple types)
#'
#' @return Vector of substituted types
#'
#' @noRd
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