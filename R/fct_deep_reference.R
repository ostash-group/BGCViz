#' simple_seg 
#'
#' @description Function, that creates simple dataframe for values to plot in deep_reference. 
#'
#' @return dataframe, with the fields: x (start of the cluster), y (arbitrary letter to change to the software further),
#' xend (End of a cluster), yend (arvitrary letter), Type (permanent), Type2 (which is temporary, if the data is renamed, 
#' used for coloring and legend), Software, ID (unique number), Start, Stop (of a cluster, to show on mouse hover)
#'
#' @noRd
simple_seg <- function(df, letter, software, soft_name ,soft_namings, inter=T, inters){
  if (inter== T){
    data <- df[df$Cluster %in% inters[[soft_namings]][[soft_name]]$from, ]
  } else{
    data <- df
  }
  
  seg_df <-  data.frame(x=as.numeric(data$Start),
                        y=rep(letter, length(data$Cluster)),
                        xend=as.numeric(  data$Stop),
                        yend=rep(letter, length(data$Cluster)),
                        Type = as.factor(data$Type),
                        Type2 = as.factor(data$Type2),
                        Software = rep(software, length(data$Cluster)),
                        ID = data$Cluster,
                        Start = data$Start,
                        Stop = data$Stop)
  return(seg_df)
}


#' add_arts 
#'
#' @description Function, that add arts-specific data to the simple dataframe
#'
#' @return dataframe with the fields, specified in simple_seg() + added Hit, Core, Count, E_value, Bitscore, Model. 
#'
#' @noRd
add_arts <- function(seg_df, soft_namings, df, inter=T,inters){
  if (inter == T){
    subset_df <- df[df$Cluster %in% inters[[soft_namings]]$arts$from, ]
  }else{
    subset_df <- df
  }
  seg_df$Hit = subset_df$Hit
  seg_df$xend = as.numeric(subset_df$Stop)
  seg_df$Core = subset_df$Core
  seg_df$Count = subset_df$Count
  seg_df$E_value = subset_df$Evalue
  seg_df$Bitscore = subset_df$Bitscore
  seg_df$Model = subset_df$Model
  return(seg_df)
}
#' add_prism_supp 
#'
#' @description  Function, that add prism-supplement-specific data to the simple dataframe
#'
#' @return dataframe with the fields, specified in simple_seg() + added Score, Name and Full_Name
#'
#' @noRd
add_prism_supp <- function(seg_df, soft_namings, df, inter=T,inters){
  if (inter == T){
    subset_df <- df[df$Cluster %in% inters[[soft_namings]]$prism_supp$from, ]
  }else{
    subset_df <- df
  }
  seg_df$xend <-  as.numeric(subset_df$Stop)
  seg_df$Score = subset_df$Score
  seg_df$Name = subset_df$Name
  seg_df$Full_name = subset_df$Full_name
  return(seg_df)
}
#' add_deep 
#'
#' @description  Function, that add deepbgc-specific data to the simple dataframe
#'
#' @return dataframe with the fields, specified in simpl_seg() + added Num_domains, deepbgc_score, activity
#'
#' @noRd
add_deep <- function(seg_df, soft_namings, df, inter=T,inters){
  if (inter == T){
    subset_df <- df[df$Cluster %in% inters[[soft_namings]]$deep$from, ]
  }else{
    subset_df <- df
  }
  seg_df$num_domains = subset_df$num_domains
  seg_df$deepbgc_score = subset_df$deepbgc_score
  seg_df$activity = subset_df$product_activity
  return(seg_df)
}
#' add_rre 
#'
#' @description  Function, that add rre-finder-specific data to the simple dataframe. The function also checks if RRE
#' data is in the long format, or in the short one. It depends on RRE-Finder run parameters
#'
#' @return dataframe with the fields, specified in simple_seg() + added E_values in short format, or + Score, E_value, 
#' P_value, RRE_Start, RRE_stop and Probability in long format
#'
#' @noRd
add_rre <- function(seg_df, soft_namings, df, inter=T, rre_more,inters){
  if (inter == T){
    subset_df <- df[df$Cluster %in% inters[[soft_namings]]$rre$from, ]
  }else{
    subset_df <- df
  }
  if (rre_more == T){
    seg_df$xend=as.numeric(subset_df$Stop)
    seg_df$Score = subset_df$Score
    seg_df$Stop = subset_df$Stop
    seg_df$E_value = subset_df$E.value
    seg_df$P_value = subset_df$P.value
    seg_df$RRE_start = subset_df$RRE.start
    seg_df$RRE_stop = subset_df$RRE.end
    seg_df$Probability = subset_df$Probability
  } else {
    seg_df$xend=subset_df$Stop
    seg_df$E_value = subset_df$E.value
  }
  
  return(seg_df)
}
#' add_gecco 
#'
#' @description  Function, that add gecco-specific data to the simple dataframe. 
#'
#' @return dataframe with the fields, specified in simple_seg() + Num_proteins, Num_domains, Average_p, Max_p
#'
#' @noRd
add_gecco <- function(seg_df, soft_namings, df, inter=T,inters){
  if (inter == T){
    subset_df <- df[df$Cluster %in% inters[[soft_namings]]$gecco$from, ]
  }else{
    subset_df <- df
  }
  seg_df$Num_proteins = subset_df$num_prot
  seg_df$Num_domains = subset_df$num_domains
  seg_df$Average_p = subset_df$average_p
  seg_df$Max_p = subset_df$max_p
  return(seg_df)
}

#' define_spec_seg_df 
#'
#' @description Master function, which invokes the appropriate add_ function, based on a label.
#'
#' @return dataframe with the fields, specified in simple_seg() + specific to the software.
#'
#' @noRd
define_spec_seg_df <- function(soft_names, index,seg_df, soft_major, df , inter=T, rre_more,inters){
  if (inter == F){
    soft_major <- "Not applicable"
  }
  if ((soft_names[index] == "prism_supp") & (soft_names[index] != soft_major)){
    seg_df <-  add_prism_supp(seg_df, soft_major, df, inter,inters)
  } else if ((soft_names[index] == "arts")& (soft_names[index] != soft_major)){
    seg_df <- add_arts(seg_df, soft_major, df, inter,inters)
  } else if ((soft_names[index] == "deep")& (soft_names[index] != soft_major)){
    seg_df <- add_deep(seg_df, soft_major, df, inter,inters)
  } else if ((soft_names[index] == "gecco")& (soft_names[index] != soft_major)){
    seg_df <- add_gecco(seg_df, soft_major, df, inter,inters)
  } else if ((soft_names[index] == "rre")& (soft_names[index] != soft_major)){
    seg_df <- add_rre(seg_df, soft_major, df, inter, rre_more,inters)
  }
  return(seg_df)
}