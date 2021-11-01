# Simple segment dataframe generation

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

# Geom generating functions
geom_anti <- function(data, rre_more){
  ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, 
                                                ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
}
geom_prism <- function(data,rre_more){
  ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,                                                                                                             ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
}
geom_deep <- function(data,rre_more){
  ggplot2::geom_segment(data=data,ggplot2::aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                               ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
                                               deepbgc_score = deepbgc_score,activity = activity ),size =3)
}
geom_rre <- function(data, rre_more){
  if (rre_more == T){
    ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type, Score = Score, Software = Software,
                                                  ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value,
                                                  P_value = P_value, RRE_start = RRE_start,RRE_stop = RRE_stop, 
                                                  Probability = Probability),size = 3)
  } else {
    ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                  ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value
    ),size = 3)
  }
}
geom_sempi <- function(data,rre_more){
  
  ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,
                                                ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
}
geom_prism_supp <- function(data,rre_more){
  ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, ID = ID,
                                                Start = Start, Stop = Stop, Type = Type, Name = Name, Full_name = Full_name,
                                                Score = Score), size = 3)
}
geom_arts <- function(data,rre_more){
  ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,                                                                                                             ID = ID, Start = Start, Stop = Stop, Type = Type, Hit = Hit, 
                                                Core = Core, E_value = E_value, Bitscore = Bitscore, Count = Count, Model = Model), size = 3)
}
geom_gecco <- function(data,rre_more){
  ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, 
                                                ID = ID, Start = Start, Stop = Stop, Type = Type, Num_proteins= Num_proteins,
                                                Num_domains = Num_domains,Average_p = Average_p, Max_p = Max_p ), size = 3)
}

# Functions to add more information to segment datadrame

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

add_more_annot <- function(seg_df, plot, soft_names, index, rre_more){
  if (dim(seg_df)[1] > 0){
    if (soft_names[index] == "anti"){
      plot <- plot + geom_anti(seg_df)
    } else if (soft_names[index] == "sempi") {
      plot <- plot + geom_sempi(seg_df)
    }
    else if (soft_names[index] == "prism") {
      plot <- plot + geom_prism(seg_df)
    }
    else if (soft_names[index] == "prism_supp") {
      plot <- plot + geom_prism_supp(seg_df)
    }
    else if (soft_names[index] == "arts") {
      plot <- plot + geom_arts(seg_df)
    }
    else if (soft_names[index] == "deep") {
      plot <- plot + geom_deep(seg_df)
    }
    else if (soft_names[index] == "rre") {
      plot <- plot + geom_rre(seg_df,rre_more)
    } else if (soft_names[index] == "gecco") {
      plot <- plot+geom_gecco(seg_df)
    }
    return(plot)
  } else{
    return(plot)
  }
}


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