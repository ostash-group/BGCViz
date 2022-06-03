#' initialize_biocircos
#'
#' @description Function, that initialized biocircos plot, setting proper chromosomes, length, arcs and their parameters
#'
#' @return Biocircos chromosomes, arcs, their start, end, labels, color in hexadecimal format
#'
#' @noRd
initialize_biocircos <- function(biocircos_anti, name, Biocircos_chromosomes, arcs_chromosomes, arcs_begin, arcs_end, arc_labels, arc_col, rename_data, chr_len, biocircos_color, coloring_datatable) {
  # Make chromosome list for Biocircos plot. Use chr_len as an input
  Biocircos_chromosomes[[name]] <- chr_len
  # Add arcs. Quantity of arcs is length of dataframes
  arcs_chromosomes <- c(arcs_chromosomes, rep(name, length(biocircos_anti$Cluster)))
  # Add arcs begin positions. (Start column)
  arcs_begin <- c(arcs_begin, biocircos_anti$Start)
  # Stop position of arcs.
  arcs_end <- c(arcs_end, biocircos_anti$Stop)
  # Add Arcs labels. Can add only one label...
  arc_labels <- c(arc_labels, biocircos_anti$Type)
  if ((biocircos_color == T)) {
    arc_colors <- sapply(biocircos_anti$Type2, function(x) {
      if (x %in% coloring_datatable$x$data$Name) {
        coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == x]
      } else {
        coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"]
      }
    })
  } else {
    arc_colors <- coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"]
  }
  arc_col <- c(arc_col, as.character(arc_colors))
  return(list(Biocircos_chromosomes, arcs_chromosomes, arcs_begin, arcs_end, arc_labels, arc_col))
}

#' add_biocircos_data
#'
#' @description Function, that adds data to the biocircos plot. It get the interception, make links between chromosomes and labels and color them (according to a chosen scheme)
#'
#' @return Chromosomes start and end, link start, end, and their color in hexadecimal format
#'
#' @noRd

add_biocircos_data <- function(data1_inter, data2_inter, data1, data2, data1_label, data2_label, rename_data, class, ref_col_biocircos, coloring_datatable) {
  inter_s_rre_n <- data1_inter
  inter_rre_s <- data2_inter
  # Add link start. Just populate certain chromosome name times the lenght of interception
  chromosomes_start <- c(rep(data2_label, length(inter_rre_s)))
  # Add link end. Just populate second output from the vectors, used above.
  chromosomes_end <- c(rep(data1_label, length(inter_s_rre_n)))
  # Add links start positions as a start from dataframe. This vector is for chromosome start
  link_pos_start <- as.numeric(c(data2$Start[match(inter_rre_s, data2$Cluster)]))
  # Add links start positions as a start from dataframe. For chromosome start variable
  link_pos_start_1 <- as.numeric(c(data2$Stop[match(inter_rre_s, data2$Cluster)]))
  # Add links start position for a chromosome stop variable
  link_pos_end <- as.numeric(c(data1$Start[match(inter_s_rre_n, data1$Cluster)]))
  # Add links start position for a chromosome stop position
  link_pos_end_2 <- as.numeric(c(data1$Stop[match(inter_s_rre_n, data1$Cluster)]))
  label_1 <- c(sapply(inter_rre_s, function(x) {
    x <- paste(paste0(data2_label, ":"), x, ",", data2$Type[data2$Cluster == x])
  }))
  label_2 <- c(sapply(inter_s_rre_n, function(x) {
    x <- paste(paste0(data1_label, ":"), x, ",", data1$Type[data1$Cluster == x])
  }))
  # browser()
  if (!is.null(inter_rre_s)) {
    if (class == "P") {
      subset_vec <- data2$Type2[match(inter_rre_s, data2$Cluster)] == data1$Type2[match(inter_s_rre_n, data1$Cluster)]
      label_color <- as.character(c(sapply(data2$Type2[match(inter_rre_s, data2$Cluster)], function(x) {
        if (x %in% coloring_datatable$x$data$Name) {
          as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == x])
        } else {
          as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"])
        }
      })))
      if (length(label_color) != 0) {
        for (t in seq(1:length(label_color))) {
          if (!is.null(subset_vec[t])) {
            if (subset_vec[t] == F) {
              label_color[t] <- as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"])
            }
          }
        }
      }
    } else if (class == "H") {
      if (grep(paste0("^", data1_label, "$"), coloring_datatable$x$data$Hierarchy) < (grep(paste0("^", data2_label, "$"), coloring_datatable$x$data$Hierarchy))) {
        label_color <- as.character(c(sapply(data1$Type2[match(inter_s_rre_n, data1$Cluster)], function(x) {
          if (x %in% coloring_datatable$x$data$Name) {
            as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == x])
          } else {
            as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"])
          }
        })))
      } else {
        label_color <- as.character(c(sapply(data2$Type2[match(inter_rre_s, data2$Cluster)], function(x) {
          if (x %in% coloring_datatable$x$data$Name) {
            as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == x])
          } else {
            as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"])
          }
        })))
      }
    } else if (class == "R") {
      if (data2_label == ref_col_biocircos) {
        label_color <- as.character(c(sapply(data1$Type2[match(inter_s_rre_n, data1$Cluster)], function(x) {
          if (x %in% coloring_datatable$x$data$Name) {
            as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == x])
          } else {
            as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"])
          }
        })))
      } else if (data1_label == ref_col_biocircos) {
        label_color <- as.character(c(sapply(data2$Type2[match(inter_rre_s, data2$Cluster)], function(x) {
          if (x %in% coloring_datatable$x$data$Name) {
            as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == x])
          } else {
            as.character(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"])
          }
        })))
      } else {
        label_color <- as.character(rep(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"], length(chromosomes_start)))
      }
    } else {
      label_color <- as.character(rep(coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"], length(chromosomes_start)))
    }
  }
  return(list(
    inter_s_rre_n, inter_s_rre_n, chromosomes_start, chromosomes_end, link_pos_start, link_pos_start_1, link_pos_end,
    link_pos_end_2, label_1, label_2, label_color
  ))
}
