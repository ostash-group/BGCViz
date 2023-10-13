read_emerald <- function(data) {
  # get rid off unneeded rows
  all <- readLines(data)
  filtered_lines <- all[!grepl("^#|^$", all)]
  data <- paste(filtered_lines, collapse = "\n")
  data_connection <- textConnection(data)
  # create  dataframe
  emerald_data <- read.table(data_connection, header = FALSE, sep = "\t", col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"))
  attributes_split <- strsplit(emerald_data$attributes, ";")
  attribute_list <- vector("list", length(attributes_split))
  key_value_pairs <- strsplit(attributes_split[[1]], "=")
  key <- lapply(key_value_pairs, function(x) x[1])
  for (i in 1:length(attributes_split)) {
    value <- lapply(key_value_pairs, function(x) x[2])
    attribute_list[[i]] <- setNames(as.list(value), unlist(key))
  }

  attributes_df <- do.call(rbind, attribute_list)
  emerald_data <- cbind(emerald_data, attributes_df)
  return (emerald_data)
}

data <- read_emerald("/Users/levchik/Downloads/streptomyces_coelicolor.gff") 
data
