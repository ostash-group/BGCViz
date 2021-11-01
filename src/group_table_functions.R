# Refining unique values for table to out in one row
refine_unique <- function(data){
  n <- tail(data, n=1)
  data <- head(data, -1)
  n_list <-  stringr::str_split(n, ",")
  out <- sapply(n_list[[1]], function(x){x %in% unlist(stringr::str_split(data, ","))})
  res <- sapply(out, function(x){
    if (x==F){
      x
    }
  })
  
  return(paste(names(Filter(Negate(is.null), res)), collapse = ","))
}
