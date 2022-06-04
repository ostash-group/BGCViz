#' refine_unique
#'
#' @description A function, that deletes identical values in a string, were numbers are separated with ','
#'
#' @return List, were all same numbers are deleted, only one instance is left
#'
#' @noRd
refine_unique <- function(data) {
    n <- utils::tail(data, n = 1)
    data <- utils::head(data, -1)
    n_list <- stringr::str_split(n, ",")
    out <- sapply(n_list[[1]], function(x) {
        x %in% unlist(stringr::str_split(data, ","))
    })
    res <- sapply(out, function(x) {
        if (x == F) {
            x
        }
    })

    return(paste(names(Filter(Negate(is.null), res)), collapse = ","))
}
