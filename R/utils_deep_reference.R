#' geom_anti
#'
#' @description A function, that returns antismash geom with the legend,
#' specific to this software (to show on mouse hover).
#'
#' @return geom_segment with specific fields
#'
#' @noRd
geom_anti <- function(data, rre_more) {
    # Silence R CMD note
    x <- y <- xend <- yend <- Type2 <-
        Software <- ID <- Start <- Stop <- Type <- NULL
    ggplot2::geom_segment(data = data, ggplot2::aes(x, y,
        xend = xend, yend = yend, color = Type2, Software = Software,
        ID = ID, Start = Start, Stop = Stop, Type = Type
    ), size = 3)
}
#' geom_prism
#'
#' @description A function, that returns prism geom with the legend, specific to this software (to show on mouse hover).
#'
#' @return geom_segment with specific fields
#'
#' @noRd
geom_prism <- function(data, rre_more) {
    # Silence R CMD note
    x <- y <- xend <- yend <- Type2 <-
        Software <- ID <- Start <- Stop <- Type <- NULL
    ggplot2::geom_segment(data = data, ggplot2::aes(x, y,
        xend = xend, yend = yend, color = Type2, Software = Software,
        ID = ID, Start = Start, Stop = Stop, Type = Type
    ), size = 3)
}
#' geom_deep
#'
#' @description A function, that returns deepbgc geom with the legend,
#' specific to this software (to show on mouse hover).
#'
#' @return geom_segment with specific fields
#'
#' @noRd
geom_deep <- function(data, rre_more) {
    # Silence R CMD note
    x <- y <- xend <- yend <- Type <-
        Software <- ID <- Start <- Stop <-
        Type <- num_domains <- deepbgc_score <-
        activity <- NULL
    ggplot2::geom_segment(data = data, ggplot2::aes(x, y,
        xend = xend, yend = yend, color = Type, Software = Software,
        ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
        deepbgc_score = deepbgc_score, activity = activity
    ), size = 3)
}
#' geom_rre
#'
#' @description A function, that returns RRE-Finder geom with the legend,
#' specific to this software (to show on mouse hover).
#'
#' @return geom_segment with specific fields
#'
#' @noRd
geom_rre <- function(data, rre_more) {
    # Silence R CMD note
    x <- y <- xend <- yend <- Type <-
        Score <- Software <- ID <- Start <-
        Stop <- Type <- E_value <- P_value <- RRE_start <-
        RRE_stop <- Probability <- NULL
    if (rre_more == T) {
        ggplot2::geom_segment(data = data, ggplot2::aes(x, y,
            xend = xend, yend = yend, color = Type, Score = Score, Software = Software,
            ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value,
            P_value = P_value, RRE_start = RRE_start, RRE_stop = RRE_stop,
            Probability = Probability
        ), size = 3)
    } else {
        ggplot2::geom_segment(data = data, ggplot2::aes(x, y,
            xend = xend, yend = yend, color = Type, Software = Software,
            ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value
        ), size = 3)
    }
}
#' geom_sempi
#'
#' @description A function, that returns SEMPI geom with the legend, specific to this software (to show on mouse hover).
#'
#' @return geom_segment with specific fields
#'
#' @noRd
geom_sempi <- function(data, rre_more) {
    # Silence R CMD note
    x <- y <- xend <- yend <- Type2 <-
        Software <- ID <- Start <- Stop <-
        Type <- NULL
    ggplot2::geom_segment(data = data, ggplot2::aes(x, y,
        xend = xend, yend = yend, color = Type2, Software = Software,
        ID = ID, Start = Start, Stop = Stop, Type = Type
    ), size = 3)
}
#' deep_reference
#'
#' @description A function, that returns Prism-Supplement geom with the legend,
#' specific to this software (to show on mouse hover).
#'
#' @return geom_segment with specific fields
#'
#' @noRd
geom_prism_supp <- function(data, rre_more) {
    # Silence R CMD note
    x <- y <- xend <- yend <- Type2 <-
        Software <- ID <- Start <- Stop <- Type <- Name <-
        Full_name <- Score <- NULL
    ggplot2::geom_segment(data = data, ggplot2::aes(x, y,
        xend = xend, yend = yend, color = Type2, Software = Software, ID = ID,
        Start = Start, Stop = Stop, Type = Type, Name = Name, Full_name = Full_name,
        Score = Score
    ), size = 3)
}
#' geom_arts
#'
#' @description A function, that returns ARTS geom with the legend,
#' specific to this software (to show on mouse hover).
#'
#' @return geom_segment with specific fields
#'
#' @noRd
geom_arts <- function(data, rre_more) {
    # Silence R CMD error
    x <- y <- xend <- yend <- Type2 <-
        Start <- Stop <- Type <- ID <- Hit <-
        Software <- Core <- E_value <-
        Bitscore <- Count <- Model <- NULL
    ggplot2::geom_segment(data = data, ggplot2::aes(x, y,
        xend = xend, yend = yend, color = Type2, Software = Software,
        ID = ID, Start = Start, Stop = Stop, Type = Type, Hit = Hit,
        Core = Core, E_value = E_value, Bitscore = Bitscore, Count = Count, Model = Model
    ), size = 3)
}
#' geom_gecco
#'
#' @description A function, that returns GECCO geom with the legend, specific to this software (to show on mouse hover).
#'
#' @return geom_segment with specific fields
#'
#' @noRd
geom_gecco <- function(data, rre_more) {
    # Silence R CMD note
    x <- y <- xend <- yend <- Type2 <- Software <-
        ID <- Start <- Stop <- Type <- Num_proteins <-
        Num_domains <- Average_p <- Max_p <- NULL
    ggplot2::geom_segment(data = data, ggplot2::aes(x, y,
        xend = xend, yend = yend, color = Type2, Software = Software,
        ID = ID, Start = Start, Stop = Stop, Type = Type, Num_proteins = Num_proteins,
        Num_domains = Num_domains, Average_p = Average_p, Max_p = Max_p
    ), size = 3)
}
#' add_more_annot
#'
#' @description A function, that that adds specific geom to the plot object
#'
#' @return plot oblect with added geom_segment
#'
#' @noRd
add_more_annot <- function(seg_df, plot, soft_names, index, rre_more) {
    if (dim(seg_df)[1] > 0) {
        if (soft_names[index] == "anti") {
            plot <- plot + geom_anti(seg_df)
        } else if (soft_names[index] == "sempi") {
            plot <- plot + geom_sempi(seg_df)
        } else if (soft_names[index] == "prism") {
            plot <- plot + geom_prism(seg_df)
        } else if (soft_names[index] == "prism_supp") {
            plot <- plot + geom_prism_supp(seg_df)
        } else if (soft_names[index] == "arts") {
            plot <- plot + geom_arts(seg_df)
        } else if (soft_names[index] == "deep") {
            plot <- plot + geom_deep(seg_df)
        } else if (soft_names[index] == "rre") {
            plot <- plot + geom_rre(seg_df, rre_more)
        } else if (soft_names[index] == "gecco") {
            plot <- plot + geom_gecco(seg_df)
        }
        return(plot)
    } else {
        return(plot)
    }
}
