#' deep_reference_2 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_deep_reference_2_ui <- function(id) {
    ns <- NS(id)
    tagList(
        div(
            id = "anno_div_1",
            shinyjqui::jqui_resizable(
                shinydashboardPlus::box(
                    title = "Annotations reference",
                    id = "annotation_reference_box",
                    height = "100%",
                    width = NULL,
                    collapsible = TRUE,
                    closable = TRUE,
                    plotly::plotlyOutput(ns("deep_reference_2")) %>%
                        shinycssloaders::withSpinner()
                ),
                options = list(handles = "w,e")
            )
        )
    )
}

#' deep_reference_2 Server Functions
#'
#' @noRd
mod_deep_reference_2_server <- function(id, vals, data_uploads, data_to_use) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        # Silence R CMD note
        x <- y <- yend <- xend <-
            ID <- Software <- Chr <-
            Type2 <- Start <- Stop <-
            Type <- num_domains <- deepbgc_score <-
            activity <- Score <- E_value <- P_value <-
            RRE_start <- RRE_stop <- Probability <-
            Name <- Full_name <- Hit <- Core <-
            Bitscore <- Count <- Model <- Num_proteins <-
            Num_domains <- Average_p <- Max_p <- NULL
        output$deep_reference_2 <- plotly::renderPlotly({
            shiny::req(vals$can_plot_deep_ref_2 == TRUE)
            vals$can_plot_deep_ref_2 == FALSE
            rename_y_axis <- shiny::isolate(vals$rename_y_axis)
            data <- NULL

            index <- 1
            for (upload in data_uploads) {
                if (is.null(data)) {
                    if (vals[[upload]] == TRUE) {
                        if (dim(vals[[data_to_use[index]]])[1] != 0) {
                            data <- vals[[data_to_use[index]]]
                        }
                    }
                }
                index <- index + 1
            }


            tooltip <- c(
                "Software", "ID", "Start", "Stop", "Type", "num_domains", "deepbgc_score", "activity", "Score", "E_value",
                "P_value", "RRE_start", "RRE_stop", "Probability", "Name", "Full_name", "Hit", "Core", "Count", "Bitscore", "Model",
                "Num_domains", "Num_proteins", "Average_p", "Max_p"
            )

            plot <- ggplot2::ggplot(data, ggplot2::aes(x = vals$chr_len, y = Chr))

            if (vals$emerald_data_input == TRUE){
              plot <- plot +
                suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_emer, ggplot2::aes(x, y,
                                                                                               xend = xend, yend = yend, color = Type2, Software = Software,
                                                                                               ID = ID, Start = Start, Stop = Stop, Type = Type
                ), size = 3))
            }
            if (vals$compare_data_input == TRUE){
              plot <- plot +
                suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_compare, ggplot2::aes(x, y,
                                                                                                    xend = xend, yend = yend, color = Type2, Software = Software,
                                                                                                    ID = ID, Start = Start, Stop = Stop, Type = Type
                ), size = 3))
            }
            if (vals$ripp_data_input == TRUE){
              plot <- plot +
                      suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_ri, ggplot2::aes(x, y,
                          xend = xend, yend = yend, color = Type2, Software = Software,
                          ID = ID, Start = Start, Stop = Stop, Type = Type
                      ), size = 3))
            }
            if (vals$anti_data_input == TRUE) {
                plot <- plot +
                    suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_a, ggplot2::aes(x, y,
                        xend = xend, yend = yend, color = Type2, Software = Software,
                        ID = ID, Start = Start, Stop = Stop, Type = Type
                    ), size = 3))
            }
            if (vals$deep_data_input == TRUE) {
                if (dim(vals$seg_df_ref_d)[1] > 0) {
                    plot <- plot +
                        suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_d, ggplot2::aes(x, y,
                            xend = xend, yend = yend, color = Type2, Software = Software,
                            ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
                            deepbgc_score = deepbgc_score, activity = activity
                        ), size = 3))
                }
            }
            if (vals$rre_data_input == TRUE) {
                if (vals$rre_more == TRUE) {
                    plot <- plot + suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_r, ggplot2::aes(x, y,
                        xend = xend, yend = yend, color = Type2, Score = Score, Software = Software,
                        ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value,
                        P_value = P_value, RRE_start = RRE_start, RRE_stop = RRE_stop,
                        Probability = Probability
                    ), size = 3))
                } else {
                    plot <- plot + suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_r, ggplot2::aes(x, y,
                        xend = xend, yend = yend, color = Type2, Software = Software,
                        ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value
                    ), size = 3))
                }
            }
            if (vals$prism_data_input == TRUE) {
                plot <- plot + suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_p, ggplot2::aes(x, y,
                    xend = xend, yend = yend, color = Type2, Software = Software,
                    ID = ID, Start = Start, Stop = Stop, Type = Type
                ), size = 3))
            }
            if (vals$sempi_data_input == TRUE) {
                plot <- plot + suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_s, ggplot2::aes(x, y,
                    xend = xend, yend = yend, color = Type2, Software = Software,
                    ID = ID, Start = Start, Stop = Stop, Type = Type
                ), size = 3))
            }
            if (vals$prism_supp_plot == TRUE) {
                plot <- plot + suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_p_s, ggplot2::aes(x, y,
                    xend = xend, yend = yend, color = Type2, Software = Software, ID = ID,
                    Start = Start, Stop = Stop, Type = Type, Name = Name, Full_name = Full_name,
                    Score = Score
                ), size = 3))
            }
            if (vals$arts_data_input == TRUE) {
                plot <- plot + suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_ar, ggplot2::aes(x, y,
                    xend = xend, yend = yend, color = Type2, Software = Software,
                    ID = ID, Start = Start, Stop = Stop, Type = Type, Hit = Hit,
                    Core = Core, E_value = E_value, Bitscore = Bitscore, Count = Count,
                    Model = Model
                ), size = 3))
            }
            if (vals$gecco_data_input == TRUE) {
                if (dim(vals$seg_df_ref_g)[1] > 0) {
                    plot <- plot + suppressWarnings(ggplot2::geom_segment(data = vals$seg_df_ref_g, ggplot2::aes(x, y,
                        xend = xend, yend = yend, color = Type2, Software = Software,
                        ID = ID, Start = Start, Stop = Stop, Type = Type, Num_proteins = Num_proteins,
                        Num_domains = Num_domains, Average_p = Average_p, Max_p = Max_p
                    ), size = 3))
                }
            }
            to_plot <- plotly::ggplotly(plot +
                ggplot2::scale_y_discrete(labels = rename_y_axis) +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10)) +
                ggplot2::ylab("") +
                ggplot2::xlab("Chromosome length") +
                ggplot2::theme(legend.title = ggplot2::element_blank()) +
                ggplot2::ggtitle("All annotations"),
            # What actually to visualize in tooltip
            tooltip = tooltip
            )
            to_plot %>% plotly::layout(
                legend = list(
                    font = list(
                        family = "sans-serif",
                        size = 12,
                        color = "#000"
                    ),
                    bordercolor = "#FFFFFF",
                    borderwidth = 2,
                    title = list(text = "<b> Cluster Types </b>")
                ),
                autosize = TRUE
            )
        }) # %>% shiny::debounce(200)
    })
}

## To be copied in the UI
# mod_deep_reference_2_ui("deep_reference_2_ui_1")

## To be copied in the server
# mod_deep_reference_2_server("deep_reference_2_ui_1", vals=vals, data_uploads = data_uploads, data_to_use = data_to_use)
