#' barplot_rank UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_barplot_rank_ui <- function(id) {
    ns <- NS(id)
    tagList(
        div(
            id = "id1",
            shinyjqui::jqui_resizable(
                shinydashboardPlus::box(
                    title = "Ranking barplot",
                    id = "ranking_barplot_box",
                    collapsible = TRUE,
                    closable = TRUE,
                    height = "100%",
                    plotly::plotlyOutput(ns("barplot_rank"), height = "600px") %>%
                        shinycssloaders::withSpinner()
                ),
                options = list(handles = "w,e")
            )
        )
    )
}

#' barplot_rank Server Functions
#'
#' @noRd
mod_barplot_rank_server <- function(id, vals, data_uploads, soft_names, soft_namings, data_to_use, abbr) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        # Silence R CMD note
        Cluster <- Count <- Type <-
            Start <- Start <- Stop <-
            Label <- NULL
        output$barplot_rank <- plotly::renderPlotly({
            shiny::req(vals$data_upload_count > 1)
            shiny::req(vals$need_filter == FALSE)
            shiny::req(vals$can_plot_barplot_rank == TRUE)

            antismash_count <- NULL
            prism_count <- NULL
            deep_count <- NULL
            rre_count <- NULL
            sempi_count <- NULL
            prism_supp_count <- NULL
            arts_count <- NULL
            gecco_count <- NULL

            if (is.null(vals$inters_filtered)) {
                inters <- vals$inters
            } else {
                inters <- vals$inters_filtered
            }
            index <- 1
            ranking_data <- NULL
            for (upload in data_uploads) {
                if (vals[[upload]] == TRUE) {
                    counts_var <- plyr::count(as.factor(unlist(sapply(inters[[soft_names[index]]], function(x) {
                        x$to
                    }))))
                    # Check if ID is in dataframe and if it is - extract all information about to the local dataframe
                    anot_var <- vals[[data_to_use[index]]][vals[[data_to_use[index]]]$Cluster %in% as.numeric(levels(counts_var$x)), ]
                    # Add prefices to the ID to plot for a barplot.
                    counts_var$x <- sapply(counts_var$x, function(x) paste0(abbr[index], ": ", x))
                    # Add label column to the dataframe, from which we will plot
                    counts_var$label <- rep(soft_namings[index], length(counts_var$x))
                    # Add type to the dataframe, from which we would plot (from annotation dataframe)
                    counts_var$Type <- anot_var$Type
                    # Add Start positions (to visualize on hover)
                    counts_var$Start <- anot_var$Start
                    # Add Stop positions (to visualize on hover)
                    counts_var$Stop <- anot_var$Stop
                    if (is.null(ranking_data)) {
                        ranking_data <- counts_var
                    } else {
                        ranking_data <- rbind(ranking_data, counts_var)
                    }
                }
                index <- index + 1
            }


            # Fix column names in the master dataframe
            colnames(ranking_data) <- c("Cluster", "Count", "Label", "Type", "Start", "Stop")
            # Plot
            plotly::ggplotly(ggplot2::ggplot(ranking_data, ggplot2::aes(x = Cluster, y = Count, Type = Type, Start = Start, Stop = Stop)) +
                ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = Label)) +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 10),
                    axis.text.y = ggplot2::element_text(size = 14)
                ) +
                ggplot2::ggtitle("Number of times cluster is annotated with other tool"),
            tooltip = c("Type", "Start", "Stop")
            )
        })
    })
}

## To be copied in the UI
# mod_barplot_rank_ui("barplot_rank_ui_1")

## To be copied in the server
# mod_barplot_rank_server("barplot_rank_ui_1", vals=vals, data_uploads = data_uploads, soft_names = soft_names, soft_namings = soft_namings, data_to_use = data_to_use, abbr = abbr)
