#' biocircos UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_biocircos_ui <- function(id) {
    ns <- NS(id)
    tagList(
        shiny::fluidRow(
            tags$div(
                id = "biocircos_data1",
                div(
                    id = "id1",
                    shinydashboardPlus::box(
                        title = "Biocircos plot",
                        id = "biocircos_plot_box",
                        collapsible = TRUE,
                        closable = TRUE,
                        width = 12,
                        shiny::checkboxInput(ns("ShowBiocircosColoring"), "Show Biocircos coloring scheme"),
                        sidebar = shinydashboardPlus::boxSidebar(
                            id = "biocircos_box_sidebar",
                            width = 25,
                            shiny::checkboxInput("biocircos_color", "Make arcs in biocircos colorful, based on the class"),
                            shiny::checkboxInput("label_color", "Make links in biocircos colorful, based on the class"),
                            shiny::selectInput("label_color_class", "Choose the mode to color the links",
                                choices = c(
                                    "Hierarchical-based" = "H",
                                    "Purity-based" = "P",
                                    "Reference column-based" = "R"
                                ),
                                selected = "H"
                            ),
                            shiny::selectInput("ref_col_biocircos", "Choose reference column to color the links", choices = c(""), selected = "")
                        ),
                        BioCircos::BioCircosOutput(ns("biocircos"), height = "900px") %>%
                            shinycssloaders::withSpinner()
                    )
                )
            )
        ),
        shiny::fluidRow(
            tags$div(
                id = "biocircos_data2",
                div(
                    id = "id1",
                    shiny::uiOutput(ns("biocircos_coloring"))
                )
            )
        )
    )
}

#' biocircos Server Functions
#'
#' @noRd
mod_biocircos_server <- function(id, vals) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        output$biocircos <- BioCircos::renderBioCircos({
            shiny::req(vals$data_upload_count > 1)

            # Plot BioCircos
            BioCircos::BioCircos(vals$tracklist, genome = vals$Biocircos_chromosomes, genomeTicksScale = 1e+6)
        })


        output$biocircos_legend <- DT::renderDataTable({
            shiny::req(vals$data_upload_count >= 1)
            rownames <- FALSE
            new_data <- vals$coloring_datatable
            color_vec <- new_data$x$data$Color
            options(DT.options = list(pageLength = 50))
            new_data %>% DT::formatStyle("Color", backgroundColor = DT::styleEqual(color_vec, color_vec))
        })

        output$biocircos_coloring <- shiny::renderUI({
            if (input$ShowBiocircosColoring == TRUE) {
                shinydashboardPlus::box(
                    title = "Biocircos coloring scheme",
                    closable = TRUE,
                    collapsible = TRUE,
                    DT::dataTableOutput(ns("biocircos_legend")) %>%
                        shinycssloaders::withSpinner()
                )
            }
        })

        # Updating values in Datatable on edit
        shiny::observeEvent(input$biocircos_legend_cell_edit, {
            if (input$biocircos_legend_cell_edit$col[1] == 0) {
                vals$coloring_datatable$x$data$Name <- input$biocircos_legend_cell_edit$value
            } else if (input$biocircos_legend_cell_edit$col[1] == 1) {
                vals$coloring_datatable$x$data$Color <- input$biocircos_legend_cell_edit$value
            } else if (input$biocircos_legend_cell_edit$col[1] == 2) {
                vals$coloring_datatable$x$data$Hierarchy <- input$biocircos_legend_cell_edit$value
            }
        })
    })
}

## To be copied in the UI
# mod_biocircos_ui("biocircos_ui_1")

## To be copied in the server
# mod_biocircos_server("biocircos_ui_1", vals = vals)
