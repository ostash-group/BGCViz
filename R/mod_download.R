#' download UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_download_ui <- function(id) {
    ns <- NS(id)
    tagList(
        div(
            id = "id6",
            shinydashboardPlus::box(
                title = "Download data",
                id = "download_data_box",
                collapsible = TRUE,
                closable = TRUE,
                width = NULL,
                shiny::downloadButton(ns("download"), "Download currently used datasets (as for Biocircos plot)")
            )
        )
    )
}

#' download Server Functions
#'
#' @noRd
mod_download_server <- function(id) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        output$download <- shiny::downloadHandler(
            filename = function() {
                paste("datasets.zip")
            },
            content = function(file) {
                flst <- c()
                # List files in directory
                files_in_dir <- list.files()
                # Iterate over those files and if found "_biocircos.csv" add to the flst vector
                for (file_names in files_in_dir) {
                    if (grepl("_biocircos.csv", file_names, fixed = TRUE)) {
                        flst <- c(flst, file_names)
                    } else if (grepl("group_by.csv", file_names, fixed = TRUE)) {
                        flst <- c(flst, file_names)
                    }
                }
                # create the zip file from flst vector
                group_by_script <- system.file("scripts", "group.py", package = "BGCViz")
                dissect_script <- system.file("scripts", "dissect.py", package = "BGCViz")
                flst <- c(flst, group_by_script, dissect_script)
                utils::zip(file, flst, flags = '-r9Xj')
            },
            contentType = "application/zip"
        )
    })
}

## To be copied in the UI
# mod_download_ui("download_ui_1")

## To be copied in the server
# mod_download_server("download_ui_1")
