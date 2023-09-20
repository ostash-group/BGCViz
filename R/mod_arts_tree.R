#' ARTS tree UI functions
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList


mod_arts_tree_ui <- function(id)
{
  ns <- NS(id)
  tagList(
    shiny::fluidRow(
      tags$div(
        id = "arts_tree_data1",
        div(
          id = "id1",
          shinydashboardPlus::box(
            title = "ARTS tree",
            id = "arts_tree_box",
            collapsible = TRUE,
            closable = TRUE,
            width = 12
          )
          
        )
      )
      
    )
      
    )
}