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
            title = "Select a tree of interest",
            id = "arts_tree_box",
            collapsible = TRUE,
            closable = TRUE,
            width = 12,
            shiny::selectInput(ns("group_by"), "Group data by", choices = c(), selected = ""),
            shiny::plotOutput(ns("arts_tree"),height = "1000px") %>%
              shinycssloaders::withSpinner()
          )
          
        )
      )
      
    )
      
    )
}

#'arts_tree server function
#'
#' @noRd
mod_arts_tree_server <- function(id,path){
  
}
