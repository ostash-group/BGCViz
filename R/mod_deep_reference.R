#' deep_reference UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_deep_reference_ui <- function(id){
  ns <- NS(id)
  shiny::tagList(
            div(
              id="id2",
              shinyjqui::jqui_resizable(
                shinydashboardPlus::box(
                  title = "Annotation comparison to the reference",
                  id = "annotation_reference_comparison_box",
                  collapsible = TRUE,                                          
                  closable = TRUE,
                  width = NULL,
                  height = "100%",
                  shiny::selectInput("ref", "Choose reference data", choices = c(""),
                                     selected = ""),
                  plotly::plotlyOutput(ns("deep_reference"))  %>%
                    shinycssloaders::withSpinner()
                ), options = list(handles="w,e"))
  )
  )
}
    
#' deep_reference Server Functions
#'
#' @noRd 
mod_deep_reference_server <- function(id, vals){
 shiny:: moduleServer( id, function(input, output, session){
    ns <- session$ns

    output$deep_reference <- plotly::renderPlotly({
      shiny::req(vals$deep_reference_to_plot)
      vals$can_plot_deep_ref_2 <- T
      vals$deep_reference_to_plot
    })
    
 })
}
    
## To be copied in the UI
# mod_deep_reference_ui("deep_reference_ui_1")
    
## To be copied in the server
# mod_deep_reference_server("deep_reference_ui_1", vals=vals)
