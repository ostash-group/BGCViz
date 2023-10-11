#' ARTS tree UI functions
#'

mod_arts_tree_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shiny::fluidRow(
      tags$div(
        id = "arts_tree_data1",
        div(
          id = "id1",
          shinydashboardPlus::box(
            title = "Phylogenetic tree",
            id = "arts_tree_box",
            collapsible = TRUE,
            closable = TRUE,
            width = 12,
            shiny::selectInput(ns("phylo_file"), "Choose a file to build a tree", choices = c(), selected = ""),
            div(
              style = "height: 600px; overflow-y: scroll; overflow-x: scroll",  # Adjust height as needed, may be needed indeed
              shiny::plotOutput(ns("arts_tree"), height = "2000px",width = "1500px") %>%
                shinycssloaders::withSpinner()
            )
          )
        )
      )
    )
  )
}

#'arts_tree server function
#'

mod_arts_tree_server <- function(id, vals) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # observing changes of input for tree plot
    observe({
      shiny::updateSelectInput(
        session,
        'phylo_file',             
        choices = paste0(vals$arts_tree_data$TreesFiles),
        selected = vals$arts_tree_data$TreesFiles[1]       
      )
    })
    
    # Define a reactive expression for the tree
    tree_data <- reactive({
      # Create the tree object
      tree$core <- vals$arts_tree_data$Trees[vals$arts_tree_data$TreesFiles == input$phylo_file][[1]]
      tree$type <- "rectangular"
      return(tree)
    })
    
    # Render the plot directly within renderPlot
    output$arts_tree <- renderPlot(res = 90,{
      req(vals$arts_data_input == TRUE)  
      
      # Create and render the plot
      tree_plot <- ggtree::ggtree(tree_data()$core, layout = tree_data()$type) +
        ggtree::geom_tree() +
        ggtree::theme_tree() +
        ggtree::geom_tiplab(size = 2.2, color = 'firebrick')
      
      return(tree_plot)
    })
  })
}
