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
            shiny::selectInput("phylo_file", "Choose a file to build a tree", choices = c(), selected = ""),
            div(
              style = "height: 600px; overflow-y: scroll",  # Adjust height as needed, may be needed indeed
              shiny::plotOutput(ns("arts_tree"), height = "2000px") %>%
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
    
    # Define a reactive expression for the tree
    tree_data <- reactive({
      tree <- vals$arts_data$Trees[[1]]
      
    ### IF YOU WANT TO SHORTER STRAINS ID UNCOMENT CODE BELOW  ###
    
      # tree$tip.label <- sapply(tree$tip.label, function(x) {
      #   split_string <- unlist(strsplit(x, "_"))
      #   return(paste(split_string[2], split_string[3], sep = " "))
      # })
      
      return(tree)
    })
    
    # Render the plot directly within renderPlot
    output$arts_tree <- renderPlot(res = 90,{
      req(vals$arts_data_input == TRUE)  
      
      # Create and render the plot
      tree_plot <- ggtree::ggtree(tree_data()) +
        ggtree::geom_tree() +
        ggtree::theme_tree() +
        ggtree::geom_tiplab(size = 2.2, color = 'firebrick')
      
      return(tree_plot)
    })
  })
}
