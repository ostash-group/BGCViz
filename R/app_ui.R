#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic 
    shinydashboardPlus::dashboardPage(
      shinydashboardPlus::dashboardHeader(title = "BGCViz"),
      shinydashboardPlus::dashboardSidebar(
        width = 350,
        shinydashboard::sidebarMenu(
          id="menu_items",
          style = "white-space: normal;",
          shinydashboard::menuItem("Upload data", tabName = "uploaddata_sidemenu", icon = icon("fas fa-upload")),
          shinydashboard::menuItem("Global options", tabName = "options_sidemenu", icon = icon("fas fa-cogs")),
          shinydashboard::menuItemOutput("deep_sidemenu_out"),
          shinydashboard::menuItemOutput("gecco_sidemenu_out"),
          shinydashboard::menuItemOutput("anno_sidemenu_out"),
          shinydashboard::menuItemOutput("biocircos_sidemenu_out"),
          shinydashboard::menuItemOutput("summarize_sidemenu_out"),
          shinydashboard::menuItem(tabName = "restore_boxes",
                                   actionButton("restore_box", "Restore all boxes", class = "bg-success"))
        )
      ),
      shinydashboard::dashboardBody(
        tags$head( 
          tags$style(HTML(".main-sidebar { font-size: 15px; }")) #change the font size to 20
        ),
        shinyjs::useShinyjs(),
        shinydisconnect::disconnectMessage(
          text = "An error occurred. Please refresh the page and try again. Also, if error persists, then you are welcome to create an issue at https://github.com/ostash-group/BGCViz/issues (:",
          refresh = "Refresh",
          background = "#FFFFFF",
          colour = "#444444",
          refreshColour = "#337AB7",
          overlayColour = "#000000",
          overlayOpacity = 0.6,
          width = 450,
          top = 50,
          size = 22,
          css = ""
        ),
        shinydashboard::tabItems(
          shinydashboard::tabItem( tabName = "deep_sidemenu",
                                   shiny::fluidRow(
                                     tags$div(
                                       id = "deep_data1",
                                       div(id = "id1",
                                           shinyjqui::jqui_resizable(shinydashboardPlus::box(
                                             title = "DeepBGC comparison",
                                             id = "deep_comparison_box",
                                             collapsible = TRUE,                                        
                                             closable = TRUE,
                                             height = "100%",
                                             shiny::plotOutput("deep_barplot", height = "500px",) %>%
                                               shinycssloaders::withSpinner()
                                           ),options = list(handles="w,e"))),
                                       div(id = "id2", 
                                           shinyjqui::jqui_resizable(shinydashboardPlus::box(
                                             title = "DeepBGC rate",
                                             id = "deep_rate_box",
                                             collapsible = TRUE,                                         
                                             height = "100%",
                                             plotly::plotlyOutput("deep_rate", height = "500px",) %>%
                                               shinycssloaders::withSpinner()
                                           ),options = list(handles="w,e"))))),
                                   shiny::fluidRow(
                                     tags$div( id = "deep_data2",
                                               div(id = "id1",
                                                   shinyjqui::jqui_resizable(shinydashboardPlus::box( 
                                                     title = "DeepBGC comparison controls",
                                                     id = "deep_comparison_controls_box",
                                                     collapsible = TRUE,                                          
                                                     closable = TRUE,
                                                     shiny::selectInput("ref_comparison", "Choose data for comparison with DeepBGC", choices = c(""), selected = ''),
                                                     # Score to use for thresholds
                                                     shiny::selectInput("score_type", "Choose score type to set threshold", choices = c("Activity score" = "Activity",
                                                                                                                                        "Cluster_type score" = "Cluster_Type",
                                                                                                                                        "DeepBGC score" = "DeepBGC"),
                                                                        selected = "Activity score"),
                                                     # Chose step for barplot (as a threshold to draw a bar)
                                                     shiny::sliderInput("plot_step", "Choose step for plots(barplot)", min = 1, max = 50,value = 10),
                                                     shiny::sliderInput("plot_start", "Chose plot start point(barplot)", min = 0, max = 99, value = 0)
                                                   ),options = list(handles="w,e")))
                                     )),
                                   sortable::sortable_js("deep_data1", options = sortable::sortable_options(swap = TRUE, group = "deep_data")),
                                   sortable::sortable_js("deep_data2", options = sortable::sortable_options(swap = TRUE, group = "deep_data"))
          ),
          shinydashboard::tabItem(
            tabName = "gecco_sidemenu",
            shiny::fluidRow(
              tags$div(
                id = "gecco_data1",
                div(
                  id = "id1",
                  shinyjqui::jqui_resizable(shinydashboardPlus::box(
                    title = "GECCO comparison",
                    id = "gecco_comparison_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    height = "100%",
                    shiny::plotOutput("gecco_barplot", height = "500px") %>%
                      shinycssloaders::withSpinner()
                  ),options = list(handles="w,e"))
                ),
                div(
                  id = "id2",
                  shinyjqui::jqui_resizable(shinydashboardPlus::box(
                    title = "GECCO rate",
                    id = "gecco_rate_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    height = "100%",
                    plotly::plotlyOutput("gecco_rate", height = "500px",)%>%
                      shinycssloaders::withSpinner()
                  ),options = list(handles="w,e"))
                ),
              )
            ),
            shiny::fluidRow(
              tags$div(
                id = "gecco_data2",
                div(
                  id = "id1",
                  shinyjqui::jqui_resizable(shinydashboardPlus::box(
                    title = "GECCO comparison controls",
                    id = "gecco_comparison_controls_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    shiny::selectInput("ref_comparison_gecco", "Choose data for comparison with Gecco", choices = c(""),selected = ''),
                    shiny::selectInput("score_type_gecco", "Choose score type to set threshold", choices = c(
                      "Average p-value" = "avg_p",
                      "Cluster_type score" = "Cluster_Type"),
                      selected = "avg_p"),
                    shiny::sliderInput("plot_step_gecco", "Choose step for plots(barplot)", min = 1, max = 50,value = 10),
                    shiny::sliderInput("plot_start_gecco", "Chose plot start point(barplot)", min = 0, max = 99, value = 0)
                  ), options = list(handles="w,e"))
                )
              )
            ),
            sortable::sortable_js("gecco_data1", options = sortable::sortable_options(swap = TRUE, group = "gecco_data")),
            sortable::sortable_js("gecco_data2", options = sortable::sortable_options(swap = TRUE, group = "gecco_data"))
          ),
          shinydashboard::tabItem(
            tabName = "anno_sidemenu",
            shiny::fluidRow(
              tags$div(
                id="anno_data1",
                shiny::column(
                  width = 12,
                  mod_deep_reference_2_ui("deep_reference_2_ui_1"),
                  mod_deep_reference_ui("deep_reference_ui_1")
                )
              )),
            sortable::sortable_js("anno_data1", options = sortable::sortable_options(swap = TRUE, group = "anno_data")),
            sortable::sortable_js("anno_data2", options = sortable::sortable_options(swap = TRUE, group = "anno_data"))
          ),
          shinydashboard::tabItem(
            tabName = "biocircos_sidemenu",
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
                    shiny::checkboxInput("ShowBiocircosColoring", "Show Biocircos coloring scheme"),
                    sidebar = shinydashboardPlus::boxSidebar(
                      id = "biocircos_box_sidebar",
                      width = 25,
                      shiny::checkboxInput("biocircos_color", "Make arcs in biocircos colorful, based on the class"),
                      shiny::checkboxInput("label_color", "Make links in biocircos colorful, based on the class"),
                      shiny::selectInput("label_color_class", "Choose the mode to color the links", choices = c("Hierarchical-based" = "H",
                                                                                                                "Purity-based" = "P",
                                                                                                                "Reference column-based" = "R"
                      ),
                      selected = 'H'),
                      shiny::selectInput("ref_col_biocircos", "Choose reference column to color the links", choices = c(""), selected = '')
                    ),
                    BioCircos::BioCircosOutput("biocircos", height = "900px")%>%
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
                  shiny::uiOutput("biocircos_coloring")
                )
              )
            ),
            sortable::sortable_js("biocircos_data1", options = sortable::sortable_options(swap = TRUE, group = "biocircos_data")),
            sortable::sortable_js("biocircos_data2", options = sortable::sortable_options(swap = TRUE, group = "biocircos_data"))
          ),
          shinydashboard::tabItem(
            tabName = "summarize_sidemenu",
            shiny::fluidRow(
              tags$div(
                id="summarize_data1",
                mod_barplot_rank_ui("barplot_rank_ui_1"),
                mod_group_table_ui("group_table_ui_1")
              )
            ),
            sortable::sortable_js("summarize_data1", options = sortable::sortable_options(swap = TRUE))
          ),
          shinydashboard::tabItem(
            tabName = "uploaddata_sidemenu",
            shiny::fluidRow(
              tags$div(
                id="upload_data1",
                div(
                  id = "id1",
                  shinydashboardPlus::box(
                    title = "Upload Antismash data",
                    id = "upload_anti_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    shiny::fileInput("anti_data",
                                     "Upload Antismash data", accept = list(".csv", ".json"))
                  )
                ),
                div(
                  id = "id2",
                  shinydashboardPlus::box(
                    title = "Upload PRISM data",
                    id = "upload_prism_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    shiny::fileInput("prism_data",
                                     "Upload PRISM data", accept = list(".csv", ".json"))
                  )
                ),
                div(
                  id = "id3",
                  shinydashboardPlus::box(
                    title = "Upload SEMPI 2.0 data",
                    id = "upload_sempi_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    shiny::fileInput("sempi_data",
                                     "Upload SEMPI 2.0 data", accept = ".csv")
                  )
                ),
                div(
                  id = "id4",
                  shinydashboardPlus::box(
                    title = "Upload DeepBGC data",
                    id = "upload_deep_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    shiny::fileInput("deep_data",
                                     "Upload DeepBGC data", accept = ".tsv")
                  )
                )
              )
            ),
            shiny::fluidRow(
              tags$div(
                id="upload_data2",
                div(
                  id = "id1",
                  shinydashboardPlus::box(
                    title = "Upload Gecco data",
                    id = "upload_gecco_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    shiny::fileInput("gecco_data",
                                     "Upload Gecco data", accept = ".tsv")
                  )
                ),
                div(
                  id = "id2",
                  shinydashboardPlus::box(
                    title = "Upload RRE-Finder data",
                    id = "upload_rre_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    shiny::fileInput("rre_data",
                                     "Upload RRE-Finder data")
                  )
                ),
                div(
                  id = "id3",
                  shinydashboardPlus::box(
                    title = "Upload ARTS data",
                    id = "upload_arts_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    shiny::fileInput("known_data",
                                     "Upload ARTS knownhits data", accept = ".csv"),
                    shiny::fileInput("dup_data",
                                     "Upload ARTS duptable data", accept = ".csv")
                  )
                ),
                div(
                  id = "id4",
                  shinydashboardPlus::box(
                    title = "Use Example data",
                    id = "use_example_data_box",
                    collapsible = TRUE,                                          
                    closable = TRUE,
                    shiny::actionButton("anti_sco", "Use Antismash example data from S.coelicolor"),
                    shiny::actionButton("prism_sco", "Use PRISM example data from S.coelicolor"),
                    shiny::actionButton("sempi_sco", "Use SEMPI example data from S.coelicolor"),
                    shiny::actionButton("deep_sco", "Use DeepBGC example data from S.coelicolor"),
                    shiny::actionButton("gecco_sco", "Use Gecco example data from S.coelicolor"),
                    shiny::actionButton("rre_sco", "Use RRE-Finder example data from S.coelicolor"),
                    shiny::actionButton("arts_sco", "Use ARTS example data from S.coelicolor"),
                    shiny::numericInput("chr_len", "Please type chr len of an organism", value = 10000000)
                  )
                )
              )
            ),
            sortable::sortable_js("upload_data1", options = sortable::sortable_options(swap = TRUE, group = "upload_data")),
            sortable::sortable_js("upload_data2", options = sortable::sortable_options(swap = TRUE, group = "upload_data"))
          ),
          shinydashboard::tabItem(
            tabName = "options_sidemenu",
            shiny::fluidRow(
              shiny::column(
                width = 6,
                tags$div(
                  id="options_data1",
                  div(
                    id = "id1",
                    shinydashboardPlus::box(
                      title = "Rename",
                      id = "rename_box",
                      collapsible = TRUE,                                          
                      closable = TRUE,
                      width = NULL,
                      shiny::checkboxInput("anti_hybrid", "Visualize AntiSMASH BGC with several types as 'Hybrid'"),
                      shiny::checkboxInput("prism_hybrid", "Visualize PRISM BGC with several types as 'Hybrid'"),
                      shiny::checkboxInput("sempi_hybrid", "Visualize SEMPI BGC with several types as 'Hybrid'"),
                      shiny::fileInput("rename_data",
                                       "Upload renaming and coloring scheme", accept = ".csv"),
                      shiny::actionButton("rename", "Rename"),
                      shiny::actionButton("reset_name", "Reset")
                    )
                  ),
                  div(
                    id = "id2",
                    shiny::uiOutput("deep_filter_box")
                  ))),
              shiny::column(
                width = 6,
                tags$div(
                  id="options_data2",
                  div(
                    id = "id3",
                    shiny::uiOutput("gecco_filter_box")
                  ),
                  div(
                    id = "id5",
                    shinydashboardPlus::box(
                      title = "Improve global visualization",
                      id = "improve_visualization_box",
                      collapsible = TRUE,                                          
                      closable = TRUE,
                      width = NULL,
                      shiny::checkboxInput("rre_width", "Add thickness to RRE results visualization"),
                      shiny::checkboxInput("prism_supp_data_input_width", "Add thickness to PRISM resistance + regulatory genes results visualization"),
                      shiny::checkboxInput("arts_width", "Add thickness to ARTS results visualization"),
                      shiny::checkboxInput("sempi_width", "Add thickness to SEMPI results visualization")
                    )
                  ),
                  div(
                    id = "id4",
                    shinydashboardPlus::box(
                      title = "Prism supplement + ARTS options",
                      id = "prism_supplement_arts_box",
                      collapsible = TRUE,                                          
                      closable = TRUE,
                      width = NULL,
                      shiny::checkboxInput("prism_supp", "Visualize PRISM resistance and regulatory genes"),
                      shiny::selectInput("dup_choice", "Choose duplicated core gene to plot only it", choices = c("All"),
                                         selected = "All")
                    )
                  ),
                  mod_download_ui("download_ui_1")
                )
              )
            ),
            sortable::sortable_js("options_data1", options = sortable::sortable_options(swap = TRUE, group = "options_data")),
            sortable::sortable_js("options_data2", options = sortable::sortable_options(swap = TRUE, group = "options_data"))
          )
        )
      )
    )
    )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'BGCViz'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

