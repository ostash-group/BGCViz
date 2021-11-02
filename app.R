# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above. 
#
# Author: Pavlo Hrab
# Made as part of Cambridge bioinformatics hackaton
# 
# This app is using bgc coordinates from DeepBGC, PRISM, ANTISMASH, RRE-Finder,
# GECCO, ARTS, SEMPI to visualized interception of those different annotations 
# in one genome
#
library(magrittr)
#options(repos = BiocManager::repositories())
options(spinner.type=6)
# Define UI 
ui <- shinydashboardPlus::dashboardPage(
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
            div(
              id="anno_div_1",
              shinyjqui::jqui_resizable(
                  shinydashboardPlus::box(
                title = "Annotations reference",
                id = "annotation_reference_box",
                height = "100%",
                width = NULL,
                collapsible = TRUE,                                          
                closable = TRUE,
                plotly::plotlyOutput("deep_reference_2")  %>%
                  shinycssloaders::withSpinner()
              ), options = list(handles="w,e"))
            ),
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
                plotly::plotlyOutput("deep_reference")  %>%
                  shinycssloaders::withSpinner()
              ), options = list(handles="w,e")),
            )
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
            div(
              id="id1",
              shinyjqui::jqui_resizable(
                  shinydashboardPlus::box(
                title = "Ranking barplot",
                id = "ranking_barplot_box",
                collapsible = TRUE,                                          
                closable = TRUE,
                height = "100%",
                plotly::plotlyOutput("barplot_rank", height = "600px")%>%
                  shinycssloaders::withSpinner()
              ),options = list(handles="w,e"))
            ),
            div(
              id="id2",
              shinyjqui::jqui_resizable( 
                  shinydashboardPlus::box(
                title = "Group table",
                id = "group_table_box",
                collapsible = TRUE,                                          
                closable = TRUE,
                style='overflow-x: scroll;height:700px;overflow-y: scroll;',
                height = "100%",
                shiny::checkboxInput("count_all", "Show all BGC for the 'group by' method (+ individually annotated BGC)"),
                shiny::selectInput("group_by", "Group data by", choices = c(""),  selected = ''),
                shiny::tableOutput("group_table")%>%
                  shinycssloaders::withSpinner()
              ),options = list(handles="w,e"))
            )
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
            div(
              id = "id6",
              shinydashboardPlus::box(
                title = "Download data",
                id = "download_data_box",
                collapsible = TRUE,                                          
                closable = TRUE,
                width = NULL,
                shiny::downloadButton("download","Download currently used datasets (as for Biocircos plot)" )
              )
            )
            )
          )
        ),
        sortable::sortable_js("options_data1", options = sortable::sortable_options(swap = TRUE, group = "options_data")),
        sortable::sortable_js("options_data2", options = sortable::sortable_options(swap = TRUE, group = "options_data"))
      )
    )
  )
)
# Define server logic
server <- function(input, output, session) {
  ##---------------------------------------------------------------
  ##        Some lists of reactive values to listen later         -
  ##---------------------------------------------------------------
  check_to_rename <- shiny::reactive({list(input$sempi_data, input$anti_data, input$prism_data,
                                           input$sempi_sco,input$anti_sco, input$prism_sco)})
  biocircos_listen <- shiny::reactive({
    list( input$biocircos_color,vals$need_filter, input$label_color, input$label_color_class, 
          input$ref_col_biocircos, vals$inters_filtered,  input$prism_supp_data_input_width, vals$prism_supp_data_input,
          input$arts_width, input$sempi_width, input$rre_width, vals$anti_data, vals$sempi_data, vals$prism_data,
          vals$coloring_datatable
          
    )
  })
  inputData <- shiny::reactive({
    list( vals$sempi_data_input, vals$rre_data_input,  vals$anti_data_input, vals$prism_data_input,
          vals$known_data_input,vals$dup_data_input,  vals$prism_supp_data_input, vals$deep_data_input, vals$gecco_data_input
    )
  })
  dynamicInput <-  shiny::reactive({
    list( input$dup_choice, vals$need_filter, input$prism_supp
    )
  }) 
  deep_reference <- shiny::reactive({
    list( vals$inters_filtered, vals$rre_more, input$ref, input$arts_width, input$sempi_width, input$rre_width,
          input$prism_supp_data_input_width, vals$anti_data, vals$prism_data, vals$sempi_data)
  })
  
  to_debounce <-  shiny::reactive({
    list(  vals$cluster_type, vals$gene_filter,vals$biodomain_filter,  vals$score_c, vals$score_d, 
           vals$score_a,  vals$score_average_gecco,vals$score_cluster_gecco, vals$domains_filter_gecco, 
           vals$prot_filter_gecco
    )
  }) %>% shiny::debounce(500)
  
  # Some dataframes that are used through the app + some vectors of untercepted values
  vals <- shiny::reactiveValues(deep_data = NULL, anti_data = NULL, rre_data=NULL, prism_data=NULL, chr_len = NULL, fullness_deep = NULL,
                                biocircos_deep = NULL, deep_data_input = FALSE,tracklist = NULL, chromosomes = NULL, fullness_gecco = NULL,
                                anti_data_input = FALSE,rre_data_input = FALSE, prism_data_input = FALSE, seg_df_ref_a = NULL,
                                seg_df_ref_d = NULL,seg_df_ref_r = NULL,seg_df_ref_p = NULL, deep_data_chromo = NULL, 
                                data_upload_count = 0, anti_type=NULL, prism_type=NULL, sempi_data = NULL, sempi_data_input= FALSE,
                                sempi_type = NULL, biocircos_color = NULL, rename_data = NULL, group_by_data = NULL, 
                                rre_interact = NULL, anti_interact = NULL, prism_interact = NULL, deep_interact = NULL,  
                                sempi_interact = NULL, df_a = NULL, df_d = NULL, df_p = NULL, df_r = NULL, prism_supp = NULL,
                                prism_json = FALSE, df_s = NULL, prism_supp_interact = NULL, known_data = NULL, dup_data = NULL,
                                known_data_input = F, dup_data_input = F, arts_data=NULL, arts_data_input = F, seg_df_ref_ar = NULL,
                                df_ps = NULL, arts_interact = NULL, rre_more = FALSE, gecco_data = NULL, gecco_data_input = FALSE,
                                gecco_data_filtered = NULL, seg_df_ref_g = NULL, prism_supp_data_input = F, computed  = NULL,
                                need_filter = F, filter_data = F, choices = list(ref=NULL, group_by=NULL, ref_col_biocircos=NULL, ref_comparison_gecco=NULL, ref_comparison = NULL),
                                renamed = NULL, renaming_notification = list(), rename_y_axis = list(), can_plot_deep_ref_2 = F, can_plot_deep_ref = F,
                                can_plot_biocircos = F, can_plot_barplot_rank = F, can_plot_group_table = F
  )
  
  vals$computed <- list(
    anti=F,deep=F, gecco=F, arts=F, prism=F, sempi=F, prism_supp=F, rre=F
  )
  # Making coloring datatable
  vals$rename_data <- read.csv("rename.csv")
  rename_data <- read.csv("rename.csv")
  coloring_datatable <-data.frame( tidyr::drop_na(data.frame(cbind(as.character(rename_data$Group_color), as.character(rename_data$Color), rename_data$Hierarchy)) ))
  coloring_datatable <- coloring_datatable[!apply(coloring_datatable == "", 1, all),]
  colnames(coloring_datatable) <- c("Name", "Color", "Hierarchy")
  vals$coloring_datatable <- DT::datatable(coloring_datatable,  rownames = F, editable = "column", options = list( dom='t',ordering=F))
  # Validation
  source("src/validate_functions.R")
  # Variables, that holds data uploads boolean (so if data is present or not)
  data_uploads <- c("anti_data_input","sempi_data_input","prism_data_input","prism_supp_data_input",
                    "arts_data_input","deep_data_input","gecco_data_input","rre_data_input")
  # Universal beginings for variables, used in the app for different data
  soft_names <- c("anti","sempi","prism","prism_supp","arts","deep","gecco","rre" )
  # The Namings, meaning how to label the data on the plots
  soft_namings <- c('Antismash', 'SEMPI','PRISM', 'PRISM-Supp', 'ARTS', 'DeepBGC', 'GECCO', 'RRE-Finder')
  # Dataframes undes vals$list, that stored the data 
  data_to_use <- c( "anti_data" ,"sempi_data" , "prism_data", "prism_supp_data","arts_data_filtered","deep_data_filtered" ,"gecco_data_filtered", "rre_data")
  # Used in barplot on summarise tab + Annotation on chromosome plots 
  abbr <- c("A", "S", "P", "P-Supp", "AR", "D", "G", "RRE")
  # Used for deep reference 2 plot
  soft_datafr <- c("seg_df_ref_a", "seg_df_ref_s" , "seg_df_ref_p", "seg_df_ref_p_s", "seg_df_ref_ar", "seg_df_ref_d", 
                   "seg_df_ref_g", "seg_df_ref_r")
  
  vals$score_a <- 50
  vals$score_d <- 50
  vals$score_c <- 50
  vals$domains_filter <- 5
  vals$biodomain_filter <- 1
  vals$gene_filter <- 1
  vals$cluster_type <- 50
  vals$score_average_gecco <- 50
  vals$score_cluster_gecco <- 50
  vals$domains_filter_gecco <- 1
  vals$prot_filter_gecco <- 1
  vals$gecco_sidebar <- FALSE
  vals$deep_sidebar <- FALSE
  vals$deep_global <- FALSE
  vals$gecco_global <- FALSE
  ##----------------------------------------------------------------
  ##                        Helper functions                       -
  ##----------------------------------------------------------------
  # Need to get them to a tidyr::separate file later
  # TODO
  files_in_dir <- list.files()
  # Iterate over those files and if found "_biocircos.csv" add remove them
  for (file_names in files_in_dir) {
    if (grepl('_biocircos.csv', file_names, fixed = TRUE)) {
      file.remove(file_names)
    } 
  }
  options(shiny.maxRequestSize=100*1024^2)
  # Small function to make integers zeros
  source("src/helper_functions.R")
  
  disable_event_logic <- function(){
    vals$can_plot_deep_ref = F
    vals$can_plot_biocircos = F
    vals$can_plot_barplot_rank = F
    vals$can_plot_group_table = F
  }
  enable_event_logic <- function(){
    vals$can_plot_deep_ref = T
    vals$can_plot_biocircos = T
    vals$can_plot_barplot_rank = T
    vals$can_plot_group_table = T
  }
  
  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                        DATA INPUT PROCESSING                        ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################
  # TODO Make tidyr::separate functions for different data types. 
  # For now you just have duplicated the code. Specifically for ARTS!
  # Reading functions:
  
  read_antismash <- function(data){
    anti_data <- data
    res_validation <- validate_basic_input(anti_data)
    if (!(res_validation[[1]])){
      anti_data <- NULL
      return(NULL)
    } else{
      anti_data <- res_validation[[2]]
    }
    # Add chromosome column
    anti_data$chromosome <-  rep("A", length(anti_data$Cluster))
    # Type magic
    anti_data$Type <- stringr::str_trim(tolower(anti_data$Type))
    anti_data['Type2'] <- stringr::str_trim(tolower(anti_data$Type))
    vals$anti_type <- anti_data$Type2
    vals$anti_data <- anti_data
    # Save file
    write.csv(vals$anti_data, "anti_data.csv", row.names = F)
    vals$anti_data_input = TRUE 
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "Antismash" = "Antismash")
    vals$choices$group_by <- c(vals$choices$group_by, "Antismash" = "Antismash")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "Antismash" = "Antismash")
    vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "Antismash" = "Antismash")
    vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "Antismash" = "Antismash")
    update_ui_with_data()
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                               selected = "Antismash" )
      shiny::updateSelectInput(session, "group_by",
                               selected = "Antismash" )
      shiny::updateSelectInput(session, "ref_comparison",
                               selected = "Antismash")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                               selected =  "Antismash")
      shiny::updateSelectInput(session, "ref_comparison_gecco",
                               selected = "Antismash")
      
    }
    return(anti_data)
  }
  read_gecco <- function(data){
    # Add chromosome column
    gecco_data <- data
    
    gecco_data$chromosome <-  rep("G", length(gecco_data$type))
    # Type magic
    gecco_data$Cluster <- seq(1:length(gecco_data$chromosome))
    gecco_data$ID <- gecco_data$Cluster
    gecco_data$Type <- stringr::str_trim(tolower(gecco_data$type))
    gecco_data$Type <- gsub("polyketide", "pks", gecco_data$Type)
    gecco_data$Type <- gsub("nrp", "nrps", gecco_data$Type)
    gecco_data$Type <- gsub("unknown", "under_threshold", gecco_data$Type)
    gecco_data['Type2'] <- stringr::str_trim(tolower(gecco_data$Type))
    drop_cols <- c("alkaloid_probability" ,  "polyketide_probability", "ripp_probability",  "saccharide_probability",
                   "terpene_probability",    "nrp_probability"  , "other_probability" )
    # Read data
    gecco_data <- gecco_data %>%
      dplyr::mutate(pks=polyketide_probability, other = other_probability, nrps = nrp_probability, alkaloid = alkaloid_probability, 
                    terpene = terpene_probability, saccharide = saccharide_probability, ripp = ripp_probability) %>%
      dplyr::select(-dplyr::one_of(drop_cols))
    gecco_data$num_prot <- sapply( stringr::str_split(as.character(gecco_data$proteins), ";"), length)
    gecco_data$num_domains <- sapply( stringr::str_split(as.character(gecco_data$domains), ";"), length)
    names(gecco_data)[names(gecco_data) == "start"] <- "Start"
    names(gecco_data)[names(gecco_data) == "end"] <-  "Stop"
    vals$gecco_data <- gecco_data
    vals$gecco_data_filtered <- filter_gecco(vals$gecco_data,vals$score_cluster_gecco,vals$score_average_gecco,vals$domains_filter_gecco,vals$prot_filter_gecco)
    # Save file
    write.csv(vals$gecco_data, "gecco_data.csv", row.names = F)
    vals$gecco_data_input = TRUE 
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "GECCO" = "GECCO")
    vals$choices$group_by <- c(vals$choices$group_by, "GECCO" = "GECCO")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "GECCO" = "GECCO")
    update_ui_with_data()
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                               selected = "GECCO" )
      shiny::updateSelectInput(session, "group_by",
                               selected = "GECCO")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                               selected =  "GECCO")
      
    }
  }
  read_prism <- function(data, json=T){
    if (json==T){
      processed_data <- process_prism_json_suppl(data)
      shiny::updateCheckboxInput(inputId = "prism_supp", value = T)
      prism_data <- processed_data[[1]]
      vals$prism_supp_data_input = T
      vals$prism_supp <- processed_data[[2]]
      vals$prism_supp_data <- processed_data[[2]]
      vals$prism_json = T 
    } else {
      prism_data <- data
    }
    res_validation <- validate_basic_input(prism_data)
    if (!(res_validation[[1]])){
      prism_data <- NULL
      return(NULL)
    } else{
      prism_data <- res_validation[[2]]
    }
    vals$choices$ref <- c(vals$choices$ref, "PRISM-Supp" = "PRISM-Supp")
    vals$choices$group_by <- c(vals$choices$group_by, "PRISM-Supp" = "PRISM-Supp")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "PRISM-Supp" = "PRISM-Supp")
    update_ui_with_data()
    prism_data$Type <- stringr::str_trim(tolower(prism_data$Type))
    prism_data['Type2'] <- stringr::str_trim(tolower(prism_data$Type))
    vals$prism_data <- prism_data
    vals$prism_type <- prism_data$Type2
    
    # Add chromosome info column
    vals$prism_data$chromosome <-  rep("P", length(vals$prism_data$Cluster))
    # Add ID column (same as Cluster)
    vals$prism_data$ID <- vals$prism_data$Cluster
    # Save file
    write.csv(vals$prism_data, "prism_data.csv", row.names = F)
    vals$prism_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "PRISM" = "PRISM")
    vals$choices$group_by <- c(vals$choices$group_by, "PRISM" = "PRISM")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "PRISM" = "PRISM")
    vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "PRISM" = "PRISM")
    vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "PRISM" = "PRISM")
    update_ui_with_data()
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                               selected = "PRISM" )
      shiny::updateSelectInput(session, "group_by",
                               selected = "PRISM" )
      shiny::updateSelectInput(session, "ref_comparison",
                               selected = "PRISM")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                               selected =  "PRISM")
      shiny::updateSelectInput(session, "ref_comparison_gecco",
                               selected = "PRISM")
    }
  }
  read_sempi <- function(data){
    sempi_data <- data
    res_validation <- validate_basic_input(sempi_data)
    if (!(res_validation[[1]])){
      sempi_data <- NULL
      return(NULL)
    } else{
      sempi_data <- res_validation[[2]]
    }
    sempi_data['Type2'] <- stringr::str_trim(tolower(sempi_data$Type))
    vals$sempi_type <- sempi_data$Type2
    vals$sempi_data <- sempi_data
    # Add chromosome info column
    vals$sempi_data$chromosome <-  rep("S", length(vals$sempi_data$Cluster))
    # Add ID column (same as Cluster)
    vals$sempi_data$ID <- vals$sempi_data$Cluster
    # Save file
    write.csv(vals$sempi_data, "sempi_data.csv", row.names = F)
    vals$sempi_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "SEMPI" = "SEMPI")
    vals$choices$group_by <- c(vals$choices$group_by, "SEMPI" = "SEMPI")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "SEMPI" = "SEMPI")
    vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "SEMPI" = "SEMPI")
    vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "SEMPI" = "SEMPI")
    update_ui_with_data()
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                               selected = "SEMPI" )
      shiny::updateSelectInput(session, "group_by",
                               selected = "SEMPI" )
      shiny::updateSelectInput(session, "ref_comparison",
                               selected = "SEMPI")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                               selected =  "SEMPI")
      shiny::updateSelectInput(session, "ref_comparison_gecco",
                               selected = "SEMPI")
    }
  }
  read_arts_knownhits <- function(data){
    locations <- sapply(data$Sequence.description, function(x){
      tail(stringr::str_split(x , "\\|")[[1]], 1)
    })
    
    start <- sapply(locations, function(x){
      stringr::str_split(x, "_")[[1]][1]
    })
    stop <- sapply(locations, function(x){
      stringr::str_split(x, "_")[[1]][2]
    })
    
    known_table <- data.frame(cbind(start, stop))
    colnames(known_table) <- c("Start", "Stop")
    rownames(known_table) <- seq(1:dim(known_table)[1])
    known_table$Start <- as.numeric(known_table$Start )
    known_table$Stop <- as.numeric(known_table$Stop)
    known_table$Description <- data$Description
    known_table$Model <- data$X.Model
    known_table$Evalue <- data$evalue
    known_table$Bitscore <- data$bitscore
    known_table$ID <- seq(1:dim(known_table)[1])
    known_table$Cluster <-known_table$ID
    known_table$Type <- 'resistance'
    known_table$Type2 <- known_table$Type
    known_table$Hit <- NA
    known_table$Core <- "Not_core"
    known_table$Count <- 1
    vals$known_data <- known_table
    vals$known_data_input <- TRUE
    write.csv(vals$known_data, "knownhits_data.csv", row.names = F)
    if ((vals$dup_data_input == T)){
      vals$choices$ref <- c(vals$choices$ref, "ARTS" = "ARTS")
      vals$choices$group_by <- c(vals$choices$group_by, "ARTS" = "ARTS")
      vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "ARTS" = "ARTS")
      update_ui_with_data()
      dup_table <- vals$dup_data
      known_table <- vals$known_data
      arts_data <- rbind(dup_table, known_table)
      arts_data <- arts_data %>%
        dplyr::arrange(Start)
      arts_data$ID <- seq(1:dim(arts_data)[1])
      arts_data$Cluster <- arts_data$ID
      vals$arts_data <- arts_data
      vals$arts_data_input <- T
      vals$data_upload_count <-  vals$data_upload_count +1
      dup_table_id <- arts_data %>%
        dplyr::filter(Core != "Not_core")
      shiny::updateSelectInput(session, "dup_choice",
                               choices = c("All", paste0("ID:",dup_table_id$ID, " ,Core:", dup_table_id$Core)),
                               selected = "All" )
      vals$upl_arts = T
      if (vals$data_upload_count == 1){
        shiny::updateSelectInput(session, "ref",
                                 selected = "ARTS" )
        shiny::updateSelectInput(session, "group_by",
                                 selected = "ARTS" )
        shiny::updateSelectInput(session, "ref_col_biocircos",
                                 selected =  "ARTS")
      }
    } 
    
  }
  read_arts_dupdata <- function(data){
    get_location_duptable <- function(x, y){
      test <- stringr::str_split(x, ";")
      test2<- sub(".*loc\\|", "", test[[1]])
      test3 <- stringr::str_split(test2, " ")
      res <- list()
      for (i in seq(1:length(test3))){
        id <- paste('hit',as.character(i), sep = "_")
        start <- test3[[i]][1]
        stop <- test3[[i]][2]
        res_1 <- list(id,start, stop)
        res <- append(res, list(res_1))
      }
      return(res)
      
    }
    
    dup_table <- data.frame()
    for (i in seq(1:dim(data)[1])){
      lst <- get_location_duptable(data$X.Hits_listed.[i])
      fin_data <- data.frame(do.call("rbind", lst))
      fin_data$Core_gene <- data$X.Core_gene[i]
      fin_data$Description <- data$Description[i]
      fin_data$Count <- data$Count[i]
      colnames(fin_data) <- c("Hit", "Start", "Stop", "Core", "Description", "Count")
      dup_table <- rbind(dup_table, fin_data)
    }
    dup_table$Hit <- unlist(dup_table$Hit)
    dup_table$Start <- unlist(dup_table$Start)
    dup_table$Stop <- unlist(dup_table$Stop)
    dup_table$Start <- as.numeric(dup_table$Start )
    dup_table$Stop <- as.numeric(dup_table$Stop)
    dup_table$ID <- seq(1: dim(dup_table)[1])
    dup_table$Cluster <- dup_table$ID
    dup_table$Type <- 'core'
    dup_table$Type2 <- dup_table$Type
    dup_table$Evalue <- NA 
    dup_table$Bitscore <- NA 
    dup_table$Model <- "Core" 
    vals$dup_data <- dup_table
    vals$dup_data_input = T
    write.csv(dup_table, "duptable_data.csv", row.names = F)
    if ((vals$known_data_input == T)){
      vals$choices$ref <- c(vals$choices$ref, "ARTS" = "ARTS")
      vals$choices$group_by <- c(vals$choices$group_by, "ARTS" = "ARTS")
      vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "ARTS" = "ARTS")
      update_ui_with_data()
      dup_table <- vals$dup_data
      known_table <- vals$known_data
      arts_data <- rbind(dup_table, known_table)
      arts_data <- arts_data %>%
        dplyr::arrange(Start)
      arts_data$ID <- seq(1:dim(arts_data)[1])
      arts_data$Cluster <- arts_data$ID
      vals$arts_data <- arts_data
      vals$data_upload_count <-  vals$data_upload_count +1
      vals$arts_data_input <- T
      dup_table_id <- arts_data %>%
        dplyr::filter(Core != "Not_core")
      shiny::updateSelectInput(session, "dup_choice",
                               choices = c("All", paste0("ID:",dup_table_id$ID, " ,Core:", dup_table_id$Core)),
                               selected = "All" )
      if (vals$data_upload_count == 1){
        shiny::updateSelectInput(session, "ref",
                                 selected = "ARTS" )
        shiny::updateSelectInput(session, "group_by",
                                 selected = "ARTS" )
        shiny::updateSelectInput(session, "ref_col_biocircos",
                                 selected =  "ARTS")
      }
    } 
  }
  read_deep <- function(data){
    res_validation <- validate_deep_input(data)
    if (!(res_validation[[1]])){
      deep_data <- NULL
      return(NULL)
    } else{
      deep_data <- res_validation[[2]]
    }
    drop_cols <- c("nrp","polyketide")
    # Read data
    deep_data <- deep_data %>%
      dplyr::mutate(pks=polyketide, nrps = nrp ) %>%
      dplyr::select(-dplyr::one_of(drop_cols))
    # Add chromosome info column
    vals$deep_data <- deep_data
    vals$deep_data$chromosome <-  rep("D", length(vals$deep_data$bgc_candidate_id))
    vals$deep_data$Start <- vals$deep_data$nucl_start
    vals$deep_data$Stop <- vals$deep_data$nucl_end
    # Add ID column as number seuquence of dataframe length
    vals$deep_data$ID <- seq(1:length(vals$deep_data$bgc_candidate_id))
    vals$deep_data$Cluster <- vals$deep_data$ID
    write.csv(vals$deep_data, "deep_data.csv", row.names = F)
    vals$deep_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
    vals$deep_data_filtered <- filter_deepbgc(vals$deep_data,vals$cluster_type,vals$score_a,vals$score_c,vals$score_d,vals$domains_filter,vals$biodomain_filter,vals$gene_filter)
    vals$choices$ref <- c(vals$choices$ref, "DeepBGC" = "DeepBGC")
    vals$choices$group_by <- c(vals$choices$group_by, "DeepBGC" = "DeepBGC")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "DeepBGC" = "DeepBGC")
    update_ui_with_data()
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                               selected = "DeepBGC" )
      shiny::updateSelectInput(session, "group_by",
                               selected = "DeepBGC" )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                               choices = "DeepBGC",
                               selected = "DeepBGC")
      
    }
  }
  read_rre <- function(data){
    res_validation <- validate_rre_input(data)
    if (!(res_validation[[1]])){
      data <- NULL
      return(NULL)
    } else{
      data <- res_validation[[2]]
    }
    # Clean RRE data. Extract coordinates and Locus tag with double underscore delimiter (__)
    vals$rre_data <- data %>%
      tidyr::separate(Gene.name, c("Sequence","Coordinates","Locus_tag"),sep = "__") %>%
      tidyr::separate(Coordinates, c("Start", "Stop"),sep = "-")
    # Add chromosome info column
    vals$rre_data$chromosome <- rep("RRE",length(vals$rre_data$Sequence))
    # Add ID column
    vals$rre_data$ID <- seq(1:length(vals$rre_data$Sequence))
    vals$rre_data$Cluster <- vals$rre_data$ID
    vals$rre_data <- data.frame(vals$rre_data)
    vals$rre_data['Type'] <- 'ripp'
    vals$rre_data['Type2'] <- 'ripp'
    vals$rre_data$Start <- as.numeric(vals$rre_data$Start) 
    vals$rre_data$Stop <- as.numeric(vals$rre_data$Stop)
    # Store rre data into local variable
    vals$rre_data <- data.frame(vals$rre_data)
    write.csv(vals$rre_data, "rre_data.csv", row.names = F)
    
    vals$rre_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "RRE-Finder" = "RRE-Finder")
    vals$choices$group_by <- c(vals$choices$group_by, "RRE-Finder" = "RRE-Finder")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "RRE-Finder" = "RRE-Finder")
    update_ui_with_data()
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                               selected = "RRE-Finder" )
      shiny::updateSelectInput(session, "group_by",
                               selected = "RRE-Finder" )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                               selected = "RRE-Finder")
      
    }
    if (!is.null(vals$rre_data$Probability)){
      vals$rre_more = T
    } else {
      vals$rre_more = F
    }
  }
  
  #----------------------------------------------------------------
  ##            Loading and processing of example data             -
  ##----------------------------------------------------------------
  shiny::observeEvent(input$anti_sco,{
    
    anti_data <- read.csv("example_data/sco_antismash.csv")
    anti_data <- read_antismash(anti_data)
    
  })
  
  shiny::observeEvent(input$gecco_sco,{
    gecco_data <- read.delim("example_data/sco_gecco.tsv")
    read_gecco(gecco_data)
    
  })
  
  shiny::observeEvent(input$prism_sco,{
    # Read data
    
    data <- rjson::fromJSON(file = "example_data/sco_prism.json")
    read_prism(data)
    
  })
  
  shiny::observeEvent(input$sempi_sco,{
    sempi_data <- read.csv("example_data/sco_sempi.csv")
    read_sempi(sempi_data)
    
  })
  
  shiny::observeEvent(input$arts_sco, {
    
    data <- read.delim("example_data/sco_duptable.tsv")
    disable_event_logic()
    read_arts_dupdata(data)
    
    data <- read.delim("example_data/sco_knownhits.tsv")
    read_arts_knownhits(data)
  })
  
  shiny::observeEvent(input$deep_sco, {
    
    data <- read.delim("example_data/sco_deep.tsv") 
    read_deep(data)
  })
  
  shiny::observeEvent(input$rre_sco, {
    
    # Read data
    data <-  read.delim("example_data/sco_rre.txt")
    read_rre(data)
    
  })
  
  ##----------------------------------------------------------------
  ##                Loading and processing user data               -
  ##----------------------------------------------------------------
  shiny::observeEvent(input$anti_data,{
    
    disable_event_logic()
    # Read data
    if (input$anti_data$type=="text/csv"){
      anti_data <- read.csv(input$anti_data$datapath)
    }else{
      data <- rjson::fromJSON(file = input$anti_data$datapath)
      types <- sapply(data$records, function(y){
        lapply(y$features, function(x){
          if (unlist(x$type == 'region')){
            tolower(x$qualifiers$product)
          }
        })
      })
      
      types <-  Filter(Negate(is.null), types)
      
      types <- sapply(types, function(x){
        if (length(unlist(x))>1){
          tmp <- stringr::str_trim(paste0(unlist(x), collapse = '', sep = " "))
          gsub(" ", "__", tmp)
        }else{
          x
        }
      })
      
      location <- sapply(data$records, function(y){
        unlist(sapply(y$features, function(x){
          if (unlist(x$type == 'region')){
            unlist(x$location)
          }
        })
        )
      })
      
      
      location <- gsub("\\[", "", location)
      location <- gsub("\\]", "", location)
      location <- data.frame(location)
      colnames(location) <- "split"
      anti_data <- location %>%
        tidyr::separate(split, c("Start", "Stop")) %>%
        dplyr::transmute(ID = rownames(location), Start, Stop)
      
      anti_data <- cbind(anti_data, types)
      colnames(anti_data) <- c("Cluster", "Start", "Stop", "Type")
      anti_data$Cluster <- as.numeric(anti_data$Cluster)
      anti_data$Start <- as.numeric(anti_data$Start)
      anti_data$Stop <- as.numeric(anti_data$Stop)
      
    }
    
    read_antismash(anti_data)
    
  })
  
  shiny::observeEvent(input$sempi_data,{
    
    
    sempi_data <- read.csv(input$sempi_data$datapath)
    read_sempi(sempi_data)
    
  })
  
  shiny::observeEvent(input$gecco_data,{
    
    gecco_data <- read.delim(input$gecco_data$datapath)
    read_gecco(gecco_data)
    
  })
  
  # These are for ARTS data processing
  # input$known_data and inoput$dup_data
  shiny::observeEvent(input$known_data, {
    disable_event_logic()
    
    data <- read.delim(input$known_data$datapath)
    read_arts_knownhits(data)
  })
  
  shiny::observeEvent(input$dup_data, {
    disable_event_logic()
    
    data <- read.delim(input$dup_data$datapath)
    
    read_arts_dupdata(data)
  })
  
  shiny::observeEvent(input$prism_data,{
    
    # Read data
    if (input$prism_data$type == "text/csv"){
      prism_data <- read.csv(input$prism_data$datapath)
      read_prism(prism_data, json=F)
    } else{
      data <- rjson::fromJSON(file = input$prism_data$datapath)
      read_prism(data)
    }
    
  })
  
  shiny::observeEvent(input$deep_data, {
    
    data <- read.delim(input$deep_data$datapath)
    read_deep(data)
  })
  
  shiny::observeEvent(input$rre_data, {
    
    # Read data
    rre_data <- read.delim(input$rre_data$datapath)
    read_rre(rre_data)
  })
  
  ############################################################################
  ############################################################################
  ###                                                                      ###
  ###                INTERFACE LOGIC: WHAT TO SHOW AND WHEN                ###
  ###                                                                      ###
  ############################################################################
  ############################################################################
  # Update choices
  update_ui_with_data <- function(){
    shiny::updateSelectInput(session, "ref",
                             choices = vals$choices$ref )
    shiny::updateSelectInput(session, "group_by",
                             choices = vals$choices$group_by )
    shiny::updateSelectInput(session, "ref_col_biocircos",
                             choices = vals$choices$ref_col_biocircos )
    shiny::updateSelectInput(session, "ref_comparison_gecco",
                             choices = vals$choices$ref_comparison_gecco )
    shiny::updateSelectInput(session, "ref_comparison",
                             choices = vals$choices$ref_comparison )
  }
  # Observe input of chromosome length
  shiny::observeEvent(input$chr_len,{
    
    vals$chr_len <- input$chr_len
  })
  ##----------------------------------------------------------------
  ##    Simple options showing/hiding logic for every data input   -
  ##----------------------------------------------------------------
  # SHOW rre_width parameter if data is available
  # and hide_viz == F
  shiny::observeEvent(vals$rre_data_input, {
    
    if (vals$rre_data_input == T){
      shinyjs::showElement(selector = "#rre_width")
    } else{
      shinyjs::hideElement(selector = "#rre_width")
    }
  })
  # Show anti_hybrid option if data is available
  # And checkbox is unchecked
  shiny::observeEvent(vals$anti_data_input, {
    
    if (vals$anti_data_input == T){
      shinyjs::showElement(selector = "#anti_hybrid")
    } else{
      shinyjs::hideElement(selector = "#anti_hybrid")
    }
  })
  # Show prism options if data is available
  # If hide anti is F (checkbox), then show them
  # Only if prism_json file, then show Prism-Supp
  # And if hide_viz == F, and prism_json, then 
  # show width
  shiny::observeEvent(vals$prism_data_input, {
    if (vals$prism_data_input == T){
      shinyjs::showElement(selector = "#prism_hybrid")
      if (vals$prism_json == T){
        shinyjs::showElement(selector = "#prism_supp")
      }
      if (vals$prism_json == T){
        shinyjs::showElement(selector = "#prism_supp_data_input_width")
      }
    } else{
      shinyjs::hideElement(selector = "#prism_header")
      shinyjs::hideElement(selector = "#prism_hybrid")
      shinyjs::hideElement(selector = "#prism_supp")
      shinyjs::hideElement(selector = "#prism_supp_data_input_width")
    }
  })
  # Show SEMPI elements on data upload
  shiny::observeEvent(vals$sempi_data_input, {
    
    if (vals$sempi_data_input == T){
      shinyjs::showElement(selector = "#sempi_hybrid")
      shinyjs::showElement(selector = "#sempi_width")
    } else{
      shinyjs::hideElement(selector = "#sempi_hybrid")
      shinyjs::hideElement(selector = "#sempi_width")
    }
  })
  # Ahow ARTS data options, if data is available
  shiny::observeEvent(vals$arts_data_input,{
    
    if (vals$arts_data_input == T){
      shinyjs::showElement(selector = "#dup_choice")
      shinyjs::showElement(selector = "#arts_width")
    } else {
      shinyjs::hideElement(selector = "#dup_choice")
      shinyjs::hideElement(selector = "#arts_width")
    }
  })
  
  shiny::observeEvent(vals$data_upload_count, {
    if ((vals$arts_data_input == T) || (vals$sempi_data_input == T) || (vals$prism_supp_data_input == T) || (vals$rre_data_input == T)){
      shinyjs::showElement(selector = "#improve_visualization_box")
    } else {
      shinyjs::hideElement(selector = "#improve_visualization_box")
    }
  })
  shiny::observeEvent(vals$data_upload_count, {
    if ((vals$arts_data_input == T)  || (vals$prism_supp_data_input == T) ){
      shinyjs::showElement(selector = "#prism_supplement_arts_box")
    } else {
      shinyjs::hideElement(selector = "#prism_supplement_arts_box")
    }
  })
  ##---------------------------------------------------------------
  ##              Data processing options show/hide               -
  ##---------------------------------------------------------------
  # Count data uploads, to show tabs and corresponding 
  # options 
  
  output$deep_sidemenu_out <- shinydashboard::renderMenu({
    if (vals$data_upload_count >=2){
      if ((vals$deep_data_input == T) & ((vals$anti_data_input == T) | (vals$prism_data_input == T) | (vals$sempi_data_input == T) )) {
        shinydashboard::menuItem("Compare data with DeepBGC", tabName = "deep_sidemenu", icon = shiny::icon("dyalog"),
                                 shinydashboard::menuItem("Compare with DeepBGC plots", tabName = "deep_sidemenu", icon = shiny::icon("chart-pie")),
                                 shinydashboard::menuItem("Filtering options", tabName = "deep_filter", icon = shiny::icon("filter"),
                                                          shiny::uiOutput("deep_filter_UI_sidemenu")
                                 )
        )
        
      }
    }
  })
  output$gecco_sidemenu_out <- shinydashboard::renderMenu({
    if (vals$data_upload_count >=2){
      if ((vals$gecco_data_input == T) & ((vals$anti_data_input == T) | (vals$prism_data_input == T) | (vals$sempi_data_input == T) )){
        shinydashboard::menuItem("Compare data with GECCO", tabName = "gecco", icon = icon("fas fa-dragon"),
                                 shinydashboard::menuItem("Compare with GECCO plots", tabName = "gecco_sidemenu", icon = shiny::icon("chart-pie")),
                                 shinydashboard::menuItem("Filtering options", tabName = "gecco_filter", icon = shiny::icon("filter"),
                                                          shiny::uiOutput("gecco_filter_UI_sidemenu")
                                 ))
      }
    }
    
  })
  output$anno_sidemenu_out <- shinydashboard::renderMenu({
    if (vals$data_upload_count >=1){
      shinydashboard::menuItem("Annotation visualization and comparison", tabName = "anno_sidemenu", icon = icon("fas fa-project-diagram"))
    }
  })
  output$biocircos_sidemenu_out <- shinydashboard::renderMenu({
    if (vals$data_upload_count >=2){
      shinydashboard::menuItem("Biocircos plot", tabName = "biocircos_sidemenu", icon = icon("fas fa-circle-notch"))
    }
  })
  output$summarize_sidemenu_out <- shinydashboard::renderMenu({
    if (vals$data_upload_count >=2){
      shinydashboard::menuItem("Summarize interception", tabName = "summarize_sidemenu", icon = icon("fas fa-chart-bar"))
    }
  })
  output$biocircos_coloring <- shiny::renderUI({
    if (input$ShowBiocircosColoring == T){
      shinydashboardPlus::box(
        title = "Biocircos coloring scheme",
        closable = TRUE,
        collapsible = TRUE,                                          
        DT::dataTableOutput("biocircos_legend") %>%
          shinycssloaders::withSpinner()
      )
    }
  })
  output$deep_filter_box <- shiny::renderUI({
    if (vals$deep_data_input == T){
      vals$deep_global <- T
      shinydashboardPlus::box(
        title = "DeepBGC filtering",
        id = "deep_filtering_box",
        collapsible = TRUE,                                          
        closable = TRUE,
        width = NULL,
        shiny::sliderInput("score_a", "Activity score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
        shiny::sliderInput("score_d", "DeepBGC score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
        shiny::sliderInput("score_c", "Cluster_type score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
        # Domains, biodomains and proteins dplyr::filter. Remain >= of set threshold
        shiny::sliderInput("domains_filter", "Domain number threshold for DeepBGC data", min = 0, max = 100, value = 5),
        shiny::sliderInput("biodomain_filter", "Biodomain number threshold for DeepBGC data", min = 0, max = 100, value = 1),
        shiny::sliderInput("gene_filter", "Protein number threshold for DeepBGC data", min = 0, max = 100, value = 1),
        shiny::sliderInput("cluster_type","Choose threshold to assign cluster type for DeepBGC data ", min = 0, max = 100, value = 50)
      )
    }
  })
  output$gecco_filter_box <- shiny::renderUI({
    if (vals$gecco_data_input == T){
      vals$gecco_global <- T
      shinydashboardPlus::box(
        title = "GECCO filtering",
        id = "gecco_filtering_box",
        collapsible = TRUE,                                          
        closable = TRUE,
        width = NULL,
        shiny::sliderInput("score_average_gecco", "Average p-value threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50 ),
        shiny::sliderInput("score_cluster_gecco", "Cluster type threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50 ),
        shiny::sliderInput("domains_filter_gecco", "Domain number threshold for Gecco data", min = 0, max = 100, value = 1),
        shiny::sliderInput("prot_filter_gecco", "Protein number threshold for Gecco data", min = 0, max = 100, value = 1)
      )
    }
  })

  output$deep_filter_UI_sidemenu <- shiny::renderUI({
    vals$deep_sidebar <- T
  shiny::tagList(
    shiny::sliderInput("score_a_sidemenu", "Activity score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
    shiny::sliderInput("score_d_sidemenu", "DeepBGC score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
    shiny::sliderInput("score_c_sidemenu", "Cluster_type score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
    # Domains, biodomains and proteins dplyr::filter. Remain >= of set threshold
    shiny::sliderInput("domains_filter_sidemenu", "Domain number threshold for DeepBGC data", min = 0, max = 100, value = 5),
    shiny::sliderInput("biodomain_filter_sidemenu", "Biodomain number threshold for DeepBGC data", min = 0, max = 100, value = 1),
    shiny::sliderInput("gene_filter_sidemenu", "Protein number threshold for DeepBGC data", min = 0, max = 100, value = 1),
    shiny::sliderInput("cluster_type_sidemenu","Choose threshold to assign cluster type for DeepBGC data ", min = 0, max = 100, value = 50)
    )
  })
  output$gecco_filter_UI_sidemenu <- shiny::renderUI({
    vals$gecco_sidebar <- T
    shiny::tagList(
      shiny::sliderInput("score_average_gecco_sidemenu", "Average p-value threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50 ),
      shiny::sliderInput("score_cluster_gecco_sidemenu", "Cluster type threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50 ),
      shiny::sliderInput("domains_filter_gecco_sidemenu", "Domain number threshold for Gecco data", min = 0, max = 100, value = 1),
      shiny::sliderInput("prot_filter_gecco_sidemenu", "Protein number threshold for Gecco data", min = 0, max = 100, value = 1)
    )
  })
  
  update_filter_values <- function(listening_value, comparing_values, updating_value, rendering_check){
    if( (as.numeric(listening_value) !=  comparing_values) && (rendering_check == F)){
      shiny::updateSliderInput(session, updating_value, NULL, listening_value)
      return(list(as.numeric(listening_value),F))
    } else{
      if (grepl("sidemenu", updating_value) == T){
        shiny::updateSliderInput(session, stringr::str_split(updating_value, "_sidemenu")[[1]][1], NULL, comparing_values)
      }else{
        shiny::updateSliderInput(session, paste0(updating_value, "_sidemenu")[[1]][1], NULL, comparing_values)
      }
      return(list(comparing_values, F))
    }
  }
  
  
  observeEvent(input$score_a,{
    res <- update_filter_values(input$score_a,vals$score_a, "score_a_sidemenu", vals$deep_sidebar)
    vals$score_a <- res[[1]]
    vals$deep_sidebar <- res[[2]]
   })
  observeEvent(input$score_d,{
    res <- update_filter_values(input$score_d,vals$score_d, "score_d_sidemenu", vals$deep_sidebar )
    vals$score_d <- res[[1]]
    vals$deep_sidebar <- res[[2]]
  })
  observeEvent(input$score_c,{
    res <- update_filter_values(input$score_c,vals$score_c, "score_c_sidemenu", vals$deep_sidebar )
    vals$score_c <- res[[1]]
    vals$deep_sidebar <- res[[2]]
  })
  observeEvent(input$domains_filter,{
    res <- update_filter_values(input$domains_filter,vals$domains_filter, "domains_filter_sidemenu", vals$deep_sidebar )
    vals$domains_filter <- res[[1]]
    vals$deep_sidebar <- res[[2]]
  })
  observeEvent(input$biodomain_filter,{
    res <- update_filter_values(input$biodomain_filter,vals$biodomain_filter, "biodomain_filter_sidemenu" , vals$deep_sidebar)
    vals$biodomain_filter <- res[[1]]
    vals$deep_sidebar <- res[[2]]
  })
  observeEvent(input$gene_filter,{
    res <- update_filter_values(input$gene_filter,vals$gene_filter, "gene_filter_sidemenu", vals$deep_sidebar )
    vals$gene_filter <- res[[1]]
    vals$deep_sidebar <- res[[2]]
  })
  observeEvent(input$cluster_type,{
    res <- update_filter_values(input$cluster_type,vals$cluster_type, "cluster_type_sidemenu", vals$deep_sidebar )
    vals$cluster_type <-res[[1]]
    vals$deep_sidebar <- res[[2]]
  })
  observeEvent(input$score_a_sidemenu,{
    res <- update_filter_values(input$score_a_sidemenu,vals$score_a, "score_a", vals$deep_global )
    vals$score_a <- res[[1]]
    vals$deep_global  <- res[[2]]
  })
  observeEvent(input$score_d_sidemenu,{
    res <- update_filter_values(input$score_d_sidemenu,vals$score_d, "score_d" , vals$deep_global )
    vals$score_d <- res[[1]]
    vals$deep_global  <- res[[2]]
  })
  observeEvent(input$score_c_sidemenu,{
    res <- update_filter_values(input$score_c_sidemenu,vals$score_c, "score_c" , vals$deep_global )
    vals$score_c <- res[[1]]
    vals$deep_global  <- res[[2]]
  })
  observeEvent(input$domains_filter_sidemenu,{
    res <- update_filter_values(input$domains_filter_sidemenu,vals$domains_filter, "domains_filter" , vals$deep_global )
    vals$domains_filter <- res[[1]]
    vals$deep_global  <- res[[2]]
  })
  observeEvent(input$biodomain_filter_sidemenu,{
    res <- update_filter_values(input$biodomain_filter_sidemenu,vals$biodomain_filter, "biodomain_filter" , vals$deep_global )
    vals$biodomain_filter <- res[[1]]
    vals$deep_global  <- res[[2]]
  })
  observeEvent(input$gene_filter_sidemenu,{
    res <- update_filter_values(input$gene_filter_sidemenu,vals$gene_filter, "gene_filter", vals$deep_global  )
    vals$gene_filter <- res[[1]]
    vals$deep_global  <- res[[2]]
  })
  observeEvent(input$cluster_type_sidemenu,{
    res <- update_filter_values(input$cluster_type_sidemenu,vals$cluster_type, "cluster_type", vals$deep_global  )
    vals$cluster_type <- res[[1]]
    vals$deep_global  <- res[[2]]
  })
  
  
  
  observeEvent(input$score_average_gecco,{
    res <- update_filter_values(input$score_average_gecco,vals$score_average_gecco, "score_average_gecco_sidemenu",vals$gecco_sidebar )
    vals$score_average_gecco <- res[[1]]
    vals$gecco_sidebar <- res[[2]]
  })
  observeEvent(input$score_cluster_gecco,{
    res <- update_filter_values(input$score_cluster_gecco,vals$score_cluster_gecco, "score_cluster_gecco_sidemenu",vals$gecco_sidebar )
    vals$score_cluster_gecco <- res[[1]]
    vals$gecco_sidebar <- res[[2]]
  })
  observeEvent(input$domains_filter_gecco,{
    res <- update_filter_values(input$domains_filter_gecco,vals$domains_filter_gecco, "domains_filter_gecco_sidemenu" ,vals$gecco_sidebar)
    vals$domains_filter_gecco <- res[[1]]
    vals$gecco_sidebar <- res[[2]]
  })
  observeEvent(input$prot_filter_gecco,{
    res <- update_filter_values(input$prot_filter_gecco,vals$prot_filter_gecco, "prot_filter_gecco_sidemenu" ,vals$gecco_sidebar)
    vals$prot_filter_gecco <- res[[1]]
    vals$gecco_sidebar <- res[[2]]
  })
  observeEvent(input$score_average_gecco_sidemenu,{
    res <- update_filter_values(input$score_average_gecco_sidemenu,vals$score_average_gecco, "score_average_gecco",vals$gecco_global )
    vals$score_average_gecco <- res[[1]]
    vals$gecco_global <- res[[2]]
  })
  observeEvent(input$score_cluster_gecco_sidemenu,{
    res <- update_filter_values(input$score_cluster_gecco_sidemenu,vals$score_cluster_gecco, "score_cluster_gecco",vals$gecco_global )
    vals$score_cluster_gecco <- res[[1]]
    vals$gecco_global <- res[[2]]
  })
  observeEvent(input$domains_filter_gecco_sidemenu,{
    res <- update_filter_values(input$domains_filter_gecco_sidemenu,vals$domains_filter_gecco, "domains_filter_gecco",vals$gecco_global )
    vals$domains_filter_gecco <- res[[1]]
    vals$gecco_global <- res[[2]]
  })
  observeEvent(input$prot_filter_gecco_sidemenu,{
    res <- update_filter_values(input$prot_filter_gecco_sidemenu,vals$prot_filter_gecco, "prot_filter_gecco",vals$gecco_global )
    vals$prot_filter_gecco <- res[[1]]
    vals$gecco_global <- res[[2]]
  })
  
  shiny::observeEvent(input$restore_box,{
    box_ids <- c("deep_comparison_box", "deep_rate_box","deep_comparison_controls_box","gecco_comparison_box",
                 "gecco_rate_box","gecco_comparison_controls_box","annotation_reference_box","annotation_reference_comparison_box",
                 "annotation_reference_comparison_controls_box","biocircos_plot_box","biocircos_controls_box",
                 "ranking_barplot_box","group_table_box","upload_anti_box","upload_prism_box",
                 "upload_sempi_box","upload_deep_box","upload_gecco_box","upload_rre_box","upload_arts_box",
                 "use_example_data_box","rename_box","prism_supplement_arts_box", "improve_visualization_box",
                 "download_data_box","gecco_filtering_box","deep_filtering_box")
    for (id in box_ids){
      shinydashboardPlus::updateBox(id, action = "restore")
    }
  })
  
  
  # Logic show/hide selectinput in Link coloring in
  # Biocircos
  shiny::observeEvent(input$label_color_class, {
    
    if (input$label_color_class == "R"){
      shinyjs::showElement(selector = "#ref_col_biocircos")
    } else {
      shinyjs::hideElement(selector = "#ref_col_biocircos")
    }
  })
  # Make hybrids from the data, if checkbox is checked   
  # TODO Put the function to the root. 
  # Tou have duplicated code
  shiny::observeEvent(input$anti_hybrid, ignoreInit=T,{
    
    if (input$anti_hybrid==T){
      vals$anti_data$Type2 <- hybrid_col(vals$anti_data)
    }else {
      vals$anti_data$Type2 <- vals$anti_type
    }
    
  })
  shiny::observeEvent(input$prism_hybrid,ignoreInit=T, {
    
    if (input$prism_hybrid==T){
      vals$prism_data$Type2 <- hybrid_col(vals$prism_data)
    }else {
      vals$prism_data$Type2 <- vals$prism_type
    }
  })
  shiny::observeEvent(input$sempi_hybrid,  ignoreInit=T,{
    
    if (input$sempi_hybrid==T){
      vals$sempi_data$Type2 <- hybrid_col(vals$sempi_data)
    }else {
      vals$sempi_data$Type2 <- vals$sempi_type
    }
  })
  # Rename the data, if button is clicked
  shiny::observeEvent(input$rename, {
    
    rename_data <- vals$rename_data
    if (vals$anti_data_input == T){
      anti_data <- read.csv("anti_data.csv")
      res <- rename_vector(anti_data, rename_data, vals$renaming_notification)
      vals$anti_type <- res[[1]]
      vals$renaming_notification <- res[[2]]
      anti_data['Type2'] <- vals$anti_type
      vals$anti_data <- anti_data
    }
    
    if (vals$sempi_data_input == T){
      sempi_data <- read.csv("sempi_data.csv")
      res <- rename_vector(sempi_data, rename_data, vals$renaming_notification)
      vals$sempi_type <- res[[1]]
      vals$renaming_notification <- res[[2]]
      sempi_data['Type2'] <- vals$sempi_type
      vals$sempi_data <- sempi_data
    }
    
    if(vals$prism_data_input == T){
      prism_data <- read.csv("prism_data.csv")
      res <- rename_vector(prism_data, rename_data, vals$renaming_notification)
      vals$prism_type <- res[[1]]
      vals$renaming_notification <- res[[2]]
      prism_data['Type2'] <-  vals$prism_type
      vals$prism_data <- prism_data
    }
    shinyjs::showElement(selector = "#reset_name")
    shinyjs::hideElement(selector = "#rename")
    vals$renamed <- T
    shiny::showNotification(paste("Please note: SEMPI, PRISM and Antismash input data will be renamed on upload"), type = "warning", duration=10)
  })
  # When the new data is uploaded and renamed
  # is T, then rename data on upload
  shiny::observeEvent(check_to_rename(), {
    
    shiny::req(vals$renamed == T)
    
    rename_data <- vals$rename_data
    if (vals$anti_data_input == T){
      anti_data <- read.csv("anti_data.csv")
      res <- rename_vector(anti_data, rename_data, vals$renaming_notification)
      vals$anti_type <- res[[1]]
      vals$renaming_notification <- res[[2]]
      anti_data['Type2'] <- vals$anti_type
      vals$anti_data <- anti_data
    }
    
    if (vals$sempi_data_input == T){
      sempi_data <- read.csv("sempi_data.csv")
      res <- rename_vector(sempi_data, rename_data, vals$renaming_notification)
      vals$sempi_type <- res[[1]]
      vals$renaming_notification <- res[[2]]      
      sempi_data['Type2'] <- vals$sempi_type
      vals$sempi_data <- sempi_data
    }
    
    if(vals$prism_data_input == T){
      prism_data <- read.csv("prism_data.csv")
      res <- rename_vector(prism_data, rename_data, vals$renaming_notification)
      vals$prism_type <- res[[1]]
      vals$renaming_notification <- res[[2]]      
      prism_data['Type2'] <-  vals$prism_type
      vals$prism_data <- prism_data
    }
  })
  # Reset the renaming. Uncheck the hybrid checkboxes
  shiny::observeEvent(input$reset_name, {
    
    vals$anti_data['Type2']  <- vals$anti_data['Type']
    vals$sempi_data['Type2'] <- vals$sempi_data['Type']
    vals$ prism_data['Type2'] <- vals$ prism_data['Type']
    if (input$anti_hybrid==T){
      shiny::showNotification(paste("Antismash cluster types are NOT visualized as hybrid anymore. You should check the option one more time"), type = "warning", duration=10)
      shiny::updateCheckboxInput(inputId = "anti_hybrid", value = F)
    }
    if (input$prism_hybrid==T){
      shiny::showNotification(paste("PRISM cluster types are NOT visualized as hybrid anymore. You should check the option one more time"), type = "warning", duration=10)
      shiny::updateCheckboxInput(inputId = "prism_hybrid", value = F)
    }
    if (input$sempi_hybrid==T){
      shiny::showNotification(paste("SEMPI cluster types are NOT visualized as hybrid anymore. You should check the option one more time"), type = "warning", duration=10)
      shiny::updateCheckboxInput(inputId = "sempi_hybrid", value = F)
    }
    shinyjs::showElement(selector = "#rename")
    shinyjs::hideElement(selector = "#reset_name")
    vals$renamed <- F
  })
  # Read the uploaded renaming scheme csv
  shiny::observeEvent(input$rename_data,{
    
    rename_data <- read.csv(input$rename_data$datapath)
    vals$rename_data <- rename_data
    coloring_datatable <-data.frame( tidyr::drop_na(data.frame(cbind(as.character(rename_data$Group_color), as.character(rename_data$Color), rename_data$Hierarchy)) ))
    coloring_datatable <- coloring_datatable[!apply(coloring_datatable == "", 1, all),]
    colnames(coloring_datatable) <- c("Name", "Color", "Hierarchy")
    vals$coloring_datatable <- DT::datatable(coloring_datatable,  rownames = F, editable = "column")
  })
  
  
  # What to do, if hide DeepBGC comparison options scheme is triggered
  
  
  ############################################################################
  ############################################################################
  ###                                                                      ###
  ###                             COMPUTATIONS                             ###
  ###                                                                      ###
  ############################################################################
  ############################################################################
  shiny::observeEvent(input$prism_supp, ignoreInit = T,priority = 3,{
    if (input$prism_supp == T){
      vals$prism_supp_data_input = T
      vals$need_filter <- T
      if (!("PRISM-Supp" %in% names(vals$choices$ref))){
        vals$choices$ref <- c(vals$choices$ref, "PRISM-Supp" = "PRISM-Supp")
        vals$choices$group_by <- c(vals$choices$group_by, "PRISM-Supp" = "PRISM-Supp")
        vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "PRISM-Supp" = "PRISM-Supp")
        update_ui_with_data()
      }
    } else {
      vals$prism_supp_data_input = F
      vals$need_filter <- T
      vals$choices$ref <- vals$choices$ref[!(names(vals$choices$ref)%in%c("PRISM-Supp"))]
      vals$choices$group_by <- vals$choices$group_by[!(names(vals$choices$group_by)%in%c("PRISM-Supp"))]
      vals$choices$ref_col_biocircos <- vals$choices$ref_col_biocircos[!(names(vals$choices$ref_col_biocircos)%in%c("PRISM-Supp"))]
      update_ui_with_data()
    }
  })
  
  # Compute all interceptions on data upload.
  # dplyr::filter while ploting then.
  shiny::observeEvent(inputData(), ignoreInit = T,priority = 5,{
    # GENERATE DATA
    if (vals$anti_data_input == TRUE){
      anti_data <-  vals$anti_data
      anti_inter <- vals$anti_data %>%
        dplyr::select(Start, Stop)
      anti_inter$seqnames <- "chr"
      
    }
    if (vals$deep_data_input == TRUE){
      deep_data <- vals$deep_data
      deep_inter <- vals$deep_data %>% 
        dplyr::select(Start, Stop)
      
      deep_inter$seqnames <- "chr"
    }
    if (vals$rre_data_input == TRUE){
      # Convert numeric columns in a dataframe as a numeric
      vals$rre_data$Start <- as.numeric(vals$rre_data$Start) 
      vals$rre_data$Stop <- as.numeric(vals$rre_data$Stop)
      # Store rre data into local variable
      rre_data <- data.frame(vals$rre_data)
      # Start/Stop columns from rre data as matrix
      rre_inter <- rre_data %>%
        dplyr::select(Start, Stop)
      rre_inter$seqnames <- "chr"
    }
    if (vals$prism_data_input == TRUE){
      # Store master prism data in local variable
      prism_data <- vals$prism_data
      # Start/Stop columns from prism data as matrix
      prism_inter <- prism_data %>%
        dplyr::select(Start,Stop)
      prism_inter$seqnames <- "chr"
    }
    if (vals$sempi_data_input == TRUE){
      # Store master prism data in local variable
      sempi_data <- vals$sempi_data
      # Start/Stop columns from prism data as matrix
      sempi_inter <- vals$sempi_data %>%
        dplyr::select(Start,Stop)
      sempi_inter$seqnames <- "chr"
    }
    if (vals$prism_supp_data_input == T){
      prism_supp_data <- vals$prism_supp_data
      prism_supp_inter <- vals$prism_supp_data %>%
        dplyr::select(Start,Stop)
      prism_supp_inter$seqnames <- "chr"
    }
    if (vals$arts_data_input == T){
      arts_data <- vals$arts_data
      arts_inter <- vals$arts_data %>%
        dplyr::select(Start,Stop) 
      arts_inter$seqnames <- "chr"
    }
    if (vals$gecco_data_input == TRUE){
      gecco_data <- vals$gecco_data
      # Start/Stop columns from prism data as matrix
      gecco_inter <- vals$gecco_data %>%
        dplyr::select(Start,Stop)
      gecco_inter$seqnames <- "chr"
    }
    
    get_inter <- function(inter1, inter2){
      query <- GenomicRanges::makeGRangesFromDataFrame(inter2)
      subject <- GenomicRanges::makeGRangesFromDataFrame(inter1)
      interseption <- GenomicRanges::findOverlaps(query,subject)
      inter_from <- interseption@from
      inter_to <- interseption@to
      return(list(from = inter_from, to = inter_to))
    }
    
    inters <- vals$inters
    index = 1
    for (i in data_uploads){
      index_2 = 1
      j = soft_names[index]
      for (p in data_uploads){
        x = soft_names[index_2]
        if ((vals[[i]] == TRUE) & (vals$computed[[j]]==F) & (j!= x)){
          if ((vals[[p]] == TRUE) & (j != soft_names[index_2])){
            res <- get_inter( eval(as.name(paste(j, '_inter', sep = ""))), eval(as.name(paste(x, '_inter', sep = "")))) 
            new_res <- list()
            new_res$from <- eval(as.name(paste(x, '_data', sep = "")))[res$from,]$Cluster
            new_res$to <- eval(as.name(paste(j, '_data', sep = "")))[res$to,]$Cluster
            inters[[j]][[x]] <- new_res
            inters[[x]][[j]] <- list(from=new_res$to, to=new_res$from)
          }
          
        }
        index_2 = index_2 +1
      }
      if (vals[[i]] == TRUE){
        vals$computed[[j]] <- TRUE
      }
      index = index +1 
    }
    
    vals$inters <- inters
    if ((vals$deep_data_input == F) & (vals$gecco_data_input == F) &(vals$arts_data_input==F)){
      vals$inters_filtered <- inters 
      enable_event_logic()
    } else{
      vals$need_filter <- T
      vals$filter_data <- T
    }
    
  })
  # dplyr::filter ARTS, DeepBGC, GECCO interception data
  # and general dataframes to plot, if data filtering 
  # options are triggered
  shiny::observeEvent({dynamicInput()
    to_debounce()
  }, ignoreInit = T, priority = 4 ,{
    shiny::req(vals$data_upload_count>=1)
    inters <- vals$inters
    if (vals$deep_data_input == TRUE){
      if (vals$need_filter == F) {
        biocircos_deep <- filter_deepbgc(vals$deep_data,vals$cluster_type,vals$score_a,vals$score_c,vals$score_d,vals$domains_filter,vals$biodomain_filter,vals$gene_filter)
        vals$deep_data_filtered <- biocircos_deep
      } else {
        biocircos_deep <-  vals$deep_data_filtered
      }
      if (vals$data_upload_count!=1){
        new_deep <- lapply(inters$deep, function(x){
          new_to <- x$to[x$to %in% biocircos_deep$Cluster]
          new_from <- x$from[x$to %in% biocircos_deep$Cluster]
          list(from=new_from, to=new_to)
        })
        new_inters <- inters
        update_list <- names(inters$deep)
        for (b in seq(1:length(update_list))){
          new_inters[[update_list[b]]]$deep$to <- new_deep[[update_list[b]]]$from
          new_inters[[update_list[b]]]$deep$from <- new_deep[[update_list[b]]]$to
        }
        new_inters$deep <- new_deep
        vals$inters_filtered <- new_inters
        inters <- new_inters
      }
    }
    if (vals$gecco_data_input == TRUE){
      if (vals$need_filter == F) {
        gecco_data <- filter_gecco(vals$gecco_data,vals$score_cluster_gecco,vals$score_average_gecco,vals$domains_filter_gecco,vals$prot_filter_gecco)
        vals$gecco_data_filtered <- gecco_data 
      } else {
        gecco_data <- vals$gecco_data_filtered
      }
      if (vals$data_upload_count!=1){
        new_gecco <- lapply(inters$gecco, function(x){
          new_to <- x$to[x$to %in% gecco_data$Cluster]
          new_from <- x$from[x$to %in% gecco_data$Cluster]
          list(from=new_from, to=new_to)
        })
        new_inters <- inters
        update_list <- names(inters$gecco)
        for (b in seq(1:length(update_list))){
          new_inters[[update_list[b]]]$gecco$to <- new_gecco[[update_list[b]]]$from
          new_inters[[update_list[b]]]$gecco$from <- new_gecco[[update_list[b]]]$to
        }
        new_inters$gecco <- new_gecco
        vals$inters_filtered <- new_inters
        inters <- new_inters
      }
    }
    if (vals$arts_data_input == TRUE){
      if (input$dup_choice != "All"){
        vals$arts_data_filtered <- data.frame(vals$arts_data) %>%
          dplyr::filter(Core == stringr::str_split(stringr::str_split(input$dup_choice, " ,")[[1]][[2]], "Core:")[[1]][[2]] | Core == "Not_core")
        if (vals$data_upload_count!=1){
          new_arts <- lapply(inters$arts, function(x){
            new_to <- x$to[x$to %in% vals$arts_data_filtered$Cluster]
            new_from <- x$from[x$to %in% vals$arts_data_filtered$Cluster]
            list(from=new_from, to=new_to)
          })
          new_inters <- inters
          update_list <- names(inters$arts)
          for (b in seq(1:length(update_list))){
            new_inters[[update_list[b]]]$arts$to <- new_arts[[update_list[b]]]$from
            new_inters[[update_list[b]]]$arts$from <- new_arts[[update_list[b]]]$to
          }
          new_inters$arts <- new_arts
          vals$inters_filtered <- new_inters
          inters <- new_inters
        }
      } else {
        vals$arts_data_filtered <- vals$arts_data
        vals$inters_filtered <- inters
        
      }
    }
    if (input$prism_supp == FALSE){
      inters$prism_supp <- NULL
      for (name in names(inters)) {
        inters[[name]][which(names(inters[[name]])%in%c("prism_supp"))]<-NULL
      }
    }
    if ((vals$gecco_data_input == F) & (vals$deep_data_input == F )& (vals$arts_data_input == F )) {
      vals$inters_filtered <- inters
    }
    vals$need_filter <- F
    vals$filter_data <- F
    vals$can_plot_deep_ref = T
    enable_event_logic()
    
  })
  # Compute the Biociros plot. Store information to plot later
  shiny::observeEvent(biocircos_listen(),  ignoreInit = T,priority = 3 ,{
    shiny::req(vals$data_upload_count >=2)
    shiny::req(vals$need_filter == F)
    shiny::req(vals$can_plot_biocircos == T)
    source("src/biocircos_functions.R")
    #BioCircos!
    Biocircos_chromosomes <- list()
    arcs_chromosomes <- c()
    arcs_begin <- c()
    arcs_end <- c()
    arc_labels <- c()
    arc_col <- c()
    
    if (is.null(vals$inters_filtered)){
      inters <- vals$inters
    } else {
      inters <- vals$inters_filtered
    }
    
    
    rename_data <- vals$rename_data
    coloring_datatable <- vals$coloring_datatable
    
    index <- 1
    # browser()
    for (upload in data_uploads){
      if (vals[[upload]] == T){
        # Store data in local variable
        corrected_data <- correct_width(vals[[data_to_use[index]]], soft_namings[index], input$sempi_width, input$prism_supp_data_input_width, input$arts_width, input$rre_width)
        init_data <- initialize_biocircos(corrected_data, soft_namings[index], Biocircos_chromosomes, arcs_chromosomes, arcs_begin , arcs_end, arc_labels, arc_col, rename_data, vals$chr_len, input$biocircos_color , coloring_datatable)
        #Make chromosome list for Biocircos plot. Use chr_len as an input
        Biocircos_chromosomes <- init_data[[1]]
        #Add arcs. Quantity of arcs is length of dataframes
        arcs_chromosomes <- init_data[[2]]
        # Add arcs begin positions. (Start column)
        arcs_begin <- init_data[[3]]
        # Stop position of arcs. 
        arcs_end <- init_data[[4]]
        # Add Arcs labels. Can add only one label...
        arc_labels <- init_data[[5]]
        
        arc_col <- init_data[[6]]
      }
      index <- index +1 
    }
    # Add to tracklist. Then it can be populated with links
    tracklist <- BioCircos::BioCircosArcTrack('myArcTrack', arcs_chromosomes, arcs_begin, arcs_end, 
                                              minRadius = 0.90, maxRadius = 0.97, labels = arc_labels,colors = arc_col )
    # Function to get interception between two matrices. Returns a list of two elements - IDs from first matrix and 
    # from second one. IDs are duplicated, if intercepted more than one time
    
    chromosomes_start <- c()
    chromosomes_end <- c()
    link_pos_start <- c()
    link_pos_start_1 <- c()
    link_pos_end <- c()
    link_pos_end_2 <- c()
    label_1 <- c()
    label_2 <- c()
    label_color <- c()
    
    #CALCULATIONS
    # -----------------------------------------
    
    
    data_uploads_2 <- data_uploads
    soft_2 <- soft_namings
    soft_names_2 <- soft_names
    data_to_use_2 <- data_to_use
    index <- 1
    for (upload in data_uploads){
      data_uploads_2 <- data_uploads_2[-1]
      soft_2 <- soft_2[-1]
      soft_names_2 <- soft_names_2[-1]
      data_to_use_2 <- data_to_use_2[-1]
      index2 <- 1
      if (vals[[upload]] == T){
        for (upload2 in data_uploads_2){
          if ((vals[[upload2]]==T) & (length(data_uploads_2) > 0) & (soft_namings[index] != soft_2[index2])){
            output <- add_biocircos_data(inters[[soft_names[index]]][[soft_names_2[index2]]]$from, inters[[soft_names[index]]][[soft_names_2[index2]]]$to, vals[[data_to_use_2[index2]]], vals[[data_to_use[index]]], soft_2[index2], soft_namings[index], rename_data, input$label_color_class, input$ref_col_biocircos, coloring_datatable)
            
            chromosomes_start <- c(chromosomes_start, output[[3]])
            # Add link end. Just populate second output from the vectors, used above. 
            chromosomes_end <- c(chromosomes_end, output[[4]] )
            # Add links start positions as a start from dataframe. This vector is for chromosome start
            link_pos_start <- as.numeric(c(link_pos_start, output[[5]]))
            # Add links start positions as a start from dataframe. For chromosome start variable
            link_pos_start_1 <- as.numeric(c(link_pos_start_1, output[[6]]))
            # Add links start position for a chromosome stop variable
            link_pos_end <- as.numeric(c(link_pos_end, output[[7]]))
            # Add links start position for a chromosome stop position
            link_pos_end_2 <- as.numeric(c(link_pos_end_2,output[[8]]))
            label_1 <- c(label_1, output[[9]])
            label_2 <- c(label_2, output[[10]])
            label_color = c(label_color, output[[11]] )
          }
          index2 <- index2 +1
        }
        write.csv(vals[[data_to_use[index]]], paste0(soft_names[index], "_biocircos.csv"), row.names = F)
      }
      index <- index +1 
    }
    
    
    
    
    # Combine labels with mapply to one list
    link_labels <- mapply(function(x,y)  paste(x, y, sep = " | "), label_1, label_2 )
    
    # Add links and labels to the track list for subsequent visualization 
    if ((input$label_color == T) & (length(chromosomes_start) > 0)){
      group_colors <- plyr::count(unlist(label_color))
      for (i in seq(1:dim(group_colors)[1])){
        subset <- unname( which(label_color %in% group_colors$x[i]))
        tracklist = tracklist + BioCircos::BioCircosLinkTrack(as.character(i), chromosomes_start[subset], link_pos_start[subset], 
                                                              link_pos_start_1[subset], chromosomes_end[subset], link_pos_end[subset], 
                                                              link_pos_end_2[subset], maxRadius = 0.85, labels = link_labels[subset],
                                                              displayLabel = FALSE, color = group_colors$x[i])
      }
    } else if ((input$label_color == F) & (length(chromosomes_start) > 0)){
      tracklist = tracklist + BioCircos::BioCircosLinkTrack('myLinkTrack_master', chromosomes_start, link_pos_start, 
                                                            link_pos_start_1, chromosomes_end, link_pos_end, 
                                                            link_pos_end_2, maxRadius = 0.85, labels = link_labels,
                                                            displayLabel = FALSE, color = coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == 'base'])
    } else{
      shiny::showNotification(paste("No interceptions are being made in the Biocircos plot. Please provide data with clusters that do have intercepting borders"), type = "warning", duration=NULL)
    }
    
    vals$tracklist <- tracklist
    vals$Biocircos_chromosomes <- Biocircos_chromosomes
  })
  
  shiny::observeEvent(deep_reference(), ignoreInit = T,{
    shiny::req(vals$data_upload_count >=1)
    shiny::req(vals$need_filter == F)
    shiny::req(vals$can_plot_deep_ref == T)
    shiny::req(input$ref != "")
    shiny::req(vals$data_upload_count >=1)
    
    if (is.null(vals$inters_filtered)){
      inters <- vals$inters
    } else {
      inters <- vals$inters_filtered
    }
    source("src/deep_reference_functions.R")
    # GENERATE DATA
    index <- 1
    for (upload in data_uploads){
      if (vals[[upload]] == T){
        data<- vals[[data_to_use[index]]]
        assign(paste0(soft_names[index], "_data"),  correct_width(data, soft_namings[index],input$sempi_width,input$prism_supp_data_input_width,input$arts_width,input$rre_width))
      }
      index <- index +1
    }
    
    
    lett <- rev(LETTERS)[1:9]
    
    tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                "P_value", "RRE_start","RRE_stop", "Probability", "Name", "Full_name",  "Hit", "Core", "Count", "Bitscore", "Model",
                "Num_domains", "Num_proteins", "Average_p", "Max_p")
    
    
    
    
    # MAKE COMPUTATIONS
    sup_index <- 1
    soft_lttrs <- lett
    rename_y_axis <- vals$rename_y_axis
    rename_y_axis <- lapply(1:(length( soft_lttrs)-1), function(x){
      soft_lttrs[x]=soft_namings[x]
    })
    names(rename_y_axis) <- soft_lttrs[-length(soft_lttrs)]
    for (upload in data_uploads){
      soft_lttr <- soft_lttrs[1]
      soft_lttrs <- soft_lttrs[-1]
      if (vals[[upload]] == T){
        soft_major <- soft_names[sup_index]
        seg_ref_g <- simple_seg(eval(as.name(paste(soft_names[sup_index], "_data", sep = ""))), "Z", soft_namings[sup_index], soft_names[sup_index],soft_major, inter = F,inters)
        seg_ref_g <- define_spec_seg_df(soft_names, sup_index,seg_ref_g, soft_major, eval(as.name(paste(soft_names[sup_index], "_data", sep = ""))), inter = F, vals$rre_more,inters)
        seg_ref <- seg_ref_g
        
        if (input$ref == soft_namings[sup_index]){
          shiny::validate(need(nrow(eval(as.name(paste(soft_names[sup_index], "_data", sep = ""))))>0,"Reference data is empty, and so, insufficient for plotting. Please select another one") )
          plot <- ggplot2::ggplot(eval(as.name(paste(soft_names[sup_index], "_data", sep = ""))), ggplot2::aes(x = vals$chr_len, y = Chr)) + 
            eval(as.name(paste0("geom_", soft_names[sup_index])))(seg_ref,vals$rre_more)
          soft_let <- abbr[sup_index]
          lettrs <- lett[2:length(lett)]
          labels_1 <- list()
          index = 1
          for (i in data_uploads){
            if ((vals[[i]] == T) & (soft_names[index] != soft_major)){
              df <- eval(as.name(paste(soft_names[index], "_data", sep = "")))
              seg_df <- simple_seg(df, lettrs[index], soft_namings[index], soft_names[index],soft_major,inter = T, inters)
              seg_df <- define_spec_seg_df(soft_names, index,seg_df, soft_major, df,inter = T, vals$rre_more, inters)
              labels_1[[lettrs[index]]] <- (paste(abbr[index], "_vs_", soft_let, sep = ""))
              plot <- add_more_annot(seg_df, plot, soft_names, index, vals$rre_more)
            }
            index = index +1
            
          }
          plot <- plot +
            ggplot2::scale_y_discrete(labels = c("Z" = input$ref, unlist(labels_1))) +
            ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10)) +
            ggplot2::ylab("")+
            ggplot2::xlab("Chromosome length")+ 
            ggplot2::theme(legend.title = ggplot2::element_blank()) +
            ggplot2::ggtitle("Annotations' comparison to the reference")
          to_plot <- plotly::ggplotly(plot, tooltip = tooltip)
          to_plot <- to_plot %>% 
            plotly::layout(legend=list(font = list(
              family = "sans-serif",
              size = 12,
              color = "#000"),
              bordercolor = "#FFFFFF",
              borderwidth = 2,
              title=list(text='<b> Cluster Types </b>')))
          
          
        }
        seg_ref$yend <- rep(soft_lttr, length(eval(as.name(paste(soft_names[sup_index], "_data", sep = "")))$Cluster))
        seg_ref$y <- rep(soft_lttr, length(eval(as.name(paste(soft_names[sup_index], "_data", sep = "")))$Cluster))
        vals[[soft_datafr[sup_index]]] <- seg_ref
      }
      sup_index <- sup_index +1
    }
    vals$rename_y_axis <- rename_y_axis
    vals$deep_reference_to_plot <- to_plot
  })
  
  
  ############################################################################
  ############################################################################
  ###                                                                      ###
  ###                             OUTPUT PLOTS                             ###
  ###                                                                      ###
  ############################################################################
  ############################################################################
  
  ##----------------------------------------------------------------
  ##                    DeepBGC Comparison tab                     -
  ##----------------------------------------------------------------
  # Render barplot
  output$deep_barplot <- shiny::renderPlot({
    shiny::req((vals$deep_data_input == T) & ((vals$anti_data_input == T) | (vals$prism_data_input == T) | (vals$sempi_data_input == T) ))
    
    
    # Create empty dataframe to populate later
    fullnes_of_annotation <- data.frame(NA, NA, NA)
    colnames(fullnes_of_annotation) <- c("Score", "Source", "Quantity")
    fullnes_of_annotation <- tidyr::drop_na(fullnes_of_annotation)
    
    deep_inter_1 <- vals$deep_data_filtered
    # Decide which score to use for basic thresholds on x axis
    if (input$score_type == "Activity") {
      score <- "score_a"
    } else if (input$score_type == "DeepBGC") {
      score <- "score_d"
    } else if (input$score_type == "Cluster_Type") {
      score <- "score_c"
    }
    deep_inter_1$score <- deep_inter_1[[score]]
    # Loop over thresholds with given step. Get the interception of antismash data with DeepBGC one at given x axis thresholds with additionsl ones
    for (dataframe_1 in seq(input$plot_start, 99, input$plot_step)){
      
      deep_inter <- deep_inter_1 %>%
        dplyr::filter(score>=dataframe_1/100) %>%
        dplyr::select(Start, Stop) 
      if (length(deep_inter$Start) > 0) {
        deep_inter$seqnames <- "chr"
      }
      
      
      # Store antismash bgc start amd atop values as matrix
      if (input$ref_comparison == 'Antismash'){
        anti_inter <- shiny::isolate(vals$anti_data) %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison == 'PRISM'){
        anti_inter <- shiny::isolate(vals$prism_data) %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison == 'SEMPI'){
        anti_inter <- shiny::isolate(vals$sempi_data) %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } 
      
      
      
      # Get the interception of two matrices
      if (length(deep_inter$Start) > 0) {
        query <- GenomicRanges::makeGRangesFromDataFrame(deep_inter)
        subject <- GenomicRanges::makeGRangesFromDataFrame(anti_inter)
        interseption <- GenomicRanges::findOverlaps(query,subject)
        inter_bgc <- length(interseption@from)
        len_new <- length(deep_inter$seqnames) - inter_bgc
      } else {
        inter_bgc <- 0
        len_new <- 0
      }
      
      if (input$ref_comparison == 'Antismash'){
        used_antismash <-  length(shiny::isolate(vals$anti_data$Cluster))-inter_bgc
        cols <-  c("Only Antismash", "DeepBGC+Antismash", "Only DeepBGC")
        title <-  ggplot2::ggtitle("Comparison of Antismash and DeepBGC annotations at given score threshold")
      } else if (input$ref_comparison == 'PRISM'){
        used_antismash <-  length(shiny::isolate(vals$prism_data$Cluster))-inter_bgc
        cols <- c("Only PRISM", "DeepBGC+PRISM", "Only DeepBGC")
        title <- ggplot2::ggtitle("Comparison of PRISM and DeepBGC annotations at given score threshold")
      } else if (input$ref_comparison == 'SEMPI') {
        used_antismash <-  length(shiny::isolate(vals$sempi_data$Cluster))-inter_bgc
        cols <- c("Only SEMPI", "DeepBGC+SEMPI", "Only DeepBGC")
        title <- ggplot2::ggtitle("Comparison of SEMPI and DeepBGC annotations at given score threshold")
      }
      
      # Combine all vectors into one dataframe
      fullnes_of_annotation_1 <- data.frame(c(rep(c(as.character(dataframe_1)),3 )), 
                                            cols, c(used_antismash, inter_bgc, len_new))
      colnames(fullnes_of_annotation_1) <- c("Score", "Source", "Quantity")
      # Combine previously created empty dataframe with this one to store results
      fullnes_of_annotation <- rbind(fullnes_of_annotation, fullnes_of_annotation_1)
      
    }
    
    # Store dataframe in reactive value for later use.
    vals$fullness_deep <- data.frame(fullnes_of_annotation)
    write.csv(fullnes_of_annotation, "fullness.csv", row.names = F)
    
    # Make text to show on a barplot to point on additional scores' thresholds
    annotateText=paste("Applied additional thresholds", paste("Activity score:", as.character(input$score_a)),
                       paste("DeepBGC score:", as.character(input$score_d)),
                       paste("Cluster type score:", as.character(input$score_c)), sep = "\n")
    
    # Plot the barplot
    ggplot2::ggplot(fullnes_of_annotation, ggplot2::aes(fill=Source, y=Quantity, x=Score)) + 
      ggplot2::geom_bar(position="dodge", stat="identity")+
      ggplot2::geom_text(ggplot2::aes(label=Quantity), position=ggplot2::position_dodge(width=0.9), vjust=-0.25) +
      ggplot2::xlab(paste(input$score_type,"Score")) +
      title +
      ggplot2::geom_label(ggplot2::aes(x=Inf,y=Inf,hjust=1,vjust=1,label=annotateText ), show.legend = F)
  })
  
  # Render interactive plot with plotly for rates of DeepBGC data in regards with antismash data
  output$deep_rate <- plotly::renderPlotly({
    shiny::req(!is.null(vals$fullness_deep))
    
    
    # Reuse stored dataframe from previous plot
    # This dataframe stores data for number of intercepted/non intercepted clusters for DeepBGC and antismash data 
    # For more information please see previous shiny::renderPlot
    fullnes_of_annotation <- data.frame(vals$fullness_deep)
    
    # Store dataframe into variable. Widen it to calculate rates
    test <- fullnes_of_annotation %>%
      tidyr::pivot_wider(names_from = Source, values_from = Quantity)
    if (input$ref_comparison == 'Antismash'){
      data <-  vals$anti_data
      title <- ggplot2::ggtitle("Rates of DeepBGC/Antismash data annotation")
      test <- test %>%
        # Calculate rates. Novelty is nummber of clusters annotated only by deepbgc/ all clusters annotated by antismash + (antismash + deepbgc)
        dplyr::mutate(Novelty_rate = test$`Only DeepBGC`/(test$`DeepBGC+Antismash` + test$`Only Antismash`), 
                      #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
                      Annotation_rate = test$`DeepBGC+Antismash`/length(data$Cluster), 
                      # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
                      Skip_rate = test$`Only Antismash`/length(data$Cluster))
    } else if (input$ref_comparison == 'PRISM'){
      data <- vals$prism_data
      title <- ggplot2::ggtitle("Rates of DeepBGC/PRISM data annotation")
      test <- test %>%
        dplyr::mutate(Novelty_rate = test$`Only DeepBGC`/(test$`DeepBGC+PRISM` + test$`Only PRISM`), 
                      #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
                      Annotation_rate = test$`DeepBGC+PRISM`/length(data$Cluster), 
                      # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
                      Skip_rate = test$`Only PRISM`/length(data$Cluster))
    } else if (input$ref_comparison == 'SEMPI'){
      data <- vals$sempi_data
      title <- ggplot2::ggtitle("Rates of DeepBGC/SEMPI data annotation")
      test <- test %>%
        dplyr::mutate(Novelty_rate = test$`Only DeepBGC`/(test$`DeepBGC+SEMPI` + test$`Only SEMPI`), 
                      #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
                      Annotation_rate = test$`DeepBGC+SEMPI`/length(data$Cluster), 
                      # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
                      Skip_rate = test$`Only SEMPI`/length(data$Cluster))
    }
    
    # Calculate rates and plot interactive plot with plotly
    plotly::ggplotly(test %>%
                       tidyr::pivot_longer(cols = c(Novelty_rate, Annotation_rate, Skip_rate), names_to = 'Rates', values_to = 'Rates_data') %>%
                       ggplot2::ggplot(ggplot2::aes(x=as.numeric(Score), y=as.numeric(Rates_data), Rate = as.numeric(Rates_data))) +
                       ggplot2::geom_line(ggplot2::aes(color=Rates)) +
                       ggplot2::geom_point(ggplot2::aes(shape=Rates), alpha = .4, size = 3) +
                       title +
                       ggplot2::ylab("Rate") +
                       ggplot2::xlab(paste(input$score_type,"Score threshold")),
                     tooltip = c("Rate"))
  })
  ##----------------------------------------------------------------
  ##                      GECCO Comparison tab                     -
  ##----------------------------------------------------------------
  # Render barplot
  output$gecco_barplot <- shiny::renderPlot({
    shiny::req((vals$gecco_data_input == T) & ((vals$anti_data_input == T) | (vals$prism_data_input == T) | (vals$sempi_data_input == T) ))
    
    # Create empty dataframe to populate later
    fullnes_of_annotation <- data.frame(NA, NA, NA)
    colnames(fullnes_of_annotation) <- c("Score", "Source", "Quantity")
    fullnes_of_annotation <- tidyr::drop_na(fullnes_of_annotation)
    
    gecco_inter_1 <- vals$gecco_data_filtered
    # Decide which score to use for basic thresholds on x axis
    if (input$score_type_gecco == "avg_p") {
      score <- "score_a"
    } else if (input$score_type_gecco == "Cluster_Type") {
      score <- "score_c"
    } 
    gecco_inter_1$score <- gecco_inter_1[[score]]
    
    # Loop over thresholds with given step. Get the interception of antismash data with DeepBGC one at given x axis thresholds with additionsl ones
    for (dataframe_1 in seq(input$plot_start_gecco, 99, input$plot_step_gecco)){
      
      # dplyr::filter dataframe. Get only rows, which >= of a given thresholds. dplyr::select only start and stop of those rows as a matrix
      gecco_inter <- gecco_inter_1 %>%
        dplyr::filter(score>=dataframe_1/100) %>%
        dplyr::select(Start, Stop) 
      if (length(gecco_inter$Start) > 0) {
        gecco_inter$seqnames <- "chr"
      }
      
      
      # Store antismash bgc start amd atop values as matrix
      if (input$ref_comparison_gecco == 'Antismash'){
        anti_inter <- vals$anti_data %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison_gecco == 'PRISM'){
        anti_inter <- vals$prism_data %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison_gecco == 'SEMPI'){
        anti_inter <- vals$sempi_data %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } 
      
      
      
      
      # Get the interception of two matrices
      if (length(gecco_inter$Start) > 0) {
        query <- GenomicRanges::makeGRangesFromDataFrame(gecco_inter)
        subject <- GenomicRanges::makeGRangesFromDataFrame(anti_inter)
        interseption <- GenomicRanges::findOverlaps(query,subject)
        inter_bgc <- length(interseption@from)
        len_new <- length(gecco_inter$seqnames) - inter_bgc
      } else {
        inter_bgc <- 0
        len_new <- 0
      }
      
      
      if (input$ref_comparison_gecco == 'Antismash'){
        used_antismash <-  length(vals$anti_data$Cluster)-inter_bgc
        cols <-  c("Only Antismash", "GECCO+Antismash", "Only GECCO")
        title <-  ggplot2::ggtitle("Comparison of Antismash and GECCO annotations at given score threshold")
      } else if (input$ref_comparison_gecco == 'PRISM'){
        used_antismash <-  length(vals$prism_data$Cluster)-inter_bgc
        cols <- c("Only PRISM", "GECCO+PRISM", "Only GECCO")
        title <- ggplot2::ggtitle("Comparison of PRISM and GECCO annotations at given score threshold")
      } else if (input$ref_comparison_gecco == 'SEMPI') {
        used_antismash <-  length(vals$sempi_data$Cluster)-inter_bgc
        cols <- c("Only SEMPI", "GECCO+SEMPI", "Only GECCO")
        title <- ggplot2::ggtitle("Comparison of SEMPI and GECCO annotations at given score threshold")
      }
      
      # Combine all vectors into one dataframe
      fullnes_of_annotation_1 <- data.frame(c(rep(c(as.character(dataframe_1)),3 )), 
                                            cols, c(used_antismash, inter_bgc, len_new))
      colnames(fullnes_of_annotation_1) <- c("Score", "Source", "Quantity")
      # Combine previously created empty dataframe with this one to store results
      fullnes_of_annotation <- rbind(fullnes_of_annotation, fullnes_of_annotation_1)
      
    }
    
    # Store dataframe in reactive value for later use.
    vals$fullness_gecco <- data.frame(fullnes_of_annotation)
    write.csv(fullnes_of_annotation, "fullness.csv", row.names = F)
    
    # Make text to show on a barplot to point on additional scores' thresholds
    annotateText=paste("Applied additional thresholds", paste("Average p-value:", as.character(input$score_average_gecco)),
                       paste("Cluster type score:", as.character(input$score_cluster_gecco)), sep = "\n")
    
    # Plot the barplot
    ggplot2::ggplot(fullnes_of_annotation, ggplot2::aes(fill=Source, y=Quantity, x=Score)) + 
      ggplot2::geom_bar(position="dodge", stat="identity")+
      ggplot2::geom_text(ggplot2::aes(label=Quantity), position=ggplot2::position_dodge(width=0.9), vjust=-0.25) +
      ggplot2::xlab(paste(input$score_type,"Score")) +
      title +
      ggplot2::geom_label(ggplot2::aes(x=Inf,y=Inf,hjust=1,vjust=1,label=annotateText ), show.legend = F)
  })
  
  # Render interactive plot with plotly for rates of DeepBGC data in regards with antismash data
  output$gecco_rate <- plotly::renderPlotly({
    shiny::req(!is.null(vals$fullness_gecco))
    
    # Reuse stored dataframe from previous plot
    # This dataframe stores data for number of intercepted/non intercepted clusters for DeepBGC and antismash data 
    # For more information please see previous shiny::renderPlot
    fullnes_of_annotation <- data.frame(vals$fullness_gecco)
    
    # Store dataframe into variable. Widen it to calculate rates
    test <- fullnes_of_annotation %>%
      tidyr::pivot_wider(names_from = Source, values_from = Quantity)
    if (input$ref_comparison_gecco == 'Antismash'){
      data <-  vals$anti_data
      title <- ggplot2::ggtitle("Rates of GECCO/Antismash data annotation")
      test <- test %>%
        # Calculate rates. Novelty is nummber of clusters annotated only by deepbgc/ all clusters annotated by antismash + (antismash + deepbgc)
        dplyr::mutate(Novelty_rate = test$`Only GECCO`/(test$`GECCO+Antismash` + test$`Only Antismash`), 
                      #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
                      Annotation_rate = test$`GECCO+Antismash`/length(data$Cluster), 
                      # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
                      Skip_rate = test$`Only Antismash`/length(data$Cluster))
    } else if (input$ref_comparison_gecco == 'PRISM'){
      data <- vals$prism_data
      title <- ggplot2::ggtitle("Rates of GECCO/PRISM data annotation")
      test <- test %>%
        dplyr::mutate(Novelty_rate = test$`Only GECCO`/(test$`GECCO+PRISM` + test$`Only PRISM`), 
                      #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
                      Annotation_rate = test$`GECCO+PRISM`/length(data$Cluster), 
                      # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
                      Skip_rate = test$`Only PRISM`/length(data$Cluster))
    } else if (input$ref_comparison_gecco == 'SEMPI'){
      data <- vals$sempi_data
      title <- ggplot2::ggtitle("Rates of GECCO/SEMPI data annotation")
      test <- test %>%
        dplyr::mutate(Novelty_rate = test$`Only GECCO`/(test$`GECCO+SEMPI` + test$`Only SEMPI`), 
                      #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
                      Annotation_rate = test$`GECCO+SEMPI`/length(data$Cluster), 
                      # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
                      Skip_rate = test$`Only SEMPI`/length(data$Cluster))
    }
    
    # Calculate rates and plot interactive plot with plotly
    plotly::ggplotly(test %>%
                       tidyr::pivot_longer(cols = c(Novelty_rate, Annotation_rate, Skip_rate), names_to = 'Rates', values_to = 'Rates_data') %>%
                       ggplot2::ggplot(ggplot2::aes(x=as.numeric(Score), y=as.numeric(Rates_data), Rate = as.numeric(Rates_data))) +
                       ggplot2::geom_line(ggplot2::aes(color=Rates)) +
                       ggplot2::geom_point(ggplot2::aes(shape=Rates), alpha = .4, size = 3) +
                       title +
                       ggplot2::ylab("Rate") +
                       ggplot2::xlab(paste(input$score_type,"Score threshold")),
                     tooltip = c("Rate"))
  })
  ##---------------------------------------------------------------
  ##              Annotation on chromosome plots' tab             -
  ##---------------------------------------------------------------
  
  # Render interactive plot, which shows bgcs of antismash, intercepted with chosen app. Also all app bgs. On hover shows all available information
  # For antismash and PRISM data showed only ID, Start, Stop, Type
  output$deep_reference <- plotly::renderPlotly({
    shiny::req(vals$deep_reference_to_plot)
    vals$can_plot_deep_ref_2 <- T
    vals$deep_reference_to_plot
  })
  
  output$deep_reference_2 <- plotly::renderPlotly({
    shiny::req(vals$can_plot_deep_ref_2 == T)
    vals$can_plot_deep_ref_2 == F
    rename_y_axis <- shiny::isolate(vals$rename_y_axis)
    data <- NULL
    
    index <- 1
    for (upload in data_uploads){
      if (is.null(data)){
        if (vals[[upload]] == T){
          if (dim(vals[[data_to_use[index]]])[1] != 0){
            data <- vals[[data_to_use[index]]]
          }
        }
      }
      index <- index+1
    }
    
    
    tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                "P_value", "RRE_start","RRE_stop", "Probability", "Name", "Full_name",  "Hit", "Core", "Count", "Bitscore", "Model",
                "Num_domains", "Num_proteins", "Average_p", "Max_p")
    
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = vals$chr_len, y = Chr))
    if (vals$anti_data_input == TRUE){
      plot <- plot + 
        ggplot2::geom_segment(data=vals$seg_df_ref_a, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,
                                                                   ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    if (vals$deep_data_input == TRUE){
      if (dim(vals$seg_df_ref_d)[1] >0) {
        plot <- plot +
          ggplot2::geom_segment(data=vals$seg_df_ref_d,ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,
                                                                    ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
                                                                    deepbgc_score = deepbgc_score,activity = activity ),size =3)
      }
    }
    if (vals$rre_data_input == TRUE){
      if (vals$rre_more == T){
        plot <- plot + ggplot2::geom_segment(data=vals$seg_df_ref_r, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Score = Score, Software = Software,
                                                                                  ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value,
                                                                                  P_value = P_value, RRE_start = RRE_start,RRE_stop = RRE_stop, 
                                                                                  Probability = Probability),size = 3)
      } else {
        plot <- plot + ggplot2::geom_segment(data=vals$seg_df_ref_r, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2,  Software = Software,
                                                                                  ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value),size = 3)
      }
    }
    if (vals$prism_data_input == TRUE){
      plot <- plot + ggplot2::geom_segment(data=vals$seg_df_ref_p, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,
                                                                                ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
      
      
    }
    if (vals$sempi_data_input == TRUE){
      plot <- plot + ggplot2::geom_segment(data=vals$seg_df_ref_s, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,
                                                                                ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
      
      
    }
    if (input$prism_supp == TRUE){
      plot <- plot + ggplot2::geom_segment(data=vals$seg_df_ref_p_s, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, ID = ID,
                                                                                  Start = Start, Stop = Stop, Type = Type, Name = Name, Full_name = Full_name,
                                                                                  Score = Score), size = 3)
    }
    if (vals$arts_data_input == TRUE){
      plot <- plot + ggplot2::geom_segment(data=vals$seg_df_ref_ar, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, 
                                                                                 ID = ID, Start = Start, Stop = Stop, Type = Type, Hit = Hit,
                                                                                 Core = Core, E_value = E_value, Bitscore = Bitscore, Count = Count,
                                                                                 Model = Model), size = 3)
    }
    if (vals$gecco_data_input == TRUE){
      if (dim(vals$seg_df_ref_g)[1] >0) {
        plot <- plot + ggplot2::geom_segment(data =  vals$seg_df_ref_g, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, 
                                                                                     ID = ID, Start = Start, Stop = Stop, Type = Type, Num_proteins= Num_proteins,
                                                                                     Num_domains = Num_domains,Average_p = Average_p, Max_p = Max_p ), size = 3) 
      }
    }
    to_plot <- plotly::ggplotly(plot +
                                  ggplot2::scale_y_discrete(labels = rename_y_axis) +
                                  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10)) +
                                  ggplot2::ylab("")+
                                  ggplot2::xlab("Chromosome length")+
                                  ggplot2::theme(legend.title = ggplot2::element_blank()) +
                                  ggplot2::ggtitle("All annotations"), 
                                # What actually to visualize in tooltip
                                tooltip = tooltip
    )
    to_plot %>% plotly::layout(legend=list(font = list(
      family = "sans-serif",
      size = 12,
      color = "#000"),
      bordercolor = "#FFFFFF",
      borderwidth = 2,
      title=list(text='<b> Cluster Types </b>')),
      autosize=TRUE)
  }) #%>% shiny::debounce(200)
  
  ##----------------------------------------------------------------
  ##                      Biocircos plot tab                       -
  ##---------------------------------------------------------------
  # Render Biocircos Plot for all-vs-all comparison
  output$biocircos <- BioCircos::renderBioCircos({
    shiny::req(vals$data_upload_count >1)
    
    # Plot BioCircos
    BioCircos::BioCircos(vals$tracklist, genome = vals$Biocircos_chromosomes, genomeTicksScale = 1e+6)
  })
  
  
  output$biocircos_legend <- DT::renderDataTable({
    shiny::req(vals$data_upload_count >=1)
    rownames = FALSE
    new_data <- vals$coloring_datatable
    color_vec <- new_data$x$data$Color
    options(DT.options = list(pageLength = 50))
    new_data %>% DT::formatStyle('Color', backgroundColor=DT::styleEqual(color_vec, color_vec))
    
    
  })
  
  # Updating values in Datatable on edit
  shiny::observeEvent(input$biocircos_legend_cell_edit, {
    if (input$biocircos_legend_cell_edit$col[1] == 0){
      vals$coloring_datatable$x$data$Name <- input$biocircos_legend_cell_edit$value
    } else if (input$biocircos_legend_cell_edit$col[1] == 1){
      vals$coloring_datatable$x$data$Color <- input$biocircos_legend_cell_edit$value
    } else if (input$biocircos_legend_cell_edit$col[1] == 2){
      vals$coloring_datatable$x$data$Hierarchy <- input$biocircos_legend_cell_edit$value
    }
  })
  ##---------------------------------------------------------------
  ##                        Summarize tab                         -
  ##---------------------------------------------------------------
  # Render barplot with number plyr::count of interception for BGC IDs
  output$barplot_rank <- plotly::renderPlotly({
    shiny::req(vals$data_upload_count >1)
    shiny::req(vals$need_filter == F)
    shiny::req(vals$can_plot_barplot_rank == T)
    
    antismash_count <-  NULL
    prism_count <- NULL
    deep_count <- NULL
    rre_count <- NULL
    sempi_count <- NULL
    prism_supp_count <- NULL
    arts_count <- NULL
    gecco_count <- NULL
    
    if (is.null(vals$inters_filtered)){
      inters <- vals$inters
    } else {
      inters <- vals$inters_filtered
    } 
    index <- 1
    ranking_data <- NULL
    for (upload in data_uploads){
      if (vals[[upload]] == T){
        counts_var <-plyr::count(as.factor(unlist(sapply(inters[[soft_names[index]]], function(x){x$to}))))
        # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
        anot_var <- vals[[data_to_use[index]]][vals[[data_to_use[index]]]$Cluster %in% as.numeric(levels(counts_var$x)),]
        # Add prefices to the ID to plot for a barplot.  
        counts_var$x <- sapply(counts_var$x, function(x) paste0(abbr[index],": ", x))
        # Add label column to the dataframe, from which we will plot  
        counts_var$label <- rep(soft_namings[index], length(counts_var$x))
        # Add type to the dataframe, from which we would plot (from annotation dataframe)  
        counts_var$Type <- anot_var$Type
        # Add Start positions (to visualize on hover)
        counts_var$Start <- anot_var$Start
        # Add Stop positions (to visualize on hover)
        counts_var$Stop <- anot_var$Stop
        if (is.null(ranking_data)){
          ranking_data <- counts_var
        } else{
          ranking_data <- rbind(ranking_data, counts_var)
        }
      }
      index <- index +1
    }
    
    
    # Fix column names in the master dataframe
    colnames(ranking_data) <- c("Cluster", "Count", "Label", "Type", "Start", "Stop")
    # Plot
    plotly::ggplotly(ggplot2::ggplot(ranking_data, ggplot2::aes(x = Cluster, y = Count, Type = Type, Start = Start, Stop = Stop)) +
                       ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = Label)) +
                       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 10),
                                      axis.text.y = ggplot2::element_text(size = 14)) +
                       ggplot2::ggtitle("Number of times cluster is annotated with other tool"),
                     tooltip=c("Type", "Start", "Stop")  
    )
    
    
  })
  
  # Render table with data
  output$group_table <- shiny::renderTable({
    shiny::req(vals$data_upload_count >1)
    shiny::req(vals$need_filter == F)
    shiny::req(vals$can_plot_group_table == T)
    source("src/group_table_functions.R")
    if (is.null(vals$inters_filtered)){
      inters <- vals$inters
    } else {
      inters <- vals$inters_filtered
    }
    df_test <- data.frame(matrix(ncol = length(abbr), nrow = 0))
    colnames(df_test) <- abbr
    added_inters <- c(soft_names[match(input$group_by, soft_namings)])
    add_inters <- list()
    if (input$count_all == F){
      df_test[nrow(df_test)+1,] <- NA
    } else{
      selected_dataframe <- data_to_use[match(input$group_by, soft_namings)]
      df_test <- data.frame(matrix(ncol = length(abbr), nrow = length(vals[[selected_dataframe]]$Cluster)))
      colnames(df_test) <- abbr
      df_test[[abbr[match(input$group_by, soft_namings)]]]<- vals[[selected_dataframe]]$Cluster
      df_test[nrow(df_test)+1,] <- NA
    }
    for (i in seq(1:length(data_uploads))){
      if (input$group_by==soft_namings[[i]]){
        exclude <- i
        soft_n <- names(inters[[soft_names[i]]])
        index <- 1
        for (d in seq(1:length(soft_n))) {
          name <- soft_n[[index]]
          df_tmp <-  data.frame(cbind(c(inters[[soft_names[i]]][[soft_n[d]]]$to), c(inters[[soft_names[i]]][[soft_n[d]]]$from)))
          for (h in seq(1:length(soft_n))){
            if (name==soft_names[match(soft_n, soft_names)][h]){
              colnames(df_tmp) <- c(abbr[i],abbr[match(soft_n, soft_names)][h])
              df_test <- merge(df_test, df_tmp,  all = T)
            }
          }
          
          index <- index +1
        }
        excluded_names <- abbr[abbr != as.name(abbr[i])]
        data <- df_test %>% dplyr::group_by_if(colnames(df_test)==abbr[i]) %>% dplyr::summarise(a = paste(eval(as.name(excluded_names[1])), collapse=","),
                                                                                                b=paste(eval(as.name(excluded_names[2])), collapse=","),
                                                                                                c=paste(eval(as.name(excluded_names[3])), collapse=","),
                                                                                                d=paste(eval(as.name(excluded_names[4])), collapse=","),
                                                                                                e=paste(eval(as.name(excluded_names[5])), collapse=","),
                                                                                                f=paste(eval(as.name(excluded_names[6])), collapse=","),
                                                                                                g=paste(eval(as.name(excluded_names[7])), collapse=","))
        colnames(data) <- c(abbr[i], excluded_names)
        for (p in abbr){
          data[[p]] <- gsub('NA,|,NA', '', data[[p]])
          data[[p]][nrow(data)] <- refine_unique(data[[p]])
          names(data)[names(data) == p] <- soft_namings[match(p, abbr)]
        }
        data["Group"] <- paste("group", rownames(data), sep = "_")
        for (f in seq(1:length(data_uploads))){
          if (vals[[data_uploads[f]]] != TRUE){
            data <- data %>%
              dplyr::select(-as.name(soft_namings[f]))
          }
        }
        
      } else {
        if ( !(soft_names[i] %in% added_inters)){
          matched_v <- match(added_inters,names(inters[[soft_names[i]]]))
          soft_n <- soft_names[-(matched_v[!is.na(matched_v)])]
          for (inter in names(inters[[soft_names[i]]])){
            if (!(inter %in% added_inters)){
              add_inters[[soft_names[i]]] <- c(add_inters[[soft_names[i]]],inters[[soft_names[i]]][[inter]]$to )
              add_inters[[inter]] <- c(add_inters[[inter]], inters[[soft_names[i]]][[inter]]$from)
            }
          }
          added_inters <- c(added_inters, soft_names[i])}
      }
    }
    
    for (name in names( add_inters) ){
      data_to_add <- sort(unique(add_inters[[name]]))
      data[nrow(data), soft_namings[match(name, soft_names)]] <-  paste(data_to_add[!(data_to_add %in% unique(unlist(c(data[soft_namings[match(name, soft_names)]]))))], collapse = ",")
    }
    write.csv(data, "group_by.csv", row.names = F)
    data
  })
  
  # Download used datasets (as for BioCircos)
  output$download <- shiny::downloadHandler(filename = function(){
    paste("datasets.zip")     
  },  
  content =  function(file){
    flst <- c()
    # List files in directory
    files_in_dir <- list.files()
    # Iterate over those files and if found "_biocircos.csv" add to the flst vector
    for (file_names in files_in_dir) {
      if (grepl('_biocircos.csv', file_names, fixed = TRUE)) {
        flst <- c(flst, file_names)
      } else if (grepl('group_by.csv', file_names, fixed = TRUE)){
        flst <- c(flst, file_names)
      } else if (grepl('group.py', file_names, fixed = TRUE)){
        flst <- c(flst, file_names)
      }
    }
    #create the zip file from flst vector
    zip(file,  flst) },
  contentType = "application/zip" )
}

# Run the application 
shiny::shinyApp(ui = ui, server = server)
