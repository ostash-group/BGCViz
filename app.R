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
library(rjson)
library(stringr)
library(GenomicRanges)
# Define UI 
ui <- shiny::fluidPage(

  # Application title
  shiny::titlePanel("BGCViz"),
  
  # Sidebar
  shinyjs::useShinyjs(),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      # Define sidebar logic. Make it fixed and not overflow
      id = "tPanel",style = "overflow-y:scroll; max-height: 90vh; position:fixed; width:inherit;",
      # Data upload
      shiny::h3("Data upload and necesary input:"),
      shiny::checkboxInput("hide_uploads", "Hide upload fields"),
      shiny::h5(id = "anti_header_upload","ANTISMASH:"),
      shiny::checkboxInput("anti_input_options", "My AntiSMASH data is a dataframe, not json results file from antismash", value = T),
      shiny::fileInput("anti_data",
                "Upload Antismash data"),
      shiny::actionButton("anti_sco", "Use Antismash example data from S.coelicolor"),
      shiny::h5(id = "prism_header_upload","PRISM:"),
      shiny::checkboxInput("prism_input_options", "My PRISM data is a dataframe, not json results file", value = T),
      shiny::fileInput("prism_data",
                "Upload PRISM data"),
      shiny::actionButton("prism_sco", "Use PRISM example data from S.coelicolor"),
      shiny::h5(id = "sempi_header_upload","SEMPI:"),
      shiny::fileInput("sempi_data",
                "Upload SEMPI 2.0 data"),
      shiny::actionButton("sempi_sco", "Use SEMPI example data from S.coelicolor"),
      shiny::h5(id = "deep_header_upload","DEEPBGC:"),
      shiny::fileInput("deep_data",
                "Upload DeepBGC data"),
      shiny::actionButton("deep_sco", "Use DeepBGC example data from S.coelicolor"),
      shiny::h5(id = "gecco_header_upload","GECCO:"),
      shiny::fileInput("gecco_data",
                "Upload Gecco data"),
      shiny::actionButton("gecco_sco", "Use Gecco example data from S.coelicolor"),
      shiny::h5(id = "rre_header_upload","RRE-FINDER:"),
      shiny::fileInput("rre_data",
                "Upload RRE-Finder data"),
      shiny::actionButton("rre_sco", "Use RRE-Finder example data from S.coelicolor"),
      shiny::h5(id = "arts_header_upload","ARTS:"),
      shiny::fileInput("known_data",
                "Upload ARTS knownhits data"),
      shiny::fileInput("dup_data",
                "Upload ARTS duptable data"),
      shiny::actionButton("arts_sco", "Use ARTS example data from S.coelicolor"),
      # Numeric input of chromosome length of analyzed sequence
      shiny::numericInput("chr_len", "Please type chr len of an organism", value = 10000000),
      shiny::h3("Data manipulation options"),
      shiny::checkboxInput("hide_anti", "Hide data manipulation fields"),
      shiny::h5(id = "anti_header","Antismash data options:"),
      shiny::checkboxInput("anti_hybrid", "Visualize AntiSMASH BGC with several types as 'Hybrid'"),
      shiny::h5(id = "prism_header","PRISM data options:"),
      shiny::checkboxInput("prism_hybrid", "Visualize PRISM BGC with several types as 'Hybrid'"),
      shiny::checkboxInput("prism_supp", "Visualize PRISM resistance and regulatory genes"),
      shiny::h5(id = "sempi_header","SEMPI data options:"),
      shiny::checkboxInput("sempi_hybrid", "Visualize SEMPI BGC with several types as 'Hybrid'"),
      shiny::h5(id = "arts_header","ARTS data options:"),
      shiny::selectInput("dup_choice", "Choose duplicated core gene to plot only it", choices = c("All"),
                  selected = "All"),
      shiny::h3(id = "genes_on_chr","Genes on chromosome plot controls:"),
      shiny::checkboxInput("hide_genes_on_chr", "Hide 'Genes on chromosome plot' fields"),
      shiny::selectInput("ref", "Choose reference data", choices = c(""),
                  selected = ""),
      shiny::h3(id = "summarize","Summarize options:"),
      shiny::checkboxInput("hide_summarize", "Hide summarize options"),
      shiny::selectInput("group_by", "Group data by", choices = c(""),  selected = ''),
      shiny::checkboxInput("count_all", "Show all BGC for the 'group by' method (+ individually annotated BGC)"),
      shiny::h3("Improve visualization:"),
      shiny::checkboxInput("hide_viz", "Hide improve visualization options"),
      shiny::fileInput("rename_data",
               "Upload renaming and coloring scheme"),
      shiny::actionButton("rename", "Rename"),
      shiny::actionButton("reset_name", "Reset"),
      shiny::checkboxInput("rre_width", "Add thickness to RRE results visualization"),
      shiny::checkboxInput("prism_supp_data_input_width", "Add thickness to PRISM resistance + regulatory genes results visualization"),
      shiny::checkboxInput("arts_width", "Add thickness to ARTS results visualization"),
      shiny::checkboxInput("sempi_width", "Add thickness to SEMPI results visualization"),
      shiny::checkboxInput("biocircos_color", "Make arcs in biocircos colorful, based on the class"),
      shiny::checkboxInput("label_color", "Make links in biocircos colorful, based on the class"),
      shiny::selectInput("label_color_class", "Choose the mode to color the links", choices = c("Hierarchical-based" = "H",
                                                                                           "Purity-based" = "P",
                                                                                           "Reference column-based" = "R"
                                                                                         ),
                  selected = 'H'),
      shiny::selectInput("ref_col_biocircos", "Choose reference column to color the links", choices = c(""), selected = ''),
      shiny::h3(id="data_comparison_header_gecco","Comparison with Gecco plots:"),
      shiny::checkboxInput("hide_data_comparison_gecco", "Hide Gecco data comparison options"),
      shiny::selectInput("ref_comparison_gecco", "Choose data for comparison with Gecco", choices = c(""),selected = ''),
      shiny::selectInput("score_type_gecco", "Choose score type to set threshold", choices = c("Average p-value" = "avg_p",
                                                                                  "Cluster_type score" = "Cluster_Type"),
                  selected = "avg_p"),
      shiny::sliderInput("plot_step_gecco", "Choose step for plots(barplot)", min = 1, max = 50,value = 10),
      shiny::sliderInput("plot_start_gecco", "Chose plot start point(barplot)", min = 0, max = 99, value = 0),
      shiny::h3(id="data_filter_header_gecco","Gecco data filtering:"),
      checkboxInput("hide_data_filter_gecco", "Hide Gecco data filtering options"),
      shiny::sliderInput("score_average_gecco", "Average p-value threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50 ),
      shiny::sliderInput("score_cluster_gecco", "Cluster type threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50 ),
      shiny::sliderInput("domains_filter_gecco", "Domain number threshold for Gecco data", min = 0, max = 100, value = 1),
      shiny::sliderInput("prot_filter_gecco", "Protein number threshold for Gecco data", min = 0, max = 100, value = 1),
      shiny::h3(id="data_comparison_header","Comparison with DeepBGC plots:"),
      shiny::checkboxInput("hide_data_comparison", "Hide DeepBGC data comparison options"),
      shiny::selectInput("ref_comparison", "Choose data for comparison with DeepBGC", choices = c(""), selected = ''),
      # Score to use for thresholds
      shiny::selectInput("score_type", "Choose score type to set threshold", choices = c("Activity score" = "Activity",
                                                                                  "Cluster_type score" = "Cluster_Type",
                                                                                  "DeepBGC score" = "DeepBGC"),
                  selected = "Activity score"),
      # Chose step for barplot (as a threshold to draw a bar)
      shiny::sliderInput("plot_step", "Choose step for plots(barplot)", min = 1, max = 50,value = 10),
      shiny::sliderInput("plot_start", "Chose plot start point(barplot)", min = 0, max = 99, value = 0),
      
      # DeepBGC data filtering 
      shiny::h3(id="data_filter_header","DeepBGC data filtering:"),
      shiny::checkboxInput("hide_data_filter", "Hide DeepBGC data filtering options"),
      # Different score filtering. Remain >= of set threshold
      shiny::sliderInput("score_a", "Activity score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      shiny::sliderInput("score_d", "DeepBGC score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      shiny::sliderInput("score_c", "Cluster_type score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      # Domains, biodomains and proteins dplyr::filter. Remain >= of set threshold
      shiny::sliderInput("domains_filter", "Domain number threshold for DeepBGC data", min = 0, max = 100, value = 5),
      shiny::sliderInput("biodomain_filter", "Biodomain number threshold for DeepBGC data", min = 0, max = 100, value = 1),
      shiny::sliderInput("gene_filter", "Protein number threshold for DeepBGC data", min = 0, max = 100, value = 1),
      shiny::sliderInput("cluster_type","Choose threshold to assign cluster type for DeepBGC data ", min = 0, max = 100, value = 50),
     
      # Donwload currently used datasets
      shiny::downloadButton("download","Download currently used datasets (as for Biocircos plot)" )
      
    ),
    
    # Show plots
    shiny::mainPanel(
      # Define tabs and their correcponding plots
      shiny::tabsetPanel(
        shiny::tabPanel(title = "Compare data with DeepBGC", value = 1 ,shiny::plotOutput("deep_barplot",height = "500px"), plotly::plotlyOutput("deep_rate")),
        shiny::tabPanel(title = "Compare data with Gecco", value = 5 ,shiny::plotOutput("gecco_barplot",height = "500px"), plotly::plotlyOutput("gecco_rate")),
        shiny::tabPanel(title = "Annotation visualization and comparison", value = 4,plotly::plotlyOutput("deep_reference_2", height = "500px"), 
                 plotly::plotlyOutput("deep_reference", height = "500px")),
        shiny::tabPanel(title = "Biocircos plot", value = 2, BioCircos::BioCircosOutput("biocircos", height = "1000px"), DT::dataTableOutput("biocircos_legend")),
        shiny::tabPanel(title = "Summarize interception", value = 3,plotly::plotlyOutput("barplot_rank", height = "600px"),shiny::tableOutput("group_table")),
        type = "tabs", id = "main"
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
          input$ref_col_biocircos, vals$inters_filtered, input$prism_supp, input$prism_supp_data_input_width,
          input$arts_width, input$sempi_width, input$rre_width
          )
  })
  inputData <- shiny::reactive({
    list( input$sempi_data,input$rre_data,  input$anti_data, input$prism_data,
          input$known_data,input$dup_data,  input$prism_supp_data, input$deep_data, input$gecco_data,
          input$sempi_sco,input$rre_sco,  input$anti_sco, input$prism_sco,
          input$arts_sco, input$prism_supp_sco, input$deep_sco, input$gecco_sco
    )
  })
  dynamicInput <-  shiny::reactive({
    list(  input$cluster_type, input$gene_filter,input$biodomain_filter,  input$score_c, input$score_d, 
          input$score_a,  input$score_average_gecco,input$score_cluster_gecco, input$domains_filter_gecco, 
          input$prot_filter_gecco, input$dup_choice, vals$need_filter
    )
  })
  


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
  vals$rename_data <- read.csv("rename.csv")
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
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }
  # Fix the duplicates in PRISM-Supp data. Therefore 1 row for 1 orf
  fix_duplicates <-  function(test_score, order_vec, regul_genes_orfs){
    dupl_names <- regul_genes_orfs[duplicated(regul_genes_orfs)]
    duplicated_values <- which(duplicated(regul_genes_orfs[order_vec]))
    test_score <- test_score[order_vec]
    to_add <- test_score[(which(duplicated(regul_genes_orfs[order_vec])))]
    test_score <- test_score[-(which(duplicated(regul_genes_orfs[order_vec])))]
    test_score[duplicated_values-1] <- paste(test_score[duplicated_values-1] , to_add, sep="/")
    return(test_score)
  }
  # PRISM JSON data processing
  process_prism_json_suppl <- function(data){
    types <- sapply(data$prism_results$clusters, function(x){
      tolower(x$type)
    })
    
    types <- sapply(types, function(x){
      if (length(unlist(x))>1){
        tmp <- str_trim(paste0(unlist(x), collapse = '', sep = " "))
        gsub(" ", "__", tmp)
      }else{
        x
      }
    })
    
    start <- sapply(data$prism_results$clusters, function(x){
      x$start
      
    })
    end <- sapply(data$prism_results$clusters, function(x){
      x$end
      
    })
    
    
    prism_data <-  data.frame(Cluster=as.numeric(seq(1:length(start))), Start=as.numeric(start), Stop = as.numeric(end), Type = types)

    regul_genes_orfs <- sapply(data$prism_results$regulatory_genes, function(x){
      x$orf
    })
    
    names <- sapply(data$prism_results$orfs[[1]]$orfs, function(y){
      y$name
    })
    coordinates <- sapply(data$prism_results$orfs[[1]]$orfs, function(y){
      y$coordinates
    })
    
    test_coords <- as.matrix(coordinates[, names %in% regul_genes_orfs ])
    
    
    reg_genes <-data.frame(t(test_coords))
    colnames(reg_genes) <- c("Start", "Stop")
    reg_genes$Type <- 'regulatory'
    reg_genes$Type2 <- reg_genes$Type
    
    test_name <- names[names %in% regul_genes_orfs]
    order_vec <- order(match(regul_genes_orfs, test_name))
    
    
    test_score <- sapply(data$prism_results$regulatory_genes, function(x){
      x$score
    })
    reg_genes$Score  <- fix_duplicates(test_score , order_vec, regul_genes_orfs)
    test_name <- sapply(data$prism_results$regulatory_genes, function(x){
      x$name
    })
    reg_genes$Name  <- fix_duplicates(test_name , order_vec, regul_genes_orfs)
    test_full_name<- sapply(data$prism_results$regulatory_genes, function(x){
      x$full_name
    })
    reg_genes$Full_name  <- fix_duplicates(test_full_name , order_vec, regul_genes_orfs)
    resist_genes_orfs <- sapply(data$prism_results$resistance_genes, function(x){
      x$orf
    })
    
    test_coords_res <- as.matrix(coordinates[, names %in% resist_genes_orfs ])
    
    
    res_genes <-data.frame(t(test_coords_res))
    
    colnames(res_genes) <- c("Start", "Stop")
    res_genes$Type <- 'resistance'
    res_genes$Type2 <- res_genes$Type
    test_name <- names[names %in% resist_genes_orfs]
    order_vec <- order(match(resist_genes_orfs, test_name))
    
    
    test_score <- sapply(data$prism_results$resistance_genes, function(x){
      x$score
    })
    res_genes$Score  <- fix_duplicates(test_score , order_vec, resist_genes_orfs)
    test_name <- sapply(data$prism_results$resistance_genes, function(x){
      x$name
    })
    res_genes$Name  <- fix_duplicates(test_name , order_vec, resist_genes_orfs)
    test_full_name<- sapply(data$prism_results$resistance_genes, function(x){
      x$full_name
    })
    res_genes$Full_name  <- fix_duplicates(test_full_name , order_vec, resist_genes_orfs)
    
    final_reg <- rbind(res_genes, reg_genes)
    final_reg$ID <- seq(1:dim(final_reg)[1])
    final_reg$Cluster <- final_reg$ID
    rownames(final_reg) <- as.numeric(seq(1:dim(final_reg)[1]))
    return(list(prism_data, final_reg))
  }
  # Filtering the DeepBGC
  filter_deepbgc <- function(){
    score_a <- apply(vals$deep_data %>% dplyr::select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
    score_d <- apply(vals$deep_data %>% dplyr::select(c("deepbgc_score")),1, function(x) max(x))
    score_c <- apply(vals$deep_data %>% dplyr::select(c("alkaloid", "nrps","other","pks","ripp","saccharide","terpene")),1, function(x) max(x))
    deep_data_chromo <- vals$deep_data %>%
      dplyr::mutate(score = apply(vals$deep_data %>%
                             dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, function(x) max(x))) 
    # Cluster_type column. Here extract colnames, and assign max value to a new column
    deep_data_chromo$Cluster_type <- colnames(deep_data_chromo %>% dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene))[apply(deep_data_chromo%>%dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, which.max) ]
    # If max score is under_threshold, print "under_threshold"
    deep_data_chromo <- deep_data_chromo%>%
      dplyr::mutate(Cluster_type = ifelse(score>as.numeric(input$cluster_type)/100, Cluster_type, "under_threshold"))
    #Finally store deepbgc data in plotting variable. Do final scores processing 
    biocircos_deep <- deep_data_chromo%>%
      dplyr::mutate( product_class = Cluster_type, score_a = score_a, score_d = score_d, score_c = score_c) %>%
      dplyr::filter(score_a >= as.numeric(input$score_a )/ 100, score_c >=as.numeric(input$score_c)/100 , 
             score_d >= as.numeric(input$score_d)/100,  num_domains >= input$domains_filter,
             num_bio_domains>=input$biodomain_filter, num_proteins>=input$gene_filter)
    biocircos_deep['Start'] <- biocircos_deep$nucl_start
    biocircos_deep['Stop'] <- biocircos_deep$nucl_end
    biocircos_deep['Type'] <- biocircos_deep$product_class
    biocircos_deep['Type2'] <- biocircos_deep$product_class
    biocircos_deep['Cluster'] <- biocircos_deep$ID
    return(biocircos_deep)
  }
  # Filtering GECCO
  filter_gecco <- function(){
    score_a_gecco <- apply(vals$gecco_data %>% dplyr::select(c("average_p")),1, function(x) max(x))
    score_c_gecco <- apply(vals$gecco_data %>% dplyr::select(c("alkaloid", "nrps","other","pks","ripp","saccharide","terpene")),1, function(x) max(x))
    # Store master prism data in local variable
    gecco_data <- vals$gecco_data %>%
      dplyr::mutate(score = apply(vals$gecco_data %>%
                             dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, function(x) max(x))) %>%
      dplyr::mutate(Cluster_type = ifelse(score>as.numeric(input$score_cluster_gecco)/100, Type2, "under_threshold")) %>%
      dplyr::mutate( Type2 = Cluster_type, score_a = score_a_gecco, score_c = score_c_gecco) %>%
      dplyr::filter(score_a >= as.numeric(input$score_average_gecco )/ 100, score_c >=as.numeric(input$score_cluster_gecco)/100 ,
             num_domains >= input$domains_filter_gecco, num_prot>=input$prot_filter_gecco)
    return(gecco_data)
  }
  # Renaming the vector for inut$rename event
  rename_vector <- function(data, renamed_dataframe){
    type <- str_split(data$Type2, "__")
    type_2 <- sapply(type, function(x){
      sapply(x, function(y){
        if (y %in% renamed_dataframe$Code){
          renamed <- as.character(renamed_dataframe$Group[renamed_dataframe$Code == y])
          if( (length(renamed) >1) & (!( as.character(y) %in% names(vals$renaming_notification)))){
            shiny::showNotification(paste("The ", as.character(y), " type have multiple renaming options: ", paste(renamed, collapse = ", ")), 
                             type = "warning", duration = NULL)
            shiny::showNotification(paste("The  ", renamed[[1]], " was chosen."), type = "warning", duration=NULL)
            vals$renaming_notification[[as.character(y)]] <- renamed[[1]] 
          }
          renamed[[1]]
        } else {
          y
        }
        
      })
      
    })
    type_3 <- sapply(type_2, function(x){
      dupl <-  x[!duplicated(x)]
      paste(dupl, collapse = "__")
    })
    type_4 <- sapply(type_3, function(y){
      if (y %in% as.character(renamed_dataframe$Code)){
        as.character(renamed_dataframe$Group[renamed_dataframe$Code == y])[[1]]
      } else {
        y
      }
    })
    return(as.character(type_4))
  }
  # Adding the thickness to the visualization for SEMPI, 
  # ARTS, RRE, PRISN-Supp data
  correct_width <- function(data, label){
    if ((label == 'SEMPI')&(input$sempi_width == T)){
      data$Stop <- data$Stop + 30000
    } else if ((label == 'PRISM-Supp')&(input$prism_supp_data_input_width == T)){
      data$Stop <- data$Stop + 20000
    } else if ((label == 'ARTS')&(input$arts_width == T)){
      data$Stop <- data$Stop + 30000
    } else if ((label == 'RRE-Finder')&(input$rre_width == T)){
      data$Stop <- data$Stop + 50000
    }
    return(data)
  }
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
  #----------------------------------------------------------------
  ##            Loading and processing of example data             -
  ##----------------------------------------------------------------
  shiny::observeEvent(input$anti_sco,{
    
    anti_data <- read.csv("example_data/sco_antismash.csv")
    # Add chromosome column
    anti_data$chromosome <-  rep("A", length(anti_data$Cluster))
    # Type magic
    anti_data$Type <- str_trim(tolower(anti_data$Type))
    anti_data['Type2'] <- str_trim(tolower(anti_data$Type))
    vals$anti_type <- anti_data$Type2
    vals$anti_data <- anti_data
    # Save file
    write.csv(vals$anti_data, "anti_data.csv", row.names = F)
    vals$anti_data_input = TRUE 
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "Antismash" = "Antismash")
    vals$choices$group_by <- c(vals$choices$group_by, "Antismash" = "A")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "Antismash" = "Antismash")
    vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "Antismash" = "A")
    vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "Antismash" = "A")
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
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "Antismash" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "A" )
      shiny::updateSelectInput(session, "ref_comparison",
                        selected = "A")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        selected =  "Antismash")
      shiny::updateSelectInput(session, "ref_comparison_gecco",
                        selected = "A")

    }
    
  })
  
  shiny::observeEvent(input$gecco_sco,{
    
    gecco_data <- read.delim("example_data/sco_gecco.tsv")
    # Add chromosome column
    gecco_data$chromosome <-  rep("G", length(gecco_data$type))
    # Type magic
    gecco_data$Cluster <- seq(1:length(gecco_data$chromosome))
    gecco_data$ID <- gecco_data$Cluster
    gecco_data$Type <- str_trim(tolower(gecco_data$type))
    gecco_data$Type <- gsub("polyketide", "pks", gecco_data$Type)
    gecco_data$Type <- gsub("nrp", "nrps", gecco_data$Type)
    gecco_data$Type <- gsub("unknown", "under_threshold", gecco_data$Type)
    gecco_data['Type2'] <- str_trim(tolower(gecco_data$Type))
    drop_cols <- c("alkaloid_probability" ,  "polyketide_probability", "ripp_probability",  "saccharide_probability",
                   "terpene_probability",    "nrp_probability"  , "other_probability" )
    # Read data
    gecco_data <- gecco_data %>%
      dplyr::mutate(pks=polyketide_probability, other = other_probability, nrps = nrp_probability, alkaloid = alkaloid_probability, 
             terpene = terpene_probability, saccharide = saccharide_probability, ripp = ripp_probability) %>%
      dplyr::select(-dplyr::one_of(drop_cols))
    gecco_data$num_prot <- sapply( str_split(as.character(gecco_data$proteins), ";"), length)
    gecco_data$num_domains <- sapply( str_split(as.character(gecco_data$domains), ";"), length)
    names(gecco_data)[names(gecco_data) == "start"] <- "Start"
    names(gecco_data)[names(gecco_data) == "end"] <-  "Stop"
    vals$gecco_data <- gecco_data
    vals$gecco_data_filtered <- filter_gecco()
    # Save file
    write.csv(vals$gecco_data, "gecco_data.csv", row.names = F)
    vals$gecco_data_input = TRUE 
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "GECCO" = "GECCO")
    vals$choices$group_by <- c(vals$choices$group_by, "GECCO" = "G")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "GECCO" = "GECCO")
    shiny::updateSelectInput(session, "ref",
                      choices = vals$choices$ref )
    shiny::updateSelectInput(session, "group_by",
                      choices = vals$choices$group_by )
    shiny::updateSelectInput(session, "ref_col_biocircos",
                      choices = vals$choices$ref_col_biocircos )
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "GECCO" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "G")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        selected =  "GECCO")
      
    }
    
  })
  
  shiny::observeEvent(input$prism_sco,{
    # Read data
    
      data <- fromJSON(file = "example_data/sco_prism.json")
      processed_data <- process_prism_json_suppl(data)
      shiny::updateCheckboxInput(inputId = "prism_supp", value = T)
      prism_data <- processed_data[[1]]
      vals$prism_supp_data_input = T
      vals$prism_supp <- processed_data[[2]]
      vals$prism_supp_data <- processed_data[[2]]
      vals$prism_json = T
      
      vals$choices$ref <- c(vals$choices$ref, "PRISM-supp" = "PRISM-supp")
      vals$choices$group_by <- c(vals$choices$group_by, "PRISM-supp" = "PS")
      vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "PRISM-Supp" = "PRISM-Supp")
      shiny::updateSelectInput(session, "ref",
                        choices = vals$choices$ref )
      shiny::updateSelectInput(session, "group_by",
                        choices = vals$choices$group_by )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        choices = vals$choices$ref_col_biocircos )
    prism_data$Type <- str_trim(tolower(prism_data$Type))
    prism_data['Type2'] <- str_trim(tolower(prism_data$Type))
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
    vals$choices$group_by <- c(vals$choices$group_by, "PRISM" = "P")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "PRISM" = "PRISM")
    vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "PRISM" = "P")
    vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "PRISM" = "P")
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
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "PRISM" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "P" )
      shiny::updateSelectInput(session, "ref_comparison",
                        selected = "P")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        selected =  "PRISM")
      shiny::updateSelectInput(session, "ref_comparison_gecco",
                        selected = "P")
    }
  })
  
  shiny::observeEvent(input$sempi_sco,{
    
    sempi_data <- read.csv("example_data/sco_sempi.csv")
    sempi_data['Type2'] <- str_trim(tolower(sempi_data$Type))
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
    vals$choices$group_by <- c(vals$choices$group_by, "SEMPI" = "S")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "SEMPI" = "SEMPI")
    vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "SEMPI" = "S")
    vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "SEMPI" = "S")
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
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "SEMPI" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "S" )
      shiny::updateSelectInput(session, "ref_comparison",
                        selected = "S")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        selected =  "SEMPI")
      shiny::updateSelectInput(session, "ref_comparison_gecco",
                        selected = "S")
    }
    
  })
  
  shiny::observeEvent(input$arts_sco, {
    
    data <- read.delim("example_data/sco_duptable.tsv")
    disable_event_logic()
    get_location_duptable <- function(x, y){
      test <- str_split(x, ";")
      test2<- sub(".*loc\\|", "", test[[1]])
      test3 <- str_split(test2, " ")
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
    
    data <- read.delim("example_data/sco_knownhits.tsv")
    locations <- sapply(data$Sequence.description, function(x){
      tail(str_split(x , "\\|")[[1]], 1)
    })
    
    start <- sapply(locations, function(x){
      str_split(x, "_")[[1]][1]
    })
    stop <- sapply(locations, function(x){
      str_split(x, "_")[[1]][2]
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
      vals$upl_arts = T
      vals$choices$ref <- c(vals$choices$ref, "ARTS" = "ARTS")
      vals$choices$group_by <- c(vals$choices$group_by, "ARTS" = "AR")
      vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "ARTS" = "ARTS")
      shiny::updateSelectInput(session, "ref",
                        choices = vals$choices$ref )
      shiny::updateSelectInput(session, "group_by",
                        choices = vals$choices$group_by )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        choices = vals$choices$ref_col_biocircos )
      if (vals$data_upload_count == 1){
        shiny::updateSelectInput(session, "ref",
                          selected = "ARTS" )
        shiny::updateSelectInput(session, "group_by",
                          selected = "AR" )
        shiny::updateSelectInput(session, "ref_col_biocircos",
                          selected =  "ARTS")
      }
  })
  
  shiny::observeEvent(input$deep_sco, {
    
    drop_cols <- c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")
    # Read data
    vals$deep_data <- read.delim("example_data/sco_deep.tsv") %>%
      dplyr::mutate(pks=Polyketide, other = Other, nrps = NRP, alkaloid = Alkaloid, 
             terpene = Terpene, saccharide = Saccharide, ripp = RiPP) %>%
      dplyr::select(-dplyr::one_of(drop_cols))
    # Add chromosome info column
    vals$deep_data$chromosome <-  rep("D", length(vals$deep_data$bgc_candidate_id))
    # Add ID column as number seuquence of dataframe length
    vals$deep_data$ID <- seq(1:length(vals$deep_data$bgc_candidate_id))
    vals$deep_data$Cluster <- vals$deep_data$ID
    write.csv(vals$deep_data, "deep_data.csv", row.names = F)
    vals$deep_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
    vals$deep_data_filtered <- filter_deepbgc()
    vals$choices$ref <- c(vals$choices$ref, "DeepBGC" = "DeepBGC")
    vals$choices$group_by <- c(vals$choices$group_by, "DeepBGC" = "D")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "DeepBGC" = "DeepBGC")
    shiny::updateSelectInput(session, "ref",
                      choices = vals$choices$ref )
    shiny::updateSelectInput(session, "group_by",
                      choices = vals$choices$group_by )
    shiny::updateSelectInput(session, "ref_col_biocircos",
                      choices = vals$choices$ref_col_biocircos )
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "DeepBGC" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "D" )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        choices = "DeepBGC",
                        selected = "DeepBGC")
      
    }
  })
  
  shiny::observeEvent(input$rre_sco, {
    
    # Read data
    vals$rre_data <- read.delim("example_data/sco_rre.txt")
    # Clean RRE data. Extract coordinates and Locus tag with double underscore delimiter (__)
    vals$rre_data <- vals$rre_data %>%
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
    write.csv(vals$rre_data, "rre_data.csv", row.names = F)
    
    vals$rre_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "RRE-Finder" = "RRE-Finder")
    vals$choices$group_by <- c(vals$choices$group_by, "RRE-Finder" = "R")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "RRE-Finder" = "RRE")
    shiny::updateSelectInput(session, "ref",
                      choices = vals$choices$ref )
    shiny::updateSelectInput(session, "group_by",
                      choices = vals$choices$group_by )
    shiny::updateSelectInput(session, "ref_col_biocircos",
                      choices = vals$choices$ref_col_biocircos )
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "RRE-Finder" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "R" )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        choices = "RRE-Finder",
                        selected = "RRE")
      
    }
    if (!is.null(vals$rre_data$Probability)){
      vals$rre_more = T
    } else {
      vals$rre_more = F
    }
  })
  
  ##----------------------------------------------------------------
  ##                Loading and processing user data               -
  ##----------------------------------------------------------------
  shiny::observeEvent(input$anti_data,{
    
    disable_event_logic()
    # Read data
    if (input$anti_input_options==T){
      anti_data <- read.csv(input$anti_data$datapath)
    }else{
       data <- fromJSON(file = input$anti_data$datapath)
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
            tmp <- str_trim(paste0(unlist(x), collapse = '', sep = " "))
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

    # Add chromosome column
    anti_data$chromosome <-  rep("A", length(anti_data$Cluster))
    # Type magic
    anti_data$Type <- str_trim(tolower(anti_data$Type))
    anti_data['Type2'] <- str_trim(tolower(anti_data$Type))
    vals$anti_type <- anti_data$Type2
    vals$anti_data <- anti_data
    # Save file
    write.csv(vals$anti_data, "anti_data.csv", row.names = F)
    vals$anti_data_input = TRUE 
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "Antismash" = "Antismash")
    vals$choices$group_by <- c(vals$choices$group_by, "Antismash" = "A")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "Antismash" = "Antismash")
    vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "Antismash" = "A")
    vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "Antismash" = "A")
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
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "Antismash" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "A" )
      shiny::updateSelectInput(session, "ref_comparison",
                        selected = "A")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        selected =  "Antismash")
      shiny::updateSelectInput(session, "ref_comparison_gecco",
                        selected = "A")
    }

  })
  
  shiny::observeEvent(input$sempi_data,{
    
    
    sempi_data <- read.csv(input$sempi_data$datapath)
    sempi_data['Type2'] <- str_trim(tolower(sempi_data$Type))
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
    vals$choices$group_by <- c(vals$choices$group_by, "SEMPI" = "S")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "SEMPI" = "SEMPI")
    vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "SEMPI" = "S")
    vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "SEMPI" = "S")
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
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "SEMPI" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "S" )
      shiny::updateSelectInput(session, "ref_comparison",
                        selected = "S")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        selected =  "SEMPI")
      shiny::updateSelectInput(session, "ref_comparison_gecco",
                        selected = "S")
    }
    
  })
  
  shiny::observeEvent(input$gecco_data,{
    
    gecco_data <- read.delim(input$gecco_data$datapath)
    gecco_data$chromosome <-  rep("G", length(gecco_data$type))
    # Type magic
    gecco_data$Cluster <- seq(1:length(gecco_data$chromosome))
    gecco_data$ID <- gecco_data$Cluster
    gecco_data$Type <- str_trim(tolower(gecco_data$type))
    gecco_data$Type <- gsub("polyketide", "pks", gecco_data$Type)
    gecco_data$Type <- gsub("nrp", "nrps", gecco_data$Type)
    gecco_data$Type <- gsub("unknown", "under_threshold", gecco_data$Type)
    gecco_data['Type2'] <- str_trim(tolower(gecco_data$Type))
    drop_cols <- c("alkaloid_probability" ,  "polyketide_probability", "ripp_probability",  "saccharide_probability",
                   "terpene_probability",    "nrp_probability"  , "other_probability" )
    # Read data
    gecco_data <- gecco_data %>%
      dplyr::mutate(pks=polyketide_probability, other = other_probability, nrps = nrp_probability, alkaloid = alkaloid_probability, 
             terpene = terpene_probability, saccharide = saccharide_probability, ripp = ripp_probability) %>%
      dplyr::select(-dplyr::one_of(drop_cols))
    gecco_data$num_prot <- sapply( str_split(as.character(gecco_data$proteins), ";"), length)
    gecco_data$num_domains <- sapply( str_split(as.character(gecco_data$domains), ";"), length)
    names(gecco_data)[names(gecco_data) == "start"] <- "Start"
    names(gecco_data)[names(gecco_data) == "end"] <-  "Stop"
    vals$gecco_data <- gecco_data
    vals$gecco_data_filtered <- filter_gecco()
    # Save file
    write.csv(vals$gecco_data, "gecco_data.csv", row.names = F)
    vals$gecco_data_input = TRUE 
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "GECCO" = "GECCO")
    vals$choices$group_by <- c(vals$choices$group_by, "GECCO" = "G")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "GECCO" = "GECCO")
    shiny::updateSelectInput(session, "ref",
                      choices = vals$choices$ref )
    shiny::updateSelectInput(session, "group_by",
                      choices = vals$choices$group_by )
    shiny::updateSelectInput(session, "ref_col_biocircos",
                      choices = vals$choices$ref_col_biocircos )
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "GECCO" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "G")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        selected =  "GECCO")
    }
    
  })
  
  # These are for ARTS data processing
  # input$known_data and inoput$dup_data
  shiny::observeEvent(input$known_data, {
    disable_event_logic()
    
    data <- read.delim(input$known_data$datapath)
    locations <- sapply(data$Sequence.description, function(x){
      tail(str_split(x , "\\|")[[1]], 1)
    })
    
    start <- sapply(locations, function(x){
      str_split(x, "_")[[1]][1]
    })
    stop <- sapply(locations, function(x){
      str_split(x, "_")[[1]][2]
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
      vals$choices$group_by <- c(vals$choices$group_by, "ARTS" = "AR")
      vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "ARTS" = "ARTS")
      shiny::updateSelectInput(session, "ref",
                        choices = vals$choices$ref )
      shiny::updateSelectInput(session, "group_by",
                        choices = vals$choices$group_by )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        choices = vals$choices$ref_col_biocircos )
      dup_table <- vals$dup_data
      known_table <- vals$known_data
      arts_data <- rbind(dup_table, known_table)
      arts_data <- arts_data %>%
        dplyr::arrange(Start)
      arts_data$ID <- seq(1:dim(arts_data)[1])
      arts_data$Cluster <- arts_data$ID
      vals$arts_data <- arts_data
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
                          selected = "AR" )
        shiny::updateSelectInput(session, "ref_col_biocircos",
                          selected =  "ARTS")
      }
    } 
  })
  
  shiny::observeEvent(input$dup_data, {
    disable_event_logic()
    
    data <- read.delim(input$dup_data$datapath)
    
    get_location_duptable <- function(x, y){
      test <- str_split(x, ";")
      test2<- sub(".*loc\\|", "", test[[1]])
      test3 <- str_split(test2, " ")
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
      vals$choices$group_by <- c(vals$choices$group_by, "ARTS" = "AR")
      vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "ARTS" = "ARTS")
      shiny::updateSelectInput(session, "ref",
                        choices = vals$choices$ref )
      shiny::updateSelectInput(session, "group_by",
                        choices = vals$choices$group_by )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        choices = vals$choices$ref_col_biocircos )
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
                          selected = "AR" )
        shiny::updateSelectInput(session, "ref_col_biocircos",
                          selected =  "ARTS")
      }
    } 
  })
  
  shiny::observeEvent(input$prism_data,{
    
    # Read data
    if (input$prism_input_options == T){
      prism_data <- read.csv(input$prism_data$datapath)
    } else{
      data <- fromJSON(file = input$prism_data$datapath)
      processed_data <- process_prism_json_suppl(data)
      shiny::updateCheckboxInput(inputId = "prism_supp", value = T)
      prism_data <- processed_data[[1]]
      vals$prism_supp <- processed_data[[2]]
      vals$prism_supp_data_input = T
      vals$prism_supp_data <- processed_data[[2]]
      vals$prism_json = T
      vals$choices$ref <- c(vals$choices$ref, "PRISM-supp" = "PRISM-supp")
      vals$choices$group_by <- c(vals$choices$group_by, "PRISM-supp" = "PS")
      vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "PRISM-Supp" = "PRISM-Supp")
      shiny::updateSelectInput(session, "ref",
                        choices = vals$choices$ref )
      shiny::updateSelectInput(session, "group_by",
                        choices = vals$choices$group_by )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        choices = vals$choices$ref_col_biocircos )
    }
    vals$choices$ref <- c(vals$choices$ref, "PRISM" = "PRISM")
    vals$choices$group_by <- c(vals$choices$group_by, "PRISM" = "P")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "PRISM" = "PRISM")
    vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "PRISM" = "P")
    vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "PRISM" = "P")
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
    prism_data$Type <- str_trim(tolower(prism_data$Type))
    prism_data['Type2'] <- str_trim(tolower(prism_data$Type))
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
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "PRISM" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "P" )
      shiny::updateSelectInput(session, "ref_comparison",
                        selected = "P")
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        selected =  "PRISM")
      shiny::updateSelectInput(session, "ref_comparison_gecco",
                        selected = "P")
    }
  })
  
  shiny::observeEvent(input$deep_data, {
    
    drop_cols <- c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")
    # Read data
    vals$deep_data <- read.delim(input$deep_data$datapath) %>%
      dplyr::mutate(pks=Polyketide, other = Other, nrps = NRP, alkaloid = Alkaloid, 
             terpene = Terpene, saccharide = Saccharide, ripp = RiPP) %>%
      dplyr::select(-dplyr::one_of(drop_cols))
    # Add chromosome info column
    vals$deep_data$chromosome <-  rep("D", length(vals$deep_data$bgc_candidate_id))
    # Add ID column as number seuquence of dataframe length
    vals$deep_data$ID <- seq(1:length(vals$deep_data$bgc_candidate_id))
    vals$deep_data$Cluster <- vals$deep_data$ID
    write.csv(vals$deep_data, "deep_data.csv", row.names = F)
    vals$deep_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
    vals$deep_data_filtered <- filter_deepbgc()
    vals$choices$ref <- c(vals$choices$ref, "DeepBGC" = "DeepBGC")
    vals$choices$group_by <- c(vals$choices$group_by, "DeepBGC" = "D")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "DeepBGC" = "DeepBGC")
    shiny::updateSelectInput(session, "ref",
                      choices = vals$choices$ref )
    shiny::updateSelectInput(session, "group_by",
                      choices = vals$choices$group_by )
    shiny::updateSelectInput(session, "ref_col_biocircos",
                      choices = vals$choices$ref_col_biocircos )
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "DeepBGC" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "D" )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        choices = "DeepBGC",
                        selected = "DeepBGC")
      
    }
  })
  
  shiny::observeEvent(input$rre_data, {
    
    # Read data
    vals$rre_data <- read.delim(input$rre_data$datapath)
    # Clean RRE data. Extract coordinates and Locus tag with double underscore delimiter (__)
    vals$rre_data <- vals$rre_data %>%
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
    write.csv(vals$rre_data, "rre_data.csv", row.names = F)
    vals$rre_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
    vals$choices$ref <- c(vals$choices$ref, "RRE-Finder" = "RRE-Finder")
    vals$choices$group_by <- c(vals$choices$group_by, "RRE-Finder" = "R")
    vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "RRE-Finder" = "RRE")
    shiny::updateSelectInput(session, "ref",
                      choices = vals$choices$ref )
    shiny::updateSelectInput(session, "group_by",
                      choices = vals$choices$group_by )
    shiny::updateSelectInput(session, "ref_col_biocircos",
                      choices = vals$choices$ref_col_biocircos )
    disable_event_logic()
    if (vals$data_upload_count == 1){
      shiny::updateSelectInput(session, "ref",
                        selected = "RRE-Finder" )
      shiny::updateSelectInput(session, "group_by",
                        selected = "R" )
      shiny::updateSelectInput(session, "ref_col_biocircos",
                        choices = "RRE-Finder",
                        selected = "RRE")
      
    }
    if (!is.null(vals$rre_data$Probability)){
      vals$rre_more = T
    } else {
      vals$rre_more = F
    }
  })
  
  ############################################################################
  ############################################################################
  ###                                                                      ###
  ###                INTERFACE LOGIC: WHAT TO SHOW AND WHEN                ###
  ###                                                                      ###
  ############################################################################
  ############################################################################
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
      if (input$hide_viz == F){
        shinyjs::showElement(selector = "#rre_width")
      }
    } else{
      shinyjs::hideElement(selector = "#rre_width")
    }
  })
  # Show anti_hybrid option if data is available
  # And checkbox is unchecked
  shiny::observeEvent(vals$anti_data_input, {
    
    if (vals$anti_data_input == T){
      if (input$hide_anti == F){
        shinyjs::showElement(selector = "#anti_header")
        shinyjs::showElement(selector = "#anti_hybrid")
      }
    } else{
      shinyjs::hideElement(selector = "#anti_header")
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
      if (input$hide_anti == F){
        shinyjs::showElement(selector = "#prism_header")
        shinyjs::showElement(selector = "#prism_hybrid")
        if (vals$prism_json == T){
          shinyjs::showElement(selector = "#prism_supp")
        }
      }
      if (input$hide_viz == F){
        if (vals$prism_json == T){
          shinyjs::showElement(selector = "#prism_supp_data_input_width")
        }
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
      if (input$hide_anti == F){
        shinyjs::showElement(selector = "#sempi_header")
        shinyjs::showElement(selector = "#sempi_hybrid")
      }
      if (input$hide_viz == F){
        shinyjs::showElement(selector = "#sempi_width")
      }
    } else{
      shinyjs::hideElement(selector = "#sempi_header")
      shinyjs::hideElement(selector = "#sempi_hybrid")
      shinyjs::hideElement(selector = "#sempi_width")
    }
  })
  # Show DeepBGC options if data is available
  shiny::observeEvent(vals$deep_data_input,{
    
    if (vals$deep_data_input == T){
      shinyjs::showElement(selector = "#ref_comparison")
      shinyjs::showElement(selector = "#hide_data_comparison")
      shinyjs::showElement(selector = "#hide_data_filter")
      shinyjs::showElement(selector = "#score_type")
      shinyjs::showElement(selector = "#plot_step")
      shinyjs::showElement(selector = "#plot_start")
      shinyjs::showElement(selector = "#score_a")
      shinyjs::showElement(selector = "#score_d")
      shinyjs::showElement(selector = "#score_c")
      shinyjs::showElement(selector = "#domains_filter")
      shinyjs::showElement(selector = "#biodomain_filter")
      shinyjs::showElement(selector = "#gene_filter")
      shinyjs::showElement(selector = "#cluster_type")
      shinyjs::showElement(selector = "#data_comparison_header")
      shinyjs::showElement(selector = "#data_filter_header")
    } else{
      shinyjs::hideElement(selector = "#ref_comparison")
      shinyjs::hideElement(selector = "#score_type")
      shinyjs::hideElement(selector = "#hide_data_comparison")
      shinyjs::hideElement(selector = "#hide_data_filter")
      shinyjs::hideElement(selector = "#plot_step")
      shinyjs::hideElement(selector = "#plot_start")
      shinyjs::hideElement(selector = "#score_a")
      shinyjs::hideElement(selector = "#score_d")
      shinyjs::hideElement(selector = "#score_c")
      shinyjs::hideElement(selector = "#domains_filter")
      shinyjs::hideElement(selector = "#biodomain_filter")
      shinyjs::hideElement(selector = "#gene_filter")
      shinyjs::hideElement(selector = "#cluster_type")
      shinyjs::hideElement(selector = "#data_comparison_header")
      shinyjs::hideElement(selector = "#data_filter_header")
    }
  })
  # Show GECCO data options, if data is uploaded
  shiny::observeEvent(vals$gecco_data_input,{
    
    if (vals$gecco_data_input == T){
      shinyjs::showElement(selector = "#data_comparison_header_gecco")
      shinyjs::showElement(selector = "#hide_data_comparison_gecco")
      shinyjs::showElement(selector = "#ref_comparison_gecco")
      shinyjs::showElement(selector = "#score_type_gecco")
      shinyjs::showElement(selector = "#plot_step_gecco")
      shinyjs::showElement(selector = "#plot_start_gecco")
      shinyjs::showElement(selector = "#data_filter_header_gecco")
      shinyjs::showElement(selector = "#hide_data_filter_gecco")
      shinyjs::showElement(selector = "#score_average_gecco")
      shinyjs::showElement(selector = "#score_cluster_gecco")
      shinyjs::showElement(selector = "#domains_filter_gecco")
      shinyjs::showElement(selector = "#prot_filter_gecco")
    } else{
      shinyjs::hideElement(selector = "#data_comparison_header_gecco")
      shinyjs::hideElement(selector = "#hide_data_comparison_gecco")
      shinyjs::hideElement(selector = "#ref_comparison_gecco")
      shinyjs::hideElement(selector = "#score_type_gecco")
      shinyjs::hideElement(selector = "#plot_step_gecco")
      shinyjs::hideElement(selector = "#plot_start_gecco")
      shinyjs::hideElement(selector = "#data_filter_header_gecco")
      shinyjs::hideElement(selector = "#hide_data_filter_gecco")
      shinyjs::hideElement(selector = "#score_average_gecco")
      shinyjs::hideElement(selector = "#score_cluster_gecco")
      shinyjs::hideElement(selector = "#domains_filter_gecco")
      shinyjs::hideElement(selector = "#prot_filter_gecco")
    }
  })
  # Ahow ARTS data options, if data is available
  shiny::observeEvent(vals$arts_data_input,{
    
    if (vals$arts_data_input == T){
      if (input$hide_anti == F){
        shinyjs::showElement(selector = "#arts_header")
        shinyjs::showElement(selector = "#dup_choice")
      }
      if (input$hide_viz == F){
        shinyjs::showElement(selector = "#arts_width")
      }
    } else {
      shinyjs::hideElement(selector = "#arts_header")
      shinyjs::hideElement(selector = "#dup_choice")
      shinyjs::hideElement(selector = "#arts_width")
    }
  })
  ##---------------------------------------------------------------
  ##              Data processing options show/hide               -
  ##---------------------------------------------------------------
  # Count data uploads, to show tabs and corresponding 
  # options 
  shiny::observeEvent(vals$data_upload_count, {
    
    if (vals$data_upload_count <2){
      shiny::hideTab("main", "2")
      shiny::hideTab("main", "3")
      
    }else if (vals$data_upload_count >=2){
      if (input$hide_summarize == F) {
        shinyjs::showElement(selector = "#summarize")
        shinyjs::showElement(selector = "#group_by")
        shinyjs::showElement(selector = "#count_all")
      }
      if (input$hide_viz == F){
        shinyjs::showElement(selector = "#biocircos_color")
        shinyjs::showElement(selector = "#label_color")
        shinyjs::showElement(selector = "#label_color_class")
      }
      shiny::showTab("main", "2")
      shiny::showTab("main", "3")
      if ((vals$gecco_data_input == T) & ((vals$anti_data_input == T) | (vals$prism_data_input == T) | (vals$sempi_data_input == T) )){
        shiny::showTab("main", "5")
      } else {
        shiny::hideTab("main", "5")
        if ((vals$gecco_data_input == T) & (vals$data_upload_count >1)){
          shiny::showNotification(paste("It seems that you would like to compare the GECCO data to the reference (in a new tab)? Please upload Antismash, SEMPI or PRISM datasets to do that."), type = "message", duration=10)
        }
      }
      if ((vals$deep_data_input == T) & ((vals$anti_data_input == T) | (vals$prism_data_input == T) | (vals$sempi_data_input == T) )) {
        shiny::showTab("main", "1")
      } else {
        shiny::hideTab("main", "1")
        if ((vals$deep_data_input == T) & (vals$data_upload_count >1)){
          shiny::showNotification(paste("It seems that you would like to compare the DeepBGC data to the reference (in a new tab)? Please upload Antismash, SEMPI or PRISM datasets to do that."), type = "message", duration=10)
        }
      }
    }
    if (vals$data_upload_count <1){
      shiny::hideTab("main", "4")
      shiny::hideTab(inputId = "main", target = "5")
      shiny::hideTab(inputId = "main", target = "1")
      shinyjs::hideElement(selector = "#genes_on_chr")
      shinyjs::hideElement(selector = "#hide_genes_on_chr")
      shinyjs::hideElement(selector = "#ref")
    }else{
      shiny::showTab("main", "4")
      if (input$hide_genes_on_chr == F){
        shinyjs::showElement(selector = "#genes_on_chr")
        shinyjs::showElement(selector = "#hide_genes_on_chr")
        shinyjs::showElement(selector = "#ref")
      }
      
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
  shiny::observeEvent(input$anti_hybrid, {
    
    hybrid_col <- function(data){
      data_split <- str_split(data$Type2, "__")
      types <- sapply(data_split, function(x){
        if (length(unlist(x))>1){
          "hybrid"
        } else{
          x
        }
      })
      return(types)
    }
    
    if (input$anti_hybrid==T){
      vals$anti_data$Type2 <- hybrid_col(vals$anti_data)
    }else {
      vals$anti_data$Type2 <- vals$anti_type
    }
    
    })
  shiny::observeEvent(input$prism_hybrid, {
    
    hybrid_col <- function(data){
      data_split <- str_split(data$Type2, "__")
      types <- sapply(data_split, function(x){
        if (length(unlist(x))>1){
          "hybrid"
        } else{
          x
        }
      })
      return(types)
    }
    
    if (input$prism_hybrid==T){
      vals$prism_data$Type2 <- hybrid_col(vals$prism_data)
    }else {
      vals$prism_data$Type2 <- vals$prism_type
    }
  })
  shiny::observeEvent(input$sempi_hybrid, {
    
    hybrid_col <- function(data){
      data_split <- str_split(data$Type2, "__")
      types <- sapply(data_split, function(x){
        if (length(unlist(x))>1){
          "hybrid"
        } else if (unlist(x) == 'nrps-pks'){
          "hybrid"
        } else {
          x
        }
      })
      return(types)
    }
    
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
      vals$anti_type <- rename_vector(anti_data, rename_data)
      anti_data['Type2'] <- vals$anti_type
      vals$anti_data <- anti_data
    }
    
    if (vals$sempi_data_input == T){
      sempi_data <- read.csv("sempi_data.csv")
      vals$sempi_type <- rename_vector(sempi_data, rename_data)
      sempi_data['Type2'] <- vals$sempi_type
      vals$sempi_data <- sempi_data
    }
    
    if(vals$prism_data_input == T){
      prism_data <- read.csv("prism_data.csv")
      vals$prism_type <- rename_vector(prism_data, rename_data)
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
      vals$anti_type <- rename_vector(anti_data, rename_data)
      anti_data['Type2'] <- vals$anti_type
      vals$anti_data <- anti_data
    }
    
    if (vals$sempi_data_input == T){
      sempi_data <- read.csv("sempi_data.csv")
      vals$sempi_type <- rename_vector(sempi_data, rename_data)
      sempi_data['Type2'] <- vals$sempi_type
      vals$sempi_data <- sempi_data
    }
    
    if(vals$prism_data_input == T){
      prism_data <- read.csv("prism_data.csv")
      vals$prism_type <- rename_vector(prism_data, rename_data)
      prism_data['Type2'] <-  vals$prism_type
      vals$prism_data <- prism_data
    }
  })
  # Reset the renaming. Uncheck the hybrid checkboxes
  shiny::observeEvent(input$reset_name, {
    
    vals$anti_data['Type2']  <- vals$anti_data['Type']
    vals$sempi_data['Type2'] <- vals$sempi_data['Type']
    vals$ prism_data['Type2'] <- vals$ prism_data['Type']
    shiny::updateCheckboxInput(inputId = "anti_hybrid", value = F)
    shiny::updateCheckboxInput(inputId = "sempi_hybrid", value =F)
    shiny::updateCheckboxInput(inputId = "prism_hybrid", value = F)
      shinyjs::showElement(selector = "#rename")
      shinyjs::hideElement(selector = "#reset_name")
    vals$renamed <- F
  })
  # Read the uploaded renaming scheme csv
  shiny::observeEvent(input$rename_data,{
    
    rename_data <- read.csv(input$rename_data$datapath)
    vals$rename_data <- rename_data
  })
  # What to do, if hide uploads scheme is triggered
  shiny::observeEvent(input$hide_uploads, {
    
    if (input$hide_uploads == T){
      shinyjs::hideElement(selector = "#anti_input_options")
      shinyjs::hideElement(selector = "#anti_data")
      shinyjs::hideElement(selector = "#prism_input_options")
      shinyjs::hideElement(selector = "#anti_header_upload")
      shinyjs::hideElement(selector = "#prism_header_upload")
      shinyjs::hideElement(selector = "#prism_data")
      shinyjs::hideElement(selector = "#sempi_header_upload")
      shinyjs::hideElement(selector = "#sempi_data")
      shinyjs::hideElement(selector = "#deep_header_upload")
      shinyjs::hideElement(selector = "#deep_data")
      shinyjs::hideElement(selector = "#gecco_header_upload")
      shinyjs::hideElement(selector = "#gecco_data")
      shinyjs::hideElement(selector = "#rre_header_upload")
      shinyjs::hideElement(selector = "#rre_data")
      shinyjs::hideElement(selector = "#chr_len")
      shinyjs::hideElement(selector = "#arts_header_upload")
      shinyjs::hideElement(selector = "#known_data")
      shinyjs::hideElement(selector = "#dup_data")
      shinyjs::hideElement(selector = "#anti_sco")
      shinyjs::hideElement(selector = "#prism_sco")
      shinyjs::hideElement(selector = "#arts_sco")
      shinyjs::hideElement(selector = "#rre_sco")
      shinyjs::hideElement(selector = "#sempi_sco")
      shinyjs::hideElement(selector = "#deep_sco")
      shinyjs::hideElement(selector = "#gecco_sco")
    }else {
      shinyjs::showElement(selector = "#anti_input_options")
      shinyjs::showElement(selector = "#anti_data")
      shinyjs::showElement(selector = "#prism_input_options")
      shinyjs::showElement(selector = "#anti_header_upload")
      shinyjs::showElement(selector = "#prism_header_upload")
      shinyjs::showElement(selector = "#prism_data")
      shinyjs::showElement(selector = "#sempi_header_upload")
      shinyjs::showElement(selector = "#sempi_data")
      shinyjs::showElement(selector = "#deep_header_upload")
      shinyjs::showElement(selector = "#deep_data")
      shinyjs::showElement(selector = "#gecco_header_upload")
      shinyjs::showElement(selector = "#gecco_data")
      shinyjs::showElement(selector = "#rre_header_upload")
      shinyjs::showElement(selector = "#rre_data")
      shinyjs::showElement(selector = "#chr_len")
      shinyjs::showElement(selector = "#arts_header_upload")
      shinyjs::showElement(selector = "#known_data")
      shinyjs::showElement(selector = "#dup_data")
      shinyjs::showElement(selector = "#anti_sco")
      shinyjs::showElement(selector = "#prism_sco")
      shinyjs::showElement(selector = "#arts_sco")
      shinyjs::showElement(selector = "#rre_sco")
      shinyjs::showElement(selector = "#sempi_sco")
      shinyjs::showElement(selector = "#deep_sco")
      shinyjs::showElement(selector = "#gecco_sco")
  }
    })
  # What to do, if hide data options scheme is triggered
  shiny::observeEvent(input$hide_anti, {
    
    if (input$hide_anti== T){
      shinyjs::hideElement(selector = "#anti_header")
      shinyjs::hideElement(selector = "#anti_hybrid")
      shinyjs::hideElement(selector = "#sempi_header")
      shinyjs::hideElement(selector = "#sempi_hybrid")
      shinyjs::hideElement(selector = "#prism_header")
      shinyjs::hideElement(selector = "#prism_hybrid")
      shinyjs::hideElement(selector = "#prism_supp")
      shinyjs::hideElement(selector = "#arts_header")
      shinyjs::hideElement(selector = "#dup_choice")
    }else{
      if (vals$anti_data_input == T){
        shinyjs::showElement(selector = "#anti_header")
        shinyjs::showElement(selector = "#anti_hybrid")
      } else{
        shinyjs::hideElement(selector = "#anti_header")
        shinyjs::hideElement(selector = "#anti_hybrid")
      }
      if (vals$prism_data_input == T){
      shinyjs::showElement(selector = "#prism_header")
      shinyjs::showElement(selector = "#prism_hybrid")
      if (vals$prism_json == T){
        shinyjs::showElement(selector = "#prism_supp")
      }
      } else {
        shinyjs::hideElement(selector = "#prism_header")
        shinyjs::hideElement(selector = "#prism_hybrid")
        shinyjs::hideElement(selector = "#prism_supp")
      }
      if (vals$sempi_data_input == T){
      shinyjs::showElement(selector = "#sempi_header")
      shinyjs::showElement(selector = "#sempi_hybrid")
      } else {
        shinyjs::hideElement(selector = "#sempi_header")
        shinyjs::hideElement(selector = "#sempi_hybrid")
      }
      if (vals$arts_data_input == T){
        shinyjs::showElement(selector = "#arts_header")
        shinyjs::showElement(selector = "#dup_choice")
      } else{
        shinyjs::hideElement(selector = "#arts_header")
        shinyjs::hideElement(selector = "#dup_choice")
      }
    }
  })
  # What to do, if hide annotation plot options scheme is triggered
  shiny::observeEvent(input$hide_genes_on_chr, {
    
    if (input$hide_genes_on_chr == T){
      shinyjs::hideElement(selector = "#ref")
    } else {
      if (vals$data_upload_count > 0){
        shinyjs::showElement(selector = "#ref")
      } else {
        shinyjs::hideElement(selector = "#genes_on_chr")
        shinyjs::hideElement(selector = "#ref")
      }
    }
  })
  # What to do, if hide summarize tab options scheme is triggered
  shiny::observeEvent(input$hide_summarize, {
    
    if (input$hide_summarize == T){
      shinyjs::hideElement(selector = "#group_by")
      shinyjs::hideElement(selector = "#count_all")
    } else {
      if (vals$data_upload_count > 1){
        shinyjs::showElement(selector = "#group_by")
        shinyjs::showElement(selector = "#count_all")
      } else {
        shinyjs::hideElement(selector = "#summarize")
        shinyjs::hideElement(selector = "#group_by")
        shinyjs::hideElement(selector = "#count_all")
      }
      
    }
  })
  # What to do, if hide improve visualization scheme is triggered
  shiny::observeEvent(input$hide_viz, {
    
    if (input$hide_viz == T){
      shinyjs::hideElement(selector = "#rename_data")
      shinyjs::hideElement(selector = "#rename")
      shinyjs::hideElement(selector = "#reset_name")
      shinyjs::hideElement(selector = "#rre_width")
      shinyjs::hideElement(selector = "#biocircos_color")
      shinyjs::hideElement(selector = "#label_color")
      shinyjs::hideElement(selector = "#label_color_class")
      shinyjs::hideElement(selector = "#ref_col_biocircos")
      shinyjs::hideElement(selector = "#arts_header")
      shinyjs::hideElement(selector = "#arts_width")
      shinyjs::hideElement(selector = "#sempi_width")
      shinyjs::hideElement(selector = "#prism_supp_data_input_width")
    } else{
      shinyjs::showElement(selector = "#rename_data")
      shinyjs::showElement(selector = "#rename")
      shinyjs::showElement(selector = "#reset_name")
      if (vals$rre_data_input == T){
        shinyjs::showElement(selector = "#rre_width")
      } else {
        shinyjs::hideElement(selector = "#rre_width")
      }
      if (vals$sempi_data_input == T){
        shinyjs::showElement(selector = "#sempi_width")
      } else {
        shinyjs::hideElement(selector = "#sempi_width")
      }
      if (vals$prism_json == T){
        shinyjs::showElement(selector = "#prism_supp_data_input_width")
      }
      else {
        shinyjs::hideElement(selector = "#prism_supp_data_input_width")
      }
      if (vals$data_upload_count > 1){
        shinyjs::showElement(selector = "#biocircos_color")
        shinyjs::showElement(selector = "#label_color")
        shinyjs::showElement(selector = "#label_color_class")
      } else {
        shinyjs::hideElement(selector = "#biocircos_color")
        shinyjs::hideElement(selector = "#label_color")
        shinyjs::hideElement(selector = "#label_color_class")
      }
      if (input$label_color_class == "R"){
        shinyjs::showElement(selector = "#ref_col_biocircos")
      } else {
        shinyjs::hideElement(selector = "#ref_col_biocircos")
      }
      if (vals$arts_data_input == T){
        shinyjs::showElement(selector = "#arts_header")
        shinyjs::showElement(selector = "#arts_width")
      } else {
        shinyjs::hideElement(selector = "#arts_header")
        shinyjs::hideElement(selector = "#arts_width")
      }
      
    }
  })
  # What to do, if hide DeepBGC comparison options scheme is triggered
  shiny::observeEvent(input$hide_data_comparison, {
    
    if ((input$hide_data_comparison == T)){
      shinyjs::hideElement(selector = "#ref_comparison")
      shinyjs::hideElement(selector = "#score_type")
      shinyjs::hideElement(selector = "#plot_step")
      shinyjs::hideElement(selector = "#plot_start")
    } else if (vals$deep_data_input == T) {
      shinyjs::showElement(selector = "#ref_comparison")
      shinyjs::showElement(selector = "#score_type")
      shinyjs::showElement(selector = "#plot_step")
      shinyjs::showElement(selector = "#plot_start")
    } else {
      shinyjs::hideElement(selector = "#ref_comparison")
      shinyjs::hideElement(selector = "#score_type")
      shinyjs::hideElement(selector = "#plot_step")
      shinyjs::hideElement(selector = "#plot_start")
    }
  })
  # What to do, if hide DeepBGC filtering options scheme is triggered
  shiny::observeEvent(input$hide_data_filter, {
    
    if ((input$hide_data_filter == T)){
      shinyjs::hideElement(selector = "#score_a")
      shinyjs::hideElement(selector = "#score_d")
      shinyjs::hideElement(selector = "#score_c")
      shinyjs::hideElement(selector = "#domains_filter")
      shinyjs::hideElement(selector = "#biodomain_filter")
      shinyjs::hideElement(selector = "#gene_filter")
      shinyjs::hideElement(selector = "#cluster_type")
    } else if  (vals$deep_data_input == T){
      shinyjs::showElement(selector = "#score_a")
      shinyjs::showElement(selector = "#score_d")
      shinyjs::showElement(selector = "#score_c")
      shinyjs::showElement(selector = "#domains_filter")
      shinyjs::showElement(selector = "#biodomain_filter")
      shinyjs::showElement(selector = "#gene_filter")
      shinyjs::showElement(selector = "#cluster_type")
    } else {
      shinyjs::hideElement(selector = "#score_a")
      shinyjs::hideElement(selector = "#score_d")
      shinyjs::hideElement(selector = "#score_c")
      shinyjs::hideElement(selector = "#domains_filter")
      shinyjs::hideElement(selector = "#biodomain_filter")
      shinyjs::hideElement(selector = "#gene_filter")
      shinyjs::hideElement(selector = "#cluster_type")
    }
  })
  # What to do, if hide GECCO comparison options scheme is triggered
  shiny::observeEvent(input$hide_data_comparison_gecco, {
    
    if ((input$hide_data_comparison_gecco == T)){
      shinyjs::hideElement(selector = "#ref_comparison_gecco")
      shinyjs::hideElement(selector = "#score_type_gecco")
      shinyjs::hideElement(selector = "#plot_step_gecco")
      shinyjs::hideElement(selector = "#plot_start_gecco")
    } else if (vals$gecco_data_input == T) {
      shinyjs::showElement(selector = "#ref_comparison_gecco")
      shinyjs::showElement(selector = "#score_type_gecco")
      shinyjs::showElement(selector = "#plot_step_gecco")
      shinyjs::showElement(selector = "#plot_start_gecco")
    } else {
      shinyjs::hideElement(selector = "#ref_comparison_gecco")
      shinyjs::hideElement(selector = "#score_type_gecco")
      shinyjs::hideElement(selector = "#plot_step_gecco")
      shinyjs::hideElement(selector = "#plot_start_gecco")
    }
  })
  # What to do, if hide GECCO filtering options scheme is triggered
  shiny::observeEvent(input$hide_data_filter_gecco, {
    
    if ((input$hide_data_filter_gecco == T)){
      shinyjs::hideElement(selector = "#score_average_gecco")
      shinyjs::hideElement(selector = "#score_average_gecco")
      shinyjs::hideElement(selector = "#domains_filter_gecco")
      shinyjs::hideElement(selector = "#prot_filter_gecco")
    } else if  (vals$gecco_data_input == T){
      shinyjs::showElement(selector = "#score_average_gecco")
      shinyjs::showElement(selector = "#score_average_gecco")
      shinyjs::showElement(selector = "#domains_filter_gecco")
      shinyjs::showElement(selector = "#prot_filter_gecco")
    } else {
      shinyjs::hideElement(selector = "#score_average_gecco")
      shinyjs::hideElement(selector = "#score_average_gecco")
      shinyjs::hideElement(selector = "#domains_filter_gecco")
      shinyjs::hideElement(selector = "#prot_filter_gecco")
    }
  })
  
  
  ############################################################################
  ############################################################################
  ###                                                                      ###
  ###                             COMPUTATIONS                             ###
  ###                                                                      ###
  ############################################################################
  ############################################################################
  # Compute all interceptions on data upload.
  # dplyr::filter while ploting then.
  # TODO make looop for data reading
  shiny::observeEvent(inputData(), {
    shiny::req(vals$data_upload_count>=1)
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
        dplyr::select(nucl_start, nucl_end)
      
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
      if (input$prism_supp_data_input_width == TRUE) {
        Stop_vals_prism_supp <- as.numeric(vals$prism_supp_data$Stop)+50000
      } else{
        Stop_vals_prism_supp <- as.numeric(vals$prism_supp_data$Stop)
      }
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
      query <- makeGRangesFromDataFrame(inter2)
      subject <- makeGRangesFromDataFrame(inter1)
      interseption <- findOverlaps(query,subject)
      inter_from <- interseption@from
      inter_to <- interseption@to
      return(list(from = inter_from, to = inter_to))
    }
    
    inters <- vals$inters
    data_uploads <- c("anti_data_input","sempi_data_input","prism_data_input","prism_supp_data_input",
                      "arts_data_input","deep_data_input","gecco_data_input","rre_data_input")
    soft_names <- c("anti","sempi","prism","prism_supp","arts","deep","gecco","rre" )
    

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
  shiny::observeEvent(dynamicInput(), {
    shiny::req(vals$data_upload_count>=1)
    inters <- vals$inters
    if (vals$deep_data_input == TRUE){
      if (vals$need_filter == F) {
        biocircos_deep <- filter_deepbgc()
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
        gecco_data <- filter_gecco()
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
          dplyr::filter(Core == str_split(str_split(input$dup_choice, " ,")[[1]][[2]], "Core:")[[1]][[2]] | Core == "Not_core")
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
    if ((vals$gecco_data_input == F) & (vals$deep_data_input == F )& (vals$arts_data_input == F )) {
      vals$inters_filtered <- inters
    }
    vals$need_filter <- F
    vals$filter_data <- F
    vals$can_plot_deep_ref = T
    enable_event_logic()
    
  })
  # Compute the Biociros plot. Store information to plot later
  shiny::observeEvent(biocircos_listen(), {
    shiny::req(vals$data_upload_count >=2)
    shiny::req(vals$need_filter == F)
    shiny::req(vals$can_plot_biocircos == T)
    initialize_biocircos <- function(biocircos_anti, name,Biocircos_chromosomes, arcs_chromosomes, arcs_begin , arcs_end, arc_labels, arc_col, rename_data ){
      #Make chromosome list for Biocircos plot. Use chr_len as an input
      Biocircos_chromosomes[[name]] <- vals$chr_len  
      #Add arcs. Quantity of arcs is length of dataframes
      arcs_chromosomes <- c(arcs_chromosomes,rep(name, length(biocircos_anti$Cluster)) )
      # Add arcs begin positions. (Start column)
      arcs_begin <- c(arcs_begin, biocircos_anti$Start)
      # Stop position of arcs. 
      biocircos_anti <- correct_width(biocircos_anti, name)
      arcs_end <- c(arcs_end, biocircos_anti$Stop )
      # Add Arcs labels. Can add only one label...
      arc_labels <- c(arc_labels, biocircos_anti$Type)
      if ((input$biocircos_color == T)){
        arc_colors <- sapply(biocircos_anti$Type2, function(x){
          if (x %in% rename_data$Group_color){
            rename_data$Color[rename_data$Group_color == x]
          } else {
            rename_data$Color[rename_data$Group_color == 'base']
          }
        })
      } else {
        arc_colors <- rename_data$Color[rename_data$Group_color == 'base']
      }
      arc_col <- c(arc_col,as.character(arc_colors) )
      return(list(Biocircos_chromosomes,arcs_chromosomes, arcs_begin , arcs_end, arc_labels, arc_col))
    }
    
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
    data_uploads <- c("anti_data_input","sempi_data_input","prism_data_input","prism_supp_data_input",
                      "arts_data_input","deep_data_input","gecco_data_input","rre_data_input")
    soft_names <- c("anti","sempi","prism","prism_supp","arts","deep","gecco","rre" )
    soft <- c("Antismash","SEMPI","PRISM","PRISM-Supp","ARTS","DeepBGC","GECCO","RRE-Finder" )
   data_to_use <- c( "anti_data" ,"sempi_data" , "prism_data", "prism_supp_data","arts_data_filtered","deep_data_filtered" ,"gecco_data_filtered", "rre_data")
    index <- 1
   # browser()
    for (upload in data_uploads){
      if (vals[[upload]] == T){
        # Store data in local variable
        init_data <- initialize_biocircos(vals[[data_to_use[index]]], soft[index], Biocircos_chromosomes, arcs_chromosomes, arcs_begin , arcs_end, arc_labels, arc_col, rename_data )
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
    
    
    add_biocircos_data <- function(data1_inter, data2_inter, data1, data2, data1_label, data2_label, rename_data, class){
      inter_s_rre_n <- data1_inter
      inter_rre_s <- data2_inter
      # Add link start. Just populate certain chromosome name times the lenght of interception 
      chromosomes_start <- c(rep(data2_label, length(inter_rre_s)))
      # Add link end. Just populate second output from the vectors, used above. 
      chromosomes_end <- c(rep(data1_label, length(inter_s_rre_n)))
      # Add links start positions as a start from dataframe. This vector is for chromosome start
      link_pos_start <- as.numeric(c(data2$Start[match(inter_rre_s,data2$Cluster)] ))
      # Add links start positions as a start from dataframe. For chromosome start variable
      link_pos_start_1 <- as.numeric(c(data2$Stop[match(inter_rre_s,data2$Cluster)]))
      # Add links start position for a chromosome stop variable
      link_pos_end <- as.numeric(c( data1$Start[match(inter_s_rre_n,data1$Cluster)]))
      # Add links start position for a chromosome stop position
      link_pos_end_2 <- as.numeric(c(data1$Stop[match(inter_s_rre_n,data1$Cluster)]))
      label_1 <- c(sapply(inter_rre_s, function(x){x = paste(paste0(data2_label,":"), x, ",", data2$Type[data2$Cluster == x])})) 
      label_2 <- c(sapply(inter_s_rre_n, function(x){x = paste(paste0(data1_label, ":"), x, ",", data1$Type[data1$Cluster == x])}))
      #browser()
      if (!is.null(inter_rre_s)){
        if (class == 'P'){
          subset_vec <- data2$Type2[match(inter_rre_s, data2$Cluster)] == data1$Type2[match(inter_s_rre_n, data1$Cluster)]
          label_color <- as.character(c(sapply(data2$Type2[match(inter_rre_s, data2$Cluster )], function (x){
            if (x %in% rename_data$Group_color) {
              as.character(rename_data$Color[rename_data$Group_color == x])
            }
            else{
              as.character(rename_data$Color[rename_data$Group_color == 'base'])
            }
          })))
          if (length(label_color) != 0){
            for (t in seq(1:length(label_color))){
              if (!is.null(subset_vec[t])){
                if (subset_vec[t] == F){
                  label_color[t] <- as.character(rename_data$Color[rename_data$Group_color == 'base'])
                }
              }
            }
          }
        } else if (class == 'H'){
          if (grep(paste0("^", data1_label, "$"), rename_data$Hierarchy) < (grep(paste0("^", data2_label, "$"), rename_data$Hierarchy))){
            label_color <- as.character(c(sapply(data1$Type2[match(inter_s_rre_n, data1$Cluster )], function (x){
              if (x %in% rename_data$Group_color) {
                as.character(rename_data$Color[rename_data$Group_color == x])
              }
              else{
                as.character(rename_data$Color[rename_data$Group_color == 'base'])
              }
            })))
          } else {
            label_color <-as.character( c(sapply(data2$Type2[ match(inter_rre_s, data2$Cluster )], function (x){
              if (x %in% rename_data$Group_color) {
                as.character(rename_data$Color[rename_data$Group_color == x])
              }
              else{
                as.character(rename_data$Color[rename_data$Group_color == 'base'])
              }
            })))
          }
        }else if (class == 'R'){
          if (data2_label == input$ref_col_biocircos){
            label_color <- as.character(c(sapply(data1$Type2[match(inter_s_rre_n, data1$Cluster)], function (x){
              if (x %in% rename_data$Group_color) {
                as.character(rename_data$Color[rename_data$Group_color == x])
              }
              else{
                as.character(rename_data$Color[rename_data$Group_color == 'base'])
              }
            })))
          } else if (data1_label == input$ref_col_biocircos){
            label_color <- as.character(c(sapply(data2$Type2[match(inter_rre_s,data2$Cluster)], function (x){
              if (x %in% rename_data$Group_color) {
                as.character(rename_data$Color[rename_data$Group_color == x])
              }
              else{
                as.character( rename_data$Color[rename_data$Group_color == 'base'])
              }
            })))
          } else{
            label_color <- as.character(rep(rename_data$Color[rename_data$Group_color == 'base'], length(chromosomes_start)))
          }
        } else {
          label_color <-as.character( rep(rename_data$Color[rename_data$Group_color == 'base'], length(chromosomes_start)))
        }
      }
      return(list( inter_s_rre_n, inter_s_rre_n, chromosomes_start, chromosomes_end, link_pos_start, link_pos_start_1, link_pos_end, 
                   link_pos_end_2, label_1, label_2, label_color))
    }
    data_uploads_2 <- data_uploads
    soft_2 <- soft
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
          if ((vals[[upload2]]==T) & (length(data_uploads_2) > 0) & (soft[index] != soft_2[index2])){
            output <- add_biocircos_data(inters[[soft_names[index]]][[soft_names_2[index2]]]$from, inters[[soft_names[index]]][[soft_names_2[index2]]]$to, vals[[data_to_use_2[index2]]], vals[[data_to_use[index]]], soft_2[index2], soft[index], rename_data, input$label_color_class)
            
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
                                                displayLabel = FALSE, color = rename_data$Color[rename_data$Group_color == 'base'])
    } else{
      shiny::showNotification(paste("No interceptions are being made in the Biocircos plot. Please provide data with clusters that do have intercepting borders"), type = "warning", duration=NULL)
    }
    
    vals$tracklist <- tracklist
    vals$Biocircos_chromosomes <- Biocircos_chromosomes
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
      if (input$ref_comparison == 'A'){
        anti_inter <- vals$anti_data %>%
        dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison == 'P'){
        anti_inter <- vals$prism_data %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison == 'S'){
        anti_inter <- vals$sempi_data %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } 
      
      
     
      
      # Get the interception of two matrices
      if (length(deep_inter$Start) > 0) {
        query <- makeGRangesFromDataFrame(deep_inter)
        subject <- makeGRangesFromDataFrame(anti_inter)
        interseption <- findOverlaps(query,subject)
        inter_bgc <- length(interseption@from)
        len_new <- length(deep_inter$seqnames) - inter_bgc
      } else {
        inter_bgc <- 0
        len_new <- 0
      }
      

      if (input$ref_comparison == 'A'){
        used_antismash <-  length(vals$anti_data$Cluster)-inter_bgc
        cols <-  c("Only Antismash", "DeepBGC+Antismash", "Only DeepBGC")
        title <-  ggplot2::ggtitle("Comparison of Antismash and DeepBGC annotations at given score threshold")
      } else if (input$ref_comparison == 'P'){
        used_antismash <-  length(vals$prism_data$Cluster)-inter_bgc
        cols <- c("Only PRISM", "DeepBGC+PRISM", "Only DeepBGC")
        title <- ggplot2::ggtitle("Comparison of PRISM and DeepBGC annotations at given score threshold")
      } else if (input$ref_comparison == 'S') {
        used_antismash <-  length(vals$sempi_data$Cluster)-inter_bgc
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
    if (input$ref_comparison == 'A'){
      data <-  vals$anti_data
      title <- ggplot2::ggtitle("Rates of DeepBGC/Antismash data annotation")
      test <- test %>%
        # Calculate rates. Novelty is nummber of clusters annotated only by deepbgc/ all clusters annotated by antismash + (antismash + deepbgc)
        dplyr::mutate(Novelty_rate = test$`Only DeepBGC`/(test$`DeepBGC+Antismash` + test$`Only Antismash`), 
               #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
               Annotation_rate = test$`DeepBGC+Antismash`/length(data$Cluster), 
               # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
               Skip_rate = test$`Only Antismash`/length(data$Cluster))
    } else if (input$ref_comparison == 'P'){
      data <- vals$prism_data
      title <- ggplot2::ggtitle("Rates of DeepBGC/PRISM data annotation")
      test <- test %>%
        dplyr::mutate(Novelty_rate = test$`Only DeepBGC`/(test$`DeepBGC+PRISM` + test$`Only PRISM`), 
               #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
               Annotation_rate = test$`DeepBGC+PRISM`/length(data$Cluster), 
               # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
               Skip_rate = test$`Only PRISM`/length(data$Cluster))
    } else if (input$ref_comparison == 'S'){
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
      if (input$ref_comparison_gecco == 'A'){
        anti_inter <- vals$anti_data %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison_gecco == 'P'){
        anti_inter <- vals$prism_data %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison_gecco == 'S'){
        anti_inter <- vals$sempi_data %>%
          dplyr::select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } 
      
      
      
      
      # Get the interception of two matrices
      if (length(gecco_inter$Start) > 0) {
        query <- makeGRangesFromDataFrame(gecco_inter)
        subject <- makeGRangesFromDataFrame(anti_inter)
        interseption <- findOverlaps(query,subject)
        inter_bgc <- length(interseption@from)
        len_new <- length(gecco_inter$seqnames) - inter_bgc
      } else {
        inter_bgc <- 0
        len_new <- 0
      }
      
      
      if (input$ref_comparison_gecco == 'A'){
        used_antismash <-  length(vals$anti_data$Cluster)-inter_bgc
        cols <-  c("Only Antismash", "GECCO+Antismash", "Only GECCO")
        title <-  ggplot2::ggtitle("Comparison of Antismash and GECCO annotations at given score threshold")
      } else if (input$ref_comparison_gecco == 'P'){
        used_antismash <-  length(vals$prism_data$Cluster)-inter_bgc
        cols <- c("Only PRISM", "GECCO+PRISM", "Only GECCO")
        title <- ggplot2::ggtitle("Comparison of PRISM and GECCO annotations at given score threshold")
      } else if (input$ref_comparison_gecco == 'S') {
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
    if (input$ref_comparison_gecco == 'A'){
      data <-  vals$anti_data
      title <- ggplot2::ggtitle("Rates of GECCO/Antismash data annotation")
      test <- test %>%
        # Calculate rates. Novelty is nummber of clusters annotated only by deepbgc/ all clusters annotated by antismash + (antismash + deepbgc)
        dplyr::mutate(Novelty_rate = test$`Only GECCO`/(test$`GECCO+Antismash` + test$`Only Antismash`), 
               #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
               Annotation_rate = test$`GECCO+Antismash`/length(data$Cluster), 
               # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
               Skip_rate = test$`Only Antismash`/length(data$Cluster))
    } else if (input$ref_comparison_gecco == 'P'){
      data <- vals$prism_data
      title <- ggplot2::ggtitle("Rates of GECCO/PRISM data annotation")
      test <- test %>%
        dplyr::mutate(Novelty_rate = test$`Only GECCO`/(test$`GECCO+PRISM` + test$`Only PRISM`), 
               #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
               Annotation_rate = test$`GECCO+PRISM`/length(data$Cluster), 
               # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
               Skip_rate = test$`Only PRISM`/length(data$Cluster))
    } else if (input$ref_comparison_gecco == 'S'){
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
    shiny::req(vals$data_upload_count >=1)
    shiny::req(vals$need_filter == F)
    shiny::req(vals$can_plot_deep_ref == T)
    
    shiny::req(vals$data_upload_count >=1)
    data_uploads <- c("anti_data_input","sempi_data_input","prism_data_input","prism_supp_data_input",
                      "arts_data_input","deep_data_input","gecco_data_input","rre_data_input")
    soft_names <- c("anti","sempi","prism","prism_supp","arts","deep","gecco","rre" )
    soft <- c("Antismash","SEMPI","PRISM","PRISM_SUPPORT","ARTS","DeepBGC","GECCO","RRE-Finder" )
    abbr <- c("A", "S", "P", "P-supp", "AR", "D", "G", "RRE")
    soft_ref <- c("Antismash","SEMPI","PRISM","PRISM-supp","ARTS","DeepBGC","GECCO","RRE-Finder" )
    soft_width <-  c("Antismash","SEMPI","PRISM","PRISM-Supp","ARTS","DeepBGC","GECCO","RRE-Finder" )
    data_to_use <- c( "anti_data" ,"sempi_data" , "prism_data", "prism_supp_data","arts_data_filtered","deep_data_filtered" ,
                      "gecco_data_filtered", "rre_data")
    soft_datafr <- c("seg_df_ref_a", "seg_df_ref_s" , "seg_df_ref_p", "seg_df_ref_p_s", "seg_df_ref_ar", "seg_df_ref_d", 
                     "seg_df_ref_g", "seg_df_ref_r")
    
    inters <- vals$inters_filtered
    # GENERATE DATA
    index <- 1
    for (upload in data_uploads){
      if (vals[[upload]] == T){
      data<- vals[[data_to_use[index]]]
      assign(paste0(soft_names[index], "_data"),  correct_width(data, soft_width[index]))
      }
      index <- index +1
    }

    
    lett <- rev(LETTERS)[1:9]
    simple_seg <- function(df, letter, software, soft_name ,soft, inter=T){
      if (inter== T){
        data <- df[df$Cluster %in% inters[[soft]][[soft_name]]$from, ]
      } else{
        data <- df
      }

      seg_df <-  data.frame(x=as.numeric(data$Start),
                         y=rep(letter, length(data$Cluster)),
                         xend=as.numeric(  data$Stop),
                         yend=rep(letter, length(data$Cluster)),
                         Type = as.factor(data$Type),
                         Type2 = as.factor(data$Type2),
                         Software = rep(software, length(data$Cluster)),
                         ID = data$Cluster,
                         Start = data$Start,
                         Stop = data$Stop)
      return(seg_df)
      }
    
    geom_anti <- function(data){
      ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, 
                                  ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    geom_prism <- function(data){
      ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,                                                                                                             ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    geom_deep <- function(data){
      ggplot2::geom_segment(data=data,ggplot2::aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                               ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
                                                               deepbgc_score = deepbgc_score,activity = activity ),size =3)
      }
    geom_rre <- function(data){
      if (vals$rre_more == T){
      ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type, Score = Score, Software = Software,
                                                           ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value,
                                                           P_value = P_value, RRE_start = RRE_start,RRE_stop = RRE_stop, 
                                                           Probability = Probability),size = 3)
      } else {
        ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                    ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value
                                    ),size = 3)
      }
      }
    geom_sempi <- function(data){

      ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,
                                  ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    geom_prism_supp <- function(data){
      ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, ID = ID,
                                  Start = Start, Stop = Stop, Type = Type, Name = Name, Full_name = Full_name,
                                  Score = Score), size = 3)
    }
    geom_arts <- function(data){
      ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,                                                                                                             ID = ID, Start = Start, Stop = Stop, Type = Type, Hit = Hit, 
                                  Core = Core, E_value = E_value, Bitscore = Bitscore, Count = Count, Model = Model), size = 3)
    }
    geom_gecco <- function(data){
      ggplot2::geom_segment(data=data, ggplot2::aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, 
                                  ID = ID, Start = Start, Stop = Stop, Type = Type, Num_proteins= Num_proteins,
                                  Num_domains = Num_domains,Average_p = Average_p, Max_p = Max_p ), size = 3)
    }
    
    tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                "P_value", "RRE_start","RRE_stop", "Probability", "Name", "Full_name",  "Hit", "Core", "Count", "Bitscore", "Model",
                "Num_domains", "Num_proteins", "Average_p", "Max_p")

    
    add_sempi <- function(seg_df, soft, df, inter = T){
      
      return(seg_df)
    }
    
    add_arts <- function(seg_df, soft, df, inter=T){
      if (inter == T){
        subset_df <- df[df$Cluster %in% inters[[soft]]$arts$from, ]
      }else{
        subset_df <- df
      }
      seg_df$Hit = subset_df$Hit
      seg_df$xend = as.numeric(subset_df$Stop)
      seg_df$Core = subset_df$Core
      seg_df$Count = subset_df$Count
      seg_df$E_value = subset_df$Evalue
      seg_df$Bitscore = subset_df$Bitscore
      seg_df$Model = subset_df$Model
      return(seg_df)
    }
    add_prism_supp <- function(seg_df, soft, df, inter=T){
      if (inter == T){
        subset_df <- df[df$Cluster %in% inters[[soft]]$prism_supp$from, ]
      }else{
        subset_df <- df
      }
      seg_df$xend <-  as.numeric(subset_df$Stop)
      seg_df$Score = subset_df$Score
      seg_df$Name = subset_df$Name
      seg_df$Full_name = subset_df$Full_name
      return(seg_df)
    }
    add_deep <- function(seg_df, soft, df, inter=T){
      if (inter == T){
        subset_df <- df[df$Cluster %in% inters[[soft]]$deep$from, ]
      }else{
        subset_df <- df
      }
      seg_df$num_domains = subset_df$num_domains
      seg_df$deepbgc_score = subset_df$deepbgc_score
      seg_df$activity = subset_df$product_activity
      return(seg_df)
    }
    add_rre <- function(seg_df, soft, df, inter=T){
      if (inter == T){
        subset_df <- df[df$Cluster %in% inters[[soft]]$rre$from, ]
      }else{
        subset_df <- df
      }
      if (vals$rre_more == T){
        seg_df$xend=as.numeric(subset_df$Stop)
        seg_df$Score = subset_df$Score
        seg_df$Stop = subset_df$Stop
        seg_df$E_value = subset_df$E.value
        seg_df$P_value = subset_df$P.value
        seg_df$RRE_start = subset_df$RRE.start
        seg_df$RRE_stop = subset_df$RRE.end
        seg_df$Probability = subset_df$Probability
      } else {
        seg_df$xend=subset_df$Stop
        seg_df$E_value = subset_df$E.value
      }
      
     return(seg_df)
    }
    add_gecco <- function(seg_df, soft, df, inter=T){
      if (inter == T){
        subset_df <- df[df$Cluster %in% inters[[soft]]$gecco$from, ]
      }else{
        subset_df <- df
      }
      seg_df$Num_proteins = subset_df$num_prot
      seg_df$Num_domains = subset_df$num_domains
      seg_df$Average_p = subset_df$average_p
      seg_df$Max_p = subset_df$max_p
      return(seg_df)
    }
    
    add_more_annot <- function(seg_df, plot, soft_names, index){
      if (dim(seg_df)[1] > 0){
        if (soft_names[index] == "anti"){
          plot <- plot + geom_anti(seg_df)
        } else if (soft_names[index] == "sempi") {
          plot <- plot + geom_sempi(seg_df)
        }
        else if (soft_names[index] == "prism") {
          plot <- plot + geom_prism(seg_df)
        }
        else if (soft_names[index] == "prism_supp") {
          plot <- plot + geom_prism_supp(seg_df)
        }
        else if (soft_names[index] == "arts") {
          plot <- plot + geom_arts(seg_df)
        }
        else if (soft_names[index] == "deep") {
          plot <- plot + geom_deep(seg_df)
        }
        else if (soft_names[index] == "rre") {
          plot <- plot + geom_rre(seg_df)
        } else if (soft_names[index] == "gecco") {
          plot <- plot+geom_gecco(seg_df)
        }
        return(plot)
      } else{
         return(plot)
      }
    }
    
    define_spec_seg_df <- function(soft_names, index,seg_df, soft_major, df , inter=T){
      if (inter == F){
        soft_major <- "Not applicable"
      }
      if ((soft_names[index] == "prism_supp") & (soft_names[index] != soft_major)){
        seg_df <-  add_prism_supp(seg_df, soft_major, df, inter)
      } else if ((soft_names[index] == "arts")& (soft_names[index] != soft_major)){
        seg_df <- add_arts(seg_df, soft_major, df, inter)
      } else if ((soft_names[index] == "deep")& (soft_names[index] != soft_major)){
        seg_df <- add_deep(seg_df, soft_major, df, inter)
      } else if ((soft_names[index] == "gecco")& (soft_names[index] != soft_major)){
        seg_df <- add_gecco(seg_df, soft_major, df, inter)
      } else if ((soft_names[index] == "rre")& (soft_names[index] != soft_major)){
        seg_df <- add_rre(seg_df, soft_major, df, inter)
      }else if ((soft_names[index] == "sempi")& (soft_names[index] != soft_major)){
        seg_df <- add_sempi(seg_df, soft_major, df, inter)
      }
      return(seg_df)
    }
    
# MAKE COMPUTATIONS
      sup_index <- 1
      soft_lttrs <- lett
      rename_y_axis <- vals$rename_y_axis
      rename_y_axis <- lapply(1:(length( soft_lttrs)-1), function(x){
        soft_lttrs[x]=soft[x]
      })
      names(rename_y_axis) <- soft_lttrs[-length(soft_lttrs)]
      for (upload in data_uploads){
        soft_lttr <- soft_lttrs[1]
        soft_lttrs <- soft_lttrs[-1]
        if (vals[[upload]] == T){
          soft_major <- soft_names[sup_index]
          seg_ref_g <- simple_seg(eval(as.name(paste(soft_names[sup_index], "_data", sep = ""))), "Z", soft[sup_index], soft_names[sup_index],soft_major, inter = F)
          seg_ref_g <- define_spec_seg_df(soft_names, sup_index,seg_ref_g, soft_major, eval(as.name(paste(soft_names[sup_index], "_data", sep = ""))), inter = F)
          seg_ref <- seg_ref_g
          if (input$ref == soft_ref[sup_index]){
            plot <- ggplot2::ggplot(eval(as.name(paste(soft_names[sup_index], "_data", sep = ""))), ggplot2::aes(x = vals$chr_len, y = Chr)) + 
              eval(as.name(paste0("geom_", soft_names[sup_index])))(seg_ref)
            soft_let <- abbr[sup_index]
            lettrs <- lett[2:length(lett)]
            labels_1 <- list()
            index = 1
            for (i in data_uploads){
              if ((vals[[i]] == T) & (soft_names[index] != soft_major)){
                df <- eval(as.name(paste(soft_names[index], "_data", sep = "")))
                seg_df <- simple_seg(df, lettrs[index], soft[index], soft_names[index],soft_major)
                seg_df <- define_spec_seg_df(soft_names, index,seg_df, soft_major, df)
                labels_1[[lettrs[index]]] <- (paste(abbr[index], "_vs_", soft_let, sep = ""))
                plot <- add_more_annot(seg_df, plot, soft_names, index)
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
      vals$can_plot_deep_ref_2 =T
    to_plot
  })
  
  output$deep_reference_2 <- plotly::renderPlotly({
    shiny::req(vals$can_plot_deep_ref_2 == T)
    vals$can_plot_deep_ref_2 == F
    shiny::req(vals$data_upload_count >=1)
    rename_y_axis <- vals$rename_y_axis
    data <- NULL
    data_uploads <- c("anti_data_input","sempi_data_input","prism_data_input","prism_supp_data_input",
                      "arts_data_input","deep_data_input","gecco_data_input","rre_data_input")
    data_to_use <- c( "anti_data" ,"sempi_data" , "prism_data", "prism_supp_data","arts_data_filtered","deep_data_filtered" ,
                      "gecco_data_filtered", "rre_data")

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
      title=list(text='<b> Cluster Types </b>')))
  })
  
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
    
    plot_data <- vals$rename_data
    new_data <- tidyr::drop_na(data.frame(cbind(as.character(plot_data$Group_color), as.character(plot_data$Color))) )
    new_data <- new_data[!apply(new_data == "", 1, all),]
    colnames(new_data) <- c("Name", "Color")
    color_vec <- new_data$Color
    options(DT.options = list(pageLength = 50))
    DT::datatable(new_data, rownames = F) %>% DT::formatStyle('Color',
                                                      backgroundColor=DT::styleEqual(color_vec, color_vec))
    
    
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

    inters <- vals$inters_filtered
    
    data_uploads <- c("anti_data_input","sempi_data_input","prism_data_input","prism_supp_data_input",
                      "arts_data_input","deep_data_input","gecco_data_input","rre_data_input")
    soft_names <- c("anti","sempi","prism","prism_supp","arts","deep","gecco","rre" )
    soft <- c("Antismash","SEMPI","PRISM","PRISM-Supp","ARTS","DeepBGC","GECCO","RRE-Finder" )
    data_to_use <- c( "anti_data" ,"sempi_data" , "prism_data", "prism_supp_data","arts_data","deep_data_filtered" ,"gecco_data_filtered", "rre_data")
    abbr <- c("A", "S", "P", "P-supp", "AR", "D", "G", "RRE")
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
        counts_var$label <- rep(soft[index], length(counts_var$x))
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
    
    refine_unique <- function(data){
      n <- tail(data, n=1)
      data <- head(data, -1)
      n_list <-  str_split(n, ",")
      out <- sapply(n_list[[1]], function(x){x %in% unlist(str_split(data, ","))})
      res <- sapply(out, function(x){
        if (x==F){
          x
        }
      })
      
      return(paste(names(Filter(Negate(is.null), res)), collapse = ","))
    }
    
    inters <- vals$inters_filtered 
    data_uploads <- c("anti_data_input","sempi_data_input","prism_data_input","prism_supp_data_input",
                      "arts_data_input","deep_data_input","gecco_data_input","rre_data_input")
    soft_names <- c("anti","sempi","prism","prism_supp","arts","deep","gecco","rre" )
    soft_let <- c("A", "S", "P", "PS", "AR", "D", "G", "R")
    soft_namings <- c('Antismash', 'SEMPI','PRISM', 'PRISM-Supp', 'ARTS', 'DeepBGC', 'GECCO', 'RRE_Finder')
    
    df_test <- data.frame(matrix(ncol = length(soft_let), nrow = 0))
    colnames(df_test) <- soft_let
    added_inters <- c(soft_names[match(input$group_by, soft_let)])
    add_inters <- list()
    if (input$count_all == F){
      df_test[nrow(df_test)+1,] <- NA
    } else{
      if( (input$group_by == "D") | (input$group_by == "G")){
        selected_dataframe <- paste0(soft_names[match(input$group_by, soft_let)], "_data_filtered")
      } else {
        selected_dataframe <- paste0(soft_names[match(input$group_by, soft_let)], "_data")
      }
      df_test <- data.frame(matrix(ncol = length(soft_let), nrow = length(vals[[selected_dataframe]]$Cluster)))
      colnames(df_test) <- soft_let
      df_test[[input$group_by]]<- vals[[selected_dataframe]]$Cluster
      df_test[nrow(df_test)+1,] <- NA
    }
    for (i in seq(1:length(data_uploads))){
      if (input$group_by==soft_let[[i]]){
        exclude <- i
        soft_n <- names(inters[[soft_names[i]]])
        index <- 1
        for (d in seq(1:length(soft_n))) {
          name <- soft_n[[index]]
          df_tmp <-  data.frame(cbind(c(inters[[soft_names[i]]][[soft_n[d]]]$to), c(inters[[soft_names[i]]][[soft_n[d]]]$from)))
          for (h in seq(1:length(soft_n))){
            if (name==soft_names[match(soft_n, soft_names)][h]){
              colnames(df_tmp) <- c(soft_let[i],soft_let[match(soft_n, soft_names)][h])
              df_test <- merge(df_test, df_tmp,  all = T)
            }
          }
          
          index <- index +1
        }
        excluded_names <- soft_let[soft_let != as.name(soft_let[i])]
        data <- df_test %>% dplyr::group_by_if(colnames(df_test)==soft_let[i]) %>% dplyr::summarise(a = paste(eval(as.name(excluded_names[1])), collapse=","),
                                                                                      b=paste(eval(as.name(excluded_names[2])), collapse=","),
                                                                                      c=paste(eval(as.name(excluded_names[3])), collapse=","),
                                                                                      d=paste(eval(as.name(excluded_names[4])), collapse=","),
                                                                                      e=paste(eval(as.name(excluded_names[5])), collapse=","),
                                                                                      f=paste(eval(as.name(excluded_names[6])), collapse=","),
                                                                                      g=paste(eval(as.name(excluded_names[7])), collapse=","))
        colnames(data) <- c(soft_let[i], excluded_names)
        for (p in soft_let){
          data[[p]] <- gsub('NA,|,NA', '', data[[p]])
          data[[p]][nrow(data)] <- refine_unique(data[[p]])
          names(data)[names(data) == p] <- soft_namings[match(p, soft_let)]
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
