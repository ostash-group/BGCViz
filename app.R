# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above. 
#
# Author: Pavlo Hrab
# Made as part of Cambridge bioinformatics hackaton
# 
# This app is using bgc coordinates from DeepBGC, PRISM, ANTISMASH and RREFinder
# to visualized interception of those different annotations in one genome
#

# Upload required libraries
library(shiny)
library(tidyverse)
library(plyr)
library(plotly)
library(BioCircos)
library(ggplot2)
library(shinyjs)
library(rjson)
library(stringr)
library(DT)
library(GenomicRanges)

# Define UI 
ui <- fluidPage(
  
  # Application title
  titlePanel("BGCViz"),
  
  # Sidebar
  useShinyjs(),
  sidebarLayout(
    sidebarPanel(
      # Data upload
      h3("Data upload and necesary input:"),
      checkboxInput("hide_uploads", "Hide upload fields"),
      h5(id = "anti_header_upload","ANTISMASH:"),
      checkboxInput("anti_input_options", "My AntiSMASH data is a dataframe, not json results file from antismash", value = T),
      fileInput("anti_data",
                "Upload Antismash data"),
      actionButton("anti_sco", "Use Antismash example data from S.coelicolor"),
      h5(id = "prism_header_upload","PRISM:"),
      checkboxInput("prism_input_options", "My PRISM data is a dataframe, not json results file", value = T),
      fileInput("prism_data",
                "Upload PRISM data"),
      actionButton("prism_sco", "Use PRISM example data from S.coelicolor"),
      h5(id = "sempi_header_upload","SEMPI:"),
      fileInput("sempi_data",
                "Upload SEMPI 2.0 data"),
      actionButton("sempi_sco", "Use SEMPI example data from S.coelicolor"),
      h5(id = "deep_header_upload","DEEPBGC:"),
      fileInput("deep_data",
                "Upload DeepBGC data"),
      actionButton("deep_sco", "Use DeepBGC example data from S.coelicolor"),
      h5(id = "gecco_header_upload","GECCO:"),
      fileInput("gecco_data",
                "Upload Gecco data"),
      actionButton("gecco_sco", "Use Gecco example data from S.coelicolor"),
      h5(id = "rre_header_upload","RRE-FINDER:"),
      fileInput("rre_data",
                "Upload RRE-Finder data"),
      actionButton("rre_sco", "Use RRE-Finder example data from S.coelicolor"),
      h5(id = "arts_header_upload","ARTS:"),
      fileInput("known_data",
                "Upload ARTS knownhits data"),
      fileInput("dup_data",
                "Upload ARTS duptable data"),
      actionButton("arts_sco", "Use ARTS example data from S.coelicolor"),
      # Numeric input of chromosome length of analyzed sequence
      numericInput("chr_len", "Please type chr len of an organism", value = 10000000),
      h3("Data manipulation options"),
      checkboxInput("hide_anti", "Hide data manipulation fields"),
      h5(id = "anti_header","Antismash data options:"),
      checkboxInput("anti_hybrid", "Visualize AntiSMASH BGC with several types as 'Hybrid'"),
      h5(id = "prism_header","PRISM data options:"),
      checkboxInput("prism_hybrid", "Visualize PRISM BGC with several types as 'Hybrid'"),
      checkboxInput("prism_supp", "Use PRISM resistance and regulatory genes information'"),
      h5(id = "sempi_header","SEMPI data options:"),
      checkboxInput("sempi_hybrid", "Visualize SEMPI BGC with several types as 'Hybrid'"),
      h5(id = "arts_header","ARTS data options:"),
      selectInput("dup_choice", "Choose duplicated core gene to plot only it", choices = c("All"),
                  selected = "All"),
      h3(id = "genes_on_chr","Genes on chromosome plot controls:"),
      checkboxInput("hide_genes_on_chr", "Hide 'Genes on chromosome plot' fields"),
      selectInput("ref", "Choose reference data", choices = c("Antismash" = "Antismash",
                                                              "DeepBGC" = "DeepBGC",
                                                              "RRE-Finder" = "RRE-Finder",
                                                              "PRISM" = "PRISM",
                                                              "SEMPI" = "SEMPI",
                                                              "PRISM-supp" = "PRISM-supp",
                                                              "ARTS" = "ARTS",
                                                              "GECCO" = "GECCO"),
                  selected = "Antismash"),
      h3(id = "summarize","Summarize options:"),
      checkboxInput("hide_summarize", "Hide summarize options"),
      selectInput("group_by", "Group data by", choices = c("Antismash" = "A",
                                                              "DeepBGC" = "D",
                                                              "RRE-Finder" =  "R",
                                                              "PRISM" = "P",
                                                           "SEMPI" = "S",
                                                           "PRISM-supp" = "PS",
                                                           "ARTS" = "AR",
                                                           "GECCO" = "G"),
                  selected = 'A'),
      checkboxInput("count_all", "Show all BGC for the 'group by' method (+ individually annotated BGC)"),
      h3("Improve visualization:"),
      checkboxInput("hide_viz", "Hide improve visualization options"),
      #Improve RREFinder annotated BCG visibility
      fileInput("rename_data",
               "Upload renaming and coloring scheme"),
      actionButton("rename", "Rename"),
      actionButton("reset_name", "Reset"),
      checkboxInput("rre_width", "Add thickness to RRE results visualization"),
      checkboxInput("prism_supp_width", "Add thickness to PRISM resistance + regulatory genes results visualization"),
      checkboxInput("arts_width", "Add thickness to ARTS results visualization"),
      checkboxInput("biocircos_color", "Make arcs in biocircos colorful, based on the class"),
      checkboxInput("label_color", "Make links in biocircos colorful, based on the class"),
      selectInput("label_color_class", "Choose the mode to color the links", choices = c("Hierarchical-based" = "H",
                                                                                           "Purity-based" = "P",
                                                                                           "Reference column-based" = "R"
                                                                                         ),
                  selected = 'H'),
      selectInput("ref_col_biocircos", "Choose reference column to color the links", choices = c("Antismash" = "Antismash",
                                                                                         "PRISM" = "PRISM",
                                                                                         "RRE-Finder" = "RRE",
                                                                                         'DeepBGC' = "DeepBGC",
                                                                                         "SEMPI" = "SEMPI",
                                                                                         "PRISM-supp" = "PRISM-supp",
                                                                                         "ARTS" = "ARTS",
                                                                                         "GECCO" = "GECCO"
      ),
      selected = 'Antismash'),
      h3(id="data_comparison_header_gecco","Comparison with Gecco plots:"),
      checkboxInput("hide_data_comparison_gecco", "Hide Gecco data comparison options"),
      selectInput("ref_comparison_gecco", "Choose data for comparison with Gecco", choices = c("Antismash" = "A",
                                                                                           "PRISM" = "P",
                                                                                           "SEMPI" = "S"),
                  selected = 'A'),
      selectInput("score_type_gecco", "Choose score type to set threshold", choices = c("Average p-value" = "avg_p",
                                                                                  "Cluster_type score" = "Cluster_Type"),
                  selected = "avg_p"),
      # Chose step for barplot (as a threshold to draw a bar)
      sliderInput("plot_step_gecco", "Choose step for plots(barplot)", min = 1, max = 50,value = 10),
      sliderInput("plot_start_gecco", "Chose plot start point(barplot)", min = 0, max = 100, value = 0),
      h3(id="data_filter_header_gecco","Gecco data filtering:"),
      checkboxInput("hide_data_filter_gecco", "Hide Gecco data filtering options"),
      sliderInput("score_average_gecco", "Average p-value threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50 ),
      sliderInput("score_cluster_gecco", "Cluster type threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50 ),
      sliderInput("domains_filter_gecco", "Domain number threshold for Gecco data", min = 0, max = 100, value = 1),
      sliderInput("prot_filter_gecco", "Protein number threshold for Gecco data", min = 0, max = 100, value = 1),
      h3(id="data_comparison_header","Comparison with DeepBGC plots:"),
      checkboxInput("hide_data_comparison", "Hide DeepBGC data comparison options"),
      selectInput("ref_comparison", "Choose data for comparison with DeepBGC", choices = c("Antismash" = "A",
                                                                     "PRISM" = "P",
                                                                     "SEMPI" = "S"),
                  selected = 'A'),
      # Score to use for thresholds
      selectInput("score_type", "Choose score type to set threshold", choices = c("Activity score" = "Activity",
                                                                                  "Cluster_type score" = "Cluster_Type",
                                                                                  "DeepBGC score" = "DeepBGC"),
                  selected = "Activity score"),
      # Chose step for barplot (as a threshold to draw a bar)
      sliderInput("plot_step", "Choose step for plots(barplot)", min = 1, max = 50,value = 10),
      sliderInput("plot_start", "Chose plot start point(barplot)", min = 0, max = 100, value = 0),
      
      # DeepBGC data filtering 
      h3(id="data_filter_header","DeepBGC data filtering:"),
      checkboxInput("hide_data_filter", "Hide DeepBGC data filtering options"),
      # Different score filtering. Remain >= of set threshold
      sliderInput("score_a", "Activity score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      sliderInput("score_d", "DeepBGC score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      sliderInput("score_c", "Cluster_type score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      # Domains, biodomains and proteins filter. Remain >= of set threshold
      sliderInput("domains_filter", "Domain number threshold for DeepBGC data", min = 0, max = 100, value = 5),
      sliderInput("biodomain_filter", "Biodomain number threshold for DeepBGC data", min = 0, max = 100, value = 1),
      sliderInput("gene_filter", "Protein number threshold for DeepBGC data", min = 0, max = 100, value = 1),
      sliderInput("cluster_type","Choose threshold to assign cluster type for DeepBGC data ", min = 0, max = 100, value = 50),
     
      # Donwload currently used datasets
      downloadButton("download","Download currently used datasets (as for Biocircos plot)" )
      
    ),
    
    # Show plots
    mainPanel(
      tabsetPanel(
        tabPanel(title = "Compare data with DeepBGC", value = 1 ,plotOutput("deep_barplot",height = "500px"), plotlyOutput("deep_rate")),
        tabPanel(title = "Compare data with Gecco", value = 5 ,plotOutput("gecco_barplot",height = "500px"), plotlyOutput("gecco_rate")),
        tabPanel(title = "Annotation visualization and comparison", value = 4,plotlyOutput("deep_reference_2", height = "500px"), 
                 plotlyOutput("deep_reference", height = "500px")),
        tabPanel(title = "Biocircos plot", value = 2, BioCircosOutput("biocircos", height = "1000px"), dataTableOutput("biocircos_legend")),
        tabPanel(title = "Summarize interception", value = 3,plotlyOutput("barplot_rank", height = "600px"),tableOutput("group_table")),
        type = "tabs", id = "main"
      )
  )
  )
)

# Define server logic
server <- function(input, output, session) {
  #
  options(shiny.maxRequestSize=100*1024^2)
  # Small function to make integers zeros
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }
  
  biocircos_listen <- reactive({
    list( input$biocircos_color, vals$arts_data_filtered, input$label_color, input$label_color_class, 
          input$ref_col_biocircos, vals$deep_data_filtered, vals$gecco_data_filtered, vals$inters_filtered
          )
  })
  inputData <- reactive({
    list( vals$sempi_data,vals$rre_data,  vals$anti_data, vals$prism_data,
          vals$arts_data, vals$prism_supp, vals$deep_data, vals$gecco_data
    )
  })
  dynamicInput <-  reactive({
    list(  input$cluster_type, input$gene_filter,input$biodomain_filter,  input$score_c, input$score_d, 
          input$score_a,  input$score_average_gecco,input$score_cluster_gecco, input$domains_filter_gecco, 
          input$prot_filter_gecco, input$dup_choice, vals$need_filter
    )
  })
  


  # Rective vals the app is using
  # Some dataframes that are used through the app + some vectors of untercepted values
  vals <- reactiveValues(deep_data = NULL, anti_data = NULL, rre_data=NULL, prism_data=NULL, chr_len = NULL, fullness = NULL,
                         biocircos_deep = NULL, deep_data_input = FALSE,tracklist = NULL, chromosomes = NULL, 
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
                         need_filter = F
                         )
  
  vals$computed <- list(
    anti=F,deep=F, gecco=F, arts=F, prism=F, sempi=F, prism_supp=F, rre=F
  )
  
  vals$rename_data <- read.csv("rename.csv")
  
  # Upload example data
  observeEvent(input$anti_sco,{
    whereami::cat_where(whereami::whereami())
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
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "Antismash" )
      updateSelectInput(session, "group_by",
                        selected = "A" )
      updateSelectInput(session, "ref_comparison",
                        selected = "A")
      updateSelectInput(session, "ref_col_biocircos",
                        selected =  "Antismash")
      updateSelectInput(session, "ref_comparison_gecco",
                        selected = "A")

    }
    
  })
  
  observeEvent(input$gecco_sco,{
    whereami::cat_where(whereami::whereami())
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
      mutate(pks=polyketide_probability, other = other_probability, nrps = nrp_probability, alkaloid = alkaloid_probability, 
             terpene = terpene_probability, saccharide = saccharide_probability, ripp = ripp_probability) %>%
      select(-one_of(drop_cols))
    gecco_data$num_prot <- sapply( str_split(as.character(gecco_data$proteins), ";"), length)
    gecco_data$num_domains <- sapply( str_split(as.character(gecco_data$domains), ";"), length)
    names(gecco_data)[names(gecco_data) == "start"] <- "Start"
    names(gecco_data)[names(gecco_data) == "end"] <-  "Stop"
    vals$gecco_data <- gecco_data
    # Save file
    write.csv(vals$gecco_data, "gecco_data.csv", row.names = F)
    vals$gecco_data_input = TRUE 
    vals$data_upload_count <- vals$data_upload_count +1
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "GECCO" )
      updateSelectInput(session, "group_by",
                        selected = "G")
      updateSelectInput(session, "ref_col_biocircos",
                        selected =  "GECCO")
      
    }
    
  })
  
  observeEvent(input$prism_sco,{
    # Read data
      data <- fromJSON(file = "example_data/sco_prism.json")
      
      
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
      vals$biocircos_color = F
      
      regul_genes_orfs <- sapply(data$prism_results$regulatory_genes, function(x){
        x$orf
      })
      
      location <-  sapply(data$prism_results$orfs[[1]]$orfs, function(y){
        sapply(regul_genes_orfs, function(x){
          if (y$name == x) {
            y$coordinates
          }
        })
      }) 
      
      location <- Filter(Negate(is.null), location)
      
      reg_genes <-data.frame(t(data.frame(sapply(location, function(x){unlist(x)}))))
      colnames(reg_genes) <- c("Start", "Stop")
      reg_genes$Type <- 'regulatory'
      reg_genes$Type2 <- reg_genes$Type
      reg_genes$Score <- sapply(data$prism_results$regulatory_genes, function(x){
        x$score
      })
      reg_genes$Name <- sapply(data$prism_results$regulatory_genes, function(x){
        x$name
      })
      reg_genes$Full_name <- sapply(data$prism_results$regulatory_genes, function(x){
        x$full_name
      })
      
      resist_genes_orfs <- sapply(data$prism_results$resistance_genes, function(x){
        x$orf
      })
      
      location <-  sapply(data$prism_results$orfs[[1]]$orfs, function(y){
        sapply(resist_genes_orfs, function(x){
          if (y$name == x) {
            y$coordinates
          }
        })
      })
      
      location <- Filter(Negate(is.null), location)
      
      res_genes <-data.frame(t(data.frame(sapply(location, function(x){unlist(x)}))))
      colnames(res_genes) <- c("Start", "Stop")
      res_genes$Type <- 'resistance'
      res_genes$Type2 <- res_genes$Type
      res_genes$Score <- sapply(data$prism_results$resistance_genes, function(x){
        x$score
      })
      res_genes$Name <- sapply(data$prism_results$resistance_genes, function(x){
        x$name
      })
      res_genes$Full_name <- sapply(data$prism_results$resistance_genes, function(x){
        x$full_name
      })
      
      final_reg <- rbind(res_genes, reg_genes)
      final_reg$ID <- seq(1:dim(final_reg)[1])
      final_reg$Cluster <- final_reg$ID
      rownames(final_reg) <- as.numeric(seq(1:dim(final_reg)[1]))
      vals$prism_supp <- final_reg
      vals$prism_supp_data <- final_reg
      vals$prism_json = T
      vals$prism_supp_data_input=T
      
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
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "PRISM" )
      updateSelectInput(session, "group_by",
                        selected = "P" )
      updateSelectInput(session, "ref_comparison",
                        selected = "P")
      updateSelectInput(session, "ref_col_biocircos",
                        selected =  "PRISM")
      updateSelectInput(session, "ref_comparison_gecco",
                        selected = "P")
    }
  })
  
  observeEvent(input$sempi_sco,{
    whereami::cat_where(whereami::whereami())
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
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "SEMPI" )
      updateSelectInput(session, "group_by",
                        selected = "S" )
      updateSelectInput(session, "ref_comparison",
                        selected = "S")
      updateSelectInput(session, "ref_col_biocircos",
                        selected =  "SEMPI")
      updateSelectInput(session, "ref_comparison_gecco",
                        selected = "S")
    }
    
  })
  
  observeEvent(input$arts_sco, {
    
    data <- read.delim("example_data/sco_duptable.tsv")
    
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
      arts_data$ID <- seq(1:dim(arts_data)[1])
      arts_data$Cluster <- arts_data$ID
      vals$arts_data <- arts_data
      vals$data_upload_count <-  vals$data_upload_count +1
      vals$arts_data_input <- T
      dup_table_id <- arts_data %>%
        filter(Core != "Not_core")
      updateSelectInput(session, "dup_choice",
                        choices = c("All", paste0("ID:",dup_table_id$ID, " ,Core:", dup_table_id$Core)),
                        selected = "All" )
      vals$upl_arts = T
      if (vals$data_upload_count == 1){
        updateSelectInput(session, "ref",
                          selected = "ARTS" )
        updateSelectInput(session, "group_by",
                          selected = "AR" )
        updateSelectInput(session, "ref_col_biocircos",
                          selected =  "ARTS")
      }
  })
  
  observeEvent(input$deep_sco, {
    drop_cols <- c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")
    # Read data
    vals$deep_data <- read.delim("example_data/sco_deep.tsv") %>%
      mutate(pks=Polyketide, other = Other, nrps = NRP, alkaloid = Alkaloid, 
             terpene = Terpene, saccharide = Saccharide, ripp = RiPP) %>%
      select(-one_of(drop_cols))
    # Add chromosome info column
    vals$deep_data$chromosome <-  rep("D", length(vals$deep_data$bgc_candidate_id))
    # Add ID column as number seuquence of dataframe length
    vals$deep_data$ID <- seq(1:length(vals$deep_data$bgc_candidate_id))
    vals$deep_data$Cluster <- vals$deep_data$ID
    write.csv(vals$deep_data, "deep_data.csv", row.names = F)
    vals$deep_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "DeepBGC" )
      updateSelectInput(session, "group_by",
                        selected = "D" )
      updateSelectInput(session, "ref_col_biocircos",
                        choices = "DeepBGC",
                        selected = "DeepBGC")
      
    }
  })
  
  observeEvent(input$rre_sco, {
    whereami::cat_where(whereami::whereami())
    # Read data
    vals$rre_data <- read.delim("example_data/sco_rre.txt")
    # Clean RRE data. Extract coordinates and Locus tag with double underscore delimiter (__)
    vals$rre_data <- vals$rre_data %>%
      separate(Gene.name, c("Sequence","Coordinates","Locus_tag"),sep = "__") %>%
      separate(Coordinates, c("Start", "Stop"),sep = "-")
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
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "RRE-Finder" )
      updateSelectInput(session, "group_by",
                        selected = "R" )
      updateSelectInput(session, "ref_col_biocircos",
                        choices = "RRE-Finder",
                        selected = "RRE")
      
    }
    if (!is.null(vals$rre_data$Probability)){
      vals$rre_more = T
    } else {
      vals$rre_more = F
    }
  })
  
  # Read the data
  observeEvent(input$anti_data,{
    # Read data
    if (input$anti_input_options==T){
      anti_data <- read.csv(input$anti_data$datapath)
      vals$biocircos_color = T
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
          separate(split, c("Start", "Stop")) %>%
          transmute(ID = rownames(location), Start, Stop)
        
        anti_data <- cbind(anti_data, types)
        colnames(anti_data) <- c("Cluster", "Start", "Stop", "Type")
        anti_data$Cluster <- as.numeric(anti_data$Cluster)
        anti_data$Start <- as.numeric(anti_data$Start)
        anti_data$Stop <- as.numeric(anti_data$Stop)
        vals$biocircos_color = F
        
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
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "Antismash" )
      updateSelectInput(session, "group_by",
                        selected = "A" )
      updateSelectInput(session, "ref_comparison",
                        selected = "A")
      updateSelectInput(session, "ref_col_biocircos",
                        selected =  "Antismash")
      updateSelectInput(session, "ref_comparison_gecco",
                        selected = "A")
    }

  })
  
  observeEvent(input$sempi_data,{
    
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
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "SEMPI" )
      updateSelectInput(session, "group_by",
                        selected = "S" )
      updateSelectInput(session, "ref_comparison",
                        selected = "S")
      updateSelectInput(session, "ref_col_biocircos",
                        selected =  "SEMPI")
      updateSelectInput(session, "ref_comparison_gecco",
                        selected = "S")
    }
    
  })
  
  observeEvent(input$gecco_data,{
    
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
      mutate(pks=polyketide_probability, other = other_probability, nrps = nrp_probability, alkaloid = alkaloid_probability, 
             terpene = terpene_probability, saccharide = saccharide_probability, ripp = ripp_probability) %>%
      select(-one_of(drop_cols))
    gecco_data$num_prot <- sapply( str_split(as.character(gecco_data$proteins), ";"), length)
    gecco_data$num_domains <- sapply( str_split(as.character(gecco_data$domains), ";"), length)
    names(gecco_data)[names(gecco_data) == "start"] <- "Start"
    names(gecco_data)[names(gecco_data) == "end"] <-  "Stop"
    vals$gecco_data <- gecco_data
    # Save file
    write.csv(vals$gecco_data, "gecco_data.csv", row.names = F)
    vals$gecco_data_input = TRUE 
    vals$data_upload_count <- vals$data_upload_count +1
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "GECCO" )
      updateSelectInput(session, "group_by",
                        selected = "G")
      updateSelectInput(session, "ref_col_biocircos",
                        selected =  "GECCO")
    }
    
  })
  
  observeEvent(input$known_data, {
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
      dup_table <- vals$dup_data
      known_table <- vals$known_data
      arts_data <- rbind(dup_table, known_table)
      arts_data$ID <- seq(1:dim(arts_data)[1])
      arts_data$Cluster <- arts_data$ID
      vals$arts_data <- arts_data
      vals$data_upload_count <-  vals$data_upload_count +1
      dup_table_id <- arts_data %>%
        filter(Core != "Not_core")
      updateSelectInput(session, "dup_choice",
                        choices = c("All", paste0("ID:",dup_table_id$ID, " ,Core:", dup_table_id$Core)),
                        selected = "All" )
      vals$upl_arts = T
      if (vals$data_upload_count == 1){
        updateSelectInput(session, "ref",
                          selected = "ARTS" )
        updateSelectInput(session, "group_by",
                          selected = "AR" )
        updateSelectInput(session, "ref_col_biocircos",
                          selected =  "ARTS")
      }
    } 
  })
  
  observeEvent(input$dup_data, {
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
      dup_table <- vals$dup_data
      known_table <- vals$known_data
      arts_data <- rbind(dup_table, known_table)
      arts_data$ID <- seq(1:dim(arts_data)[1])
      arts_data$Cluster <- arts_data$ID
      vals$arts_data <- arts_data
      vals$data_upload_count <-  vals$data_upload_count +1
      vals$arts_data_input <- T
      dup_table_id <- arts_data %>%
        filter(Core != "Not_core")
      updateSelectInput(session, "dup_choice",
                        choices = c("All", paste0("ID:",dup_table_id$ID, " ,Core:", dup_table_id$Core)),
                        selected = "All" )
      if (vals$data_upload_count == 1){
        updateSelectInput(session, "ref",
                          selected = "ARTS" )
        updateSelectInput(session, "group_by",
                          selected = "AR" )
        updateSelectInput(session, "ref_col_biocircos",
                          selected =  "ARTS")
      }
    } 
  })
  
  observeEvent(input$prism_data,{
    # Read data
    if (input$prism_input_options == T){
      prism_data <- read.csv(input$prism_data$datapath)
      vals$biocircos_color = T
    } else{
      data <- fromJSON(file = input$prism_data$datapath)
      
      
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
      vals$biocircos_color = F
      
      regul_genes_orfs <- sapply(data$prism_results$regulatory_genes, function(x){
        x$orf
      })
      
      location <-  sapply(data$prism_results$orfs[[1]]$orfs, function(y){
        sapply(regul_genes_orfs, function(x){
          if (y$name == x) {
            y$coordinates
          }
        })
      }) 
      
      location <- Filter(Negate(is.null), location)
      
      reg_genes <-data.frame(t(data.frame(sapply(location, function(x){unlist(x)}))))
      colnames(reg_genes) <- c("Start", "Stop")
      reg_genes$Type <- 'regulatory'
      reg_genes$Type2 <- reg_genes$Type
      reg_genes$Score <- sapply(data$prism_results$regulatory_genes, function(x){
        x$score
      })
      reg_genes$Name <- sapply(data$prism_results$regulatory_genes, function(x){
        x$name
      })
      reg_genes$Full_name <- sapply(data$prism_results$regulatory_genes, function(x){
        x$full_name
      })
      
      resist_genes_orfs <- sapply(data$prism_results$resistance_genes, function(x){
        x$orf
      })
      
      location <-  sapply(data$prism_results$orfs[[1]]$orfs, function(y){
        sapply(resist_genes_orfs, function(x){
          if (y$name == x) {
            y$coordinates
          }
        })
      })
      
      location <- Filter(Negate(is.null), location)
      
      res_genes <-data.frame(t(data.frame(sapply(location, function(x){unlist(x)}))))
      colnames(res_genes) <- c("Start", "Stop")
      res_genes$Type <- 'resistance'
      res_genes$Type2 <- res_genes$Type
      res_genes$Score <- sapply(data$prism_results$resistance_genes, function(x){
        x$score
      })
      res_genes$Name <- sapply(data$prism_results$resistance_genes, function(x){
        x$name
      })
      res_genes$Full_name <- sapply(data$prism_results$resistance_genes, function(x){
        x$full_name
      })
      
      final_reg <- rbind(res_genes, reg_genes)
      final_reg$ID <- seq(1:dim(final_reg)[1])
      final_reg$Cluster <- final_reg$ID
      rownames(final_reg) <- as.numeric(seq(1:dim(final_reg)[1]))
      vals$prism_supp <- final_reg
      vals$prism_supp_data <- final_reg
      vals$prism_json = T
      vals$prism_supp_data_input = T
      
    }
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
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "PRISM" )
      updateSelectInput(session, "group_by",
                        selected = "P" )
      updateSelectInput(session, "ref_comparison",
                        selected = "P")
      updateSelectInput(session, "ref_col_biocircos",
                        selected =  "PRISM")
      updateSelectInput(session, "ref_comparison_gecco",
                        selected = "P")
    }
  })
  
  observeEvent(input$deep_data, {
    drop_cols <- c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")
    # Read data
    vals$deep_data <- read.delim(input$deep_data$datapath) %>%
      mutate(pks=Polyketide, other = Other, nrps = NRP, alkaloid = Alkaloid, 
             terpene = Terpene, saccharide = Saccharide, ripp = RiPP) %>%
      select(-one_of(drop_cols))
    # Add chromosome info column
    vals$deep_data$chromosome <-  rep("D", length(vals$deep_data$bgc_candidate_id))
    # Add ID column as number seuquence of dataframe length
    vals$deep_data$ID <- seq(1:length(vals$deep_data$bgc_candidate_id))
    vals$deep_data$Cluster <- vals$deep_data$ID
    write.csv(vals$deep_data, "deep_data.csv", row.names = F)
    vals$deep_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "DeepBGC" )
      updateSelectInput(session, "group_by",
                        selected = "D" )
      updateSelectInput(session, "ref_col_biocircos",
                        choices = "DeepBGC",
                        selected = "DeepBGC")
      
    }
  })
  
  observeEvent(input$rre_data, {
    # Read data
    vals$rre_data <- read.delim(input$rre_data$datapath)
    # Clean RRE data. Extract coordinates and Locus tag with double underscore delimiter (__)
    vals$rre_data <- vals$rre_data %>%
      separate(Gene.name, c("Sequence","Coordinates","Locus_tag"),sep = "__") %>%
      separate(Coordinates, c("Start", "Stop"),sep = "-")
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
    if (vals$data_upload_count == 1){
      updateSelectInput(session, "ref",
                        selected = "RRE-Finder" )
      updateSelectInput(session, "group_by",
                        selected = "R" )
      updateSelectInput(session, "ref_col_biocircos",
                        choices = "RRE-Finder",
                        selected = "RRE")
      
    }
    if (!is.null(vals$rre_data$Probability)){
      vals$rre_more = T
    } else {
      vals$rre_more = F
    }
  })
  
  # Observe input of chromosome length
  observeEvent(input$chr_len,{
    vals$chr_len <- input$chr_len
  })
  
  # Logic for showing/hiding UI when input 
  observeEvent(vals$rre_data_input, {
    whereami::cat_where(whereami::whereami())
    if (vals$rre_data_input == T){
      if (input$hide_viz == F){
        showElement(selector = "#rre_width")
      }
    } else{
      hideElement(selector = "#rre_width")
    }
  })

  observeEvent(vals$anti_data_input, {
    whereami::cat_where(whereami::whereami())
    if (vals$anti_data_input == T){
      if (input$hide_anti == F){
        showElement(selector = "#anti_header")
        showElement(selector = "#anti_hybrid")
      }
    } else{
      hideElement(selector = "#anti_header")
      hideElement(selector = "#anti_hybrid")
    }
  })
  
  observeEvent(vals$prism_data_input, {
    if (vals$prism_data_input == T){
      if (input$hide_anti == F){
        showElement(selector = "#prism_header")
        showElement(selector = "#prism_hybrid")
        if (vals$prism_json == T){
          showElement(selector = "#prism_supp")
        }
      }
      if (input$hide_viz == F){
        if (vals$prism_json == T){
          showElement(selector = "#prism_supp_width")
        }
      }
    } else{
      hideElement(selector = "#prism_header")
      hideElement(selector = "#prism_hybrid")
      hideElement(selector = "#prism_supp")
      hideElement(selector = "#prism_supp_width")
    }
  })
  
  observeEvent(vals$sempi_data_input, {
    whereami::cat_where(whereami::whereami())
    if (vals$sempi_data_input == T){
      if (input$hide_anti == F){
        showElement(selector = "#sempi_header")
        showElement(selector = "#sempi_hybrid")
      }
    } else{
      hideElement(selector = "#sempi_header")
      hideElement(selector = "#sempi_hybrid")
    }
  })

  observeEvent(vals$deep_data_input,{
    if (vals$deep_data_input == T){
      showElement(selector = "#ref_comparison")
      showElement(selector = "#hide_data_comparison")
      showElement(selector = "#hide_data_filter")
      showElement(selector = "#score_type")
      showElement(selector = "#plot_step")
      showElement(selector = "#plot_start")
      showElement(selector = "#score_a")
      showElement(selector = "#score_d")
      showElement(selector = "#score_c")
      showElement(selector = "#domains_filter")
      showElement(selector = "#biodomain_filter")
      showElement(selector = "#gene_filter")
      showElement(selector = "#cluster_type")
      showElement(selector = "#data_comparison_header")
      showElement(selector = "#data_filter_header")
    } else{
      hideElement(selector = "#ref_comparison")
      hideElement(selector = "#score_type")
      hideElement(selector = "#hide_data_comparison")
      hideElement(selector = "#hide_data_filter")
      hideElement(selector = "#plot_step")
      hideElement(selector = "#plot_start")
      hideElement(selector = "#score_a")
      hideElement(selector = "#score_d")
      hideElement(selector = "#score_c")
      hideElement(selector = "#domains_filter")
      hideElement(selector = "#biodomain_filter")
      hideElement(selector = "#gene_filter")
      hideElement(selector = "#cluster_type")
      hideElement(selector = "#data_comparison_header")
      hideElement(selector = "#data_filter_header")
    }
  })
  
  observeEvent(vals$gecco_data_input,{
    if (vals$gecco_data_input == T){
      showElement(selector = "#data_comparison_header_gecco")
      showElement(selector = "#hide_data_comparison_gecco")
      showElement(selector = "#ref_comparison_gecco")
      showElement(selector = "#score_type_gecco")
      showElement(selector = "#plot_step_gecco")
      showElement(selector = "#plot_start_gecco")
      showElement(selector = "#data_filter_header_gecco")
      showElement(selector = "#hide_data_filter_gecco")
      showElement(selector = "#score_average_gecco")
      showElement(selector = "#score_cluster_gecco")
      showElement(selector = "#domains_filter_gecco")
      showElement(selector = "#prot_filter_gecco")
    } else{
      hideElement(selector = "#data_comparison_header_gecco")
      hideElement(selector = "#hide_data_comparison_gecco")
      hideElement(selector = "#ref_comparison_gecco")
      hideElement(selector = "#score_type_gecco")
      hideElement(selector = "#plot_step_gecco")
      hideElement(selector = "#plot_start_gecco")
      hideElement(selector = "#data_filter_header_gecco")
      hideElement(selector = "#hide_data_filter_gecco")
      hideElement(selector = "#score_average_gecco")
      hideElement(selector = "#score_cluster_gecco")
      hideElement(selector = "#domains_filter_gecco")
      hideElement(selector = "#prot_filter_gecco")
    }
  })
  
  observeEvent(vals$data_upload_count, {
    whereami::cat_where(whereami::whereami())
    if (vals$data_upload_count <2){
      hideTab("main", "2")
      hideTab("main", "3")
      
    }else if (vals$data_upload_count >=2){
      if (input$hide_summarize == F) {
        showElement(selector = "#summarize")
        showElement(selector = "#group_by")
        showElement(selector = "#count_all")
      }
      if (input$hide_viz == F){
        showElement(selector = "#biocircos_color")
        showElement(selector = "#label_color")
        showElement(selector = "#label_color_class")
      }
      showTab("main", "2")
      showTab("main", "3")
      if (vals$gecco_data_input == T) {
        showTab("main", "5")
      } else {
        hideTab("main", "5")
      }
      if (vals$deep_data_input == T) {
        showTab("main", "1")
      } else {
        hideTab("main", "1")
      }
    }
    if (vals$data_upload_count <1){
      hideTab("main", "4")
      hideTab(inputId = "main", target = "5")
      hideTab(inputId = "main", target = "1")
      hideElement(selector = "#genes_on_chr")
      hideElement(selector = "#hide_genes_on_chr")
      hideElement(selector = "#ref")
    }else{
      showTab("main", "4")
      if (input$hide_genes_on_chr == F){
        showElement(selector = "#genes_on_chr")
        showElement(selector = "#hide_genes_on_chr")
        showElement(selector = "#ref")
      }
      
    }
  })
  
  observeEvent(input$label_color_class, {
    if (input$label_color_class == "R"){
      showElement(selector = "#ref_col_biocircos")
    } else {
      hideElement(selector = "#ref_col_biocircos")
    }
  })
  
  observeEvent(input$anti_hybrid, {
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

  observeEvent(input$prism_hybrid, {
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
  
  observeEvent(input$sempi_hybrid, {
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
    
    if (input$sempi_hybrid==T){
      vals$sempi_data$Type2 <- hybrid_col(vals$sempi_data)
    }else {
      vals$sempi_data$Type2 <- vals$sempi_type
    }
  })
  
  observeEvent(vals$arts_data_input,{
    if (vals$arts_data_input == T){
      if (input$hide_anti == F){
        showElement(selector = "#arts_header")
        showElement(selector = "#dup_choice")
      }
      if (input$hide_viz == F){
        showElement(selector = "#arts_width")
      }
    } else {
      hideElement(selector = "#arts_header")
      hideElement(selector = "#dup_choice")
      hideElement(selector = "#arts_width")
    }
  })
  
  observeEvent(input$rename, {
    rename_vector <- function(data, renamed_dataframe){
      type <- str_split(data$Type, "__")
      type_2 <- sapply(type, function(x){
        sapply(x, function(y){
          if (y %in% renamed_dataframe$Code){
            renamed <- as.character(renamed_dataframe$Group[renamed_dataframe$Code == y])
            if (length(renamed) >1){
              showNotification(paste("The ", as.character(y), " type have multiple renaming options: ", paste(renamed, collapse = ", ")), 
                               type = "warning", duration = NULL)
              showNotification(paste("The  ", renamed[[1]], " was chosen."), type = "warning", duration=NULL)
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
    vals$biocircos_color = T
      })
  
  observeEvent(input$reset_name, {

    vals$anti_data['Type2']  <- vals$anti_data['Type']
    vals$sempi_data['Type2'] <- vals$sempi_data['Type']
    vals$ prism_data['Type2'] <- vals$ prism_data['Type']
    vals$biocircos_color = FALSE
  })
  
  observeEvent(input$rename_data,{
    rename_data <- read.csv(input$rename_data$datapath)
    vals$rename_data <- rename_data
  })
  
  observeEvent(input$hide_uploads, {
    if (input$hide_uploads == T){
      hideElement(selector = "#anti_input_options")
      hideElement(selector = "#anti_data")
      hideElement(selector = "#prism_input_options")
      hideElement(selector = "#anti_header_upload")
      hideElement(selector = "#prism_header_upload")
      hideElement(selector = "#prism_data")
      hideElement(selector = "#sempi_header_upload")
      hideElement(selector = "#sempi_data")
      hideElement(selector = "#deep_header_upload")
      hideElement(selector = "#deep_data")
      hideElement(selector = "#gecco_header_upload")
      hideElement(selector = "#gecco_data")
      hideElement(selector = "#rre_header_upload")
      hideElement(selector = "#rre_data")
      hideElement(selector = "#chr_len")
      hideElement(selector = "#arts_header_upload")
      hideElement(selector = "#known_data")
      hideElement(selector = "#dup_data")
      hideElement(selector = "#anti_sco")
      hideElement(selector = "#prism_sco")
      hideElement(selector = "#arts_sco")
      hideElement(selector = "#rre_sco")
      hideElement(selector = "#sempi_sco")
      hideElement(selector = "#deep_sco")
      hideElement(selector = "#gecco_sco")
    }else {
      showElement(selector = "#anti_input_options")
      showElement(selector = "#anti_data")
      showElement(selector = "#prism_input_options")
      showElement(selector = "#anti_header_upload")
      showElement(selector = "#prism_header_upload")
      showElement(selector = "#prism_data")
      showElement(selector = "#sempi_header_upload")
      showElement(selector = "#sempi_data")
      showElement(selector = "#deep_header_upload")
      showElement(selector = "#deep_data")
      showElement(selector = "#gecco_header_upload")
      showElement(selector = "#gecco_data")
      showElement(selector = "#rre_header_upload")
      showElement(selector = "#rre_data")
      showElement(selector = "#chr_len")
      showElement(selector = "#arts_header_upload")
      showElement(selector = "#known_data")
      showElement(selector = "#dup_data")
      showElement(selector = "#anti_sco")
      showElement(selector = "#prism_sco")
      showElement(selector = "#arts_sco")
      showElement(selector = "#rre_sco")
      showElement(selector = "#sempi_sco")
      showElement(selector = "#deep_sco")
      showElement(selector = "#gecco_sco")
  }
    })
  
  observeEvent(input$hide_anti, {
    if (input$hide_anti== T){
      hideElement(selector = "#anti_header")
      hideElement(selector = "#anti_hybrid")
      hideElement(selector = "#sempi_header")
      hideElement(selector = "#sempi_hybrid")
      hideElement(selector = "#prism_header")
      hideElement(selector = "#prism_hybrid")
      hideElement(selector = "#prism_supp")
      hideElement(selector = "#arts_header")
      hideElement(selector = "#dup_choice")
    }else{
      if (vals$anti_data_input == T){
        showElement(selector = "#anti_header")
        showElement(selector = "#anti_hybrid")
      } else{
        hideElement(selector = "#anti_header")
        hideElement(selector = "#anti_hybrid")
      }
      if (vals$prism_data_input == T){
      showElement(selector = "#prism_header")
      showElement(selector = "#prism_hybrid")
      if (vals$prism_json == T){
        showElement(selector = "#prism_supp")
      }
      } else {
        hideElement(selector = "#prism_header")
        hideElement(selector = "#prism_hybrid")
        hideElement(selector = "#prism_supp")
      }
      if (vals$sempi_data_input == T){
      showElement(selector = "#sempi_header")
      showElement(selector = "#sempi_hybrid")
      } else {
        hideElement(selector = "#sempi_header")
        hideElement(selector = "#sempi_hybrid")
      }
      if (vals$arts_data_input == T){
        showElement(selector = "#arts_header")
        showElement(selector = "#dup_choice")
      } else{
        hideElement(selector = "#arts_header")
        hideElement(selector = "#dup_choice")
      }
    }
  })
  
  observeEvent(input$hide_genes_on_chr, {
    if (input$hide_genes_on_chr == T){
      hideElement(selector = "#ref")
    } else {
      if (vals$data_upload_count > 0){
        showElement(selector = "#ref")
      } else {
        hideElement(selector = "#genes_on_chr")
        hideElement(selector = "#ref")
      }
    }
  })
  
  observeEvent(input$hide_summarize, {
    if (input$hide_summarize == T){
      hideElement(selector = "#group_by")
      hideElement(selector = "#count_all")
    } else {
      if (vals$data_upload_count > 1){
        showElement(selector = "#group_by")
        showElement(selector = "#count_all")
      } else {
        hideElement(selector = "#summarize")
        hideElement(selector = "#group_by")
        hideElement(selector = "#count_all")
      }
      
    }
  })
  
  observeEvent(input$hide_viz, {
    if (input$hide_viz == T){
      hideElement(selector = "#rename_data")
      hideElement(selector = "#rename")
      hideElement(selector = "#reset_name")
      hideElement(selector = "#rre_width")
      hideElement(selector = "#biocircos_color")
      hideElement(selector = "#label_color")
      hideElement(selector = "#label_color_class")
      hideElement(selector = "#ref_col_biocircos")
      hideElement(selector = "#arts_header")
    } else{
      showElement(selector = "#rename_data")
      showElement(selector = "#rename")
      showElement(selector = "#reset_name")
      if (vals$rre_data_input == T){
        showElement(selector = "#rre_width")
      } else {
        hideElement(selector = "#rre_width")
      }
      if (vals$prism_json == T){
        showElement(selector = "#prism_supp_width")
      }
      else {
        hideElement(selector = "#prism_supp_width")
      }
      if (vals$data_upload_count > 1){
        showElement(selector = "#biocircos_color")
        showElement(selector = "#label_color")
        showElement(selector = "#label_color_class")
      } else {
        hideElement(selector = "#biocircos_color")
        hideElement(selector = "#label_color")
        hideElement(selector = "#label_color_class")
      }
      if (input$label_color_class == "R"){
        showElement(selector = "#ref_col_biocircos")
      } else {
        hideElement(selector = "#ref_col_biocircos")
      }
      if (vals$arts_data_input == T){
        showElement(selector = "#arts_header")
      } else {
        hideElement(selector = "#arts_header")
      }
      
    }
  })
  
  observeEvent(input$hide_data_comparison, {
    if ((input$hide_data_comparison == T)){
      hideElement(selector = "#ref_comparison")
      hideElement(selector = "#score_type")
      hideElement(selector = "#plot_step")
      hideElement(selector = "#plot_start")
    } else if (vals$deep_data_input == T) {
      showElement(selector = "#ref_comparison")
      showElement(selector = "#score_type")
      showElement(selector = "#plot_step")
      showElement(selector = "#plot_start")
    } else {
      hideElement(selector = "#ref_comparison")
      hideElement(selector = "#score_type")
      hideElement(selector = "#plot_step")
      hideElement(selector = "#plot_start")
    }
  })
  
  observeEvent(input$hide_data_filter, {
    if ((input$hide_data_filter == T)){
      hideElement(selector = "#score_a")
      hideElement(selector = "#score_d")
      hideElement(selector = "#score_c")
      hideElement(selector = "#domains_filter")
      hideElement(selector = "#biodomain_filter")
      hideElement(selector = "#gene_filter")
      hideElement(selector = "#cluster_type")
    } else if  (vals$deep_data_input == T){
      showElement(selector = "#score_a")
      showElement(selector = "#score_d")
      showElement(selector = "#score_c")
      showElement(selector = "#domains_filter")
      showElement(selector = "#biodomain_filter")
      showElement(selector = "#gene_filter")
      showElement(selector = "#cluster_type")
    } else {
      hideElement(selector = "#score_a")
      hideElement(selector = "#score_d")
      hideElement(selector = "#score_c")
      hideElement(selector = "#domains_filter")
      hideElement(selector = "#biodomain_filter")
      hideElement(selector = "#gene_filter")
      hideElement(selector = "#cluster_type")
    }
  })
  
  observeEvent(input$hide_data_comparison_gecco, {
    if ((input$hide_data_comparison_gecco == T)){
      hideElement(selector = "#ref_comparison_gecco")
      hideElement(selector = "#score_type_gecco")
      hideElement(selector = "#plot_step_gecco")
      hideElement(selector = "#plot_start_gecco")
    } else if (vals$gecco_data_input == T) {
      showElement(selector = "#ref_comparison_gecco")
      showElement(selector = "#score_type_gecco")
      showElement(selector = "#plot_step_gecco")
      showElement(selector = "#plot_start_gecco")
    } else {
      hideElement(selector = "#ref_comparison_gecco")
      hideElement(selector = "#score_type_gecco")
      hideElement(selector = "#plot_step_gecco")
      hideElement(selector = "#plot_start_gecco")
    }
  })
  
  observeEvent(input$hide_data_filter_gecco, {
    if ((input$hide_data_filter_gecco == T)){
      hideElement(selector = "#score_average_gecco")
      hideElement(selector = "#score_average_gecco")
      hideElement(selector = "#domains_filter_gecco")
      hideElement(selector = "#prot_filter_gecco")
    } else if  (vals$gecco_data_input == T){
      showElement(selector = "#score_average_gecco")
      showElement(selector = "#score_average_gecco")
      showElement(selector = "#domains_filter_gecco")
      showElement(selector = "#prot_filter_gecco")
    } else {
      hideElement(selector = "#score_average_gecco")
      hideElement(selector = "#score_average_gecco")
      hideElement(selector = "#domains_filter_gecco")
      hideElement(selector = "#prot_filter_gecco")
    }
  })
  
  observeEvent(inputData(), {
    req(vals$data_upload_count >1)
    whereami::cat_where(whereami::whereami())
    # GENERATE DATA
    if (vals$anti_data_input == TRUE){
      anti_data <-  vals$anti_data
      anti_inter <- vals$anti_data %>%
        select(Start, Stop)
      anti_inter$seqnames <- "chr"
      
    }
    if (vals$deep_data_input == TRUE){
      deep_data <- vals$deep_data
      deep_inter <- vals$deep_data %>% 
        select(nucl_start, nucl_end)
      
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
        select(Start, Stop)
      rre_inter$seqnames <- "chr"
    }
    if (vals$prism_data_input == TRUE){
      # Store master prism data in local variable
      prism_data <- vals$prism_data
      # Start/Stop columns from prism data as matrix
      prism_inter <- prism_data %>%
        select(Start,Stop)
      prism_inter$seqnames <- "chr"
    }
    if (vals$sempi_data_input == TRUE){
      # Store master prism data in local variable
      sempi_data <- vals$sempi_data
      # Start/Stop columns from prism data as matrix
      sempi_inter <- vals$sempi_data %>%
        select(Start,Stop)
      sempi_inter$seqnames <- "chr"
    }
    if (vals$prism_supp_data_input == T){
      prism_supp_data <- vals$prism_supp_data
      prism_supp_inter <- vals$prism_supp_data %>%
        select(Start,Stop)
      prism_supp_inter$seqnames <- "chr"
      if (vals$prism_supp_data_input_width == TRUE) {
        Stop_vals_prism_supp <- as.numeric(vals$prism_supp$Stop)+50000
      } else{
        Stop_vals_prism_supp <- as.numeric(vals$prism_supp$Stop)
      }
    }
    if (vals$arts_data_input == T){
      arts_data <- vals$arts_data
      arts_inter <- vals$arts_data %>%
        select(Start,Stop) 
      arts_inter$seqnames <- "chr"
    }
    if (vals$gecco_data_input == TRUE){
      gecco_data <- vals$gecco_data
      # Start/Stop columns from prism data as matrix
      gecco_inter <- vals$gecco_data %>%
        select(Start,Stop)
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
    #inters_filtered <- vals$inters_filtered
    data_uploads <- c("anti_data_input","sempi_data_input","prism_data_input","prism_supp_data_input",
                      "arts_data_input","deep_data_input","gecco_data_input","rre_data_input")
    soft_names <- c("anti","sempi","prism","prism_supp","arts","deep","gecco","rre" )
    
   #TESTING
#     computed <- list(
#      anti=F,deep=F, gecco=F, arts=F, prism=F, sempi=F, prism_supp=F, rre=F
#    )
    
    
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
            #inters_filtered[[j]][[x]] <- new_res
            #inters_filtered[[x]][[j]] <- list(from=new_res$to, to=new_res$from)
            
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
   } else{
     vals$need_filter <- T
   }
   
   
#### TESTING
#   score_a <- apply(vals$deep_data %>% select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
#   score_d <- apply(vals$deep_data %>% select(c("deepbgc_score")),1, function(x) max(x))
#   score_c <- apply(vals$deep_data %>% select(c("alkaloid", "nrps","other","pks","ripp","saccharide","terpene")),1, function(x) max(x))
#   
#   deep_data <- vals$deep_data%>%
#     mutate(score_a = score_a, score_d = score_d, score_c = score_c) %>%
#     filter(score_a >= 0.5, score_c >=0.5 , score_d >= 0.8,  num_domains >= 5,
#            num_bio_domains>=3, num_proteins>=3)
#   vals <- list(
#     anti_data_input = T,
#     sempi_data_input = T,
#     prism_data_input = T,
#     prism_supp = T,
#     arts_data_input = T,
#     deep_data_input = T,
#     gecco_data_input = T,
#     rre_data_input = T,
#     anti_data = anti_data,
#     deep_data = deep_data,
#     prism_data = prism_data,
#     rre_data = rre_data,
#     arts_data = arts_data,
#     sempi_data = sempi_data,
#     gecco_data = gecco_data,
#     prism_supp_data = prism_supp,
#     rre_more = F,
#     biocircos_gecco = gecco_data,
#   computed=computed
#   )

  })
  
  observeEvent(dynamicInput(), {
    req(vals$data_upload_count>1)
    whereami::cat_where(whereami::whereami())
    inters <- vals$inters
    if (vals$deep_data_input == TRUE){
      # Get vector of max values from chosen columns from deepbgc data
      score_a <- apply(vals$deep_data %>% select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
      score_d <- apply(vals$deep_data %>% select(c("deepbgc_score")),1, function(x) max(x))
      score_c <- apply(vals$deep_data %>% select(c("alkaloid", "nrps","other","pks","ripp","saccharide","terpene")),1, function(x) max(x))
      deep_data_chromo <- vals$deep_data %>%
        mutate(score = apply(vals$deep_data %>%
                               select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, function(x) max(x))) 
      # Cluster_type column. Here extract colnames, and assign max value to a new column
      deep_data_chromo$Cluster_type <- colnames(deep_data_chromo %>% select(alkaloid, nrps, other, pks, ripp, saccharide, terpene))[apply(deep_data_chromo%>%select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, which.max) ]
      # If max score is under_threshold, print "under_threshold"
      deep_data_chromo <- deep_data_chromo%>%
        mutate(Cluster_type = ifelse(score>as.numeric(input$cluster_type)/100, Cluster_type, "under_threshold"))
      #Finally store deepbgc data in plotting variable. Do final scores processing 
      biocircos_deep <- deep_data_chromo%>%
        mutate( product_class = Cluster_type, score_a = score_a, score_d = score_d, score_c = score_c) %>%
        filter(score_a >= as.numeric(input$score_a )/ 100, score_c >=as.numeric(input$score_c)/100 , 
               score_d >= as.numeric(input$score_d)/100,  num_domains >= input$domains_filter,
               num_bio_domains>=input$biodomain_filter, num_proteins>=input$gene_filter)
      biocircos_deep['Start'] <- biocircos_deep$nucl_start
      biocircos_deep['Stop'] <- biocircos_deep$nucl_end
      biocircos_deep['Type'] <- biocircos_deep$product_class
      biocircos_deep['Type2'] <- biocircos_deep$product_class
      biocircos_deep['Cluster'] <- biocircos_deep$ID
      vals$deep_data_filtered <- biocircos_deep
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
    if (vals$gecco_data_input == TRUE){
      score_a_gecco <- apply(vals$gecco_data %>% select(c("average_p")),1, function(x) max(x))
      score_c_gecco <- apply(vals$gecco_data %>% select(c("alkaloid", "nrps","other","pks","ripp","saccharide","terpene")),1, function(x) max(x))
      # Store master prism data in local variable
      gecco_data <- vals$gecco_data %>%
        mutate(score = apply(vals$gecco_data %>%
                               dplyr::select(alkaloid, nrps, other, pks, ripp, saccharide, terpene),1, function(x) max(x))) %>%
        mutate(Cluster_type = ifelse(score>as.numeric(input$score_cluster_gecco)/100, Type2, "under_threshold")) %>%
        mutate( Type2 = Cluster_type, score_a = score_a_gecco, score_c = score_c_gecco) %>%
        filter(score_a >= as.numeric(input$score_average_gecco )/ 100, score_c >=as.numeric(input$score_cluster_gecco)/100 ,
               num_domains >= input$domains_filter_gecco, num_prot>=input$prot_filter_gecco)
      
      vals$gecco_data_filtered <- gecco_data 
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
    if (vals$arts_data_input == TRUE){
      if (input$dup_choice != "All"){
        vals$arts_data_filtered <- data.frame(vals$arts_data) %>%
          filter(Core == input$dup_choice | Core == "Not_core")
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
      } else {
        vals$arts_data_filtered <- vals$arts_data
      }
    }
    if ((vals$gecco_data_input == F) & (vals$deep_data_input == F )& (vals$arts_data_input == F )) {
      vals$inters_filtered <- inters
    }
    vals$need_filter <- F
  })
  
  observeEvent(biocircos_listen(), {
    req(vals$data_upload_count >=2)
    whereami::cat_where(whereami::whereami())
    #BioCircos!
    Biocircos_chromosomes <- list()
    arcs_chromosomes <- c()
    arcs_begin <- c()
    arcs_end <- c()
    arc_labels <- c()
    arc_col <- c()
    
    inters <- vals$inters_filtered
    
    rename_data <- vals$rename_data
    
    # ANTISMASH
    if (vals$anti_data_input == TRUE){
      # Store data in local variable
      biocircos_anti <- vals$anti_data
      #Make chromosome list for Biocircos plot. Use chr_len as an input
      Biocircos_chromosomes[["Antismash"]] <- vals$chr_len  
      #Add arcs. Quantity of arcs is length of dataframes
      arcs_chromosomes <- c(arcs_chromosomes,rep("Antismash", length(biocircos_anti$Cluster)) )
      # Add arcs begin positions. (Start column)
      arcs_begin <- c(arcs_begin, biocircos_anti$Start)
      # Stop position of arcs. 
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
    }
    
    #DEEPBGC
    if (vals$deep_data_input == TRUE){
      biocircos_deep <- vals$deep_data_filtered
      #Make chromosome list for Biocircos plot. Use chr_len as an input
      Biocircos_chromosomes[["DeepBGC"]] <- vals$chr_len
      #Add arcs. Quantity of arcs is length of dataframes
      arcs_chromosomes <- c(arcs_chromosomes, rep("DeepBGC", length(biocircos_deep$bgc_candidate_id)))
      # Add arcs begin positions. (Start column)
      arcs_begin <- c(arcs_begin, biocircos_deep$nucl_start )
      # Stop position of arcs. 
      arcs_end <- c(arcs_end, biocircos_deep$nucl_end)
      # Add Arcs labels. Can add only one label...
      arc_labels <- c(arc_labels, biocircos_deep$product_class)
      if ((input$biocircos_color == T)){
        arc_colors <- sapply(biocircos_deep$Type, function(x){
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
    }
    
    #RRE-FINDER
    if (vals$rre_data_input == TRUE){
      biocircos_rre <- data.frame(vals$rre_data)
      biocircos_rre$Start <- as.numeric(biocircos_rre$Start)
      biocircos_rre$Stop <- as.numeric(biocircos_rre$Stop)
      #Make chromosome list for Biocircos plot. Use chr_len as an input
      Biocircos_chromosomes[["RRE"]] <- vals$chr_len
      #Add arcs. Quantity of arcs is length of dataframes
      arcs_chromosomes <- c(arcs_chromosomes, rep("RRE", length(biocircos_rre$Locus_tag)))
      # Add arcs begin positions. (Start column)
      arcs_begin <- c(arcs_begin, biocircos_rre$Start)
      # Stop position of arcs. 
      if (input$rre_width == TRUE) {
        arcs_end <- c(arcs_end,  as.numeric(biocircos_rre$Stop)+50000)
      }else{
        arcs_end <- c(arcs_end,  as.numeric(biocircos_rre$Stop))
      }
      # Add Arcs labels. Can add only one label...
      arc_labels <- c(arc_labels,  biocircos_rre$E.value)
      if ((input$biocircos_color == T)){
        arc_colors <- sapply(biocircos_rre$Type, function(x){
          if (x %in% rename_data$Group_color){
            rename_data$Color[rename_data$Group_color == x]
          } else {
            rename_data$Color[rename_data$Group_color == 'base']
          }
        })
      } else {
        arc_colors <-rename_data$Color[rename_data$Group_color == 'base']
      }
      arc_col <- c(arc_col,as.character(arc_colors) )
    }
    
    # PRISM
    if (vals$prism_data_input == TRUE){
      # Store data in local variable
      biocircos_prism <- vals$prism_data
      #Make chromosome list for Biocircos plot. Use chr_len as an input
      Biocircos_chromosomes[["PRISM"]] <- vals$chr_len
      #Add arcs. Quantity of arcs is length of dataframes
      arcs_chromosomes <- c(arcs_chromosomes, rep("PRISM", length(biocircos_prism$Cluster)))
      # Add arcs begin positions. (Start column)
      arcs_begin <- c(arcs_begin, biocircos_prism$Start )
      # Stop position of arcs.
      arcs_end <- c(arcs_end, biocircos_prism$Stop )
      # Add Arcs labels. Can add only one label...
      arc_labels <- c(arc_labels,biocircos_prism$Type )
      if ((input$biocircos_color == T)){
        arc_colors <- sapply(biocircos_prism$Type2, function(x){
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
    }
    
    #SEMPI
    if (vals$sempi_data_input == TRUE){
      # Store data in local variable
      biocircos_sempi <- vals$sempi_data
      #Make chromosome list for Biocircos plot. Use chr_len as an input
      Biocircos_chromosomes[["SEMPI"]] <- vals$chr_len
      #Add arcs. Quantity of arcs is length of dataframes
      arcs_chromosomes <- c(arcs_chromosomes, rep("SEMPI", length(biocircos_sempi$Cluster)))
      
      # Add arcs begin positions. (Start column)
      arcs_begin <- c(arcs_begin, biocircos_sempi$Start )
      # Stop position of arcs.
      arcs_end <- c(arcs_end, biocircos_sempi$Stop )
      # Add Arcs labels. Can add only one label...
      arc_labels <- c(arc_labels,biocircos_sempi$Type )
      if ((input$biocircos_color == T)){
        arc_colors <- sapply(biocircos_sempi$Type2, function(x){
          if (x %in% rename_data$Group_color){
            rename_data$Color[rename_data$Group_color == x]
          } else {
            rename_data$Color[rename_data$Group_color == 'base']
          }
        })
      } else {
        arc_colors <-  rename_data$Color[rename_data$Group_color == 'base']
      }
      arc_col <- c(arc_col,as.character(arc_colors) )
    }
    
    if (vals$prism_supp_data_input == TRUE){
      # Store data in local variable
      biocircos_prism_supp <- vals$prism_supp
      #Make chromosome list for Biocircos plot. Use chr_len as an input
      Biocircos_chromosomes[["PRISM-supp"]] <- vals$chr_len
      #Add arcs. Quantity of arcs is length of dataframes
      arcs_chromosomes <- c(arcs_chromosomes, rep("PRISM-supp", length(biocircos_prism_supp$Cluster)))
      # Add arcs begin positions. (Start column)
      arcs_begin <- c(arcs_begin, biocircos_prism_supp$Start )
      if (vals$prism_supp_data_input_width == TRUE) {
        arcs_end <- c(arcs_end,  as.numeric(biocircos_prism_supp$Stop)+50000)
      }else{
        arcs_end <- c(arcs_end,  as.numeric(biocircos_prism_supp$Stop))
      }
      # Add Arcs labels. Can add only one label...
      arc_labels <- c(arc_labels,biocircos_prism_supp$Type )
      if ((input$biocircos_color == T)){
        arc_colors <- sapply(biocircos_prism_supp$Type2, function(x){
          if (x %in% rename_data$Group_color){
            rename_data$Color[rename_data$Group_color == x]
          } else {
            rename_data$Color[rename_data$Group_color == 'base']
          }
        })
      } else {
        arc_colors <-  rename_data$Color[rename_data$Group_color == 'base']
      }
      arc_col <- c(arc_col,as.character(arc_colors) )
    }
    
    
    if (vals$arts_data_input == TRUE){
      biocircos_arts <- data.frame(vals$arts_data_filtered)
      #Make chromosome list for Biocircos plot. Use chr_len as an input
      Biocircos_chromosomes[["ARTS"]] <- vals$chr_len
      #Add arcs. Quantity of arcs is length of dataframes
      arcs_chromosomes <- c(arcs_chromosomes, rep("ARTS", length(biocircos_arts$Cluster)))
      # Add arcs begin positions. (Start column)
      arcs_begin <- c(arcs_begin, biocircos_arts$Start)
      # Stop position of arcs. 
      if (input$arts_width == TRUE) {
        arcs_end <- c(arcs_end,  as.numeric(biocircos_arts$Stop)+20000)
      }else{
        arcs_end <- c(arcs_end,  as.numeric(biocircos_arts$Stop))
      }
      # Add Arcs labels. Can add only one label...
      arc_labels <- c(arc_labels,  biocircos_arts$Type2)
      if ((input$biocircos_color == T)){
        arc_colors <- sapply(biocircos_arts$Type2, function(x){
          if (x %in% rename_data$Group_color){
            rename_data$Color[rename_data$Group_color == x]
          } else {
            rename_data$Color[rename_data$Group_color == 'base']
          }
        })
      } else {
        arc_colors <-rename_data$Color[rename_data$Group_color == 'base']
      }
      arc_col <- c(arc_col,as.character(arc_colors) )
    }
    
    if (vals$gecco_data_input == TRUE){
      biocircos_gecco <- vals$gecco_data_filtered
      #Make chromosome list for Biocircos plot. Use chr_len as an input
      Biocircos_chromosomes[["GECCO"]] <- vals$chr_len  
      #Add arcs. Quantity of arcs is length of dataframes
      arcs_chromosomes <- c(arcs_chromosomes,rep("GECCO", length(biocircos_gecco$Cluster)) )
      # Add arcs begin positions. (Start column)
      arcs_begin <- c(arcs_begin, biocircos_gecco$Start)
      # Stop position of arcs. 
      arcs_end <- c(arcs_end, biocircos_gecco$Stop )
      # Add Arcs labels. Can add only one label...
      arc_labels <- c(arc_labels, biocircos_gecco$Type)
      if ((input$biocircos_color == T)){
        arc_colors <- sapply(biocircos_gecco$Type2, function(x){
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
    }
    
    # Add to tracklist. Then it can be populated with links
    tracklist <- BioCircosArcTrack('myArcTrack', arcs_chromosomes, arcs_begin, arcs_end, 
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
      link_pos_start <- as.numeric(c(data2$Start[data2$Cluster %in% inter_rre_s] ))
      # Add links start positions as a start from dataframe. For chromosome start variable
      link_pos_start_1 <- as.numeric(c(data2$Stop[data2$Cluster %in% inter_rre_s]))
      # Add links start position for a chromosome stop variable
      link_pos_end <- as.numeric(c( data1$Start[data1$Cluster %in% inter_s_rre_n]))
      # Add links start position for a chromosome stop position
      link_pos_end_2 <- as.numeric(c(data1$Stop[data1$Cluster %in% inter_s_rre_n]))
      label_1 <- c(sapply(inter_rre_s, function(x){x = paste(paste0(data2_label,":"), x, ",", data2$Type[data2$Cluster == x])})) 
      label_2 <- c(sapply(inter_s_rre_n, function(x){x = paste(paste0(data1_label, ":"), x, ",", data1$Type[data1$Cluster == x])}))
      if (!is.null(inter_rre_s)){
        if (class == 'P'){
          subset_vec <- data2$Type2[data2$Cluster %in% inter_rre_s] == data1$Type2[data1$Cluster %in% inter_s_rre_n]
          label_color <- as.character(c(sapply(data2$Type2[data2$Cluster %in% inter_rre_s], function (x){
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
          if (grep(data1_label, rename_data$Hierarchy) < (grep(data2_label, rename_data$Hierarchy))){
            label_color <- as.character(c(sapply(data1$Type2[data1$Cluster %in% inter_s_rre_n], function (x){
              if (x %in% rename_data$Group_color) {
                as.character(rename_data$Color[rename_data$Group_color == x])
              }
              else{
                as.character(rename_data$Color[rename_data$Group_color == 'base'])
              }
            })))
          } else {
            label_color <-as.character( c(sapply(data2$Type2[data2$Cluster %in% inter_rre_s], function (x){
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
            label_color <- as.character(c(sapply(data1$Type2[data1$Cluster %in% inter_s_rre_n], function (x){
              if (x %in% rename_data$Group_color) {
                as.character(rename_data$Color[rename_data$Group_color == x])
              }
              else{
                as.character(rename_data$Color[rename_data$Group_color == 'base'])
              }
            })))
          } else if (data1_label == input$ref_col_biocircos){
            label_color <- as.character(c(sapply(data2$Type2[data2$Cluster %in% inter_rre_s], function (x){
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
    
    
    
    if (vals$gecco_data_input == TRUE){
      if (vals$anti_data_input == TRUE){
        output <- add_biocircos_data(inters$gecco$anti$from, inters$gecco$anti$to, biocircos_anti, biocircos_gecco, "Antismash", "GECCO", rename_data, input$label_color_class)
        
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
      # Get interception of antismash with PRISM
      if (vals$prism_data_input == TRUE){
        output <- add_biocircos_data(inters$gecco$prism$from, inters$gecco$prism$to, biocircos_prism, biocircos_gecco, "PRISM", "GECCO", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Get interception of antismash with deepbgc
      if (vals$deep_data_input == TRUE){
        output <- add_biocircos_data(inters$gecco$deep$from, inters$gecco$deep$to, biocircos_deep, biocircos_gecco, "DeepBGC", "GECCO", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
        # Safe used local variables to the reactive ones
      } 
      # Get interception of antismash with RREFinder
      if (vals$rre_data_input == TRUE){
        output <- add_biocircos_data(inters$gecco$rre$from, inters$gecco$rre$to, biocircos_rre, biocircos_gecco, "RRE", "GECCO", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Get interception of antismash with SEMPI
      if (vals$sempi_data_input == TRUE){
        output <- add_biocircos_data(inters$gecco$sempi$from, inters$gecco$sempi$to, biocircos_sempi, biocircos_gecco, "SEMPI", "GECCO", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$prism_supp_data_input == TRUE){
        output <- add_biocircos_data(inters$gecco$prism_supp$from, inters$gecco$prism_supp$to, biocircos_prism_supp, biocircos_gecco, "PRISM-supp", "GECCO", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$arts_data_input == TRUE){
        output <- add_biocircos_data(inters$gecco$arts$from, inters$gecco$arts$to, biocircos_arts, biocircos_gecco, "ARTS", "GECCO", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Write csvs with locally used variables
      vals$biocircos_gecco <- biocircos_gecco
      write.csv(biocircos_gecco, "gecco_biocircos.csv", row.names = F)
    }
    
    # ANTISMASH
    if (vals$anti_data_input == TRUE){
      
      # Get interception of antismash with PRISM
      if (vals$prism_data_input == TRUE){
        output <- add_biocircos_data(inters$anti$prism$from, inters$anti$prism$to, biocircos_prism, biocircos_anti, "PRISM", "Antismash", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Get interception of antismash with deepbgc
      if (vals$deep_data_input == TRUE){
        output <- add_biocircos_data(inters$anti$deep$from, inters$anti$deep$to, biocircos_deep, biocircos_anti, "DeepBGC", "Antismash", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
        # Safe used local variables to the reactive ones
      } 
      # Get interception of antismash with RREFinder
      if (vals$rre_data_input == TRUE){
        output <- add_biocircos_data(inters$anti$rre$from, inters$anti$rre$to, biocircos_rre, biocircos_anti, "RRE", "Antismash", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Get interception of antismash with SEMPI
      if (vals$sempi_data_input == TRUE){
        output <- add_biocircos_data(inters$anti$sempi$from, inters$anti$sempi$to, biocircos_sempi, biocircos_anti, "SEMPI", "Antismash", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$prism_supp_data_input == TRUE){
        output <- add_biocircos_data(inters$anti$prism_supp$from, inters$anti$prism_supp$to, biocircos_prism_supp, biocircos_anti, "PRISM-supp", "Antismash", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$arts_data_input == TRUE){
        output <- add_biocircos_data(inters$anti$arts$from, inters$anti$arts$to, biocircos_arts, biocircos_anti, "ARTS", "Antismash", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Write csvs with locally used variables
      write.csv(biocircos_anti, "antismash_biocircos.csv", row.names = F)
    }
    
    # DEEPBGC 
    if (vals$deep_data_input == TRUE){
      
      # Get interception of DeepBGC with rrefinder
      if (vals$rre_data_input == TRUE){
        output <- add_biocircos_data(inters$deep$rre$from, inters$deep$rre$to, biocircos_rre, biocircos_deep, "RRE", "DeepBGC", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
        # Safe used local variables to the reactive ones
      }
      # Get interception of DeepBGC with PRISM
      if (vals$prism_data_input == TRUE){
        output <- add_biocircos_data(inters$deep$prism$from, inters$deep$prism$to, biocircos_prism, biocircos_deep, "PRISM", "DeepBGC", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
        # Safe used local variables to the reactive ones
      }
      # Get interception of DeepBGC with SEMPI
      if (vals$sempi_data_input == TRUE){
        output <- add_biocircos_data(inters$deep$sempi$from, inters$deep$sempi$to, biocircos_sempi, biocircos_deep, "SEMPI", "DeepBGC", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$prism_supp_data_input == TRUE){
        output <- add_biocircos_data(inters$deep$prism_supp$from, inters$deep$prism_supp$to, biocircos_prism_supp, biocircos_deep, "PRISM-supp", "DeepBGC", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$arts_data_input == TRUE){
        output <- add_biocircos_data(inters$deep$arts$from, inters$deep$arts$to, biocircos_arts, biocircos_deep, "ARTS", "DeepBGC", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Safe used local variables to the reactive ones
      vals$biocircos_deep <- biocircos_deep
      # Write csvs with locally used variables
      write.csv(biocircos_deep, "deepbgc_biocircos.csv", row.names = F)
    }
    
    # PRISM
    if (vals$prism_data_input == TRUE){
      
      # Get interception of PRISM with rrefinder
      if (vals$rre_data_input == TRUE){
        output <- add_biocircos_data(inters$prism$rre$from, inters$prism$rre$to, biocircos_rre, biocircos_prism, "RRE", "PRISM", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Get interception of PRISM with SEMPI
      if (vals$sempi_data_input == TRUE){
        output <- add_biocircos_data(inters$prism$sempi$from, inters$prism$sempi$to, biocircos_sempi, biocircos_prism, "SEMPI", "PRISM", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$prism_supp_data_input == TRUE){
        output <- add_biocircos_data(inters$prism$prism_supp$from, inters$prism$prism_supp$to, biocircos_prism_supp, biocircos_prism, "PRISM-supp", "PRISM", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$arts_data_input == TRUE){
        output <- add_biocircos_data(inters$prism$arts$from, inters$prism$arts$to, biocircos_arts, biocircos_prism, "ARTS", "PRISM", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Write csvs with locally used variables
      write.csv(biocircos_prism, "prism_biocircos.csv", row.names = F)
    }
    
    # RRE-FINDER 
    if (vals$rre_data_input == TRUE){
      
      # Get interception of RRE with SEMPI
      if (vals$sempi_data_input == TRUE){
        output <- add_biocircos_data(inters$rre$sempi$from, inters$rre$sempi$to, biocircos_sempi, biocircos_rre, "SEMPI", "RRE", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$prism_supp_data_input == TRUE){
        output <- add_biocircos_data(inters$rre$prism_supp$from, inters$rre$prism_supp$to, biocircos_prism_supp, biocircos_rre, "PRISM-supp", "RRE", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$arts_data_input == TRUE){
        output <- add_biocircos_data(inters$rre$arts$from, inters$rre$arts$to, biocircos_arts, biocircos_rre, "ARTS", "RRE", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Write csvs with locally used variables
      write.csv(biocircos_rre, "rre_biocircos.csv", row.names = F)
    }
    
    #SEMPI 
    if (vals$sempi_data_input == TRUE){
      if (vals$prism_supp_data_input == TRUE){
        output <- add_biocircos_data(inters$sempi$prism_supp$from, inters$sempi$prism_supp$to,biocircos_prism_supp, biocircos_sempi, "PRISM-supp", "SEMPI", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      if (vals$arts_data_input == TRUE){
        output <- add_biocircos_data(inters$sempi$arts$from, inters$sempi$arts$to,biocircos_arts, biocircos_sempi, "ARTS", "SEMPI", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      # Write csvs with locally used variables
      write.csv(biocircos_sempi, "sempi_biocircos.csv", row.names = F)
    }
    
    if (vals$prism_supp_data_input == TRUE){
      if (vals$arts_data_input == TRUE){
        output <- add_biocircos_data(inters$prism_supp$arts$from, inters$prism_supp$arts$to,biocircos_arts, biocircos_prism_supp, "ARTS", "PRISM-supp", rename_data, input$label_color_class)
        # Add link start. Just populate certain chromosome name times the lenght of interception 
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
      #Write csvs with locally used variables
      write.csv(biocircos_prism_supp, "prism_supp_biocircos.csv", row.names = F)
    }
    
    if (vals$arts_data_input == TRUE){
      #Write csvs with locally used variables
      write.csv(biocircos_arts, "arts_biocircos.csv", row.names = F)
    }
    
    
    
    
    # Combine labels with mapply to one list
    link_labels <- mapply(function(x,y)  paste(x, y, sep = " | "), label_1, label_2 )
    
    # Add links and labels to the track list for subsequent visualization 
    if (input$label_color == T){
      group_colors <- count(unlist(label_color))
      for (i in seq(1:dim(group_colors)[1])){
        subset <- unname( which(label_color %in% group_colors$x[i]))
        tracklist = tracklist + BioCircosLinkTrack(as.character(i), chromosomes_start[subset], link_pos_start[subset], 
                                                   link_pos_start_1[subset], chromosomes_end[subset], link_pos_end[subset], 
                                                   link_pos_end_2[subset], maxRadius = 0.85, labels = link_labels[subset],
                                                   displayLabel = FALSE, color = group_colors$x[i])
      }
    } else{
      tracklist = tracklist + BioCircosLinkTrack('myLinkTrack_master', chromosomes_start, link_pos_start, 
                                                 link_pos_start_1, chromosomes_end, link_pos_end, 
                                                 link_pos_end_2, maxRadius = 0.85, labels = link_labels,
                                                 displayLabel = FALSE, color = rename_data$Color[rename_data$Group_color == 'base'])
    }
    
    vals$tracklist <- tracklist
    vals$Biocircos_chromosomes <- Biocircos_chromosomes
  })
  #Render output plots

  # Render barplot
  output$deep_barplot <- renderPlot({
    whereami::cat_where(whereami::whereami())
    
    # Create empty dataframe to populate later
    fullnes_of_annotation <- data.frame(NA, NA, NA)
    colnames(fullnes_of_annotation) <- c("Score", "Source", "Quantity")
    fullnes_of_annotation <- drop_na(fullnes_of_annotation)
    
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
        filter(score>=dataframe_1/100) %>%
        select(Start, Stop) 
      if (length(deep_inter$Start) > 0) {
        deep_inter$seqnames <- "chr"
      }
      
      
      # Store antismash bgc start amd atop values as matrix
      if (input$ref_comparison == 'A'){
        anti_inter <- vals$anti_data %>%
        select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison == 'P'){
        anti_inter <- vals$prism_data %>%
          select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison == 'S'){
        anti_inter <- vals$sempi_data %>%
          select(Start, Stop) 
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
        title <-  ggtitle("Comparison of Antismash and DeepBGC annotations at given score threshold")
      } else if (input$ref_comparison == 'P'){
        used_antismash <-  length(vals$prism_data$Cluster)-inter_bgc
        cols <- c("Only PRISM", "DeepBGC+PRISM", "Only DeepBGC")
        title <- ggtitle("Comparison of PRISM and DeepBGC annotations at given score threshold")
      } else if (input$ref_comparison == 'S') {
        used_antismash <-  length(vals$sempi_data$Cluster)-inter_bgc
        cols <- c("Only SEMPI", "DeepBGC+SEMPI", "Only DeepBGC")
        title <- ggtitle("Comparison of SEMPI and DeepBGC annotations at given score threshold")
      }

      # Combine all vectors into one dataframe
      fullnes_of_annotation_1 <- data.frame(c(rep(c(as.character(dataframe_1)),3 )), 
                                            cols, c(used_antismash, inter_bgc, len_new))
      colnames(fullnes_of_annotation_1) <- c("Score", "Source", "Quantity")
      # Combine previously created empty dataframe with this one to store results
      fullnes_of_annotation <- rbind(fullnes_of_annotation, fullnes_of_annotation_1)
      
    }
    
    # Store dataframe in reactive value for later use.
    vals$fullness <- data.frame(fullnes_of_annotation)
    write.csv(fullnes_of_annotation, "fullness.csv", row.names = F)
    
    # Make text to show on a barplot to point on additional scores' thresholds
    annotateText=paste("Applied additional thresholds", paste("Activity score:", as.character(input$score_a)),
                       paste("DeepBGC score:", as.character(input$score_d)),
                       paste("Cluster type score:", as.character(input$score_c)), sep = "\n")
    
    # Plot the barplot
    ggplot(fullnes_of_annotation, aes(fill=Source, y=Quantity, x=Score)) + 
      geom_bar(position="dodge", stat="identity")+
      geom_text(aes(label=Quantity), position=position_dodge(width=0.9), vjust=-0.25) +
      xlab(paste(input$score_type,"Score")) +
      title +
      geom_label(aes(x=Inf,y=Inf,hjust=1,vjust=1,label=annotateText ), show.legend = F)
  })
  
  # Render interactive plot with plotly for rates of DeepBGC data in regards with antismash data
  output$deep_rate <- renderPlotly({
    whereami::cat_where(whereami::whereami())
    
    # Reuse stored dataframe from previous plot
    # This dataframe stores data for number of intercepted/non intercepted clusters for DeepBGC and antismash data 
    # For more information please see previous renderPlot
    fullnes_of_annotation <- data.frame(vals$fullness)
    
    # Store dataframe into variable. Widen it to calculate rates
    test <- fullnes_of_annotation %>%
      pivot_wider(names_from = Source, values_from = Quantity)
    if (input$ref_comparison == 'A'){
      data <-  vals$anti_data
      title <- ggtitle("Rates of DeepBGC/Antismash data annotation")
      test <- test %>%
        # Calculate rates. Novelty is nummber of clusters annotated only by deepbgc/ all clusters annotated by antismash + (antismash + deepbgc)
        mutate(Novelty_rate = test$`Only DeepBGC`/(test$`DeepBGC+Antismash` + test$`Only Antismash`), 
               #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
               Annotation_rate = test$`DeepBGC+Antismash`/length(data$Cluster), 
               # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
               Skip_rate = test$`Only Antismash`/length(data$Cluster))
    } else if (input$ref_comparison == 'P'){
      data <- vals$prism_data
      title <- ggtitle("Rates of DeepBGC/PRISM data annotation")
      test <- test %>%
        mutate(Novelty_rate = test$`Only DeepBGC`/(test$`DeepBGC+PRISM` + test$`Only PRISM`), 
               #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
               Annotation_rate = test$`DeepBGC+PRISM`/length(data$Cluster), 
               # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
               Skip_rate = test$`Only PRISM`/length(data$Cluster))
    } else if (input$ref_comparison == 'S'){
      data <- vals$sempi_data
      title <- ggtitle("Rates of DeepBGC/SEMPI data annotation")
      test <- test %>%
        mutate(Novelty_rate = test$`Only DeepBGC`/(test$`DeepBGC+SEMPI` + test$`Only SEMPI`), 
               #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
               Annotation_rate = test$`DeepBGC+SEMPI`/length(data$Cluster), 
               # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
               Skip_rate = test$`Only SEMPI`/length(data$Cluster))
    }
    
    # Calculate rates and plot interactive plot with plotly
    ggplotly(test %>%
               pivot_longer(cols = c(Novelty_rate, Annotation_rate, Skip_rate), names_to = 'Rates', values_to = 'Rates_data') %>%
               ggplot(aes(x=as.numeric(Score), y=as.numeric(Rates_data), Rate = as.numeric(Rates_data))) +
               geom_line(aes(color=Rates)) +
               geom_point(aes(shape=Rates), alpha = .4, size = 3) +
               title +
               ylab("Rate") +
               xlab(paste(input$score_type,"Score threshold")),
             tooltip = c("Rate"))
  })
  
  # Render barplot
  output$gecco_barplot <- renderPlot({
    whereami::cat_where(whereami::whereami())
    # Create empty dataframe to populate later
    fullnes_of_annotation <- data.frame(NA, NA, NA)
    colnames(fullnes_of_annotation) <- c("Score", "Source", "Quantity")
    fullnes_of_annotation <- drop_na(fullnes_of_annotation)
    
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

      # Filter dataframe. Get only rows, which >= of a given thresholds. Select only start and stop of those rows as a matrix
      gecco_inter <- gecco_inter_1 %>%
        filter(score>=dataframe_1/100) %>%
        select(Start, Stop) 
      if (length(gecco_inter$Start) > 0) {
        gecco_inter$seqnames <- "chr"
      }
      
      
      # Store antismash bgc start amd atop values as matrix
      if (input$ref_comparison_gecco == 'A'){
        anti_inter <- vals$anti_data %>%
          select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison_gecco == 'P'){
        anti_inter <- vals$prism_data %>%
          select(Start, Stop) 
        anti_inter$seqnames <- "chr"
      } else if (input$ref_comparison_gecco == 'S'){
        anti_inter <- vals$sempi_data %>%
          select(Start, Stop) 
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
        title <-  ggtitle("Comparison of Antismash and GECCO annotations at given score threshold")
      } else if (input$ref_comparison_gecco == 'P'){
        used_antismash <-  length(vals$prism_data$Cluster)-inter_bgc
        cols <- c("Only PRISM", "GECCO+PRISM", "Only GECCO")
        title <- ggtitle("Comparison of PRISM and GECCO annotations at given score threshold")
      } else if (input$ref_comparison_gecco == 'S') {
        used_antismash <-  length(vals$sempi_data$Cluster)-inter_bgc
        cols <- c("Only SEMPI", "GECCO+SEMPI", "Only GECCO")
        title <- ggtitle("Comparison of SEMPI and GECCO annotations at given score threshold")
      }
      
      # Combine all vectors into one dataframe
      fullnes_of_annotation_1 <- data.frame(c(rep(c(as.character(dataframe_1)),3 )), 
                                            cols, c(used_antismash, inter_bgc, len_new))
      colnames(fullnes_of_annotation_1) <- c("Score", "Source", "Quantity")
      # Combine previously created empty dataframe with this one to store results
      fullnes_of_annotation <- rbind(fullnes_of_annotation, fullnes_of_annotation_1)
      
    }
    
    # Store dataframe in reactive value for later use.
    vals$fullness <- data.frame(fullnes_of_annotation)
    write.csv(fullnes_of_annotation, "fullness.csv", row.names = F)
    
    # Make text to show on a barplot to point on additional scores' thresholds
    annotateText=paste("Applied additional thresholds", paste("Average p-value:", as.character(input$score_average_gecco)),
                       paste("Cluster type score:", as.character(input$score_cluster_gecco)), sep = "\n")
    
    # Plot the barplot
    ggplot(fullnes_of_annotation, aes(fill=Source, y=Quantity, x=Score)) + 
      geom_bar(position="dodge", stat="identity")+
      geom_text(aes(label=Quantity), position=position_dodge(width=0.9), vjust=-0.25) +
      xlab(paste(input$score_type,"Score")) +
      title +
      geom_label(aes(x=Inf,y=Inf,hjust=1,vjust=1,label=annotateText ), show.legend = F)
  })
  
  # Render interactive plot with plotly for rates of DeepBGC data in regards with antismash data
  output$gecco_rate <- renderPlotly({
    whereami::cat_where(whereami::whereami())
    # Reuse stored dataframe from previous plot
    # This dataframe stores data for number of intercepted/non intercepted clusters for DeepBGC and antismash data 
    # For more information please see previous renderPlot
    fullnes_of_annotation <- data.frame(vals$fullness)
    
    # Store dataframe into variable. Widen it to calculate rates
    test <- fullnes_of_annotation %>%
      pivot_wider(names_from = Source, values_from = Quantity)
    if (input$ref_comparison_gecco == 'A'){
      data <-  vals$anti_data
      title <- ggtitle("Rates of GECCO/Antismash data annotation")
      test <- test %>%
        # Calculate rates. Novelty is nummber of clusters annotated only by deepbgc/ all clusters annotated by antismash + (antismash + deepbgc)
        mutate(Novelty_rate = test$`Only GECCO`/(test$`GECCO+Antismash` + test$`Only Antismash`), 
               #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
               Annotation_rate = test$`GECCO+Antismash`/length(data$Cluster), 
               # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
               Skip_rate = test$`Only Antismash`/length(data$Cluster))
    } else if (input$ref_comparison_gecco == 'P'){
      data <- vals$prism_data
      title <- ggtitle("Rates of GECCO/PRISM data annotation")
      test <- test %>%
        mutate(Novelty_rate = test$`Only GECCO`/(test$`GECCO+PRISM` + test$`Only PRISM`), 
               #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
               Annotation_rate = test$`GECCO+PRISM`/length(data$Cluster), 
               # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
               Skip_rate = test$`Only PRISM`/length(data$Cluster))
    } else if (input$ref_comparison_gecco == 'S'){
      data <- vals$sempi_data
      title <- ggtitle("Rates of GECCO/SEMPI data annotation")
      test <- test %>%
        mutate(Novelty_rate = test$`Only GECCO`/(test$`GECCO+SEMPI` + test$`Only SEMPI`), 
               #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
               Annotation_rate = test$`GECCO+SEMPI`/length(data$Cluster), 
               # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
               Skip_rate = test$`Only SEMPI`/length(data$Cluster))
    }
    
    # Calculate rates and plot interactive plot with plotly
    ggplotly(test %>%
               pivot_longer(cols = c(Novelty_rate, Annotation_rate, Skip_rate), names_to = 'Rates', values_to = 'Rates_data') %>%
               ggplot(aes(x=as.numeric(Score), y=as.numeric(Rates_data), Rate = as.numeric(Rates_data))) +
               geom_line(aes(color=Rates)) +
               geom_point(aes(shape=Rates), alpha = .4, size = 3) +
               title +
               ylab("Rate") +
               xlab(paste(input$score_type,"Score threshold")),
             tooltip = c("Rate"))
  })
  
  # Render interactive plot, which shows bgcs of antismash, intercepted with chosen app. Also all app bgs. On hover shows all available information
  # For antismash and PRISM data showed only ID, Start, Stop, Type
  output$deep_reference <- renderPlotly({

    whereami::cat_where(whereami::whereami())
   # req(vals$data_upload_count >=1)
    data_uploads <- c("anti_data_input","sempi_data_input","prism_data_input","prism_supp_data_input",
                      "arts_data_input","deep_data_input","gecco_data_input","rre_data_input")
    soft_names <- c("anti","sempi","prism","prism_supp","arts","deep","gecco","rre" )
    soft <- c("Antismash","SEMPI","PRISM","PRISM_SUPPORT","ARTS","DeepBGC","GECCO","RRE-Finder" )
    abbr <- c("A", "S", "P", "P-supp", "AR", "D", "G", "RRE")
    
    inters <- vals$inters_filtered
    
    # GENERATE DATA
    if (vals$anti_data_input == TRUE){
      # Store antismash data in local variable, with column renaming
      anti_data <-  vals$anti_data 
      # Extract only Start and Stop from antismash data into matrix
      anti_inter <- vals$anti_data %>%
        select(Start, Stop)
      anti_inter$seqnames <- "chr"
      
      }
    if (vals$deep_data_input == TRUE){
      deep_data <- vals$deep_data_filtered
      deep_inter <- deep_data %>% 
        select(nucl_start, nucl_end)
      
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
        select(Start, Stop)
      rre_inter$seqnames <- "chr"
      if (input$rre_width == TRUE) {
        Stop_vals_RRE <- as.numeric(vals$rre_data$Stop)+50000
      } else{
        Stop_vals_RRE <- as.numeric(vals$rre_data$Stop)
      }
    }
    if (vals$prism_data_input == TRUE){
      # Store master prism data in local variable
      prism_data <- vals$prism_data
      # Start/Stop columns from prism data as matrix
      prism_inter <- prism_data %>%
        select(Start,Stop)
      prism_inter$seqnames <- "chr"
    }
    if (vals$sempi_data_input == TRUE){
      # Store master prism data in local variable
      sempi_data <- vals$sempi_data
      # Start/Stop columns from prism data as matrix
      sempi_inter <- vals$sempi_data %>%
        select(Start,Stop)
      sempi_inter$seqnames <- "chr"
    }
    if (vals$prism_supp_data_input == T){
      prism_supp_data <- vals$prism_supp
      prism_supp_inter <- prism_supp %>%
        select(Start,Stop)
      prism_supp_inter$seqnames <- "chr"
      if (vals$prism_supp_data_input_width == TRUE) {
        Stop_vals_prism_supp <- as.numeric(vals$prism_supp$Stop)+50000
      } else{
        Stop_vals_prism_supp <- as.numeric(vals$prism_supp$Stop)
      }
    }
    if (vals$arts_data_input == T){
      arts_data <- vals$arts_data_filtered
      arts_inter <- arts_data %>%
        select(Start,Stop) 
      arts_inter$seqnames <- "chr"
      if (input$arts_width == TRUE) {
        Stop_vals_arts <- as.numeric(arts_data$Stop)+20000
      } else{
        Stop_vals_arts<- as.numeric(arts_data$Stop)
      }
    }
    if (vals$gecco_data_input == TRUE){
     gecco_data <- vals$gecco_data_filtered
      # Start/Stop columns from prism data as matrix
      gecco_inter <- gecco_data %>%
        select(Start,Stop)
      gecco_inter$seqnames <- "chr"
    }
    lett <- rev(LETTERS)[1:9]
    simple_seg <- function(df, letter, software, soft_name ,soft){
      data <- df[df$Cluster %in% inters[[soft]][[soft_name]]$from, ]
      
      seg_df <-  data.frame(x=as.numeric(data$Start),
                         y=rep(letter, length(data$Cluster)),
                         xend=as.numeric(  data$Stop),
                         yend=rep(letter, length(data$Cluster)),
                         Type = as.factor(data$Type),
                         Type2 = as.factor(data$Type),
                         Software = rep(software, length(data$Cluster)),
                         ID = data$Cluster,
                         Start = data$Start,
                         Stop = data$Stop)
      return(seg_df)
      }
    
    geom_anti <- function(data){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, 
                                  ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    geom_prism <- function(data){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,                                                                                                             ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    geom_deep <- function(data){
      geom_segment(data=data,aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                               ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
                                                               deepbgc_score = deepbgc_score,activity = activity ),size =3)
      }
    geom_rre <- function(data){
      if (vals$rre_more == T){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type, Score = Score, Software = Software,
                                                           ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value,
                                                           P_value = P_value, RRE_start = RRE_start,RRE_stop = RRE_stop, 
                                                           Probability = Probability),size = 3)
      } else {
        geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                    ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value
                                    ),size = 3)
      }
      }
    geom_sempi <- function(data){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,
                                  ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    geom_prism_supp <- function(data){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, ID = ID,
                                  Start = Start, Stop = Stop, Type = Type, Name = Name, Full_name = Full_name,
                                  Score = Score), size = 3)
    }
    geom_arts <- function(data){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,                                                                                                             ID = ID, Start = Start, Stop = Stop, Type = Type, Hit = Hit, 
                                  Core = Core, E_value = E_value, Bitscore = Bitscore, Count = Count, Model = Model), size = 3)
    }
    geom_gecco <- function(data){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, 
                                  ID = ID, Start = Start, Stop = Stop, Type = Type, Num_proteins= Num_proteins,
                                  Num_domains = Num_domains,Average_p = Average_p, Max_p = Max_p ), size = 3)
    }
    
    tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                "P_value", "RRE_start","RRE_stop", "Probability", "Name", "Full_name",  "Hit", "Core", "Count", "Bitscore", "Model",
                "Num_domains", "Num_proteins", "Average_p", "Max_p")

    
    
    add_arts <- function(seg_df, soft, df){
      if (input$arts_width == TRUE) {
        Stop_vals_arts_in <- as.numeric(df[df$Cluster %in% inters[[soft]]$arts$from, ]$Stop)+20000
      } else{
        Stop_vals_arts_in <- as.numeric(df[df$Cluster %in% inters[[soft]]$arts$from, ]$Stop)
      }
      seg_df$Hit = df[df$Cluster %in% inters[[soft]]$arts$from, ]$Hit
      seg_df$xend = Stop_vals_arts_in
      seg_df$Core = df[df$Cluster %in% inters[[soft]]$arts$from, ]$Core
      seg_df$Count = df[df$Cluster %in% inters[[soft]]$arts$from, ]$Count
      seg_df$E_value = df[df$Cluster %in% inters[[soft]]$arts$from, ]$Evalue
      seg_df$Bitscore = df[df$Cluster %in% inters[[soft]]$arts$from, ]$Bitscore
      seg_df$Model = df[df$Cluster %in% inters[[soft]]$arts$from, ]$Model
      return(seg_df)
    }
    add_prism_supp <- function(seg_df, soft, df){
      if (vals$prism_supp_data_input_width == TRUE) {
        Stop_vals_prism_supp_in <- as.numeric(df[df$Cluster %in% inters[[soft]]$prism_supp$from, ]$Stop)+50000
      } else{
        Stop_vals_prism_supp_in <- as.numeric(df[df$Cluster %in% inters[[soft]]$prism_supp$from, ]$Stop)
      }
      seg_df$xend <- Stop_vals_prism_supp_in
      seg_df$Score = df[df$Cluster %in% inters[[soft]]$prism_supp$from, ]$Score
      seg_df$Name = df[df$Cluster %in% inters[[soft]]$prism_supp$from, ]$Name
      seg_df$Full_name = df[df$Cluster %in% inters[[soft]]$prism_supp$from, ]$Full_name
      return(seg_df)
    }
    add_deep <- function(seg_df, soft, df){
      seg_df$num_domains = df[df$Cluster %in% inters[[soft]]$deep$from, ]$num_domains
      seg_df$deepbgc_score = df[df$Cluster %in% inters[[soft]]$deep$from, ]$deepbgc_score
      seg_df$activity = df[df$Cluster %in% inters[[soft]]$deep$from, ]$product_activity
      return(seg_df)
    }
    add_rre <- function(seg_df, soft, df){
      if (input$rre_width == TRUE) {
        Stop_vals_RRE_in <- as.numeric(df[df$Cluster %in% inters[[soft]]$rre$from, ]$Stop)+50000
      } else{
        Stop_vals_RRE_in <- as.numeric(df[df$Cluster %in% inters[[soft]]$rre$from, ]$Stop)
      }
      if (vals$rre_more == T){
        seg_df$xend=Stop_vals_RRE_in
        seg_df$Score = df[df$Cluster %in% inters[[soft]]$rre$from, ]$Score
        seg_df$Stop = df[df$Cluster %in% inters[[soft]]$rre$from, ]$Stop
        seg_df$E_value = df[df$Cluster %in% inters[[soft]]$rre$from, ]$E.value
        seg_df$P_value = df[df$Cluster %in% inters[[soft]]$rre$from, ]$P.value
        seg_df$RRE_start = df[df$Cluster %in% inters[[soft]]$rre$from, ]$RRE.start
        seg_df$RRE_stop = df[df$Cluster %in% inters[[soft]]$rre$from, ]$RRE.end
        seg_df$Probability = df[df$Cluster %in% inters[[soft]]$rre$from, ]$Probability
      } else {
        seg_df$xend=Stop_vals_RRE_in
        seg_df$E_value = df[df$Cluster %in% inters[[soft]]$rre$from, ]$E.value
      }
      
     return(seg_df)
    }
    add_gecco <- function(seg_df, soft, df){
      seg_df$Num_proteins = df[df$Cluster %in% inters[[soft]]$gecco$from, ]$num_prot
      seg_df$Num_domains = df[df$Cluster %in% inters[[soft]]$gecco$from, ]$num_domains
      seg_df$Average_p = df[df$Cluster %in% inters[[soft]]$gecco$from, ]$average_p
      seg_df$Max_p = df[df$Cluster %in% inters[[soft]]$gecco$from, ]$max_p
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
    
    define_spec_seg_df <- function(soft_names, index,seg_df, soft_major, df ){
      if ((soft_names[index] == "prism_supp") & (soft_names[index] != soft_major)){
        seg_df <-  add_prism_supp(seg_df, soft_major, df)
      } else if ((soft_names[index] == "arts")& (soft_names[index] != soft_major)){
        seg_df <- add_arts(seg_df, soft_major, df)
      } else if ((soft_names[index] == "deep")& (soft_names[index] != soft_major)){
        seg_df <- add_deep(seg_df, soft_major, df)
      } else if ((soft_names[index] == "gecco")& (soft_names[index] != soft_major)){
        seg_df <- add_gecco(seg_df, soft_major, df)
      } else if ((soft_names[index] == "rre")& (soft_names[index] != soft_major)){
        seg_df <- add_rre(seg_df, soft_major, df)
      }
      return(seg_df)
    }
    
    # MAKE COMPUTATIONS
    if (vals$gecco_data_input == TRUE){
      seg_ref_g <- data.frame(x=as.numeric(vals$biocircos_gecco$Start),
                              y=rep("Z", length(vals$biocircos_gecco$Cluster)),
                              xend=as.numeric(vals$biocircos_gecco$Stop),
                              yend=rep("Z", length(vals$biocircos_gecco$Cluster)),
                              Type = as.factor(vals$biocircos_gecco$Type),
                              Type2 = as.factor(vals$biocircos_gecco$Type2),
                              Software = rep("GECCO", length(vals$biocircos_gecco$Cluster)),
                              ID = vals$biocircos_gecco$Cluster,
                              Start = vals$biocircos_gecco$Start,
                              Stop = vals$biocircos_gecco$Stop,
                              Num_proteins = vals$biocircos_gecco$num_prot,
                              Num_domains = vals$biocircos_gecco$num_domains,
                              Average_p = vals$biocircos_gecco$average_p,
                              Max_p = vals$biocircos_gecco$max_p)
      seg_ref <- seg_ref_g
      
      if (input$ref == "GECCO") {
        plot <- ggplot(vals$biocircos_gecco, aes(x = vals$chr_len, y = Chr)) + geom_gecco(seg_ref)
        soft_major <- "gecco"
        soft_let <- "G"
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
          scale_y_discrete(labels = c("Z" = input$ref, unlist(labels_1))) +
          theme(axis.text.y = element_text(size = 10)) +
          ylab("")+
          xlab("Chromosome length")+ 
          theme(legend.title = element_blank()) +
          ggtitle("Annotations' comparison to the reference")
        to_plot <- ggplotly(plot, tooltip = tooltip)
        to_plot <- to_plot %>% 
          layout(legend=list(font = list(
            family = "sans-serif",
            size = 12,
            color = "#000"),
            bordercolor = "#FFFFFF",
            borderwidth = 2,
            title=list(text='<b> Cluster Types </b>')))
        
        
      }
      vals$seg_df_ref_g <- data.frame(x=as.numeric(vals$biocircos_gecco$Start),
                                      y=rep("S", length(vals$biocircos_gecco$Cluster)),
                                      xend=as.numeric(vals$biocircos_gecco$Stop),
                                      yend=rep("S", length(vals$biocircos_gecco$Cluster)),
                                      Type = as.factor(vals$biocircos_gecco$Type),
                                      Type2 = as.factor(vals$biocircos_gecco$Type2),
                                      Software = rep("GECCO", length(vals$biocircos_gecco$Cluster)),
                                      ID = vals$biocircos_gecco$Cluster,
                                      Start = vals$biocircos_gecco$Start,
                                      Stop = vals$biocircos_gecco$Stop,
                                      Num_proteins = vals$biocircos_gecco$num_prot,
                                      Num_domains = vals$biocircos_gecco$num_domains,
                                      Average_p = vals$biocircos_gecco$average_p,
                                      Max_p = vals$biocircos_gecco$max_p)
    }
    if (vals$anti_data_input == TRUE){
      seg_ref_a <- data.frame(x=as.numeric(anti_data$Start),
                              y=rep("Z", length(anti_data$Cluster)),
                              xend=as.numeric(anti_data$Stop),
                              yend=rep("Z", length(anti_data$Cluster)),
                              Type = as.factor(anti_data$Type),
                              Type2 = as.factor(anti_data$Type2),
                              Software = rep("Antismash", length(anti_data$Cluster)),
                              ID = anti_data$Cluster,
                              Start = anti_data$Start,
                              Stop = anti_data$Stop)
      seg_ref <- seg_ref_a
      if (input$ref == "Antismash") {
        plot <- ggplot(anti_data, aes(x = vals$chr_len, y = Chr)) + geom_anti(seg_ref)
        soft_major <- "anti"
        soft_let <- "A"
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
          scale_y_discrete(labels = c("Z" = input$ref, unlist(labels_1))) +
          theme(axis.text.y = element_text(size = 10)) +
          ylab("")+
          xlab("Chromosome length")+ 
          theme(legend.title = element_blank()) +
          ggtitle("Annotations' comparison to the reference")
        to_plot <- ggplotly(plot, tooltip = tooltip)
        to_plot <- to_plot %>% 
          layout(legend=list(font = list(
            family = "sans-serif",
            size = 12,
            color = "#000"),
            bordercolor = "#FFFFFF",
            borderwidth = 2,
            title=list(text='<b> Cluster Types </b>')))
        
      }
      vals$seg_df_ref_a <- seg_ref_a
    }
    if (vals$deep_data_input == TRUE){
      # Create a dataframe with all deepbgc data + the additional info to visualize on hover
      seg_ref_d <- data.frame(x=as.numeric(  deep_data$nucl_start),
                              y=rep("Z", length(deep_data$ID)),
                              xend=as.numeric(  deep_data$nucl_end),
                              yend=rep("Z", length(deep_data$ID)),
                              Type = as.factor(deep_data$product_class),
                              Software = rep("DeepBGC", length(deep_data$ID)),
                              ID = deep_data$ID,
                              Start = deep_data$Start,
                              Stop = deep_data$Stop,
                              num_domains = deep_data$num_domains,
                              deepbgc_score = deep_data$deepbgc_score,
                              activity = deep_data$product_activity)
      seg_ref <- seg_ref_d
      
      if (input$ref == "DeepBGC") {
        plot <- ggplot(deep_data_chromo, aes(x = vals$chr_len, y = Chr)) + geom_deep(seg_ref)
        soft_major <- "deep"
        soft_let <- "D"
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
          scale_y_discrete(labels = c("Z" = input$ref, unlist(labels_1))) +
          theme(axis.text.y = element_text(size = 10)) +
          ylab("")+
          xlab("Chromosome length")+ 
          theme(legend.title = element_blank()) +
          ggtitle("Annotations' comparison to the reference")
        to_plot <- ggplotly(plot, tooltip = tooltip)
        to_plot <- to_plot %>% 
          layout(legend=list(font = list(
            family = "sans-serif",
            size = 12,
            color = "#000"),
            bordercolor = "#FFFFFF",
            borderwidth = 2,
            title=list(text='<b> Cluster Types </b>')))
      }
     
      vals$seg_df_ref_d <- data.frame(x=as.numeric(  deep_data$nucl_start),
                                      y=rep("X", length(deep_data$ID)),
                                      xend=as.numeric(  deep_data$nucl_end),
                                      yend=rep("X", length(deep_data$ID)),
                                      Type = as.factor(deep_data$product_class),
                                      Software = rep("DeepBGC", length(deep_data$ID)),
                                      ID = deep_data$ID,
                                      Start = deep_data$Start,
                                      Stop = deep_data$Stop,
                                      num_domains = deep_data$num_domains,
                                      deepbgc_score = deep_data$deepbgc_score,
                                      activity = deep_data$product_activity)
    }
    
    if (vals$rre_data_input == TRUE){
      if (vals$rre_more == T){
      seg_ref_r <- data.frame(x=vals$rre_data$Start,
                              y=rep("Z", length(vals$rre_data$Locus_tag)),
                              xend=Stop_vals_RRE,
                              yend=rep("Z", length(vals$rre_data$Locus_tag)),
                              Type = rep("ripp", length(vals$rre_data$Locus_tag)),
                              Score = vals$rre_data$Score,
                              Software = rep("RREFinder", length(vals$rre_data$Locus_tag)),
                              ID = vals$rre_data$Locus_tag,
                              Start = vals$rre_data$Start,
                              Stop = vals$rre_data$Stop,
                              E_value = vals$rre_data$E.value,
                              P_value = vals$rre_data$P.value,
                              RRE_start = vals$rre_data$RRE.start,
                              RRE_stop = vals$rre_data$RRE.end,
                              Probability = vals$rre_data$Probability)
      } else {
        seg_ref_r <- data.frame(x=vals$rre_data$Start,
                                y=rep("Z", length(vals$rre_data$Locus_tag)),
                                xend=Stop_vals_RRE,
                                yend=rep("Z", length(vals$rre_data$Locus_tag)),
                                Type = rep("ripp", length(vals$rre_data$Locus_tag)),
                                Software = rep("RREFinder", length(vals$rre_data$Locus_tag)),
                                ID = vals$rre_data$Locus_tag,
                                Start = vals$rre_data$Start,
                                Stop = vals$rre_data$Stop,
                                E_value = vals$rre_data$E.value)
      }
      seg_ref <- seg_ref_r
      if (input$ref == "RRE-Finder") {
        plot <- ggplot(rre_data, aes(x = vals$chr_len, y = Chr)) + geom_rre(seg_ref)
        soft_major <- "rre"
        soft_let <- "RRE"
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
          scale_y_discrete(labels = c("Z" = input$ref, unlist(labels_1))) +
          theme(axis.text.y = element_text(size = 10)) +
          ylab("")+
          xlab("Chromosome length")+ 
          theme(legend.title = element_blank()) +
          ggtitle("Annotations' comparison to the reference")
        to_plot <- ggplotly(plot, tooltip = tooltip)
        to_plot <- to_plot %>% 
          layout(legend=list(font = list(
            family = "sans-serif",
            size = 12,
            color = "#000"),
            bordercolor = "#FFFFFF",
            borderwidth = 2,
            title=list(text='<b> Cluster Types </b>')))
      }
      
      if (vals$rre_more == T){
      vals$seg_df_ref_r <- data.frame(x=vals$rre_data$Start,
                                      y=rep("Y", length(vals$rre_data$Locus_tag)),
                                      xend=Stop_vals_RRE,
                                      yend=rep("Y", length(vals$rre_data$Locus_tag)),
                                      Type = rep("ripp", length(vals$rre_data$Locus_tag)),
                                      Score = vals$rre_data$Score,
                                      Software = rep("RREFinder", length(vals$rre_data$Locus_tag)),
                                      ID = vals$rre_data$Locus_tag,
                                      Start = vals$rre_data$Start,
                                      Stop = vals$rre_data$Stop,
                                      E_value = vals$rre_data$E.value,
                                      P_value = vals$rre_data$P.value,
                                      RRE_start = vals$rre_data$RRE.start,
                                      RRE_stop = vals$rre_data$RRE.end,
                                      Probability = vals$rre_data$Probability)
      } else {
        vals$seg_df_ref_r <- data.frame(x=vals$rre_data$Start,
                                        y=rep("Y", length(vals$rre_data$Locus_tag)),
                                        xend=Stop_vals_RRE,
                                        yend=rep("Y", length(vals$rre_data$Locus_tag)),
                                        Type = rep("ripp", length(vals$rre_data$Locus_tag)),
                                        Software = rep("RREFinder", length(vals$rre_data$Locus_tag)),
                                        ID = vals$rre_data$Locus_tag,
                                        Start = vals$rre_data$Start,
                                        Stop = vals$rre_data$Stop,
                                        E_value = vals$rre_data$E.value)
      }
    }
    if (vals$prism_data_input == TRUE){
      # Create a dataframe with PRISM data with all the additional info to visualize on hover
      seg_ref_p <- data.frame(x=as.numeric(prism_data$Start),
                                y=rep("Z", length(prism_data$Cluster)),
                                xend=as.numeric(prism_data$Stop),
                                yend=rep("Z", length(prism_data$Cluster)),
                                Type = as.factor(prism_data$Type),
                                Type2 = as.factor(prism_data$Type2),
                                Software = rep("PRISM", length(prism_data$Cluster)),
                                ID = prism_data$Cluster,
                                Start = prism_data$Start,
                                Stop = prism_data$Stop)
      seg_ref <- seg_ref_p 
      
      if (input$ref == "PRISM") {
      plot <- ggplot(prism_data, aes(x = vals$chr_len, y = Chr)) + geom_prism(seg_ref)
      soft_major <- "prism"
      soft_let <- "P"
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
        scale_y_discrete(labels = c("Z" = input$ref, unlist(labels_1))) +
        theme(axis.text.y = element_text(size = 10)) +
        ylab("")+
        xlab("Chromosome length")+ 
        theme(legend.title = element_blank()) +
        ggtitle("Annotations' comparison to the reference")
      to_plot <- ggplotly(plot, tooltip = tooltip)
      to_plot <- to_plot %>% 
        layout(legend=list(font = list(
          family = "sans-serif",
          size = 12,
          color = "#000"),
          bordercolor = "#FFFFFF",
          borderwidth = 2,
          title=list(text='<b> Cluster Types </b>')))
      }
      vals$seg_df_ref_p <- data.frame(x=as.numeric(prism_data$Start),
                                      y=rep("W", length(prism_data$Cluster)),
                                      xend=as.numeric(prism_data$Stop),
                                      yend=rep("W", length(prism_data$Cluster)),
                                      Type = as.factor(prism_data$Type),
                                      Type2 = as.factor(prism_data$Type2),
                                      Software = rep("PRISM", length(prism_data$Cluster)),
                                      ID = prism_data$Cluster,
                                      Start = prism_data$Start,
                                      Stop = prism_data$Stop)
    }
    if (vals$sempi_data_input == TRUE){
      # Create a dataframe with sempi data with all the additional info to visualize on hover
      seg_ref_s <- data.frame(x=as.numeric(sempi_data$Start),
                              y=rep("Z", length(sempi_data$Cluster)),
                              xend=as.numeric(sempi_data$Stop),
                              yend=rep("Z", length(sempi_data$Cluster)),
                              Type = as.factor(sempi_data$Type),
                              Type2 = as.factor(sempi_data$Type2),
                              Software = rep("SEMPI", length(sempi_data$Cluster)),
                              ID = sempi_data$Cluster,
                              Start = sempi_data$Start,
                              Stop = sempi_data$Stop)
      seg_ref <- seg_ref_s
      
      if (input$ref == "SEMPI") {
        plot <- ggplot(sempi_data, aes(x = vals$chr_len, y = Chr)) + geom_sempi(seg_ref)
        soft_major <- "sempi"
        soft_let <- "S"
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
          scale_y_discrete(labels = c("Z" = input$ref, unlist(labels_1))) +
          theme(axis.text.y = element_text(size = 10)) +
          ylab("")+
          xlab("Chromosome length")+ 
          theme(legend.title = element_blank()) +
          ggtitle("Annotations' comparison to the reference")
        to_plot <- ggplotly(plot, tooltip = tooltip)
        to_plot <- to_plot %>% 
          layout(legend=list(font = list(
            family = "sans-serif",
            size = 12,
            color = "#000"),
            bordercolor = "#FFFFFF",
            borderwidth = 2,
            title=list(text='<b> Cluster Types </b>')))
      }
      vals$seg_df_ref_s <- data.frame(x=as.numeric(sempi_data$Start),
                                      y=rep("V", length(sempi_data$Cluster)),
                                      xend=as.numeric(sempi_data$Stop),
                                      yend=rep("V", length(sempi_data$Cluster)),
                                      Type = as.factor(sempi_data$Type),
                                      Type2 = as.factor(sempi_data$Type2),
                                      Software = rep("SEMPI", length(sempi_data$Cluster)),
                                      ID = sempi_data$Cluster,
                                      Start = sempi_data$Start,
                                      Stop = sempi_data$Stop)
    }
    if (vals$prism_supp_data_input == TRUE){
      # Create a dataframe with sempi data with all the additional info to visualize on hover
      seg_ref_s <- data.frame(x=prism_supp$Start,
                              y=rep("Z", length(prism_supp$Cluster)),
                              xend=Stop_vals_prism_supp,
                              yend=rep("Z", length(prism_supp$Cluster)),
                              Type = as.factor(prism_supp$Type),
                              Type2 = as.factor(prism_supp$Type2),
                              Score = prism_supp$Score,
                              Software = rep("PRISM", length(prism_supp$Cluster)),
                              ID = prism_supp$Cluster,
                              Start = prism_supp$Start,
                              Stop = prism_supp$Stop,
                              Name = prism_supp$Name,
                              Full_name = prism_supp$Full_name
      )
      seg_ref <- seg_ref_s
      
      if (input$ref == "PRISM-supp") {
        plot <- ggplot(prism_supp, aes(x = vals$chr_len, y = Chr)) + geom_prism_supp(seg_ref)
        soft_major <- "prism_supp"
        soft_let <- "PS"
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
          scale_y_discrete(labels = c("Z" = input$ref, unlist(labels_1))) +
          theme(axis.text.y = element_text(size = 10)) +
          ylab("")+
          xlab("Chromosome length")+ 
          theme(legend.title = element_blank()) +
          ggtitle("Annotations' comparison to the reference")
        to_plot <- ggplotly(plot, tooltip = tooltip)
        to_plot <- to_plot %>% 
          layout(legend=list(font = list(
            family = "sans-serif",
            size = 12,
            color = "#000"),
            bordercolor = "#FFFFFF",
            borderwidth = 2,
            title=list(text='<b> Cluster Types </b>')))
      }
      vals$seg_df_ref_p_s <- data.frame(x=prism_supp$Start,
                                      y=rep("U", length(prism_supp$Cluster)),
                                      xend=Stop_vals_prism_supp,
                                      yend=rep("U", length(prism_supp$Cluster)),
                                      Type = as.factor(prism_supp$Type),
                                      Type2 = as.factor(prism_supp$Type2),
                                      Score = prism_supp$Score,
                                      Software = rep("PRISM", length(prism_supp$Cluster)),
                                      ID = prism_supp$Cluster,
                                      Start = prism_supp$Start,
                                      Stop = prism_supp$Stop,
                                      Name = prism_supp$Name,
                                      Full_name = prism_supp$Full_name
      )
    }
    if (vals$arts_data_input == TRUE){
      # Create a dataframe with sempi data with all the additional info to visualize on hover
      seg_ref_s <- data.frame(x=arts_data$Start,
                              y=rep("Z", length(arts_data$Cluster)),
                              xend=Stop_vals_arts,
                              yend=rep("Z", length(arts_data$Cluster)),
                              Type = as.factor(arts_data$Type),
                              Type2 = as.factor(arts_data$Type2),
                              Hit = arts_data$Hit,
                              Software = rep("ARTS", length(arts_data$Cluster)),
                              ID = arts_data$Cluster,
                              Start = arts_data$Start,
                              Stop = arts_data$Stop,
                              Core = arts_data$Core,
                              Count = arts_data$Count,
                              E_value = arts_data$Evalue,
                              Bitscore = arts_data$Bitscore,
                              Model = arts_data$Model
      )
      seg_ref <- seg_ref_s
      
      if (input$ref == "ARTS") {
        plot <- ggplot(arts_data, aes(x = vals$chr_len, y = Chr)) + geom_arts(seg_ref)
        soft_major <- "arts"
        soft_let <- "AR"
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
          scale_y_discrete(labels = c("Z" = input$ref, unlist(labels_1))) +
          theme(axis.text.y = element_text(size = 10)) +
          ylab("")+
          xlab("Chromosome length")+ 
          theme(legend.title = element_blank()) +
          ggtitle("Annotations' comparison to the reference")
        to_plot <- ggplotly(plot, tooltip = tooltip)
        to_plot <- to_plot %>% 
          layout(legend=list(font = list(
          family = "sans-serif",
          size = 12,
          color = "#000"),
          bordercolor = "#FFFFFF",
          borderwidth = 2,
          title=list(text='<b> Cluster Types </b>')))
      }
      vals$seg_df_ref_ar <- data.frame(x=arts_data$Start,
                                      y=rep("T", length(arts_data$Cluster)),
                                      xend=Stop_vals_arts,
                                      yend=rep("T", length(arts_data$Cluster)),
                                      Type = as.factor(arts_data$Type),
                                      Type2 = as.factor(arts_data$Type2),
                                      Hit = arts_data$Hit,
                                      Software = rep("ARTS", length(arts_data$Cluster)),
                                      ID = arts_data$Cluster,
                                      Start = arts_data$Start,
                                      Stop = arts_data$Stop,
                                      Core = arts_data$Core,
                                      Count = arts_data$Count,
                                      E_value = arts_data$Evalue,
                                      Bitscore = arts_data$Bitscore,
                                      Model = arts_data$Model
      )
    }

    to_plot
  })
  
  output$deep_reference_2 <- renderPlotly({
    whereami::cat_where(whereami::whereami())
    if (vals$rre_data_input == TRUE){
      data <- data.frame(vals$rre_data)
    } else if (vals$arts_data_input == TRUE){
      data <- vals$arts_data_filtered
    } else if (vals$anti_data_input == TRUE){
      data <-  vals$anti_data %>%
        mutate(ID = Cluster, Chr = chromosome) %>%
        dplyr::select(ID,Chr ,Start, Stop, Type, Type2)
    }else if (vals$deep_data_input == TRUE){
      data <- vals$deep_data_filtered
    }else if (vals$prism_data_input == TRUE){
      data <- vals$prism_data
    }else if (vals$sempi_data_input == TRUE){
      data <- vals$sempi_data
    }else if (vals$gecco_data_input == TRUE){
      data <- vals$gecco_data_filtered
    }
    
    tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                "P_value", "RRE_start","RRE_stop", "Probability", "Name", "Full_name",  "Hit", "Core", "Count", "Bitscore", "Model",
                "Num_domains", "Num_proteins", "Average_p", "Max_p")
    
    plot <- ggplot(data, aes(x = vals$chr_len, y = Chr))
    if (vals$anti_data_input == TRUE){
      plot <- plot + 
        geom_segment(data=vals$seg_df_ref_a, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,
                                                 ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    if (vals$deep_data_input == TRUE){
      plot <- plot +
        geom_segment(data=vals$seg_df_ref_d,aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                      ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
                                      deepbgc_score = deepbgc_score,activity = activity ),size =3)
    }
    if (vals$rre_data_input == TRUE){
      if (vals$rre_more == T){
      plot <- plot + geom_segment(data=vals$seg_df_ref_r, aes(x, y, xend=xend, yend=yend, color = Type, Score = Score, Software = Software,
                                                    ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value,
                                                    P_value = P_value, RRE_start = RRE_start,RRE_stop = RRE_stop, 
                                                    Probability = Probability),size = 3)
      } else {
        plot <- plot + geom_segment(data=vals$seg_df_ref_r, aes(x, y, xend=xend, yend=yend, color = Type,  Software = Software,
                                                                ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value),size = 3)
      }
    }
    if (vals$prism_data_input == TRUE){
      plot <- plot + geom_segment(data=vals$seg_df_ref_p, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,
                                                    ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
      
      
    }
    if (vals$sempi_data_input == TRUE){
      plot <- plot + geom_segment(data=vals$seg_df_ref_s, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software,
                                                              ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
      
      
    }
    if (vals$prism_supp_data_input == TRUE){
      plot <- plot + geom_segment(data=vals$seg_df_ref_p_s, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, ID = ID,
                                                            Start = Start, Stop = Stop, Type = Type, Name = Name, Full_name = Full_name,
                                                            Score = Score), size = 3)
    }
    if (vals$arts_data_input == TRUE){
      plot <- plot + geom_segment(data=vals$seg_df_ref_ar, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, 
                                                                ID = ID, Start = Start, Stop = Stop, Type = Type, Hit = Hit,
                                                                Core = Core, E_value = E_value, Bitscore = Bitscore, Count = Count,
                                                                Model = Model), size = 3)
    }
    if (vals$gecco_data_input == TRUE){
      plot <- plot + geom_segment(data =  vals$seg_df_ref_g, aes(x, y, xend=xend, yend=yend, color = Type2, Software = Software, 
                                                                 ID = ID, Start = Start, Stop = Stop, Type = Type, Num_proteins= Num_proteins,
                                                                 Num_domains = Num_domains,Average_p = Average_p, Max_p = Max_p ), size = 3) 
    }
    to_plot <- ggplotly(plot +
                          scale_y_discrete(labels = c("Z" = "Antismash","X" = "DeepBGC", "Y" = "RRE", "W" = "PRISM", 
                                                      "V" = "SEMPI", "U" = "P-supp", "T" = "ARTS", "S" = "GECCO")) +
                          theme(axis.text.y = element_text(size = 10)) +
                          ylab("")+
                          xlab("Chromosome length")+
                          theme(legend.title = element_blank()) +
                          ggtitle("All annotations"), 
                        # What actually to visualize in tooltip
                        tooltip = tooltip
    )
    to_plot %>% layout(legend=list(font = list(
      family = "sans-serif",
      size = 12,
      color = "#000"),
      bordercolor = "#FFFFFF",
      borderwidth = 2,
      title=list(text='<b> Cluster Types </b>')))
  })
  
  # Render Biocircos Plot for all-vs-all comparison
  output$biocircos <- renderBioCircos({
    whereami::cat_where(whereami::whereami())
    # Plot BioCircos
    BioCircos(vals$tracklist, genome = vals$Biocircos_chromosomes, genomeTicksScale = 1e+6)
  })
  
  
  output$biocircos_legend <- renderDataTable({
    whereami::cat_where(whereami::whereami())
    plot_data <- vals$rename_data
    new_data <- drop_na(data.frame(cbind(as.character(plot_data$Group_color), as.character(plot_data$Color))) )
    new_data <- new_data[!apply(new_data == "", 1, all),]
    colnames(new_data) <- c("Name", "Color")
    color_vec <- new_data$Color
    options(DT.options = list(pageLength = 50))
    datatable(new_data, rownames = F) %>% formatStyle('Color',
                                                      backgroundColor=styleEqual(color_vec, color_vec))
    
    
  })
  
  # Render barplot with number count of interception for BGC IDs
  output$barplot_rank <- renderPlotly({
    whereami::cat_where(whereami::whereami())
    antismash_count <-  NULL
    prism_count <- NULL
    deep_count <- NULL
    rre_count <- NULL
    sempi_count <- NULL
    prism_supp_count <- NULL
    arts_count <- NULL
    gecco_count <- NULL

    inters <- vals$inters_filtered
    
    if (vals$anti_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      antismash_count <-count(as.factor(unlist(sapply(inters$anti, function(x){x$to}))))
      # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
      anti_anot <- vals$anti_data[vals$anti_data$Cluster %in% as.numeric(levels(antismash_count$x)),]
      # Add prefices to the ID to plot for a barplot.  
      antismash_count$x <- sapply(antismash_count$x, function(x) paste("A: ", x))
      # Add label column to the dataframe, from which we will plot  
      antismash_count$label <- rep("Antismash", length(antismash_count$x))
      # Add type to the dataframe, from which we would plot (from annotation dataframe)  
      antismash_count$Type <- anti_anot$Type
      # Add Start positions (to visualize on hover)
      antismash_count$Start <- anti_anot$Start
      # Add Stop positions (to visualize on hover)
      antismash_count$Stop <- anti_anot$Stop
    }
    if (vals$gecco_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      gecco_count <- count(as.factor(unlist(sapply(inters$gecco, function(x){x$to}))))
      # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
      gecco_anot <- vals$gecco_data_filtered[vals$gecco_data_filtered$Cluster %in% as.numeric(levels(gecco_count$x)),]
      # Add prefices to the ID to plot for a barplot.  
      gecco_count$x <- sapply(gecco_count$x, function(x) paste("G: ", x))
      # Add label column to the dataframe, from which we will plot  
      gecco_count$label <- rep("GECCO", length(gecco_count$x))
      # Add type to the dataframe, from which we would plot (from annotation dataframe)  
      gecco_count$Type <- gecco_anot$Type
      # Add Start positions (to visualize on hover)
      gecco_count$Start <- gecco_anot$Start
      # Add Stop positions (to visualize on hover)
      gecco_count$Stop <- gecco_anot$Stop
    }
    if (vals$deep_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      deep_count <- count(as.factor(unlist(sapply(inters$deep, function(x){x$to}))))
      # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
      deep_anot <- vals$deep_data_filtered[vals$deep_data_filtered$ID %in% as.numeric(levels(deep_count$x)),]
      # Add prefices to the ID to plot for a barplot.  
      deep_count$x <- sapply(deep_count$x, function(x) paste("D: ", x))
      # Add label column to the dataframe, from which we will plot  
      deep_count$label <- rep("DeepBGC", length(deep_count$x))
      # Add type to the dataframe, from which we would plot (from annotation dataframe)  
      deep_count$Type <- deep_anot$Type
      deep_count$Start <- deep_anot$Start
      # Add Stop positions (to visualize on hover)
      deep_count$Stop <- deep_anot$Stop
    }
    if (vals$rre_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      rre_count <- count(as.factor(unlist(sapply(inters$rre, function(x){x$to}))))
      # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
      rre_anot <- vals$rre_data[vals$rre_data$ID %in% as.numeric(levels(rre_count$x)),]
      # Add prefices to the ID to plot for a barplot.  
      rre_count$x <- sapply(rre_count$x, function(x) paste("RRE: ", x))
      # Add label column to the dataframe, from which we will plot  
      rre_count$label <- rep("RRE", length(rre_count$x))
      # Add type to the dataframe, from which we would plot (from annotation dataframe)  
      rre_count$Type <- rep("ripp", length(rre_anot$Sequence))
      # Add Start positions (to visualize on hover)
      rre_count$Start <- rre_anot$Start
      # Add Stop positions (to visualize on hover)
      rre_count$Stop <- rre_anot$Stop
    }
    if (vals$prism_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      prism_count <- count(as.factor(unlist(sapply(inters$prism, function(x){x$to}))))
      # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
      prism_anot <- vals$prism_data[vals$prism_data$Cluster %in% as.numeric(levels(prism_count$x)),]
      # Add prefices to the ID to plot for a barplot.  
      prism_count$x <- sapply(prism_count$x, function(x) paste("P: ", x))
      # Add label column to the dataframe, from which we will plot  
      prism_count$label <- rep("PRISM", length(prism_count$x))
      # Add type to the dataframe, from which we would plot (from annotation dataframe)  
      prism_count$Type <- prism_anot$Type
      # Add Start positions (to visualize on hover)
      prism_count$Start <- prism_anot$Start
      # Add Stop positions (to visualize on hover)
      prism_count$Stop <- prism_anot$Stop
    }
    if (vals$sempi_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      sempi_count <- count(as.factor(unlist(sapply(inters$sempi, function(x){x$to}))))
      # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
      sempi_anot <- vals$sempi_data[vals$sempi_data$Cluster %in% as.numeric(levels(sempi_count$x)),]
      # Add prefices to the ID to plot for a barplot.  
      sempi_count$x <- sapply(sempi_count$x, function(x) paste("S: ", x))
      # Add label column to the dataframe, from which we will plot  
      sempi_count$label <- rep("SEMPI", length(sempi_count$x))
      # Add type to the dataframe, from which we would plot (from annotation dataframe)  
      sempi_count$Type <- sempi_anot$Type
      # Add Start positions (to visualize on hover)
      sempi_count$Start <- sempi_anot$Start
      # Add Stop positions (to visualize on hover)
      sempi_count$Stop <- sempi_anot$Stop
    }
    if (vals$prism_supp_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      prism_supp_count <- count(as.factor(unlist(sapply(inters$prism_supp, function(x){x$to}))))
      # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
      prism_supp_anot <- vals$prism_supp[vals$prism_supp$Cluster %in% as.numeric(levels(prism_supp_count$x)),]
      # Add prefices to the ID to plot for a barplot.  
      prism_supp_count$x <- sapply(prism_supp_count$x, function(x) paste("PS: ", x))
      # Add label column to the dataframe, from which we will plot  
      prism_supp_count$label <- rep("PRISM-supp", length(prism_supp_count$x))
      # Add type to the dataframe, from which we would plot (from annotation dataframe)  
      prism_supp_count$Type <- prism_supp_anot$Type
      # Add Start positions (to visualize on hover)
      prism_supp_count$Start <- prism_supp_anot$Start
      # Add Stop positions (to visualize on hover)
      prism_supp_count$Stop <- prism_supp_anot$Stop
    }
    if (vals$arts_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      arts_count <- count(as.factor(unlist(sapply(inters$arts, function(x){x$to}))))
      # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
      arts_anot <- vals$arts_data_filtered[vals$arts_data_filtered$Cluster %in% as.numeric(levels(arts_count$x)),]
      # Add prefices to the ID to plot for a barplot.  
      arts_count$x <- sapply(arts_count$x, function(x) paste("AR: ", x))
      # Add label column to the dataframe, from which we will plot  
      arts_count$label <- rep("ARTS", length(arts_count$x))
      # Add type to the dataframe, from which we would plot (from annotation dataframe)  
      arts_count$Type <- arts_anot$Type
      # Add Start positions (to visualize on hover)
      arts_count$Start <- arts_anot$Start
      # Add Stop positions (to visualize on hover)
      arts_count$Stop <- arts_anot$Stop
    }

      
    # Integrate all those dataframe to the master one 
    ranking_data <- rbind(gecco_count, antismash_count,prism_count, deep_count,rre_count, sempi_count, prism_supp_count, arts_count)
    # Fix column names in the master dataframe
    colnames(ranking_data) <- c("Cluster", "Count", "Label", "Type", "Start", "Stop")
    # Plot
    ggplotly(ggplot(ranking_data, aes(x = Cluster, y = Count, Type = Type, Start = Start, Stop = Stop)) +
               geom_bar(stat = "identity", aes(fill = Label)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
                     axis.text.y = element_text(size = 14)) +
               ggtitle("Number of times cluster is annotated with other tool"),
             tooltip=c("Type", "Start", "Stop")  
             )

    
  })
  
  # Render table with data
  output$group_table <- renderTable({
    whereami::cat_where(whereami::whereami())
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
        data <- df_test %>% group_by_if(colnames(df_test)==soft_let[i]) %>% summarise(a = paste(eval(as.name(excluded_names[1])), collapse=","),
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
              select(-as.name(soft_namings[f]))
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
  output$download <- downloadHandler(filename = function(){
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
shinyApp(ui = ui, server = server)
