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
library(IntervalSurgeon)
library(plotly)
library(BioCircos)
library(ggplot2)
library(reticulate)
library(shinyjs)
library(rjson)
library(stringr)
library(RSQLite)



use_virtualenv("clinker")
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
      h5("ANTISMASH:"),
      checkboxInput("anti_input_options", "My AnriSMASH data is a dataframe, not json results file from antismash"),
      fileInput("anti_data",
                "Upload antismash data"),
      h5("PRISM:"),
      checkboxInput("prism_input_options", "My PRISM data is a dataframe, not json results file"),
      fileInput("prism_data",
                "Upload PRISM data"),
      h5("SEMPI:"),
      fileInput("sempi_data",
                "Upload SEMPI 2.0 data"),
      h5("DEEPBGC:"),
      fileInput("deep_data",
                "Upload DeepBGC data"),
      h5("RRE-FINDER:"),
      fileInput("rre_data",
                "Upload RREFinder data"),
      # Numeric input of chromosome length of analyzed sequence
      numericInput("chr_len", "Please type chr len of an organism", value = 8773899),
      h3(id = "anti_header","Antismash data options:"),
      checkboxInput("anti_hybrid", "Visualize AntiSMASH BGC with several types as 'Hybrid'"),
      h3(id = "prism_header","PRISM data options:"),
      checkboxInput("prism_hybrid", "Visualize PRISM BGC with several types as 'Hybrid'"),
      h3(id = "sempi_header","SEMPI data options:"),
      checkboxInput("sempi_hybrid", "Visualize SEMPI BGC with several types as 'Hybrid'"),
      h3(id = "genes_on_chr","Genes on chromosome plot controls:"),
      selectInput("ref", "Choose reference data", choices = c("Antismash" = "Antismash",
                                                              "DeepBGC" = "DeepBGC",
                                                              "RRE-Finder" = "RRE-Finder",
                                                              "PRISM" = "PRISM",
                                                              "SEMPI" = "SEMPI"),
                  selected = "Antismash"),
      h3(id = "summarize","Summarize options:"),
      selectInput("group_by", "Group data by", choices = c("Antismash" = "A",
                                                              "DeepBGC" = "D",
                                                              "RRE-Finder" =  "R",
                                                              "PRISM" = "P",
                                                           "SEMPI" = "S"),
                  selected = 'A'),
      checkboxInput("count_all", "Show all BGC for the 'group by' method (+ individually annotated BGC)"),
      h3("Improve visualization:"),
      #Improve RREFinder annotated BCG visibility
      checkboxInput("rre_width", "Add thickness (+50000) to RRE results visualization (do not alter any results)"),
      h3(id="data_comparison_header","Comparison with DeepBGC plots:"),
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
        tabPanel(title = "Annotation visualization and comparison", value = 4,plotlyOutput("deep_reference_2", height = "500px"), 
                 plotlyOutput("deep_reference", height = "500px")),
        tabPanel(title = "Biocircos plot", value = 2, BioCircosOutput("biocircos", height = "1000px")),
        tabPanel(title = "Summarize interception", value = 3,plotlyOutput("barplot_rank", height = "600px"),tableOutput("group_table")),
        type = "tabs", id = "main"
      )
  )
  )
)

# Define server logic
server <- function(input, output) {
  options(shiny.maxRequestSize=100*1024^2)
  # Small function to make integers zeros
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }
  
  anti_listen <- reactive({
    list( input$anti_hybrid, input$anti_input_options)
  })
    
  # Rective vals the app is using
  # Some dataframes that are used through the app + some vectors of untercepted values
  vals <- reactiveValues(deep_data = NULL, anti_data = NULL, rre_data=NULL, prism_data=NULL, chr_len = NULL, fullness = NULL,
                         inter_a1 = NULL, inter_a2 = NULL, inter_d_ref_n = NULL,inter_d_rre=NULL,
                         inter_rre_ref_n = NULL,inter_rre_d_n = NULL, inter_a3 = NULL, inter_p_ref_n=NULL,
                         inter_d_p=NULL, inter_p_d_n = NULL, inter_p_rre = NULL, inter_rre_p_n = NULL, inter_d_rre_ID = NULL,
                         inter_d_p_ID = NULL,inter_d_ref_n_ID = NULL , biocircos_deep = NULL, deep_data_input = FALSE,
                         anti_data_input = FALSE,rre_data_input = FALSE, prism_data_input = FALSE, seg_df_ref_a = NULL,
                         seg_df_ref_d = NULL,seg_df_ref_r = NULL,seg_df_ref_p = NULL, deep_data_chromo = NULL, 
                         data_upload_count = 0, anti_type=NULL, prism_type=NULL, sempi_data = NULL, sempi_data_input= FALSE,
                         sempi_type = NULL, inter_s_ref_n = NULL, inter_a4 = NULL, inter_d_s = NULL, inter_s_d_n = NULL,
                         inter_d_s_ID = NULL, inter_s_p_n = NULL, inter_p_s = NULL,  inter_s_rre_n = NULL, inter_ree_s = NULL
                         )
  
  # Observe antismash data input and save as reactive value
  observeEvent(input$anti_data,{
    # Read data
    if (input$anti_input_options==T){
      vals$anti_data <- read.csv(input$anti_data$datapath)
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
        
        vals$anti_type <- types
        types <- sapply(types, function(x){
          if (length(unlist(x))>1){
            tmp <- str_trim(paste0(unlist(x), collapse = '', sep = " "))
            gsub(" ", "_", tmp)
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
        vals$anti_data <- anti_data
    }

    # Add chromosome column
    vals$anti_data$chromosome <-  rep("A", length(vals$anti_data$Cluster))
    # Save file
    write.csv(vals$anti_data, "anti_data.csv", row.names = F)
    vals$anti_data_input = TRUE 
    vals$data_upload_count <- vals$data_upload_count +1
  })
  
  observeEvent(input$sempi_data,{
    conn <- dbConnect(RSQLite::SQLite(), input$sempi_data$datapath)
    
    dbListTables(conn)
    
    data <- dbGetQuery(conn, "SELECT * FROM tbl_segments")
    
    
    data <- data %>%
      filter(trackid==6)
    
    vals$sempi_type <- data$name

      types <- sapply(data$name, function(x){
        tmp <- str_trim(x)
        tmp <- gsub(", ", "", tmp)
        gsub(" ", "_", tmp)
      })

    sempi_data <- data.frame(cbind(seq(1:length(data$trackid)),data$start, data$end,as.character(types)))
    colnames(sempi_data) <- c("Cluster", "Start", "Stop", "Type")
    sempi_data$Cluster <- as.numeric(sempi_data$Cluster)
    sempi_data$Start <- as.numeric(sempi_data$Start)
    sempi_data$Stop <- as.numeric(sempi_data$Stop)
    vals$sempi_data <- sempi_data
    # Add chromosome info column
    vals$sempi_data$chromosome <-  rep("S", length(vals$sempi_data$Cluster))
    # Add ID column (same as Cluster)
    vals$sempi_data$ID <- vals$sempi_data$Cluster
    # Save file
    write.csv(vals$sempi_data, "sempi_data.csv", row.names = F)
    vals$sempi_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
  })
  
  # Observe PRISM data input and save in reactive dataframe
  observeEvent(input$prism_data,{
    # Read data
    if (input$prism_input_options == T){
      vals$prism_data <- read.csv(input$prism_data$datapath)
    } else{
      data <- fromJSON(file = input$prism_data$datapath)
      
      
      types <- sapply(data$prism_results$clusters, function(x){
        tolower(x$type)
      })
      
      vals$prism_type <- types
      types <- sapply(types, function(x){
        if (length(unlist(x))>1){
          tmp <- str_trim(paste0(unlist(x), collapse = '', sep = " "))
          gsub(" ", "_", tmp)
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
      
      
      prism_data <- data.frame(cbind(start, end, types))
      prism_data <- prism_data %>%
        transmute(Cluster=as.numeric(rownames(prism_data)), Start=as.numeric(start), Stop = as.numeric(end), Type = types)
      vals$prism_data <- prism_data
      
    }
    
    # Add chromosome info column
    vals$prism_data$chromosome <-  rep("P", length(vals$prism_data$Cluster))
    # Add ID column (same as Cluster)
    vals$prism_data$ID <- vals$prism$Cluster
    # Save file
    write.csv(vals$prism_data, "prism_data.csv", row.names = F)
    vals$prism_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
  })
  
# Read and clean DeepBGC data
  observeEvent(input$deep_data, {
    # Read data
    vals$deep_data <- read.delim(input$deep_data$datapath)
    # Add chromosome info column
    vals$deep_data$chromosome <-  rep("D", length(vals$deep_data$bgc_candidate_id))
    # Add ID column as number seuquence of dataframe length
    vals$deep_data$ID <- seq(1:length(vals$deep_data$bgc_candidate_id))
    write.csv(vals$deep_data, "deep_data.csv", row.names = F)
    vals$deep_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
  })
  
  # Read RREFinder data
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
    write.csv(vals$rre_data, "rre_data.csv", row.names = F)
    vals$rre_data_input = TRUE
    vals$data_upload_count <- vals$data_upload_count +1
  })
  # Observe input of chromosome length
  observeEvent(input$chr_len,{
    vals$chr_len <- input$chr_len
  })
  
  observeEvent(vals$rre_data_input, {
    if (vals$rre_data_input == T){
      showElement(selector = "#rre_width")
    } else{
      hideElement(selector = "#rre_width")
    }
  })
  
  observeEvent(vals$deep_data_input,{
    if (vals$deep_data_input == T){
      showElement(selector = "#ref_comparison")
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
      showTab(inputId = "main", target = "1")
    } else{
      hideElement(selector = "#ref_comparison")
      hideElement(selector = "#score_type")
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
      hideTab(inputId = "main", target = "1")
    }
  })
  
  observeEvent(vals$data_upload_count, {
    if (vals$data_upload_count <2){
      hideTab("main", "2")
      hideTab("main", "3")
      hideElement(selector = "#summarize")
      hideElement(selector = "#group_by")
      hideElement(selector = "#count_all")
    }else{
      showTab("main", "2")
      showTab("main", "3")
      showElement(selector = "#summarize")
      showElement(selector = "#group_by")
      showElement(selector = "#count_all")
    }
    if (vals$data_upload_count <1){
      hideTab("main", "4")
      hideElement(selector = "#genes_on_chr")
      hideElement(selector = "#ref")
    }else{
      showTab("main", "4")
      showElement(selector = "#genes_on_chr")
      showElement(selector = "#ref")
    }
  })
  
  observeEvent(input$anti_hybrid, {
    if (input$anti_hybrid==T){
    types <- sapply(vals$anti_type, function(x){
      if (length(unlist(x))>1){
        "Hybrid"
      } else{
        x
      }
    })
    types <- sapply(types, function(x){
      if (length(unlist(x))>1){
        tmp <- str_trim(paste0(unlist(x), collapse = '', sep = " "))
        gsub(" ", "_", tmp)
      }else{
        x
      }
    })
    vals$anti_data$Type <- types
    } else {
      types <- sapply(vals$anti_type, function(x){
        if (length(unlist(x))>1){
          tmp <- str_trim(paste0(unlist(x), collapse = '', sep = " "))
          gsub(" ", "_", tmp)
        }else{
          x
        }
      })
    vals$anti_data$Type <- types
  }
  })
  
  observeEvent(vals$anti_data_input,{
    if ((vals$anti_data_input == T) & (input$anti_input_options==F)){
      showElement(selector = "#anti_header")
      showElement(selector = "#anti_hybrid")
    } else{
      hideElement(selector = "#anti_header")
      hideElement(selector = "#anti_hybrid")
    }
  })
  
  observeEvent(vals$prism_data_input,{
    if ((vals$prism_data_input == T) & (input$prism_input_options==F)){
      showElement(selector = "#prism_header")
      showElement(selector = "#prism_hybrid")
    } else{
      hideElement(selector = "#prism_header")
      hideElement(selector = "#prism_hybrid")
    }
  })
  
  observeEvent(input$prism_hybrid, {
    if (input$prism_hybrid==T){
    types <- sapply(vals$prism_type, function(x){
      if (length(unlist(x))>1){
        "Hybrid"
      } else{
        x
      }
    })
    types <- sapply(types, function(x){
      if (length(unlist(x))>1){
        tmp <- str_trim(paste0(unlist(x), collapse = '', sep = " "))
        gsub(" ", "_", tmp)
      }else{
        x
      }
    })
    vals$prism_data$Type <- types
    } else{
      types <- sapply(vals$prism_type, function(x){
        if (length(unlist(x))>1){
          tmp <- str_trim(paste0(unlist(x), collapse = '', sep = " "))
          gsub(" ", "_", tmp)
        }else{
          x
        }
      })
      vals$prism_data$Type <- types
    }
  })
  
  observeEvent(input$sempi_hybrid, {
    if (input$sempi_hybrid==T){
      types <- str_split(vals$sempi_type, " ")
      types <- str_split(types, "-")
      types <- sapply(types, function(x){
        if (length(unlist(x))>1){
          "Hybrid"
        } else{
          x
        }
      })
      
      vals$sempi_data$Type <- types
    } else {
      types <- sapply(vals$sempi_type, function(x){
        tmp <- str_trim(x)
        tmp <- gsub(", ", "", tmp)
        gsub(" ", "_", tmp)
      })
    vals$sempi_data$Type <- types
  }
  })
  
  
  #Render output plots

  # Render barplot
  output$deep_barplot <- renderPlot({
    # Require deepBGC data and sntismash data to visualize plot
    req(input$deep_data)
    
    # Create empty dataframe to populate later
    fullnes_of_annotation <- data.frame(NA, NA, NA)
    colnames(fullnes_of_annotation) <- c("Score", "Source", "Quantity")
    fullnes_of_annotation <- drop_na(fullnes_of_annotation)
    
    # Vectors of columns of score values in DeepBGC data for later subset
    score_activity <- c("antibacterial", "cytotoxic","inhibitor","antifungal")
    score_bgc <- c("deepbgc_score")
    score_cluster_type <- c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")
    
    # Subset dataframe with scores' vectors. Get max value vectors
    score_a <- apply(vals$deep_data %>% select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
    score_d <- apply(vals$deep_data %>% select(c("deepbgc_score")),1, function(x) max(x))
    score_c <- apply(vals$deep_data %>% select(c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")),1, function(x) max(x))
    
    # Decide which score to use for basic thresholds on x axis
    if (input$score_type == "Activity") {
      score_type <- score_activity
    } else if (input$score_type == "DeepBGC") {
      score_type <- score_bgc
    } else if (input$score_type == "Cluster_Type") {
      score_type <- score_cluster_type
    }
    
    # Get max value vector for chosen score
    chosen_score_vector <- apply(vals$deep_data %>% select(score_type),1, function(x) max(x))
    
    # Loop over thresholds with given step. Get the interception of antismash data with DeepBGC one at given x axis thresholds with additionsl ones
    for (dataframe_1 in seq(input$plot_start, 99, input$plot_step)){
      # Store DeepBGC dataframe in variable
      deep_inter <- vals$deep_data
      #Store max value in separate column
      deep_inter$score <- chosen_score_vector
      # Filter dataframe. Get only rows, which >= of a given thresholds. Select only start and stop of those rows as a matrix
      deep_inter <- deep_inter %>% 
        mutate(score_a = score_a, score_d = score_d, score_c = score_c) %>%
        filter(num_domains>=input$domains_filter, score>=dataframe_1/100, num_bio_domains>=input$biodomain_filter, num_proteins>=input$gene_filter,
               score_a >= as.numeric(input$score_a )/ 100, score_c >=as.numeric(input$score_c)/100 , score_d >= as.numeric(input$score_d)/100) %>%
        select(nucl_start, nucl_end) %>%
        as.matrix()
      
      # Store antismash bgc start amd atop values as matrix
      if (input$ref_comparison == 'A'){
        anti_inter <- vals$anti_data %>%
        select(Start, Stop) %>%
        as.matrix()
      } else if (input$ref_comparison == 'P'){
        anti_inter <- vals$prism_data %>%
          select(Start, Stop) %>%
          as.matrix()
      } else if (input$ref_comparison == 'S'){
    anti_inter <- vals$sempi_data %>%
      select(Start, Stop) %>%
      as.matrix()
  } 
      
      
     
      
      # Get the interception of two matrices
      interseption <- annotate(deep_inter, anti_inter) #Here we use IntervalSurgeon to get intersection list
      
      # Add zero elements to the new vector. Those are the values that have no interception and annotated only by deepBGC
      new_vect <- c()
      if (length(interseption) >0 ) {
        for (i in 1:length(interseption )){
          if (is.integer0(interseption[[i]])) {
            new_vect = c(new_vect , i)
          }
        }
      }
      
      # Get number for interception bgc
      inter_bgc <-  length(unlist(interseption))
      # Get nember of non intercepted antismash data
      if (input$ref_comparison == 'A'){
        used_antismash <-  length(vals$anti_data$Cluster)-length(unlist(interseption))
        cols <-  c("Only Antismash", "DeepBGC+Antismash", "Only DeepBGC")
        title <-  ggtitle("Comparison of Antismash and DeepBGC annotations at given score threshold")
      } else if (input$ref_comparison == 'P'){
        used_antismash <-  length(vals$prism_data$Cluster)-length(unlist(interseption))
        cols <- c("Only PRISM", "DeepBGC+PRISM", "Only DeepBGC")
        title <- ggtitle("Comparison of PRISM and DeepBGC annotations at given score threshold")
      } else if (input$ref_comparison == 'S') {
        used_antismash <-  length(vals$sempi_data$Cluster)-length(unlist(interseption))
        cols <- c("Only SEMPI", "DeepBGC+SEMPI", "Only DeepBGC")
        title <- ggtitle("Comparison of SEMPI and DeepBGC annotations at given score threshold")
      }
      
      # Number of only DeepBGC annotated clusters
      len_new <- length(new_vect)
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
    # Require DeepBGC and antismash data to begin plotting
    req(input$deep_data)
    
    # Store scores columns in vectors for DeepBGC data
    score_activity <- c("antibacterial", "cytotoxic","inhibitor","antifungal")
    score_bgc <- c("deepbgc_score")
    score_cluster_type <- c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")
    
    # Logic of what plot to draw on x axis 
    if (input$score_type == "Activity") {
      score_type <- score_activity
    } else if (input$score_type == "DeepBGC") {
      score_type <- score_bgc
    } else if (input$score_type == "Cluster_Type") {
      score_type <- score_cluster_type
    }
    
    # Store max values for chosen score in a vector
    chosen_score_vector <- apply(vals$deep_data %>% select(score_type),1, function(x) max(x))
    
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
  
  # Render interactive plot, which shows bgcs of antismash, intercepted with chosen app. Also all app bgs. On hover shows all available information
  # For antismash and PRISM data showed only ID, Start, Stop, Type
  output$deep_reference <- renderPlotly({
#if (input$ref == "Antismash") {}
#if (vals$anti_data_input == TRUE){}
#if (vals$deep_data_input == TRUE){}
#if (vals$prism_data_input == TRUE){}
#if (vals$rre_data_input == TRUE){}
#if (vals$sempi_data_input == TRUE){}
    # GENERATE DATA
    if (vals$anti_data_input == TRUE){
      # Store antismash data in local variable, with column renaming
      anti_data_chromo <-  vals$anti_data %>%
        mutate(ID = Cluster, Chr = chromosome) %>%
        dplyr::select(ID,Chr ,Start, Stop, Type)
      # Extract only Start and Stop from antismash data into matrix
      anti_inter <- vals$anti_data %>%
        select(Start, Stop) %>%
        as.matrix()
      }
    if (vals$deep_data_input == TRUE){
      # Store deepbgc max score in a vector for chosen columns
      score_a <- apply(vals$deep_data %>% select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
      score_d <- apply(vals$deep_data %>% select(c("deepbgc_score")),1, function(x) max(x))
      score_c <- apply(vals$deep_data %>% select(c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")),1, function(x) max(x))
      # Store DeepBGC data in local variable.
      deep_data_chromo <- vals$deep_data %>%
        mutate(score = apply(vals$deep_data %>%
                               dplyr::select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene),1, function(x) max(x))) 
      
      # Add Cluster_type column, which store only the max name of the cluster type
      deep_data_chromo$Cluster_type <- colnames(deep_data_chromo %>% dplyr::select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene))[apply(deep_data_chromo%>%dplyr::select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene),1, which.max) ]
      
      # Clean data, using, thesholds
      deep_data_chromo <- deep_data_chromo%>%
        # Change to "Under threshold"  Cluster_type column values, if they are under chosen theshold (no cluster type will be visualised)
        mutate(Cluster_type = ifelse(score>as.numeric(input$cluster_type)/100, Cluster_type, "Under threshold")) %>%
        # Add new columns and change product_class to Cluster_type
        mutate( product_class = Cluster_type, score_a = score_a, score_d = score_d, score_c = score_c) %>%
        filter(score_a >= as.numeric(input$score_a )/ 100, score_c >=as.numeric(input$score_c)/100 , score_d >= as.numeric(input$score_d)/100,  num_domains >= input$domains_filter,
               num_bio_domains>=input$biodomain_filter, num_proteins>=input$gene_filter)
      # Store cleaned DeepBGC data into other variable to subset only Start and Stop
      deep_inter <- deep_data_chromo
      vals$deep_data_chromo <- deep_data_chromo
      deep_inter <- deep_inter %>% 
        select(nucl_start, nucl_end) %>%
        as.matrix()
    }
    if (vals$rre_data_input == TRUE){
      # Convert numeric columns in a dataframe as a numeric
      vals$rre_data$Start <- as.numeric(vals$rre_data$Start) 
      vals$rre_data$Stop <- as.numeric(vals$rre_data$Stop)
      # Store rre data into local variable
      rre_data <- data.frame(vals$rre_data)
      # Start/Stop columns from rre data as matrix
      rre_inter <- rre_data %>%
        select(Start, Stop) %>%
        as.matrix()
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
      prism_inter <- vals$prism_data %>%
        select(Start,Stop) %>%
        as.matrix()
    }
    if (vals$sempi_data_input == TRUE){
      # Store master prism data in local variable
      sempi_data <- vals$sempi_data
      # Start/Stop columns from prism data as matrix
      sempi_inter <- vals$sempi_data %>%
        select(Start,Stop) %>%
        as.matrix()
    }
    
    get_prism_inter <- function(x,prism_inter,prism_data, letter ){
      # Get an interception of prism and smth
      interseption <- annotate(x,prism_inter)
      inter <- unlist(interseption, use.names=FALSE)
      data <- prism_data[inter,]
      # Create a dataframe with antismash data, intercepted with PRISM, with all the additional info to visualize on hover
      seg_df_3 <- data.frame(x=as.numeric(data$Start),
                             y=rep(letter, length(data$Cluster)),
                             xend=as.numeric(data$Stop),
                             yend=rep(letter, length(data$Cluster)),
                             Type = as.factor(data$Type),
                             Software = rep("PRISM", length(data$Cluster)),
                             ID = data$Cluster,
                             Start = data$Start,
                             Stop = data$Stop)
      return(seg_df_3)
    }
    get_sempi_inter <- function(x,sempi_inter,sempi_data, letter ){
      # Get an interception of sempi and smth
      interseption <- annotate(x,sempi_inter)
      inter <- unlist(interseption, use.names=FALSE)
      data <- sempi_data[inter,]
      # Create a dataframe with antismash data, intercepted with PRISM, with all the additional info to visualize on hover
      seg_df <- data.frame(x=as.numeric(data$Start),
                             y=rep(letter, length(data$Cluster)),
                             xend=as.numeric(data$Stop),
                             yend=rep(letter, length(data$Cluster)),
                             Type = as.factor(data$Type),
                             Software = rep("SEMPI", length(data$Cluster)),
                             ID = data$Cluster,
                             Start = data$Start,
                             Stop = data$Stop)
      return(seg_df)
    }
    get_anti_inter <- function(x,anti_inter,anti_data_chromo, letter ){
      # Get an interception of deepBGC and antismash data 
      interseption <- annotate( x,anti_inter)
      inter <- unlist(interseption, use.names=FALSE)
      data <- anti_data_chromo[inter,]
      seg_df_1 <- data.frame(x=as.numeric(  data$Start),
                             y=rep(letter, length(data$ID)),
                             xend=as.numeric(  data$Stop),
                             yend=rep(letter, length(data$ID)),
                             Type = as.factor(data$Type),
                             Software = rep("Antismash", length(data$ID)),
                             ID = data$ID,
                             Start = data$Start,
                             Stop = data$Stop)
      return(seg_df_1)
    }
    get_deep_inter <- function(x,deep_inter,deep_data_chromo, letter ){
      # Get an interception of deepBGC and antismash data 
      interseption <- annotate(x, deep_inter)
      inter <- unlist(interseption, use.names=FALSE)
      data <- deep_data_chromo[inter,]
      # Create a dataframe with antismash data, interce[ted with deepbgc with all the additional info to visualize on hover
      seg_df_1 <-  data.frame(x=as.numeric(  data$nucl_start),
                              y=rep(letter, length(data$ID)),
                              xend=as.numeric(  data$nucl_end),
                              yend=rep(letter, length(data$ID)),
                              Type = as.factor(data$product_class),
                              Software = rep("DeepBGC", length(data$ID)),
                              ID = data$ID,
                              Start = data$nucl_start,
                              Stop = data$nucl_end,
                              num_domains = data$num_domains,
                              deepbgc_score = data$deepbgc_score,
                              activity = data$product_activity)
      return(seg_df_1)
    }
    get_RRE_inter <- function(x,rre_inter,rre_data, letter ){
      # Get an interception of RREFinder and antismash 
      interseption <- annotate(x, rre_inter)
      inter <- unlist(interseption, use.names=FALSE)
      data <- rre_data[inter,]
      if (input$rre_width == TRUE) {
        Stop_vals_RRE_in <- as.numeric(data$Stop)+50000
      } else{
        Stop_vals_RRE_in <- as.numeric(data$Stop)
      }
      # Create a dataframe with antismash data, intercepted with RREFinder, with all the additional info to visualize on hover
      seg_df_2 <- data.frame(x=data$Start,
                             y=rep(letter, length(data$Locus_tag)),
                             xend=Stop_vals_RRE_in,
                             yend=rep(letter, length(data$Locus_tag)),
                             Type = rep("RiPP", length(data$Locus_tag)),
                             Score = data$Score,
                             Software = rep("RREFinder", length(data$Locus_tag)),
                             ID = data$Locus_tag,
                             Start = data$Start,
                             Stop = data$Stop,
                             E_value = data$E.value,
                             P_value = data$P.value,
                             RRE_start = data$RRE.start,
                             RRE_stop = data$RRE.end,
                             Probability = data$Probability)
      return(seg_df_2)
    }
    
    geom_anti <- function(data){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,                                                                                                             ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    geom_prism <- function(data){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,                                                                                                             ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    geom_deep <- function(data){
      geom_segment(data=data,aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                               ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
                                                               deepbgc_score = deepbgc_score,activity = activity ),size =3)
      }
    geom_rre <- function(data){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type, Score = Score, Software = Software,
                                                           ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value,
                                                           P_value = P_value, RRE_start = RRE_start,RRE_stop = RRE_stop, 
                                                           Probability = Probability),size = 3)
      }
    geom_sempi <- function(data){
      geom_segment(data=data, aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,                                                                                                             ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    
    
    # MAKE COMPUTATIONS
    if (vals$anti_data_input == TRUE){
      seg_ref_a <- data.frame(x=as.numeric(  anti_data_chromo$Start),
                                 y=rep("Z", length(anti_data_chromo$ID)),
                                 xend=as.numeric(  anti_data_chromo$Stop),
                                 yend=rep("Z", length(anti_data_chromo$ID)),
                                 Type = as.factor(anti_data_chromo$Type),
                                 Software = rep("Antismash", length(anti_data_chromo$ID)),
                                 ID = anti_data_chromo$ID,
                                 Start = anti_data_chromo$Start,
                                 Stop = anti_data_chromo$Stop)
      seg_ref <- seg_ref_a
      
      if (input$ref == "Antismash") {
        plot <- ggplot(anti_data_chromo, aes(x = vals$chr_len, y = Chr)) + geom_anti(seg_ref)
        if (vals$deep_data_input == TRUE){
          seg_df <-  get_deep_inter(anti_inter,deep_inter,deep_data_chromo, "Y" )
          if (dim(seg_df[1]) > 0){
          plot <- plot + geom_anti(seg_df)
          }
        }
        if (vals$rre_data_input == TRUE){
          seg_df <- get_RRE_inter(anti_inter, rre_inter, rre_data, "X")
          if (dim(seg_df[1]) > 0){
          plot <- plot + geom_rre(seg_df)
          }
        }
        if (vals$prism_data_input == TRUE){
          seg_df <- get_prism_inter(anti_inter, prism_inter, prism_data, "W")
          if (dim(seg_df[1]) > 0){
          plot <- plot +  geom_prism(seg_df)
          }
        }
        if (vals$sempi_data_input == TRUE){
          seg_df <- get_sempi_inter(anti_inter, sempi_inter, sempi_data, "V")
          if (dim(seg_df[1]) > 0){
          plot <- plot +  geom_sempi(seg_df)
          }
        }
        plot <- plot +
                   scale_y_discrete(labels = c("Z" = "Antismash","Y" = "D_vs_A", "X" = "RRE_vs_A", "W" = "P_vs_A", "V" = "S_vs_A")) +
                   theme(axis.text.y = element_text(size = 10)) +
                   ylab("")+
                   xlab("Chromosome length")+
                   ggtitle("Annotations' comparison to the reference")
        to_plot <- ggplotly(plot, tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                                   "P_value", "RRE_start","RRE_stop", "Probability" ))
        
      }
      vals$seg_df_ref_a <- seg_ref_a
    }
    if (vals$deep_data_input == TRUE){
      # Create a dataframe with all deepbgc data + the additional info to visualize on hover
      seg_ref_d <- data.frame(x=as.numeric(  deep_data_chromo$nucl_start),
                              y=rep("Z", length(deep_data_chromo$ID)),
                              xend=as.numeric(  deep_data_chromo$nucl_end),
                              yend=rep("Z", length(deep_data_chromo$ID)),
                              Type = as.factor(deep_data_chromo$product_class),
                              Software = rep("DeepBGC", length(deep_data_chromo$ID)),
                              ID = deep_data_chromo$ID,
                              Start = deep_data_chromo$nucl_start,
                              Stop = deep_data_chromo$nucl_end,
                              num_domains = deep_data_chromo$num_domains,
                              deepbgc_score = deep_data_chromo$deepbgc_score,
                              activity = deep_data_chromo$product_activity)
      seg_ref <- seg_ref_d
      
      if (input$ref == "DeepBGC") {
        plot <- ggplot(deep_data_chromo, aes(x = vals$chr_len, y = Chr)) + geom_deep(seg_ref)
      if (vals$anti_data_input == TRUE){
        seg_df <- get_anti_inter(deep_inter, anti_inter, anti_data_chromo, "X")
        if (dim(seg_df[1]) > 0){
        plot <- plot + geom_anti(seg_df)
        }
      }
      if (vals$rre_data_input == TRUE){
        seg_df <- get_RRE_inter(deep_inter, rre_inter, rre_data, "Y")
        if (dim(seg_df[1]) > 0){
        plot <- plot + geom_rre(seg_df)
        }
      }
      if (vals$prism_data_input == TRUE){
        seg_df<- get_prism_inter(deep_inter, prism_inter, prism_data, "W")
        if (dim(seg_df[1]) > 0){
        plot <- plot + geom_prism(seg_df)
        }
      }
      if (vals$sempi_data_input == TRUE){
        seg_df<- get_sempi_inter(deep_inter, sempi_inter, sempi_data, "V")
          if (dim(seg_df[1]) > 0){
          plot <- plot + geom_sempi(seg_df)
          }
        }
       to_plot <- ggplotly(plot +
                   scale_y_discrete(labels = c("Z" = "DeepBGC","X" = "A_vs_D", "Y" = "RRE_vs_D", "W" = "P_vs_D", "V" = "S_vs_D")) +
                   theme(axis.text.y = element_text(size = 10)) +
                   ylab("")+
                   xlab("Chromosome length")+
                   ggtitle("Annotations' comparison to the reference"), 
                 # What actually to visualize in tooltip
                 tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                             "P_value", "RRE_start","RRE_stop", "Probability" )
        )
      }
     
      vals$seg_df_ref_d <- data.frame(x=as.numeric(  deep_data_chromo$nucl_start),
                                      y=rep("X", length(deep_data_chromo$ID)),
                                      xend=as.numeric(  deep_data_chromo$nucl_end),
                                      yend=rep("X", length(deep_data_chromo$ID)),
                                      Type = as.factor(deep_data_chromo$product_class),
                                      Software = rep("DeepBGC", length(deep_data_chromo$ID)),
                                      ID = deep_data_chromo$ID,
                                      Start = deep_data_chromo$nucl_start,
                                      Stop = deep_data_chromo$nucl_end,
                                      num_domains = deep_data_chromo$num_domains,
                                      deepbgc_score = deep_data_chromo$deepbgc_score,
                                      activity = deep_data_chromo$product_activity)
    }
    if (vals$rre_data_input == TRUE){
      seg_ref_r <- data.frame(x=vals$rre_data$Start,
                              y=rep("Z", length(vals$rre_data$Locus_tag)),
                              xend=Stop_vals_RRE,
                              yend=rep("Z", length(vals$rre_data$Locus_tag)),
                              Type = rep("RiPP", length(vals$rre_data$Locus_tag)),
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
      seg_ref <- seg_ref_r
      if (input$ref == "RRE-Finder") {
        plot <- ggplot(rre_data, aes(x = vals$chr_len, y = Chr)) + geom_rre(seg_ref)
      if (vals$anti_data_input == TRUE){
        seg_df <- get_anti_inter(rre_inter, anti_inter, anti_data_chromo, "X")
        if (dim(seg_df[1]) > 0){
        plot <- plot + geom_anti(seg_df)
        }
      }
      if (vals$deep_data_input == TRUE){
        seg_df <- get_deep_inter(rre_inter, deep_inter, deep_data_chromo, "Y")
        if (dim(seg_df[1]) > 0){
        plot <- plot +  geom_deep(seg_df)
        }
      }
      if (vals$prism_data_input == TRUE){
        seg_df <- get_prism_inter(rre_inter, prism_inter, prism_data, "W")
        if (dim(seg_df[1]) > 0){
        plot <- plot + geom_prism(seg_df)
        }
      }
        if (vals$sempi_data_input == TRUE){
          seg_df <- get_sempi_inter(rre_inter, sempi_inter, sempi_data, "V")
          if (dim(seg_df[1]) > 0){
          plot <- plot + geom_sempi(seg_df)
          }
        }
        to_plot <- ggplotly(plot +
                   scale_y_discrete(labels = c("Z" = "RRE","X" = "A_vs_RRE", "Y" = "D_vs_RRE", "W" = "P_vs_RRE", "V" = "S_vs_RRE")) +
                   theme(axis.text.y = element_text(size = 10)) +
                   ylab("")+
                   xlab("Chromosome length")+
                   ggtitle("Annotations' comparison to the reference"), 
                 # What actually to visualize in tooltip
                 tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                             "P_value", "RRE_start","RRE_stop", "Probability" )
        )
      }
      
      vals$seg_df_ref_r <- data.frame(x=vals$rre_data$Start,
                                      y=rep("Y", length(vals$rre_data$Locus_tag)),
                                      xend=Stop_vals_RRE,
                                      yend=rep("Y", length(vals$rre_data$Locus_tag)),
                                      Type = rep("RiPP", length(vals$rre_data$Locus_tag)),
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
    }
    if (vals$prism_data_input == TRUE){
      # Create a dataframe with PRISM data with all the additional info to visualize on hover
      seg_ref_p <- data.frame(x=as.numeric(prism_data$Start),
                                y=rep("Z", length(prism_data$Cluster)),
                                xend=as.numeric(prism_data$Stop),
                                yend=rep("Z", length(prism_data$Cluster)),
                                Type = as.factor(prism_data$Type),
                                Software = rep("PRISM", length(prism_data$Cluster)),
                                ID = prism_data$Cluster,
                                Start = prism_data$Start,
                                Stop = prism_data$Stop)
      seg_ref <- seg_ref_p 
      
      if (input$ref == "PRISM") {
      plot <- ggplot(prism_data, aes(x = vals$chr_len, y = Chr)) + geom_prism(seg_ref)
      if (vals$anti_data_input == TRUE){
        seg_df <- get_anti_inter(prism_inter, anti_inter, anti_data_chromo, "X")
        if (dim(seg_df[1]) > 0){
        plot <- plot + geom_anti(seg_df)
        }
      }
      if (vals$deep_data_input == TRUE){
        seg_df <- get_deep_inter(prism_inter, deep_inter, deep_data_chromo, "Y")
        if (dim(seg_df[1]) > 0){
        plot <- plot + geom_deep(seg_df)
        }
      }
      if (vals$rre_data_input == TRUE){
        seg_df <- get_RRE_inter(prism_inter, rre_inter, rre_data, "W")
        if (dim(seg_df[1]) > 0){
        plot <- plot +geom_rre(seg_df)
        }
      }
      if (vals$sempi_data_input == TRUE){
        seg_df <- get_sempi_inter(prism_inter, sempi_inter, sempi_data, "V")
        if (dim(seg_df[1]) > 0){
        plot <- plot + geom_sempi(seg_df)
        }
      }
        # Create a dataframe with PRISM data with all the additional info to visualize on hover
        to_plot <- ggplotly(plot +
                   scale_y_discrete(labels = c("Z" = "PRISM","X" = "A_vs_P", "Y" = "D_vs_P", "W" = "RRE_vs_P", "V" = "S_vs_P")) +
                   theme(axis.text.y = element_text(size = 10)) +
                   ylab("")+
                   xlab("Chromosome length")+
                   ggtitle("Annotations' comparison to the reference"), 
                 # What actually to visualize in tooltip
                 tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                             "P_value", "RRE_start","RRE_stop", "Probability" )
                 )
      }
      vals$seg_df_ref_p <- data.frame(x=as.numeric(prism_data$Start),
                                      y=rep("W", length(prism_data$Cluster)),
                                      xend=as.numeric(prism_data$Stop),
                                      yend=rep("W", length(prism_data$Cluster)),
                                      Type = as.factor(prism_data$Type),
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
                              Software = rep("SEMPI", length(sempi_data$Cluster)),
                              ID = sempi_data$Cluster,
                              Start = sempi_data$Start,
                              Stop = sempi_data$Stop)
      seg_ref <- seg_ref_s
      
      if (input$ref == "SEMPI") {
        plot <- ggplot(sempi_data, aes(x = vals$chr_len, y = Chr)) + geom_sempi(seg_ref)
        if (vals$anti_data_input == TRUE){
          seg_df <- get_anti_inter(sempi_inter, anti_inter, anti_data_chromo, "X")
          if (dim(seg_df[1]) > 0){
          plot <- plot + geom_anti(seg_df)
          }
        }
        if (vals$deep_data_input == TRUE){
          seg_df <- get_deep_inter(sempi_inter, deep_inter, deep_data_chromo, "Y")
          if (dim(seg_df[1]) > 0){
          plot <- plot + geom_deep(seg_df)
          }
        }
        if (vals$rre_data_input == TRUE){
          seg_df <- get_RRE_inter(sempi_inter, rre_inter, rre_data, "W")
          if (dim(seg_df[1]) > 0){
          plot <- plot + geom_deep(seg_df)
          }
        }
        if (vals$prism_data_input == TRUE){
          seg_df <- get_prism_inter(sempi_inter, prism_inter, prism_data, "V")
          if (dim(seg_df[1]) > 0){
          plot <- plot + geom_prism(seg_df)
          }
        }
        # Create a dataframe with PRISM data with all the additional info to visualize on hover
        to_plot <- ggplotly(plot +
                              scale_y_discrete(labels = c("Z" = "SEMPI","X" = "A_vs_S", "Y" = "D_vs_S", "W" = "RRE_vs_S", "V" = "P_vs_S")) +
                              theme(axis.text.y = element_text(size = 10)) +
                              ylab("")+
                              xlab("Chromosome length")+
                              ggtitle("Annotations' comparison to the reference"), 
                            # What actually to visualize in tooltip
                            tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                                        "P_value", "RRE_start","RRE_stop", "Probability" )
        )
      }
      vals$seg_df_ref_s <- data.frame(x=as.numeric(sempi_data$Start),
                                      y=rep("V", length(sempi_data$Cluster)),
                                      xend=as.numeric(sempi_data$Stop),
                                      yend=rep("V", length(sempi_data$Cluster)),
                                      Type = as.factor(sempi_data$Type),
                                      Software = rep("SEMPI", length(sempi_data$Cluster)),
                                      ID = sempi_data$Cluster,
                                      Start = sempi_data$Start,
                                      Stop = sempi_data$Stop)
    }

    to_plot
  })
  
  output$deep_reference_2 <- renderPlotly({
    if (vals$anti_data_input == TRUE){
      data <-  vals$anti_data %>%
        mutate(ID = Cluster, Chr = chromosome) %>%
        dplyr::select(ID,Chr ,Start, Stop, Type)
    }
    if (vals$deep_data_input == TRUE){
      data <- vals$deep_data_chromo
    }
    if (vals$rre_data_input == TRUE){
      data <- data.frame(vals$rre_data)
    }
    if (vals$prism_data_input == TRUE){
      data <- vals$prism_data
    }
    if (vals$sempi_data_input == TRUE){
      data <- vals$sempi_data
    }
    
    plot <- ggplot(data, aes(x = vals$chr_len, y = Chr))
    if (vals$anti_data_input == TRUE){
      plot <- plot + 
        geom_segment(data=vals$seg_df_ref_a, aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,                                                                                                             ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
    }
    if (vals$deep_data_input == TRUE){
      plot <- plot +
        geom_segment(data=vals$seg_df_ref_d,aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                      ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
                                      deepbgc_score = deepbgc_score,activity = activity ),size =3)
    }
    if (vals$rre_data_input == TRUE){
      plot <- plot + geom_segment(data=vals$seg_df_ref_r, aes(x, y, xend=xend, yend=yend, color = Type, Score = Score, Software = Software,
                                                    ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value,
                                                    P_value = P_value, RRE_start = RRE_start,RRE_stop = RRE_stop, 
                                                    Probability = Probability),size = 3)
    }
    if (vals$prism_data_input == TRUE){
      plot <- plot + geom_segment(data=vals$seg_df_ref_p, aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                    ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
      
      
    }
    if (vals$sempi_data_input == TRUE){
      plot <- plot + geom_segment(data=vals$seg_df_ref_s, aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                              ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3)
      
      
    }
    to_plot <- ggplotly(plot +
                          scale_y_discrete(labels = c("Z" = "Antismash","X" = "DeepBGC", "Y" = "RRE", "W" = "PRISM", "V" = "SEMPI")) +
                          theme(axis.text.y = element_text(size = 10)) +
                          ylab("")+
                          xlab("Chromosome length")+
                          ggtitle("All annotations"), 
                        # What actually to visualize in tooltip
                        tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                                    "P_value", "RRE_start","RRE_stop", "Probability" )
    )
    to_plot
  })
  
  # Render Biocircos Plot for all-vs-all comparison
  output$biocircos <- renderBioCircos({
    #BioCircos!
    Biocircos_chromosomes <- list()
    arcs_chromosomes <- c()
    arcs_begin <- c()
    arcs_end <- c()
    arc_labels <- c()
    
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
    }
    
    #DEEPBGC
    if (vals$deep_data_input == TRUE){
      # Get vector of max values from chosen columns from deepbgc data
      score_a <- apply(vals$deep_data %>% select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
      score_d <- apply(vals$deep_data %>% select(c("deepbgc_score")),1, function(x) max(x))
      score_c <- apply(vals$deep_data %>% select(c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")),1, function(x) max(x))
      deep_data_chromo <- vals$deep_data %>%
        mutate(score = apply(vals$deep_data %>%
                               select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene),1, function(x) max(x))) 
      # Cluster_type column. Here extract colnames, and assign max value to a new column
      deep_data_chromo$Cluster_type <- colnames(deep_data_chromo %>% select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene))[apply(deep_data_chromo%>%select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene),1, which.max) ]
      # If max score is under threshold, print "Under threshold"
      deep_data_chromo <- deep_data_chromo%>%
        mutate(Cluster_type = ifelse(score>as.numeric(input$cluster_type)/100, Cluster_type, "Under threshold"))
      #Finally store deepbgc data in plotting variable. Do final scores processing 
      biocircos_deep <- deep_data_chromo%>%
        mutate( product_class = Cluster_type, score_a = score_a, score_d = score_d, score_c = score_c) %>%
        filter(score_a >= as.numeric(input$score_a )/ 100, score_c >=as.numeric(input$score_c)/100 , 
               score_d >= as.numeric(input$score_d)/100,  num_domains >= input$domains_filter,
               num_bio_domains>=input$biodomain_filter, num_proteins>=input$gene_filter)
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
    }
    
    #RRE-FINDER
    if (vals$rre_data_input == TRUE){
      biocircos_rre <- vals$rre_data
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
    }



    # Add to tracklist. Then it can be populated with links
    tracklist <- BioCircosArcTrack('myArcTrack', arcs_chromosomes, arcs_begin, arcs_end, 
                                   minRadius = 0.90, maxRadius = 0.97, labels = arc_labels)
    # Function to get interception between two matrices. Returns a list of two elements - IDs from first matrix and 
    # from second one. IDs are duplicated, if intercepted more than one time
    get_interception <- function(x,y) {
      interseption <- annotate(x, y)
      inter_x <- unlist(interseption, use.names=FALSE)
      inter_tmp <- which(sapply(interseption,length )!=0)
      inter_y <- c()
      if (length(inter_tmp) != 0) {
        tmp <- sapply(interseption,length)
        for (i in seq(1:length(tmp[which(tmp != 0)]))) {
          inter_y <- c(inter_y,rep(inter_tmp[i],tmp[which(tmp != 0)][i]))
        }}
      return(list(inter_x, inter_y))
    }
    
    #    TEMPLATE TO COPY    
    #    # Add link start. Just populate certain chromosome name times the lenght of interception 
    #    chromosomes_start <- c(chromosomes_start, )
    #    # Add link end. Just populate second output from the vectors, used above. 
    #    chromosomes_end <- c(chromosomes_end, )
    #    # Add links start positions as a start from dataframe. This vector is for chromosome start
    #    link_pos_start <- as.numeric(c(link_pos_start, ))
    #    # Add links start positions as a start from dataframe. For chromosome start variable
    #   link_pos_start_1 <- as.numeric(c(link_pos_start_1, ))
    #   # Add links start position for a chromosome stop variable
    #    link_pos_end <- as.numeric(c(link_pos_end, ))
    #    # Add links start position for a chromosome stop position
    #    link_pos_end_2 <- as.numeric(c(link_pos_end_2, )) 
    #    label_1 <- c(label_1, ) 
    #    label_2 <- c(label_2, )
    
    
    chromosomes_start <- c()
    chromosomes_end <- c()
    link_pos_start <- c()
    link_pos_start_1 <- c()
    link_pos_end <- c()
    link_pos_end_2 <- c()
    label_1 <- c()
    label_2 <- c()
    # REVERSE THE ORDER, ACCORDING TO THE QUANTITY OF THE LINKS FOR _INTER COMPUTATION?
    
    # ANTISMASH
    if (vals$anti_data_input == TRUE){
      anti_inter <- biocircos_anti %>%
        select(Start, Stop) %>%
        as.matrix()
    }
    
    # PRISM
    if (vals$prism_data_input == TRUE){
      # Store PRISM data start/stop as matrix  for further interception calculation
      prism_inter <- biocircos_prism %>%
        select(Start, Stop) %>%
        as.matrix()
    }
    
    #SEMPI
    if (vals$sempi_data_input == TRUE){
      # Store SEMPI data start/stop as matrix  for further interception calculation
      sempi_inter <- biocircos_sempi %>%
        select(Start, Stop) %>%
        as.matrix()
    }
    #RRE-FINDER
    if (vals$rre_data_input == TRUE){
      # Store RREFinder data start/stop as matrix  for futher interception calculation
      rre_inter <- biocircos_rre%>%
        select(Start, Stop) %>%
        as.matrix()
    }
    
    #DEEPBGC
    if (vals$deep_data_input == TRUE){
      # Store deepbgc data start/stop as matrix  for futher interception calculation
      deep_inter <- biocircos_deep %>% 
        select(nucl_start, nucl_end) %>%
        as.matrix()
    }
    
    #CALCULATIONS
    # -----------------------------------------
    
    # ANTISMASH
    if (vals$anti_data_input == TRUE){
      # Get interception of antismash with PRISM
      if (vals$prism_data_input == TRUE){
        inter_a1_t<- get_interception(prism_inter, anti_inter)
        inter_p_ref_n <- unlist(inter_a1_t[2])
        inter_a3 <- unlist(inter_a1_t[1])
        #Store values as reactive ones in order to use later
        vals$inter_p_ref_n <- inter_p_ref_n
        vals$inter_a3 <- inter_a3
        # Add link start. Just populate certain chromosome name times the lenght of interception 
        chromosomes_start <- c(chromosomes_start, rep("Antismash",length(c(inter_a3))))
        # Add link end. Just populate second output from the vectors, used above. 
        chromosomes_end <- c(chromosomes_end, rep("PRISM", length(inter_p_ref_n)) )
        # Add links start positions as a start from dataframe. This vector is for chromosome start
        link_pos_start <- as.numeric(c(link_pos_start, biocircos_anti$Start[c(inter_a3)]))
        # Add links start positions as a start from dataframe. For chromosome start variable
        link_pos_start_1 <- as.numeric(c(link_pos_start_1, biocircos_anti$Stop[c(inter_a3)]))
        # Add links start position for a chromosome stop variable
        link_pos_end <- as.numeric(c(link_pos_end, biocircos_prism$Start[inter_p_ref_n]))
        # Add links start position for a chromosome stop position
        link_pos_end_2 <- as.numeric(c(link_pos_end_2, biocircos_prism$Stop[inter_p_ref_n] ))
        label_1 <- c(label_1, sapply(inter_a3, function(x){x = paste("Antismash:", x, ",",biocircos_anti$Type[x])}))
        label_2 <- c(label_2, sapply(inter_p_ref_n, function(x){x = paste("PRISM:", x, ",", biocircos_prism$Type[x])}))
      }
      # Get interception of antismash with deepbgc
      if (vals$deep_data_input == TRUE){
        inter_a1_t<- get_interception(deep_inter,anti_inter )
        inter_d_ref_n <- unlist(inter_a1_t[2])
        inter_a1 <- unlist(inter_a1_t[1])
        # Store values as reactive ones in order to use later
        vals$inter_a1 <- inter_a1
        vals$inter_d_ref_n <- inter_d_ref_n
        # Add link start. Just populate certain chromosome name times the lenght of interception 
        chromosomes_start <- c(chromosomes_start, rep("Antismash",length(c(inter_a1))) )
        # Add link end. Just populate second output from the vectors, used above. 
        chromosomes_end <- c(chromosomes_end,rep("DeepBGC", length(inter_d_ref_n)) )
        # Add links start positions as a start from dataframe. This vector is for chromosome start
        link_pos_start <- as.numeric(c(link_pos_start, biocircos_anti$Start[c(inter_a1)]))
        # Add links start positions as a start from dataframe. For chromosome start variable
        link_pos_start_1 <- as.numeric(c(link_pos_start_1, biocircos_anti$Stop[c(inter_a1)]))
        # Add links start position for a chromosome stop variable
        link_pos_end <- as.numeric(c(link_pos_end, biocircos_deep$nucl_start[inter_d_ref_n]))
        # Add links start position for a chromosome stop position
        link_pos_end_2 <- as.numeric(c(link_pos_end_2, biocircos_deep$nucl_end[inter_d_ref_n]))
        # Add labels
        label_1 <- c(label_1, sapply(inter_a1, function(x){x = paste("Antismash:", x, ",",biocircos_anti$Type[x])}))
        label_2 <- c(label_2, sapply(inter_d_ref_n, function(x){x = paste("DeepBGC:", x,",", biocircos_deep$product_class[x])}))
        # Safe used local variables to the reactive ones
        vals$inter_d_ref_n_ID <- biocircos_deep$ID[inter_d_ref_n]
      } 
      # Get interception of antismash with RREFinder
      if (vals$rre_data_input == TRUE){
        inter_a1_t<- get_interception(rre_inter, anti_inter )
        inter_rre_ref_n <- unlist(inter_a1_t[2])
        inter_a2 <- unlist(inter_a1_t[1])
        #Store values as reactive ones in order to use later
        vals$inter_a2 <- inter_a2
        vals$inter_rre_ref_n <- inter_rre_ref_n
        # Add link start. Just populate certain chromosome name times the lenght of interception 
        chromosomes_start <- c(chromosomes_start, rep("Antismash",length(c(inter_a2))))
        # Add link end. Just populate second output from the vectors, used above. 
        chromosomes_end <- c(chromosomes_end,rep("RRE", length(inter_rre_ref_n)) )
        # Add links start positions as a start from dataframe. This vector is for chromosome start
        link_pos_start <- as.numeric(c(link_pos_start, biocircos_anti$Start[c(inter_a2)]))
        # Add links start positions as a start from dataframe. For chromosome start variable
        link_pos_start_1 <- as.numeric(c(link_pos_start_1, biocircos_anti$Stop[c(inter_a2)]))
        # Add links start position for a chromosome stop variable
        link_pos_end <- as.numeric(c(link_pos_end, biocircos_rre$Start[inter_rre_ref_n] ))
        # Add links start position for a chromosome stop position
        link_pos_end_2 <- as.numeric(c(link_pos_end_2, biocircos_rre$Start[inter_rre_ref_n]+50000))
        label_1 <- c(label_1, sapply(inter_a2, function(x){x = paste("Antismash:", x, ",",biocircos_anti$Type[x])}))
        label_2 <- c(label_2, sapply(inter_rre_ref_n, function(x){x = paste("RRE:", x, ",", "RiPP")}))
      }
      # Get interception of antismash with SEMPI
      if (vals$sempi_data_input == TRUE){
        inter_a1_t<- get_interception(sempi_inter, anti_inter)
        inter_s_ref_n <- unlist(inter_a1_t[2])
        inter_a4 <- unlist(inter_a1_t[1])
        #Store values as reactive ones in order to use later
        vals$inter_s_ref_n <- inter_s_ref_n
        vals$inter_a4 <- inter_a4
        # Add link start. Just populate certain chromosome name times the lenght of interception 
        chromosomes_start <- c(chromosomes_start, rep("Antismash",length(c(inter_a4))))
        # Add link end. Just populate second output from the vectors, used above. 
        chromosomes_end <- c(chromosomes_end, rep("SEMPI", length(inter_s_ref_n)) )
        # Add links start positions as a start from dataframe. This vector is for chromosome start
        link_pos_start <- as.numeric(c(link_pos_start, biocircos_anti$Start[c(inter_a4)]))
        # Add links start positions as a start from dataframe. For chromosome start variable
        link_pos_start_1 <- as.numeric(c(link_pos_start_1, biocircos_anti$Stop[c(inter_a4)]))
        # Add links start position for a chromosome stop variable
        link_pos_end <- as.numeric(c(link_pos_end, biocircos_sempi$Start[inter_s_ref_n]))
        # Add links start position for a chromosome stop position
        link_pos_end_2 <- as.numeric(c(link_pos_end_2, biocircos_sempi$Stop[inter_s_ref_n] ))
        label_1 <- c(label_1, sapply(inter_a4, function(x){x = paste("Antismash:", x, ",",biocircos_anti$Type[x])}))
        label_2 <- c(label_2, sapply(inter_s_ref_n, function(x){x = paste("SEMPI:", x, ",", biocircos_sempi$Type[x])}))
      }
      # Write csvs with locally used variables
      write.csv(biocircos_anti, "antismash_biocircos.csv", row.names = F)
    }
    
    # DEEPBGC 
    if (vals$deep_data_input == TRUE){
      # Get interception of DeepBGC with rrefinder
      if (vals$rre_data_input == TRUE){
        inter_a1_t<- get_interception(rre_inter, deep_inter )
        inter_rre_d_n <- unlist(inter_a1_t[2])
        inter_d_rre <- unlist(inter_a1_t[1])
        #Store values as reactive ones in order to use later
        vals$inter_d_rre <- inter_d_rre
        vals$inter_rre_d_n <- inter_rre_d_n
        # Add link start. Just populate certain chromosome name times the lenght of interception 
        chromosomes_start <- c(chromosomes_start, rep("DeepBGC",length(c(inter_d_rre))) )
        # Add link end. Just populate second output from the vectors, used above. 
        chromosomes_end <- c(chromosomes_end, rep("RRE", length(inter_rre_d_n)))
        # Add links start positions as a start from dataframe. This vector is for chromosome start
        link_pos_start <- as.numeric(c(link_pos_start, biocircos_deep$nucl_start[c(inter_d_rre)]))
        # Add links start positions as a start from dataframe. For chromosome start variable
        link_pos_start_1 <- as.numeric(c(link_pos_start_1,biocircos_deep$nucl_end[c(inter_d_rre)]))
        # Add links start position for a chromosome stop variable
        link_pos_end <- as.numeric(c(link_pos_end, biocircos_rre$Start[inter_rre_d_n]))
        # Add links start position for a chromosome stop position
        link_pos_end_2 <- as.numeric(c(link_pos_end_2, biocircos_rre$Start[inter_rre_d_n]+50000 ))
        label_1 <- c(label_1, sapply(c(inter_d_rre), function(x){x = paste("DeepBGC:", x, ",",biocircos_deep$product_class[x])}))
        label_2 <- c(label_2, sapply(inter_rre_d_n, function(x){x = paste("RRE:", x, ",", "RiPP")}))
        # Safe used local variables to the reactive ones
        vals$inter_d_rre_ID <- biocircos_deep$ID[inter_d_rre]
      }
      # Get interception of DeepBGC with PRISM
      if (vals$prism_data_input == TRUE){
        inter_a1_t<- get_interception(prism_inter, deep_inter)
        inter_p_d_n <- unlist(inter_a1_t[2])
        inter_d_p <- unlist(inter_a1_t[1])
        #Store values as reactive ones in order to use later
        vals$inter_d_p <- inter_d_p
        vals$inter_p_d_n <- inter_p_d_n
        # Add link start. Just populate certain chromosome name times the lenght of interception 
        chromosomes_start <- c(chromosomes_start, rep("DeepBGC",length(c(inter_d_p))))
        # Add link end. Just populate second output from the vectors, used above. 
        chromosomes_end <- c(chromosomes_end, rep("PRISM", length(inter_p_d_n)))
        # Add links start positions as a start from dataframe. This vector is for chromosome start
        link_pos_start <- as.numeric(c(link_pos_start, biocircos_deep$nucl_start[c(inter_d_p)]))
        # Add links start positions as a start from dataframe. For chromosome start variable
        link_pos_start_1 <- as.numeric(c(link_pos_start_1, biocircos_deep$nucl_end[c(inter_d_p)]))
        # Add links start position for a chromosome stop variable
        link_pos_end <- as.numeric(c(link_pos_end,biocircos_prism$Start[inter_p_d_n] ))
        # Add links start position for a chromosome stop position
        link_pos_end_2 <- as.numeric(c(link_pos_end_2,biocircos_prism$Stop[inter_p_d_n] ))
        label_1 <- c(label_1, sapply(c(inter_d_p), function(x){x = paste("DeepBGC:", x, ",",biocircos_deep$product_class[x])}))
        label_2 <- c(label_2, sapply(inter_p_d_n, function(x){x = paste("PRISM:", x, ",", biocircos_prism$Type[x])}))
        # Safe used local variables to the reactive ones
        vals$inter_d_p_ID <- biocircos_deep$ID[inter_d_p]
      }
      # Get interception of DeepBGC with SEMPI
      if (vals$sempi_data_input == TRUE){
        inter_a1_t<- get_interception(sempi_inter, deep_inter)
        inter_s_d_n <- unlist(inter_a1_t[2])
        inter_d_s <- unlist(inter_a1_t[1])
        #Store values as reactive ones in order to use later
        vals$inter_d_s <- inter_d_s
        vals$inter_s_d_n <- inter_s_d_n
        # Add link start. Just populate certain chromosome name times the lenght of interception 
        chromosomes_start <- c(chromosomes_start, rep("DeepBGC",length(c(inter_d_s))))
        # Add link end. Just populate second output from the vectors, used above. 
        chromosomes_end <- c(chromosomes_end, rep("SEMPI", length(inter_s_d_n)))
        # Add links start positions as a start from dataframe. This vector is for chromosome start
        link_pos_start <- as.numeric(c(link_pos_start, biocircos_deep$nucl_start[c(inter_d_s)]))
        # Add links start positions as a start from dataframe. For chromosome start variable
        link_pos_start_1 <- as.numeric(c(link_pos_start_1, biocircos_deep$nucl_end[c(inter_d_s)]))
        # Add links start position for a chromosome stop variable
        link_pos_end <- as.numeric(c(link_pos_end,biocircos_sempi$Start[inter_s_d_n] ))
        # Add links start position for a chromosome stop position
        link_pos_end_2 <- as.numeric(c(link_pos_end_2,biocircos_sempi$Stop[inter_s_d_n] ))
        label_1 <- c(label_1, sapply(c(inter_d_s), function(x){x = paste("DeepBGC:", x, ",",biocircos_deep$product_class[x])}))
        label_2 <- c(label_2, sapply(inter_s_d_n, function(x){x = paste("SEMPI:", x, ",", biocircos_sempi$Type[x])}))
        # Safe used local variables to the reactive ones
        vals$inter_d_s_ID <- biocircos_deep$ID[inter_d_s]
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
        inter_a1_t<- get_interception(rre_inter, prism_inter)
        inter_rre_p_n <- unlist(inter_a1_t[2])
        inter_p_rre <- unlist(inter_a1_t[1])
        #Store values as reactive ones in order to use later
        vals$inter_p_rre <- inter_p_rre
        vals$inter_rre_p_n <- inter_rre_p_n
        # Add link start. Just populate certain chromosome name times the lenght of interception 
        chromosomes_start <- c(chromosomes_start, rep("PRISM", length(inter_p_rre)))
        # Add link end. Just populate second output from the vectors, used above. 
        chromosomes_end <- c(chromosomes_end, rep("RRE", length(inter_rre_p_n)))
        # Add links start positions as a start from dataframe. This vector is for chromosome start
        link_pos_start <- as.numeric(c(link_pos_start, biocircos_prism$Start[inter_p_rre] ))
        # Add links start positions as a start from dataframe. For chromosome start variable
        link_pos_start_1 <- as.numeric(c(link_pos_start_1, biocircos_prism$Stop[inter_p_rre] ))
        # Add links start position for a chromosome stop variable
        link_pos_end <- as.numeric(c(link_pos_end,  biocircos_rre$Start[inter_rre_p_n]))
        # Add links start position for a chromosome stop position
        link_pos_end_2 <- as.numeric(c(link_pos_end_2,biocircos_rre$Start[inter_rre_p_n]+50000 ))
        label_1 <- c(label_1, sapply(inter_p_rre, function(x){x = paste("PRISM:", x, ",",biocircos_prism$Type[x])})) 
        label_2 <- c(label_2, sapply(inter_rre_p_n, function(x){x = paste("RRE", x, ",", "RiPP")}))
      }
      # Get interception of PRISM with SEMPI
      if (vals$sempi_data_input == TRUE){
        inter_a1_t<- get_interception(sempi_inter, prism_inter)
        inter_s_p_n <- unlist(inter_a1_t[2])
        inter_p_s <- unlist(inter_a1_t[1])
        #Store values as reactive ones in order to use later
        vals$inter_p_s <- inter_p_s
        vals$inter_s_p_n <- inter_s_p_n
        # Add link start. Just populate certain chromosome name times the lenght of interception 
        chromosomes_start <- c(chromosomes_start, rep("PRISM", length(inter_p_s)))
        # Add link end. Just populate second output from the vectors, used above. 
        chromosomes_end <- c(chromosomes_end, rep("SEMPI", length(inter_s_p_n)))
        # Add links start positions as a start from dataframe. This vector is for chromosome start
        link_pos_start <- as.numeric(c(link_pos_start, biocircos_prism$Start[inter_p_s] ))
        # Add links start positions as a start from dataframe. For chromosome start variable
        link_pos_start_1 <- as.numeric(c(link_pos_start_1, biocircos_prism$Stop[inter_p_s] ))
        # Add links start position for a chromosome stop variable
        link_pos_end <- as.numeric(c(link_pos_end,  biocircos_sempi$Start[inter_s_p_n]))
        # Add links start position for a chromosome stop position
        link_pos_end_2 <- as.numeric(c(link_pos_end_2,biocircos_sempi$Stop[inter_s_p_n] ))
        label_1 <- c(label_1, sapply(inter_p_s, function(x){x = paste("PRISM:", x, ",",biocircos_prism$Type[x])})) 
        label_2 <- c(label_2, sapply(inter_s_p_n, function(x){x = paste("SEMPI", x, ",", biocircos_sempi$Type[x])}))
      }
      # Write csvs with locally used variables
      write.csv(biocircos_prism, "prism_biocircos.csv", row.names = F)
    }
    
    # RRE-FINDER 
    if (vals$rre_data_input == TRUE){
      # Get interception of RRE with SEMPI
      if (vals$sempi_data_input == TRUE){
        inter_a1_t<- get_interception(sempi_inter, rre_inter)
        inter_s_rre_n <- unlist(inter_a1_t[2])
        inter_rre_s <- unlist(inter_a1_t[1])
        #Store values as reactive ones in order to use later
        vals$inter_rre_s <- inter_rre_s
        vals$inter_s_rre_n <- inter_s_rre_n
        # Add link start. Just populate certain chromosome name times the lenght of interception 
        chromosomes_start <- c(chromosomes_start, rep("RRE", length(inter_rre_s)))
        # Add link end. Just populate second output from the vectors, used above. 
        chromosomes_end <- c(chromosomes_end, rep("SEMPI", length(inter_s_rre_n)))
        # Add links start positions as a start from dataframe. This vector is for chromosome start
        link_pos_start <- as.numeric(c(link_pos_start, biocircos_rre$Start[inter_rre_s] ))
        # Add links start positions as a start from dataframe. For chromosome start variable
        link_pos_start_1 <- as.numeric(c(link_pos_start_1, biocircos_rre$Stop[inter_rre_s] + 50000 ))
        # Add links start position for a chromosome stop variable
        link_pos_end <- as.numeric(c(link_pos_end,  biocircos_sempi$Start[inter_s_rre_n]))
        # Add links start position for a chromosome stop position
        link_pos_end_2 <- as.numeric(c(link_pos_end_2,biocircos_sempi$Stop[inter_s_rre_n]))
        label_1 <- c(label_1, sapply(inter_rre_s, function(x){x = paste("RRE:", x, ",","RiPP")})) 
        label_2 <- c(label_2, sapply(inter_s_rre_n, function(x){x = paste("SEMPI", x, ",", biocircos_sempi$Type[x])}))
      }
      # Write csvs with locally used variables
      write.csv(biocircos_rre, "rre_biocircos.csv", row.names = F)
    }
    #SEMPI (NO VALUE)
    if (vals$sempi_data_input == TRUE){
      # Write csvs with locally used variables
      write.csv(biocircos_sempi, "sempi_biocircos.csv", row.names = F)
    }
    
 
    
   
    # Combine labels with mapply to one list
    link_labels <- mapply(function(x,y)  paste(x, y, sep = " | "), label_1, label_2 )
    # Add links and labels to the track list for subsequent visualization 
    tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', chromosomes_start, link_pos_start, 
                                               link_pos_start_1, chromosomes_end, link_pos_end, 
                                               link_pos_end_2, maxRadius = 0.85, labels = link_labels,
                                               displayLabel = FALSE)

    # Plot BioCircos
    BioCircos(tracklist, genome = Biocircos_chromosomes, genomeTicksScale = 1e+6)
  })
  
  # Render barplot with number count of interception for BGC IDs
  output$barplot_rank <- renderPlotly({
    antismash_count <-  NULL
    prism_count <- NULL
    deep_count <- NULL
    rre_count <- NULL
    sempi_count <- NULL
    
    if (vals$anti_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      antismash_count <- count(as.factor(c(vals$inter_a1, vals$inter_a2, vals$inter_a3, vals$inter_a4)))
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
    if (vals$deep_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      deep_count <- count(as.factor(c(vals$inter_d_ref_n_ID, vals$inter_d_rre_ID, vals$inter_d_p_ID, vals$inter_d_s_ID)))
      # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
      deep_anot <- vals$biocircos_deep[vals$biocircos_deep$ID %in% as.numeric(levels(deep_count$x)),]
      # Add prefices to the ID to plot for a barplot.  
      deep_count$x <- sapply(deep_count$x, function(x) paste("D: ", x))
      # Add label column to the dataframe, from which we will plot  
      deep_count$label <- rep("DeepBGC", length(deep_count$x))
      # Add type to the dataframe, from which we would plot (from annotation dataframe)  
      deep_count$Type <- deep_anot$product_class
      # Add Start positions (to visualize on hover)
      deep_count$Start <- deep_anot$nucl_start
      # Add Stop positions (to visualize on hover)
      deep_count$Stop <- deep_anot$nucl_end
    }
    if (vals$rre_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      rre_count <- count(as.factor(c(vals$inter_rre_ref_n, vals$inter_rre_d_n, vals$inter_rre_p_n, vals$inter_rre_s)))
      # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
      rre_anot <- vals$rre_data[vals$rre_data$ID %in% as.numeric(levels(rre_count$x)),]
      # Add prefices to the ID to plot for a barplot.  
      rre_count$x <- sapply(rre_count$x, function(x) paste("RRE: ", x))
      # Add label column to the dataframe, from which we will plot  
      rre_count$label <- rep("RRE", length(rre_count$x))
      # Add type to the dataframe, from which we would plot (from annotation dataframe)  
      rre_count$Type <- rep("RiPP", length(rre_anot$Sequence))
      # Add Start positions (to visualize on hover)
      rre_count$Start <- rre_anot$Start
      # Add Stop positions (to visualize on hover)
      rre_count$Stop <- rre_anot$Stop
    }
    if (vals$prism_data_input == TRUE){
      # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
      prism_count <- count(as.factor(c(vals$inter_p_ref_n,vals$inter_p_d_n,vals$inter_p_rre, vals$inter_p_s )))
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
      sempi_count <- count(as.factor(c(vals$inter_s_ref_n,vals$inter_s_d_n,vals$inter_s_rre_n, vals$inter_s_p_n )))
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

      
    # Integrate all those dataframe to the master one 
    ranking_data <- rbind(antismash_count,prism_count, deep_count,rre_count, sempi_count)
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
    if (vals$anti_data_input == TRUE){
      if ((input$group_by=="A") & (input$count_all == T)){
        df_f_a <- data.frame(seq(1:length(vals$anti_data$chromosome)))
        colnames(df_f_a) <- c("A")
      } else{
        df_f_a <- data.frame(A=NA)
      }
    }
    if (vals$deep_data_input == TRUE){
      if ((input$group_by=="D") & (input$count_all == T)){
        df_f_d <- data.frame(seq(1:length(vals$deep_data_chromo$ID)))
        colnames(df_f_d) <- c("D")
      } else{
        df_f_d <- data.frame(D=NA)
      }
    }
    df_f_r <- data.frame(R=NA)
    df_f_s <- data.frame(S=NA)
    if (vals$rre_data_input == TRUE){
      if ((input$group_by=="R") & (input$count_all == T)){
        df_f_r <- data.frame(seq(1:length(vals$rre_data$ID)))
        colnames(df_f_r) <- c("R")
      } else{
        df_f_r <- data.frame(R=NA)
      }
    }
    if (vals$prism_data_input == TRUE){
        if ((input$group_by=="P") & (input$count_all == T)){
          df_f_p <- data.frame(seq(1:length(vals$prism_data$Cluster)))
          colnames(df_f_p) <- c("P")
        } else{
          df_f_p <- data.frame(P=NA)
        }
    }
    if (vals$sempi_data_input == TRUE){
      if ((input$group_by=="S") & (input$count_all == T)){
        df_f_s <- data.frame(seq(1:length(vals$sempi_data$Cluster)))
        colnames(df_f_s) <- c("S")
      } else{
        df_f_s <- data.frame(S=NA)
      }
    }
    
    df_a <- data.frame(A=NA, D=NA, P=NA, R=NA, S=NA)
    if (vals$anti_data_input == TRUE){
        df_tmp <- data.frame(A=NA, D=NA)
        df_tmp1 <- data.frame(A=NA, P=NA)
        df_tmp2 <- data.frame(A=NA, R=NA)
        df_tmp3 <- data.frame(A=NA, S=NA)
        if (!is.null(vals$inter_d_ref_n)){
          df_tmp <- data.frame(cbind(vals$inter_a1,vals$inter_d_ref_n))
          colnames(df_tmp) <- c("A", "D")
        }
        if (!is.null(vals$inter_p_ref_n)){
          df_tmp1 <- data.frame(cbind(vals$inter_a3,vals$inter_p_ref_n))
          colnames(df_tmp1) <- c("A", "P")
        }
        
        if(!is.null(vals$inter_rre_ref_n)){
          df_tmp2 <- data.frame(cbind(vals$inter_a2,vals$inter_rre_ref_n))
        colnames(df_tmp2) <- c("A", "R")
        }
        if(!is.null(vals$inter_s_ref_n)){
          df_tmp3 <- data.frame(cbind(vals$inter_a4,vals$inter_s_ref_n))
          colnames(df_tmp3) <- c("A", "S")
        }
        
        
     df_a <- merge(df_f_a, df_tmp, all =T)
     df_a <- merge(df_a, df_tmp1, all = T)
     df_a <- merge(df_a, df_tmp2, all = T )
     df_a <- merge(df_a, df_tmp3, all = T )
     
    }
    
    df_d <- data.frame(D=NA, P=NA, R=NA, S=NA)
    if (vals$deep_data_input == TRUE){
      df_tmp <- data_frame(D=NA, R=NA)
      df_tmp1 <- data_frame(D=NA, P=NA)
      df_tmp2 <- data_frame(D=NA, S=NA)
      if (!is.null(vals$inter_rre_d_n)){
        df_tmp <- data.frame(cbind(vals$inter_d_rre, vals$inter_rre_d_n ))
        colnames(df_tmp) <- c("D", "R")
      }
      if (!is.null(vals$inter_p_d_n)){
        df_tmp1 <- data.frame(cbind(vals$inter_d_p ,vals$inter_p_d_n ))
        colnames(df_tmp1) <- c("D", "P")
      }
      if (!is.null(vals$inter_s_d_n)){
        df_tmp2 <- data.frame(cbind(vals$inter_d_s ,vals$inter_s_d_n ))
        colnames(df_tmp2) <- c("D", "S")
      }
      df_d <- merge(df_f_d, df_tmp, by.x="D", by.y="D", all =T)
      df_d <- merge(df_d, df_tmp1, by.x="D", by.y="D", all =T)
      df_d <- merge(df_d, df_tmp2, by.x="D", by.y="D", all =T)
    }
    
    df_p <- data_frame(P=NA, R=NA, S=NA)
    if (vals$prism_data_input == TRUE){
      df_p <- data_frame(P=NA, R=NA, S=NA)
      df_tmp <- data_frame(P=NA, R=NA)
      df_tmp2 <- data_frame(P=NA, S=NA)
      if (!is.null(vals$inter_rre_p_n)){
      df_tmp<- data.frame(cbind(vals$inter_p_rre,vals$inter_rre_p_n ))
      colnames(df_tmp) <- c("P", "R")
      }
      if (!is.null(vals$inter_s_p_n)){
        df_tmp2<- data.frame(cbind(vals$inter_p_s,vals$inter_s_p_n ))
        colnames(df_tmp2) <- c("P", "S")
      }
      df_p <- merge(df_f_p, df_tmp, by.x="P", by.y="P", all =T)
      df_p <- merge(df_p, df_tmp2, by.x="P", by.y="P", all =T)
    }
    
    df_r <- data_frame( R=NA, S=NA)
    if (vals$rre_data_input == TRUE){
      df_tmp <- data_frame(R=NA, S=NA)
      if (!is.null(vals$inter_s_rre_n)){
        df_tmp<- data.frame(cbind(vals$inter_rre_s, vals$inter_s_rre_n ))
        colnames(df_tmp) <- c("R", "S")
      }
      df_r <- merge(df_f_r, df_tmp, by.x="R", by.y="R", all =T)
    }
    if (vals$sempi_data_input == TRUE){
    }
    df_1 <- merge(df_d,df_a, all=T)
    df_2 <- merge(df_1, df_p, all=T)
    df_3 <- merge(df_2, df_r, all=T)
    df_fin <- merge(df_3, df_f_s, all=T)
    if (input$group_by=="A"){
      data <- df_fin %>% group_by(A) %>% summarise(D=paste(D, collapse=","),
                                                              R=paste(R, collapse=","),
                                                              P=paste(P, collapse=","),
                                                   S=paste(S, collapse=","))
    }else if (input$group_by=="P"){
      data <- df_fin %>% group_by(P) %>% summarise(D=paste(D, collapse=","),
                                                              R=paste(R, collapse=","),
                                                              A=paste(A, collapse=","),
                                                   S=paste(S, collapse=","))
    }else if (input$group_by=="R"){
      data <- df_fin %>% group_by(R) %>% summarise(D=paste(D, collapse=","),
                                                              A=paste(A, collapse=","),
                                                              P=paste(P, collapse=","),
                                                   S=paste(S, collapse=","))
    }else if (input$group_by=="D"){
      data <- df_fin %>% group_by(D) %>% summarise(A=paste(A, collapse=","),
                                                              R=paste(R, collapse=","),
                                                              P=paste(P, collapse=","),
                                                   S=paste(S, collapse=","))
    } else if (input$group_by=="S"){
      data <- df_fin %>% group_by(S) %>% summarise(A=paste(A, collapse=","),
                                                   R=paste(R, collapse=","),
                                                   P=paste(P, collapse=","),
                                                   D=paste(D, collapse=","))
    }
    
    if (vals$anti_data_input != TRUE){
      data <- data %>%
        select(-A)
    }
    if (vals$deep_data_input != TRUE){
      data <- data %>%
        select(-D)
    }
    if (vals$rre_data_input != TRUE){
      data <- data %>%
        select(-R)
    }
    if (vals$prism_data_input != TRUE){
      data <- data %>%
        select(-P)
    }
    if (vals$sempi_data_input != TRUE){
      data <- data %>%
        select(-S)
    }
    
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
      }
    }
    #create the zip file from flst vector
    zip(file,  flst) },
  contentType = "application/zip" )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
