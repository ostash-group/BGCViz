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

# Define UI 
ui <- fluidPage(
  
  # Application title
  titlePanel("BGCViz"),
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      # Data upload
      h3("Data upload and necesary input:"),
      fileInput("anti_data",
                "Upload antismash data"),
      fileInput("prism_data",
                "Upload PRISM data"),
      fileInput("deep_data",
                "Upload DeepBGC data"),
      fileInput("rre_data",
                "Upload RREFinder data"),
      # Numeric input of chromosome length of analyzed sequence
      numericInput("chr_len", "Please type chr len of an organism", value = 8773899),
      # Some controls for first two plots.
      h3("Controls for DeepBGC data exploration:"),
      # Score to use for thresholds
      selectInput("score_type", "Choose score type to set threshold", choices = c("Activity score" = "Activity",
                                                                                  "Cluster_type score" = "Cluster_Type",
                                                                                  "DeepBGC score" = "DeepBGC"),
                  selected = "Activity score"),
      # Chose step for barplot (as a threshold to draw a bar)
      sliderInput("plot_step", "Choose step for plots(barplot)", min = 1, max = 50,value = 10),
      sliderInput("plot_start", "Chose plot start point(barplot)", min = 0, max = 100, value = 0),
      # DeepBGC data filtering 
      h3("DeepBGC data filtering:"),
      # Different score filtering. Remain >= of set threshold
      sliderInput("score_a", "Activity score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      sliderInput("score_d", "DeepBGC score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      sliderInput("score_c", "Cluster_type score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      # Domains, biodomains and proteins filter. Remain >= of set threshold
      sliderInput("domains_filter", "Domain number threshold for DeepBGC data", min = 0, max = 100, value = 5),
      sliderInput("biodomain_filter", "Biodomain number threshold for DeepBGC data", min = 0, max = 100, value = 1),
      sliderInput("gene_filter", "Protein number threshold for DeepBGC data", min = 0, max = 100, value = 1),
      sliderInput("cluster_type","Choose threshold to assign cluster type for DeepBGC data ", min = 0, max = 100, value = 50),
      h3("Improve visualization:"),
      #Improve RREFinder annotated BCG visibility
      checkboxInput("rre_width", "Add thickness (+50000) to RRE results visualization (can alter interception results)"),
      # Donwload currently used datasets
      downloadButton("download","Download currently used datasets (as for Biocircos plot)" )
      
    ),
    
    # Show plots
    mainPanel(
      plotOutput("deep_barplot",height = "500px"),
      plotlyOutput("deep_rate"),
      plotlyOutput("deep_reference", height = "500px"),
      BioCircosOutput("biocircos", height = "1000px"),
      plotlyOutput("barplot_rank", height = "600px")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Small function to make integers zeros
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }
  
  
  # Rective vals the app is using
  # Some dataframes that are used through the app + some vectors of untercepted values
  vals <- reactiveValues(deep_data = NULL, anti_data = NULL, rre_data=NULL, prism_data=NULL, chr_len = NULL, fullness = NULL,
                         inter_a1 = NULL, inter_a2 = NULL, inter_d_ref_n = NULL,inter_d_rre=NULL,
                         inter_rre_ref_n = NULL,inter_rre_d_n = NULL, inter_a3 = NULL, inter_p_ref_n=NULL,
                         inter_d_p=NULL, inter_p_d_n = NULL, inter_p_rre = NULL, inter_rre_p_n = NULL, inter_d_rre_ID = NULL,
                         inter_d_p_ID = NULL,inter_d_ref_n_ID = NULL , biocircos_deep = NULL, deep_data_input = FALSE,
                         anti_data_input = FALSE,rre_data_input = FALSE, prism_data_input = FALSE, deep_data_biocircos = FALSE
                         )
  
  # Observe antismash data input and save as reactive value
  observeEvent(input$anti_data,{
    # Read data
    vals$anti_data <- read.csv(input$anti_data$datapath)
    # Add chromosome column
    vals$anti_data$chromosome <-  rep("A", length(vals$anti_data$Cluster))
    # Save file
    write.csv(vals$anti_data, "anti_data.csv", row.names = F)
    vals$anti_data_input = TRUE
  })
  
  # Observe PRISM data input and save in rective dataframe
  observeEvent(input$prism_data,{
    # Read data
    vals$prism_data <- read.csv(input$prism_data$datapath)
    # Add chromosome info column
    vals$prism_data$chromosome <-  rep("P", length(vals$prism_data$Cluster))
    # Add ID column (same as Cluster)
    vals$prism_data$ID <- vals$prism$Cluster
    # Save file
    write.csv(vals$prism_data, "prism_data.csv", row.names = F)
    vals$prism_data_input = TRUE
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
  })
  # Observe input of chromosome length
  observeEvent(input$chr_len,{
    vals$chr_len <- input$chr_len
  })
  
  #Render output plots

  # Render barplot
  output$deep_barplot <- renderPlot({
    # Require deepBGC data and sntismash data to visualize plot
    req(input$anti_data)
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
      anti_inter <- vals$anti_data %>%
        select(Start, Stop) %>%
        as.matrix()
      
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
      used_antismash <-  length(vals$anti_data$Cluster)-length(unlist(interseption))
      # Number of only DeepBGC annotated clusters
      len_new <- length(new_vect)
      # Combine all vectors into one dataframe
      fullnes_of_annotation_1 <- data.frame(c(rep(c(as.character(dataframe_1)),3 )), 
                                            c("Only Antismash", "DeepBGC+Antismash", "Only DeepBGC"), c(used_antismash, inter_bgc, len_new))
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
      ggtitle("Comparison of Antismash and DeepBGC annotations at given score threshold") +
      geom_label(aes(x=Inf,y=Inf,hjust=1,vjust=1,label=annotateText ), show.legend = F)
  })
  
  # Render interactive plot with plotly for rates of DeepBGC data in regards with antismash data
  output$deep_rate <- renderPlotly({
    # Require DeepBGC and antismash data to begin plotting
    req(input$anti_data)
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
    
    # Calculate rates and plot interactive plot with plotly
    ggplotly(test %>%
                # Calculate rates. Novelty is nummber of clusters annotated only by deepbgc/ all clusters annotated by antismash + (antismash + deepbgc)
               mutate(Novelty_rate = test$`Only DeepBGC`/(test$`DeepBGC+Antismash` + test$`Only Antismash`), 
                      #Annotation rate = clusters, annotated by antismash+deepBGC/ clusters annotated only by antismash (We assume that antismash annotation is full and reference)
                      Annotation_rate = test$`DeepBGC+Antismash`/length(vals$anti_data$Cluster), 
                      # Skip rate = clusters, annotated only by antismash/ all antismash clusters. Points to how much clusters DeepBGC missed
                      Skip_rate = test$`Only Antismash`/length(vals$anti_data$Cluster)) %>%
               pivot_longer(cols = c(Novelty_rate, Annotation_rate, Skip_rate), names_to = 'Rates', values_to = 'Rates_data') %>%
               ggplot(aes(x=as.numeric(Score), y=as.numeric(Rates_data), Rate = as.numeric(Rates_data))) +
               geom_line(aes(color=Rates)) +
               geom_point(aes(shape=Rates), alpha = .4, size = 3) +
               ggtitle("Rates of DeepBGC/Antismash data annotation") +
               ylab("Rate") +
               xlab(paste(input$score_type,"Score threshold")),
             tooltip = c("Rate"))
  })
  
  # Render interactive plot, which shows bgcs of antismash, intercepted with chosen app. Also all app bgs. On hover shows all available information
  # For antismash and PRISM data showed only ID, Start, Stop, Type
  output$deep_reference <- renderPlotly({
    # Require antismash, PRISM, deepbgc and RREFinder data to start
    req(input$anti_data)
    req(input$prism_data)
    req(input$deep_data)
    req(input$rre_data)

    # Store master prism data in local variable
    prism_data <- vals$prism_data
    
    # Store deepbgc max score in a vector for chosen columns
    score_a <- apply(vals$deep_data %>% select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
    score_d <- apply(vals$deep_data %>% select(c("deepbgc_score")),1, function(x) max(x))
    score_c <- apply(vals$deep_data %>% select(c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")),1, function(x) max(x))
    
    # Store antismash data in local variable, with column renaming
    anti_data_chromo <-  vals$anti_data %>%
      mutate(ID = Cluster, Chr = chromosome) %>%
      dplyr::select(ID,Chr ,Start, Stop, Type)
    
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
    deep_inter <- deep_inter %>% 
      select(nucl_start, nucl_end) %>%
      as.matrix()
    
    # Extract only Start and Stop from antismash data into matrix
    anti_inter <- vals$anti_data %>%
      select(Start, Stop) %>%
      as.matrix()
    
    # Convert numeric columns in a dataframe as a numeric
    vals$rre_data$Start <- as.numeric(vals$rre_data$Start) 
    vals$rre_data$Stop <- as.numeric(vals$rre_data$Stop)
    # Store rre data into local variable
    rre_data <- data.frame(vals$rre_data)
    # Start/Stop columns from rre data as matrix
    rre_inter <- rre_data %>%
      select(Start, Stop) %>%
      as.matrix()
    # Start/Stop columns from prism data as matrix
    prism_inter <- vals$prism_data %>%
      select(Start,Stop) %>%
      as.matrix()
    # Get an interception of deepBGC and antismash data 
    interseption <- annotate(deep_inter, anti_inter)
    inter <- unlist(interseption, use.names=FALSE)
    anti_data_d <- anti_data_chromo[inter,]
    # Get an interception of RREFinder and antismash 
    interseption <- annotate(rre_inter, anti_inter)
    inter <- unlist(interseption, use.names=FALSE)
    anti_data_r <- anti_data_chromo[inter,]
    # GEt an interception of prism and antismash
    interseption <- annotate(prism_inter, anti_inter)
    inter <- unlist(interseption, use.names=FALSE)
    anti_data_p <- anti_data_chromo[inter,]
    # Get an interception of prism and DeepBGC
    interseption <- annotate(prism_inter, deep_inter)
    inter <- unlist(interseption, use.names=FALSE)
    deep_data_p <- deep_data_chromo[inter,]

    # Create a dataframe with antismash data with all the additional info to visualize on hover    
    seg_df_ref <- data.frame(x=as.numeric(  anti_data_chromo$Start),
                              y=rep("Z", length(anti_data_chromo$ID)),
                              xend=as.numeric(  anti_data_chromo$Stop),
                              yend=rep("Z", length(anti_data_chromo$ID)),
                              Type = as.factor(anti_data_chromo$Type),
                              Software = rep("Antismash", length(anti_data_chromo$ID)),
                              ID = anti_data_chromo$ID,
                              Start = anti_data_chromo$Start,
                              Stop = anti_data_chromo$Stop)
    # Create a dataframe with antismash data, interce[ted with deepbgc with all the additional info to visualize on hover
    seg_df_anti <- data.frame(x=as.numeric(  anti_data_d$Start),
                              y=rep("Y", length(anti_data_d$ID)),
                              xend=as.numeric(  anti_data_d$Stop),
                              yend=rep("Y", length(anti_data_d$ID)),
                              Type = as.factor(anti_data_d$Type),
                              Software = rep("Antismash", length(anti_data_d$ID)),
                              ID = anti_data_d$ID,
                              Start = anti_data_d$Start,
                              Stop = anti_data_d$Stop)
    # Create a dataframe with all deepbgc data + the additional info to visualize on hover
    seg_df_deep <- data.frame(x=as.numeric(  deep_data_chromo$nucl_start),
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
    # Create a dataframe with antismash data, intercepted with RREFinder, with all the additional info to visualize on hover
    seg_df_anti_2 <- data.frame(x=as.numeric(  anti_data_r$Start),
                                y=rep("W", length(anti_data_r$ID)),
                                xend=as.numeric(anti_data_r$Stop),
                                yend=rep("W", length(anti_data_r$ID)),
                                Type = as.factor(anti_data_r$Type),
                                Software = rep("Antismash", length(anti_data_r$ID)),
                                ID = anti_data_r$ID,
                                Start = anti_data_r$Start,
                                Stop = anti_data_r$Stop)
    
    # Some logic of how wide RREFinder bgcs should be
    if (input$rre_width == TRUE) {
      # Create a dataframe with RREFinder data with all the additional info to visualize on hover + 50000 to Stop to make it look thick
      seg_df_rre <- data.frame(x=vals$rre_data$Start,
                               y=rep("V", length(vals$rre_data$Locus_tag)),
                               xend=as.numeric(vals$rre_data$Stop)+50000,
                               yend=rep("V", length(vals$rre_data$Locus_tag)),
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
      
    } else {
      # Create a dataframe with antismash data with all the additional info to visualize on hover. Plot as it is
       seg_df_rre <- data.frame(x=vals$rre_data$Start,
                             y=rep("V", length(vals$rre_data$Locus_tag)),
                             xend=as.numeric(vals$rre_data$Stop),
                             yend=rep("V", length(vals$rre_data$Locus_tag)),
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
    # Create a dataframe with antismash data, intercepted with PRISM, with all the additional info to visualize on hover
    seg_df_anti_3 <- data.frame(x=as.numeric( anti_data_p$Start),
                                y=rep("U", length(anti_data_p$ID)),
                                xend=as.numeric(anti_data_p$Stop),
                                yend=rep("U", length(anti_data_p$ID)),
                                Type = as.factor(anti_data_p$Type),
                                Software = rep("Antismash", length(anti_data_p$ID)),
                                ID = anti_data_p$ID,
                                Start = anti_data_p$Start,
                                Stop = anti_data_p$Stop)

    # Create a dataframe with deepbgc data, intercepted with PRISM, with all the additional info to visualize on hover
    seg_df_deep_2 <- data.frame(x=as.numeric(  deep_data_p$nucl_start),
                                y=rep("T", length(deep_data_p$ID)),
                                xend=as.numeric(  deep_data_p$nucl_end),
                                yend=rep("T", length(deep_data_p$ID)),
                                Type = as.factor(deep_data_p$product_class),
                                Software = rep("DeepBGC", length(deep_data_p$ID)),
                                ID = deep_data_p$ID,
                                Start = deep_data_p$nucl_start,
                                Stop = deep_data_p$nucl_end,
                                num_domains = deep_data_p$num_domains,
                                deepbgc_score = deep_data_p$deepbgc_score,
                                activity = deep_data_p$product_activity)
    # Add ID column to PRISM (Not sure, but I think that was done in the preprocessing step. Delete?)
    prism_data$ID <- prism_data$Cluster
    # Create a dataframe with PRISM data with all the additional info to visualize on hover
    seg_df_prism <- data.frame(x=as.numeric(prism_data$Start),
                               y=rep("S", length(prism_data$Cluster)),
                               xend=as.numeric(prism_data$Stop),
                               yend=rep("S", length(prism_data$Cluster)),
                               Type = as.factor(prism_data$Type),
                               Software = rep("PRISM", length(prism_data$Cluster)),
                               ID = prism_data$ID,
                               Start = prism_data$Start,
                               Stop = prism_data$Stop)
    
    # Convert chromosome info in antismash data to as factor for proper visualization
    anti_data_chromo$Chr <- as.factor(anti_data_chromo$Chr)
    #Add levels
    levels(anti_data_chromo$Chr) <- c("A_vs_D", "D","A_vs_RRE" ,"RRE" )
    
    
    ggplotly(ggplot(anti_data_chromo, aes(x = as.numeric(Start), y = Chr)) +
              # Add all the dataframe info with all the metadata we want to display
               geom_segment(data=seg_df_ref, aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                  ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3) +
               geom_segment(data=seg_df_anti, aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                  ID = ID, Start = Start, Stop = Stop, Type = Type ), size = 3) + 
               geom_segment(data=seg_df_deep, aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                  ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
                                                  deepbgc_score = deepbgc_score,activity = activity ), size = 3) +
               geom_segment(data=seg_df_anti_2,aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                   ID = ID, Start = Start, Stop = Stop, Type = Type),size =3) +
               geom_segment(data=seg_df_rre, aes(x, y, xend=xend, yend=yend, color = Type, Score = Score, Software = Software,
                                                 ID = ID, Start = Start, Stop = Stop, Type = Type, E_value = E_value,
                                                 P_value = P_value, RRE_start = RRE_start,RRE_stop = RRE_stop, 
                                                 Probability = Probability),size = 3)+
               geom_segment(data=seg_df_anti_3,aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                   ID = ID, Start = Start, Stop = Stop, Type = Type),size =3) +
               geom_segment(data=seg_df_deep_2,aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                   ID = ID, Start = Start, Stop = Stop, Type = Type, num_domains = num_domains,
                                                   deepbgc_score = deepbgc_score,activity = activity ),size =3) +
               geom_segment(data=seg_df_prism,aes(x, y, xend=xend, yend=yend, color = Type, Software = Software,
                                                  ID = ID, Start = Start, Stop = Stop, Type = Type),size =3) +
              # Actually rename chromosomes. The initial name sare here for proper "top_to_bottom" (because names are sorted alphabetically from bottom_to_top) visualization.
               scale_y_discrete(labels = c("Z" = "Antismash","Y" = "A_vs_D", "X" = "D", "W" = "A_vs_RRE","V"= "RRE",
                                           "U" = "A_vs_PRISM", "T" = "D_vs_PRISM", "S" = "PRISM")) +
               theme(axis.text.y = element_text(size = 10)) +
               ylab("")+
               xlab("Chromosome length")+
               ggtitle("Annotations' comparison to the reference"), 
               # What actually to visualize in tooltip
             tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                         "P_value", "RRE_start","RRE_stop", "Probability" ))
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
      Biocircos_chromosomes[["Antismash"]] <- input$chr_len  
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
      Biocircos_chromosomes[["DeepBGC"]] <- input$chr_len
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
      Biocircos_chromosomes[["RRE"]] <- input$chr_len
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
      Biocircos_chromosomes[["PRISM"]] <- input$chr_len
      #Add arcs. Quantity of arcs is length of dataframes
      arcs_chromosomes <- c(arcs_chromosomes, rep("PRISM", length(biocircos_prism$Cluster)))
      # Add arcs begin positions. (Start column)
      arcs_begin <- c(arcs_begin, biocircos_prism$Start )
      # Stop position of arcs.
      arcs_end <- c(arcs_end, biocircos_prism$Stop )
      # Add Arcs labels. Can add only one label...
      arc_labels <- c(arc_labels,biocircos_prism$Type )
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
      # Write csvs with locally used variables
      write.csv(biocircos_prism, "prism_biocircos.csv", row.names = F)
    }
    
    # RRE-FINDER (NO VALUE HERE)
    if (vals$rre_data_input == TRUE){
      # Write csvs with locally used variables
      write.csv(biocircos_rre, "rre_biocircos.csv", row.names = F)
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
    # Begin to plot only if all data is uploaded
    
    req(input$anti_data)
    req(input$prism_data)
    req(input$deep_data)
    req(input$rre_data)
    # Count every interception for every tool. If count is >=1, this means, that given cluster at least intercepted with 1 other from another tool annotation
    antismash_count <- count(as.factor(c(vals$inter_a1, vals$inter_a2, vals$inter_a3)))
    prism_count <- count(as.factor(c(vals$inter_p_ref_n,vals$inter_p_d_n,vals$inter_p_rre )))
    deep_count <- count(as.factor(c(vals$inter_d_ref_n_ID, vals$inter_d_rre_ID, vals$inter_d_p_ID)))
    rre_count <- count(as.factor(c(vals$inter_rre_ref_n, vals$inter_rre_d_n, vals$inter_rre_p_n)))
    # Check if ID is in dataframe and if it is - extract all information about to the local dataframe  
    anti_anot <- vals$anti_data[vals$anti_data$Cluster %in% as.numeric(levels(antismash_count$x)),]
    prism_anot <- vals$prism_data[vals$prism_data$Cluster %in% as.numeric(levels(prism_count$x)),]
    rre_anot <- vals$rre_data[vals$rre_data$ID %in% as.numeric(levels(rre_count$x)),]
    deep_anot <- vals$biocircos_deep[vals$biocircos_deep$ID %in% as.numeric(levels(deep_count$x)),]
    # Add prefices to the ID to plot for a barplot.  
    antismash_count$x <- sapply(antismash_count$x, function(x) paste("A: ", x))
    prism_count$x <- sapply(prism_count$x, function(x) paste("P: ", x))
    deep_count$x <- sapply(deep_count$x, function(x) paste("D: ", x))
    rre_count$x <- sapply(rre_count$x, function(x) paste("RRE: ", x))
    # Add label column to the dataframe, from which we will plot  
    antismash_count$label <- rep("Antismash", length(antismash_count$x))
    prism_count$label <- rep("PRISM", length(prism_count$x))
    deep_count$label <- rep("DeepBGC", length(deep_count$x))
    rre_count$label <- rep("RRE", length(rre_count$x))
    # Add type to the dataframe, from which we would plot (from annotation dataframe)  
    antismash_count$Type <- anti_anot$Type
    prism_count$Type <- prism_anot$Type
    rre_count$Type <- rep("RiPP", length(rre_anot$Sequence))
    deep_count$Type <- deep_anot$product_class
    # Add Start positions (to visualize on hover)
    antismash_count$Start <- anti_anot$Start
    prism_count$Start <- prism_anot$Start
    rre_count$Start <- rre_anot$Start
    deep_count$Start <- deep_anot$nucl_start
    # Add Stop positions (to visualize on hover)
    antismash_count$Stop <- anti_anot$Stop
    prism_count$Stop <- prism_anot$Stop
    rre_count$Stop <- rre_anot$Stop
    deep_count$Stop <- deep_anot$nucl_end
    # Integrate all those dataframe to the master one 
    ranking_data <- rbind(antismash_count,prism_count, deep_count,rre_count)
    # Fix column names in the master dataframe
    colnames(ranking_data) <- c("Cluster", "Count", "Label", "Type", "Start", "Stop")
    # Plot
    ggplotly(ggplot(ranking_data, aes(x = Cluster, y = Count, Type = Type, Start = Start, Stop = Stop)) +
               geom_bar(stat = "identity", aes(fill = Label)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
                     axis.text.y = element_text(size = 14)) +
               ggtitle("Number of times cluster is annotated with other tool"),
             
             
             
             tooltip=c("Type", "Start", "Stop")  )
    
    
    
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
