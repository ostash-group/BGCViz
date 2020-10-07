#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(chromoMap)
library(tidyverse)
library(plyr)
library(IntervalSurgeon)
library(plotly)
library(BioCircos)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("BGCViz"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h3("Data upload and necesary input:"),
      fileInput("anti_data",
                "Upload antismash data"),
      fileInput("prism_data",
                "Upload PRISM data"),
      fileInput("deep_data",
                "Upload DeepBGC data"),
      fileInput("rre_data",
                "Upload RREFinder data"),
      numericInput("chr_len", "Please type chr len of an organism", value = 8773899),
      h3("Controls for DeepBGC data exploration:"),
      selectInput("score_type", "Choose score type to set threshold", choices = c("Activity score" = "Activity",
                                                                                  "Cluster_type score" = "Cluster_Type",
                                                                                  "DeepBGC score" = "DeepBGC"),
                  selected = "Activity score"),
      sliderInput("plot_step", "Choose step for plots(barplot)", min = 1, max = 50,value = 10),
      sliderInput("plot_start", "Chose plot start point(barplot)", min = 0, max = 100, value = 0),
      h3("DeepBGC data filtering:"),
      sliderInput("score_a", "Activity score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      sliderInput("score_d", "DeepBGC score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      sliderInput("score_c", "Cluster_type score threshold for DeepBGC data", min = 0, max = 100, value = 50 ),
      sliderInput("domains_filter", "Domain number threshold for DeepBGC data", min = 0, max = 100, value = 5),
      sliderInput("biodomain_filter", "Biodomain number threshold for DeepBGC data", min = 0, max = 100, value = 1),
      sliderInput("gene_filter", "Protein number threshold for DeepBGC data", min = 0, max = 100, value = 1),
      sliderInput("cluster_type","Choose threshold to assign cluster type for DeepBGC data ", min = 0, max = 100, value = 50),
      h3("Improve visualization:"),
      checkboxInput("rre_width", "Add thickness (+50000) to RRE results visualization (can alter interception results)"),
      downloadButton("download","Download currently used datasets (as for Biocircos plot)" )
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("deep_barplot",height = "500px"),
      plotlyOutput("deep_rate"),
      plotlyOutput("deep_reference", height = "500px"),
      BioCircosOutput("biocircos", height = "1000px"),
      plotlyOutput("barplot_rank", height = "600px")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }
  
  
  
  vals <- reactiveValues(deep_data = NULL, anti_data = NULL, rre_data=NULL, prism_data=NULL, chr_len = NULL, fullness = NULL,
                         inter_a1 = NULL, inter_a2 = NULL, inter_d_ref_n = NULL,inter_d_rre=NULL,
                         inter_rre_ref_n = NULL,inter_rre_d_n = NULL, inter_a3 = NULL, inter_p_ref_n=NULL,
                         inter_d_p=NULL, inter_p_d_n = NULL, inter_p_rre = NULL, inter_rre_p_n = NULL, inter_d_rre_ID = NULL,
                         inter_d_p_ID = NULL,inter_d_ref_n_ID = NULL , biocircos_deep = NULL)
  
  observeEvent(input$anti_data,{
    vals$anti_data <- read.csv(input$anti_data$datapath)
    vals$anti_data$chromosome <-  rep("A", length(vals$anti_data$Cluster))
    write.csv(vals$anti_data, "anti_data.csv", row.names = F)
  })
  
  observeEvent(input$prism_data,{
    vals$prism_data <- read.csv(input$prism_data$datapath)
    vals$prism_data$chromosome <-  rep("A", length(vals$prism_data$Cluster))
    vals$prism_data$ID <- vals$prism$Cluster
    write.csv(vals$prism_data, "prism_data.csv", row.names = F)
  })
  
  observeEvent(input$deep_data, {
    vals$deep_data <- read.delim(input$deep_data$datapath)
    vals$deep_data$chromosome <-  rep("D", length(vals$deep_data$bgc_candidate_id))
    vals$deep_data$ID <- seq(1:length(vals$deep_data$bgc_candidate_id))
    write.csv(vals$deep_data, "deep_data.csv", row.names = F)
  })
  
  observeEvent(input$rre_data, {
    vals$rre_data <- read.delim(input$rre_data$datapath)
    vals$rre_data <- vals$rre_data %>%
      separate(Gene.name, c("Sequence","Coordinates","Locus_tag"),sep = "__") %>%
      separate(Coordinates, c("Start", "Stop"),sep = "-")
    vals$rre_data$chromosome <- rep("RRE",length(vals$rre_data$Sequence))
    vals$rre_data$ID <- seq(1:length(vals$rre_data$Sequence))
    write.csv(vals$rre_data, "rre_data.csv", row.names = F)
  })
  
  observeEvent(input$chr_len,{
    vals$chr_len <- input$chr_len
  })
  
  output$deep_barplot <- renderPlot({
    req(input$anti_data)
    req(input$deep_data)
    
    fullnes_of_annotation <- data.frame(NA, NA, NA)
    colnames(fullnes_of_annotation) <- c("Score", "Source", "Quantity")
    fullnes_of_annotation <- drop_na(fullnes_of_annotation)
    
    
    score_activity <- c("antibacterial", "cytotoxic","inhibitor","antifungal")
    score_bgc <- c("deepbgc_score")
    score_cluster_type <- c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")
    
    score_a <- apply(vals$deep_data %>% select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
    score_d <- apply(vals$deep_data %>% select(c("deepbgc_score")),1, function(x) max(x))
    score_c <- apply(vals$deep_data %>% select(c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")),1, function(x) max(x))
    
    if (input$score_type == "Activity") {
      score_type <- score_activity
    } else if (input$score_type == "DeepBGC") {
      score_type <- score_bgc
    } else if (input$score_type == "Cluster_Type") {
      score_type <- score_cluster_type
    }
    
    chosen_score_vector <- apply(vals$deep_data %>% select(score_type),1, function(x) max(x))
    
    for (dataframe_1 in seq(input$plot_start, 99, input$plot_step)){
      deep_inter <- vals$deep_data
      deep_inter$score <- chosen_score_vector
      deep_inter <- deep_inter %>% 
        mutate(score_a = score_a, score_d = score_d, score_c = score_c) %>%
        filter(num_domains>=input$domains_filter, score>=dataframe_1/100, num_bio_domains>=input$biodomain_filter, num_proteins>=input$gene_filter,
               score_a >= as.numeric(input$score_a )/ 100, score_c >=as.numeric(input$score_c)/100 , score_d >= as.numeric(input$score_d)/100) %>%
        select(nucl_start, nucl_end) %>%
        as.matrix()
      
      anti_inter <- vals$anti_data %>%
        select(Start, Stop) %>%
        as.matrix()
      
      interseption <- annotate(deep_inter, anti_inter) #Here we use IntervalSurgeon to get intersection list
      
      new_vect <- c()
      if (length(interseption) >0 ) {
        for (i in 1:length(interseption )){
          if (is.integer0(interseption[[i]])) {
            new_vect = c(new_vect , i)
          }
        }
      }
      
      inter_bgc <-  length(unlist(interseption))
      used_antismash <-  length(vals$anti_data$Cluster)-length(unlist(interseption))
      len_new <- length(new_vect)
      fullnes_of_annotation_1 <- data.frame(c(rep(c(as.character(dataframe_1)),3 )), 
                                            c("Only Antismash", "DeepBGC+Antismash", "Only DeepBGC"), c(used_antismash, inter_bgc, len_new))
      colnames(fullnes_of_annotation_1) <- c("Score", "Source", "Quantity")
      fullnes_of_annotation <- rbind(fullnes_of_annotation, fullnes_of_annotation_1)
      
    }
    
    vals$fullness <- data.frame(fullnes_of_annotation)
    write.csv(fullnes_of_annotation, "fullness.csv", row.names = F)
    
    annotateText=paste("Applied additional thresholds", paste("Activity score:", as.character(input$score_a)),
                       paste("DeepBGC score:", as.character(input$score_d)),
                       paste("Cluster type score:", as.character(input$score_c)), sep = "\n")
    
    
    ggplot(fullnes_of_annotation, aes(fill=Source, y=Quantity, x=Score)) + 
      geom_bar(position="dodge", stat="identity")+
      geom_text(aes(label=Quantity), position=position_dodge(width=0.9), vjust=-0.25) +
      xlab(paste(input$score_type,"Score")) +
      ggtitle("Comparison of Antismash and DeepBGC annotations at given score threshold") +
      geom_label(aes(x=Inf,y=Inf,hjust=1,vjust=1,label=annotateText ), show.legend = F)
  })
  
  output$deep_rate <- renderPlotly({
    req(input$anti_data)
    req(input$deep_data)
    
    score_activity <- c("antibacterial", "cytotoxic","inhibitor","antifungal")
    score_bgc <- c("deepbgc_score")
    score_cluster_type <- c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")
    
    if (input$score_type == "Activity") {
      score_type <- score_activity
    } else if (input$score_type == "DeepBGC") {
      score_type <- score_bgc
    } else if (input$score_type == "Cluster_Type") {
      score_type <- score_cluster_type
    }
    
    chosen_score_vector <- apply(vals$deep_data %>% select(score_type),1, function(x) max(x))
    
    fullnes_of_annotation <- data.frame(vals$fullness)
    
    test <- fullnes_of_annotation %>%
      pivot_wider(names_from = Source, values_from = Quantity)
    
    ggplotly(test %>%
               mutate(Novelty_rate = test$`Only DeepBGC`/(test$`DeepBGC+Antismash` + test$`Only Antismash`), 
                      Annotation_rate = test$`DeepBGC+Antismash`/length(vals$anti_data$Cluster), 
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
  
  output$deep_reference <- renderPlotly({
    req(input$anti_data)
    req(input$prism_data)
    req(input$deep_data)
    req(input$rre_data)
    prism_data <- vals$prism_data
    
    score_a <- apply(vals$deep_data %>% select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
    score_d <- apply(vals$deep_data %>% select(c("deepbgc_score")),1, function(x) max(x))
    score_c <- apply(vals$deep_data %>% select(c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")),1, function(x) max(x))
    
    anti_data_chromo <-  vals$anti_data %>%
      mutate(ID = Cluster, Chr = chromosome) %>%
      dplyr::select(ID,Chr ,Start, Stop, Type)
    
    deep_data_chromo <- vals$deep_data %>%
      mutate(score = apply(vals$deep_data %>%
                             dplyr::select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene),1, function(x) max(x))) 
    
    deep_data_chromo$Cluster_type <- colnames(deep_data_chromo %>% dplyr::select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene))[apply(deep_data_chromo%>%dplyr::select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene),1, which.max) ]
    
    deep_data_chromo <- deep_data_chromo%>%
      mutate(Cluster_type = ifelse(score>as.numeric(input$cluster_type)/100, Cluster_type, "Under threshold")) %>%
      mutate( product_class = Cluster_type, score_a = score_a, score_d = score_d, score_c = score_c) %>%
      filter(score_a >= as.numeric(input$score_a )/ 100, score_c >=as.numeric(input$score_c)/100 , score_d >= as.numeric(input$score_d)/100,  num_domains >= input$domains_filter,
             num_bio_domains>=input$biodomain_filter, num_proteins>=input$gene_filter)
    
    deep_inter <- deep_data_chromo
    deep_inter <- deep_inter %>% 
      select(nucl_start, nucl_end) %>%
      as.matrix()
    
    anti_inter <- vals$anti_data %>%
      select(Start, Stop) %>%
      as.matrix()
    
    rre_data <- data.frame(vals$rre_data)
    vals$rre_data$Start <- as.numeric(vals$rre_data$Start) 
    vals$rre_data$Stop <- as.numeric(vals$rre_data$Stop)
    rre_data$Start <- as.numeric(rre_data$Start) 
    rre_data$Stop <- as.numeric(rre_data$Stop)
    rre_inter <- rre_data %>%
      select(Start, Stop) %>%
      as.matrix()
    
    prism_inter <- vals$prism_data %>%
      select(Start,Stop) %>%
      as.matrix()
    
    interseption <- annotate(deep_inter, anti_inter)
    inter <- unlist(interseption, use.names=FALSE)
    anti_data_d <- anti_data_chromo[inter,]
    
    interseption <- annotate(rre_inter, anti_inter)
    inter <- unlist(interseption, use.names=FALSE)
    anti_data_r <- anti_data_chromo[inter,]
    
    interseption <- annotate(prism_inter, anti_inter)
    inter <- unlist(interseption, use.names=FALSE)
    anti_data_p <- anti_data_chromo[inter,]
    
    interseption <- annotate(prism_inter, deep_inter)
    inter <- unlist(interseption, use.names=FALSE)
    deep_data_p <- deep_data_chromo[inter,]
    
    
    seg_df_ref <- data.frame(x=as.numeric(  anti_data_chromo$Start),
                              y=rep("Z", length(anti_data_chromo$ID)),
                              xend=as.numeric(  anti_data_chromo$Stop),
                              yend=rep("Z", length(anti_data_chromo$ID)),
                              Type = as.factor(anti_data_chromo$Type),
                              Software = rep("Antismash", length(anti_data_chromo$ID)),
                              ID = anti_data_chromo$ID,
                              Start = anti_data_chromo$Start,
                              Stop = anti_data_chromo$Stop)
    
    seg_df_anti <- data.frame(x=as.numeric(  anti_data_d$Start),
                              y=rep("Y", length(anti_data_d$ID)),
                              xend=as.numeric(  anti_data_d$Stop),
                              yend=rep("Y", length(anti_data_d$ID)),
                              Type = as.factor(anti_data_d$Type),
                              Software = rep("Antismash", length(anti_data_d$ID)),
                              ID = anti_data_d$ID,
                              Start = anti_data_d$Start,
                              Stop = anti_data_d$Stop)
    
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
    
    seg_df_anti_2 <- data.frame(x=as.numeric(  anti_data_r$Start),
                                y=rep("W", length(anti_data_r$ID)),
                                xend=as.numeric(anti_data_r$Stop),
                                yend=rep("W", length(anti_data_r$ID)),
                                Type = as.factor(anti_data_r$Type),
                                Software = rep("Antismash", length(anti_data_r$ID)),
                                ID = anti_data_r$ID,
                                Start = anti_data_r$Start,
                                Stop = anti_data_r$Stop)
    
    
    if (input$rre_width == TRUE) {
      
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
   
    seg_df_anti_3 <- data.frame(x=as.numeric( anti_data_p$Start),
                                y=rep("U", length(anti_data_p$ID)),
                                xend=as.numeric(anti_data_p$Stop),
                                yend=rep("U", length(anti_data_p$ID)),
                                Type = as.factor(anti_data_p$Type),
                                Software = rep("Antismash", length(anti_data_p$ID)),
                                ID = anti_data_p$ID,
                                Start = anti_data_p$Start,
                                Stop = anti_data_p$Stop)
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
    prism_data$ID <- prism_data$Cluster
    seg_df_prism <- data.frame(x=as.numeric(prism_data$Start),
                               y=rep("S", length(prism_data$Cluster)),
                               xend=as.numeric(prism_data$Stop),
                               yend=rep("S", length(prism_data$Cluster)),
                               Type = as.factor(prism_data$Type),
                               Software = rep("PRISM", length(prism_data$Cluster)),
                               ID = prism_data$ID,
                               Start = prism_data$Start,
                               Stop = prism_data$Stop)
    
    anti_data_chromo$Chr <- as.factor(anti_data_chromo$Chr)
    levels(anti_data_chromo$Chr) <- c("A_vs_D", "D","A_vs_RRE" ,"RRE" )
    
    
    ggplotly(ggplot(anti_data_chromo, aes(x = as.numeric(Start), y = Chr)) +
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
               scale_y_discrete(labels = c("Z" = "Antismash","Y" = "A_vs_D", "X" = "D", "W" = "A_vs_RRE","V"= "RRE",
                                           "U" = "A_vs_PRISM", "T" = "D_vs_PRISM", "S" = "PRISM")) +
               theme(axis.text.y = element_text(size = 10)) +
               ylab("")+
               xlab("Chromosome length")+
               ggtitle("Annotations' comparison to the reference"), 
             tooltip = c("Software", "ID", "Start", "Stop", "Type","num_domains",  "deepbgc_score", "activity","Score","E_value",
                         "P_value", "RRE_start","RRE_stop", "Probability" ))
  })
  
  output$biocircos <- renderBioCircos({
    req(input$anti_data)
    req(input$prism_data)
    req(input$deep_data)
    req(input$rre_data)
    
    #BioCircos!
    
    score_a <- apply(vals$deep_data %>% select(c("antibacterial", "cytotoxic","inhibitor","antifungal")),1, function(x) max(x))
    score_d <- apply(vals$deep_data %>% select(c("deepbgc_score")),1, function(x) max(x))
    score_c <- apply(vals$deep_data %>% select(c("Alkaloid", "NRP","Other","Polyketide","RiPP","Saccharide","Terpene")),1, function(x) max(x))
    #Upload and filter data
    biocircos_prism <- vals$prism_data
    biocircos_anti <- vals$anti_data
    
    deep_data_chromo <- vals$deep_data %>%
      mutate(score = apply(vals$deep_data %>%
                             select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene),1, function(x) max(x))) 
    
    deep_data_chromo$Cluster_type <- colnames(deep_data_chromo %>% select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene))[apply(deep_data_chromo%>%select(Alkaloid, NRP, Other, Polyketide, RiPP, Saccharide, Terpene),1, which.max) ]
    
    deep_data_chromo <- deep_data_chromo%>%
      mutate(Cluster_type = ifelse(score>as.numeric(input$cluster_type)/100, Cluster_type, "Under threshold"))
    
    biocircos_deep <- deep_data_chromo%>%
      mutate( product_class = Cluster_type, score_a = score_a, score_d = score_d, score_c = score_c) %>%
      filter(score_a >= as.numeric(input$score_a )/ 100, score_c >=as.numeric(input$score_c)/100 , score_d >= as.numeric(input$score_d)/100,  num_domains >= input$domains_filter, 
             num_bio_domains>=input$biodomain_filter, num_proteins>=input$gene_filter)
    
    biocircos_rre <- vals$rre_data
    biocircos_rre$Start <- as.numeric(biocircos_rre$Start)
    biocircos_rre$Stop <- as.numeric(biocircos_rre$Stop) 
    
    #Make chromosome list
    Biocircos_chromosomes <- list(
      "Antismash"  = input$chr_len,
      "PRISM" = input$chr_len,
      "DeepBGC" = input$chr_len,
      "RRE" = input$chr_len
    )
    
    
    #Add arcs
    arcs_chromosomes <- c(rep("Antismash", length(biocircos_anti$Cluster)), 
                          rep("PRISM", length(biocircos_prism$Cluster)),
                          rep("DeepBGC", length(biocircos_deep$bgc_candidate_id)), 
                          rep("RRE", length(biocircos_rre$Locus_tag)))
    
    
    arcs_begin <- c(biocircos_anti$Start,
                    biocircos_prism$Start,
                    biocircos_deep$nucl_start,
                    biocircos_rre$Start)
    
    if (input$rre_width == TRUE) {
      arcs_end <- c(biocircos_anti$Stop,
                  biocircos_prism$Stop,
                  biocircos_deep$nucl_end,
                  as.numeric(biocircos_rre$Stop)+50000)
      
    } else {
      
      arcs_end <- c(biocircos_anti$Stop,
                    biocircos_prism$Stop,
                    biocircos_deep$nucl_end,
                    as.numeric(biocircos_rre$Stop))
    }
    
    
    arc_labels <- c(biocircos_anti$Type,
                    biocircos_prism$Type,
                    biocircos_deep$product_class,
                    biocircos_rre$E.value)
    
    tracklist <- BioCircosArcTrack('myArcTrack', arcs_chromosomes, arcs_begin, arcs_end, 
                                   minRadius = 0.90, maxRadius = 0.97, labels = arc_labels)
    #Make links
    
    deep_inter <- biocircos_deep %>% 
      select(nucl_start, nucl_end) %>%
      as.matrix()
    
    anti_inter <- biocircos_anti %>%
      select(Start, Stop) %>%
      as.matrix()
    
    prism_inter <- biocircos_prism %>%
      select(Start, Stop) %>%
      as.matrix()
    
    rre_inter <- biocircos_rre%>%
      select(Start, Stop) %>%
      as.matrix() 
    
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
    
    
    inter_a1_t<- get_interception(deep_inter,anti_inter )
    inter_d_ref_n <- unlist(inter_a1_t[2])
    inter_a1 <- unlist(inter_a1_t[1])
    
    
    inter_a1_t<- get_interception(rre_inter, anti_inter )
    inter_rre_ref_n <- unlist(inter_a1_t[2])
    inter_a2 <- unlist(inter_a1_t[1])
    
    inter_a1_t<- get_interception(rre_inter, deep_inter )
    inter_rre_d_n <- unlist(inter_a1_t[2])
    inter_d_rre <- unlist(inter_a1_t[1])
    
    
    inter_a1_t<- get_interception(prism_inter, anti_inter)
    inter_p_ref_n <- unlist(inter_a1_t[2])
    inter_a3 <- unlist(inter_a1_t[1])
    
    
    inter_a1_t<- get_interception(prism_inter, deep_inter)
    inter_p_d_n <- unlist(inter_a1_t[2])
    inter_d_p <- unlist(inter_a1_t[1])
    
    
    inter_a1_t<- get_interception(prism_inter, rre_inter)
    inter_rre_p_n <- unlist(inter_a1_t[2])
    inter_p_rre <- unlist(inter_a1_t[1])
    
    
    #Store values in order to use
    vals$inter_a1 <- inter_a1
    vals$inter_d_ref_n <- inter_d_ref_n
    vals$inter_a2 <- inter_a2
    vals$inter_rre_ref_n <- inter_rre_ref_n
    vals$inter_d_rre <- inter_d_rre
    vals$inter_rre_d_n <- inter_rre_d_n
    vals$inter_a3 <- inter_a3
    vals$inter_p_ref_n <- inter_p_ref_n
    vals$inter_d_p <- inter_d_p
    vals$inter_p_d_n <- inter_p_d_n
    vals$inter_p_rre <- inter_p_rre
    vals$inter_rre_p_n <- inter_rre_p_n
    
    
    
    chromosomes_start <- c(rep("Antismash",length(c(inter_a1, inter_a2, inter_a3))),
                           rep("DeepBGC",length(c(inter_d_rre, inter_d_p))), rep("PRISM", length(inter_p_rre)))
    chromosomes_end <- c(rep("DeepBGC", length(inter_d_ref_n)),rep("RRE", length(inter_rre_ref_n)), 
                         rep("PRISM", length(inter_p_ref_n)),rep("RRE", length(inter_rre_d_n)), 
                         rep("PRISM", length(inter_p_d_n)), rep("RRE", length(inter_rre_p_n)))
    
    link_pos_start <- as.numeric(c(biocircos_anti$Start[c(inter_a1, inter_a2,inter_a3)], 
                                   biocircos_deep$nucl_start[c(inter_d_rre,inter_d_p)],
                                   biocircos_prism$Start[inter_p_rre]))
    link_pos_start_1 <- as.numeric(c(biocircos_anti$Stop[c(inter_a1, inter_a2,inter_a3)], 
                                     biocircos_deep$nucl_end[c(inter_d_rre,inter_d_p)],
                                     biocircos_prism$Stop[inter_p_rre ]))
    link_pos_end <- as.numeric(c(biocircos_deep$nucl_start[inter_d_ref_n], biocircos_rre$Start[inter_rre_ref_n],
                                 biocircos_prism$Start[inter_p_ref_n], biocircos_rre$Start[inter_rre_d_n],
                                 biocircos_prism$Start[inter_p_d_n], biocircos_rre$Start[inter_rre_p_n]))
    link_pos_end_2 <- as.numeric(c(biocircos_deep$nucl_end[inter_d_ref_n], biocircos_rre$Start[inter_rre_ref_n]+50000,
                                   biocircos_prism$Stop[inter_p_ref_n], biocircos_rre$Start[inter_rre_d_n]+50000,
                                   biocircos_prism$Stop[inter_p_d_n], biocircos_rre$Start[inter_rre_p_n]+50000))
    
    for (p in seq(1:length(c(inter_a1, inter_a2, inter_a3)))) {
      tmp_1 = sapply(c(inter_a1, inter_a2, inter_a3), function(x){x = paste("Antismash:", x, ",",biocircos_anti$Type[x])})
    }
    for (p in seq(1:length(c(inter_d_rre, inter_d_p)))) {
      tmp_2 = c(sapply(c(inter_d_rre, inter_d_p), function(x){x = paste("DeepBGC:", x, ",",biocircos_deep$product_class[x])}))
    }
    for (p in seq(1:length(inter_p_rre))) {
      tmp_3 = c(sapply(inter_p_rre, function(x){x = paste("PRISM:", x, ",",biocircos_prism$Type[x])}))
    }
    tmp_label_1 <- c(tmp_1,tmp_2, tmp_3)
    
    for (p in seq(1:length(inter_d_ref_n))) {
      tmp_4 = c(sapply(inter_d_ref_n, function(x){x = paste("DeepBGC:", x,",", biocircos_deep$product_class[x])}))
    }
    for (p in seq(1:length(inter_rre_ref_n))) {
      tmp_5 = c(sapply(inter_rre_ref_n, function(x){x = paste("RRE:", x, ",", "RiPP")}))
    }
    for (p in seq(1:length(inter_p_ref_n))) {
      tmp_6 = c(sapply(inter_p_ref_n, function(x){x = paste("PRISM:", x, ",", biocircos_prism$Type[x])}))
    }
    for (p in seq(1:length(inter_rre_d_n))) {
      tmp_7 = c(sapply(inter_rre_d_n, function(x){x = paste("RRE:", x, ",", "RiPP")}))
    }
    for (p in seq(1:length(inter_p_d_n))) {
      tmp_8 = c(sapply(inter_p_d_n, function(x){x = paste("PRISM:", x, ",", biocircos_prism$Type[x])}))
    }
    for (p in seq(1:length(inter_rre_p_n))) {
      tmp_9 = c(sapply(inter_rre_p_n, function(x){x = paste("RRE", x, ",", "RiPP")}))
    }
    
    tmp_label_3 <- c(tmp_4, tmp_5, tmp_6, tmp_7, tmp_8, tmp_9)
    
    link_labels <- mapply(function(x,y)  paste(x, y, sep = " | "), tmp_label_1, tmp_label_3 )
    
    tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', chromosomes_start, link_pos_start, 
                                               link_pos_start_1, chromosomes_end, link_pos_end, 
                                               link_pos_end_2, maxRadius = 0.85, labels = link_labels,
                                               displayLabel = FALSE)
    
    vals$inter_d_rre_ID <- biocircos_deep$ID[inter_d_rre]
    vals$inter_d_p_ID <- biocircos_deep$ID[inter_d_p]
    vals$inter_d_ref_n_ID <- biocircos_deep$ID[inter_d_ref_n]
    vals$biocircos_deep <- biocircos_deep
    
    write.csv(biocircos_anti, "antismash_biocircos.csv", row.names = F)
    write.csv(biocircos_deep, "deepbgc_biocircos.csv", row.names = F)
    write.csv(biocircos_prism, "prism_biocircos.csv", row.names = F)
    write.csv(biocircos_rre, "rre_biocircos.csv", row.names = F)
    
    BioCircos(tracklist, genome = Biocircos_chromosomes, genomeTicksScale = 1e+6)
  })
  
  output$barplot_rank <- renderPlotly({
    req(input$anti_data)
    req(input$prism_data)
    req(input$deep_data)
    req(input$rre_data)
    
    
    antismash_count <- count(as.factor(c(vals$inter_a1, vals$inter_a2, vals$inter_a3)))
    prism_count <- count(as.factor(c(vals$inter_p_ref_n,vals$inter_p_d_n,vals$inter_p_rre )))
    deep_count <- count(as.factor(c(vals$inter_d_ref_n_ID, vals$inter_d_rre_ID, vals$inter_d_p_ID)))
    rre_count <- count(as.factor(c(vals$inter_rre_ref_n, vals$inter_rre_d_n, vals$inter_rre_p_n)))
    
    
    anti_anot <- vals$anti_data[vals$anti_data$Cluster %in% as.numeric(levels(antismash_count$x)),]
    prism_anot <- vals$prism_data[vals$prism_data$Cluster %in% as.numeric(levels(prism_count$x)),]
    rre_anot <- vals$rre_data[vals$rre_data$ID %in% as.numeric(levels(rre_count$x)),]
    deep_anot <- vals$biocircos_deep[vals$biocircos_deep$ID %in% as.numeric(levels(deep_count$x)),]
    
    
    antismash_count$x <- sapply(antismash_count$x, function(x) paste("A: ", x))
    prism_count$x <- sapply(prism_count$x, function(x) paste("P: ", x))
    deep_count$x <- sapply(deep_count$x, function(x) paste("D: ", x))
    rre_count$x <- sapply(rre_count$x, function(x) paste("RRE: ", x))
    
    
    antismash_count$label <- rep("Antismash", length(antismash_count$x))
    prism_count$label <- rep("PRISM", length(prism_count$x))
    deep_count$label <- rep("DeepBGC", length(deep_count$x))
    rre_count$label <- rep("RRE", length(rre_count$x))
    
    
    antismash_count$Type <- anti_anot$Type
    prism_count$Type <- prism_anot$Type
    rre_count$Type <- rep("RiPP", length(rre_anot$Sequence))
    deep_count$Type <- deep_anot$product_class
    
    antismash_count$Start <- anti_anot$Start
    prism_count$Start <- prism_anot$Start
    rre_count$Start <- rre_anot$Start
    deep_count$Start <- deep_anot$nucl_start
    
    antismash_count$Stop <- anti_anot$Stop
    prism_count$Stop <- prism_anot$Stop
    rre_count$Stop <- rre_anot$Stop
    deep_count$Stop <- deep_anot$nucl_end
    
    ranking_data <- rbind(antismash_count,prism_count, deep_count,rre_count)
    
    colnames(ranking_data) <- c("Cluster", "Count", "Label", "Type", "Start", "Stop")
    
    ggplotly(ggplot(ranking_data, aes(x = Cluster, y = Count, Type = Type, Start = Start, Stop = Stop)) +
               geom_bar(stat = "identity", aes(fill = Label)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
                     axis.text.y = element_text(size = 14)) +
               ggtitle("Number of times cluster is annotated with other tool"),
             
             
             
             tooltip=c("Type", "Start", "Stop")  )
    
    
    
  })
  
  
  output$download <- downloadHandler(filename = function(){
    paste("datasets.zip")     
  },  
  content =  function(file){
    flst <- c()
    files_in_dir <- list.files()
    for (file_names in files_in_dir) {
      if (grepl('_biocircos.csv', file_names, fixed = TRUE)) {
        flst <- c(flst, file_names)
      }
    }
    #create the zip file
    zip(file,  flst) },
  contentType = "application/zip" )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
