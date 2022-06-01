#' deep_barplot UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_deepbgc_plots_ui <- function(id){
  ns <- NS(id)
  tagList(
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
          shiny::plotOutput(ns("deep_barplot"), height = "500px",) %>%
            shinycssloaders::withSpinner()
        ),options = list(handles="w,e"))),
    div(id = "id2",
        shinyjqui::jqui_resizable(shinydashboardPlus::box( 
          title = "DeepBGC comparison controls",
          id = "deep_comparison_controls_box",
          collapsible = TRUE,                                          
          closable = TRUE,
          shiny::selectInput(ns("ref_comparison"), "Choose data for comparison with DeepBGC", choices = c(""), selected = ''),
          # Score to use for thresholds
          shiny::selectInput(ns("score_type"), "Choose score type to set threshold", choices = c("Activity score" = "Activity",
                                                                                             "Cluster_type score" = "Cluster_Type",
                                                                                             "DeepBGC score" = "DeepBGC"),
                             selected = "Activity score"),
          # Chose step for barplot (as a threshold to draw a bar)
          shiny::sliderInput(ns("plot_step"), "Choose step for plots(barplot)", min = 1, max = 50,value = 10),
          shiny::sliderInput(ns("plot_start"), "Chose plot start point(barplot)", min = 0, max = 99, value = 0)
        ),options = list(handles="w,e")))
      )
    ),
    shiny::fluidRow(
      tags$div( id = "deep_data2",
                div(id = "id2", 
                    shinyjqui::jqui_resizable(shinydashboardPlus::box(
                      title = "DeepBGC rate",
                      id = "deep_rate_box",
                      collapsible = TRUE,                                         
                      height = "100%",
                      plotly::plotlyOutput(ns("deep_rate"), height = "500px",) %>%
                        shinycssloaders::withSpinner()
                    ),options = list(handles="w,e")))    
      ))
  )
}
    
#' deep_barplot Server Functions
#'
#' @noRd 
mod_deepbgc_plots_server <- function(id, vals, score_a, score_d, score_c){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    # Silence R CMD note
    Start <- Stop <- Source <- 
      Quantity <- Score <- Novelty_rate <- 
      Annotation_rate <- Skip_rate <- Rates_data <- 
      Rates <- NULL
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
      #write.csv(fullnes_of_annotation, "fullness.csv", row.names = F)
      
      # Make text to show on a barplot to point on additional scores' thresholds
      annotateText=paste("Applied additional thresholds", paste("Activity score:", as.character(score_a)),
                         paste("DeepBGC score:", as.character(score_d)),
                         paste("Cluster type score:", as.character(score_c)), sep = "\n")
      
      # Plot the barplot
      ggplot2::ggplot(fullnes_of_annotation, ggplot2::aes(fill=Source, y=Quantity, x=Score)) + 
        ggplot2::geom_bar(position="dodge", stat="identity")+
        ggplot2::geom_text(ggplot2::aes(label=Quantity), position=ggplot2::position_dodge(width=0.9), vjust=-0.25) +
        ggplot2::xlab(paste(input$score_type,"Score")) +
        title +
        ggplot2::geom_label(ggplot2::aes(x=Inf,y=Inf,hjust=1,vjust=1,label=annotateText ), show.legend = F)
    })
    
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
  })
}
## To be copied in the UI
# mod_deepbgc_plots_ui("deep_barplot_ui_1")
    
## To be copied in the server
# mod_deepbgc_plots_server("deep_barplot_ui_1", vals, score_a, score_d, score_c)
