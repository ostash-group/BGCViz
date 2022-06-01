#' gecco_plots UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_gecco_plots_ui <- function(id){
  ns <- NS(id)
  tagList(
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
            shiny::plotOutput(ns("gecco_barplot"), height = "500px") %>%
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
            plotly::plotlyOutput(ns("gecco_rate"), height = "500px",)%>%
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
            shiny::selectInput(ns("ref_comparison_gecco"), "Choose data for comparison with Gecco", choices = c(""),selected = ''),
            shiny::selectInput(ns("score_type_gecco"), "Choose score type to set threshold", choices = c(
              "Average p-value" = "avg_p",
              "Cluster_type score" = "Cluster_Type"),
              selected = "avg_p"),
            shiny::sliderInput(ns("plot_step_gecco"), "Choose step for plots(barplot)", min = 1, max = 50,value = 10),
            shiny::sliderInput(ns("plot_start_gecco"), "Chose plot start point(barplot)", min = 0, max = 99, value = 0)
          ), options = list(handles="w,e"))
        )
      )
    )
  )
}
    
#' gecco_plots Server Functions
#'
#' @noRd 
mod_gecco_plots_server <- function(id, vals,score_average_gecco,score_cluster_gecco){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    # Silence R CMD note
    Start <- Stop <- Source <- Quantity <- 
      Score <- Novelty_rate <- Annotation_rate <- 
      Skip_rate <- Skip_rate <- Rates_data <- 
      Rates <- NULL
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
      #write.csv(fullnes_of_annotation, "fullness.csv", row.names = F)
      
      # Make text to show on a barplot to point on additional scores' thresholds
      annotateText=paste("Applied additional thresholds", paste("Average p-value:", as.character(score_average_gecco)),
                         paste("Cluster type score:", as.character(score_cluster_gecco)), sep = "\n")
      
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
  })
}
    
## To be copied in the UI
# mod_gecco_plots_ui("gecco_plots_ui_1")
    
## To be copied in the server
# mod_gecco_plots_server("gecco_plots_ui_1", vals,score_average_gecco,score_cluster_gecco )
