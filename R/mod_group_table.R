#' group_table UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_group_table_ui <- function(id){
  ns <- NS(id)
  tagList(
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
          shiny::checkboxInput(ns("count_all"), "Show all BGC for the 'group by' method (+ individually annotated BGC)"),
          shiny::selectInput(ns("group_by"), "Group data by", choices = c(""),  selected = ''),
          shiny::tableOutput(ns("group_table"))%>%
            shinycssloaders::withSpinner()
        ),options = list(handles="w,e"))
    )
  )
}
    
#' group_table Server Functions
#'
#' @noRd 
mod_group_table_server <- function(id,vals,data_uploads,soft_names,soft_namings,data_to_use,abbr){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    output$group_table <- shiny::renderTable({
      shiny::req(vals$data_upload_count >1)
      shiny::req(vals$need_filter == F)
      shiny::req(vals$can_plot_group_table == T)
      print(paste0("input$count_all: ", input$count_all, ", input$group_by: ", input$group_by, ", "))
      #source("src/group_table_functions.R")
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
  })
}
    
## To be copied in the UI
# mod_group_table_ui("group_table_ui_1")
    
## To be copied in the server
# mod_group_table_server("group_table_ui_1", vals=vals, data_uploads = data_uploads, soft_names = soft_names, soft_namings = soft_namings, data_to_use = data_to_use, abbr = abbr)
