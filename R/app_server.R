#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
    # Your application server logic
    # Silence R CMD note
    Start <- Stop <- Core <- Chr <- NULL
    ## ---------------------------------------------------------------
    ##        Some lists of reactive values to listen later         -
    ## ---------------------------------------------------------------
    check_to_rename <- shiny::reactive({
        list(
            input$sempi_data, input$anti_data, input$prism_data,
            input$sempi_sco, input$anti_sco, input$prism_sco,
            input$ripp_sco, input$ripp_data
        )
    })
    biocircos_listen <- shiny::reactive({
        list(
            input$biocircos_color, vals$need_filter, input$label_color, input$label_color_class,
            input$ref_col_biocircos, vals$inters_filtered, input$prism_supp_data_input_width, vals$prism_supp_data_input,
            input$arts_width, input$sempi_width, input$rre_width, vals$anti_data, vals$sempi_data, vals$prism_data,
            vals$coloring_datatable,
            vals$ripp_data
        )
    })
    inputData <- shiny::reactive({
        list(
            vals$sempi_data_input, vals$rre_data_input, vals$anti_data_input, vals$prism_data_input,
            vals$prism_supp_data_input, vals$deep_data_input, vals$gecco_data_input, vals$arts_data_input,
            vals$ripp_data_input
        )
    })
    dynamicInput <- shiny::reactive({
        list(input$dup_choice, vals$need_filter, input$prism_supp, input$phylo_file)
    })
    deep_reference <- shiny::reactive({
        list(
            vals$inters_filtered, vals$rre_more, input$ref, input$arts_width, input$sempi_width, input$rre_width,
            input$prism_supp_data_input_width, vals$anti_data, vals$prism_data, vals$sempi_data, vals$arts_data,
            vals$ripp_data, vals$arts_tree_data
        )
    })

    to_debounce <- shiny::reactive({
        list(
            vals$cluster_type, vals$gene_filter, vals$biodomain_filter, vals$score_c, vals$score_d,
            vals$score_a, vals$score_average_gecco, vals$score_cluster_gecco, vals$domains_filter_gecco,
            vals$prot_filter_gecco
        )
    }) %>% shiny::debounce(500)

    # Some dataframes that are used through the app + some vectors of untercepted values
    vals <- shiny::reactiveValues(
        deep_data = NULL, anti_data = NULL, rre_data = NULL, prism_data = NULL, chr_len = NULL, fullness_deep = NULL,
        biocircos_deep = NULL, deep_data_input = FALSE, tracklist = NULL, chromosomes = NULL, fullness_gecco = NULL,
        anti_data_input = FALSE, rre_data_input = FALSE, prism_data_input = FALSE, seg_df_ref_a = NULL,
        seg_df_ref_d = NULL, seg_df_ref_r = NULL, seg_df_ref_p = NULL, deep_data_chromo = NULL,
        data_upload_count = 0, anti_type = NULL, prism_type = NULL, sempi_data = NULL, sempi_data_input = FALSE,
        sempi_type = NULL, biocircos_color = NULL, rename_data = NULL, group_by_data = NULL,
        rre_interact = NULL, anti_interact = NULL, prism_interact = NULL, deep_interact = NULL,
        sempi_interact = NULL, df_a = NULL, df_d = NULL, df_p = NULL, df_r = NULL, prism_supp = NULL,
        prism_json = FALSE, df_s = NULL, prism_supp_interact = NULL, known_data = NULL, dup_data = NULL,
        known_data_input = FALSE, dup_data_input = FALSE, arts_data = NULL, arts_tree_data = NULL, arts_data_input = FALSE, seg_df_ref_ar = NULL,
        df_ps = NULL, arts_interact = NULL, rre_more = FALSE, gecco_data = NULL, gecco_data_input = FALSE,
        gecco_data_filtered = NULL, seg_df_ref_g = NULL, prism_supp_data_input = FALSE, computed = NULL,
        need_filter = FALSE, filter_data = FALSE, choices = list(ref = NULL, group_by = NULL, ref_col_biocircos = NULL, ref_comparison_gecco = NULL, ref_comparison = NULL),
        renamed = NULL, renaming_notification = list(), rename_y_axis = list(), can_plot_deep_ref_2 = FALSE, can_plot_deep_ref = FALSE,
        can_plot_biocircos = FALSE, can_plot_barplot_rank = FALSE, can_plot_group_table = FALSE, prism_supp_plot = FALSE,
        ripp_data = NULL, ripp_data_input = FALSE, ripp_type = NULL, ripp_interact = NULL, seg_df_ref_ri = NULL
    )

    vals$computed <- list(
        anti = FALSE, deep = FALSE, gecco = FALSE, arts = FALSE, prism = FALSE, sempi = FALSE, prism_supp = FALSE, rre = FALSE, ripp = FALSE
    )
    # Making coloring datatable
    rename_file <- system.file("extdata", "rename.csv", package = "BGCViz")
    vals$rename_data <- utils::read.csv(rename_file)
    rename_data <- utils::read.csv(rename_file)
    coloring_datatable <- data.frame(tidyr::drop_na(data.frame(cbind(as.character(rename_data$Group_color), as.character(rename_data$Color), rename_data$Hierarchy))))
    coloring_datatable <- coloring_datatable[!apply(coloring_datatable == "", 1, all), ]
    colnames(coloring_datatable) <- c("Name", "Color", "Hierarchy")
    vals$coloring_datatable <- DT::datatable(coloring_datatable, rownames = FALSE, editable = "column", options = list(dom = "t", ordering = FALSE))
    # Variables, that holds data uploads boolean (so if data is present or not)
    data_uploads <- c(
        "anti_data_input", "sempi_data_input", "prism_data_input", "prism_supp_data_input",
        "arts_data_input", "deep_data_input", "gecco_data_input", "rre_data_input",
        "ripp_data_input"
    )
    data_uploads_inter <- c(
        "anti_data_input", "sempi_data_input", "prism_data_input", "prism_json",
        "arts_data_input", "deep_data_input", "gecco_data_input", "rre_data_input",
        "ripp_data_input"
    )
    # Universal beginings for variables, used in the app for different data
    soft_names <- c("anti", "sempi", "prism", "prism_supp", "arts", "deep", "gecco", "rre", "ripp")
    # The Namings, meaning how to label the data on the plots
    soft_namings <- c("Antismash", "SEMPI", "PRISM", "PRISM-Supp", "ARTS", "DeepBGC", "GECCO", "RRE-Finder", "RippMiner")
    # Dataframes undes vals$list, that stored the data
    data_to_use <- c("anti_data", "sempi_data", "prism_data", "prism_supp_data", "arts_data_filtered", "deep_data_filtered", "gecco_data_filtered", "rre_data","ripp_data")
    # Used in barplot on summarise tab + Annotation on chromosome plots
    abbr <- c("A", "S", "P", "P-Supp", "AR", "D", "G", "RRE", "Ripp")
    # Used for deep reference 2 plot
    soft_datafr <- c(
        "seg_df_ref_a", "seg_df_ref_s", "seg_df_ref_p", "seg_df_ref_p_s", "seg_df_ref_ar", "seg_df_ref_d",
        "seg_df_ref_g", "seg_df_ref_r", "seg_df_ref_ri"
    )

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
    ## ----------------------------------------------------------------
    ##                        Helper functions                       -
    ## ----------------------------------------------------------------
    # Need to get them to a tidyr::separate file later
    # TODO
    files_in_dir <- list.files()
    # Iterate over those files and if found "_biocircos.csv" add remove them
    for (file_names in files_in_dir) {
        if (grepl("_biocircos.csv", file_names, fixed = TRUE)) {
            file.remove(file_names)
        }
    }
    options(shiny.maxRequestSize = 100 * 1024^2)

    disable_event_logic <- function() {
        vals$can_plot_deep_ref <- FALSE
        vals$can_plot_biocircos <- FALSE
        vals$can_plot_barplot_rank <- FALSE
        vals$can_plot_group_table <- FALSE
    }
    enable_event_logic <- function() {
        vals$can_plot_deep_ref <- TRUE
        vals$can_plot_biocircos <- TRUE
        vals$can_plot_barplot_rank <- TRUE
        vals$can_plot_group_table <- TRUE
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
    process_rippminer <- function(data, example_data = FALSE) {
      if (example_data == TRUE) {
        ripp_data <- data
      } else {
        ripp_data <- read_ripp(data)
      }
      vals$ripp_type <- ripp_data$Type2
      vals$ripp_data <- ripp_data
      vals$ripp_data_input <- TRUE
      vals$data_upload_count <- vals$data_upload_count + 1
      vals$choices$ref <- c(vals$choices$ref, "RippMiner" = "RippMiner")
      vals$choices$group_by <- c(vals$choices$group_by, "RippMiner" = "RippMiner")
      vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "RippMiner" = "RippMiner")
      vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "RippMiner" = "RippMiner")
      vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "RippMiner" = "RippMiner")
      update_ui_with_data()
      disable_event_logic()
      if (vals$data_upload_count == 1) {
        shiny::updateSelectInput(session, "ref",
                                 selected = "RippMiner"
        )
        shiny::updateSelectInput(session, "group_table_ui_1-group_by",
                                 selected = "RippMiner"
        )
        shiny::updateSelectInput(session, "deep_barplot_ui_1-ref_comparison",
                                 selected = "RippMiner"
        )
        shiny::updateSelectInput(session, "ref_col_biocircos",
                                 selected = "RippMiner"
        )
        shiny::updateSelectInput(session, "gecco_plots_ui_1-ref_comparison_gecco",
                                 selected = "RippMiner"
        )
      }
    }
    
    process_antismash <- function(data, example_data = FALSE) {
        if (example_data == TRUE) {
            anti_data <- data
        } else {
            anti_data <- read_anti(data)
        }
        vals$anti_type <- anti_data$Type2
        vals$anti_data <- anti_data
        vals$anti_data_input <- TRUE
        vals$data_upload_count <- vals$data_upload_count + 1
        vals$choices$ref <- c(vals$choices$ref, "Antismash" = "Antismash")
        vals$choices$group_by <- c(vals$choices$group_by, "Antismash" = "Antismash")
        vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "Antismash" = "Antismash")
        vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "Antismash" = "Antismash")
        vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "Antismash" = "Antismash")
        update_ui_with_data()
        disable_event_logic()
        if (vals$data_upload_count == 1) {
            shiny::updateSelectInput(session, "ref",
                selected = "Antismash"
            )
            shiny::updateSelectInput(session, "group_table_ui_1-group_by",
                selected = "Antismash"
            )
            shiny::updateSelectInput(session, "deep_barplot_ui_1-ref_comparison",
                selected = "Antismash"
            )
            shiny::updateSelectInput(session, "ref_col_biocircos",
                selected = "Antismash"
            )
            shiny::updateSelectInput(session, "gecco_plots_ui_1-ref_comparison_gecco",
                selected = "Antismash"
            )
        }
    }
    process_gecco <- function(data, example_data = FALSE) {
        if (example_data == TRUE) {
            gecco_data <- data
        } else {
            gecco_data <- read_gecco(data)
        }
        vals$gecco_data <- gecco_data
        vals$gecco_data_filtered <- filter_gecco(vals$gecco_data, vals$score_cluster_gecco, vals$score_average_gecco, vals$domains_filter_gecco, vals$prot_filter_gecco)
        vals$gecco_data_input <- TRUE
        vals$data_upload_count <- vals$data_upload_count + 1
        vals$choices$ref <- c(vals$choices$ref, "GECCO" = "GECCO")
        vals$choices$group_by <- c(vals$choices$group_by, "GECCO" = "GECCO")
        vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "GECCO" = "GECCO")
        update_ui_with_data()
        disable_event_logic()
        if (vals$data_upload_count == 1) {
            shiny::updateSelectInput(session, "ref",
                selected = "GECCO"
            )
            shiny::updateSelectInput(session, "group_table_ui_1-group_by",
                selected = "GECCO"
            )
            shiny::updateSelectInput(session, "ref_col_biocircos",
                selected = "GECCO"
            )
        }
    }
    process_prism <- function(data, json = TRUE, example_data = FALSE) {
        if (example_data == TRUE) {
            prism_data <- data
            prism_data$Type <- stringr::str_trim(tolower(prism_data$Type))
            prism_data["Type2"] <- stringr::str_trim(tolower(prism_data$Type))
            vals$prism_supp_data_input <- TRUE
            vals$prism_supp <- BGCViz:::prism_supp_data
            vals$prism_supp_data <- BGCViz:::prism_supp_data
            vals$prism_supp_plot <- TRUE
            vals$prism_json <- TRUE
            shiny::updateCheckboxInput(inputId = "prism_supp", value = TRUE)
        } else {
            if (json == TRUE) {
                processed <- read_prism(data, json = TRUE)
                prism_data <- processed[[1]]
                vals$prism_supp_data_input <- TRUE
                vals$prism_supp <- processed[[2]]
                vals$prism_supp_data <- processed[[2]]
                vals$prism_json <- TRUE
                vals$prism_supp_plot <- TRUE
                shiny::updateCheckboxInput(inputId = "prism_supp", value = TRUE)
            } else {
                processed <- read_prism(data, json = FALSE)
                prism_data <- processed[[1]]
            }
        }
        vals$choices$ref <- c(vals$choices$ref, "PRISM-Supp" = "PRISM-Supp")
        vals$choices$group_by <- c(vals$choices$group_by, "PRISM-Supp" = "PRISM-Supp")
        vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "PRISM-Supp" = "PRISM-Supp")
        update_ui_with_data()
        vals$prism_data <- prism_data
        vals$prism_type <- prism_data$Type2

        # Add chromosome info column
        vals$prism_data$chromosome <- rep("P", length(vals$prism_data$Cluster))
        # Add ID column (same as Cluster)
        vals$prism_data$ID <- vals$prism_data$Cluster
        vals$prism_data_input <- TRUE
        vals$data_upload_count <- vals$data_upload_count + 1
        vals$choices$ref <- c(vals$choices$ref, "PRISM" = "PRISM")
        vals$choices$group_by <- c(vals$choices$group_by, "PRISM" = "PRISM")
        vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "PRISM" = "PRISM")
        vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "PRISM" = "PRISM")
        vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "PRISM" = "PRISM")
        update_ui_with_data()
        disable_event_logic()
        if (vals$data_upload_count == 1) {
            shiny::updateSelectInput(session, "ref",
                selected = "PRISM"
            )
            shiny::updateSelectInput(session, "group_table_ui_1-group_by",
                selected = "PRISM"
            )
            shiny::updateSelectInput(session, "deep_barplot_ui_1-ref_comparison",
                selected = "PRISM"
            )
            shiny::updateSelectInput(session, "ref_col_biocircos",
                selected = "PRISM"
            )
            shiny::updateSelectInput(session, "gecco_plots_ui_1-ref_comparison_gecco",
                selected = "PRISM"
            )
        }
    }
    process_sempi <- function(data, zip = TRUE, example_data = FALSE) {
        if (example_data == TRUE) {
            sempi_data <- data
        } else {
            if (zip == TRUE) {
                sempi_data <- read_sempi(data, zip = TRUE)
            } else {
                sempi_data <- read_sempi(data, zip = FALSE)
            }
        }
        vals$sempi_type <- sempi_data$Type2
        vals$sempi_data <- sempi_data
        # Add chromosome info column
        vals$sempi_data$chromosome <- rep("S", length(vals$sempi_data$Cluster))
        # Add ID column (same as Cluster)
        vals$sempi_data$ID <- vals$sempi_data$Cluster
        vals$sempi_data_input <- TRUE
        vals$data_upload_count <- vals$data_upload_count + 1
        vals$choices$ref <- c(vals$choices$ref, "SEMPI" = "SEMPI")
        vals$choices$group_by <- c(vals$choices$group_by, "SEMPI" = "SEMPI")
        vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "SEMPI" = "SEMPI")
        vals$choices$ref_comparison_gecco <- c(vals$choices$ref_comparison_gecco, "SEMPI" = "SEMPI")
        vals$choices$ref_comparison <- c(vals$choices$ref_comparison, "SEMPI" = "SEMPI")
        update_ui_with_data()
        disable_event_logic()
        if (vals$data_upload_count == 1) {
            shiny::updateSelectInput(session, "ref",
                selected = "SEMPI"
            )
            shiny::updateSelectInput(session, "group_table_ui_1-group_by",
                selected = "SEMPI"
            )
            shiny::updateSelectInput(session, "deep_barplot_ui_1-ref_comparison",
                selected = "SEMPI"
            )
            shiny::updateSelectInput(session, "ref_col_biocircos",
                selected = "SEMPI"
            )
            shiny::updateSelectInput(session, "gecco_plots_ui_1-ref_comparison_gecco",
                selected = "SEMPI"
            )
        }
    }
    process_arts_archive <- function(archive, zip = TRUE, example_data = FALSE) {
        if (example_data == TRUE) {
            arts_data <- BGCViz:::arts_data
        } else {
            if (zip == TRUE) {
                arts_data <- read_arts_archive(archive, zip = TRUE)
            } else {
                arts_data <- utils::read.csv(archive)
            }
        }
        vals$arts_tree_data <- arts_data
        vals$arts_data <- arts_data[,!(names(arts_data) %in% c("Trees", "TreesFiles"))]
        vals$choices$ref <- c(vals$choices$ref, "ARTS" = "ARTS")
        vals$choices$group_by <- c(vals$choices$group_by, "ARTS" = "ARTS")
        vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "ARTS" = "ARTS")
        update_ui_with_data()
        vals$data_upload_count <- vals$data_upload_count + 1
        vals$arts_data_input <- TRUE
        dup_table_id <- vals$arts_data %>%
            dplyr::filter(Core != "Not_core")
        shiny::updateSelectInput(session, "dup_choice",
            choices = c("All", paste0("ID:", dup_table_id$ID, " ,Core:", dup_table_id$Core)),
            selected = "All"
        )
        

        
        if (vals$data_upload_count == 1) {
            shiny::updateSelectInput(session, "ref",
                selected = "ARTS"
            )
            shiny::updateSelectInput(session, "group_table_ui_1-group_by",
                selected = "ARTS"
            )
            shiny::updateSelectInput(session, "ref_col_biocircos",
                selected = "ARTS"
            )
        }
    }
    process_deep <- function(data, example_data = FALSE) {
        if (example_data == TRUE) {
            deep_data <- data
        } else {
            deep_data <- read_deep(data)
        }
        vals$deep_data <- deep_data
        vals$deep_data$chromosome <- rep("D", length(vals$deep_data$bgc_candidate_id))
        vals$deep_data$Start <- vals$deep_data$nucl_start
        vals$deep_data$Stop <- vals$deep_data$nucl_end
        # Add ID column as number seuquence of dataframe length
        vals$deep_data$ID <- seq(1:length(vals$deep_data$bgc_candidate_id))
        vals$deep_data$Cluster <- vals$deep_data$ID
        vals$deep_data_input <- TRUE
        vals$data_upload_count <- vals$data_upload_count + 1
        vals$deep_data_filtered <- filter_deepbgc(vals$deep_data, vals$cluster_type, vals$score_a, vals$score_c, vals$score_d, vals$domains_filter, vals$biodomain_filter, vals$gene_filter)
        vals$choices$ref <- c(vals$choices$ref, "DeepBGC" = "DeepBGC")
        vals$choices$group_by <- c(vals$choices$group_by, "DeepBGC" = "DeepBGC")
        vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "DeepBGC" = "DeepBGC")
        update_ui_with_data()
        disable_event_logic()
        if (vals$data_upload_count == 1) {
            shiny::updateSelectInput(session, "ref",
                selected = "DeepBGC"
            )
            shiny::updateSelectInput(session, "group_table_ui_1-group_by",
                selected = "DeepBGC"
            )
            shiny::updateSelectInput(session, "ref_col_biocircos",
                choices = "DeepBGC",
                selected = "DeepBGC"
            )
        }
    }
    process_rre <- function(data, example_data = FALSE) {
        if (example_data == TRUE) {
            rre_data <- data
        } else {
            rre_data <- read_rre(data)
        }
        vals$rre_data <- rre_data
        # write.csv(vals$rre_data, "rre_data.csv", row.names = FALSE)

        vals$rre_data_input <- TRUE
        vals$data_upload_count <- vals$data_upload_count + 1
        vals$choices$ref <- c(vals$choices$ref, "RRE-Finder" = "RRE-Finder")
        vals$choices$group_by <- c(vals$choices$group_by, "RRE-Finder" = "RRE-Finder")
        vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "RRE-Finder" = "RRE-Finder")
        update_ui_with_data()
        disable_event_logic()
        if (vals$data_upload_count == 1) {
            shiny::updateSelectInput(session, "ref",
                selected = "RRE-Finder"
            )
            shiny::updateSelectInput(session, "group_table_ui_1-group_by",
                selected = "RRE-Finder"
            )
            shiny::updateSelectInput(session, "ref_col_biocircos",
                selected = "RRE-Finder"
            )
        }
        if (!is.null(vals$rre_data$Probability)) {
            vals$rre_more <- TRUE
        } else {
            vals$rre_more <- FALSE
        }
    }


    #----------------------------------------------------------------
    ##            Loading and processing of example data             -
    ## ----------------------------------------------------------------
    shiny::observeEvent(input$ripp_sco, {
      process_rippminer(BGCViz:::ripp_data, example_data = TRUE)
    })
    
    shiny::observeEvent(input$anti_sco, {
        process_antismash(BGCViz:::anti_data, example_data = TRUE)
    })

    shiny::observeEvent(input$gecco_sco, {
        process_gecco(BGCViz:::gecco_data, example_data = TRUE)
    })

    shiny::observeEvent(input$prism_sco, {
        process_prism(BGCViz:::prism_data, example_data = TRUE)
    })

    shiny::observeEvent(input$sempi_sco, {
        process_sempi(BGCViz:::sempi_data, example_data = TRUE)
    })

    shiny::observeEvent(input$arts_sco, {
        process_arts_archive(BGCViz:::arts_data, example_data = TRUE)
    })

    shiny::observeEvent(input$deep_sco, {
        process_deep(BGCViz:::deep_data, example_data = TRUE)
    })

    shiny::observeEvent(input$rre_sco, {
        process_rre(BGCViz:::rre_data, example_data = TRUE)
    })

    ## ----------------------------------------------------------------
    ##                Loading and processing user data               -
    ## ----------------------------------------------------------------
    
    shiny::observeEvent(input$ripp_data, {
      
      # Read data
      ripp_data <- utils::read.delim(input$ripp_data$datapath)
      process_rippminer(ripp_data)
    })
    
    shiny::observeEvent(input$anti_data, {
        disable_event_logic()
        # Read data
        if (input$anti_data$type == "text/csv") {
            anti_data <- utils::read.csv(input$anti_data$datapath)
        } else {
            data <- rjson::fromJSON(file = input$anti_data$datapath)
            types <- sapply(data$records, function(y) {
                lapply(y$features, function(x) {
                    if (unlist(x$type == "region")) {
                        tolower(x$qualifiers$product)
                    }
                })
            })

            types <- Filter(Negate(is.null), types)

            types <- sapply(types, function(x) {
                if (length(unlist(x)) > 1) {
                    tmp <- stringr::str_trim(paste0(unlist(x), collapse = "", sep = " "))
                    gsub(" ", "__", tmp)
                } else {
                    x
                }
            })

            location <- sapply(data$records, function(y) {
                unlist(sapply(y$features, function(x) {
                    if (unlist(x$type == "region")) {
                        unlist(x$location)
                    }
                }))
            })


            location <- gsub("\\[", "", location)
            location <- gsub("\\]", "", location)
            location <- gsub("<", "", location)
            location <- gsub(">", "", location)
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

        process_antismash(anti_data)
    })

    shiny::observeEvent(input$sempi_data, {
        if (input$sempi_data$type == "text/csv") {
            sempi_data <- utils::read.csv(input$sempi_data$datapath)
            process_sempi(sempi_data, zip = FALSE)
        } else {
            process_sempi(input$sempi_data$datapath, zip = TRUE)
        }
    })

    shiny::observeEvent(input$gecco_data, {
        gecco_data <- utils::read.delim(input$gecco_data$datapath)
        process_gecco(gecco_data)
    })

    # These are for ARTS data processing
    # input$known_data and inoput$dup_data

    shiny::observeEvent(input$arts_data, {
        disable_event_logic()
        if (input$arts_data$type == "text/csv") {
            process_arts_archive(input$arts_data$datapath, zip = FALSE)
        } else {
            process_arts_archive(input$arts_data$datapath, zip = TRUE)
        }
    })


    shiny::observeEvent(input$prism_data, {

        # Read data
        if (input$prism_data$type == "text/csv") {
            prism_data <- utils::read.csv(input$prism_data$datapath)
            process_prism(prism_data, json = FALSE)
        } else {
            data <- rjson::fromJSON(file = input$prism_data$datapath)
            process_prism(data)
        }
    })

    shiny::observeEvent(input$deep_data, {
        data <- utils::read.delim(input$deep_data$datapath)
        process_deep(data)
    })

    shiny::observeEvent(input$rre_data, {

        # Read data
        rre_data <- utils::read.delim(input$rre_data$datapath)
        process_rre(rre_data)
    })

    ############################################################################
    ############################################################################
    ###                                                                      ###
    ###                INTERFACE LOGIC: WHAT TO SHOW AND WHEN                ###
    ###                                                                      ###
    ############################################################################
    ############################################################################
    # Update choices
    update_ui_with_data <- function() {
        shiny::updateSelectInput(session, "ref",
            choices = vals$choices$ref
        )
        shiny::updateSelectInput(session, "group_table_ui_1-group_by",
            choices = vals$choices$group_by
        )
        shiny::updateSelectInput(session, "ref_col_biocircos",
            choices = vals$choices$ref_col_biocircos
        )
        shiny::updateSelectInput(session, "gecco_plots_ui_1-ref_comparison_gecco",
            choices = vals$choices$ref_comparison_gecco
        )
        shiny::updateSelectInput(session, "deep_barplot_ui_1-ref_comparison",
            choices = vals$choices$ref_comparison
        )
    }
    # Observe input of chromosome length
    shiny::observeEvent(input$chr_len, {
        vals$chr_len <- input$chr_len
    })
    ## ----------------------------------------------------------------
    ##    Simple options showing/hiding logic for every data input   -
    ## ----------------------------------------------------------------
    # SHOW rre_width parameter if data is available
    # and hide_viz == FALSE
    shiny::observeEvent(vals$rre_data_input, {
        if (vals$rre_data_input == TRUE) {
            shinyjs::showElement(selector = "#rre_width")
        } else {
            shinyjs::hideElement(selector = "#rre_width")
        }
    })
    # Show anti_hybrid option if data is available
    # And checkbox is unchecked
    shiny::observeEvent(vals$anti_data_input, {
        if (vals$anti_data_input == TRUE) {
            shinyjs::showElement(selector = "#anti_hybrid")
        } else {
            shinyjs::hideElement(selector = "#anti_hybrid")
        }
    })
    
    # Show ripp_hybrid options
    shiny::observeEvent(vals$ripp_data_input, {
      if (vals$ripp_data_input == TRUE){
        shinyjs::showElement(selector = "#ripp_hybrid")
      } else {
        shinyjs::hideElement(selector = "#ripp_hybrid")
      }
    })
    # Show prism options if data is available
    # If hide anti is FALSE (checkbox), then show them
    # Only if prism_json file, then show Prism-Supp
    # And if hide_viz == FALSE, and prism_json, then
    # show width
    shiny::observeEvent(vals$prism_data_input, {
        if (vals$prism_data_input == TRUE) {
            shinyjs::showElement(selector = "#prism_hybrid")
            if (vals$prism_json == TRUE) {
                shinyjs::showElement(selector = "#prism_supp")
            }
            if (vals$prism_json == TRUE) {
                shinyjs::showElement(selector = "#prism_supp_data_input_width")
            }
        } else {
            shinyjs::hideElement(selector = "#prism_header")
            shinyjs::hideElement(selector = "#prism_hybrid")
            shinyjs::hideElement(selector = "#prism_supp")
            shinyjs::hideElement(selector = "#prism_supp_data_input_width")
        }
    })
    # Show SEMPI elements on data upload
    shiny::observeEvent(vals$sempi_data_input, {
        if (vals$sempi_data_input == TRUE) {
            shinyjs::showElement(selector = "#sempi_hybrid")
            shinyjs::showElement(selector = "#sempi_width")
        } else {
            shinyjs::hideElement(selector = "#sempi_hybrid")
            shinyjs::hideElement(selector = "#sempi_width")
        }
    })
    # Ahow ARTS data options, if data is available
    shiny::observeEvent(vals$arts_data_input, {
        if (vals$arts_data_input == TRUE) {
            shinyjs::showElement(selector = "#dup_choice")
            shinyjs::showElement(selector = "#arts_width")
            shinyjs::showElement(selector = "#phylo_file")
        } else {
            shinyjs::hideElement(selector = "#dup_choice")
            shinyjs::hideElement(selector = "#arts_width")
            shinyjs::hideElement(selector = "#phylo_file")
        }
    })

    shiny::observeEvent(vals$data_upload_count, {
        if ((vals$arts_data_input == TRUE) || (vals$sempi_data_input == TRUE) || (vals$prism_supp_data_input == TRUE) || (vals$rre_data_input == TRUE)) {
            shinyjs::showElement(selector = "#improve_visualization_box")
        } else {
            shinyjs::hideElement(selector = "#improve_visualization_box")
        }
    })
    shiny::observeEvent(vals$data_upload_count, {
        if ((vals$arts_data_input == TRUE) || (vals$prism_json == TRUE)) {
            shinyjs::showElement(selector = "#prism_supplement_arts_box")
        } else {
            shinyjs::hideElement(selector = "#prism_supplement_arts_box")
        }
    })
    ## ---------------------------------------------------------------
    ##              Data processing options show/hide               -
    ## ---------------------------------------------------------------
    # Count data uploads, to show tabs and corresponding
    # options

    output$deep_sidemenu_out <- shinydashboard::renderMenu({
        if (vals$data_upload_count >= 2) {
            if ((vals$deep_data_input == TRUE) & ((vals$anti_data_input == TRUE) | (vals$prism_data_input == TRUE) | (vals$sempi_data_input == TRUE))) {
                shinydashboard::menuItem("Compare data with DeepBGC",
                    tabName = "deep_sidemenu", icon = shiny::icon("dyalog"),
                    shinydashboard::menuItem("Compare with DeepBGC plots", tabName = "deep_sidemenu", icon = shiny::icon("chart-pie")),
                    shinydashboard::menuItem("Filtering options",
                        tabName = "deep_filter", icon = shiny::icon("filter"),
                        shiny::uiOutput("deep_filter_UI_sidemenu")
                    )
                )
            }
        }
    })
    output$gecco_sidemenu_out <- shinydashboard::renderMenu({
        if (vals$data_upload_count >= 2) {
            if ((vals$gecco_data_input == TRUE) & ((vals$anti_data_input == TRUE) | (vals$prism_data_input == TRUE) | (vals$sempi_data_input == TRUE))) {
                shinydashboard::menuItem("Compare data with GECCO",
                    tabName = "gecco", icon = icon("fas fa-dragon"),
                    shinydashboard::menuItem("Compare with GECCO plots", tabName = "gecco_sidemenu", icon = shiny::icon("chart-pie")),
                    shinydashboard::menuItem("Filtering options",
                        tabName = "gecco_filter", icon = shiny::icon("filter"),
                        shiny::uiOutput("gecco_filter_UI_sidemenu")
                    )
                )
            }
        }
    })
    output$anno_sidemenu_out <- shinydashboard::renderMenu({
        if (vals$data_upload_count >= 1) {
            shinydashboard::menuItem("Annotation visualization and comparison", tabName = "anno_sidemenu", icon = icon("fas fa-project-diagram"))
        }
    })
    output$biocircos_sidemenu_out <- shinydashboard::renderMenu({
        if (vals$data_upload_count >= 2) {
            shinydashboard::menuItem("Biocircos plot", tabName = "biocircos_sidemenu", icon = icon("fas fa-circle-notch"))
        }
    })
    output$summarize_sidemenu_out <- shinydashboard::renderMenu({
        if (vals$data_upload_count >= 2) {
            shinydashboard::menuItem("Summarize interception", tabName = "summarize_sidemenu", icon = icon("fas fa-chart-bar"))
        }
    })
    output$arts_tree_sidemenu_out <- shinydashboard::renderMenu(
      {
        if (vals$arts_data_input == TRUE){
            shinydashboard::menuItem('ARTS phylogeny', tabName = "arts_tree_sidemenu", icon = icon("tree")) 
        }
      }
    )

    output$deep_filter_box <- shiny::renderUI({
        if (vals$deep_data_input == TRUE) {
            vals$deep_global <- TRUE
            shinydashboardPlus::box(
                title = "DeepBGC filtering",
                id = "deep_filtering_box",
                collapsible = TRUE,
                closable = TRUE,
                width = NULL,
                shiny::sliderInput("score_a", "Activity score threshold for DeepBGC data", min = 0, max = 100, value = 50),
                shiny::sliderInput("score_d", "DeepBGC score threshold for DeepBGC data", min = 0, max = 100, value = 50),
                shiny::sliderInput("score_c", "Cluster_type score threshold for DeepBGC data", min = 0, max = 100, value = 50),
                # Domains, biodomains and proteins dplyr::filter. Remain >= of set threshold
                shiny::sliderInput("domains_filter", "Domain number threshold for DeepBGC data", min = 0, max = 100, value = 5),
                shiny::sliderInput("biodomain_filter", "Biodomain number threshold for DeepBGC data", min = 0, max = 100, value = 1),
                shiny::sliderInput("gene_filter", "Protein number threshold for DeepBGC data", min = 0, max = 100, value = 1),
                shiny::sliderInput("cluster_type", "Choose threshold to assign cluster type for DeepBGC data ", min = 0, max = 100, value = 50)
            )
        }
    })
    output$gecco_filter_box <- shiny::renderUI({
        if (vals$gecco_data_input == TRUE) {
            vals$gecco_global <- TRUE
            shinydashboardPlus::box(
                title = "GECCO filtering",
                id = "gecco_filtering_box",
                collapsible = TRUE,
                closable = TRUE,
                width = NULL,
                shiny::sliderInput("score_average_gecco", "Average p-value threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50),
                shiny::sliderInput("score_cluster_gecco", "Cluster type threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50),
                shiny::sliderInput("domains_filter_gecco", "Domain number threshold for Gecco data", min = 0, max = 100, value = 1),
                shiny::sliderInput("prot_filter_gecco", "Protein number threshold for Gecco data", min = 0, max = 100, value = 1)
            )
        }
    })

    output$deep_filter_UI_sidemenu <- shiny::renderUI({
        vals$deep_sidebar <- TRUE
        shiny::tagList(
            shiny::sliderInput("score_a_sidemenu", "Activity score threshold for DeepBGC data", min = 0, max = 100, value = 50),
            shiny::sliderInput("score_d_sidemenu", "DeepBGC score threshold for DeepBGC data", min = 0, max = 100, value = 50),
            shiny::sliderInput("score_c_sidemenu", "Cluster_type score threshold for DeepBGC data", min = 0, max = 100, value = 50),
            # Domains, biodomains and proteins dplyr::filter. Remain >= of set threshold
            shiny::sliderInput("domains_filter_sidemenu", "Domain number threshold for DeepBGC data", min = 0, max = 100, value = 5),
            shiny::sliderInput("biodomain_filter_sidemenu", "Biodomain number threshold for DeepBGC data", min = 0, max = 100, value = 1),
            shiny::sliderInput("gene_filter_sidemenu", "Protein number threshold for DeepBGC data", min = 0, max = 100, value = 1),
            shiny::sliderInput("cluster_type_sidemenu", "Choose threshold to assign cluster type for DeepBGC data ", min = 0, max = 100, value = 50)
        )
    })
    output$gecco_filter_UI_sidemenu <- shiny::renderUI({
        vals$gecco_sidebar <- TRUE
        shiny::tagList(
            shiny::sliderInput("score_average_gecco_sidemenu", "Average p-value threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50),
            shiny::sliderInput("score_cluster_gecco_sidemenu", "Cluster type threshold for Gecco data (%, mapped from 0 to 1)", min = 0, max = 100, value = 50),
            shiny::sliderInput("domains_filter_gecco_sidemenu", "Domain number threshold for Gecco data", min = 0, max = 100, value = 1),
            shiny::sliderInput("prot_filter_gecco_sidemenu", "Protein number threshold for Gecco data", min = 0, max = 100, value = 1)
        )
    })

    update_filter_values <- function(listening_value, comparing_values, updating_value, rendering_check) {
        if ((as.numeric(listening_value) != comparing_values) && (rendering_check == FALSE)) {
            shiny::updateSliderInput(session, updating_value, NULL, listening_value)
            return(list(as.numeric(listening_value), FALSE))
        } else {
            if (grepl("sidemenu", updating_value) == TRUE) {
                shiny::updateSliderInput(session, stringr::str_split(updating_value, "_sidemenu")[[1]][1], NULL, comparing_values)
            } else {
                shiny::updateSliderInput(session, paste0(updating_value, "_sidemenu")[[1]][1], NULL, comparing_values)
            }
            return(list(comparing_values, FALSE))
        }
    }


    shiny::observeEvent(input$score_a, {
        res <- update_filter_values(input$score_a, vals$score_a, "score_a_sidemenu", vals$deep_sidebar)
        vals$score_a <- res[[1]]
        vals$deep_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$score_d, {
        res <- update_filter_values(input$score_d, vals$score_d, "score_d_sidemenu", vals$deep_sidebar)
        vals$score_d <- res[[1]]
        vals$deep_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$score_c, {
        res <- update_filter_values(input$score_c, vals$score_c, "score_c_sidemenu", vals$deep_sidebar)
        vals$score_c <- res[[1]]
        vals$deep_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$domains_filter, {
        res <- update_filter_values(input$domains_filter, vals$domains_filter, "domains_filter_sidemenu", vals$deep_sidebar)
        vals$domains_filter <- res[[1]]
        vals$deep_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$biodomain_filter, {
        res <- update_filter_values(input$biodomain_filter, vals$biodomain_filter, "biodomain_filter_sidemenu", vals$deep_sidebar)
        vals$biodomain_filter <- res[[1]]
        vals$deep_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$gene_filter, {
        res <- update_filter_values(input$gene_filter, vals$gene_filter, "gene_filter_sidemenu", vals$deep_sidebar)
        vals$gene_filter <- res[[1]]
        vals$deep_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$cluster_type, {
        res <- update_filter_values(input$cluster_type, vals$cluster_type, "cluster_type_sidemenu", vals$deep_sidebar)
        vals$cluster_type <- res[[1]]
        vals$deep_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$score_a_sidemenu, {
        res <- update_filter_values(input$score_a_sidemenu, vals$score_a, "score_a", vals$deep_global)
        vals$score_a <- res[[1]]
        vals$deep_global <- res[[2]]
    })
    shiny::observeEvent(input$score_d_sidemenu, {
        res <- update_filter_values(input$score_d_sidemenu, vals$score_d, "score_d", vals$deep_global)
        vals$score_d <- res[[1]]
        vals$deep_global <- res[[2]]
    })
    shiny::observeEvent(input$score_c_sidemenu, {
        res <- update_filter_values(input$score_c_sidemenu, vals$score_c, "score_c", vals$deep_global)
        vals$score_c <- res[[1]]
        vals$deep_global <- res[[2]]
    })
    shiny::observeEvent(input$domains_filter_sidemenu, {
        res <- update_filter_values(input$domains_filter_sidemenu, vals$domains_filter, "domains_filter", vals$deep_global)
        vals$domains_filter <- res[[1]]
        vals$deep_global <- res[[2]]
    })
    shiny::observeEvent(input$biodomain_filter_sidemenu, {
        res <- update_filter_values(input$biodomain_filter_sidemenu, vals$biodomain_filter, "biodomain_filter", vals$deep_global)
        vals$biodomain_filter <- res[[1]]
        vals$deep_global <- res[[2]]
    })
    shiny::observeEvent(input$gene_filter_sidemenu, {
        res <- update_filter_values(input$gene_filter_sidemenu, vals$gene_filter, "gene_filter", vals$deep_global)
        vals$gene_filter <- res[[1]]
        vals$deep_global <- res[[2]]
    })
    shiny::observeEvent(input$cluster_type_sidemenu, {
        res <- update_filter_values(input$cluster_type_sidemenu, vals$cluster_type, "cluster_type", vals$deep_global)
        vals$cluster_type <- res[[1]]
        vals$deep_global <- res[[2]]
    })



    shiny::observeEvent(input$score_average_gecco, {
        res <- update_filter_values(input$score_average_gecco, vals$score_average_gecco, "score_average_gecco_sidemenu", vals$gecco_sidebar)
        vals$score_average_gecco <- res[[1]]
        vals$gecco_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$score_cluster_gecco, {
        res <- update_filter_values(input$score_cluster_gecco, vals$score_cluster_gecco, "score_cluster_gecco_sidemenu", vals$gecco_sidebar)
        vals$score_cluster_gecco <- res[[1]]
        vals$gecco_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$domains_filter_gecco, {
        res <- update_filter_values(input$domains_filter_gecco, vals$domains_filter_gecco, "domains_filter_gecco_sidemenu", vals$gecco_sidebar)
        vals$domains_filter_gecco <- res[[1]]
        vals$gecco_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$prot_filter_gecco, {
        res <- update_filter_values(input$prot_filter_gecco, vals$prot_filter_gecco, "prot_filter_gecco_sidemenu", vals$gecco_sidebar)
        vals$prot_filter_gecco <- res[[1]]
        vals$gecco_sidebar <- res[[2]]
    })
    shiny::observeEvent(input$score_average_gecco_sidemenu, {
        res <- update_filter_values(input$score_average_gecco_sidemenu, vals$score_average_gecco, "score_average_gecco", vals$gecco_global)
        vals$score_average_gecco <- res[[1]]
        vals$gecco_global <- res[[2]]
    })
    shiny::observeEvent(input$score_cluster_gecco_sidemenu, {
        res <- update_filter_values(input$score_cluster_gecco_sidemenu, vals$score_cluster_gecco, "score_cluster_gecco", vals$gecco_global)
        vals$score_cluster_gecco <- res[[1]]
        vals$gecco_global <- res[[2]]
    })
    shiny::observeEvent(input$domains_filter_gecco_sidemenu, {
        res <- update_filter_values(input$domains_filter_gecco_sidemenu, vals$domains_filter_gecco, "domains_filter_gecco", vals$gecco_global)
        vals$domains_filter_gecco <- res[[1]]
        vals$gecco_global <- res[[2]]
    })
    shiny::observeEvent(input$prot_filter_gecco_sidemenu, {
        res <- update_filter_values(input$prot_filter_gecco_sidemenu, vals$prot_filter_gecco, "prot_filter_gecco", vals$gecco_global)
        vals$prot_filter_gecco <- res[[1]]
        vals$gecco_global <- res[[2]]
    })

    shiny::observeEvent(input$restore_box, {
        box_ids <- c(
            "deep_comparison_box", "deep_rate_box", "deep_comparison_controls_box", "gecco_comparison_box",
            "gecco_rate_box", "gecco_comparison_controls_box", "annotation_reference_box", "annotation_reference_comparison_box",
            "annotation_reference_comparison_controls_box", "biocircos_plot_box", "biocircos_controls_box",
            "ranking_barplot_box", "group_table_box", "upload_anti_box","upload_ripp_box", "upload_prism_box",
            "upload_sempi_box", "upload_deep_box", "upload_gecco_box", "upload_rre_box", "upload_arts_box",
            "use_example_data_box", "rename_box", "prism_supplement_arts_box", "improve_visualization_box",
            "download_data_box", "gecco_filtering_box", "deep_filtering_box", "arts_tree_box"
        )
        for (id in box_ids) {
            shinydashboardPlus::updateBox(id, action = "restore")
        }
    })


    # Logic show/hide selectinput in Link coloring in
    # Biocircos
    shiny::observeEvent(input$label_color_class, {
        if (input$label_color_class == "R") {
            shinyjs::showElement(selector = "#ref_col_biocircos")
        } else {
            shinyjs::hideElement(selector = "#ref_col_biocircos")
        }
    })
    # Make hybrids from the data, if checkbox is checked
    # TODO Put the function to the root.
    # Tou have duplicated code
    shiny::observeEvent(input$ripp_hybrid, ignoreInit = TRUE, {
      if (input$ripp_hybrid == TRUE) {
        vals$ripp_data$Type2 <- hybrid_col(vals$ripp_data)
      } else {
        vals$ripp_data$Type2 <- vals$ripp_type
      }
    })
    shiny::observeEvent(input$anti_hybrid, ignoreInit = TRUE, {
        if (input$anti_hybrid == TRUE) {
            vals$anti_data$Type2 <- hybrid_col(vals$anti_data)
        } else {
            vals$anti_data$Type2 <- vals$anti_type
        }
    })
    shiny::observeEvent(input$prism_hybrid, ignoreInit = TRUE, {
        if (input$prism_hybrid == TRUE) {
            vals$prism_data$Type2 <- hybrid_col(vals$prism_data)
        } else {
            vals$prism_data$Type2 <- vals$prism_type
        }
    })
    shiny::observeEvent(input$sempi_hybrid, ignoreInit = TRUE, {
        if (input$sempi_hybrid == TRUE) {
            vals$sempi_data$Type2 <- hybrid_col(vals$sempi_data)
        } else {
            vals$sempi_data$Type2 <- vals$sempi_type
        }
    })
    # Rename the data, if button is clicked
    shiny::observeEvent(input$rename, {
        rename_data <- vals$rename_data
        if (vals$anti_data_input == TRUE) {
            anti_data <- vals$anti_data
            res <- rename_vector(anti_data, rename_data, vals$renaming_notification)
            vals$anti_type <- res[[1]]
            vals$renaming_notification <- res[[2]]
            anti_data["Type2"] <- vals$anti_type
            vals$anti_data <- anti_data
        }

        if (vals$sempi_data_input == TRUE) {
            sempi_data <- vals$sempi_data
            res <- rename_vector(sempi_data, rename_data, vals$renaming_notification)
            vals$sempi_type <- res[[1]]
            vals$renaming_notification <- res[[2]]
            sempi_data["Type2"] <- vals$sempi_type
            vals$sempi_data <- sempi_data
        }

        if (vals$prism_data_input == TRUE) {
            prism_data <- vals$prism_data
            res <- rename_vector(prism_data, rename_data, vals$renaming_notification)
            vals$prism_type <- res[[1]]
            vals$renaming_notification <- res[[2]]
            prism_data["Type2"] <- vals$prism_type
            vals$prism_data <- prism_data
        }
        if (vals$ripp_data_input == TRUE) {
            ripp_data <- vals$ripp_data
            res <- rename_vector(ripp_data, rename_data, vals$renaming_notification)
            vals$ripp_type <- res[[1]]
            vals$renaming_notification <-res[[2]]
            ripp_data["Type2"] <- vals$ripp_data
            vals$ripp_data <- ripp_data
        }
        shinyjs::showElement(selector = "#reset_name")
        shinyjs::hideElement(selector = "#rename")
        vals$renamed <- TRUE
        shiny::showNotification(paste("Please note: SEMPI, PRISM and Antismash input data will be renamed on upload"), type = "warning", duration = 10)
    })
    # When the new data is uploaded and renamed
    # is TRUE, then rename data on upload
    shiny::observeEvent(check_to_rename(), {
        shiny::req(vals$renamed == TRUE)

        rename_data <- vals$rename_data
        if (vals$anti_data_input == TRUE) {
            anti_data <- vals$anti_data
            res <- rename_vector(anti_data, rename_data, vals$renaming_notification)
            vals$anti_type <- res[[1]]
            vals$renaming_notification <- res[[2]]
            anti_data["Type2"] <- vals$anti_type
            vals$anti_data <- anti_data
        }

        if (vals$sempi_data_input == TRUE) {
            sempi_data <- vals$sempi_data
            res <- rename_vector(sempi_data, rename_data, vals$renaming_notification)
            vals$sempi_type <- res[[1]]
            vals$renaming_notification <- res[[2]]
            sempi_data["Type2"] <- vals$sempi_type
            vals$sempi_data <- sempi_data
        }

        if (vals$prism_data_input == TRUE) {
            prism_data <- vals$prism_data
            res <- rename_vector(prism_data, rename_data, vals$renaming_notification)
            vals$prism_type <- res[[1]]
            vals$renaming_notification <- res[[2]]
            prism_data["Type2"] <- vals$prism_type
            vals$prism_data <- prism_data
        }
        if (vals$ripp_data_input == TRUE) {
            ripp_data <- vals$ripp_data
            res <- rename_vector(ripp_data, rename_data, vals$renaming_notification)
            vals$ripp_type <- res[[1]]
            vals$renaming_notification <-res[[2]]
            ripp_data["Type2"] <- vals$ripp_data
            vals$ripp_data <- ripp_data
        }
    })
    # Reset the renaming. Uncheck the hybrid checkboxes
    shiny::observeEvent(input$reset_name, {
        vals$anti_data["Type2"] <- vals$anti_data["Type"]
        vals$sempi_data["Type2"] <- vals$sempi_data["Type"]
        vals$ prism_data["Type2"] <- vals$ prism_data["Type"]
        vals$ripp_data["Type2"] <- vals$ripp_data["Type"]
        if (input$ripp_hybrid == TRUE) {
            shiny::showNotification(paste("RippMiner cluster types are NOT visualized as hybrid anymore. You should check the option one more time"), type = "warning", duration = 10 )
            shiny::showNotification(inputId ="ripp_hybrid", value = FALSE)
        }
        if (input$anti_hybrid == TRUE) {
            shiny::showNotification(paste("Antismash cluster types are NOT visualized as hybrid anymore. You should check the option one more time"), type = "warning", duration = 10)
            shiny::updateCheckboxInput(inputId = "anti_hybrid", value = FALSE)
        }
        if (input$prism_hybrid == TRUE) {
            shiny::showNotification(paste("PRISM cluster types are NOT visualized as hybrid anymore. You should check the option one more time"), type = "warning", duration = 10)
            shiny::updateCheckboxInput(inputId = "prism_hybrid", value = FALSE)
        }
        if (input$sempi_hybrid == TRUE) {
            shiny::showNotification(paste("SEMPI cluster types are NOT visualized as hybrid anymore. You should check the option one more time"), type = "warning", duration = 10)
            shiny::updateCheckboxInput(inputId = "sempi_hybrid", value = FALSE)
        }
        shinyjs::showElement(selector = "#rename")
        shinyjs::hideElement(selector = "#reset_name")
        vals$renamed <- FALSE
    })
    # Read the uploaded renaming scheme csv
    shiny::observeEvent(input$rename_data, {
        rename_data <- utils::read.csv(input$rename_data$datapath)
        vals$rename_data <- rename_data
        coloring_datatable <- data.frame(tidyr::drop_na(data.frame(cbind(as.character(rename_data$Group_color), as.character(rename_data$Color), rename_data$Hierarchy))))
        coloring_datatable <- coloring_datatable[!apply(coloring_datatable == "", 1, all), ]
        colnames(coloring_datatable) <- c("Name", "Color", "Hierarchy")
        vals$coloring_datatable <- DT::datatable(coloring_datatable, rownames = FALSE, editable = "column")
    })


    # What to do, if hide DeepBGC comparison options scheme is triggered


    ############################################################################
    ############################################################################
    ###                                                                      ###
    ###                             COMPUTATIONS                             ###
    ###                                                                      ###
    ############################################################################
    ############################################################################
    shiny::observeEvent(input$prism_supp, ignoreInit = TRUE, priority = 3, {
        if (input$prism_supp == TRUE) {
            vals$need_filter <- TRUE
            vals$prism_supp_data_input <- TRUE
            vals$prism_supp_plot <- TRUE
            if (!("PRISM-Supp" %in% names(vals$choices$ref))) {
                vals$choices$ref <- c(vals$choices$ref, "PRISM-Supp" = "PRISM-Supp")
                vals$choices$group_by <- c(vals$choices$group_by, "PRISM-Supp" = "PRISM-Supp")
                vals$choices$ref_col_biocircos <- c(vals$choices$ref_col_biocircos, "PRISM-Supp" = "PRISM-Supp")
                update_ui_with_data()
            }
        } else {
            vals$prism_supp_data_input <- FALSE
            vals$need_filter <- TRUE
            vals$prism_supp_plot <- FALSE
            vals$choices$ref <- vals$choices$ref[!(names(vals$choices$ref) %in% c("PRISM-Supp"))]
            vals$choices$group_by <- vals$choices$group_by[!(names(vals$choices$group_by) %in% c("PRISM-Supp"))]
            vals$choices$ref_col_biocircos <- vals$choices$ref_col_biocircos[!(names(vals$choices$ref_col_biocircos) %in% c("PRISM-Supp"))]
            update_ui_with_data()
        }
    })

    # Compute all interceptions on data upload.
    # dplyr::filter while ploting then.
    shiny::observeEvent(inputData(), ignoreInit = TRUE, priority = 5, {
        # GENERATE DATA
        if (vals$ripp_data_input == TRUE) {
          ripp_data <- vals$ripp_data
          ripp_inter <- vals$ripp_data %>%
            dplyr::select(Start, Stop)
          ripp_inter$seqnames <- "chr"
        }
      
        if (vals$anti_data_input == TRUE) {
            anti_data <- vals$anti_data
            anti_inter <- vals$anti_data %>%
                dplyr::select(Start, Stop)
            anti_inter$seqnames <- "chr"
        }
        if (vals$deep_data_input == TRUE) {
            deep_data <- vals$deep_data
            deep_inter <- vals$deep_data %>%
                dplyr::select(Start, Stop)
            deep_inter$seqnames <- "chr"
        }
        if (vals$rre_data_input == TRUE) {
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
        if (vals$prism_data_input == TRUE) {
            # Store master prism data in local variable
            prism_data <- vals$prism_data
            # Start/Stop columns from prism data as matrix
            prism_inter <- prism_data %>%
                dplyr::select(Start, Stop)
            prism_inter$seqnames <- "chr"
        }
        if (vals$sempi_data_input == TRUE) {
            # Store master prism data in local variable
            sempi_data <- vals$sempi_data
            # Start/Stop columns from prism data as matrix
            sempi_inter <- vals$sempi_data %>%
                dplyr::select(Start, Stop)
            sempi_inter$seqnames <- "chr"
        }
        if (vals$prism_json == TRUE) {
            prism_supp_data <- vals$prism_supp_data
            prism_supp_inter <- vals$prism_supp_data %>%
                dplyr::select(Start, Stop)
            prism_supp_inter$seqnames <- "chr"
        }
        if (vals$arts_data_input == TRUE) {
            arts_data <- vals$arts_data
            arts_inter <- vals$arts_data %>%
                dplyr::select(Start, Stop)
            arts_inter$seqnames <- "chr"
        }
        if (vals$gecco_data_input == TRUE) {
            gecco_data <- vals$gecco_data
            # Start/Stop columns from prism data as matrix
            gecco_inter <- vals$gecco_data %>%
                dplyr::select(Start, Stop)
            gecco_inter$seqnames <- "chr"
        }

        get_inter <- function(inter1, inter2) {
            query <- GenomicRanges::makeGRangesFromDataFrame(inter2)
            subject <- GenomicRanges::makeGRangesFromDataFrame(inter1)
            interseption <- GenomicRanges::findOverlaps(query, subject)
            inter_from <- interseption@from
            inter_to <- interseption@to
            return(list(from = inter_from, to = inter_to))
        }

        inters <- vals$inters
        index <- 1
        for (i in data_uploads_inter) {
            index_2 <- 1
            j <- soft_names[index]
            for (p in data_uploads_inter) {
                x <- soft_names[index_2]
                if ((vals[[i]] == TRUE) & (vals$computed[[j]] == FALSE) & (j != x)) {
                    if ((vals[[p]] == TRUE) & (j != soft_names[index_2])) {
                        res <- get_inter(eval(as.name(paste(j, "_inter", sep = ""))), eval(as.name(paste(x, "_inter", sep = ""))))
                        new_res <- list()
                        new_res$from <- eval(as.name(paste(x, "_data", sep = "")))[res$from, ]$Cluster
                        new_res$to <- eval(as.name(paste(j, "_data", sep = "")))[res$to, ]$Cluster
                        inters[[j]][[x]] <- new_res
                        inters[[x]][[j]] <- list(from = new_res$to, to = new_res$from)
                    }
                }
                index_2 <- index_2 + 1
            }
            if (vals[[i]] == TRUE) {
                vals$computed[[j]] <- TRUE
            }
            index <- index + 1
        }

        vals$inters <- inters
        if ((vals$deep_data_input == FALSE) & (vals$gecco_data_input == FALSE) & (vals$arts_data_input == FALSE)) {
            vals$inters_filtered <- inters
            enable_event_logic()
        } else {
            vals$need_filter <- TRUE
            vals$filter_data <- TRUE
        }
    })
    # dplyr::filter ARTS, DeepBGC, GECCO interception data
    # and general dataframes to plot, if data filtering
    # options are triggered
    shiny::observeEvent(
        {
            dynamicInput()
            to_debounce()
        },
        ignoreInit = TRUE,
        priority = 4,
        {
            shiny::req(vals$data_upload_count >= 1)
            inters <- vals$inters
            if (vals$deep_data_input == TRUE) {
                if (vals$need_filter == FALSE) {
                    biocircos_deep <- filter_deepbgc(vals$deep_data, vals$cluster_type, vals$score_a, vals$score_c, vals$score_d, vals$domains_filter, vals$biodomain_filter, vals$gene_filter)
                    vals$deep_data_filtered <- biocircos_deep
                } else {
                    biocircos_deep <- vals$deep_data_filtered
                }
                if (vals$data_upload_count != 1) {
                    new_deep <- lapply(inters$deep, function(x) {
                        new_to <- x$to[x$to %in% biocircos_deep$Cluster]
                        new_from <- x$from[x$to %in% biocircos_deep$Cluster]
                        list(from = new_from, to = new_to)
                    })
                    new_inters <- inters
                    update_list <- names(inters$deep)
                    for (b in seq(1:length(update_list))) {
                        new_inters[[update_list[b]]]$deep$to <- new_deep[[update_list[b]]]$from
                        new_inters[[update_list[b]]]$deep$from <- new_deep[[update_list[b]]]$to
                    }
                    new_inters$deep <- new_deep
                    vals$inters_filtered <- new_inters
                    inters <- new_inters
                }
            }
            if (vals$gecco_data_input == TRUE) {
                if (vals$need_filter == FALSE) {
                    gecco_data <- filter_gecco(vals$gecco_data, vals$score_cluster_gecco, vals$score_average_gecco, vals$domains_filter_gecco, vals$prot_filter_gecco)
                    vals$gecco_data_filtered <- gecco_data
                } else {
                    gecco_data <- vals$gecco_data_filtered
                }
                if (vals$data_upload_count != 1) {
                    new_gecco <- lapply(inters$gecco, function(x) {
                        new_to <- x$to[x$to %in% gecco_data$Cluster]
                        new_from <- x$from[x$to %in% gecco_data$Cluster]
                        list(from = new_from, to = new_to)
                    })
                    new_inters <- inters
                    update_list <- names(inters$gecco)
                    for (b in seq(1:length(update_list))) {
                        new_inters[[update_list[b]]]$gecco$to <- new_gecco[[update_list[b]]]$from
                        new_inters[[update_list[b]]]$gecco$from <- new_gecco[[update_list[b]]]$to
                    }
                    new_inters$gecco <- new_gecco
                    vals$inters_filtered <- new_inters
                    inters <- new_inters
                }
            }
            if (vals$arts_data_input == TRUE) {
                if (input$dup_choice != "All") {
                    vals$arts_data_filtered <- data.frame(vals$arts_data) %>%
                        dplyr::filter(Core == stringr::str_split(stringr::str_split(input$dup_choice, " ,")[[1]][[2]], "Core:")[[1]][[2]] | Core == "Not_core")
                    if (vals$data_upload_count != 1) {
                        new_arts <- lapply(inters$arts, function(x) {
                            new_to <- x$to[x$to %in% vals$arts_data_filtered$Cluster]
                            new_from <- x$from[x$to %in% vals$arts_data_filtered$Cluster]
                            list(from = new_from, to = new_to)
                        })
                        new_inters <- inters
                        update_list <- names(inters$arts)
                        for (b in seq(1:length(update_list))) {
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
            if (input$prism_supp == FALSE) {
                inters$prism_supp <- NULL
                for (name in names(inters)) {
                    inters[[name]][which(names(inters[[name]]) %in% c("prism_supp"))] <- NULL
                }
            }
            if ((vals$gecco_data_input == FALSE) & (vals$deep_data_input == FALSE) & (vals$arts_data_input == FALSE)) {
                vals$inters_filtered <- inters
            }
            vals$need_filter <- FALSE
            vals$filter_data <- FALSE
            vals$can_plot_deep_ref <- TRUE
            enable_event_logic()
        }
    )
    # Compute the Biociros plot. Store information to plot later
    shiny::observeEvent(biocircos_listen(), ignoreInit = TRUE, priority = 3, {
        shiny::req(vals$data_upload_count >= 2)
        shiny::req(vals$need_filter == FALSE)
        shiny::req(vals$can_plot_biocircos == TRUE)
        ## source("src/biocircos_functions.R")
        # BioCircos!
        Biocircos_chromosomes <- list()
        arcs_chromosomes <- c()
        arcs_begin <- c()
        arcs_end <- c()
        arc_labels <- c()
        arc_col <- c()

        if (is.null(vals$inters_filtered)) {
            inters <- vals$inters
        } else {
            inters <- vals$inters_filtered
        }


        rename_data <- vals$rename_data
        coloring_datatable <- vals$coloring_datatable

        index <- 1
        # browser()
        for (upload in data_uploads) {
            if (vals[[upload]] == TRUE) {
                # Store data in local variable
                corrected_data <- correct_width(vals[[data_to_use[index]]], soft_namings[index], input$sempi_width, input$prism_supp_data_input_width, input$arts_width, input$rre_width)
                init_data <- initialize_biocircos(corrected_data, soft_namings[index], Biocircos_chromosomes, arcs_chromosomes, arcs_begin, arcs_end, arc_labels, arc_col, rename_data, vals$chr_len, input$biocircos_color, coloring_datatable)
                # Make chromosome list for Biocircos plot. Use chr_len as an input
                Biocircos_chromosomes <- init_data[[1]]
                # Add arcs. Quantity of arcs is length of dataframes
                arcs_chromosomes <- init_data[[2]]
                # Add arcs begin positions. (Start column)
                arcs_begin <- init_data[[3]]
                # Stop position of arcs.
                arcs_end <- init_data[[4]]
                # Add Arcs labels. Can add only one label...
                arc_labels <- init_data[[5]]

                arc_col <- init_data[[6]]
            }
            index <- index + 1
        }
        # Add to tracklist. Then it can be populated with links
        tracklist <- BioCircos::BioCircosArcTrack("myArcTrack", arcs_chromosomes, arcs_begin, arcs_end,
            minRadius = 0.90, maxRadius = 0.97, labels = arc_labels, colors = arc_col
        )
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

        # CALCULATIONS
        # -----------------------------------------


        data_uploads_2 <- data_uploads
        soft_2 <- soft_namings
        soft_names_2 <- soft_names
        data_to_use_2 <- data_to_use
        index <- 1
        for (upload in data_uploads) {
            data_uploads_2 <- data_uploads_2[-1]
            soft_2 <- soft_2[-1]
            soft_names_2 <- soft_names_2[-1]
            data_to_use_2 <- data_to_use_2[-1]
            index2 <- 1
            if (vals[[upload]] == TRUE) {
                for (upload2 in data_uploads_2) {
                    if ((vals[[upload2]] == TRUE) & (length(data_uploads_2) > 0) & (soft_namings[index] != soft_2[index2])) {
                        output <- add_biocircos_data(inters[[soft_names[index]]][[soft_names_2[index2]]]$from, inters[[soft_names[index]]][[soft_names_2[index2]]]$to, vals[[data_to_use_2[index2]]], vals[[data_to_use[index]]], soft_2[index2], soft_namings[index], rename_data, input$label_color_class, input$ref_col_biocircos, coloring_datatable)
                        chromosomes_start <- c(chromosomes_start, output[[3]])
                        # Add link end. Just populate second output from the vectors, used above.
                        chromosomes_end <- c(chromosomes_end, output[[4]])
                        # Add links start positions as a start from dataframe. This vector is for chromosome start
                        link_pos_start <- as.numeric(c(link_pos_start, output[[5]]))
                        # Add links start positions as a start from dataframe. For chromosome start variable
                        link_pos_start_1 <- as.numeric(c(link_pos_start_1, output[[6]]))
                        # Add links start position for a chromosome stop variable
                        link_pos_end <- as.numeric(c(link_pos_end, output[[7]]))
                        # Add links start position for a chromosome stop position
                        link_pos_end_2 <- as.numeric(c(link_pos_end_2, output[[8]]))
                        label_1 <- c(label_1, output[[9]])
                        label_2 <- c(label_2, output[[10]])
                        label_color <- c(label_color, output[[11]])
                    }
                    index2 <- index2 + 1
                }
                utils::write.csv(vals[[data_to_use[index]]], paste0(soft_names[index], "_biocircos.csv"), row.names = FALSE)
            }
            index <- index + 1
        }




        # Combine labels with mapply to one list
        link_labels <- mapply(function(x, y) paste(x, y, sep = " | "), label_1, label_2)

        # Add links and labels to the track list for subsequent visualization
        if ((input$label_color == TRUE) & (length(chromosomes_start) > 0)) {
            group_colors <- plyr::count(unlist(label_color))
            for (i in seq(1:dim(group_colors)[1])) {
                subset <- unname(which(label_color %in% group_colors$x[i]))
                tracklist <- tracklist + BioCircos::BioCircosLinkTrack(as.character(i), chromosomes_start[subset], link_pos_start[subset],
                    link_pos_start_1[subset], chromosomes_end[subset], link_pos_end[subset],
                    link_pos_end_2[subset],
                    maxRadius = 0.85, labels = link_labels[subset],
                    displayLabel = FALSE, color = group_colors$x[i]
                )
            }
        } else if ((input$label_color == FALSE) & (length(chromosomes_start) > 0)) {
            tracklist <- tracklist + BioCircos::BioCircosLinkTrack("myLinkTrack_master", chromosomes_start, link_pos_start,
                link_pos_start_1, chromosomes_end, link_pos_end,
                link_pos_end_2,
                maxRadius = 0.85, labels = link_labels,
                displayLabel = FALSE, color = coloring_datatable$x$data$Color[coloring_datatable$x$data$Name == "base"]
            )
        } else {
            shiny::showNotification(paste("No interceptions are being made in the Biocircos plot. Please provide data with clusters that do have intercepting borders"), type = "warning", duration = NULL)
        }

        vals$tracklist <- tracklist
        vals$Biocircos_chromosomes <- Biocircos_chromosomes
    })

    shiny::observeEvent(deep_reference(), ignoreInit = TRUE, {
        shiny::req(vals$data_upload_count >= 1)
        shiny::req(vals$need_filter == FALSE)
        shiny::req(vals$can_plot_deep_ref == TRUE)
        shiny::req(input$ref != "")
        shiny::req(vals$data_upload_count >= 1)

        if (is.null(vals$inters_filtered)) {
            inters <- vals$inters
        } else {
            inters <- vals$inters_filtered
        }
        ## source("src/deep_reference_functions.R")
        # GENERATE DATA
        index <- 1
        for (upload in data_uploads) {
            if (vals[[upload]] == TRUE) {
                data <- vals[[data_to_use[index]]]
                assign(paste0(soft_names[index], "_data"), correct_width(data, soft_namings[index], input$sempi_width, input$prism_supp_data_input_width, input$arts_width, input$rre_width))
            }
            index <- index + 1
        }

        lett <- rev(LETTERS)[1:(length(data_uploads)+1)]
        

        tooltip <- c(
            "Software", "ID", "Start", "Stop", "Type", "num_domains", "deepbgc_score", "activity", "Score", "E_value",
            "P_value", "RRE_start", "RRE_stop", "Probability", "Name", "Full_name", "Hit", "Core", "Count", "Bitscore", "Model",
            "Num_domains", "Num_proteins", "Average_p", "Max_p"
        )




        # MAKE COMPUTATIONS
        sup_index <- 1
        soft_lttrs <- lett
        rename_y_axis <- vals$rename_y_axis
        rename_y_axis <- lapply(1:(length(soft_lttrs) - 1), function(x) {
            soft_lttrs[x] <- soft_namings[x]
        })
        names(rename_y_axis) <- soft_lttrs[-length(soft_lttrs)]
        for (upload in data_uploads) {
            soft_lttr <- soft_lttrs[1]
            soft_lttrs <- soft_lttrs[-1]
            if (vals[[upload]] == TRUE) {
                soft_major <- soft_names[sup_index]
                seg_ref_g <- simple_seg(eval(as.name(paste(soft_names[sup_index], "_data", sep = ""))), "Z", soft_namings[sup_index], soft_names[sup_index], soft_major, inter = FALSE, inters)
                seg_ref_g <- define_spec_seg_df(soft_names, sup_index, seg_ref_g, soft_major, eval(as.name(paste(soft_names[sup_index], "_data", sep = ""))), inter = FALSE, vals$rre_more, inters)
                seg_ref <- seg_ref_g

                if (input$ref == soft_namings[sup_index]) {
                    shiny::validate(need(nrow(eval(as.name(paste(soft_names[sup_index], "_data", sep = "")))) > 0, "Reference data is empty, and so, insufficient for plotting. Please select another one"))
                    plot <- ggplot2::ggplot(eval(as.name(paste(soft_names[sup_index], "_data", sep = ""))), ggplot2::aes(x = vals$chr_len, y = Chr)) +
                        suppressWarnings(eval(as.name(paste0("geom_", soft_names[sup_index])))(seg_ref, vals$rre_more))
                    soft_let <- abbr[sup_index]
                    lettrs <- lett[2:length(lett)]
                    labels_1 <- list()
                    index <- 1
                    for (i in data_uploads) {
                        if ((vals[[i]] == TRUE) & (soft_names[index] != soft_major)) {
                            df <- eval(as.name(paste(soft_names[index], "_data", sep = "")))
                            seg_df <- simple_seg(df, lettrs[index], soft_namings[index], soft_names[index], soft_major, inter = TRUE, inters)
                            seg_df <- define_spec_seg_df(soft_names, index, seg_df, soft_major, df, inter = TRUE, vals$rre_more, inters)
                            labels_1[[lettrs[index]]] <- (paste(abbr[index], "_vs_", soft_let, sep = ""))
                            plot <- suppressWarnings(add_more_annot(seg_df, plot, soft_names, index, vals$rre_more))
                        }
                        index <- index + 1
                    }
                    plot <- plot +
                        ggplot2::scale_y_discrete(labels = c("Z" = input$ref, unlist(labels_1))) +
                        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10)) +
                        ggplot2::ylab("") +
                        ggplot2::xlab("Chromosome length") +
                        ggplot2::theme(legend.title = ggplot2::element_blank()) +
                        ggplot2::ggtitle("Annotations' comparison to the reference")
                    to_plot <- plotly::ggplotly(plot, tooltip = tooltip)
                    to_plot <- to_plot %>%
                        plotly::layout(legend = list(
                            font = list(
                                family = "sans-serif",
                                size = 12,
                                color = "#000"
                            ),
                            bordercolor = "#FFFFFF",
                            borderwidth = 2,
                            title = list(text = "<b> Cluster Types </b>")
                        ))
                }
                seg_ref$yend <- rep(soft_lttr, length(eval(as.name(paste(soft_names[sup_index], "_data", sep = "")))$Cluster))
                seg_ref$y <- rep(soft_lttr, length(eval(as.name(paste(soft_names[sup_index], "_data", sep = "")))$Cluster))
                vals[[soft_datafr[sup_index]]] <- seg_ref
            }
            sup_index <- sup_index + 1
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

    ## ----------------------------------------------------------------
    ##                          ARTS phylogenetic tree                  -
    ## ----------------------------------------------------------------
    # Plot tree
    
    mod_arts_tree_server("arts_tree_1",vals = vals)
    
    ## ----------------------------------------------------------------
    ##                    DeepBGC Comparison tab                     -
    ## ----------------------------------------------------------------
    # Render barplot
    mod_deepbgc_plots_server("deep_barplot_ui_1", vals = vals, score_a = vals$score_a, score_d = vals$score_d, score_c = vals$score_c)

    # Render interactive plot with plotly for rates of DeepBGC data in regards with antismash data

    ## ----------------------------------------------------------------
    ##                      GECCO Comparison tab                     -
    ## ----------------------------------------------------------------
    # Render barplot
    mod_gecco_plots_server("gecco_plots_ui_1",
        vals = vals, score_average_gecco = vals$score_average_gecco,
        score_cluster_gecco = vals$score_cluster_gecco
    )
    ## ---------------------------------------------------------------
    ##              Annotation on chromosome plots' tab             -
    ## ---------------------------------------------------------------

    # Render interactive plot, which shows bgcs of antismash, intercepted with chosen app. Also all app bgs. On hover shows all available information
    # For antismash and PRISM data showed only ID, Start, Stop, Type
    mod_deep_reference_server("deep_reference_ui_1", vals = vals)

    mod_deep_reference_2_server("deep_reference_2_ui_1", vals = vals, data_uploads = data_uploads, data_to_use = data_to_use)
    ## ----------------------------------------------------------------
    ##                      Biocircos plot tab                       -
    ## ---------------------------------------------------------------
    # Render Biocircos Plot for all-vs-all comparison
    mod_biocircos_server("biocircos_ui_1", vals = vals)
    ## ---------------------------------------------------------------
    ##                        Summarize tab                         -
    ## ---------------------------------------------------------------
    # Render barplot with number plyr::count of interception for BGC IDs
    mod_barplot_rank_server("barplot_rank_ui_1", vals = vals, data_uploads = data_uploads, soft_names = soft_names, soft_namings = soft_namings, data_to_use = data_to_use, abbr = abbr)


    # Render table with data
    mod_group_table_server("group_table_ui_1", vals = vals, data_uploads = data_uploads, soft_names = soft_names, soft_namings = soft_namings, data_to_use = data_to_use, abbr = abbr)

    # Download used datasets (as for BioCircos)
    mod_download_server("download_ui_1")

    shiny::onSessionEnded(function() {
        # List files in directory
        files_in_dir <- list.files()
        # Iterate over those files and if found "_biocircos.csv" add to the flst vector
        for (file_names in files_in_dir) {
            if (grepl("_biocircos.csv", file_names, fixed = TRUE)) {
                file.remove(file_names)
            } else if (grepl("group_by.csv", file_names, fixed = TRUE)) {
                file.remove(file_names)
            }
        }

        shiny::stopApp()
    })
}
