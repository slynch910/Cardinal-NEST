# Cardinal NEST - metadata module: pixel/feature metadata explorer

library(shiny)

metadataUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        h3("Metadata Explorer"),
        p("View and filter pixel and feature metadata"),
        uiOutput(ns("active_dataset"))
      )
    ),
    
    fluidRow(
      column(6,
        bslib::card(
          bslib::card_header("Pixel Metadata"),
          bslib::card_body(
          p("View and filter pixel coordinates and run information."),
          
          selectInput(ns("pixel_view"), "View Mode",
                      choices = c("All Pixels" = "all",
                                  "Filtered" = "filtered"),
                      selected = "all"),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'filtered'", ns("pixel_view")),
            numericInput(ns("x_min"), "X Min", value = NA),
            numericInput(ns("x_max"), "X Max", value = NA),
            numericInput(ns("y_min"), "Y Min", value = NA),
            numericInput(ns("y_max"), "Y Max", value = NA),
            actionButton(ns("apply_pixel_filter"), "Apply Filter", class = "btn-primary")
          ),
          
          DT::dataTableOutput(ns("pixel_table"))
          )
        )
      ),
      
      column(6,
        bslib::card(
          bslib::card_header("Feature Metadata"),
          bslib::card_body(
          p("View and filter m/z values and feature information."),
          
          selectInput(ns("feature_view"), "View Mode",
                      choices = c("All Features" = "all",
                                  "By m/z Range" = "mz_range"),
                      selected = "all"),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'mz_range'", ns("feature_view")),
            numericInput(ns("mz_min"), "m/z Min", value = NA),
            numericInput(ns("mz_max"), "m/z Max", value = NA),
            actionButton(ns("apply_mz_filter"), "Apply Filter", class = "btn-primary")
          ),
          
          DT::dataTableOutput(ns("feature_table"))
          )
        )
      )
    ),
    
    fluidRow(
      column(12,
        bslib::card(
          bslib::card_header("Memory Management"),
          bslib::card_body(
          p("View the active dataset's storage mode and load file-backed data into memory for faster access."),
          
          uiOutput(ns("memory_status")),
          
          uiOutput(ns("load_memory_ui"))
          )
        )
      )
    )
  )
}

metadataServer <- function(id, datasets, active_id) {
  moduleServer(id, function(input, output, session) {
    
    library(Cardinal)
    library(DT)

    mse_data <- active_dataset(datasets, active_id)
    output$active_dataset <- active_dataset_ui(datasets, active_id)

    # Filter params update only on Apply, not on every keystroke.
    pixel_filter_params <- reactiveVal(NULL)
    mz_filter_params <- reactiveVal(NULL)

    observeEvent(input$apply_pixel_filter, {
      pixel_filter_params(list(
        x_min = input$x_min,
        x_max = input$x_max,
        y_min = input$y_min,
        y_max = input$y_max
      ))
    })

    observeEvent(input$apply_mz_filter, {
      mz_filter_params(list(
        mz_min = input$mz_min,
        mz_max = input$mz_max
      ))
    })

    # Reset stored params when switching back to "all".
    observeEvent(input$pixel_view, {
      if (input$pixel_view == "all") {
        pixel_filter_params(NULL)
      }
    })
    observeEvent(input$feature_view, {
      if (input$feature_view == "all") {
        mz_filter_params(NULL)
      }
    })

    pixel_data_filtered <- reactive({
      data <- mse_data()
      if (is.null(data)) return(NULL)
      
      pd <- Cardinal::pixelData(data)

      params <- pixel_filter_params()
      if (input$pixel_view == "all" || is.null(params)) {
        return(pd)
      }
      
      coord <- Cardinal::coord(data)
      filter_mask <- rep(TRUE, nrow(pd))

      if (!is.na(params$x_min)) {
        filter_mask <- filter_mask & coord$x >= params$x_min
      }
      if (!is.na(params$x_max)) {
        filter_mask <- filter_mask & coord$x <= params$x_max
      }

      if (!is.na(params$y_min)) {
        filter_mask <- filter_mask & coord$y >= params$y_min
      }
      if (!is.na(params$y_max)) {
        filter_mask <- filter_mask & coord$y <= params$y_max
      }
      
      pd[filter_mask, , drop = FALSE]
    })

    feature_data_filtered <- reactive({
      data <- mse_data()
      if (is.null(data)) return(NULL)
      
      fd <- Cardinal::featureData(data)

      params <- mz_filter_params()
      if (input$feature_view == "all" || is.null(params)) {
        return(fd)
      }
      
      mz_vals <- Cardinal::mz(data)
      filter_mask <- rep(TRUE, nrow(fd))

      if (!is.na(params$mz_min)) {
        filter_mask <- filter_mask & mz_vals >= params$mz_min
      }
      if (!is.na(params$mz_max)) {
        filter_mask <- filter_mask & mz_vals <= params$mz_max
      }
      
      fd[filter_mask, , drop = FALSE]
    })

    output$pixel_table <- DT::renderDataTable({
      pd <- pixel_data_filtered()
      if (is.null(pd)) return(NULL)
      pd_df <- as.data.frame(pd)
      
      DT::datatable(
        pd_df,
        options = list(
          pageLength = 10
        ),
        filter = 'top'
      )
    })

    output$feature_table <- DT::renderDataTable({
      fd <- feature_data_filtered()
      if (is.null(fd)) return(NULL)
      fd_df <- as.data.frame(fd)
      
      DT::datatable(
        fd_df,
        options = list(
          pageLength = 10
        ),
        filter = 'top'
      )
    })

    output$memory_status <- renderUI({
      data <- mse_data()
      if (is.null(data)) {
        return(div(class = "alert alert-info", "No data loaded"))
      }
      
      if (isTRUE(is_file_backed(data))) {
        div(class = "alert alert-info",
            strong("Status: File-backed"), br(),
            "Data is stored on disk for memory efficiency")
      } else {
        div(class = "alert alert-success", 
            strong("Status: In-memory"), br(),
            "Data is loaded into RAM for fast access")
      }
    })

    # Only offer the load button while the dataset is still file-backed.
    output$load_memory_ui <- renderUI({
      data <- mse_data()
      if (is.null(data) || !isTRUE(is_file_backed(data))) return(NULL)
      actionButton(session$ns("load_memory"), "Load into Memory", class = "btn-primary")
    })

    observeEvent(input$load_memory, {
      ds <- datasets()
      aid <- active_id()
      if (is.null(aid) || is.null(ds[[aid]])) {
        showNotification("No data loaded", type = "error")
        return()
      }
      
      tryCatch({
        data_in_memory <- withProgress(
          message = "Loading dataset into memory", value = 0.2,
          ensure_in_memory(ds[[aid]]$data)
        )
        ds[[aid]]$data <- data_in_memory
        datasets(ds)
        showNotification("Data loaded into memory", type = "message")
      }, error = function(e) {
        showNotification(paste("Failed to load into memory:", e$message), type = "error")
      })
    })
  })
}
