#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
if (!exists("setup_conda_environment")) {
  source("global.R")
}
# A reusable function to render a styled error message
# --- Server Logic ---
server <- function(input, output, session) {
  dir.create("logs", showWarnings = FALSE, recursive = TRUE)
  dir.create("www", showWarnings = FALSE, recursive = TRUE)
  track_usage(
    storage_mode = store_json(path = "logs/")
  )
  neighbour_results <- NULL
  # Per-session flag (was a module-level global shared across all users)
  hierarch_tabs_check <- reactiveVal(FALSE)
  conda_path <- Sys.which("python")
  render_error_ui("", output = output)
  
  output$csvtest <- reactive({
    !is.null(input$file1)
  })
  
  outputOptions(output, "csvtest", suspendWhenHidden = FALSE)
  #Trigger to reset heatmap
  reset_trigger <- reactiveVal(0)
  
  # Reactive value to store the TCR data
  tcr_data <- reactive({
    render_error_ui("", output = output) #Clear any error messages
    req(input$file1)
    # Determine separator based on file extension
    sep <- ifelse(tools::file_ext(input$file1$name) == "tsv", "\t", ",")
    df <- read.csv(input$file1$datapath, sep = sep)
    if (!"count" %in% colnames(df)) {
      df$count <- 1L
    }
    df
  })
  
  chains<-reactive({
    req(input$chains)
    chains <- strsplit(input$chains, ",")[[1]]
    return(chains)
  })
  
  # Run tcrdist3 analysis when the run button is clicked
  analysis_results <- eventReactive(input$run, {
    render_error_ui("", output = output) #Clear any error messages
    req(tcr_data())
    if (!("epitope" %in% colnames(tcr_data())) ||
        !("subject" %in% colnames(tcr_data())))
    {
      shinyalert(title = "Missing columns",
                 text = "Column names: 'epitope' or 'subject' missing in input data",
                 type = "error")
      return(NULL)
    }
    
    hierarch_tabs_check(FALSE)
    render_error_ui("", output = output)
    withProgress(message = 'Running TCRdist3 analysis...', value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Importing Python modules...")
        #Import necessary python libraries
        pd <- import("pandas")
        inspect <- import("inspect")
        incProgress(0.3, detail = "Initialising...")
        # Convert R data frame to pandas DataFrame
        cell_df <- r_to_py(tcr_data())
        # Get chains for analysis
        py <- import_main()
        
        chains <- chains()
        message(paste0(chains, collapse = ":"))
        
        
        incProgress(0.5, detail = "Computing distances...")
        
        # Initialize TCRrep object and compute distances
        tr <- get_epitope_dist_matrix(
          tcr_data(),
          organism = input$organism,
          chains = chains,
          dbfile = input$dbfile
        )

        #message(paste0(names(tr),collapse = "\n"))
        incProgress(0.8, detail = "Finalizing...")
        return(tr)
      }, error = function(e)
      {
        message(conda_path)
        message.error <- paste0("An error occurred: ",
                                e$message,
                                " with this python:",
                                conda_path)
        render_error_ui(message.error, output)
        message(message.error)
      })
    })
  })
  
  
  #######Update selectisize options with values based on all results###########
  observe({
    tr<-analysis_results()
    req(tr)
    #Get epitopes
    epitopes <- unique(tcr_data()$epitope)
    #Get chains from the distance compute
    all.tr_names <- names(tr)
    chains_tr <- all.tr_names[grepl("^pw_", all.tr_names)]
    chains_tr_only <- gsub("pw_", "", chains_tr[!grepl("cdr|pmhc", chains_tr)])
    genes <- names(tr$clone_df)[grepl("gene", names(tr$clone_df))]
    #Load selection input for epitopes
    output$col_select_ui <- renderUI({
      tagList(
        selectizeInput("epitope_sel", "Select Epitope:", choices = epitopes),
        selectizeInput(
          "columns_sel",
          "Select Columns:",
          choices = genes,
          multiple = T
        )
      )
    })
    
    #Load selection input for epitopes in KPCA
    output$kpca_repertoire_option <- renderUI({
      tagList(
        selectizeInput("epitope_sel_kpca", "Select Epitope:", choices = epitopes)
      )
    })
    
    #Include a select option for the chains
    output$matrixoptions <- renderUI({
      tagList(
        selectInput(
          "chain_select_heatmap",
          "Select a chain:",
          selected = chains_tr_only[1],
          choices = chains_tr_only
        ),
        selectizeInput("epitope_sel_heatmap", 
                       "Select Epitope:", 
                       choices = epitopes)
      )
    })
    
    #Include a select option for the chains in KPCA tab
    output$kpca_chain_option <- renderUI({
      tagList(
        selectInput(
          "chain_select_kpca",
          "Select a chain:",
          selected = chains_tr_only[1],
          choices = chains_tr_only
        ))
    })
    
    output$kpca_colorby_option<-renderUI({
      tagList(
        selectizeInput("colorby_sel_kpca", 
                       "Colour by:", 
                       choices = genes))
    })
    
    #Load epitopes in Neighbourhood
    output$epitope_select_ui <- renderUI({
      tagList(
        selectizeInput(
          "epitope_sel_n",
          "Select up to 2 Epitopes:",
          choices = epitopes,
          multiple = T,
          options = list(placeholder = "Select Epitope:", maxItems = 2)
        )
      )
    })
    
    output$epitope_select_ui_motif<-renderUI({
      tagList(
        selectizeInput(
          "epitope_sel_motif",
          "Select up to 2 Epitopes:",
          choices = epitopes,
          multiple = T,
          options = list(placeholder = "Select Epitope:")#, maxItems = 2)
        )
      )
    })
    #Load epitopes in trees tab
    output$epitope_select_ui_trees <- renderUI({
      tagList(
        selectizeInput(
          "epitope_sel_trees",
          "Select Epitopes:",
          selected= epitopes,
          choices = epitopes,
          multiple = T,
          options = list(placeholder = "Select Epitope:")
        )
      )
    })
    
    output$epitope_select_hierarch_ui <- renderUI({
      tagList(
        selectizeInput(
          "epitope_sel_h",
          "Select up to 2 Epitopes:",
          choices = epitopes,
          multiple = T,
          options = list(placeholder = "Select up to 2 Epitopes", 
                         maxItems =2)
        )
      )
    })
    
    #####################################################
    #Load network graph after computing distance matrix #
    #####################################################
    default_tcr_distance <- 100
    chains_net <- strsplit(input$chains, ",")[[1]]
    pw_mats <- lapply(chains_net, function(ch) tr[[paste0("pw_", ch)]])
    pw_mats <- Filter(Negate(is.null), pw_mats)
    pw_alpha_beta_max <- if (length(pw_mats) > 0) max(Reduce(`+`, pw_mats)) else 200
    
    all_color_choices <- c(
      "Epitope" = "epitope",
      "Subject" = "subject",
      "V-alpha-gene" = "v_a_gene",
      "V-beta-gene" = "v_b_gene",
      "CDR3-alpha-aa" = "cdr3_a_aa",
      "CDR3-beta-aa" = "cdr3_b_aa",
      "V-gamma-gene" = "v_g_gene",
      "V-delta-gene" = "v_d_gene",
      "CDR3-gamma-aa" = "cdr3_g_aa",
      "CDR3-delta-aa" = "cdr3_d_aa"
    )
    avail_color_choices <- all_color_choices[all_color_choices %in% names(tr$clone_df)]
    
    output$networkcard <- renderUI({
      tagList(
        hr(),
        card(
          layout_sidebar(
            sidebar = sidebar(
              h4("Network options"),
              hr(),
              helpText("Filter edges by radius, epitope, and choose how to color the nodes."),
              hr(),
              selectizeInput(
                "epitope_sel_network",
                "Select Epitopes:",
                selected = epitopes,
                choices = epitopes,
                multiple = T,
                options = list(placeholder = "Select Epitope:")
              ),
              selectizeInput(
                "chain_select_network",
                "Select chain(s):",
                choices = chains_net,
                selected = chains_net,
                multiple = TRUE,
                options = list(placeholder = "Select chain(s)")
              ),
              # Slider to filter maximum distance
              sliderInput("max_dist",
                          "TCR distance:",
                          min = 0,
                          max = pw_alpha_beta_max,
                          value = min(default_tcr_distance, pw_alpha_beta_max),
                          step = 1),
              # Dropdown to color nodes
              selectInput("color_var",
                          "Color Nodes By:",
                          choices = avail_color_choices),
              hr(),
              checkboxInput("check_net_selfloops", "Remove CDR3 self-loops"),
              hr(),
            ),
            card(
              p(strong("Network Stats:")),
              textOutput("node_count"),
              textOutput("edge_count"),
              hr(),
              withSpinner(visNetworkOutput("network_plot", height = "750px"), )
            ))
        ))
    })
    
    #Load chains in Neighbourhood
    output$chain_select_ui <- renderUI({
      tagList(
        selectizeInput(
          "chain_select_neighbours",
          "Select chain(s):",
          choices = chains_tr_only,
          selected= chains_tr_only[1:length(chains_tr_only)],
          multiple = TRUE,
          options=list(placeholder = "Select chain(s)",maxItems=2)
        )
      )
    })
    
    #Load chains in motif selection
    output$chain_select_ui_motif <- renderUI({
      tagList(
        selectizeInput(
          "chain_select_motif",
          "Select chain:",
          choices = c("",chains_tr_only),
          options=list( placeholder = "Select chain:")
        )
      )
    })
    
    output$chain_select_ui_trees <- renderUI({
      tagList(
        selectInput(
          "chain_select_tree",
          "Select chain(s):",
          selected = chains_tr_only[1],
          multiple = T,
          choices = chains_tr_only
        )
      )
    })
    
    output$chain_select_ui_knn <- renderUI({
      tagList(
        selectizeInput(
          "chain_select_knn",
          "Select a chain:",
          selected = chains_tr_only,
          choices = chains_tr_only,
          multiple = TRUE,
          options = list(placeholder = "Select a chain:")
        )
      )
    })
    
  })
  
  # Render the network visualization (top-level observer, not nested inside observe())
  observeEvent(c(input$max_dist, input$color_var, input$check_net_selfloops, input$epitope_sel_network, input$chain_select_network), {
    req(analysis_results())
    req(input$chain_select_network, input$epitope_sel_network, input$color_var, !is.null(input$max_dist))
    tryCatch({
      tr <- analysis_results()
      chains_sel <- input$chain_select_network
      tr <- subset_tr_by_epitope(tr, epitope_name = input$epitope_sel_network, chain = chains_sel)
      
      net_data <- generate_network(tr, input$max_dist, input$check_net_selfloops)
      
      if (nrow(net_data) == 0) {
        # No edges at this threshold
        output$network_plot <- renderVisNetwork({
          visNetwork(
            nodes = data.frame(id = 1, label = "No edges at this distance threshold"),
            edges = data.frame(from = numeric(0), to = numeric(0))
          )
        })
        output$node_count <- renderText("Nodes: 0")
        output$edge_count <- renderText("Edges: 0")
        return()
      }
      
      # Extract unique nodes and their metadata
      clone_df <- as.data.frame(tr$clone_df)
      all_node_ids <- unique(c(net_data$node_1, net_data$node_2))
      row_idx <- all_node_ids + 1
      
      # Create nodes data frame with color variable
      color_var <- input$color_var
      if (!color_var %in% colnames(clone_df)) {
        color_var <- "epitope"
      }
      safe_val <- function(col) {
        if (col %in% colnames(clone_df)) as.character(clone_df[[col]][row_idx]) else "N/A"
      }
      nodes <- data.frame(
        id = all_node_ids,
        label = all_node_ids,
        title = paste0(
          "Clone ID: ", all_node_ids, "<br>",
          "Subject: ", safe_val("subject"), "<br>",
          "Epitope: ", safe_val("epitope"), "<br>",
          "V-Alpha: ", safe_val("v_a_gene"), "<br>",
          "V-Beta: ", safe_val("v_b_gene"), "<br>",
          "CDR3-A: ", safe_val("cdr3_a_aa"), "<br>",
          "CDR3-B: ", safe_val("cdr3_b_aa"), "<br>",
          "Col by: ", color_var),
        group = as.character(clone_df[[color_var]][row_idx]),
        size = 20
      )
      
      max_edge_dist <- max(net_data$dist, na.rm = TRUE)
      edge_widths <- if (max_edge_dist > 0) 5 - (net_data$dist / max_edge_dist) * 4 else rep(3, nrow(net_data))
      
      # Create edges data frame with weights inversely proportional to distance
      # (closer clones = thicker edges)
      edges <- data.frame(
        from = net_data$node_1,
        to = net_data$node_2,
        width = edge_widths,  # Width between 1-5
        title = paste0(
          "Distance: ", round(net_data$dist, 2), "<br>",
          "Epitope 1: ", net_data$epitope_1, "<br>",
          "Epitope 2: ", net_data$epitope_2
        )
      ) %>% filter(from != to)
      
      # Render the network
      output$network_plot <- renderVisNetwork({
        visNetwork(nodes, edges) %>%
          visOptions(
            highlightNearest = TRUE,
            nodesIdSelection = TRUE,
            selectedBy = "group"
          ) %>%
          visLayout(randomSeed = 123) %>%
          visPhysics(
            enabled = TRUE,
            stabilization = list(iterations = 200)
          ) %>%
          visInteraction(navigationButtons = TRUE)
      })
      
      # Update statistics
      output$node_count <- renderText(paste("Nodes:", nrow(nodes)))
      output$edge_count <- renderText(paste("Edges:", nrow(edges)))
      
    }, error = function(e) {
      if (!inherits(e, "shiny.silent.error")) {
        message.error <- paste0("Error rendering network: ", e$message)
        render_error_ui(message.error, output)
        message(message.error)
      }
    })
  })
  
  
  #######################Plot Diversity plots###################################
  output$diversityPlot <- renderPlot({
    req(analysis_results())
    counts <- analysis_results()$clone_df$count
    # Generate a simple Lorenz Curve for Diversity
    df_plot <- data.frame(counts = sort(counts, decreasing = TRUE)) %>%
      mutate(
        cum_counts = cumsum(counts) / sum(counts),
        index = row_number() / n()
      )
    
    ggplot(df_plot, aes(x = index, y = cum_counts)) +
      geom_line(color = "steelblue", linewidth = 1.2) +
      geom_abline(lty = 2) +
      labs(title = "TCR Repertoire Evenness (Lorenz Curve)",
           x = "Fraction of Unique Clonotypes",
           y = "Cumulative Fraction of Total Template Count") +
      theme_minimal()
  })
  
  #####################kPCA plot########################################
  # Reactive distance matrix calculation for KPCA plot
  ######################################################################
  kpcaresults<-eventReactive(c(input$epitope_sel_kpca,input$run_kpca), {
    render_error_ui("", output = output) #Clear any error messages
    req(analysis_results())
    req(input$epitope_sel_kpca)
    results <- subset_tr_by_epitope(analysis_results(), input$epitope_sel_kpca)
    render_error_ui("", output)
    return(results)
  })
  
  observeEvent(c(input$run_kpca,input$kpca_facet,input$epitope_sel_kpca),{
    req(kpcaresults())
    tryCatch({
      tr <- kpcaresults()
      dist = tr[[paste0("pw_",input$chain_select_kpca)]]
      meta = as.data.frame(tr$clone_df)
      output$kpcaPlot <- renderPlotly({
        # renderPlotly()'s expression is evaluated lazily by Shiny, outside
        # this tryCatch's call stack, so it needs its own error handling to
        # stay consistent with the rest of the app.
        tryCatch({
          p<-generate_kpca_plot(dist_matrix = dist,
                             meta_df = meta,
                             color_by =input$colorby_sel_kpca,
                             title = paste0("kPCA of TCR distances for ", input$epitope_sel_kpca, " epitope"),
                             enable.facet = input$kpca_facet
                             )
          ggplotly(p) %>% layout(autosize = TRUE,legend = list(orientation = "h", y = -0.2)) %>% 
            config(responsive = TRUE)
        }, error = function(e) {
          if (!inherits(e, "shiny.silent.error")) {
            message.error <- paste0("Error rendering kPCA plot: ", e$message)
            render_error_ui(message.error, output)
            message(message.error)
          }
          NULL
        })
      })
    }, error = function(e)
    {
      if (!inherits(e, "shiny.silent.error"))
      {
        message(conda_path)
        message.error <- paste0("An error occurred: ",
                                e$message)
        render_error_ui(message.error, output)
        message(message.error)
      }
    })
  })
  
  
  observeEvent(input$chain_select_kpca,{
    output$kpcaPlot <- renderPlotly({
      return(NULL)
    })
  })
  
  #Update distance matrix options based on selected chain when updating heatmap
  observeEvent(input$chain_select_heatmap,{
    req(analysis_results())
    tr<-analysis_results()
    #Get chains from the distance compute
    all.tr_names <- names(tr)
    message(paste0(all.tr_names, collapse = "\n"))
    chain.selected<-input$chain_select_heatmap
    chains_tr <- all.tr_names[grepl("^pw", all.tr_names)]
    chain.start<-substr(chain.selected,1,1)
    pattern <- paste0("^[^_]+_([",chain.start,"cp]([^_]+_",chain.start,"|[^_]*$))")
    chains_tr <- chains_tr[grepl(pattern,chains_tr)]
    message(chains_tr)
    #Include a select option for the chains
    output$matrixoptions2 <- renderUI({
      tagList(
        selectInput(
          "matrix_select",
          "Select distance matrix:",
          selected = chains_tr[1],
          choices = chains_tr
        )
      )
    })
  })
  
  #######################Heatmap Tab###########################
  #Update heatmap and matrix table on the following events
  #############################################################
  heat_data<-eventReactive(input$run_heatmap, {
    render_error_ui("", output = output) #Clear any error messages
    tryCatch({
      req(analysis_results())
      # Distances were already computed in analysis_results(); just subset
      # to the selected epitope rather than re-running tcrdist3.
      results <- subset_tr_by_epitope(analysis_results(), input$epitope_sel_heatmap)
      render_error_ui("", output)
      chain2check <- gsub("pw_", "", (input$matrix_select))
      if (sum(chain2check %in% colnames(results$clone_df)) == 0)
        chain2check <- paste0("cdr3_", substring(chain2check, 1, 1), "_aa")
      message(chain2check)
      pw_matrix <- results[[input$matrix_select]]
      clone_df = results$clone_df
      mat <- pw_matrix
      rownames(mat) <- clone_df[[chain2check]]
      colnames(mat) <- clone_df[[chain2check]]
      # Display the distance matrix
      session.id <- session$token
      output$dist_matrix <- DT::renderDataTable({
        # renderDataTable()'s expression is evaluated lazily by Shiny,
        # outside this tryCatch's call stack, so it needs its own error
        # handling to stay consistent with the rest of the app.
        tryCatch({
          DT::datatable(
            mat,
            extensions = 'Buttons',
            options = list(
              scrollX = TRUE,
              pageLength = 10,
              rownames = F,
              dom = 'lfrtipB',
              #'flrtBip',
              buttons  = list(
                'copy',
                list(
                  extend = 'csv',
                  text = "csv",
                  filename = paste0(
                    'TCR-distance-matrix-',
                    input$matrix_select,
                    "-",
                    session.id
                  )
                ),
                list(
                  extend = "excel",
                  text = 'excel',
                  filename = paste0(
                    'TCR-distance-matrix-',
                    input$matrix_select,
                    "-",
                    session.id
                  )
                ),
                list(
                  extend = 'pdf',
                  text = 'pdf',
                  filename = paste0(
                    'TCR-distance-matrix-',
                    input$matrix_select,
                    "-",
                    session.id
                  )
                ),
                list(
                  extend = 'print',
                  text = "print",
                  filename = paste0(
                    'TCR-distance-matrix-',
                    input$matrix_select,
                    "-",
                    session.id
                  )
                )
              )
            )
          )
        }, error = function(e) {
          if (!inherits(e, "shiny.silent.error")) {
            message.error <- paste0("Error rendering distance matrix table: ", e$message)
            render_error_ui(message.error, output)
            message(message.error)
          }
          NULL
        })
      }, server = T)
      return(mat)
      # Display the heatmap
      # output$heatmap <- renderPlot({
      #   pheatmap(
      #     mat,
      #     clustering_distance_rows = "euclidean",
      #     clustering_distance_cols = "euclidean",
      #     clustering_method = "ward.D2",
      #     main = "TCR Distance Heatmap"
      #   )
      #})
    }, error = function(e)
    {
      if (!inherits(e, "shiny.silent.error"))
      {
        message.error <- paste0("Error rendering heatmap: ", e$message)
        render_error_ui(message.error, output)
      }
      
    })
  })
  
  observeEvent(c(input$run_heatmap,input$show_row_names,input$show_col_names), {
    mat<-heat_data()
    req(mat)
    ht <- draw(Heatmap(mat,
                  name = "TCRdist",
                  column_title = paste("Clonotypes for", input$epitope_sel_heatmap),
                  show_row_dend = TRUE,
                  show_column_dend = TRUE,
                  show_row_names = input$show_row_names,
                  show_column_names = input$show_col_names))
    # Make it interactive in Shiny
    makeInteractiveComplexHeatmap(
      input, output, session, 
      ht, 
      heatmap_id = "heatmap"
    )
    reset_trigger(0) # Reset the trigger after rendering the heatmap
  })
  

  ######################clone tab##############################################
  #Load clonal information into a tab
  #W###########################################################################
  observe({
    req(analysis_results())
    clone_df<-analysis_results()$clone_df
    
    # For each epitope (and chain), find the CDR3 sequence shared by the
    # largest number of distinct clones -- i.e. the most "public"/common
    # clonotype for that epitope.
    output$topCdr3Plot <- renderPlot({
      tryCatch({
        cdr3_cols <- grep("^cdr3_[ab]_aa$", names(clone_df), value = TRUE)
        req(length(cdr3_cols) > 0)
        
        top_cdr3 <- clone_df %>%
          select(epitope, clone_id, all_of(cdr3_cols)) %>%
          pivot_longer(cols = all_of(cdr3_cols), names_to = "chain", values_to = "cdr3") %>%
          filter(!is.na(cdr3), cdr3 != "") %>%
          mutate(chain = ifelse(grepl("_a_", chain), "alpha", "beta")) %>%
          summarise(n_clones = n_distinct(clone_id), .by = c(epitope, chain, cdr3)) %>%
          slice_max(n_clones, n = 1, with_ties = FALSE, by = c(epitope, chain))
        req(nrow(top_cdr3) > 0)
        
        ggplot(top_cdr3, aes(x = epitope, y = n_clones, fill = chain)) +
          geom_col(position = position_dodge(width = 0.7), width = 0.6) +
          geom_text(aes(label = cdr3),
                    position = position_dodge(width = 0.7),
                    vjust = -0.4, size = 3, fontface = "bold") +
          scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
          labs(title = "Most Common CDR3 per Epitope",
               subtitle = "Bar height = number of clones sharing that CDR3; label = sequence",
               x = "Epitope", y = "Number of Clones", fill = "Chain") +
          theme_minimal(base_size = 13) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
            panel.grid.minor = element_blank()
          )
      }, error = function(e) {
        if (!inherits(e, "shiny.silent.error")) {
          message.error <- paste0("Error rendering top CDR3 plot: ", e$message)
          render_error_ui(message.error, output)
          message(message.error)
        }
        NULL
      })
    })
    
    output$clone_matrix <- DT::renderDataTable({
      DT::datatable(
        clone_df,
        caption = "Summary of TCRdist clone df",
        extensions = 'Buttons',
        rownames = F,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'flrtip',
          buttons = list(
            list(
              extend = 'csv',
              text = 'Download CSV',
              filename = paste0("TCRDist3-clone-", session$token)
            ),
            list(
              extend = 'excel',
              text = 'Download Excel',
              filename = paste0("TCRDist3-clone-", session$token)
            )
          )
        )
      )
    }, server = T)
  })
  
  #Reset heatmap window
  observeEvent(c(input$matrix_select,input$epitope_sel_heatmap,input$chain_select_heatmap),{
    output$dist_matrix <- DT::renderDataTable({NULL})
  })
  

  # 
  # observeEvent(c(input$epitope_sel_heatmap,input$chain_select_heatmap), {
  #   output$heatmap_container <- renderUI({
  #       return(p("Updating values... Please click 'Load Heatmap' to update the heatmap and distance matrix."))
  #   })
  # })
  
  
  
  #Add sankey gen button to the dashboard
  observeEvent(c(input$epitope_sel, input$columns_sel), {
    render_error_ui("", output = output) #Clear any error messages
    req(input$epitope_sel)
    output$btn_gensankey_ui <- renderUI(actionButton(
      "btn_gensankey",
      HTML(
        "<span class='glyphicon glyphicon-play-circle'></span> Generate Sankey plot"
      )
    ))
  })
  
  #Display summary stats
  observeEvent(input$run, {
    render_error_ui("", output = output) #Clear any error messages
    tr <- analysis_results()
    tryCatch({
      req(tr)
      clone_df <- py_to_r(tr$clone_df)
      summary.stats <- clone_df |>
        summarise(
          num_individuals = n_distinct(subject),
          num_clones = n_distinct(clone_id),
          .by = epitope
        )
      
      output$summary_matrix <- DT::renderDataTable({
        # renderDataTable()'s expression is evaluated lazily by Shiny,
        # outside this tryCatch's call stack, so it needs its own error
        # handling to stay consistent with the rest of the app.
        tryCatch({
          DT::datatable(summary.stats, options = list(scrollX = TRUE, pageLength = 10))
        }, error = function(e) {
          if (!inherits(e, "shiny.silent.error")) {
            message.error <- paste0("Error rendering summary table: ", e$message)
            render_error_ui(message.error, output)
            message(message.error)
          }
          NULL
        })
      })
      
      output$summaryPlot <- renderPlot({
        # renderPlot()'s expression is evaluated lazily by Shiny, outside
        # the surrounding tryCatch's call stack, so it needs its own error
        # handling to stay consistent with the rest of the app.
        tryCatch({
          ggplot(summary.stats, aes(x = epitope, y = num_clones, fill = epitope)) +
            geom_col(width = 0.6, show.legend = FALSE) +
            geom_text(aes(label = num_clones), vjust = -0.4, size = 4, fontface = "bold") +
            scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
            labs(title = "Number of Clones per Epitope",
                 x = "Epitope", y = "Number of Clones") +
            theme_minimal(base_size = 13) +
            theme(
              plot.title = element_text(face = "bold", hjust = 0.5),
              axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
              panel.grid.minor = element_blank()
            )
        }, error = function(e) {
          if (!inherits(e, "shiny.silent.error")) {
            message.error <- paste0("Error rendering summary plot: ", e$message)
            render_error_ui(message.error, output)
            message(message.error)
          }
          NULL
        })
      })
      
      clone_plt_data<- clone_df %>% 
        select("epitope",matches("cdr3"),"clone_id","subject") %>% 
        pivot_longer(cols = matches("cdr3"),
                            names_to = c(".value", "chain"),
                            names_sep = "_")
      
      
    }, error = function(e)
    {
      if (!inherits(e, "shiny.silent.error"))
      {
        message.error <- paste0("Error loading summary results: ", e$message)
        output$genepairsplots <- renderUI(NULL)
        render_error_ui(message = message.error, output)
        message(message.error)
      }
    })
  })
  
  #Generate sankey plot
  observeEvent(c(input$epitope_sel, input$btn_gensankey, input$columns_sel),
               {
                 render_error_ui("", output = output) #Clear any error messages
                 tryCatch({
                   tr <- analysis_results()
                   req(tr)
                   req(input$columns_sel)
                   req(length(input$columns_sel) > 1)
                   render_error_ui(NULL, output)
                   plotting <- import("tcrdist.plotting")
                   svg_PA  = plotting$plot_pairings(tr$clone_df[tr$clone_df$epitope == input$epitope_sel, ],
                                                    cols = r_to_py(input$columns_sel),
                                                    count_col = 'count')
                   output$genepairsplots <- renderUI(HTML(svg_PA))
                 }, error = function(e)
                 {
                   if (!inherits(e, "shiny.silent.error"))
                   {
                     message.error <- paste0("Error rendering sankey plot: ", e$message)
                     output$genepairsplots <- renderUI(NULL)
                     render_error_ui(message = message.error, output)
                     #output$errorreport <- renderText(message.error)
                     message(message.error)
                   }
                 })
               })
  
  epitope_dist_mat_KNN<-eventReactive(input$epitope_sel_n, {
    render_error_ui("", output = output) #Clear any error messages
    tr<-analysis_results()
    req(tr)
    epitope_knn_sel<-input$epitope_sel_n
    #Check if only single epitope is selected then return all
    
    if(length(epitope_knn_sel)==1)
    {
      return(subset_tr_by_epitope(tr, NULL))
    }
    else{
      return(subset_tr_by_epitope(tr,epitope_knn_sel))
    }
  })
  
  #Run Fixed Radius Neighbourhoods
  observeEvent(c(input$run_fixn), {
    render_error_ui("", output = output) #Clear any error messages
    tryCatch({
      withProgress(
        message = paste0("Run KNN with ", paste0(input$epitope_sel_n, collapse = ", "), ":"),
        value = 0,
        {
          tr <- epitope_dist_mat_KNN()
          req(tr)
          req(input$epitope_sel_n)
          message("Import tcrdist.rep_diff")
          rep_diff <- import("tcrdist.rep_diff")
          incProgress(0.1, "Preparing data...")
          np <- import("numpy")
          chain <- input$chain_select_knn
          message("Computing for chain", chain)
          # diff testing is based on binary comparison, so all epitopes not matching epi_1 are set to 'Other' or epitope 2
          epitope_n <- epi_1 <- input$epitope_sel_n
          if(length(epi_1)>1)
            epi_1<-epi_1[1]
          
          message("Selected epitope is ",epi_1,"\n")
          clone_df <- as.data.frame(tr$clone_df)
          rownames(clone_df) <- NULL
          
          if(length(epitope_n)==1)
          {
            # Code to handle multiple epitopes
            clone_df[, epi_1] <- ifelse(
              clone_df$epitope == epi_1, epi_1, "Other")
            epitope_n <- list(epitope_n)
          }
          else{
            clone_df[, epi_1] <- ifelse(
              clone_df$epitope == epi_1, epi_1, epitope_n[2])
          }
          
          if(length(chain)==1)
          {
            pwmat <- tr[[paste0("pw_", chain)]]
          }
          else if(length(chain) > 1)
          {
            pwmat <- tr[[paste0("pw_", chain[1])]]+tr[[paste0("pw_", chain[2])]]
          }

          incProgress(0.5, "Running neighbourhood differential...")
          # Larger Radius
          nn_df = rep_diff$neighborhood_diff(
            clone_df = clone_df,
            pwmat = pwmat,
            count_col = 'count',
            x_cols = r_to_py(list(c(epi_1))),
            knn_radius = r_to_py(input$knn_radius)
          )
          nn_df1 <- nn_df %>% select(
            c(
              'K_neighbors',
              'val_0',
              'ct_0',
              'val_2',
              'ct_2',
              'RR',
              'OR',
              'pvalue',
              'FWERp',
              'FDRq'
            )
          ) %>% arrange(FDRq, desc(OR), desc(ct_0))
          output$fixedradius_output_ui <- renderUI({
            withSpinner(DT::dataTableOutput("fixedradius_output"))
          })
          incProgress(0.9, "Displaying results table...")
          
          session.id <- session$token
          file.name = paste0(
            'TCRdist3-KNN-clustering-',
            paste0(input$epitope_sel_n, collapse = "_"),
            "-",
            paste0(chain, collapse = "_"),
            "-",
            session.id
          )
          output$fixedradius_output <- DT::renderDT({
            # renderDT()'s expression is evaluated lazily by Shiny, outside
            # this tryCatch's call stack, so it needs its own error
            # handling to stay consistent with the rest of the app.
            tryCatch({
              DT::datatable(
                nn_df1,
                caption = "KNN neighbours",
                rownames = F,
                extensions = 'Buttons',
                options = list(
                  scrollX = TRUE,
                  pageLength = 10,
                  dom = 'lfrtipB',
                  buttons  = list(
                    'copy',
                    list(
                      extend = 'csv',
                      text = "csv",
                      filename = file.name
                    ),
                    list(
                      extend = "excel",
                      text = 'excel',
                      filename = file.name
                    ),
                    list(
                      extend = 'pdf',
                      text = 'pdf',
                      filename = file.name
                    ),
                    list(
                      extend = 'print',
                      text = "print",
                      filename = file.name
                    )
                  )
                )
              )
            }, error = function(e) {
              if (!inherits(e, "shiny.silent.error")) {
                message.error <- paste0("Error rendering fixed radius results table: ", e$message)
                render_error_ui(message.error, output)
                message(message.error)
              }
              NULL
            })
          }, server = F)
        }
      )
    }, error = function(e)
    {
      if (!inherits(e, "shiny.silent.error"))
      {
        message.error <- paste0("Error computing fixed radius neighbourhoods: ",
                                e$message)
        output$genepairsplots <- renderUI(NULL)
        render_error_ui(message = message.error, output)
        #output$errorreport <- renderText(message.error)
        message(message.error)
      }
    })
  })
  
  epitope_results_hierarch <- eventReactive(input$epitope_sel_h, {
    render_error_ui("", output = output) #Clear any error messages
    req(tcr_data())
    epitope_h_sel<-input$epitope_sel_h
    tr<-analysis_results()
    req(tr)
    if(length(epitope_h_sel)==1)
    {
      return(subset_tr_by_epitope(tr, NULL))
    }
    else{
      return(subset_tr_by_epitope(tr, epitope_h_sel))
    }
  })
  
  #Listen to hierarchical neighbourhood clustering events
  #neighbour_results <- eventReactive(
  observeEvent(input$run_Hierarch, {
    render_error_ui("", output = output) #Clear any error messages
    hierarch_tabs_check(TRUE)
    withProgress(message = 'Running hierarchical neighbourhood clustering...', value = 0, {
      tryCatch({
        tr.this<-epitope_results_hierarch()
        req(tr.this)
        req(input$epitope_sel_h)
        message("Running neighbouring")
        epitope_h <-epi_1<-input$epitope_sel_h
        if(length(epi_1)>1)
          epi_1<-epi_1[1]
        message("Selected epitope is ",epi_1,"\n")
        clone_df <- as.data.frame(tr.this$clone_df)
        rownames(clone_df) <- NULL
        
        if(length(epitope_h)==1)
        {
          # Code to handle multiple epitopes
          clone_df[, epi_1] <- ifelse(
            clone_df$epitope == epi_1, epi_1, "Other")
          epitope_h <- list(epitope_h)
        }
        else{
          clone_df[, epi_1] <- ifelse(
            clone_df$epitope == epi_1, epi_1, epitope_h[2])
        }
        rep_diff <- import("tcrdist.rep_diff")
        hierdiff <- import("hierdiff")
        pd <- import("pandas")
        incProgress(0.1, "Run hierarchical neighbourhoods")
        chain <- input$chain_select_neighbours
        
        if(length(chain)==1)
        {
          pwmat <- tr.this[[paste0("pw_", chain)]]
        }
        else if(length(chain) > 1)
        {
          pwmat <- tr.this[[paste0("pw_", chain[1])]]+tr.this[[paste0("pw_", chain[2])]]
        }
        
        res_Z <- rep_diff$hcluster_diff(clone_df, pwmat, x_cols = r_to_py(as.list(c(epi_1))), count_col = 'count')
        res <- res_Z[[1]]
        Z <- res_Z[[2]]
        
        res_summary <- rep_diff$member_summ(
          res_df = res,
          clone_df = clone_df,
          addl_cols = r_to_py(list('epitope'))
        )
        res_detailed = r_to_py(cbind(res, res_summary))
        incProgress(0.3, detail = "Preparing results...")
        res_detailed.df <- py_to_r(res_detailed)
        columns.tooltip <- colnames(res_detailed.df)[(grep("FDRq", colnames(res_detailed.df)) +
                                                        1):ncol(res_detailed.df)]
        message(paste0(columns.tooltip, ";"))
        incProgress(0.5, "Generating cluster diagram...")
        html = hierdiff$plot_hclust_props(
          Z,
          title ='',
          res = res_detailed,
          tooltip_cols = r_to_py(columns.tooltip),
          alpha = 0.00001,
          colors = r_to_py(c('orange', 'lightgray')),
          alpha_col = 'pvalue'
        )
        incProgress(0.9, detail = "Rendering results table...")
        session.id <- session$token
        file.name = paste0('TCRdist3-hierarchical-clustering-',
                           paste0(input$epitope_sel_h,collapse = "_"),
                           "-",
                           paste0(chain,collapse = "_"),
                           "-",
                           session.id)
        output$hierachy_output <- DT::renderDataTable({
          # renderDataTable()'s expression is evaluated lazily by Shiny,
          # outside this tryCatch's call stack, so it needs its own error
          # handling to stay consistent with the rest of the app.
          tryCatch({
            DT::datatable(
              py_to_r(res_detailed),
              caption = "Hierarchical neighbours",
              extensions = 'Buttons',
              rownames = F,
              options = list(
                scrollX = TRUE,
                pageLength = 10,
                dom = 'lfrtipB',
                #'flrtBip',
                buttons  = list(
                  'copy',
                  list(
                    extend = 'csv',
                    text = "csv",
                    filename = file.name
                  ),
                  list(
                    extend = "excel",
                    text = 'excel',
                    filename = file.name
                  ),
                  list(
                    extend = 'pdf',
                    text = 'pdf',
                    filename = file.name
                  ),
                  list(
                    extend = 'print',
                    text = "print",
                    filename = file.name
                  )
                )
              )
            )
          }, error = function(e) {
            if (!inherits(e, "shiny.silent.error")) {
              message.error <- paste0("Error rendering hierarchical results table: ", e$message)
              render_error_ui(message.error, output)
              message(message.error)
            }
            NULL
          })
        }, server = F)
        
        #message(html)
        incProgress(0.95, detail = "Save html file...")
        temp_path <- tempfile()
        remove_www_tempfiles(session, "hierachical.html")
        html.file <- paste0(session$token, "_hierachical.html")
        write(html, temp_path)
        file.copy(temp_path, file.path("www", html.file), overwrite = TRUE)
        
        output$hierachy_output_plot <- renderUI(
          tags$iframe(
            src = html.file,
            width = "100%",
            style = "height: 100vh;",
            # Use CSS for full viewport height
            frameBorder = "0"
          )
        )
        neighbour_results <<- res_detailed.df
      }, error = function(e)
      {
        if (!inherits(e, "shiny.silent.error"))
        {
          message.error <- paste0("Error computing hierarchical neighbourhood clusters: ",
                                  e$message)
          output$genepairsplots <- renderUI(NULL)
          render_error_ui(message = message.error, output)
          #output$errorreport <- renderText(message.error)
          message(message.error)
        }
      })
    })
  })
  
  
  #Generate Network Graph options when chain selection changes
  observeEvent(input$chain_select_neighbours, {
    tr <- analysis_results()
    
    tryCatch({
      req(tr)
      chain <- input$chain_select_neighbours
      req(length(chain) > 0)
      first_chain <- chain[1]
      genes <- tr$clone_df[[paste0("v_", substr(first_chain, 1, 1), "_gene")]]
      chain2check <- paste0("cdr3_", substring(first_chain, 1, 1), "_aa")
      pw_matrix <- tr[[paste0("pw_", chain2check)]]
      req(!is.null(pw_matrix))
      updateSelectizeInput(
        session,
        "v_gene_filter",
        choices = c("All", genes),
        server = TRUE
      ) # This is the magic argument
      output$dist_threshold_ui <- renderUI({
        tagList(
          sliderInput(
            "dist_threshold",
            "TCR distance:",
            min = floor(min(pw_matrix)),
            max = ceiling(max(pw_matrix)),
            value = ceiling(median(pw_matrix))
          )
        )
      })
    }, error = function(e) {
      if (!inherits(e, "shiny.silent.error"))
      {
        message.error <- paste0("Error updating TRVs: ", e$message)
        output$genepairsplots <- renderUI(NULL)
        render_error_ui(message = message.error, output)
        message(message.error)
      }
    })
    
  })
  
  
  # Reactive: Filter Data based on UI inputs
  filtered_data <- reactive({
    tr <- analysis_results()
    nb <- neighbour_results
    tryCatch({
      req(tr)
      req(nb)
      # Generate synthetic metadata for 50 TCR clones
      tcr_metadata <- as.data.frame(nb)
      # Generate synthetic KNN results (Edge List)
      # In reality, this comes from your tcrdist3 'neighbors' output
      chain.char <- substr(input$chain_select_neighbours, 1, 1)
      #chain2check <- paste0("cdr3_",chain.char, "_aa")
      pw_matrix = tr[[paste0("pw_", "cdr3_", chain.char, "_aa")]]
      message("Chain starting character:", chain.char)
      knn_edges <- tcr_metadata[, c(
        "cid",
        paste0("cdr3_", chain.char, "_aa"),
        "epitope",
        paste0("v_", chain.char, "_gene"),
        "neighbors"
      )] %>% unnest(neighbors)
      colnames(knn_edges) <- c("from", "cdr3_aa", "epitope", "v_gene", "to")
      rownames(knn_edges) <- NULL
      knn_edges <- knn_edges[!duplicated(knn_edges), ]
      knn_edges <- knn_edges[order(knn_edges$from, decreasing = F), ]
      tcr.data <- tcr_data()
      clone.df <- tr$clone_df
      knn_edges$distance <- as.numeric(apply(knn_edges, 1, function(x)
        pw_matrix[x[["from"]], x[["to"]]]))
      tcr_metadata <- tcr_metadata[, c(
        "cid",
        paste0("cdr3_", chain.char, "_aa"),
        paste0("v_", chain.char, "_gene"),
        "epitope"
      )]
      colnames(tcr_metadata) <- c("clone_id", "cdr3_aa", "v_gene", "epitope")
      nodes <- tcr_metadata
      edges <- knn_edges
      
      # Filter nodes by V-Gene if selected
      if (input$v_gene_filter != "All") {
        valid_clones <- nodes %>% filter(grepl(input$v_gene_filter, v_gene)) %>% pull(clone_id)
        
        # Keep nodes matching filter OR nodes connected to them
        edges <- edges %>%
          filter(from %in% valid_clones | to %in% valid_clones)
        
        active_nodes <- unique(c(edges$from, edges$to))
        nodes <- nodes %>% filter(clone_id %in% active_nodes)
      }
      
      # Filter edges by distance
      edges <- edges %>% filter(distance <= input$dist_threshold)
      list(nodes = nodes, edges = edges)
      
    }, error = function(e) {
      if (!inherits(e, "shiny.silent.error"))
      {
        message.error <- paste0("Error updating TRVs", e$message)
        output$genepairsplots <- renderUI(NULL)
        render_error_ui(message = message.error, output)
        message(message.error)
      }
    })
  })
  
  # Render Network
  output$knnGraph <- renderVisNetwork({
    data <- filtered_data()
    
    # Prepare Nodes for visNetwork
    nodes_viz <- data$nodes %>%
      rename(id = clone_id) %>%
      mutate(
        # Dynamic coloring based on selection
        group = .data[[input$epitope_color]],
        # Tooltip info
        title = paste0(
          "<p><b>ID:</b> ",
          id,
          "<br>",
          "<b>CDR3:</b> ",
          cdr3_aa,
          "<br>",
          "<b>V-Gene:</b> ",
          v_gene,
          "</p>"
        ),
        label = id # Or use NA to hide labels if too crowded
      )
    
    # Prepare Edges for visNetwork
    edges_viz <- data$edges %>%
      mutate(
        width = (60 - distance) / 5,
        # Thicker lines for closer neighbors
        title = paste("Dist:", round(distance, 2))
      )
    
    visNetwork(nodes_viz, edges_viz) %>%
      visOptions(highlightNearest = TRUE,
                 nodesIdSelection = TRUE) %>%
      visLayout(randomSeed = 123) %>%
      visPhysics(stabilization = FALSE) %>%
      visInteraction(navigationButtons = TRUE) %>%
      visLegend()
  })
  output$showtrees <- reactive({
    return("donotshow")
  })
  outputOptions(output, "showtrees", suspendWhenHidden = FALSE)
  
  ##############################################################
  #.          Tab used to generate trees                       #
  ##############################################################
  
  observeEvent(input$run_trees, {
    render_error_ui("", output = output) #Clear any error messages
    df_tcr <- tcr_data()
    req(df_tcr)
    withProgress(message = 'Running trees...', value = 0, {
      tryCatch({
        #If more than one chain has been selected, then recompute distance matrix based on chain
        chain_trees <- input$chain_select_tree
        output$showtrees <- reactive({
          return("show")
        })
        epitopes_trees <- input$epitope_sel_trees
        if (length(chain_trees) == length(chains()) && length(unique(df_tcr$epitope)) == length(epitopes_trees))
        {
          message("Using already computed results")
          tcr <- analysis_results()
        }
        else{
          incProgress(0.1, detail = "Computing distance matrix..")
          tcr <- get_epitope_dist_matrix(
            df_tcr,
            epitope_name = epitopes_trees,
            organism = input$organism,
            chains = chain_trees,
            dbfile = input$dbfile
          )  
        }
        outputOptions(output, "showtrees", suspendWhenHidden = FALSE)
        output$tree_output_plot <- renderUI(NULL)
        incProgress(0.3, detail = "Importing required library...")
        TCRtree <- import("tcrdist.tree")
        incProgress(0.4, detail = "Preparing data...")
        file.name <- paste0(session$token, "_tree.html")
        temp.file <- tempfile()
        message("Temp file for trees: ", temp.file)
        remove_www_tempfiles(session, "tree.html")
        message("Running trees")
        tcrtree = TCRtree$TCRtree(tcrrep = tcr, html_name = temp.file)
        incProgress(0.5, detail = "Building tree...")
        tcrtree$build_tree()
        incProgress(0.9, detail = "Finalizing...")
        file.copy(temp.file, file.path("www", file.name), overwrite = TRUE)
        output$tree_output_plot <- renderUI({
          withSpinner(tags$iframe(
            src = file.name,
            width = "95%",
            style = "height: 90vh;",
            # Use CSS for full view port height
            frameBorder = "0"
          ))
        })
      }, error = function(e)
      {
        if (!inherits(e, "shiny.silent.error"))
        {
          message.error <- paste0("Error generating tree: ", e$message)
          output$genepairsplots <- renderUI(NULL)
          render_error_ui(message = message.error, output)
          #output$errorreport <- renderText(message.error)
          message(message.error)
        }
      })
    })
  })
  
  # Download handler for the download button
  output$download_data <- downloadHandler(
    filename = function() {
      paste("my-clone-data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      tr = analysis_results()
      req(tr)
      write.csv(tr$clone_df, file, row.names = FALSE)
    }
  )
  
  #Clean-up temp files when session ends
  onStop(function() {
    remove_www_tempfiles(session)
  })
  
  #####Test diversity plot
  
  # Reactive data processing
  processed_data <- reactive({
    req(tcr_data())
    tr <- analysis_results()
    req(tr)
    get_repertoir_diversity(tr, chains())
  })
  
  # Top Plot: 1 - Simpson's (Clonality)
  output$simpsonPlot <- renderPlot({
    data <- processed_data()
    ggplot(data.frame(val = data$one_minus_simpson, label = data$epitope), 
           aes(x = label, y = val, color = label)) +
      geom_point(size = 4, show.legend = FALSE) +
      ylim(0, 1) +
      labs(title = "1-Simpson's Diversity (Clonality)", x = "Epitope", y = "1 - Simpson") +
      theme_minimal()
  })
  
  # Bottom Plot: Inverse Simpson's (Diversity)
  output$invSimpsonPlot <- renderPlot({
    data <- processed_data()
    ggplot(data.frame(val = data$inv_simpson, label = data$epitope), 
           aes(x = label, y = val, fill = label)) +
      geom_col(show.legend = FALSE) +
      scale_y_log10() +
      labs(title = "Inverse Simpson's Diversity", x = "Epitope", y = "Index (Log Scale)") +
      theme_minimal()
  })
  
  
  
  ##############################################################################
  #                     Motif Enrichment Tab                                   #
  ##############################################################################
  
  
  output$cluster_selector <- renderUI({
    req(neighbors(), cluster_order())
    ord   <- cluster_order()
    sizes <- sapply(neighbors()[ord], length)
    
    choice_labels <- paste0("Cluster ", ord, " (n=", sizes, ")")
    choices <- setNames(ord, choice_labels)
    
    selectInput("cluster_id", "Choose cluster (largest first)",
                choices = choices)
  })
  motif_filtered_tcr<-reactive({
    tr <- analysis_results()
    req(tr)
    chain_motif <- input$chain_select_motif
    epitopes_motif <- input$epitope_sel_motif
    tryCatch({
      tr <- subset_tr_by_epitope(tr, epitopes_motif, chain_motif)
      req(tr)
      return(tr)
    }, error = function(e) {
      print(e$message)
      showModal(modalDialog(
        title = "Error",
        message = paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  
  motif_pwmat <- reactive({
    chain_motif <- input$chain_select_motif
    req(chain_motif)
    tryCatch({
      tr <- motif_filtered_tcr()
      req(tr)
      message("Chain selected:", chain_motif)
      # Get pairwise matrix dynamically for alpha, beta, gamma, or delta
      dist_mat <- tr[[paste0("pw_", chain_motif)]]
      req(!is.null(dist_mat))
      return(dist_mat)
    }, error = function(e) {
      print(e$message)
      showModal(modalDialog(
        title = "Error",
        message = paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  observeEvent(c(input$chain_select_motif, input$epitope_sel_motif), {
    req(motif_pwmat())
    dist_mat <- motif_pwmat()
    output$motif_radius_ui <- renderUI({
      sliderInput("motif_radius",
                  "TCR distance:",
                  min = 0,
                  max = max(dist_mat, na.rm = TRUE),
                  value = min(50L, max(dist_mat, na.rm = TRUE)),
                  step = 1)
    })
  })
  
  neighbors<-reactive({
    req(motif_pwmat(), !is.null(input$motif_radius))
    dist_mat <- motif_pwmat()
    radius <- input$motif_radius
    # Find neighborhoods (clusters) around a radius (1-indexed for R subsetting)
    neighbors <- lapply(seq_len(nrow(dist_mat)), function(i) {
      which(dist_mat[i, ] <= radius)
    })
    return(neighbors)
  })
  
  # ---- Deduplicated neighbor list -----------------------------------
  neighbors_dedup <- reactive({
    req(neighbors())
    nb <- neighbors()
    # build a signature for each cluster: sorted, comma-joined member indices
    # (two clusters are "the same" if they contain exactly the same members)
    signatures <- sapply(nb, function(idx) paste(sort(idx), collapse = ","))
    
    # keep only the first occurrence of each unique signature
    keep <- !duplicated(signatures)
    kept_indices <- which(keep)
    setNames(kept_indices, kept_indices)
  })
  
  # ---- Cluster ordering: largest first, using deduplicated list -----
  cluster_order <- reactive({
    req(neighbors(), neighbors_dedup())
    keep_idx <- neighbors_dedup()
    sizes <- sapply(neighbors()[keep_idx], length)
    keep_idx[order(sizes, decreasing = TRUE)]
  })
  
  # cluster_order <- reactive({
  #   req(neighbors())
  #   sizes <- sapply(neighbors(), length)
  #   # order of original cluster indices, largest cluster first
  #   order(sizes, decreasing = TRUE)
  # })
  
  observeEvent(c(input$chain_select_motif,input$epitope_sel_motif,input$motif_radius),{
    # Code to run motif enrichment analysis
    print("Running chain event")
    req(input$chain_select_motif)
    req(input$epitope_sel_motif)
    tryCatch({
      req(neighbors())
      neighbors<-neighbors()
      output$cluster_select_ui <- renderUI({
        req(neighbors(), cluster_order())
        ord   <- cluster_order()
        sizes <- sapply(neighbors(), length)[ord]
        
        # named choices: label shows rank + size, value is the ORIGINAL cluster index
        choice_labels <- paste0("Cluster ", ord, " (n=", sizes, ")")
        choices <- setNames(ord, choice_labels)
        
        selectizeInput("cluster_id","Choose cluster (largest first)",
                       choices = choices)
      })
      # output$cluster_select_ui <- renderUI({
      #   selectizeInput("cluster_id", "Choose cluster (by seed index)",
      #               choices = seq_along(neighbors))
      # })
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
    })
  
  # ---- Subset of the original data belonging to the selected cluster ----
  cluster_data <- reactive({
    req(neighbors(),motif_filtered_tcr(), input$cluster_id)
    idx <- neighbors()[[as.integer(input$cluster_id)]]
    motif_filtered_tcr()$clone_df[idx, , drop = FALSE]
  })
  
 observeEvent(c(input$run_motifs, input$cluster_id), {
      req(neighbors(), motif_filtered_tcr(), input$cluster_id, input$chain_select_motif)
      idx <- neighbors()[[as.integer(input$cluster_id)]]
      palmotif <- import("palmotif")
      chain_motif <- input$chain_select_motif
      cdr3 <- paste0("cdr3_", substr(chain_motif, 1, 1), "_aa")
      cluster_seqs <- motif_filtered_tcr()$clone_df[[cdr3]][idx]
      if (length(cluster_seqs) < 2) return(NULL)
      # palmotif needs a DataFrame with sequences (and optionally counts/weights)
      motif_input <- data.frame(seq = cluster_seqs, stringsAsFactors = FALSE)
      py_motif_input <- r_to_py(motif_input)
      motif_result <- palmotif$compute_pal_motif(
        seqs = py_motif_input$seq$tolist(),
        centroid = cluster_seqs[1]   # reference/anchor sequence for alignment
      )

    # ---- Render the SVG logo directly in the app -----------------
  prob_matrix<-py_to_r(motif_result[[1]]) %>% as.matrix() %>% t()
  output$motif_plot<-renderPlot({
      ggseqlogo(prob_matrix, method = "bits", seq_type = "aa") +
        ggtitle(paste0(cdr3," Motif")) 
  })
  
  
    # output$motif_svg <- renderUI({
    #   req(motif_result)
    #   svg_string <- palmotif$svg_logo(motif_result, return_str = TRUE) #SVG text returned by palmotif
    #   tagList(
    #     tags$div(
    #       style = "border:1px solid #ddd; padding:10px; background:white;",
    #       HTML(svg_string)
    #     )
    #   )
    # })
    
    output$motif_table <- renderTable({
      req(motif_result)
      py_to_r(motif_result[[1]])
    })
  })
 
 
 # ---- Motif Sidebar quick summary --------------------------------------
 output$cluster_summary <- renderUI({
   req(cluster_data(), input$chain_select_motif)
   cd <- cluster_data()
   n_seqs   <- nrow(cd)
   chain_motif <- input$chain_select_motif
   chain_prefix <- substr(chain_motif, 1, 1)
   cdr3  <- paste0("cdr3_", chain_prefix, "_aa")
   vgene <- paste0("v_", chain_prefix, "_gene")
   jgene <- paste0("j_", chain_prefix, "_gene")
   n_unique <- length(unique(cd[[cdr3]]))
   total_ct <- if ("count" %in% names(cd)) sum(cd$count, na.rm = TRUE) else NA
   dom_v <- if (vgene %in% names(cd) && nrow(cd) > 0) names(sort(table(cd[[vgene]]), decreasing = TRUE))[1] else "N/A"
   dom_j <- if (jgene %in% names(cd) && nrow(cd) > 0) names(sort(table(cd[[jgene]]), decreasing = TRUE))[1] else "N/A"
   
   tagList(
     h5("Cluster summary"),
     p(strong("Sequences in cluster: "), n_seqs),
     p(strong("Unique CDR3s: "), n_unique),
     if (!is.na(total_ct)) p(strong("Total clone count: "), total_ct),
     p(strong("Dominant V gene: "), dom_v),
     p(strong("Dominant J gene: "), dom_j)
   )
 })
 
 # ---- Cluster Details sub-tab -----------------------------------------
 output$cluster_seq_table <- renderDT({
   req(cluster_data())
   datatable(cluster_data(), options = list(pageLength = 10), rownames = FALSE)
 })
 
 output$v_gene_plot <- renderPlot({
   req(cluster_data())
   cd <- cluster_data()
   chain_motif <- input$chain_select_motif
   vgene<-paste0("v_",substr(chain_motif,1,1),"_gene")
   barplot(
     table(cd[[vgene]]),
     las = 2, col = "steelblue",
     main = NULL, ylab = "Count"
   )
 })
 
 output$j_gene_plot <- renderPlot({
   req(cluster_data())
   cd <- cluster_data()
   chain_motif <- input$chain_select_motif
   jgene<-paste0("j_",substr(chain_motif,1,1),"_gene")
   barplot(
     table(cd[[jgene]]),
     las = 2, col = "darkorange",
     main = NULL, ylab = "Count"
   )
 })
 
 ##############################################################################
 #                     Explore Results Tab                                    #
 ##############################################################################
 
 # Populate Explorer sidebar UI controls once tcrdist3 analysis is complete
 observe({
   tr <- analysis_results()
   req(tr)
   clone_df <- as.data.frame(tr$clone_df)
   epitopes <- unique(as.character(clone_df$epitope))
   subjects <- unique(as.character(clone_df$subject))
   
   all_names <- names(tr)
   pw_names <- all_names[grepl("^pw_", all_names)]
   chains_avail <- gsub("^pw_", "", pw_names[!grepl("cdr|pmhc", pw_names)])
   
   # Include a combined distance option when more than 1 chain was analyzed
   pw_choices <- setNames(pw_names, pw_names)
   if (length(chains_avail) > 1) {
     pw_choices <- c("Combined (Sum of Chains)" = "pw_combined", pw_choices)
   }
   
   color_candidates <- c(
     "Epitope" = "epitope",
     "Subject" = "subject",
     "V-Alpha Gene" = "v_a_gene",
     "J-Alpha Gene" = "j_a_gene",
     "V-Beta Gene" = "v_b_gene",
     "J-Beta Gene" = "j_b_gene",
     "V-Gamma Gene" = "v_g_gene",
     "J-Gamma Gene" = "j_g_gene",
     "V-Delta Gene" = "v_d_gene",
     "J-Delta Gene" = "j_d_gene",
     "CDR3 Length" = "cdr3_length"
   )
   avail_colors <- color_candidates[
     color_candidates %in% c(colnames(clone_df), "cdr3_length")
   ]
   
   output$explore_epitope_ui <- renderUI({
     selectizeInput(
       "explore_epitopes",
       "Filter Epitope(s):",
       choices = epitopes,
       selected = epitopes,
       multiple = TRUE,
       options = list(placeholder = "Select Epitope(s)")
     )
   })
   
   output$explore_subject_ui <- renderUI({
     selectizeInput(
       "explore_subjects",
       "Filter Subject(s):",
       choices = subjects,
       selected = subjects,
       multiple = TRUE,
       options = list(placeholder = "Select Subject(s)")
     )
   })
   
   output$explore_chain_ui <- renderUI({
     selectInput(
       "explore_chain",
       "Active Chain (for CDR3 / V-J):",
       choices = chains_avail,
       selected = chains_avail[1]
     )
   })
   
   output$explore_pw_matrix_ui <- renderUI({
     selectInput(
       "explore_pw_matrix",
       "Distance Matrix:",
       choices = pw_choices,
       selected = pw_choices[1]
     )
   })
   
   output$explore_color_by_ui <- renderUI({
     selectInput(
       "explore_color_by",
       "Color / Group By:",
       choices = avail_colors,
       selected = "epitope"
     )
   })
 })
 
 # Filtered clone_df and distance matrices for the Explorer tab
 explore_data <- reactive({
   tr <- analysis_results()
   req(tr, input$explore_epitopes, input$explore_subjects, input$explore_chain, input$explore_pw_matrix)
   
   clone_df <- as.data.frame(tr$clone_df)
   keep_idx <- which(
     clone_df$epitope %in% input$explore_epitopes &
       clone_df$subject %in% input$explore_subjects
   )
   req(length(keep_idx) >= 2)
   
   sub_df <- clone_df[keep_idx, , drop = FALSE]
   rownames(sub_df) <- NULL
   
   # Strip allele suffix (*01) on gene columns if requested
   gene_cols <- grep("_gene$", colnames(sub_df), value = TRUE)
   if (isTRUE(input$explore_strip_allele)) {
     for (gc in gene_cols) {
       sub_df[[gc]] <- gsub("\\*.*", "", as.character(sub_df[[gc]]))
     }
   }
   
   chain_prefix <- substr(input$explore_chain, 1, 1)
   cdr3_col <- paste0("cdr3_", chain_prefix, "_aa")
   v_col <- paste0("v_", chain_prefix, "_gene")
   j_col <- paste0("j_", chain_prefix, "_gene")
   
   if (cdr3_col %in% colnames(sub_df)) {
     sub_df$cdr3_active <- as.character(sub_df[[cdr3_col]])
     sub_df$cdr3_length <- nchar(sub_df$cdr3_active)
   } else {
     sub_df$cdr3_active <- NA_character_
     sub_df$cdr3_length <- NA_integer_
   }
   
   # Compute selected pairwise distance matrix
   pw_names <- names(tr)[grepl("^pw_", names(tr))]
   chains_avail <- gsub("^pw_", "", pw_names[!grepl("cdr|pmhc", pw_names)])
   
   if (input$explore_pw_matrix == "pw_combined") {
     mats <- lapply(paste0("pw_", chains_avail), function(nm) tr[[nm]][keep_idx, keep_idx, drop = FALSE])
     dist_mat <- Reduce(`+`, mats)
   } else {
     req(input$explore_pw_matrix %in% names(tr))
     dist_mat <- tr[[input$explore_pw_matrix]][keep_idx, keep_idx, drop = FALSE]
   }
   
   # Compute nearest-neighbor distance (excluding self) for each clone
   diag_saved <- diag(dist_mat)
   diag(dist_mat) <- Inf
   nn_idx <- apply(dist_mat, 1, which.min)
   nn_dist <- apply(dist_mat, 1, min)
   diag(dist_mat) <- diag_saved
   
   sub_df$nearest_neighbor_dist <- nn_dist
   sub_df$nearest_neighbor_epitope <- sub_df$epitope[nn_idx]
   if ("clone_id" %in% colnames(sub_df)) {
     sub_df$nearest_neighbor_clone <- sub_df$clone_id[nn_idx]
   } else {
     sub_df$nearest_neighbor_clone <- nn_idx
   }
   
   # Publicness: number of distinct subjects sharing the active CDR3 sequence
   if (cdr3_col %in% colnames(sub_df)) {
     pub_counts <- sub_df %>%
       dplyr::group_by(cdr3_active) %>%
       dplyr::summarise(
         n_subjects_sharing_cdr3 = dplyr::n_distinct(subject),
         n_epitopes_sharing_cdr3 = dplyr::n_distinct(epitope),
         .groups = "drop"
       )
     sub_df <- sub_df %>%
       dplyr::left_join(pub_counts, by = "cdr3_active") %>%
       dplyr::mutate(
         clonotype_sharing = ifelse(n_subjects_sharing_cdr3 > 1, "Public (>1 Subject)", "Private (1 Subject)")
       )
   } else {
     sub_df$n_subjects_sharing_cdr3 <- 1L
     sub_df$n_epitopes_sharing_cdr3 <- 1L
     sub_df$clonotype_sharing <- "Private (1 Subject)"
   }
   
   list(
     clone_df = sub_df,
     dist_mat = dist_mat,
     keep_idx = keep_idx,
     cdr3_col = cdr3_col,
     v_col = v_col,
     j_col = j_col,
     chains_avail = chains_avail
   )
 })
 
 # 1A. Distance Distribution Plot
 output$explore_dist_plot <- renderPlotly({
   ed <- explore_data()
   req(ed)
   dist_mat <- ed$dist_mat
   sub_df <- ed$clone_df
   
   upper_idx <- which(row(dist_mat) < col(dist_mat), arr.ind = TRUE)
   req(nrow(upper_idx) > 0)
   
   r_idx <- upper_idx[, 1]
   c_idx <- upper_idx[, 2]
   
   pairs_df <- data.frame(
     dist = dist_mat[upper_idx],
     epi_1 = as.character(sub_df$epitope[r_idx]),
     epi_2 = as.character(sub_df$epitope[c_idx]),
     subj_1 = as.character(sub_df$subject[r_idx]),
     subj_2 = as.character(sub_df$subject[c_idx]),
     stringsAsFactors = FALSE
   ) %>%
     dplyr::mutate(
       epitope_pair_type = ifelse(epi_1 == epi_2, "Within Same Epitope", "Between Different Epitopes"),
       subject_pair_type = ifelse(subj_1 == subj_2, "Within Same Subject", "Between Different Subjects"),
       within_epitope = ifelse(epi_1 == epi_2, epi_1, NA_character_)
     )
   
   group_mode <- input$explore_dist_group_mode
   if (group_mode == "within_epitope") {
     pairs_df <- pairs_df %>% dplyr::filter(!is.na(within_epitope))
     req(nrow(pairs_df) > 0)
     pairs_df$group <- pairs_df$within_epitope
   } else {
     pairs_df$group <- pairs_df[[group_mode]]
   }
   
   # Cap pair sample size for smooth interactive plotting on very large repertoires
   if (nrow(pairs_df) > 25000) {
     set.seed(42)
     pairs_df <- pairs_df[sample(nrow(pairs_df), 25000), , drop = FALSE]
   }
   
   plot_type <- input$explore_dist_plot_type
   if (plot_type == "violin") {
     p <- ggplot(pairs_df, aes(x = group, y = dist, fill = group)) +
       geom_violin(alpha = 0.6, trim = FALSE) +
       geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
       labs(
         title = paste("Pairwise TCR Distance Distribution (", input$explore_pw_matrix, ")"),
         x = NULL,
         y = "Pairwise TCR Distance",
         fill = "Pair Group"
       ) +
       theme_minimal(base_size = 12)
   } else if (plot_type == "density") {
     p <- ggplot(pairs_df, aes(x = dist, color = group, fill = group)) +
       geom_density(alpha = 0.3, linewidth = 1) +
       labs(
         title = paste("Pairwise TCR Distance Density (", input$explore_pw_matrix, ")"),
         x = "Pairwise TCR Distance",
         y = "Density",
         color = "Pair Group",
         fill = "Pair Group"
       ) +
       theme_minimal(base_size = 12)
   } else if (plot_type == "ecdf") {
     p <- ggplot(pairs_df, aes(x = dist, color = group)) +
       stat_ecdf(geom = "step", linewidth = 1.1) +
       labs(
         title = paste("Empirical Cumulative Distribution of TCR Distances (", input$explore_pw_matrix, ")"),
         x = "Pairwise TCR Distance",
         y = "Cumulative Probability F(d)",
         color = "Pair Group"
       ) +
       theme_minimal(base_size = 12)
   } else {
     p <- ggplot(pairs_df, aes(x = dist, fill = group)) +
       geom_histogram(bins = 35, alpha = 0.65, position = "identity") +
       labs(
         title = paste("Pairwise TCR Distance Histogram (", input$explore_pw_matrix, ")"),
         x = "Pairwise TCR Distance",
         y = "Pair Count",
         fill = "Pair Group"
       ) +
       theme_minimal(base_size = 12)
   }
   
   ggplotly(p) %>% config(responsive = TRUE)
 })
 
 # 1B. Chain vs. Chain (or CDR3 vs. full chain) Distance Comparison Plot
 output$explore_chain_corr_plot <- renderPlotly({
   ed <- explore_data()
   tr <- analysis_results()
   req(ed, tr)
   keep_idx <- ed$keep_idx
   sub_df <- ed$clone_df
   chains_avail <- ed$chains_avail
   
   upper_idx <- which(row(ed$dist_mat) < col(ed$dist_mat), arr.ind = TRUE)
   req(nrow(upper_idx) > 0)
   if (nrow(upper_idx) > 8000) {
     set.seed(42)
     upper_idx <- upper_idx[sample(nrow(upper_idx), 8000), , drop = FALSE]
   }
   
   if (length(chains_avail) >= 2) {
     ch1 <- chains_avail[1]
     ch2 <- chains_avail[2]
     mat1 <- tr[[paste0("pw_", ch1)]][keep_idx, keep_idx, drop = FALSE]
     mat2 <- tr[[paste0("pw_", ch2)]][keep_idx, keep_idx, drop = FALSE]
     x_lab <- paste0("pw_", ch1, " Distance")
     y_lab <- paste0("pw_", ch2, " Distance")
   } else {
     ch1 <- chains_avail[1]
     ch_prefix <- substr(ch1, 1, 1)
     cdr3_pw <- paste0("pw_cdr3_", ch_prefix, "_aa")
     mat1 <- tr[[paste0("pw_", ch1)]][keep_idx, keep_idx, drop = FALSE]
     mat2 <- if (cdr3_pw %in% names(tr)) tr[[cdr3_pw]][keep_idx, keep_idx, drop = FALSE] else mat1
     x_lab <- paste0("pw_", ch1, " (Full Chain Distance)")
     y_lab <- paste0(cdr3_pw, " (CDR3 Distance)")
   }
   
   r_idx <- upper_idx[, 1]
   c_idx <- upper_idx[, 2]
   df_corr <- data.frame(
     d1 = mat1[upper_idx],
     d2 = mat2[upper_idx],
     pair_type = ifelse(
       sub_df$epitope[r_idx] == sub_df$epitope[c_idx],
       paste0("Same Epitope (", sub_df$epitope[r_idx], ")"),
       "Different Epitopes"
     ),
     stringsAsFactors = FALSE
   )
   
   p <- ggplot(df_corr, aes(x = d1, y = d2, color = pair_type)) +
     geom_point(alpha = 0.55, size = 1.8) +
     labs(
       title = paste(x_lab, "vs.", y_lab),
       x = x_lab,
       y = y_lab,
       color = "Pair Type"
     ) +
     theme_minimal(base_size = 12)
   
   ggplotly(p) %>% config(responsive = TRUE)
 })
 
 # 2. 2D / 3D Classical MDS Embedding Plot
 output$explore_mds_plot <- renderPlotly({
   ed <- explore_data()
   req(ed, input$explore_color_by)
   dist_mat <- ed$dist_mat
   sub_df <- ed$clone_df
   
   k_dims <- ifelse(input$explore_mds_dims == "3d", 3, 2)
   req(nrow(dist_mat) > k_dims)
   
   mds_coords <- tryCatch({
     stats::cmdscale(as.dist(dist_mat), k = k_dims, eig = FALSE)
   }, error = function(e) {
     matrix(0, nrow = nrow(dist_mat), ncol = k_dims)
   })
   
   if (ncol(mds_coords) < k_dims) {
     mds_coords <- cbind(mds_coords, matrix(0, nrow = nrow(dist_mat), ncol = k_dims - ncol(mds_coords)))
   }
   
   color_col <- input$explore_color_by
   if (!color_col %in% colnames(sub_df)) {
     color_col <- "epitope"
   }
   
   sub_df$MDS1 <- mds_coords[, 1]
   sub_df$MDS2 <- mds_coords[, 2]
   sub_df$color_group <- sub_df[[color_col]]
   sub_df$pt_size <- if (isTRUE(input$explore_size_by_count) && "count" %in% colnames(sub_df)) {
     pmax(6, pmin(22, 6 + 3 * log2(sub_df$count + 1)))
   } else {
     rep(9, nrow(sub_df))
   }
   
   hover_txt <- paste0(
     "<b>Epitope:</b> ", sub_df$epitope,
     "<br><b>Subject:</b> ", sub_df$subject,
     "<br><b>Count:</b> ", sub_df$count,
     "<br><b>Active CDR3:</b> ", sub_df$cdr3_active,
     "<br><b>Nearest Neighbor Dist:</b> ", round(sub_df$nearest_neighbor_dist, 1),
     "<br><b>", color_col, ":</b> ", sub_df$color_group
   )
   
   if (input$explore_mds_dims == "3d") {
     sub_df$MDS3 <- mds_coords[, 3]
     plot_ly(
       data = sub_df,
       x = ~MDS1,
       y = ~MDS2,
       z = ~MDS3,
       color = ~as.factor(color_group),
       text = hover_txt,
       hoverinfo = "text",
       type = "scatter3d",
       mode = "markers",
       marker = list(size = ~pt_size / 1.5, opacity = 0.8)
     ) %>%
       layout(
         title = paste("3D MDS of TCR Distances (", input$explore_pw_matrix, ") - Colored by", color_col),
         scene = list(
           xaxis = list(title = "MDS 1"),
           yaxis = list(title = "MDS 2"),
           zaxis = list(title = "MDS 3")
         )
       )
   } else {
     plot_ly(
       data = sub_df,
       x = ~MDS1,
       y = ~MDS2,
       color = ~as.factor(color_group),
       text = hover_txt,
       hoverinfo = "text",
       type = "scatter",
       mode = "markers",
       marker = list(size = ~pt_size, opacity = 0.78)
     ) %>%
       layout(
         title = paste("2D MDS of TCR Distances (", input$explore_pw_matrix, ") - Colored by", color_col),
         xaxis = list(title = "MDS 1"),
         yaxis = list(title = "MDS 2")
       )
   }
 })
 
 # 3A. V-J Gene Pairing Heatmap
 output$explore_vj_heatmap <- renderPlotly({
   ed <- explore_data()
   req(ed)
   sub_df <- ed$clone_df
   v_col <- ed$v_col
   j_col <- ed$j_col
   req(v_col %in% colnames(sub_df), j_col %in% colnames(sub_df))
   
   top_n <- max(5, as.integer(input$explore_top_n_genes %||% 15))
   metric <- input$explore_vj_metric %||% "clones"
   
   top_v <- sub_df %>%
     dplyr::count(.data[[v_col]], wt = if (metric == "count") count else NULL, sort = TRUE) %>%
     dplyr::slice_head(n = top_n) %>%
     dplyr::pull(1)
   top_j <- sub_df %>%
     dplyr::count(.data[[j_col]], wt = if (metric == "count") count else NULL, sort = TRUE) %>%
     dplyr::slice_head(n = top_n) %>%
     dplyr::pull(1)
   
   vj_df <- sub_df %>%
     dplyr::filter(.data[[v_col]] %in% top_v, .data[[j_col]] %in% top_j) %>%
     dplyr::group_by(V_Gene = .data[[v_col]], J_Gene = .data[[j_col]]) %>%
     dplyr::summarise(
       Abundance = if (metric == "count") sum(count, na.rm = TRUE) else dplyr::n(),
       .groups = "drop"
     )
   
   p <- ggplot(vj_df, aes(x = J_Gene, y = V_Gene, fill = Abundance)) +
     geom_tile(color = "white") +
     scale_fill_viridis_c(option = "plasma") +
     labs(
       title = paste("V-J Gene Pairing (", toupper(input$explore_chain), "Chain)"),
       x = "J Gene",
       y = "V Gene",
       fill = ifelse(metric == "count", "Total Count", "Clones")
     ) +
     theme_minimal(base_size = 11) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   ggplotly(p) %>% config(responsive = TRUE)
 })
 
 # 3B. Top V & J Gene Barplot
 output$explore_vj_barplot <- renderPlotly({
   ed <- explore_data()
   req(ed)
   sub_df <- ed$clone_df
   v_col <- ed$v_col
   j_col <- ed$j_col
   req(v_col %in% colnames(sub_df), j_col %in% colnames(sub_df))
   
   top_n <- max(5, as.integer(input$explore_top_n_genes %||% 15))
   metric <- input$explore_vj_metric %||% "clones"
   
   v_long <- sub_df %>%
     dplyr::transmute(epitope = epitope, segment = "V Gene", gene = .data[[v_col]], count = count)
   j_long <- sub_df %>%
     dplyr::transmute(epitope = epitope, segment = "J Gene", gene = .data[[j_col]], count = count)
   
   combined_genes <- dplyr::bind_rows(v_long, j_long)
   top_genes <- combined_genes %>%
     dplyr::group_by(segment, gene) %>%
     dplyr::summarise(total = if (metric == "count") sum(count, na.rm = TRUE) else dplyr::n(), .groups = "drop") %>%
     dplyr::group_by(segment) %>%
     dplyr::slice_max(total, n = top_n, with_ties = FALSE)
   
   plot_df <- combined_genes %>%
     dplyr::inner_join(top_genes %>% dplyr::select(segment, gene), by = c("segment", "gene")) %>%
     dplyr::group_by(segment, gene, epitope) %>%
     dplyr::summarise(
       Abundance = if (metric == "count") sum(count, na.rm = TRUE) else dplyr::n(),
       .groups = "drop"
     )
   
   p <- ggplot(plot_df, aes(x = reorder(gene, Abundance), y = Abundance, fill = epitope)) +
     geom_col() +
     coord_flip() +
     facet_wrap(~segment, scales = "free") +
     labs(
       title = paste("Top V & J Gene Usage by Epitope (", toupper(input$explore_chain), "Chain)"),
       x = NULL,
       y = ifelse(metric == "count", "Total Template Count", "Number of Clones"),
       fill = "Epitope"
     ) +
     theme_minimal(base_size = 11)
   
   ggplotly(p) %>% config(responsive = TRUE)
 })
 
 # 4A. CDR3 Amino Acid Length Distribution (Spectratype)
 output$explore_cdr3_len_plot <- renderPlotly({
   ed <- explore_data()
   req(ed)
   sub_df <- ed$clone_df %>% dplyr::filter(!is.na(cdr3_length))
   req(nrow(sub_df) > 0)
   
   len_df <- sub_df %>%
     dplyr::group_by(epitope, cdr3_length) %>%
     dplyr::summarise(n_clones = dplyr::n(), total_count = sum(count, na.rm = TRUE), .groups = "drop")
   
   p <- ggplot(len_df, aes(x = factor(cdr3_length), y = n_clones, fill = epitope)) +
     geom_col(position = "dodge") +
     labs(
       title = paste("CDR3 Length Spectratype (", toupper(input$explore_chain), "Chain)"),
       x = "CDR3 Amino Acid Length (aa)",
       y = "Number of Clones",
       fill = "Epitope"
     ) +
     theme_minimal(base_size = 12)
   
   ggplotly(p) %>% config(responsive = TRUE)
 })
 
 # 4B. Dynamic selector for CDR3 length logo
 output$explore_logo_length_ui <- renderUI({
   ed <- explore_data()
   req(ed)
   sub_df <- ed$clone_df %>% dplyr::filter(!is.na(cdr3_length), cdr3_length > 0)
   req(nrow(sub_df) > 0)
   
   len_tab <- sort(table(sub_df$cdr3_length), decreasing = TRUE)
   avail_lens <- as.integer(names(len_tab))
   choices_labeled <- setNames(avail_lens, paste0(avail_lens, " aa (n=", as.integer(len_tab), ")"))
   
   selectInput(
     "explore_logo_len",
     "Select CDR3 Length for Sequence Logo:",
     choices = choices_labeled,
     selected = avail_lens[1]
   )
 })
 
 # 4C. Length-Stratified Sequence Logo
 output$explore_cdr3_logo <- renderPlot({
   ed <- explore_data()
   req(ed, input$explore_logo_len)
   target_len <- as.integer(input$explore_logo_len)
   seqs <- ed$clone_df %>%
     dplyr::filter(cdr3_length == target_len, !is.na(cdr3_active)) %>%
     dplyr::pull(cdr3_active)
   req(length(seqs) >= 1)
   
   ggseqlogo(seqs, method = input$explore_logo_method %||% "bits", seq_type = "aa") +
     labs(
       title = paste0(
         toupper(input$explore_chain),
         " Chain CDR3 Sequence Logo (Length = ",
         target_len,
         " aa, n = ",
         length(seqs),
         " clones)"
       )
     ) +
     theme_minimal(base_size = 13)
 })
 
 # 5A. Nearest-Neighbor Distance vs. Clone Count Plot
 output$explore_nn_plot <- renderPlotly({
   ed <- explore_data()
   req(ed)
   sub_df <- ed$clone_df
   
   hover_txt <- paste0(
     "<b>Active CDR3:</b> ", sub_df$cdr3_active,
     "<br><b>Epitope:</b> ", sub_df$epitope,
     "<br><b>Subject:</b> ", sub_df$subject,
     "<br><b>Clone Count:</b> ", sub_df$count,
     "<br><b>Nearest Neighbor Dist:</b> ", round(sub_df$nearest_neighbor_dist, 1),
     "<br><b>Nearest Neighbor Epitope:</b> ", sub_df$nearest_neighbor_epitope,
     "<br><b>Subjects Sharing CDR3:</b> ", sub_df$n_subjects_sharing_cdr3
   )
   
   p <- ggplot(
     sub_df,
     aes(
       x = nearest_neighbor_dist,
       y = count,
       color = epitope,
       shape = clonotype_sharing,
       text = hover_txt
     )
   ) +
     geom_point(alpha = 0.8, size = 3) +
     scale_y_log10() +
     labs(
       title = paste("Nearest-Neighbor TCR Distance vs. Clone Count (", input$explore_pw_matrix, ")"),
       x = "Distance to Closest Neighbor TCR",
       y = "Clone Count (log10 scale)",
       color = "Epitope",
       shape = "CDR3 Sharing"
     ) +
     theme_minimal(base_size = 12)
   
   ggplotly(p, tooltip = "text") %>% config(responsive = TRUE)
 })
 
 # 5B. Enriched Clonotype Table
 output$explore_clones_table <- DT::renderDataTable({
   ed <- explore_data()
   req(ed)
   sub_df <- ed$clone_df
   
   DT::datatable(
     sub_df,
     caption = "Filtered & Enriched TCRdist3 Clonotype Explorer Table",
     extensions = "Buttons",
     rownames = FALSE,
     filter = "top",
     options = list(
       scrollX = TRUE,
       pageLength = 10,
       dom = "lfrtipB",
       buttons = list(
         "copy",
         list(extend = "csv", text = "CSV", filename = paste0("TCRdist3-Explorer-", session$token)),
         list(extend = "excel", text = "Excel", filename = paste0("TCRdist3-Explorer-", session$token))
       )
     )
   )
 }, server = TRUE)
 
}
