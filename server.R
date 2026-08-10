#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
source("global.R")
# A reusable function to render a styled error message
# --- Server Logic ---
server <- function(input, output, session) {
  track_usage(
    storage_mode = store_json(path = "logs/")
  )
  neighbour_results <- NULL
  # Per-session flag (was a module-level global shared across all users)
  hierarch_tabs_check <- reactiveVal(FALSE)
  conda_path <- Sys.which("python")
  render_error_ui("", output = output)
  render_error_ui(conda_path, output = output)
  
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
    read.csv(input$file1$datapath, sep = sep)
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
    observe({
      tr<-analysis_results()
      req(tr)
      chains<-strsplit(input$chains,",")[[1]]
      if(length(chains)==2)
        pw_alpha_beta_max <- max(tr$pw_beta + tr$pw_alpha)
      else
        pw_alpha_beta_max<-max(tr[[paste0("pw_",chains)]])
      net_data<-generate_network(tr,120)
      output$networkcard<-renderUI({
        tagList(
          card(
            layout_sidebar(
              sidebar = sidebar(
                helpText("Filter edges by distance and choose how to color the nodes."),
                # Slider to filter maximum distance
                sliderInput("max_dist", 
                            "Maximum TCR Distance (Edge Threshold):", 
                            min = 0, 
                            max = pw_alpha_beta_max,
                            value = 120, 
                            step = 1),
                # Dropdown to color nodes
                selectInput("color_var", 
                            "Color Nodes By:", 
                            choices = c("Epitope" = "epitope", 
                                        "Subject" = "subject", 
                                        "V Beta Gene" = "v_b_gene",
                                        "V Alpha Gene" = "v_a_gene",
                                        "CDR3 Alpha AA" = "cdr3_a_aa",
                                        "CDR3 Beta AA" ="cdr3_b_aa")),
                hr()
              ),
              card(
                p(strong("Network Stats:")),
                textOutput("node_count"),
                textOutput("edge_count"),
                hr(),
                visNetworkOutput("network_plot", height = "750px")
              ))
          ))
      })
    })
    
    # Render the network visualization
    observeEvent(c(input$max_dist, input$color_var), {
      req(analysis_results())
      tryCatch({
        tr <- analysis_results()
        chains <- strsplit(input$chains, ",")[[1]]
        
        # Recompute network with current edge threshold
        if (length(chains) == 2) {
          pw_alpha_beta_max <- max(tr$pw_beta + tr$pw_alpha)
        } else {
          pw_alpha_beta_max <- max(tr[[paste0("pw_", chains)]])
        }
        
        net_data <- generate_network(tr, input$max_dist)
        
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
        
        # Create nodes data frame with color variable
        color_var <- input$color_var
        nodes <- data.frame(
          id = all_node_ids,
          label = all_node_ids,
          title = paste0(
            "Clone ID: ", all_node_ids, "<br>",
            "Epitope: ", clone_df[[color_var]][all_node_ids + 1],"<br>",
            "V-Alpha: ", clone_df$v_a_gene[all_node_ids + 1],"<br>",
            "V-Beta: ", clone_df$v_b_gene[all_node_ids + 1],"<br>",
            "CDR3-A: ", clone_df$cdr3_a_aa[all_node_ids + 1],"<br>",
            "CDR3-B: ", clone_df$cdr3_b_aa[all_node_ids + 1]),
          group = as.character(clone_df[[color_var]][all_node_ids + 1]),
          size = 20
        )
        
        # Create edges data frame with weights inversely proportional to distance
        # (closer clones = thicker edges)
        edges <- data.frame(
          from = net_data$node_1,
          to = net_data$node_2,
          width = 5 - (net_data$dist / max(net_data$dist)) * 4,  # Width between 1-5
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
    output$chain_select_ui_trees <- renderUI({
      tagList(
        selectInput(
          "chain_select_tree",
          "Select a chain:",
          selected = chains_tr_only[1],
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
  
  
  #######################Plot Diversity plots###################################
  output$diversityPlot <- renderPlot({
    req(analysis_results())
    counts <- analysis_results()$cell_df$count
    # Generate a simple Lorenz Curve for Diversity
    df_plot <- data.frame(counts = sort(counts, decreasing = TRUE)) %>%
      mutate(
        cum_counts = cumsum(counts) / sum(counts),
        index = row_number() / n()
      )
    
    ggplot(df_plot, aes(x = index, y = cum_counts)) +
      geom_line(color = "steelblue", size = 1.2) +
      geom_abline(lty = 2) +
      labs(title = "TCR Repertoire Evenness (Lorenz Curve)",
           x = "Fraction of Unique Clotypes",
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
      return(tr)
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
        message = paste0("Run KNN with ", input$epitope_sel_n, ":"),
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
          # diff testing is pasted on binary comparison, so all epitope not 'PA' are set to 'X'
          epitope_n <-epi_1<-input$epitope_sel_n
          if(length(epi_1)>1)
            epi_1<-epi_1[1]
          
          message("Selected epitope is ",epi_1,"\n")
          
          if(length(epitope_n)==1)
          {
            # Code to handle mu]ltiple epitopes
            tr$clone_df[,epi_1]= ifelse(
              tr$clone_df$epitope==epi_1,epi_1,"Other")
            epitope_n <- list(epitope_n)
          }
          else{
            tr$clone_df[, epi_1] = ifelse(tr$clone_df$epitope ==
                                                 epi_1, epi_1,epitope_n[2])
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
            clone_df = tr$clone_df,
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
          ) %>% arrange('FDRq') %>% arrange(desc(OR), desc(ct_0))
          output$fixedradius_output_ui <- renderUI({
            withSpinner(DT::dataTableOutput("fixedradius_output"))
          })
          incProgress(0.9, "Displaying results table...")
          
          session.id <- session$token
          file.name = paste0(
            'TCRdist3-KNN-clustering-',
            input$epitope_sel_n,
            "-",
            chain,
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
  
  
  observeEvent(input$run_Hierarch, {
    hierarch_tabs_check(FALSE)
    #neighbour_results()
  })
  
  epitope_results_hierarch <- eventReactive(input$epitope_sel_h, {
    render_error_ui("", output = output) #Clear any error messages
    req(tcr_data())
    epitope_h_sel<-input$epitope_sel_h
    tr<-analysis_results()
    if(length(epitope_h_sel)==1)
    {
      return(tr)
    }
    else{
      return(subset_tr_by_epitope(analysis_results(),epitope_h_sel))
    }
  })
  
  #Listen to hierarchical neighbourhood clustering events
  #neighbour_results <- eventReactive(
  observeEvent(input$run_Hierarch, {
    render_error_ui("", output = output) #Clear any error messages
    req(!hierarch_tabs_check())
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
        
        if(length(epitope_h)==1)
        {
          # Code to handle mu]ltiple epitopes
          tr.this$clone_df[,epi_1]= ifelse(
            tr.this$clone_df$epitope==epi_1,epi_1,"Other")
          epitope_h <- list(epitope_h)
        }
        else{
          tr.this$clone_df[, epi_1] = ifelse(tr.this$clone_df$epitope ==
                                               epi_1, epi_1,epitope_h[2])
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
        
        res_Z <- rep_diff$hcluster_diff(tr.this$clone_df,pwmat,x_cols = r_to_py(as.list(c(epi_1))),count_col = 'count')
        res <- res_Z[[1]]
        Z <- res_Z[[2]]
        
        res_summary <- rep_diff$member_summ(
          res_df = res,
          clone_df = tr.this$clone_df,
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
        file.copy(temp_path, file.path("www/", html.file))
        
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
  
  
  #Generate Network Graph that will
  observeEvent(input$chain_select_neighbours, {
    tr <- analysis_results()
    
    tryCatch({
      req(tr)
      chain <- input$chain_select_neighbours
      genes <- tr$clone_df[[paste0("v_", substr(chain, 1, 1), "_gene")]]
      chain2check <- paste0("cdr3_", substring(chain, 1, 1), "_aa")
      pw_matrix = tr[[paste0("pw_", chain2check)]]
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
            "Max Distance to Show:",
            min = floor(min(pw_matrix)),
            max = ceiling(max(pw_matrix)),
            value = ceiling(median(pw_matrix))
          )
        )
      })
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
  #Populate trees
  observeEvent(input$run_trees, {
    render_error_ui("", output = output) #Clear any error messages
    req(tcr_data())
    withProgress(message = 'Running trees...', value = 0, {
      tryCatch({
        #If more than one chain has been selected, then recompute distance matrix based on chain
        chains <- strsplit(input$chains, ",")[[1]]
        output$showtrees <- reactive({
          return("show")
        })
        outputOptions(output, "showtrees", suspendWhenHidden = FALSE)
        output$tree_output_plot<-renderUI(NULL)
        if(length(chains)>1)
        {
          pd <- import("pandas")
          inspect <- import("inspect")
          message("Finished importing pandas")
          TCRrep <- import("tcrdist.repertoire")$TCRrep
          
          message("Finished importing TCRrep")
          incProgress(0.3, detail = "Preparing data...")
          
          # Convert R data frame to pandas DataFrame
          cell_df <- r_to_py(tcr_data())
          # Get chains for analysis
          py <- import_main()
          
          message(paste0(chains, collapse = ":"))
          chains <- list(input$chain_select_tree)
          message(chains)
          incProgress(0.5, detail = "Computing distances...")
          # Initialize TCRrep object and compute distances
          tr <- TCRrep(
            cell_df = cell_df,
            organism = r_to_py(input$organism),
            chains = r_to_py(chains),
            db_file = r_to_py(input$dbfile)
          )  
        } else{
          tr <- analysis_results()
        }
        req(tr)
        incProgress(0.1, detail = "Importing required library...")
        
        TCRtree <- import("tcrdist.tree")
        incProgress(0.3, detail = "Preparing data...")
        file.name <- paste0(session$token, "_tree.html")
        temp.file <- tempfile()
        message("Temp file for trees: ", temp.file)
        remove_www_tempfiles(session, "tree.html")
        message("Running trees")
        tcrtree = TCRtree$TCRtree(tcrrep = tr, html_name = temp.file)
        incProgress(0.5, detail = "Building tree...")
        tcrtree$build_tree()
        incProgress(0.9, detail = "Finalizing...")
        file.copy(temp.file, paste0("www/", file.name))
        output$tree_output_plot <- renderUI({
          withSpinner(tags$iframe(
            src = file.name,
            width = "95%",
            style = "height: 90vh;",
            # Use CSS for full view port height
            frameBorder = "0",
            onload = "iframeLoaded()"
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
    ggplot(data.frame(val = data$one_minus_simpson, label = "Repertoire"), 
           aes(x = label, y = val)) +
      geom_point(color = "red", size = 4) +
      ylim(0, 1) +
      labs(title = "1-Simpson's Diversity (Clonality)", y = "1 - Simpson") +
      theme_minimal()
  })
  
  # Bottom Plot: Inverse Simpson's (Diversity)
  output$invSimpsonPlot <- renderPlot({
    data <- processed_data()
    ggplot(data.frame(val = data$inv_simpson, label = "Repertoire"), 
           aes(x = label, y = val)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      scale_y_log10() +
      labs(title = "Inverse Simpson's Diversity", y = "Index (Log Scale)") +
      theme_minimal()
  })
  
}
