#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
source("global.R")
# A reusable function to render a styled error message

#Declare global variable
hierarch_tabs.check <- F

render_error_ui <- function(message, output) {
  # if(message!="")
  #   shinyalert(
  #     title = "Error!",
  #     text = message,
  #     type = "error"
  #   )
  output$errorreport <- renderUI(HTML(ifelse(
    message == "",
    "",
    paste0(
      '<div style="color: red; font-weight: bold; padding: 10px; border: 1px solid red; border-radius: 5px; background-color: #ffe6e6;">',
      message,
      '</div>'
    )
  )))
  #feedbackDanger("my_text_input", TRUE, messa)
  #showNotification(message)
  
}

remove_www_tempfiles <- function(session, pattern = "*")
{
  temp.files <- list.files("www",
                           pattern = paste0(session$token, ".", pattern),
                           full.names = T)
  message(paste0(temp.files, collapse = ";"))
  for (temp_file in temp.files)
  {
    if (file.exists(temp_file)) {
      file.remove(temp_file)
      message("Temporary file removed: ", temp_file)
    }
  }
}

# --- Server Logic ---
server <- function(input, output, session) {
  neighbour_results<-NULL
  conda_path <- Sys.which("python")
  render_error_ui("", output = output)
  render_error_ui(conda_path, output = output)
  
  output$csvtest <- reactive({
    !is.null(input$file1)
  })
  
  outputOptions(output, "csvtest", suspendWhenHidden = FALSE)
  # Reactive value to store the TCR data
  tcr_data <- reactive({
    render_error_ui("", output = output) #Clear any error messages
    req(input$file1)
    # Determine separator based on file extension
    sep <- ifelse(tools::file_ext(input$file1$name) == "tsv", "\t", ",")
    read.csv(input$file1$datapath, sep = sep)
  })
  
  
  # Run tcrdist3 analysis when the run button is clicked
  analysis_results <- eventReactive(input$run, {
    render_error_ui("", output = output) #Clear any error messages
    req(tcr_data())
    if(!("epitope" %in% colnames(tcr_data())) || !("subject" %in% colnames(tcr_data())))
    {
      shinyalert(
        title = "Missing columns",
        text = "Column names: 'epitope' or 'subject' missing in input data",
        type = "error"
      )
    }
    
    hierarch_tabs.check <<- F
    render_error_ui("", output = output)
    withProgress(message = 'Running TCRdist3 analysis...', value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Importing Python modules...")
        #message(py_config())
        #Import necessary python libraries
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
        
        chains <- strsplit(input$chains, ",")[[1]]
        message(paste0(chains, collapse = ":"))
        
        if (length(chains) == 1)
          chains <- list(chains)
        
        incProgress(0.5, detail = "Computing distances...")
        # Initialize TCRrep object and compute distances
        tr <- TCRrep(
          cell_df = cell_df,
          organism = r_to_py(input$organism),
          chains = r_to_py(chains),
          db_file = r_to_py(input$dbfile)
        )
        #message(paste0(names(tr),collapse = "\n"))
        incProgress(0.8, detail = "Finalizing...")
        #Get chains from the distance compute
        all.tr_names <- names(tr)
        chains_tr <- all.tr_names[grepl("^pw_", all.tr_names)]
        chains_tr_only <- gsub("pw_", "", chains_tr[!grepl("cdr|pmhc", chains_tr)])
        #Include a select option for the chains
        output$matrixoptions <- renderUI({
          tagList(
            hr(),
            selectInput(
              "matrix_select",
              "Select distance matrix:",
              selected = chains_tr[1],
              choices = chains_tr
            ),
            hr()
          )
        })
        
        epitopes <- unique(tr$clone_df$epitope)
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
        #Load epitopes in Neighbourhood
        output$epitope_select_ui <- renderUI({
          tagList(selectizeInput("epitope_sel_n", "Select Epitope:", choices = epitopes))
        })
        
        output$epitope_select_hierarch_ui <- renderUI({
          tagList(selectizeInput("epitope_sel_h", "Select Epitope:", choices = epitopes))
        })
        
        #Load chains in Neighbourhood
        output$chain_select_ui <- renderUI({
          tagList(
            selectInput(
              "chain_select_neighbours",
              "Select a chain:",
              selected = chains_tr_only[1],
              choices = chains_tr_only
            )
          )
        })
          output$chain_select_ui_knn <- renderUI({
            tagList(
              selectInput(
                "chain_select_knn",
                "Select a chain:",
                selected = chains_tr_only[1],
                choices = chains_tr_only
              )
            )
        })
        
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
  
  #Update heatmap and matrix table on the following events
  observeEvent(c(input$run, input$matrix_select), {
    render_error_ui("", output = output) #Clear any error messages
    tryCatch({
      results <- analysis_results()
      req(results)
      render_error_ui("", output)
      chain2check <- gsub("pw_","", (input$matrix_select))
      if (sum(chain2check %in% colnames(results$clone_df))==0)
        chain2check <- paste0("cdr3_", substring(chain2check,1, 1), "_aa")
      message(chain2check)
      #pw_matrix = results[[paste0("pw_", chain2check)]]
      pw_matrix <- results[[input$matrix_select]]
      clone_df = results$clone_df
      mat <- pw_matrix
      rownames(mat) <- clone_df[[chain2check]]
      colnames(mat) <- clone_df[[chain2check]]
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
              list(extend = 'csv', text = 'Download CSV',filename=paste0("TCRDist3-clone-",session$token)),
              list(extend = 'excel', text = 'Download Excel',filename=paste0("TCRDist3-clone-",session$token))
            )
          )
        )
      },server = F)
      
      # Display the distance matrix
      session.id<-session$token
      output$dist_matrix <- DT::renderDataTable({
        DT::datatable(
          mat,
          extensions = 'Buttons',
          options = list(
            scrollX = TRUE,
            pageLength = 10,
            rownames = F,
            dom = 'lfrtipB', #'flrtBip',
            buttons  = list('copy', list(extend='csv',text="csv",filename=paste0('TCR-distance-matrix-',input$matrix_select,"-",session.id)), 
                         list(extend="excel",text='excel',filename=paste0('TCR-distance-matrix-',input$matrix_select,"-",session.id)), 
                         list(extend='pdf',text='pdf',filename=paste0('TCR-distance-matrix-',input$matrix_select,"-",session.id)), 
                         list(extend='print',text="print",filename=paste0('TCR-distance-matrix-',input$matrix_select,"-",session.id)))
          )
        )
      },server = F)
      
      # Display the heatmap
      output$heatmap <- renderPlot({
        pheatmap(
          mat,
          clustering_distance_rows = "euclidean",
          clustering_distance_cols = "euclidean",
          clustering_method = "ward.D2",
          main = "TCR Distance Heatmap"
        )
      })
    }, error = function(e)
    {
      if (e$message != "")
      {
        message.error <- paste0("Error rendering heatmap: ", e$message)
        render_error_ui(message.error, output)
      }
      
    })
    
  })
  
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
      summary.stats <- do.call(rbind, lapply(unique(clone_df$epitope), function(x)
      {
        data.frame(
          epitope = x,
          num_individuals = length(unique(
            subset(clone_df, epitope == x)$subject
          )),
          num_clones = length(unique(
            subset(clone_df, epitope == x)$clone_id
          ))
        )
      }))
      
      output$summary_matrix <- DT::renderDataTable({
        DT::datatable(summary.stats, options = list(scrollX = TRUE, pageLength = 10))
      })
    }, error = function(e)
    {
      if (e$message != "")
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
                   if (e$message != "")
                   {
                     message.error <- paste0("Error rendering sankey plot: ", e$message)
                     output$genepairsplots <- renderUI(NULL)
                     render_error_ui(message = message.error, output)
                     #output$errorreport <- renderText(message.error)
                     message(message.error)
                   }
                 })
               })
  
  
  #Run Fixed Radius Neighbourhoods
  observeEvent(c(input$run_fixn,input$hierarch_tabs), {
    render_error_ui("", output = output) #Clear any error messages
    tryCatch({
      withProgress(message = paste0("Run fixed neighbourhoods with ", input$epitope_sel_n, ":"), value = 0, {
      tr <- analysis_results()
      req(tr)
      req(input$epitope_sel_n)
      rep_diff <- import("tcrdist.rep_diff")
      incProgress(0.1,"Preparing data...")
      np <- import("numpy")
      chain <- input$chain_select_knn
      # diff testing is pasted on binary comparison, so all epitope not 'PA' are set to 'X'
      tr$clone_df[, input$epitope_sel_n] = ifelse(tr$clone_df$epitope ==
                                                    input$epitope_sel_n,
                                                  input$epitope_sel_n,
                                                  "X")
      pwmat <- r_to_py(tr[[paste0("pw_", chain)]])
      
      incProgress(0.5,"Running neighbourhood differential...")
      # Larger Radius
      nn_df = rep_diff$neighborhood_diff(
        clone_df = tr$clone_df,
        pwmat = pwmat,
        count_col = 'count',
        x_cols = r_to_py(list(input$epitope_sel_n)),
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
      output$fixedradius_output_ui <-renderUI({
        withSpinner(DT::dataTableOutput("fixedradius_output"))
      })
      incProgress(0.9,"Displaying results table...")
      
      session.id<-session$token
      file.name=paste0('TCRdist3-KNN-clustering-',input$epitope_sel_n,"-",chain,"-",session.id)
      output$fixedradius_output<-DT::renderDT({
        DT::datatable(
          nn_df1,
          caption = "KNN neighbours",
          rownames = F,
          extensions = 'Buttons',
          options = list(
            scrollX = TRUE,
            pageLength = 10,
            dom = 'lfrtipB',
            buttons  = list('copy', list(extend='csv',text="csv",filename= file.name),
                            list(extend="excel",text='excel',filename=file.name), 
                            list(extend='pdf',text='pdf',filename=file.name), 
                            list(extend='print',text="print",filename=file.name))
          )
          )
      }, server = F)
      })
    }, error = function(e)
    {
      if (e$message != "")
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
    hierarch_tabs.check <<- F
    #neighbour_results()
  })
  
  #Listen to hierarchical neighbourhood clustering events
  #neighbour_results <- eventReactive(
  observeEvent(
      input$run_Hierarch,
    {
      render_error_ui("", output = output) #Clear any error messages
      message("Running neighbouring")
      req(!hierarch_tabs.check)
      hierarch_tabs.check <<- T
      withProgress(message = 'Running hierarchical neighbourhood clustering...', value = 0, {
        tryCatch({
          tr <- analysis_results()
          req(tr)
          req(input$epitope_sel_h)
          epitope_h<-input$epitope_sel_h
          rep_diff <- import("tcrdist.rep_diff")
          hierdiff <- import("hierdiff")
          pd <- import("pandas")
          chains <- strsplit(input$chains, ",")[[1]]
          incProgress(0.1,"Run hierarchical neighbourhoods")
          # diff testing is pasted on binary comparison, so all epitope not 'PA' are set to 'X'
          tr$clone_df[, epitope_h] = ifelse(tr$clone_df$epitope ==
                                                        epitope_h,
                                                      epitope_h,
                                                      "X")
          pwmat <- r_to_py(tr[[paste0("pw_", chains[1])]])
          chain <- input$chain_select_neighbours
          #tr$pw_beta
          pwmat <- tr[[paste0("pw_", chain)]]
          res_Z <- rep_diff$hcluster_diff(tr$clone_df,
                                          pwmat,
                                          x_cols = r_to_py(list(epitope_h)),
                                          count_col = 'count')
          res <- res_Z[[1]]
          Z <- res_Z[[2]]
          
          res_summary <- rep_diff$member_summ(
            res_df = res,
            clone_df = tr$clone_df,
            addl_cols = r_to_py(list('epitope'))
          )
          res_detailed = r_to_py(cbind(res, res_summary))
          incProgress(0.3, detail = "Preparing results...")
          res_detailed.df <- py_to_r(res_detailed)
          columns.tooltip <- colnames(res_detailed.df)[(grep("FDRq", colnames(res_detailed.df)) +
                                                          1):ncol(res_detailed.df)]
          message(paste0(columns.tooltip, ";"))
          incProgress(0.5,"Generating cluster diagram...")
          html = hierdiff$plot_hclust_props(
            Z,
            title = paste0(epitope_h, ' Epitope & ', chain, "-chain"),
            res = res_detailed,
            tooltip_cols = r_to_py(columns.tooltip),
            alpha = 0.00001,
            colors = r_to_py(c('blue', 'gray')),
            alpha_col = 'pvalue'
          )
          incProgress(0.9, detail = "Rendering results table...")
          session.id<-session$token
          file.name=paste0('TCRdist3-hierarchical-clustering-',epitope_h,"-",chain,"-",session.id)
          output$hierachy_output <- DT::renderDataTable({
            DT::datatable(
              py_to_r(res_detailed),
              caption = "Hierarchical neighbours",
              extensions = 'Buttons',
              rownames = F,
              options = list(
                scrollX = TRUE,
                pageLength = 10,
                dom = 'lfrtipB', #'flrtBip',
                buttons  = list('copy', list(extend='csv',text="csv",filename= file.name),
                                list(extend="excel",text='excel',filename=file.name), 
                                list(extend='pdf',text='pdf',filename=file.name), 
                                list(extend='print',text="print",filename=file.name))
              )
            )
            
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
          neighbour_results<<-res_detailed.df
        }, error = function(e)
        {
          if (e$message != "")
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
    }
  )
  
  
  #Generate Network Graph that will
  observeEvent(input$chain_select_neighbours, {
    tr <- analysis_results()
    
    tryCatch({
      req(tr)
      chain <- input$chain_select_neighbours
      genes <- tr$clone_df[[paste0("v_",substr(chain,1,1),"_gene")]]
      chain2check <- paste0("cdr3_", substring(chain, 1, 1), "_aa")
      pw_matrix = tr[[paste0("pw_", chain2check)]]
      updateSelectizeInput(session, "v_gene_filter", 
                           choices = c("All",genes), 
                           server = TRUE) # This is the magic argument
      output$dist_threshold_ui<-renderUI({
        tagList(
          sliderInput(
            "dist_threshold",
            "Max Distance to Show:",
            min = floor(min(pw_matrix)),
            max = ceiling(max(pw_matrix)),
            value = ceiling(median(pw_matrix))
          ),  
        )
      })
    }, error = function(e) {
      if (e$message != "")
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
      chain.char<-substr(input$chain_select_neighbours,1,1)
      #chain2check <- paste0("cdr3_",chain.char, "_aa")
      pw_matrix = tr[[paste0("pw_","cdr3_",chain.char, "_aa")]]
      message("Chain starting character:", chain.char)
      knn_edges<-tcr_metadata[,c("cid",paste0("cdr3_",chain.char,"_aa"),"epitope",paste0("v_",chain.char,"_gene"),"neighbors")] %>% unnest(neighbors)
      colnames(knn_edges)<-c("from","cdr3_aa","epitope","v_gene","to")
      rownames(knn_edges)<-NULL
      knn_edges<-knn_edges[!duplicated(knn_edges),]
      knn_edges<-knn_edges[order(knn_edges$from,decreasing = F),]
      View(knn_edges)
      tcr.data<-tcr_data()
      clone.df<-tr$clone_df
      save(clone.df,knn_edges,tcr.data,nb,pw_matrix,file="Test.RData")
      knn_edges$distance <- as.numeric(apply(knn_edges,1,function(x) pw_matrix[x[["from"]],x[["to"]]]))
      View(knn_edges)
      stop()
      tcr_metadata<-tcr_metadata[,c("cid",paste0("cdr3_",chain.char,"_aa"),paste0("v_",chain.char,"_gene"),"epitope")]
      colnames(tcr_metadata)<-c("clone_id","cdr3_aa","v_gene","epitope")
      nodes <- tcr_metadata
      edges <- knn_edges
      
      # Filter nodes by V-Gene if selected
      if (input$v_gene_filter != "All") {
        valid_clones <- nodes %>% filter(grepl(input$v_gene_filter,v_gene)) %>% pull(clone_id)
        
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
      if (e$message != "")
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
  
  #Populate trees
  observeEvent(input$run_trees, {
    render_error_ui("", output = output) #Clear any error messages
    withProgress(message = 'Running trees...', value = 0, {
      tryCatch({
        tr <- analysis_results()
        req(tr)
        incProgress(0.1, detail = "Importing required library...")
        TCRtree <- import("tcrdist.tree")
        incProgress(0.3, detail = "Preparing data...")
        file.name <- paste0(session$token, "_tree.html")
        temp.file <- tempfile()
        message("Temp file for trees: ", temp.file)
        remove_www_tempfiles(session, "tree.html")
        tcrtree = TCRtree$TCRtree(tcrrep = tr, html_name = temp.file)
        incProgress(0.5, detail = "Building tree...")
        tcrtree$build_tree()
        incProgress(0.9, detail = "Finalizing...")
        file.copy(temp.file, paste0("www/", file.name))
        output$tree_output_plot <- renderUI({
          tags$iframe(
            src = file.name,
            width = "95%",
            style = "height: 90vh;",
            # Use CSS for full view port height
            frameBorder = "0"
          )
        })
      }, error = function(e)
      {
        message.error <- paste0("Error computing trees: ", e$message)
        output$genepairsplots <- renderUI(NULL)
        render_error_ui(message = message.error, output)
        #output$errorreport <- renderText(message.error)
        message(message.error)
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
  
  
  
}
