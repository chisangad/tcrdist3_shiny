#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
source("global.R")
# A reusable function to render a styled error message

#Declare global variable
hierarch_tabs.check<-F

render_error_ui <- function(message,output) {
  output$errorreport <- renderUI(HTML(
    ifelse(message=="","",paste0(
      '<div style="color: red; font-weight: bold; padding: 10px; border: 1px solid red; border-radius: 5px; background-color: #ffe6e6;">',
      message,
      '</div>'
    ))
  ))
  #feedbackDanger("my_text_input", TRUE, messa)
  #showNotification(message)
}

# Required column names
REQUIRED_COLS <- c(
  "v_a_gene", "j_a_gene", "cdr3_a_aa", 
  "v_b_gene", "j_b_gene", "cdr3_b_aa"
)

# --- Server Logic ---
server <- function(input, output, session) {
  
  output$csvtest <- reactive({
    !is.null(input$file1)
  })
  outputOptions(output, "csvtest", suspendWhenHidden = FALSE)
  # Reactive value to store the TCR data
  tcr_data <- reactive({
    render_error_ui("",output = output) #Clear any error messages
    req(input$file1)
    message("Reading file")
    # Determine separator based on file extension
    sep <- ifelse(tools::file_ext(input$file1$name) == "tsv", "\t", ",")
    df<-read.csv(input$file1$datapath, sep = sep)
    message(all(REQUIRED_COLS %in% colnames(df)))
    # Validation checks
    if(!all(REQUIRED_COLS %in% colnames(df)))
    {
      show_toast(title = "Error!",
                 text = paste0("Required columns not found in the file. Should include: ",paste0(REQUIRED_COLS,collapse = ", ")),
                 position = "center",
                 timer = NULL,
                 type = "error")
      return()
    }
    # If all checks pass, return the data frame
    return(df)
    
  })
  
  
  # Run tcrdist3 analysis when the run button is clicked
  analysis_results <- eventReactive(input$run, {
    render_error_ui("",output = output) #Clear any error messages
    req(tcr_data())
    hierarch_tabs.check<<-F
    render_error_ui("",output = output)
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
        message(paste0(chains,collapse = ":"))
        
        if(length(chains)==1) chains <-list(chains)
        
        incProgress(0.5, detail = "Computing distances...")
        message(r_to_py(chains))
        message(input$organism)
        message(input$dbfile)
        # Initialize TCRrep object and compute distances
        tr <- TCRrep(cell_df = cell_df,organism = r_to_py(input$organism),chains = r_to_py(chains),db_file = r_to_py(input$dbfile))
       #message(paste0(names(tr),collapse = "\n"))
        incProgress(0.8, detail = "Finalizing...")
        #Get chains from the distance compute
        all.tr_names<-names(tr)
        print(all.tr_names)
        chains_tr<-all.tr_names[grepl("^pw_",all.tr_names)]
        #chains_tr<-gsub("pw_","",chains_tr)#[!grepl("cdr|pmhc",chains_tr)])
        
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
        epitopes<-unique(tr$clone_df$epitope)
        genes<-names(tr$clone_df)[grepl("gene",names(tr$clone_df))]
        
        #Load selection input for epitopes
        output$col_select_ui <- renderUI({
          tagList(
            selectizeInput("epitope_sel", "Select Epitope:",
                             choices = epitopes),
            selectizeInput("columns_sel", "Select Columns:",
                             choices = genes,
                              multiple=T)
          )
        })
        
        #Load epitopes in Neighbourhood
        output$epitope_select_ui<-renderUI({
          tagList(
          selectizeInput("epitope_sel_n", "Select Epitope:",
                         choices = epitopes)
          )
        })
        
        return(tr)
      }, error=function(e)
      {
        if(!is.null(e$message))
        {
          message.error<-paste0("An error occurred: ", e$message)
          render_error_ui(message.error,output)
          message(message.error)
        }
      })})
  })
  
  #Update heatmap and matrix table on the following events
  observeEvent(c(input$run,input$matrix_select), {
    render_error_ui("",output = output) #Clear any error messages
    tryCatch({
      results <- analysis_results()
      req(results)
      render_error_ui("",output)
      chain2check<-input$matrix_select #paste0("cdr3_",substring(input$matrix_select,1,1),"_aa")
      if(is.null(input$matrix_select))
        chain2check<-
      message(chain2check)
      pw_matrix = results[[chain2check]]
      clone_df = results$clone_df
      mat <- pw_matrix
      chain.col_name<-gsub("^pw_","",chain2check)
      if(chain2check %in% c("pw_alpha","pw_beta","pw_gamma"))
        chain.col_name<-paste0("cdr3_",substring(gsub("pw_","",chain2check),1,1),"_aa")
      rownames(mat) <- clone_df[[chain.col_name]]
      colnames(mat) <- clone_df[[chain.col_name]]
      
      output$clone_matrix <- DT::renderDataTable({
        DT::datatable(clone_df, 
                      caption = "Summary of TCRdist clone df",
                      extensions = 'Buttons',
                      options = list(
                        scrollX = TRUE, 
                        pageLength = 10,
                       dom = 'flrtip',
                       buttons =list(
                         list(extend = 'csv', text = 'Download CSV'),
                         list(extend = 'excel', text = 'Download Excel')
                       )))
      }, server = F)
      # Display the distance matrix
      output$dist_matrix <- DT::renderDataTable({
        DT::datatable(mat, options = list(scrollX = TRUE, pageLength = 10))
      })
      
      # Display the heatmap
      output$heatmap <- renderPlot({
        pheatmap(mat,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "ward.D2",
                 main = "TCRdist Heatmap")
      })
    },error=function(e)
    {
      if(e$message!="")
      {
        message.error<-paste0("Error rendering heatmap: ", e$message)
        render_error_ui(message.error,output)
      }

    })
  })
  #Add sankey gen button to the dashboard
  observeEvent(c(input$epitope_sel,input$columns_sel),{
    render_error_ui("",output = output) #Clear any error messages
    req(input$epitope_sel)
    output$btn_gensankey_ui<-renderUI(
      actionButton("btn_gensankey",
                   HTML(
                     "<span class='glyphicon glyphicon-play-circle'></span> Generate Sankey plot"
                   ))
    )
  })

  #Display summary stats  
  observeEvent(input$run,{
    render_error_ui("",output = output) #Clear any error messages
    tr<-analysis_results()
    tryCatch({
      req(tr)
      clone_df<-py_to_r(tr$clone_df)
      summary.stats<-do.call(rbind, lapply(unique(clone_df$epitope), function(x)
      {
         data.frame(epitope=x,
                    num_individuals=length(unique(subset(clone_df, epitope==x)$subject)),
                    num_clones=length(unique(subset(clone_df, epitope==x)$clone_id)))
      }))
      
      output$summary_matrix <- DT::renderDataTable({
        DT::datatable(summary.stats, options = list(scrollX = TRUE, pageLength = 10))
      })
    }, error=function(e)
      {
      if(e$message!="")
      {
        message.error<-paste0("Error loading summary results: ", e$message)
        output$genepairsplots<-renderUI(NULL)
        render_error_ui(message = message.error,output)
        message(message.error)
      }
    })
  })
  
  #Generate sankey plot
  observeEvent(c(input$epitope_sel,input$btn_gensankey,input$columns_sel), {
    render_error_ui("",output = output) #Clear any error messages
    tryCatch({
      tr <- analysis_results()
      req(tr)
      req(input$columns_sel)
      req(length(input$columns_sel)>1)
      render_error_ui(NULL,output)
      plotting <- import("tcrdist.plotting")
      svg_PA  = plotting$plot_pairings(tr$clone_df[tr$clone_df$epitope == input$epitope_sel,],cols = r_to_py(input$columns_sel),count_col='count')
      output$genepairsplots<-renderUI(HTML(svg_PA))
    },error=function(e)
    {
      if(e$message!="")
      {
        message.error<-paste0("Error rendering sankey plot: ", e$message)
        output$genepairsplots<-renderUI(NULL)
        render_error_ui(message = message.error,output)
        #output$errorreport <- renderText(message.error)
        message(message.error)
      }
    })
  })
  
  #Run Fixed Radius Neighbourhoods
  observeEvent(input$run_fixn,{
    render_error_ui("",output = output) #Clear any error messages
    tryCatch({
      tr<-analysis_results()
      req(tr)
      req(input$epitope_sel_n)
      rep_diff <- import("tcrdist.rep_diff")
      message("Run fixed neighbourhoods:")
      np<-import("numpy")
      chains <- strsplit(input$chains, ",")[[1]]
      # diff testing is pasted on binary comparison, so all epitope not 'PA' are set to 'X'
      tr$clone_df[,input$epitope_sel_n] = ifelse(tr$clone_df$epitope==input$epitope_sel_n,input$epitope_sel_n,"X")
      pwmat<-r_to_py(tr[[paste0("pw_",chains[1])]])

      for(chain in chains[-1])
      {
        message(chain)
        pwmat<-pwmat+tr[[paste0("pw_",chain)]]
      }
      # Larger Radius
      nn_df = rep_diff$neighborhood_diff(clone_df = tr$clone_df, pwmat =pwmat, count_col = 'count', x_cols = r_to_py(list(input$epitope_sel_n)),knn_radius = r_to_py(input$knn_radius))
      nn_df1<-nn_df %>% select(c('K_neighbors', 'val_0', 'ct_0', 'val_2', 'ct_2', 'RR','OR', 'pvalue', 'FWERp','FDRq')) %>% arrange('FDRq') %>% arrange(desc(OR),desc(ct_0))
      output$fixedradius_output <- DT::renderDT({
        DT::datatable(nn_df1,caption = "KNN neighbours",
                      extensions = 'Buttons',
                      options = list(
                        scrollX = TRUE, 
                        pageLength = 10,
                        dom = 'Bflrtip',
                        buttons =list(
                          list(extend = 'csv', text = 'Download CSV'),
                          list(extend = 'excel', text = 'Download Excel')
                        )))
      },server = F)
      
    }, error=function(e)
    {
      if(e$message!="")
      {
        message.error<-paste0("Error computing fixed radius neighbourhoods: ", e$message)
        output$genepairsplots<-renderUI(NULL)
        render_error_ui(message = message.error,output)
        #output$errorreport <- renderText(message.error)
        message(message.error)
      }
    })
  })
  
  observeEvent(input$run_Hierarch, {
    hierarch_tabs.check<<-F
  })
  
  #Listen to hierarchical neighbourhood clustering events
  observeEvent(c(input$hierarch_tabs, input$epitope_sel_n, input$run_Hierarch), {
    render_error_ui("",output = output) #Clear any error messages
    req(!hierarch_tabs.check)
    hierarch_tabs.check<<-T
    withProgress(message = 'Running hierarchical neighbourhood clustering...', value = 0, {
    tryCatch(
      {
        tr<-analysis_results()
        req(tr)
        req(input$epitope_sel_n)
        rep_diff <- import("tcrdist.rep_diff")
        hierdiff <- import("hierdiff")
        pd <- import("pandas")
        chains <- strsplit(input$chains, ",")[[1]]
        message("Run hierarchical neighbourhoods:")
        # diff testing is pasted on binary comparison, so all epitope not 'PA' are set to 'X'
        tr$clone_df[,input$epitope_sel_n] = ifelse(tr$clone_df$epitope==input$epitope_sel_n,input$epitope_sel_n,"X")
        pwmat<-r_to_py(tr[[paste0("pw_",chains[1])]])
        
        for(chain in chains[-1])
        {
          message(chain)
          pwmat<-pwmat+tr[[paste0("pw_",chain)]]
        }
        res_Z<-rep_diff$hcluster_diff(tr$clone_df, tr$pw_beta, x_cols = r_to_py(list(input$epitope_sel_n)), count_col = 'count')
        res<-res_Z[[1]]
        Z<-res_Z[[2]]

        res_summary <- rep_diff$member_summ(res_df = res, clone_df = tr$clone_df, addl_cols=r_to_py(list('epitope')))
        res_detailed = r_to_py(cbind(res, res_summary))
        incProgress(0.3, detail = "Preparing results...")
        res_detailed.df<-py_to_r(res_detailed)
        columns.tooltip<-colnames(res_detailed.df)[(grep("FDRq",colnames(res_detailed.df))+1):ncol(res_detailed.df)]
        message(paste0(columns.tooltip,";"))     
        html = hierdiff$plot_hclust_props(Z,title=paste0(input$epitope_sel_n,' Epitope'),res=res_detailed,tooltip_cols=r_to_py(columns.tooltip),alpha=0.00001, colors = r_to_py(c('blue','gray')),alpha_col='pvalue')
        incProgress(0.9, detail = "Rendering results table...")
        output$hierachy_output <- DT::renderDataTable({
          DT::datatable(py_to_r(res_detailed),caption = "Hierarchical neighbours",
                        extensions = 'Buttons',
                        options = list(
                          scrollX = TRUE, 
                          pageLength = 10,
                          dom = 'Bflrtip',
                          buttons =list(
                            list(extend = 'csv', text = 'Download CSV'),
                            list(extend = 'excel', text = 'Download Excel')
                          )))
        
        }, server = F)
        #message(html)
        incProgress(0.95, detail = "Save html file...")
        if(!dir.exists("www"))
          dir.create("www")
        html.file<-paste0(session$token,"_Test_hierachical.html")
        write(html,file.path("www/",html.file))
        output$hierachy_output_plot<-renderUI(
          tags$iframe(
            src = html.file,
            width = "95%",
            style = "height: 90vh;", # Use CSS for full viewport height
            frameBorder = "0"
          )
        )
        
      }, error=function(e)
      {
        if(e$message!="")
        {
          message.error<-paste0("Error computing hierarchical neighbourhood clusters: ", e$message)
          output$genepairsplots<-renderUI(NULL)
          render_error_ui(message = message.error,output)
          #output$errorreport <- renderText(message.error)
          message(message.error)
        }
      }
    )})
  })
  

  #Populate trees
  observeEvent(input$run_trees,{
    render_error_ui("",output = output) #Clear any error messages
    withProgress(message = 'Running trees...', value = 0, {
      tryCatch({
        tr<-analysis_results()
        req(tr)
        incProgress(0.1, detail = "Importing libraries...")
        TCRtree<-import("tcrdist.tree")
        incProgress(0.3, detail = "Preparing data...")
        file.name<-paste0(session$token,"dash.mouse.b.tree.html")
        tcrtree = TCRtree$TCRtree(tcrrep = tr, html_name = paste0("www/",file.name))
        incProgress(0.5, detail = "Building tree...")
        tcrtree$build_tree()
        incProgress(0.9, detail = "Finalizing...")
        output$tree_output_plot<-renderUI({
          tags$iframe(
            src = file.name,
            width = "95%",
            style = "height: 90vh;", # Use CSS for full viewport height
            frameBorder = "0"
          )
        })
      }, error=function(e)
        {
        message.error<-paste0("Error computing trees: ", e$message)
        output$genepairsplots<-renderUI(NULL)
        render_error_ui(message = message.error,output)
        #output$errorreport <- renderText(message.error)
        message(message.error)
    })})
  })
  # Download handler for the download button
  output$download_data <- downloadHandler(
    filename = function() {
      paste("my-clone-data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      tr=analysis_results()
      req(tr)
      write.csv(tr$clone_df, file, row.names = FALSE)
    }
  )

  #Clean-up temp files when session ends
  onStop(function() {
    temp.files<-list.files("www",pattern = paste0(session$token,".*"),full.names = T)
    message(paste0(temp.files,collapse = ";"))
    for(temp_file in temp.files)
    {
      if (file.exists(temp_file)) {
        file.remove(temp_file)
        message("Temporary file removed: ", temp_file)
      }
    }

  })
}
