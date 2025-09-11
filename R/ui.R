#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

source("global.R")

# --- UI Definition ---
ui <- fluidPage(
  titlePanel("TCRdist3 Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Upload Data"),
      fileInput("file1",NULL,
                placeholder = "Choose CSV/TSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",
                           ".tsv")),
      
      hr(),
      
      h4("2. Set Parameters"),
      selectInput("organism", "Organism:",
                  choices = c("human", "mouse"),
                  selected = "human"),
      
      selectInput("chains", "TCR Chains:",
                  choices = c('alpha,beta','alpha', 'beta', 'gamma,delta','gamma','delta'),
                  selected = "alpha,beta"),
      
      selectInput("dbfile", "DB file:",
                  choices = c('alphabeta_gammadelta_db.tsv','alphabeta_db.tsv','gammadelta_db.tsv'),
                  selected = "alphabeta_gammadelta_db.tsv"),
      
      conditionalPanel(
        condition = "output.csvtest",
        hr(),
        h4("3. Run Analysis"),
        actionButton("run", "Run TCRdist3"),
      ),
      hr(),
      h4("Input File Format"),
      p("The input file must be a CSV or TSV with the following columns for single-chain analysis: 'v_b_gene', 'j_b_gene', 'cdr3_b_aa'."),
      p("For paired-chain analysis ('alpha,beta'), it must also include: 'v_a_gene', 'j_a_gene', 'cdr3_a_aa'."),
      p("A 'count' column is optional.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("TCR Distance",
                 conditionalPanel(condition = "output.csvtest",
                 tabsetPanel(
                   tabPanel(
                     "Summary",
                     tags$br(),
                     DT::dataTableOutput("summary_matrix")
                   ),
                   tabPanel(
                     "Clones",
                     tags$br(),
                     DT::dataTableOutput("clone_matrix"),
                     hr(),
                     downloadButton("download_data", "Download Clone Matrix")
                   ),
                   tabPanel(
                     "Heatmap", 
                      uiOutput("matrixoptions"),
                      hr(),
                      plotOutput("heatmap", height = "800px"),
                      hr(),
                      DT::dataTableOutput("dist_matrix")
                      )
                 ))),
        tabPanel("Gene Level Analyses",
                 hr(),
                 p("Loading a TCR dataset will immediately enable gene-level analyses of repertoires, for example Sankey plots of V/D/J gene frequency or statistical tests for differential gene enrichment in two or more conditions:"),
                 hr(),
                 uiOutput("col_select_ui"),
                 uiOutput("btn_gensankey_ui"),
                 hr(),
                 htmlOutput("genepairsplots")
        ),
        tabPanel("Neighbourhoods",
                 hr(),
                 uiOutput("epitope_select_ui"),
                 hr(),
                tabsetPanel(
                   tabPanel(
                     "Fixed Radius Neighbourhoods",
                     hr(),
                     numericInput("knn_radius",label = "KNN radius",value = 150),
                     br(),
                     actionButton("run_fixn","Run KNN"),
                     hr(),
                     DT::dataTableOutput("fixedradius_output"),
                     hr()
                   ),
                  tabPanel(
                     "Hierarchical Neighbourhoods",
                     hr(),
                     actionButton("run_Hierarch","Run Hierarchical"),
                     hr(),
                     htmlOutput("hierachy_output_plot"),
                     DT::dataTableOutput("hierachy_output")
                   )
                  )
                 ),
        tabPanel(
          "Trees",
          hr(),
          actionButton("run_trees","Get Trees"),
          hr(),
        ),
        tabPanel("Export Results",
                 hr(),
                 p("")
        ),
        # tabPanel("About", 
        #          h4("About this App"),
        #          p("This Shiny app provides an interface to the `tcrdist3` Python package for TCR sequence analysis."),
        #          p("`tcrdist3` is a powerful tool for calculating distances between TCRs, which can be used to identify clusters of functionally related receptors."),
        #          p("For more information on `tcrdist3`, please see the official documentation: ", 
        #            a("https://tcrdist3.readthedocs.io/", href = "https://tcrdist3.readthedocs.io/")),
        #          h4("How to Cite"),
        #          p("If you use this tool in your research, please cite the `tcrdist3` publication:"),
        #          p("Mayer-Blackwell, K., Schattgen, S., Cohen-Lavi, L., Crawford, J. C., Souquette, A., Gaevert, J., ... & Bradley, P. (2021). 
        #            Flexible distance-based T-cell receptor analysis in Python with tcrdist3. bioRxiv.")
        # )
      ),
      htmlOutput("errorreport")
    )
  )
)

