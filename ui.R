#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

#source("global.R")

# --- UI Definition ---
ui <- page_sidebar(
  title = "TCRdist3 Shiny",
  useShinyalert(),
  theme = bs_theme(
    bootswatch = "cosmo",
    #"darkly", #cosmo","flatly", "darkly", "journal"
    base_font = font_google("Roboto")
  ),
  #fluidPage(titlePanel("TCRdist3 Analysis"),
  sidebar = sidebar(
    # sidebarLayout(
    #   sidebarPanel(
    h4("1. Upload Data"),
    fileInput(
      "file1",
      NULL,
      placeholder = "Choose CSV/TSV File",
      multiple = FALSE,
      accept = c(
        "text/csv",
        "text/comma-separated-values,text/plain",
        ".csv",
        ".tsv"
      )
    ),
    
    hr(style = "border-top: 1px solid #7f8c8d;"),
    
    h4("2. Set Parameters"),
    selectInput(
      "organism",
      "Organism:",
      choices = c("human", "mouse"),
      selected = "human"
    ),
    
    selectInput(
      "chains",
      "TCR Chains:",
      choices = c('alpha,beta', 'alpha', 'beta', 'gamma,delta', 'gamma', 'delta'),
      selected = "alpha,beta"
    ),
    
    selectInput(
      "dbfile",
      "DB file:",
      choices = c(
        'alphabeta_gammadelta_db.tsv',
        'alphabeta_db.tsv',
        'gammadelta_db.tsv'
      ),
      selected = "alphabeta_gammadelta_db.tsv"
    ),
    
    conditionalPanel(
      condition = "output.csvtest",
      hr(),
      h4("3. Run Analysis"),
      actionButton(
        "run",
        HTML(
          "<span class='glyphicon glyphicon-play-circle'></span> Run TCRdist3"
        )
      ),
    ),
    hr(),
    h4("Input File Format"),
    p(
      "The input file must be a CSV or TSV with the following columns for single-chain analysis: 'v_b_gene', 'j_b_gene', 'cdr3_b_aa'."
    ),
    p(
      "For paired-chain analysis ('alpha,beta'), it must also include: 'v_a_gene', 'j_a_gene', 'cdr3_a_aa'."
    ),
    p("A 'count' column is optional."),
    width = 350,
    # Wider than default (default is ~250px)
    #bg = "#2c3e50",       # Dark blue-grey background
    #fg = "white",         # White text
  ),
  
  card(
    tabsetPanel(
      tabPanel(
        "TCR Distance",
        conditionalPanel(condition = "output.csvtest", tabsetPanel(
          tabPanel("Summary", tags$br(), 
                   DT::dataTableOutput("summary_matrix")),
          tabPanel(
            "Heatmaps",
            uiOutput("matrixoptions"),
            hr(),
            withSpinner(
              plotOutput("heatmap", height = "900px"),
              type = 6,
              color = "#75AADB"
            ),
            hr(),
            DT::dataTableOutput("dist_matrix")
          ),
          tabPanel(
            "Clones",
            tags$br(),
            withSpinner(DT::dataTableOutput("clone_matrix")),
            hr(),
            downloadButton("download_data", "Download Clone Matrix")
          )
        ))
      ),
      tabPanel(
        "Gene Level Analysis",
        hr(),
        uiOutput("col_select_ui"),
        uiOutput("btn_gensankey_ui"),
        hr(),
        htmlOutput("genepairsplots")
      ),
      tabPanel(
        "Neighbourhoods",
        hr(),
        # The Main Panel
        tabsetPanel(
          tabPanel(
            "Fixed Radius",
            hr(),
            layout_sidebar(
              # The Sidebar
              sidebar = sidebar(
                h5("Fixed Radius"),
                uiOutput("epitope_select_ui"),
                uiOutput("chain_select_ui_knn"),
                numericInput("knn_radius", label = "KNN radius", value = 150),
                br(),
                actionButton(
                  "run_fixn",
                  HTML(
                    "<span class='glyphicon glyphicon-play-circle'></span> Run KNN"
                  )
                )
              ),
              card(
                uiOutput("fixedradius_output_ui"),
                full_screen = TRUE 
                #withSpinner(DT::dataTableOutput("fixedradius_output")),
              )
            )
          ),
          tabPanel(
            "Hierarchical",
            hr(),
            layout_sidebar(
              sidebar = sidebar(
                h5("Hierarchical"),
                uiOutput("epitope_select_hierarch_ui"),
                uiOutput("chain_select_ui"),
                actionButton(
                  "run_Hierarch",
                  HTML(
                    "<span class='glyphicon glyphicon-play-circle'></span> Run Hierarchical"
                  )
                )
                # ,hr(),
                # h5("Graph options"),
                # selectizeInput("v_gene_filter", "Select V gene:", choices = NULL),
                # selectInput(
                #   "epitope_color",
                #   "Color Nodes By:",
                #   choices = c("v_gene", "epitope")
                # ),
                # uiOutput("dist_threshold_ui"),
                # selectInput(
                #   "color",
                #   "Choose Color:",
                #   choices = c(
                #     "Blue" = "steelblue",
                #     "Red" = "firebrick",
                #     "Green" = "forestgreen"
                #   )
                # )
              ),
              card(
                tabsetPanel(
                  id = "hierarch_tabs",
                  tabPanel("Cluster diagram", uiOutput("hierachy_output_plot")),
                  tabPanel("Table", DT::dataTableOutput("hierachy_output"))
                  # tabPanel(
                  #   "Network Graph",
                  #   helpText("Select a node on the graph to see details below."),
                  #   width = 10,
                  #   visNetworkOutput("knnGraph", height = "600px")
                  # )
                ),
                full_screen = TRUE 
              )
            )
          )
        )
      ),
      tabPanel(
        "Trees",
        hr(),
        actionButton(
          "run_trees",
          HTML(
            "<span class='glyphicon glyphicon-play-circle'></span> Get Trees"
          )
        ),
        hr(),
        uiOutput("tree_output_plot"),
        hr(),
      ),
      tabPanel("Export Results", hr(), p("")),
      tabPanel(
        "About",
        h4("About this App"),
        p(
          "This Shiny app provides an interface to the `tcrdist3` Python package for TCR sequence analysis."
        ),
        p(
          "`tcrdist3` is a powerful tool for calculating distances between TCRs, which can be used to identify clusters of functionally related receptors."
        ),
        p(
          "For more information on `tcrdist3`, please see the official documentation: ",
          a("https://tcrdist3.readthedocs.io/", href = "https://tcrdist3.readthedocs.io/")
        ),
        #h4("How to Cite"),
        # p("If you use this tool in your research, please cite the `tcrdist3` publication:"),
        # p("Mayer-Blackwell, K., Schattgen, S., Cohen-Lavi, L., Crawford, J. C., Souquette, A., Gaevert, J., ... & Bradley, P. (2021).
        #   Flexible distance-based T-cell receptor analysis in Python with tcrdist3. bioRxiv.")
      )
    ),
    htmlOutput("errorreport"),
    #useShinyFeedback(),
    #width = 9
  )
)#)
