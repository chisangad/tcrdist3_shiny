#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

#source("global.R")

# --- UI Definition ---
ui <- page_sidebar(
  #use_tracking(),
  title = "TCRdist3 Shiny",
  theme = bs_theme(
    bootswatch = "cosmo",
    #"darkly", #cosmo","flatly", "darkly", "journal"
    base_font = font_google("Roboto")
  ),
  # tags$head(
  #   tags$style(HTML("
  #     /* Target the selected items (tags) inside the selectize input */
  #     .selectize-input .item {
  #       background-color: orange !important;
  #       background-image: none !important; /* Removes the default blue gradient */
  #       color: black !important;           /* Sets text color */
  #       border: 1px solid darkorange !important; /* Optional: adds a matching border */
  #     }
  #   "))),
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
      selected = "mouse"
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
        "alpha beta & gammadelta"='alphabeta_gammadelta_db.tsv',
        adaptive_imgt_mapping='adaptive_imgt_mapping.csv',
        combo_db='combo_db.tsv',
        "combo_xcr_2023-12-30"='combo_xcr_2023-12-30.tsv',
        "combo_xcr_2024-03-05"='combo_xcr_2024-03-05.tsv'
        # 'alphabeta_gammadelta_db.tsv',
        # 'alphabeta_db.tsv',
        # 'gammadelta_db.tsv'
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
      ##############################################################
      #.                         Summary Panel                     #
      ##############################################################
      tabPanel(
        "Summary",
        tags$br(),
        layout_columns(
          col_widths = c(4, 8),
          plotOutput("summaryPlot"),
          DT::dataTableOutput("summary_matrix")
        ),
        
        uiOutput("networkcard")
      ),
      ###########################################################
      #                        TCR Distance Tab                 #
      ###########################################################
      tabPanel("TCR Distance",
               hr(),
               conditionalPanel(condition = "output.csvtest", 
                 card(tabsetPanel(
                   tabPanel(
                     "Heatmap",
                     layout_sidebar(
                       sidebar = sidebar(
                         h5("Options"),
                         hr(),
                         uiOutput("matrixoptions"),
                         uiOutput("matrixoptions2"),
                         actionButton(
                           "run_heatmap",
                           HTML(
                             "<span class='glyphicon glyphicon-play-circle'></span> Load heatmap"
                           )
                         ),
                         hr(),
                         h6("Heatmap options"),
                         checkboxInput("show_row_names", "Show Row Labels", value = F),
                         checkboxInput("show_col_names", "Show Column Labels", value = F)
                       ),
                       card(
                     withSpinner(
                       #uiOutput("heatmap_container"),
                       InteractiveComplexHeatmapOutput(
                         "heatmap",
                         layout = "1|(2-3)",
                         compact = T,
                         width1 = "800px",
                         height1 = "700px"
                       ),
                       type = 6,
                       color = "#75AADB"
                     ),
                     hr(),
                     withSpinner(DT::dataTableOutput("dist_matrix"))
                   ))),
                   #tabPanel("Distance Matrix", withSpinner(DT::dataTableOutput("dist_matrix"))),
                   tabPanel("kPCA", layout_sidebar(
                     sidebar = sidebar(
                       h5("Options"),
                       hr(),
                       uiOutput("kpca_chain_option"),
                       uiOutput("kpca_repertoire_option"),
                       actionButton(
                         "run_kpca",
                         HTML(
                           "<span class='glyphicon glyphicon-play-circle'></span> Load kPCA"
                         )
                       ),
                       hr(),
                       h6("KPCA options"),
                       checkboxInput("kpca_facet", "Show all", value = F),
                       conditionalPanel(
                         condition = "!input.kpca_facet",
                         uiOutput("kpca_colorby_option")
                       )
                     ),
                     card(
                          card_header("Dataset Overview"),
                          plotlyOutput("kpcaPlot"),
                          full_screen = TRUE, # bslib magic: adds an 'expand' icon automatically
                          card_footer("Click the expand icon in the bottom right to go full screen."))
                   )), 
                   tabPanel("Diversity", card(
                        plotOutput("diversityPlot"),
                        plotOutput("simpsonPlot"),
                        plotOutput("invSimpsonPlot")
                   )),
                   tabPanel(
                     "Clones",
                     tags$br(),
                     plotOutput("topCdr3Plot"),
                     hr(),
                     withSpinner(DT::dataTableOutput("clone_matrix")),
                     hr(),
                     downloadButton("download_data", "Download Clone Matrix")
                   )
                 ))
               )),
      ###########################################################
      #                        Gene Level Analysis              #
      ###########################################################
      tabPanel(
        "Gene Level Analysis",
        hr(),
        uiOutput("col_select_ui"),
        uiOutput("btn_gensankey_ui"),
        hr(),
        htmlOutput("genepairsplots")
      ),
      ###########################################################
      #                        Neighbourhoods                   #
      ###########################################################
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
          ###########################################################
          #                        Hierarchical Clustering Tab      #
          ###########################################################          
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
      ###########################################################
      #                        Motifs                           #
      ###########################################################      
      tabPanel(
        "Motif Enrichment",
        hr(),
        layout_sidebar(
          # The Sidebar
          sidebar = sidebar(
            h5("Options"),
            uiOutput("epitope_select_ui_motif"),
            uiOutput("chain_select_ui_motif"),
            uiOutput("motif_radius_ui"),
            uiOutput("cluster_selector_select_ui"),
            uiOutput("cluster_select_ui"),
            hr(),
            uiOutput("cluster_summary"),  
            br(),
            actionButton(
              "run_motifs",
              HTML(
                "<span class='glyphicon glyphicon-play-circle'></span> Get Motifs"
              )
            )
          ),
          card(
            tabsetPanel(
              tabPanel("Motif Logo",
                       plotOutput("motif_plot"),
                       #uiOutput("motif_svg"),
                       tableOutput("motif_table")
              ),
              tabPanel("Cluster Details",
                       h4("Sequences in this cluster"),
                       DTOutput("cluster_seq_table"),
                       h4("V gene usage"),
                       plotOutput("v_gene_plot", height = "250px"),
                       h4("J gene usage"),
                       plotOutput("j_gene_plot", height = "250px")
              )
              ),
            full_screen = TRUE 
            #withSpinner(DT::dataTableOutput("fixedradius_output")),
          )
        )
      ),
      ###########################################################
      #                        Trees                            #
      ###########################################################      
      tabPanel(
        "Trees",
        hr(),
        layout_sidebar(
          # The Sidebar
          sidebar = sidebar(
            h5("Options"),
            uiOutput("epitope_select_ui_trees"),
            uiOutput("chain_select_ui_trees"),
            br(),
            actionButton(
              "run_trees",
              HTML(
                "<span class='glyphicon glyphicon-play-circle'></span> Get Trees"
              )
            )
          ),
          card(
            conditionalPanel(
            condition="output.showtrees=='show'",
            withSpinner(uiOutput("tree_output_plot"))
            ),
            #uiOutput("tree_output_plot"),
            full_screen = TRUE 
            #withSpinner(DT::dataTableOutput("fixedradius_output")),
          )
        )
        
        # 
        # actionButton(
        #   "run_trees",
        #   HTML(
        #     "<span class='glyphicon glyphicon-play-circle'></span> Get Trees"
        #   )
        # ),
        # hr(),
        # uiOutput("tree_output_plot"),
        # hr(),
      ),
      ###########################################################
      #                        Explore Results                  #
      ###########################################################
      tabPanel(
        "Explore Results",
        hr(),
        conditionalPanel(
          condition = "output.csvtest",
          layout_sidebar(
            sidebar = sidebar(
              h5("Explorer Filters"),
              uiOutput("explore_epitope_ui"),
              uiOutput("explore_subject_ui"),
              uiOutput("explore_chain_ui"),
              uiOutput("explore_pw_matrix_ui"),
              uiOutput("explore_color_by_ui"),
              checkboxInput("explore_strip_allele", "Strip allele suffix (*01)", value = TRUE),
              checkboxInput("explore_size_by_count", "Scale point size by clone count", value = TRUE),
              width = 300
            ),
            card(
              tabsetPanel(
                id = "explore_subtabs",
                tabPanel(
                  "Distance Distributions",
                  tags$br(),
                  layout_columns(
                    col_widths = c(4, 8),
                    selectInput(
                      "explore_dist_plot_type",
                      "Distribution Plot Style:",
                      choices = c(
                        "Violin + Boxplot" = "violin",
                        "Density Curves" = "density",
                        "Empirical CDF (ECDF)" = "ecdf",
                        "Histogram" = "histogram"
                      ),
                      selected = "violin"
                    ),
                    selectInput(
                      "explore_dist_group_mode",
                      "Compare Pairs By:",
                      choices = c(
                        "Within vs. Between Epitope" = "epitope_pair_type",
                        "Within vs. Between Subject" = "subject_pair_type",
                        "By Epitope (Within-Epitope Pairs)" = "within_epitope"
                      ),
                      selected = "epitope_pair_type"
                    )
                  ),
                  withSpinner(plotlyOutput("explore_dist_plot", height = "420px")),
                  hr(),
                  h5("Component Distance Comparison (Chain / Region Correlation)"),
                  withSpinner(plotlyOutput("explore_chain_corr_plot", height = "420px"))
                ),
                tabPanel(
                  "2D / 3D Embedding (MDS)",
                  tags$br(),
                  layout_columns(
                    col_widths = c(4, 8),
                    radioButtons(
                      "explore_mds_dims",
                      "Embedding Dimensions:",
                      choices = c("2D MDS" = "2d", "3D MDS" = "3d"),
                      selected = "2d",
                      inline = TRUE
                    ),
                    helpText("Classical Multidimensional Scaling (PCoA) projects clones into 2D/3D space preserving pairwise TCR distances.")
                  ),
                  withSpinner(plotlyOutput("explore_mds_plot", height = "650px"))
                ),
                tabPanel(
                  "V-J Gene Usage",
                  tags$br(),
                  layout_columns(
                    col_widths = c(4, 4, 4),
                    selectInput(
                      "explore_vj_metric",
                      "Abundance Metric:",
                      choices = c("Unique Clones" = "clones", "Total Template Count" = "count"),
                      selected = "clones"
                    ),
                    numericInput(
                      "explore_top_n_genes",
                      "Max top genes to display:",
                      value = 15,
                      min = 5,
                      max = 50,
                      step = 5
                    ),
                    div()
                  ),
                  h5("V-J Gene Pairing Heatmap"),
                  withSpinner(plotlyOutput("explore_vj_heatmap", height = "450px")),
                  hr(),
                  h5("Top V & J Gene Distribution"),
                  withSpinner(plotlyOutput("explore_vj_barplot", height = "400px"))
                ),
                tabPanel(
                  "CDR3 Spectratype & Logo",
                  tags$br(),
                  h5("CDR3 Amino Acid Length Distribution"),
                  withSpinner(plotlyOutput("explore_cdr3_len_plot", height = "380px")),
                  hr(),
                  h5("Length-Stratified CDR3 Sequence Logo"),
                  layout_columns(
                    col_widths = c(6, 6),
                    uiOutput("explore_logo_length_ui"),
                    radioButtons(
                      "explore_logo_method",
                      "Logo Y-Axis:",
                      choices = c("Information Content (Bits)" = "bits", "Probability" = "prob"),
                      selected = "bits",
                      inline = TRUE
                    )
                  ),
                  withSpinner(plotOutput("explore_cdr3_logo", height = "280px"))
                ),
                tabPanel(
                  "Clonotype & Public TCR Explorer",
                  tags$br(),
                  h5("Nearest-Neighbor TCR Distance vs. Clonal Expansion"),
                  helpText("Clones with small nearest-neighbor distance belong to tight sequence clusters; points are shaped by whether the CDR3 is shared across multiple subjects (public) or private."),
                  withSpinner(plotlyOutput("explore_nn_plot", height = "420px")),
                  hr(),
                  h5("Enriched Clonotype Table"),
                  withSpinner(DT::dataTableOutput("explore_clones_table"))
                )
              ),
              full_screen = TRUE
            )
          )
        )
      ),
      ###########################################################
      #                        Export Results                   #
      ###########################################################
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
