# Load required packages with automatic installation
#
# Package list organized by purpose:
# - Core Shiny & web: shiny, shinyjs, shinyalert, bslib, shinycssloaders
# - Data & visualization: tidyverse (includes ggplot2, dplyr), plotly
# - Tables: DT
# - TCR-specific: ComplexHeatmap, InteractiveComplexHeatmap, kernlab
# - Python/reticulate: reticulate
# - Other: visNetwork, openxlsx

path2conda=normalizePath("~/miniconda3/envs/tcrdist3_env/bin/python")
Sys.setenv(RETICULATE_PYTHON = path2conda)

message("Checking if packages exist, install and load")
packages <- list(
  # Core dependencies
  core = c("reticulate", "shiny", "shinyjs", "shinyalert", "bslib", "shinycssloaders","shinylogs","ggseqlogo"),
  # Data wrangling & visualization
  data = c("tidyverse", "plotly", "DT"),
  # Heatmap & kernel methods
  heatmap = c("ComplexHeatmap", "kernlab","InteractiveComplexHeatmap"),
  # Network visualization
  network = c("visNetwork"),
  # Utilities
  utils = c("openxlsx")
)

# Flatten and load all packages
all_packages <- unlist(packages)
for (pkg in all_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    # Special handling for BioConductor packages
    if (pkg %in% c("InteractiveComplexHeatmap","ComplexHeatmap")) {
      message(paste0("Installing ", pkg, " from BioConductor..."))
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

message("Packages loaded!")

py_run_string("
def process_list(my_list):
        return my_list
")

# --- Python Environment Setup ---
env_name <- "tcrdist3_env"
message("Check if Conda environment exists")
setup_conda_environment <- function(env_name = "tcrdist3_env") {
  tryCatch({
    if (!reticulate::condaenv_exists(env_name)) {
      message(paste0("Conda environment '", env_name, "' not found. Creating..."))
      
      reticulate::conda_create(
        envname = env_name,
        packages = "python=3.8",
        channel = c("conda-forge", "bioconda", "defaults")
      )
      
      message(paste0("Installing tcrdist3 in '", env_name, "'..."))
      reticulate::conda_install(env_name, packages = "tcrdist3", pip = TRUE)
    }
    
    reticulate::use_condaenv(env_name)
    message(paste0("Conda environment '", env_name, "' ready."))
  }, error = function(e) {
    stop(paste0("Failed to setup conda environment: ", conditionMessage(e)))
  })
}

# Initialize conda environment
setup_conda_environment(env_name)
message("Conda environment loaded")

#' Generate kPCA Plot from tcrdist3 data
#'
#' @param dist_matrix A square distance matrix (from tcrdist3 $pw_beta or $pw_alpha)
#' @param meta_df A dataframe containing TCR metadata (e.g., V-genes, J-genes)
#' @param color_by The column name in meta_df to color the points by (default: "v_a_gene")
#' @param title Plot title (default: "TCRdist kPCA")
#' @param features Number of kernel PCA components to extract (default: 2)
#' @param enable.facet If TRUE, create faceted plots by multiple genes; if FALSE, single plot (default: TRUE)
#' @param cols.to.use Columns to facet by when enable.facet = TRUE
#'
#' @return A ggplot2 object
#'
generate_kpca_plot <- function(dist_matrix,
                               meta_df,
                               color_by = "v_a_gene",
                               title = "TCRdist kPCA",
                               features = 2,
                               enable.facet = TRUE,
                               cols.to.use = c("v_a_gene", "v_b_gene", "j_a_gene", "j_b_gene")) {
  # Validate inputs
  if (!color_by %in% colnames(meta_df)) {
    stop(paste0("Column '", color_by, "' not found in meta_df"))
  }
  
  # 1. Convert distances to a Gram Matrix (Similarity Matrix)
  # gram = 1 - (D / Dmax); high similarity -> values near 1, high distance -> values near 0
  d_max <- max(dist_matrix)
  gram_matrix <- 1 - (dist_matrix / d_max)
  
  # 2. Perform Kernel PCA using precomputed Gram matrix
  kpca_res <- kernlab::kpca(kernlab::as.kernelMatrix(as.matrix(gram_matrix)), features = features)
  if (!enable.facet)
  {
    # 3. Prepare plotting data
    plot_data <- data.frame(
      PC1 = rotated(kpca_res)[, 1],
      PC2 = rotated(kpca_res)[, 2],
      label = meta_df[[color_by]]
    )
    
    plot_data$label <- gsub("\\*.*", "", plot_data$label)
    
    # Calculate cluster centroids (median position of each gene/label)
    label_data <- plot_data %>%
      dplyr::group_by(label) %>%
      dplyr::summarise(PC1 = median(PC1), PC2 = median(PC2), .groups = "drop")
    
    # Create ggplot with direct labels (no legend for plotly compatibility)
    p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = label)) +
      geom_point(alpha = 0.6, size = 2.5) +
      geom_text(
        data = label_data,
        aes(label = label),
        size = 4,
        fontface = "bold",
        vjust = -1.5,
        check_overlap = TRUE,  # Respects overlaps (ggplotly-compatible alternative to ggrepel)
        show.legend = FALSE
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = title,
        subtitle = paste("Colored by:", color_by),
        x = "kPC1 (Kernel PC 1)",
        y = "kPC2 (Kernel PC 2)"
      ) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )
    return(p)
  }
  else
  {
    # 3. Prepare plotting data with facet variable
    plot_data <- cbind(data.frame(PC1 = rotated(kpca_res)[, 1], PC2 = rotated(kpca_res)[, 2]), meta_df[, cols.to.use]) %>% 
      pivot_longer(
        cols = all_of(cols.to.use),
        names_to = "facet_var",
        values_to = "label"
      )
    plot_data$label <- gsub("\\*.*", "", plot_data$label)
    
    # Calculate cluster centroids per facet for label placement
    label_data <- plot_data %>%
      dplyr::group_by(facet_var, label) %>%
      dplyr::summarise(PC1 = median(PC1), PC2 = median(PC2), .groups = "drop")
    
    # Create faceted ggplot with direct labels (no ggrepel for plotly compatibility)
    p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = label)) +
      geom_point(alpha = 0.6, size = 2.5) +
      geom_text(
        data = label_data,
        aes(label = label),
        size = 3.5,
        fontface = "bold",
        vjust = -1.2,
        check_overlap = TRUE,
        show.legend = FALSE
      ) +
      ggplot2::labs(
        title = title,
        subtitle = paste("Colored and Faceted by:", color_by),
        x = "kPC1 (Kernel PC 1)",
        y = "kPC2 (Kernel PC 2)"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      ) +
      ggplot2::facet_wrap(~ facet_var, scales = "free")
    return(p)
  }
  
}


#' Function to get distance matrix for a specific epitope
#' @param df Imported dataframe in tcrdist3 format (with columns like 'epitope', 'v_a_gene', 'v_b_gene', etc.)
#' @param epitope_name name of epitope to run tcrdist on (only one epitope at a time, if NULL runs on all data)
#' @param organism name of organism for which repertoire is being analyzed (e.g., "human", "mouse")
#' @param chains name of chain(s) to analyze (e.g., "alpha", "beta", or c("alpha","beta"))
#' @param dbfile path to tcrdist3 database file (e.g.,)
get_epitope_dist_matrix <- function(df,
                                    epitope_name = NULL,
                                    organism,
                                    chains,
                                    dbfile) {
  # Subset for specific epitope
  if (!is.null(epitope_name))
    subset_df <- df[df$epitope %in% epitope_name, ]
  else
    subset_df<-df
  
  pd <- import("pandas")
  message("Finished importing pandas")
  TCRrep <- import("tcrdist.repertoire")$TCRrep
  message("Finished importing TCRrep")
  # Convert R data frame to pandas DataFrame
  cell_df <- r_to_py(subset_df)
  # Get chains for analysis
  py <- import_main()
  
  if (length(chains) == 1)
    chains <- list(chains)
  
  # Initialize TCRrep object and compute distances
  tr <- TCRrep(
    cell_df = cell_df,
    organism = r_to_py(organism),
    chains = r_to_py(chains),
    db_file = r_to_py(dbfile)
  )
  cat("Done computing distances\n")
  # Return the tr object
  return(tr)
}


#' Subset an already-computed TCRrep result to one or more epitopes
#'
#' Distance computation (`tr$compute_distances()`) is by far the most
#' expensive step in a tcrdist3 analysis. Once it has been run once (e.g. in
#' the main "Run TCRdist3" analysis), narrowing down to a specific epitope
#' (and/or chain) is just an indexing operation on the clone table and the
#' pairwise distance matrices -- it does not require re-running tcrdist3.
#'
#' @param tr A TCRrep result (as returned by get_epitope_dist_matrix()) that
#'   distances have already been computed for.
#' @param epitope_name Optional vector of epitope(s) to subset to. If NULL,
#'   all clones are kept.
#' @return A list with a `clone_df` element and one element per pairwise
#'   distance matrix (named e.g. `pw_alpha`, `pw_cdr3_b_aa`, ...), each
#'   subset to the rows/columns matching `epitope_name`. Mirrors the parts
#'   of a TCRrep object that the rest of the app relies on, so it can be
#'   used as a drop-in replacement wherever `tr$clone_df` / `tr[[pw_name]]`
#'   are accessed.
subset_tr_by_epitope <- function(tr, epitope_name = NULL,chain=NULL) {
  clone_df <- as.data.frame(tr$clone_df)
  if (!is.null(epitope_name)) {
    idx <- which(clone_df$epitope %in% epitope_name)
  } else {
    idx <- seq_len(nrow(clone_df))
  }
  result <- list(clone_df = clone_df[idx, , drop = FALSE])
  pw_names <- names(tr)[grepl("^pw_", names(tr))]
  for (pw in pw_names) {
    mat <- tr[[pw]]
    result[[pw]] <- mat[idx, idx, drop = FALSE]
  }
  write.table(data.frame(Names=names(tr)),file="Columns.txt")
  chains.keep<-tr[["chains"]]
  if(!is.null(chain))
  {
    chains.keep<-list(chains.keep[chains.keep %in% chain])  
  }
  result[["chains"]]<-chains.keep
  return(result)
}


#' Render error message
#' @param message error message to be rendered
#' @param output shiny output object to render the error message in

#' Render a styled error message in the Shiny UI
#'
#' @param message Error message text; if empty string, clears the error display
#' @param output Shiny output object where the error will be rendered
#'
render_error_ui <- function(message, output) {
  if (identical(message, "")) {
    # Clear error display
    error_html <- ""
  } else {
    # Render error in a red bordered box
    error_html <- paste0(
      '<div style="color: red; font-weight: bold; padding: 10px; ',
      'border: 1px solid red; border-radius: 5px; background-color: #ffe6e6;">',
      message,
      '</div>'
    )
  }
  
  output$errorreport <- shiny::renderUI(HTML(error_html))
}


#' Remove temporary files created by the current Shiny session
#'
#' Cleans up HTML/plot files written to www/ directory during a session.
#'
#' @param session Shiny session object (provides unique session token)
#' @param pattern File name pattern to match (default: "*" for all session files)
#'
remove_www_tempfiles <- function(session, pattern = "*") {
  temp_files <- list.files(
    "www",
    pattern = paste0(session$token, ".", pattern),
    full.names = TRUE
  )
  
  if (length(temp_files) == 0) {
    return(invisible(NULL))
  }
  
  for (temp_file in temp_files) {
    if (file.exists(temp_file)) {
      tryCatch(
        file.remove(temp_file),
        error = function(e) message("Failed to remove ", temp_file, ": ", e$message)
      )
    }
  }
}



#' Compute repertoire diversity metrics per epitope
#'
#' @param tr An already-computed TCRrep result (e.g. from analysis_results()
#'   / get_epitope_dist_matrix()). Distances must already be computed.
#' @param chains Character vector of chain(s) used to compute `tr`, e.g.
#'   c("alpha","beta"). Their pairwise distance matrices (`pw_alpha`,
#'   `pw_beta`, ...) are summed to obtain the combined distance used for
#'   fuzzy diversity, mirroring how tcrdist3 combines chains elsewhere.
get_repertoir_diversity <- function(tr, chains)
{
  pd <- import("pandas")
  generalized_simpsons_entropy <- import("tcrdist.diversity")$generalized_simpsons_entropy
  fuzzy_diversity <- import("tcrdist.diversity")$fuzzy_diversity
  
  clone_df <- as.data.frame(tr$clone_df)
  epitopes <- unique(clone_df$epitope)
  message(paste0(epitopes, collapse = "\t"))
  
  # Combine whichever chain-level distance matrices were actually computed
  # (previously this was hard-coded to pw_alpha + pw_beta, which fails
  # whenever a single-chain or gamma/delta analysis was run).
  pw_names <- paste0("pw_", chains)
  pw_names <- pw_names[pw_names %in% names(tr)]
  if (length(pw_names) == 0)
    stop("None of the requested chain distance matrices (",
         paste0("pw_", chains, collapse = ", "),
         ") were found on the TCRrep result.")
  full_pwmat <- Reduce(`+`, lapply(pw_names, function(nm) tr[[nm]]))
  
  ###Compute the Simpson's diversity index (SDI) with 95% CI and the "effective number" of TCRs for a TCR repertoire.
  do.call(rbind, lapply(epitopes, function(epitope) {
    idx <- which(clone_df$epitope == epitope)
    sub_clone_df <- clone_df[idx, , drop = FALSE]
    sub_pwmat <- full_pwmat[idx, idx, drop = FALSE]
    
    # generalized_simpsons_entropy() only returns its named-column result
    # (with 'Z' and 'D') when given a pandas.Series -- a plain array/numeric
    # vector makes it return an unnamed (Z, lcl, ucl) tuple instead, which
    # is what the previous version passed in.
    counts_series <- pd$Series(sub_clone_df$count)
    simpson_res <- py_to_r(generalized_simpsons_entropy(counts_series))
    
    fuzzy_div <- fuzzy_diversity(r_to_py(sub_clone_df$count),
                                 sub_pwmat,
                                 order = 2,
                                 threshold = 75)
    data.frame(
      epitope = epitope,
      # Z: unbiased estimate of 1 - Simpson's index (probability two
      # randomly drawn clones differ)
      one_minus_simpson = simpson_res$Z[1],
      # D: "effective number" of clones, i.e. inverse Simpson
      inv_simpson = simpson_res$D[1],
      fuzzy_diversity = fuzzy_div
    )
  }))
}


generate_clustered_heatmap<-function(tr)
{
  
  annotation_df <- data.frame(
    Epitope = rep(c("Epitope_A", "Epitope_B", "Epitope_C"), each = 5)
  )
  rownames(annotation_df) <- rownames(mock_data)
  
  # Specify colors for your epitopes
  ann_colors <- list(
    Epitope = c(Epitope_A = "#E41A1C", Epitope_B = "#377EB8", Epitope_C = "#4DAF4A")
  )
  
  # ---------------------------------------------------------
  # 3. Plot the Clustered Heatmap
  # ---------------------------------------------------------
  pheatmap(mock_data,
           clustering_method = "ward.D2", # Ward's method minimizes variance within clusters
           annotation_col = annotation_df,
           annotation_row = annotation_df,
           annotation_colors = ann_colors,
           color = colorRampPalette(c("firebrick", "white", "navy"))(50),
           main = "TCRdist3 Distance Matrix Clustering",
           fontsize_row = 8,
           fontsize_col = 8)
  
}



#' Generate a network edge list from TCRrep distance matrices
#'
#' Constructs a network where edges connect TCR clones with distance <= threshold.
#' Distance is computed as the sum of alpha and beta chain pairwise distances.
#'
#' @param tr A TCRrep result (as returned by get_epitope_dist_matrix())
#' @param edge_threshold Maximum pairwise distance for an edge to exist in the network
#'
#' @return A data.frame with columns:
#'   - node_1, node_2: Clone indices (0-indexed)
#'   - dist: Combined alpha+beta distance between clones
#'   - v_b_gene_1, v_b_gene_2: V-gene assignments
#'   - cdr3_b_aa_1, cdr3_b_aa_2: CDR3 amino acid sequences
#'   - subject_1, subject_2: Subject/individual IDs
#'   - epitope_1, epitope_2: Epitope targets
#'   - K_neighbors: Number of neighbors (degree) for node_1
#'
#' @examples
#' # Create network with edge threshold at 50% of max distance
#' tr <- get_epitope_dist_matrix(data_df, organism = "human", chains = c("alpha", "beta"))
#' max_dist <- max(tr$pw_alpha + tr$pw_beta)
#' net <- generate_network(tr, edge_threshold = max_dist / 2)
#'
generate_network <- function(tr, edge_threshold, removeselfloops) {
  clone_df <- as.data.frame(tr$clone_df)
  # Sum alpha and beta pairwise distance matrices to get combined distance
  # (this allows measuring TCR similarity accounting for both chains)
  chains<-tr[["chains"]][[1]]
  print(chains)
  pw_alpha_beta<-NULL
  for(i in chains)
  {
    if(is.null(pw_alpha_beta))
      pw_alpha_beta<-tr[[paste0("pw_",i)]]
    else
      pw_alpha_beta <- pw_alpha_beta+tr[[paste0("pw_",i)]]
  }
    
  
  # Find all pairs of clones with combined distance <= threshold
  # (excluding self-edges via row(...) != col(...))
  edges <- which(
    pw_alpha_beta <= edge_threshold & row(pw_alpha_beta) != col(pw_alpha_beta),
    arr.ind = TRUE
  )
  
  # Extract 1-indexed row/col positions (R uses 1-indexing)
  node_1_idx <- edges[, "row"]
  node_2_idx <- edges[, "col"]
  
  # Calculate node degree (K_neighbors): number of neighbors per clone
  # (counting edges <= threshold, minus 1 for the self-loop)
  k_neighbors <- rowSums(pw_alpha_beta <= edge_threshold) - 1
  
  # Build edge list as a data.frame
  # Note: Indices are converted to 0-indexed (node - 1) for consistency with Python/Cytoscape
  df_net <- data.frame(
    node_1      = node_1_idx - 1,                      # 0-indexed
    node_2      = node_2_idx - 1,                      # 0-indexed
    dist        = pw_alpha_beta[edges],
    v_a_gene_1  = clone_df$v_a_gene[node_1_idx],
    v_a_gene_2  = clone_df$v_a_gene[node_2_idx],
    v_b_gene_1  = clone_df$v_b_gene[node_1_idx],
    v_b_gene_2  = clone_df$v_b_gene[node_2_idx],
    cdr3_b_aa_1 = clone_df$cdr3_b_aa[node_1_idx],
    cdr3_b_aa_2 = clone_df$cdr3_b_aa[node_2_idx],
    cdr3_a_aa_1 = clone_df$cdr3_a_aa[node_1_idx],
    cdr3_a_aa_2 = clone_df$cdr3_a_aa[node_2_idx],
    subject_1   = clone_df$subject[node_1_idx],
    subject_2   = clone_df$subject[node_2_idx],
    epitope_1   = clone_df$epitope[node_1_idx],
    epitope_2   = clone_df$epitope[node_2_idx],
    K_neighbors = k_neighbors[node_1_idx]
  ) %>% mutate(
    # Create temporary sorted columns to identify matching pairs
    min_col = pmin(node_1, node_2),
    max_col = pmax(node_1, node_2)
  ) %>% 
    # Keep only the first occurrence of each pair
    distinct(min_col, max_col, .keep_all = T) %>%
    # Remove the temporary columns
    select(-min_col, -max_col)
  
  #find self loops
  if(removeselfloops)
  {
    #match all self loops and remove them from the rest of the analysis
    if(length(chains)>1){
        df_net<-df_net[df_net$cdr3_b_aa_1!=df_net$cdr3_b_aa_2&df_net$cdr3_a_aa_1!=df_net$cdr3_a_aa_2,]
        
        df_net<- df_net %>% mutate(
          # Create temporary sorted columns to identify matching pairs
          min_col = pmin(cdr3_a_aa_1, cdr3_a_aa_2),
          max_col = pmax(cdr3_a_aa_2, cdr3_a_aa_1)
        ) %>% 
          # Keep only the first occurrence of each pair
          distinct(min_col, max_col, .keep_all = T) %>%
          # Remove the temporary columns
          select(-min_col, -max_col) %>%
          mutate(
            # Create temporary sorted columns to identify matching pairs
            min_col = pmin(cdr3_b_aa_1, cdr3_b_aa_2),
            max_col = pmax(cdr3_b_aa_2, cdr3_b_aa_1)
          ) %>% 
          # Keep only the first occurrence of each pair
          distinct(min_col, max_col, .keep_all = T) %>%
          # Remove the temporary columns
          select(-min_col, -max_col) %>%
          filter(!(cdr3_a_aa_1 == cdr3_a_aa_2 & cdr3_b_aa_1 == cdr3_b_aa_2)) %>%
          mutate(
            # Create temporary sorted columns to identify matching pairs
            min_col = pmin(cdr3_a_aa_1, cdr3_b_aa_1),
            max_col = pmax(cdr3_a_aa_2, cdr3_b_aa_2)
          ) %>% 
          # Keep only the first occurrence of each pair
          distinct(min_col, max_col, .keep_all = T) 
    }
    else
    {
      cdr<-paste0("cdr3_",substr(chains,1,1),"_aa")
      df_net<-df_net[df_net[[paste0(cdr,"_1")]]!=df_net[[paste0(cdr,"_2")]],]
      
      if(chains=="alpha")
        df_net<- df_net %>% mutate(
          # Create temporary sorted columns to identify matching pairs
          min_col = pmin(cdr3_a_aa_1, cdr3_a_aa_2),
          max_col = pmax(cdr3_a_aa_2, cdr3_a_aa_1)
        ) %>% 
        # Keep only the first occurrence of each pair
        distinct(min_col, max_col, .keep_all = T) %>%
        # Remove the temporary columns
        select(-min_col, -max_col) 
      else{
        df_net<- df_net %>% mutate(
          # Create temporary sorted columns to identify matching pairs
          min_col = pmin(cdr3_b_aa_1, cdr3_b_aa_2),
          max_col = pmax(cdr3_b_aa_2, cdr3_b_aa_1)
        ) %>% 
        # Keep only the first occurrence of each pair
        distinct(min_col, max_col, .keep_all = T) %>%
        # Remove the temporary columns
        select(-min_col, -max_col) 
      }
    }
  }
  write.table(df_net,file="network_edges.txt",sep="\t",row.names = F,quote=F)
  return(df_net)
}


#' R Function to generate motif enrichment
#' @param tcr_data A data.frame containing your TCR sequence data
#' @return A data.frame with the tcrdist3 results
generate_tcrdist_motifs <- function(tr,tcr_df,chain,radius=50L) {
  
  
  # Generate a sequence motif (logo) for a cluster -----------------
  # palmotif builds a probability matrix + SVG logo from an alignment
  # of CDR3 sequences within a neighborhood/cluster.
  
  cluster_idx <- neighbors[[1]] + 1  # back to 1-indexed for R
  cluster_seqs <- tcr_df[[paste0("cdr3_",substr(chain,1,1),"_aa")]][cluster_idx]
  
  palmotif <- import("palmotif")
  
  # palmotif needs a DataFrame with sequences (and optionally counts/weights)
  motif_input <- data.frame(seq = cluster_seqs, stringsAsFactors = FALSE)
  py_motif_input <- r_to_py(motif_input)
  motif_result <- palmotif$compute_pal_motif(
    seqs = py_motif_input$seq$tolist(),
    centroid = cluster_seqs[1]   # reference/anchor sequence for alignment
  )
  
  # motif_result typically returns (motif_df, stats) — a matrix of
  # position-specific amino acid probabilities plus an SVG logo string
  motif_matrix <- motif_result[[1]]
  svg_logo     <- motif_result[[2]]
  print(svg_logo)
  # ---- 6. Save the motif logo as an SVG file ------------------------------
  writeLines(svg_logo, con = "cluster1_motif.svg")
  
  # ---- 7. Inspect the motif matrix in R -----------------------------------
  motif_df_r <- py_to_r(motif_matrix)
  
  View(svg_logo)
}

message("Global file loaded!")
