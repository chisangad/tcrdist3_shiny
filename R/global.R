source("checkConda.R")
if(!require("reticulate"))
{
  install.packages("reticulate",dependencies = TRUE)
  library("reticulate")
}

if(!require("shiny"))
{
  install.packages("shiny")
  library("shiny")
}

if(!require(shinyjs))
{
  install.packages("shinyjs")
  library(shinyjs)
}

if(!require("pheatmap"))
{
  install.packages("pheatmap")
  library("pheatmap")
}

if(!require("DT"))
{
  install.packages("DT")
  library("DT")
}
if(!require("networkD3"))
{
  install.packages("networkD3")
  library(networkD3)
}
if(!require('openxlsx'))
{
  install.packages("openxlsx")
  library("openxlsx")
} 

if(!require('dplyr'))
{
  install.packages("dplyr")
  library("dplyr")
} 
if(!require(shinyFeedback))
{
  install.packages("shinyFeedback")
  library(shinyFeedback)
}
  
# py_run_string("
# def process_list(my_list):
#         return my_list
# ")

env_name <- "tcrdist3_env"

# --- Python Environment Setup ---
tryCatch({
  if(!condaenv_exists(env_name))
  {
    message(paste0("Conda environment ",env_name," not found."))
    message(paste0("Creating conda environment with name ",env_name))
    pkgs <- c("python=3.8")
    conda_create(envname = env_name, packages = pkgs,
                 channel= c("conda-forge","bioconda","defaults"))
    use_condaenv(env_name, required = TRUE)
    #conda_install(env_name,packages = c("make","autoconf", "automake", "libtool","glib-tools"))
    conda_install(env_name,packages = "tcrdist3",pip = T)
  }
  use_condaenv(env_name)
}, error = function(e) {
    message(paste0("Error initialising conda environment - ",e))
    #conda_remove(env_name)
    stop("Quiting")
})

py_run_string("
def process_r_vector(input_vector):
    print('Received object type in Python:', type(input_vector))
    print('Received Python list:', input_vector)
    return input_vector
")