# Define a function to check and install Miniconda
check_and_install_miniconda <- function() {
  # Check if Miniconda is already installed
  if (nchar(Sys.which("conda"))>0) {
    message("Miniconda is already installed. ✅")
    return(invisible(TRUE))
  } else {
    message("Miniconda not found. Attempting to download and install... ⏳") 
    # Define installation paths and filenames
    if (Sys.info()["sysname"] == "Windows") {
      # Windows
      installer_url <- "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe"
      installer_name <- "Miniconda_installer.exe"
      miniconda_path <- file.path(Sys.getenv("LOCALAPPDATA"), "Continuum", "miniconda3")
    } else if (Sys.info()["sysname"] == "Darwin") {
      # macOS
      installer_url <- "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
      installer_name <- "Miniconda_installer.sh"
      miniconda_path <- file.path("~", "miniconda3")
    } else {
      # Linux
      installer_url <- "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
      installer_name <- "Miniconda_installer.sh"
      miniconda_path <- file.path("~", "miniconda3")
    }
    
    # Download the installer
    tryCatch({
      download.file(installer_url, installer_name, mode = "wb")
    }, error = function(e) {
      stop("Failed to download Miniconda installer. Please check your internet connection or URL. ❌")
    })
    
    # Install Miniconda
    if (Sys.info()["sysname"] == "Windows") {
      # For Windows, run the silent installer
      system2(installer_name, args = c("/S", paste("/D=", miniconda_path, sep = "")))
    } else {
      # For macOS and Linux, run the shell script with silent installation options
      system2("bash", args = c(installer_name, "-b", "-p", miniconda_path))
    }
    system(paste0(miniconda_path,"/bin/conda tos accept"),intern = T)
    # Clean up the installer file
    file.remove(installer_name)
    
    # Check if installation was successful by verifying the 'conda' command
    if (nzchar(Sys.which("conda"))) {
      message("Miniconda installed successfully! 🎉")
      return(invisible(TRUE))
    } else {
      message("Installation complete, but 'conda' command not found in PATH. You may need to restart your R session or add Miniconda to your PATH manually. ⚠️")
      return(invisible(FALSE))
    }
  }
}

# Run the check and install function
check_and_install_miniconda()
