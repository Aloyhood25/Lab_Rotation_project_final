#create a working directory with folders for raw data, clean data, and plots
#check if project directory already exists, otherwise create one
init_project <- function(path){
  #path <- '/Desktop/sperm_motility_analysis_project'  # Change this to your desired path
  subfolders <- c("01 Raw data", "02 Clean data", "03 R scripts",
                  "04 R plots", "05 Results")
  
  if (!dir.exists(path)) {dir.create(path, recursive = TRUE)} else {stop('project folder already exists, no files overwritten')}
  
  for(folder in subfolders){
    dir.create(file.path(path, folder), showWarnings = FALSE, recursive = TRUE)
  }
  
  setwd(path)
  message("📁 Project initialized at: ", path)
}

init_project("~/Desktop/SpermMotility_Project")
