# R script to create a project directory structure for 
# a reproducible data analysis

#Before you start working in the lab set up a project directory structure
# with the following folders:
# 01_raw_data
# 02_clean_data
# 03_R_scripts
# 04_R_plots
# 05_R_reports
# 06_manuscript
# 07_references

# 1. define where you want to create the project directory structure, for this
#you need to create a project folder and specify its path here:
existing_path <- '/Users/oliverotti/Documents/various stuff'

# 2. now run the following code to create the subfolders
init_project <- function(path){
  path <- existing_path
  # Change this to your desired path
  subfolders <- c("01_raw_data", "02_clean_data", "03_R_scripts",
                  "04_R_plots", "05_R_reports", "06_manuscript",
                  "07_references")
  
  if(!dir.exists(path)) dir.create(path, recursive = TRUE)
  
  for(folder in subfolders){
    dir.create(file.path(path, folder), showWarnings = FALSE, recursive = TRUE)
  }
  
  setwd(path)
  message("📁 Project initialized at: ", path)
}

# Change the following path to your desired location
init_project()#if you have not set the path in the function, it will use 
#the existing_path object defined above, otherwise you can pass a different path
#as an argument in the init_project() function

# Now you have a reproducible project directory structure set up!