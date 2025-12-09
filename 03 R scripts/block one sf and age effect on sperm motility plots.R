#age effect on sperm motility in Drosophila melanogaster block one
#Oliver Otti
#10.05.22

#BEFORE WE START SOME THINGS TO REMEMBER
# always comment code
# indent and spacing, use control+i to tell R to indent your code
# keep yourself to 80 characters per line, if possible: split comments and 
# function calls if necessary
# never save the environment!

# Load packages ----------------------------------------------------------

#clear R's brain
rm(list=ls())

###############################################################
# Generic Automated Mixed-Effects Analysis for Sperm Motility
# Description:
# - Automatically detects response & predictors
# - Treats ID as a factor for plotting & summary
# - Fits mixed-effect model with random intercept for ID
# - Uses ggplot2 for all plots
###############################################################

#----------------------------#
# 0. Load Required Packages
#----------------------------#
packages <- c("tidyverse", "broom.mixed", "lme4")
invisible(lapply(packages, function(p) {
  if (!require(p, character.only = TRUE)) install.packages(p, dependencies = TRUE)
}))
library(tidyverse)
library(lme4)
library(broom.mixed)

#----------------------------#
# 1. USER INPUT
#----------------------------#
file_path <- "01 Raw Data /sperm_age vs SF_age.csv"  # Change this line if needed

#----------------------------#
# 2. Load and Clean Data
#----------------------------#
my_data <- read.csv(file_path)
names(my_data) <- make.names(names(my_data))  # Clean column names
cat("Dataset loaded successfully: ", file_path, "\n")

cat("Variables detected:\n")
print(names(my_data))

#----------------------------#
# 3. Identify Response Variable
#----------------------------#
response_candidates <- grep("sperm.*motil", names(my_data), ignore.case = TRUE, value = TRUE)
if (length(response_candidates) == 0) stop("No sperm motility variable found.")
response_var <- response_candidates[1]
cat("Detected response variable:", response_var, "\n")

#----------------------------#
# 4. Identify Predictors and ID
#----------------------------#
predictors <- setdiff(names(my_data), response_var)

# Automatically detect ID column (e.g., male_ID, ID, bird_ID)
id_col <- grep("id", names(my_data), ignore.case = TRUE, value = TRUE)
if (length(id_col) == 0) {
  warning("No ID column detected. Mixed model will run without random effect.")
  id_col <- NULL
} else {
  id_col <- id_col[1]
  cat("Detected ID column:", id_col, "\n")
  my_data[[id_col]] <- as.factor(my_data[[id_col]])  # Treat ID as factor
}

# Exclude ID column from predictors for fixed effects
predictors <- predictors[predictors != id_col]

# Detect variable types
categorical_vars <- predictors[sapply(my_data[predictors], function(x) is.character(x) || is.factor(x))]
numeric_vars <- predictors[sapply(my_data[predictors], is.numeric)]

# Convert character columns to factors
my_data[categorical_vars] <- lapply(my_data[categorical_vars], as.factor)

cat("Detected categorical predictors:", categorical_vars, "\n")
cat("Detected numeric predictors:", numeric_vars, "\n")

#----------------------------#
# 5. Create Results Folder
#----------------------------#
dir.create("results", showWarnings = FALSE)

#----------------------------#
# 6. Visualizations (ggplot2)
#----------------------------#
cat("Generating exploratory plots...\n")

## ---- 6.1 Scatter plots for numeric predictors ---- ##
for (v in numeric_vars) {
  p <- ggplot(my_data, aes_string(x = v, y = response_var, color = id_col)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(aes_string(group = id_col), alpha = 0.4) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Sperm motility vs", v, "(colored by ID)"),
      x = v, y = response_var, color = "ID"
    )
  ggsave(paste0("results/scatter_", v, "_byID.png"), p, width = 7, height = 5)
}

## ---- 6.2 Boxplots for categorical predictors ---- ##
for (v in categorical_vars) {
  p <- ggplot(my_data, aes_string(x = v, y = response_var, fill = id_col)) +
    geom_boxplot(alpha = 0.7, outlier.color = "black", position = position_dodge2()) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.5, size = 2) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Sperm motility by", v, "(grouped by ID)"),
      x = v, y = response_var, fill = "ID"
    )
  ggsave(paste0("results/boxplot_", v, "_byID.png"), p, width = 7, height = 5)
}

## ---- 6.3 Per-ID Means ± SE ---- ##
if (!is.null(id_col)) {
  cat("Calculating per-ID means...\n")
  
  id_means <- my_data %>%
    group_by(across(all_of(id_col))) %>%
    summarise(
      mean_response = mean(.data[[response_var]], na.rm = TRUE),
      se = sd(.data[[response_var]], na.rm = TRUE) / sqrt(n())
    )
  
  # Plot mean ± SE per ID
  p <- ggplot(id_means, aes_string(x = id_col, y = "mean_response")) +
    geom_col(fill = "skyblue", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_response - se, ymax = mean_response + se), width = 0.3) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Mean ± SE of", response_var, "per ID"),
      x = "ID", y = paste("Mean", response_var)
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("results/mean_per_ID.png", p, width = 8, height = 5)
}

## ---- 6.4 Group Means by Treatment & ID ---- ##
if (length(categorical_vars) > 0) {
  for (cat_var in categorical_vars) {
    group_means <- my_data %>%
      group_by(across(all_of(c(cat_var, id_col)))) %>%
      summarise(
        mean_response = mean(.data[[response_var]], na.rm = TRUE),
        se = sd(.data[[response_var]], na.rm = TRUE) / sqrt(n())
      )
    
    p <- ggplot(my_data, aes_string(x = cat_var, y = response_var, color = id_col, group = id_col)) +
      geom_point(alpha = 0.5, size = 2) +
      geom_line(alpha = 0.4) +
      geom_point(data = group_means, aes_string(y = "mean_response"), size = 3, shape = 18) +
      geom_errorbar(data = group_means, aes(ymin = mean_response - se, ymax = mean_response + se), width = 0.15) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste("Group means ± SE of", response_var, "by", cat_var, "and ID"),
        x = cat_var, y = response_var, color = "ID"
      )
    ggsave(paste0("results/means_", cat_var, "_byID.png"), p, width = 8, height = 5)
  }
}

#----------------------------#
# 7. Mixed-Effect Model
#----------------------------#
cat("Running mixed-effects model...\n")

if (length(predictors) == 1) {
  fixed_formula <- predictors[1]
} else {
  fixed_formula <- paste(predictors[1:min(2, length(predictors))], collapse = " * ")
}

if (!is.null(id_col)) {
  formula <- as.formula(paste(response_var, "~", fixed_formula, "+ (1|", id_col, ")"))
} else {
  formula <- as.formula(paste(response_var, "~", fixed_formula))
}

if (!is.null(id_col)) {
  model <- lmer(formula, data = my_data, REML = FALSE)
  model_type <- "Linear Mixed-Effects Model"
} else {
  model <- lm(as.formula(paste(response_var, "~", fixed_formula)), data = my_data)
  model_type <- "Linear Model (no random effects)"
}

cat("Model type:", model_type, "\n")
print(summary(model))

summary_df <- tidy(model, effects = "fixed")
write.csv(summary_df, "results/model_summary.csv", row.names = FALSE)

#----------------------------#
# 8. Wrap up
#----------------------------#
cat("Analysis complete.\nResults and plots saved in 'results/' folder.\n")

###############################################################
# End of Script
###############################################################