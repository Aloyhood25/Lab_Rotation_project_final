#Generic script to analyse sperm motility data ----------------
#project name:
#authors:
#start date:
#last updated:

#remove all objects from global R environment
rm(list = ls())

# 1. Load required packages----------------------------
library(tidyverse)   # Data manipulation and visualization, version 2.0.0
library(broom)       # Tidy model outputs,version 1.0.2
library(swirl)       # Enhanced console output, version 2.4.5
library(lme4)        # Linear mixed-effects models, version 1.1-31
library(lmerTest)    # Mixed models with p-values, version 3.1-3
library(DHARMa)      # Model diagnostics, version 0.4.6
library(ggfortify)   # autoplot for models, version 0.4.14
library(emmeans)     # Estimated marginal means, version 1.8.5-1

# 2. Load and inspect my_data----------------------------
my_data <- read.csv('01 Raw data/seminal_fluid vs male age interaction.csv')
#my_data <- read.csv('01 Raw data/sperm_age vs SF_age.csv')

# Display variable names
cat("Variables detected:\n", names(my_data), sep = "\n")

# Inspect data
str(my_data)
head(my_data)

# Identify response variables
response_candidates <- grep("sperm.*motil", names(my_data), ignore.case = TRUE, value = TRUE)
if (length(response_candidates) == 0) stop("No sperm motility variable found.")
response_var <- response_candidates[1]
cat("Detected response variable:\n", response_var, sep = "\n")

# Identify predictor variables
predictors <- setdiff(names(my_data), response_var)
cat("Detected predictor variables:\n", predictors, sep = "\n")

# Automatically detect ID column
id_col <- grep("id", names(my_data), ignore.case = TRUE, value = TRUE)
if (length(id_col) > 0) {
  id_col <- id_col[1]
  my_data[[id_col]] <- as.factor(my_data[[id_col]])
  cat("Detected ID column:\n", id_col, sep = "\n")
} else {
  warning("No ID column detected — random effects will not be used.")
  id_col <- NULL
}

# 3. Keyword-based variable detection -----------
keyword_list <- list(
  treatment     = c("treatment_combination", "medium", "condition", "treatment"),
  age           = c("male_age", "sperm_age", "age"),
  species       = c("species"),
  temperature   = c("temp", "temperature"),
  time          = c("time_factor", "time"),
  concentration = c("conc", "concentration", "density"),
  id            = c("id", "male_ID", "subject")
)

# Function to search for variable name based on keyword list
detect_var <- function(keywords, data_names) {
  unique(unlist(lapply(keywords, function(k) 
    grep(k, data_names, ignore.case = TRUE, value = TRUE))))
}

# Detect variables
detected_vars <- lapply(keyword_list, detect_var, data_names = names(my_data))

# Helper to print detected variables
print_var <- function(label, var) {
  if (length(var) == 0) cat(label, ": no variable found\n") 
  else cat(label, ":", paste(var, collapse = ", "), "\n")
}

cat("Keyword-based variable detection:\n")
for (n in names(detected_vars)) {
  label <- paste0(toupper(substr(n, 1, 1)), substr(n, 2, nchar(n)), 
                  " variable(s)")
  print_var(label, detected_vars[[n]])
}

# Apply detection for each variable type
treatment_var <- detect_var(keyword_list$treatment, names(my_data))
age_var <- detect_var(keyword_list$age, names(my_data))
species_var <- detect_var(keyword_list$species, names(my_data))
temp_var <- detect_var(keyword_list$temperature, names(my_data))
time_var <- detect_var(keyword_list$time, names(my_data))
conc_var <- detect_var(keyword_list$concentration, names(my_data))
id_var <- detect_var(keyword_list$id, names(my_data))

# 4. Combine detected variables-----------
categorical_vars <- unique(c(treatment_var, species_var, age_var, time_var,
                             id_var))
numeric_vars <- unique(c(temp_var, conc_var))


# Fallback if nothing detected by keywords    
if (length(categorical_vars) == 0) {
  categorical_vars <- predictors[sapply(my_data[predictors],
                                        function(x) is.character(x) || is.factor(x))]
}
if (length(numeric_vars) == 0) {
  numeric_vars <- predictors[sapply(my_data[predictors], is.numeric)]
}

# Convert categorical variables to factors
my_data[categorical_vars] <- lapply(my_data[categorical_vars], as.factor)

cat("Final variable classification:\n", 
    "Categorical predictors:\n", categorical_vars, "\n",
    "Numeric predictors:\n", numeric_vars, "\n")


# 5. Create results directory----------------------------
dir.create("04 R plots", showWarnings = FALSE)


# 6. Visualization with ggplot2----------------------------
cat("Generating exploratory plots...\n")

# Only use the fixed factors for finding the correct plot
fixed.predictors <- c(treatment_var, age_var, time_var)

# Display the appropriate plot based on number of predictors detected
if (length(fixed.predictors) == 1) {
    ggplot(my_data, aes_string(x = time_var, y = response_var)) +
    geom_point(aes(group = treatment_var), alpha = 0.6, size = 2) +
    theme_classic() +
    labs(x = "Time (min)", y = "Sperm motility", color = "Treatments"
    )+
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size=16, face="bold"),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          strip.text = element_text(size=12, face="bold"),
          strip.background = element_rect(fill="lightgrey"))
} else if (length(fixed.predictors) == 2) {
    ggplot(my_data, aes_string(x = time_var, y = response_var, 
                                  color = treatment_var)) +
    geom_point(aes(group = treatment_var), alpha = 0.6, size = 2) +
    theme_classic() +
    labs(x = "Time (min)", y = "Sperm motility", color = "Treatments"
    )+
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size=16, face="bold"),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          strip.text = element_text(size=12, face="bold"),
          strip.background = element_rect(fill="lightgrey"))
} else if (length(fixed.predictors) == 3) {
    ggplot(my_data, aes_string(x = time_var, y = response_var, 
                                  color = treatment_var)) +
    geom_point(aes(group = treatment_var), alpha = 0.6, size = 2) +
    theme_classic() +
    facet_wrap(as.formula(paste("~", age_var[1])), nrow = 1) +
    labs(x = "Time (min)", y = "Sperm motility", color = "Treatments"
    )+
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size=16, face="bold"),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          strip.text = element_text(size=12, face="bold"),
          strip.background = element_rect(fill="lightgrey"))
  
} else {
  cat("More than three predictors detected. Think of how you want to plot this.\n")
}

# Save the exploratory plot
ggsave("04 R plots/exploratory_plot.png", width = 10, height = 6)

#create group means using detected predictor variables
# Identify which group of variables are present
group_vars <- c(
  if (length(age_var) >= 1) age_var[1],
  if (length(treatment_var) >= 1) treatment_var[1],
  if (length(time_var) >= 1) time_var[1]
)

# Stop if none exist
if (length(group_vars) == 0) {
  stop("No age, treatment, or time variable detected. Cannot proceed with 
       group means calculation.")
}

# Group and calculate means
general_mean_SM <- my_data %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    mean_SM = mean(.data[[response_var]], na.rm = TRUE),
    se_SM = sd(.data[[response_var]], na.rm = TRUE) / sqrt(n()),
    max_SM = max(.data[[response_var]], na.rm = TRUE),
    min_SM = min(.data[[response_var]], na.rm = TRUE)
  )

glimpse(general_mean_SM)

#order levels of time variable
levels(general_mean_SM[[time_var[1]]]) <- c("0", "2", "4", "6")

#Plot sperm motility over time variable, age variable, and treatment variables 
#use a conditional statement to plot depending on the number of predictor 
#variables detected

ggplot() +
  geom_point(data = my_data, aes_string(x = time_var[1], 
                                        y = response_var,
           color = treatment_var[1]), col = "grey", size = 0.5, alpha = 0.8) +
  geom_point(data = general_mean_SM, aes_string(x = time_var[1], y = "mean_SM",
           color = treatment_var[1]), size = 3) +
  facet_wrap(as.formula(paste("~", age_var[1])),nrow=1)+
  geom_line(data = general_mean_SM, aes_string(x = time_var[1], y = "mean_SM",
           group = treatment_var[1], color = treatment_var[1]),
           linewidth = 0.5,
           alpha = 0.5) +
  theme_classic() +
  geom_errorbar(data = general_mean_SM,
                aes_string(x = time_var[1], y = "mean_SM",
                           ymin = "mean_SM - se_SM",
                           ymax = "mean_SM + se_SM",
                           color = treatment_var[1]), width = 0.1) +
  labs(x = "Time (min)",
       y = "Sperm motility",
       color = "Treatment")+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size=16, face="bold"),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12, face="bold"),
        strip.background = element_rect(fill="lightgrey"))

# Save the sperm motility over time plot
ggsave("04 R plots/sperm_motility_over_time.png", width = 10, height = 6)

#Calculate overall mean sperm motility across all groups
overall_mean_SM <- my_data %>%
  group_by(across(all_of(id_var[1])),across(all_of(group_vars[1:2]))) %>%
  summarise(
    overall_mean_SM = mean(.data[[response_var]], na.rm = TRUE)
  )
group_vars
overall_mean_SM

#Plot overall mean sperm motility
ggplot() +
  geom_boxplot(data = overall_mean_SM,
               aes_string(x = age_var[1],
                          y = "overall_mean_SM",
                          fill = treatment_var[1]),
               position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(data = overall_mean_SM, aes_string(x = age_var[1], 
                                                y = "overall_mean_SM",
                                                fill = treatment_var[1]),
             position = position_dodge(width = 1),
             size = 2, alpha = 0.6) +
  theme_classic() +
  labs(x = "Age in days",
       y = "Mean individual sperm motility over time",
       fill = "Treatment")+
  theme(legend.position = "bottom",
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))

# Save the mean sperm motility plot
ggsave("04 R plots/mean_sperm_motility.png", width = 10, height = 6)

#calculate max sperm motility per individual across group variables 
    overall_max <- my_data %>%
      group_by(across(all_of(id_var[1])), across(all_of(group_vars[1:2]))) %>%
      summarise(
        max_SM = max(.data[[response_var]], na.rm = TRUE),
        se_SM = sd(.data[[response_var]], na.rm = TRUE) / sqrt(n()),
      )
  
overall_max


#Plot overall maximum sperm motility per individual across treatments 
ggplot() +
  geom_point(data = overall_max, aes(x = .data[[treatment_var[1]]], y = max_SM,
      color = .data[[age_var[1]]]),size = 2) +
  geom_errorbar(data = overall_max,aes(x = .data[[treatment_var[1]]],
      ymin = max_SM - se_SM,
      ymax = max_SM + se_SM,
      color = .data[[age_var[1]]]), width = 0.12, linewidth = 0.7) +
  geom_line(
    data = overall_max, aes(x = .data[[treatment_var[1]]], y = max_SM,
      group = .data[[id_var[1]]],
      color = .data[[age_var[1]]]), linewidth = 1.1, alpha = 0.7) +
  facet_wrap(as.formula(paste("~", age_var[1])), nrow=1) +
  labs(x = "Treatment",
    y = "Maximum individual sperm motility") +
  guides(color = "none") +
  theme_classic() +
  theme(
    legend.position = "right",,
    axis.title = element_text(size=14),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(size=12,angle=45, hjust=1),
    strip.text = element_text(size=12, face="bold"),
    strip.background = element_rect(fill="lightgrey")
  )

# Save the maximum sperm motility plot
ggsave("04 R plots/max_sperm_motility.png", width = 10, height = 6)

#calculate min sperm motility per individual across group variables 
overall_min <- my_data %>%
  group_by(across(all_of(id_var[1])), across(all_of(group_vars[1:2]))) %>%
  summarise(
    min_SM = min(.data[[response_var]], na.rm = TRUE),
    se_SM = sd(.data[[response_var]], na.rm = TRUE) / sqrt(n()),
  )
overall_min

#Plot overall minimum sperm motility per individual across treatments 
ggplot() +
  geom_point(data = overall_min, aes(x = .data[[treatment_var[1]]], y = min_SM,
                                         color = .data[[age_var[1]]]),size = 2) +
  geom_errorbar(data = overall_min,aes(x = .data[[treatment_var[1]]],
                                           ymin = min_SM - se_SM,
                                           ymax = min_SM + se_SM,
                                           color = .data[[age_var[1]]]), 
                width = 0.12, linewidth = 0.7) +
  geom_line(data = overall_min, aes(x = .data[[treatment_var[1]]], y = min_SM,
                                group = .data[[id_var[1]]],
                                color = .data[[age_var[1]]]), linewidth = 1, 
            alpha = 0.7) +
  facet_wrap(as.formula(paste("~", age_var[1])), nrow=1) +
  labs(x = "Treatment",
       y = "Minimum individual sperm motility") +
  guides(color = "none") +
  theme_classic() +
  theme(
    legend.position = "right",,
    axis.title = element_text(size=14),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(size=12,angle=45, hjust=1),
    strip.text = element_text(size=12, face="bold"),
    strip.background = element_rect(fill="lightgrey")
  )

# Save the minimum sperm motility plot
ggsave("04 R plots/min_sperm_motility.png", width = 10, height = 6)


cat("Running statistical model...\n")

# 7. Model Fitting----------------------------
#first run a linear model without random effects
linear_model <- lm(as.formula(paste('mean_SM', "~", paste(group_vars,
                                                          collapse = " + "))), 
                   data = general_mean_SM)
summary(linear_model)
anova(linear_model)
autoplot(linear_model, smooth.color = NA)


#determine the model type based on the data variables 
#run a mixed effects model if ID column is detected

#create model formula based on the number of predictors detected
if (length(predictors) == 1) {
  fixed_formula <- group_vars[1]
} else {
  fixed_formula <- paste(group_vars[1:min(3, length(group_vars))], 
                         collapse = " * ")
}

#presence or absence of random effects based on ID column detection
if (!is.null(id_col)) {
  formula <- as.formula(paste(response_var, "~", fixed_formula, "+ 
                              (1|", id_col, ")"))
} else {
  formula <- as.formula(paste(response_var, "~", fixed_formula))
}

#Run linear mixed effects model or linear model based on ID column detection
if (!is.null(id_col)) {
  model <- lmer(formula, data = my_data, REML = FALSE)
  model_type <- "Linear Mixed-Effects Model"
} else {
  model <- lm(as.formula(paste(response_var, "~", fixed_formula)), 
              data = my_data)
  model_type <- "Linear Model (no random effects)"
}
group_vars
head(my_data)
cat("Model type:", model_type, "\n")
print(summary(model))
vcov(model)

#8. Find the emmeans ----------------------------
library(emmeans)
emm <- emmeans(model, specs = as.formula(paste("~", paste(group_vars, 
                                                          collapse = " * "))))
emm
emm_plot <- as.data.frame(emm)
emm_plot
#pairwise comparisons within each time of measurement
simp <- pairs(emm, simple = treatment_var)
simp

#emmeans plot 
levels(emm_plot[[time_var[1]]]) <- c("0", "2", "4","6")
ggplot()+
  geom_line(data=emm_plot, aes_string(x = time_var[1], y = "emmean", 
                               group = treatment_var[1], 
                               color = treatment_var[1]),linewidth=1, 
            position=position_dodge(0.4))+
  geom_point(data=emm_plot, aes_string(x = time_var[1], y = "emmean",
                                group = treatment_var[1], 
                                color =  treatment_var[1]),size=3, 
             position=position_dodge(0.4))+
  geom_errorbar(data=emm_plot, aes_string(x = time_var[1], y = "emmean", 
                                   group = treatment_var[1], 
                                   color = treatment_var[1], 
                                   ymin  =  "emmean-SE",
                                   ymax  =   "emmean+SE"), width =  0.1, 
                linewidth  =  0.5, position=position_dodge(0.4))+
  geom_point(data=my_data, aes_string(x = time_var[1], y = response_var, 
                               group = treatment_var[1], 
                               color = treatment_var[1]),alpha=0.1,
             position=position_dodge(0.4))+
  facet_wrap(as.formula(paste("~", age_var[1])), nrow=1)+
  labs(x = "Time in minutes", 
       y = expression(paste("Sperm motility (temporal noise ",sigma,")")),
       color = "Treatment")+
  theme_classic()+
  theme(legend.position="bottom",
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12, face="bold"),
        strip.background = element_rect(fill="lightgrey"))


# Check model diagnostics using DHARMa
simulation_output <- simulateResiduals(fittedModel = model, n = 1000)

# Plot residuals
plot(simulation_output)
testResiduals(simulation_output)
testCategorical(simulation_output, catPred = my_data[[treatment_var[1]]])

if(testOutliers(simulation_output)$p.value < 0.05){
  cat("Warning: Outliers detected in the model residuals (p < 0.05).\n")
} else {
  cat("No significant outliers detected in the model residuals.\n")
}
if(testUniformity(simulation_output)$p.value < 0.05){
  cat("Warning: Residuals are not uniformly distributed (p < 0.05).\n")
} else {
  cat("Residuals are uniformly distributed.\n")
}
if(testDispersion(simulation_output)$p.value < 0.05){
  cat("Warning: Dispersion issues detected in the model (p < 0.05).\n")
} else {
  cat("No significant dispersion issues detected in the model.\n")
}
if(testQuantiles(simulation_output)$p.value < 0.05){
  cat("Warning: Quantile issues detected in the model (p < 0.05).\n")
} else {
  cat("No significant quantile issues detected in the model.\n")
}
if(testZeroInflation(simulation_output)$p.value < 0.05){
  cat("Warning: Zero-inflation issues detected in the model (p < 0.05).\n")
} else {
  cat("No significant zero-inflation issues detected in the model.\n")
}



run_dharma_tests <- function(simulation_output) {
  
  # List of tests to apply, with labels
  tests <- list(
    "Outliers"      = testOutliers,
    "Uniformity"    = testUniformity,
    "Dispersion"    = testDispersion,
    "Quantiles"     = testQuantiles,
    "Zero-inflation"= testZeroInflation
  )
  
  # Loop through each test
  for (test_name in names(tests)) {
    test_fun <- tests[[test_name]]
    
    # Run test
    result <- test_fun(simulation_output)
    
    # Extract p-value safely
    p <- result$p.value
    
    # Print results
    if (p < 0.05) {
      cat(sprintf("Warning: %s issues detected (p = %.4f)\n", test_name, p))
    } else {
      cat(sprintf("No significant %s issues detected (p = %.4f)\n", test_name, p))
    }
  }
}

run_dharma_tests(simulation_output)

#if tests indicate issues (p < 0.05), consider model refinement or consult 
#a statistician.


# 9. Model Summary and Visualization of Fixed Effects---------------------------
summary_df <- tidy(model, effects = "fixed")
summary_df
ggplot(summary_df, aes(x = term, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2) +
  labs(title = "Fixed Effects Estimates", x = "Terms", y = "Estimates") +
  theme_minimal()


# Save model summary to CSV
write.csv(summary_df, file = "04 R plots2/mixed_effects_model_summary.csv", 
          row.names = FALSE)

