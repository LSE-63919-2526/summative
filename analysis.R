#####################################################################################################
#
# Codebook and Associated Graphics
# MY456 | WT 2026 | Summative Assigment
#
#####################################################################################################

# Global options for qmd
knitr::opts_chunk$set(warning = FALSE, message = FALSE, linewidth = 60, echo = T)

# install.packages(c("survey", "tidyverse"))
library(survey)
library(tidyverse)

#####################################################################################################
# Part C - Preliminary Analysis
#####################################################################################################

## Loading Data
df <- read.csv("data/silicon_sample_20260406_173053.csv", stringsAsFactors = FALSE)

# Quick check
nrow(df)       # Should be 2,204
names(df)      # Check variable names are as expected

#####################################################################################################

## Clean and Recode Cogbot Responses

# Cogbot encodes responses as letters. Inspect unique values first:
table(df$response, useNA = "always")

# Recode to meaningful labels
# A / A. = Yes  |  B / B. = No  |  C = Not sure  |  D = (treat as NA)
df <- df %>%
  mutate(
    env_steps = case_when(
      response %in% c("A", "A.") ~ "Yes",
      response %in% c("B", "B.") ~ "No",
      response %in% c("C", "C.") ~ "Not sure",
      TRUE ~ NA_character_          # D and anything unexpected -> NA
    ),
    env_steps = factor(env_steps, levels = c("Yes", "No", "Not sure"))
  )

# Verify recode
table(df$env_steps, useNA = "always")

#####################################################################################################

## Complex Survey Design
# stratum and psu: Sampling structure
# anweight: Combines design, calibration, and analysis weights

design <- svydesign(
  ids      = ~psu,        # Primary sampling units (clusters)
  strata   = ~stratum,    # Stratification variable
  weights  = ~anweight,   # Combined analysis weight
  data     = df,
  nest     = TRUE         # PSU IDs are nested within strata
)

# Summary of the design object
summary(design)

#####################################################################################################

## Weighted Frequency Table

# Estimated population proportions with standard errors
props <- svymean(~env_steps, design, na.rm = TRUE)
print(props)

# Convert to a tidy data frame with confidence intervals
props_df <- as.data.frame(confint(props)) %>%
  rownames_to_column("response") %>%
  mutate(
    response   = str_remove(response, "env_steps"),   # Tidy up label
    proportion = coef(props),
    pct        = proportion * 100,
    ci_lower   = `2.5 %` * 100,
    ci_upper   = `97.5 %` * 100
  ) %>%
  select(response, pct, ci_lower, ci_upper)

print(props_df)

#####################################################################################################

## Weighted Counts
# Estimated number of UK adults in each category

counts <- svytotal(~env_steps, design, na.rm = TRUE)
print(counts)

#####################################################################################################

## Bar Chart of Weighted Proportions

ggplot(props_df, aes(x = response, y = pct, fill = response)) +
  geom_col(width = 0.5, colour = "white") +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct)),
    vjust = -0.8,
    size  = 4
  ) +
  scale_fill_manual(values = c("Yes" = "#2E86AB", "No" = "#E84855", "Not sure" = "#A8A8A8")) +
  scale_y_continuous(limits = c(0, 80), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "In the past month, have you taken any steps to\nreduce your personal impact on the environment?",
    subtitle = "Weighted estimates, UK adult population (n = 2,204)",
    x        = NULL,
    y        = "Weighted percentage (%)",
    caption  = "Error bars show 95% confidence intervals.\nWeights: anweight (design + calibration + analysis). Source: Cogbot silicon sample."
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("environment_initial_design_plot.png", width = 7, height = 5, dpi = 300)

#####################################################################################################

## Breakdown by Sex

# Weighted cross-tabulation: response by sex
sex_tab <- svyby(~env_steps, ~sex, design, svymean, na.rm = TRUE)
print(sex_tab)

# Clean up for plotting
sex_df <- sex_tab %>%
  as.data.frame() %>%
  select(sex, env_stepsYes, env_stepsNo, `env_stepsNot sure`) %>%
  pivot_longer(
    cols      = starts_with("env_steps"),
    names_to  = "response",
    values_to = "proportion"
  ) %>%
  mutate(
    response = str_remove(response, "env_steps"),
    pct      = proportion * 100
  )

ggplot(sex_df, aes(x = sex, y = pct, fill = response)) +
  geom_col(position = "dodge", width = 0.6, colour = "white") +
  scale_fill_manual(values = c("Yes" = "#2E86AB", "No" = "#E84855", "Not sure" = "#A8A8A8")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Environmental steps taken in past month, by sex",
    subtitle = "Weighted estimates, UK adult population",
    x        = NULL,
    y        = "Weighted percentage (%)",
    fill     = "Response",
    caption  = "Source: Cogbot silicon sample."
  ) +
  theme_minimal(base_size = 13)

ggsave("environment_by_sex.png", width = 7, height = 5, dpi = 300)

#####################################################################################################

## Breakdown by Age Group

# Create age bands from continuous age variable
df <- df %>%
  mutate(
    age_group = cut(
      agea,
      breaks = c(17, 29, 44, 59, 74, Inf),
      labels = c("18-29", "30-44", "45-59", "60-74", "75+")
    )
  )

# Redefine design with age_group added
design2 <- svydesign(
  ids     = ~psu,
  strata  = ~stratum,
  weights = ~anweight,
  data    = df,
  nest    = TRUE
)

age_tab <- svyby(~env_steps, ~age_group, design2, svymean, na.rm = TRUE)
print(age_tab)

#####################################################################################################

## Significance Testing

```{r}
# Chi-squared test: Is is the distribution of responses independent of sex?
sex_test <- svychisq(~env_steps + sex, design)
print(sex_test)

```

Text goes here about preliminary analysis

#####################################################################################################

# Helper Functions (following seminar 5 style)

# Proportions: SE, 95% CI, and optionally DEFF / DEFT / Neff
calc_prop_table <- function(design, formula, include_deff = TRUE) {
  prop  <- svymean(formula, design = design, na.rm = TRUE)
  ci    <- confint(prop)
  result <- data.frame(
    Proportion = coef(prop),
    SE         = SE(prop),
    CI_Lower   = ci[, 1],
    CI_Upper   = ci[, 2]
  )
  if (include_deff) {
    var_name <- all.vars(formula)
    n        <- sum(!is.na(design$variables[[var_name]]))
    p        <- coef(prop)
    se       <- SE(prop)
    srs_var  <- p * (1 - p) / n
    result$DEFF <- (se^2) / srs_var
    result$DEFT <- se / sqrt(srs_var)
    result$Neff <- n / result$DEFF
  }
  result
}

#####################################################################################################
# PART D - Testing Estimation Options
#####################################################################################################

## Levels of Complex Survey Design


# Design 1: SRS  (no weight, no strata, no PSU)
# Design 2: Weight only
# Design 3: Weight + stratification
# Design 4: Weight + stratification + PSU  <-- correct full design
 
 
# ---- Design 1: Simple random sample (naive) ----------------
 
design_srs <- svydesign(ids = ~1, data = df)
 
cat("\n===== DESIGN 1: SRS =====\n", "\n-- env_steps --\n")
print(calc_prop_table(design_srs, ~env_steps, include_deff = FALSE))
 
cat("\n-- citizen --\n")
print(calc_prop_table(design_srs, ~citizen, include_deff = FALSE))
 
 
# ---- Design 2: Weight only ---------------------------------
 
design_wt <- svydesign(ids = ~1, weights = ~anweight, data = df)
 
cat("\n===== DESIGN 2: Weight only =====\n", "\n-- env_steps --\n")
print(calc_prop_table(design_wt, ~env_steps))
 
cat("\n-- citizen --\n")
print(calc_prop_table(design_wt, ~citizen))
 
 
# ---- Design 3: Weight + stratification ---------------------
 
design_wt_str <- svydesign(
  ids     = ~1,
  weights = ~anweight,
  strata  = ~stratum,
  data    = df,
  nest    = TRUE
)
 
cat("\n===== DESIGN 3: Weight + Stratification =====\n", "\n-- env_steps --\n")
print(calc_prop_table(design_wt_str, ~env_steps))
 
cat("\n-- citizen --\n")
print(calc_prop_table(design_wt_str, ~citizen))
 
 
# ---- Design 4: Weight + stratification + PSU (full design) -
 
design_full <- svydesign(
  ids     = ~psu,
  weights = ~anweight,
  strata  = ~stratum,
  data    = df,
  nest    = TRUE
)
 
cat("\n===== DESIGN 4: Full design (Weight + Strata + PSU) =====\n", "\n-- env_steps --\n")
print(calc_prop_table(design_full, ~env_steps))
 
cat("\n-- citizen --\n")
print(calc_prop_table(design_full, ~citizen))

#####################################################################################################

## Clean Side-by-Side Comparison Tables

# Helper: extract one named category row across the four designs
make_comparison <- function(cat_label, d1, d2, d3, d4, formula) {
  extract <- function(d, deff) {
    res <- calc_prop_table(d, formula, include_deff = deff)
    res[grepl(cat_label, rownames(res), fixed = TRUE), , drop = FALSE]
  }
  bind_rows(
    mutate(extract(d1, FALSE), Design = "1. SRS"),
    mutate(extract(d2, TRUE),  Design = "2. Weight only"),
    mutate(extract(d3, TRUE),  Design = "3. Weight + Strata"),
    mutate(extract(d4, TRUE),  Design = "4. Weight + Strata + PSU")
  ) %>%
    select(Design, everything()) %>%
    mutate(across(where(is.numeric), ~round(.x, 4)))
}
 
cat("\n===== COMPARISON TABLE: env_steps — 'Yes' category =====\n")
print(make_comparison("Yes", design_srs, design_wt, design_wt_str, design_full, ~env_steps))
 
cat("\n===== COMPARISON TABLE: citizen — 'UK citizen' category =====\n")
print(make_comparison("UK citizen", design_srs, design_wt, design_wt_str, design_full, ~citizen))

##################################################################################################### 

## Weight Distribution Plot

ggplot(df, aes(x = anweight)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Distribution of analysis weights (anweight)",
    subtitle = "Cogbot silicon sample, n = 2,204",
    x        = "anweight",
    y        = "Frequency",
    caption  = "Right skew indicates some respondents up-weighted to represent under-sampled groups."
  )
 
ggsave("weight_distribution.png", width = 7, height = 4, dpi = 300)

#####################################################################################################

## Bar Charts: Full Design Estimates with 95% CIs
 
# env_steps
props_env <- calc_prop_table(design_full, ~env_steps) %>%
  rownames_to_column("category") %>%
  mutate(category = str_remove(category, "env_steps"))
 
ggplot(props_env, aes(x = category, y = Proportion * 100, fill = category)) +
  geom_col(width = 0.5, colour = "white") +
  geom_errorbar(aes(ymin = CI_Lower * 100, ymax = CI_Upper * 100), width = 0.15) +
  geom_text(aes(label = sprintf("%.1f%%", Proportion * 100)), vjust = -0.8, size = 4) +
  scale_fill_manual(values = c("Yes" = "#2E86AB", "No" = "#E84855", "Not sure" = "#A8A8A8")) +
  scale_y_continuous(limits = c(0, 80), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Steps taken to reduce environmental impact (past month)",
    subtitle = "Full complex design estimates, UK adult population (n = 2,204)",
    x = NULL, y = "Weighted proportion (%)",
    caption  = "Error bars = 95% CI. Design: anweight + stratum + PSU."
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
 
ggsave("env_steps_full_design.png", width = 7, height = 5, dpi = 300)

# citizen
props_cit <- calc_prop_table(design_full, ~citizen) %>%
  rownames_to_column("category") %>%
  mutate(category = str_remove(category, "citizen"))
 
ggplot(props_cit, aes(x = category, y = Proportion * 100, fill = category)) +
  geom_col(width = 0.4, colour = "white") +
  geom_errorbar(aes(ymin = CI_Lower * 100, ymax = CI_Upper * 100), width = 0.12) +
  geom_text(aes(label = sprintf("%.1f%%", Proportion * 100)), vjust = -0.8, size = 4) +
  scale_fill_manual(values = c("UK citizen" = "#3D9970", "Not UK citizen" = "#FF851B")) +
  scale_y_continuous(limits = c(0, 110), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "UK citizenship status",
    subtitle = "Full complex design estimates, UK adult population (n = 2,204)",
    x = NULL, y = "Weighted proportion (%)",
    caption  = "Error bars = 95% CI. Design: anweight + stratum + PSU."
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
 
ggsave("citizen_full_design.png", width = 6, height = 5, dpi = 300)

#####################################################################################################

# PARTS E & F — Key Quantities Printed for Commentary

res_env_full <- calc_prop_table(design_full, ~env_steps)
cat("\nFull design estimates for env_steps:\n")
print(round(res_env_full[, c("Proportion", "SE", "DEFF", "DEFT", "Neff")], 4))
 
res_cit_full <- calc_prop_table(design_full, ~citizen)
cat("\nFull design estimates for citizen:\n")
print(round(res_cit_full[, c("Proportion", "SE", "DEFF", "DEFT", "Neff")], 4))
 
res_env_wt <- calc_prop_table(design_wt, ~env_steps)
res_cit_wt <- calc_prop_table(design_wt, ~citizen)
cat("\nWeight-only Neff for env_steps:\n")
print(round(res_env_wt[, c("Proportion", "SE", "DEFF", "Neff")], 4))
cat("\nWeight-only Neff for citizen:\n")
print(round(res_cit_wt[, c("Proportion", "SE", "DEFF", "Neff")], 4))
