# Propensity score matching: pregnant vs non-pregnant MS cohorts
# Aim: 1:1 matching of early postpartum women with MS to non-pregnant women with MS

#### CONTENTS #####
# 1. WORKSPACE SETUP
# 2. PROPENSITY SCORE MATCHING

#----------------------------------------------------------------------------------------------------------------------------------
#### 1. WORKSPACE SETUP ####
rm(list = ls())
library(tidyverse)
library(MatchIt)

# Input: one row per woman, both cohorts combined
# preg_status = 1 (early postpartum) or 0 (non-pregnant)
# Covariates: age_at_mri, phenotype_at_mri, disease_duration_at_mri, edss_closest,
#             relapse_preceding_3mo, year, het_met_preceding_yr
combined_data <- read_csv("combined_cohorts_for_psm.csv")

# Set variable types
str(combined_data)

combined_data <- combined_data %>%
  mutate(
    phenotype_at_mri      = factor(phenotype_at_mri),
    centre_code           = factor(centre_code),
    relapse_preceding_3mo = factor(relapse_preceding_3mo, levels = c(FALSE, TRUE)),
    preg_status           = factor(preg_status, levels = c(0, 1)),
    het_met_preceding_yr  = factor(het_met_preceding_yr, levels = c(FALSE, TRUE))
  )

str(combined_data)

# Check numbers
combined_data %>% count(preg_status) # n=74 pregnant, n=129 non-pregnant

#----------------------------------------------------------------------------------------------------------------------------------
#### 2. PROPENSITY SCORE MATCHING ####
# Nearest-neighbour matching without replacement
psm_1to1 <- matchit(
  preg_status ~ age_at_mri +
    phenotype_at_mri +
    disease_duration_at_mri +
    edss_closest +
    relapse_preceding_3mo +
    year +
    het_met_preceding_yr,
  data    = combined_data,
  method  = "nearest",
  ratio   = 1,
  caliper = 0.2,
  replace = FALSE
)

summary(psm_1to1) # all absolute SMDs < 0.15

# Extract matched cohort
matched_1to1 <- match.data(psm_1to1)
matched_1to1 %>% count(preg_status) # n=70 per group
