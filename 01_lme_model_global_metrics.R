# Longitudinal analysis of global brain volumes across the peripregnancy period
# Aim: Estimate change in global brain volumes from baseline to each postpartum 
#      timepoint in women with MS, using linear mixed-effects models

#### CONTENTS #####
# 1. WORKSPACE SETUP
# 2. LINEAR MIXED-EFFECTS MODELS
# 3. SUMMARY TABLE

#----------------------------------------------------------------------------------------------------------------------------------
#### 1. WORKSPACE SETUP ####
rm(list = ls())
library(tidyverse)
library(lmerTest)
library(emmeans)

# Input: one row per scan (345 scans from 127 women)
# Columns: subject, timepoint, gm, wm, csf, bpf (volumes in ml; BPF a ratio)
# Covariates are of two kinds:
#   Time-invariant (one value per woman, constant across all her scans):
#     age and disease duration at baseline, MS phenotype, MSBase centre,
#     baseline total intracranial volume, pregnancy duration, year of baseline MRI
#   Time-varying (a separate value for each scan, so the covariate can change
#   within a woman across the peripregnancy period):
#     scanner (composite of manufacturer and field strength), relapse within
#     3 months before the scan, exposure to moderate- or high-efficacy DMT
brain_volumes <- read_csv("peripregnancy_global_brain_volumes.csv")

# Set categorical covariates as factors
brain_volumes <- brain_volumes %>%
  mutate(across(c(scanner, relapse, dmt, centre, phenotype), factor))

# Baseline as reference level
brain_volumes <- brain_volumes %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("baseline", "preg", "pp0-3", "pp3-12", "pp12-24")))

# Check numbers
nrow(brain_volumes)              # n=345 scans
n_distinct(brain_volumes$subject) # n=127 women
table(brain_volumes$timepoint)   # n=127, 8, 74, 80, 56

#----------------------------------------------------------------------------------------------------------------------------------
#### 2. LINEAR MIXED-EFFECTS MODELS ####
lmm_terms <- c("timepoint", "age", "disease_dur", "phenotype", "centre", "tiv_cov",
               "year", "scanner", "relapse", "dmt", "preg_dur", "(1|subject)")

run_lmm <- function(metric) { # fitted to "gm", "wm", "csf"
  m   <- lmer(reformulate(lmm_terms, response = metric), data = brain_volumes)
  con <- contrast(emmeans(m, ~ timepoint), "trt.vs.ctrl") # postpartum timepoints against baseline
  con <- con[!grepl("preg", summary(con)$contrast)]       # drop pregnancy before adjusting
  ci  <- as.data.frame(confint(con, adjust = "none"))     # unadjusted 95% CIs
  p   <- as.data.frame(summary(con, adjust = "holm"))     # estimates, t, df, Holm-adjusted P
  
  tibble(
    metric     = metric,
    comparison = p$contrast,
    estimate   = p$estimate,
    ci_low     = ci$lower.CL,
    ci_high    = ci$upper.CL,
    t          = p$t.ratio,
    df         = p$df,
    p_holm     = p$p.value
  )
}

lme_results_table <- map_dfr(c("gm", "wm", "csf"), run_lmm) # see below for BPF

print(lme_results_table %>%
        mutate(across(c(estimate, ci_low, ci_high, t), ~ round(.x, 2)),
               df = round(df, 1)), n = Inf)

# Estimated marginal means
emm_all <- map_dfr(c("gm", "wm", "csf"), function(metric) {
  m <- lmer(reformulate(lmm_terms, response = metric), data = brain_volumes)
  as.data.frame(emmeans(m, ~ timepoint)) %>% mutate(metric = metric)
})

print(emm_all, digits = 5)

# BPF model, excluding TIV as BPF is already normalised to it
m_bpf <- lmer(reformulate(setdiff(lmm_terms, "tiv_cov"), response = "bpf"),
              data = brain_volumes)

con_bpf <- contrast(emmeans(m_bpf, ~ timepoint), "trt.vs.ctrl")
con_bpf <- con_bpf[!grepl("preg", summary(con_bpf)$contrast)] # drop pregnancy before adjusting

bpf_out <- cbind(
  as.data.frame(summary(con_bpf, adjust = "holm"))[, c("contrast", "estimate", "t.ratio", "p.value")],
  as.data.frame(confint(con_bpf, adjust = "none"))[, c("lower.CL", "upper.CL")]
)

print(bpf_out, digits = 4)

emm_bpf <- as.data.frame(emmeans(m_bpf, ~ timepoint))
emm_bpf

#----------------------------------------------------------------------------------------------------------------------------------
#### 3. SUMMARY TABLE ####
# Combine contrasts from both models
contrasts_all <- bind_rows(
  lme_results_table %>%
    select(metric, comparison, estimate, ci_low, ci_high, t, p_holm),
  bpf_out %>%
    transmute(metric     = "bpf",
              comparison = contrast,
              estimate,
              ci_low     = lower.CL,
              ci_high    = upper.CL,
              t          = t.ratio,
              p_holm     = p.value)
)

# Combine estimated marginal means from both models
emm_combined <- bind_rows(emm_all, emm_bpf %>% mutate(metric = "bpf"))

# Relative change, as a percentage of the estimated baseline value
base_all <- emm_combined %>%
  filter(timepoint == "baseline") %>%
  select(metric, base = emmean)

contrasts_all <- contrasts_all %>%
  left_join(base_all, by = "metric") %>%
  mutate(pct         = round(estimate / base * 100, 2),
         pct_ci_low  = round(ci_low   / base * 100, 2),
         pct_ci_high = round(ci_high  / base * 100, 2)) %>%
  select(-base)

# Summary table: estimated means, change from baseline, and relative change
summary_table <- emm_combined %>%
  select(metric, timepoint, emmean, emm_low = lower.CL, emm_high = upper.CL) %>%
  filter(timepoint != "preg") %>%
  left_join(
    contrasts_all %>%
      mutate(timepoint = sub("^\\(?([^)]*)\\)? - baseline$", "\\1", comparison)) %>%
      select(-comparison),
    by = c("metric", "timepoint")
  ) %>%
  arrange(factor(metric, levels = c("bpf", "gm", "wm", "csf")),
          factor(timepoint, levels = c("baseline", "pp0-3", "pp3-12", "pp12-24"))) %>%
  as_tibble()

print(summary_table, n = Inf)
