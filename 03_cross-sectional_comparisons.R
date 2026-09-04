# Cross-sectional comparison of global brain volumes: early postpartum vs non-pregnant MS cohorts
# Aim: Compare global brain volumes between early postpartum women with MS and
#      propensity score-matched non-pregnant women with MS

#### CONTENTS #####
# 1. WORKSPACE SETUP
# 2. CHECK NORMALITY
# 3. COMPARE GLOBAL BRAIN VOLUMES

#----------------------------------------------------------------------------------------------------------------------------------
#### 1. WORKSPACE SETUP ####
rm(list = ls())
library(tidyverse)

# Input: one row per woman, one scan each (first postpartum MRI for the pregnant cohort)
# Columns: subject, gm, wm, csf (ml), bpf (ratio), tiv (ml)
postpartum_volumes  <- read_csv("ms_early_postpartum_volumes.csv")
nonpregnant_volumes <- read_csv("ms_nonpregnant_volumes.csv")

# Check numbers
nrow(postpartum_volumes)  # n=70
nrow(nonpregnant_volumes) # n=70

#----------------------------------------------------------------------------------------------------------------------------------
#### 2. CHECK NORMALITY ####
metrics <- c("gm", "wm", "csf", "bpf")

long_cross <- bind_rows(
  postpartum_volumes  %>% mutate(group = "Postpartum"),
  nonpregnant_volumes %>% mutate(group = "Non-pregnant")
) %>%
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value")

# Shapiro-Wilk tests
long_cross %>%
  group_by(group, metric) %>%
  summarise(w = shapiro.test(value)$statistic,
            p = shapiro.test(value)$p.value,
            .groups = "drop")

# Histograms
ggplot(long_cross, aes(x = value)) +
  geom_histogram(bins = 15, fill = "#1f3a5f", colour = "white") +
  facet_wrap(~ metric + group, scales = "free", ncol = 2) +
  theme_classic()

# QQ plots
ggplot(long_cross, aes(sample = value)) +
  stat_qq(size = 0.8) +
  stat_qq_line(colour = "grey50") +
  facet_grid(metric ~ group, scales = "free") +
  theme_classic() # all normal

#----------------------------------------------------------------------------------------------------------------------------------
#### 3. COMPARE GLOBAL BRAIN VOLUMES ####
# Mean and SD for each metric
cross_stats <- long_cross %>%
  group_by(group, metric) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            sd   = sd(value, na.rm = TRUE),
            .groups = "drop")

print(as.data.frame(cross_stats))

# Between-subject SD, for comparison with within-woman longitudinal change
between_sd <- long_cross %>%
  group_by(metric) %>%
  summarise(sd_pooled = sd(value, na.rm = TRUE), .groups = "drop")

print(as.data.frame(between_sd))

# Welch's t-tests for all metrics
for (m in metrics) {
  cat("\n---", m, "---\n") # header for each output
  print(t.test(postpartum_volumes[[m]], nonpregnant_volumes[[m]]))
}
