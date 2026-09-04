# Longitudinal structural brain changes associated with pregnancy in women with multiple sclerosis

Statistical analysis code for Shipley et al., Brain Communications, 2026.

A multicentre cohort study of structural brain change across the peripregnancy
period in 127 women with multiple sclerosis, using data from the MSBase Registry
and MSBase Imaging Repository.

## Data

The data are not included here. Access to the deidentified study data is subject
to the approval of the individual Principal Investigators at each MSBase centre.
The scripts therefore cannot be run as they stand, and are provided to document
the analyses reported in the manuscript.

## Scripts

**`01_lme_model_global_metrics.R`**
Linear mixed-effects models of global brain volumes across four timepoints, with
a random intercept for participant. Produces Figure 1 and Table 2.

**`02_propensity_score_matching.R`**
1:1 nearest-neighbour propensity score matching of early postpartum to
non-pregnant women. Produces the matched cohort, the supplementary balance table
and the Love plot.

**`03_cross-sectional_comparisons.R`**
Welch's t-tests comparing global brain volumes between the matched cohorts.
Produces Table 5.

## Voxel-based and surface-based morphometry

Regional analyses were run in CAT12 (SPM25) using a second-level flexible
factorial design with subject and time as factors. The model specification files
include participant identifiers and are therefore not deposited here; the design,
covariates and thresholds are described in full in the manuscript Methods.

## Software

R v4.2.2 (R Foundation for Statistical Computing, Vienna, Austria) with
`tidyverse`, `lmerTest`, `emmeans`, `MatchIt` and `cobalt`.

CAT12 v12.9, implemented in SPM25 v25.01 (Wellcome Centre for Human
Neuroimaging, UCL, London, UK) and run in MATLAB R2024b (MathWorks, Natick, MA).

## Contact

Corresponding author: Prof Vilija G. Jokubaitis
[Vilija.jokubaitis@monash.edu](mailto:Vilija.jokubaitis@monash.edu)
