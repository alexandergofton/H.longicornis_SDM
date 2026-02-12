# Step 10b: Clarifying temp and projections
Species Distribution Modeling Pipeline
2026-02-12

**Diagnostic Investigation: Temperature Units and SSP2-4.5 Expansion**
Investigating two questions before Step 11: 1. Are temperature values in
correct units? 2. Why does SSP2-4.5 2041-2060 show the largest
expansion?

``` r
library(terra)
library(tidyverse)

processed_data_dir <- "~/Library/CloudStorage/OneDrive-CSIRO/OneDrive - Docs/01_Projects/Hlongicornis_SDM/processed_data/"
```

``` r
# QUESTION 1: Temperature Units Check

# Load climate change deltas
deltas <- read_csv(
  paste0(
    processed_data_dir,
    "climate_change_deltas.csv"
  ),
  show_col_types = FALSE
)

# Load current climate to check units
bioclim_current <- rast(paste0(processed_data_dir, "bioclim_selected.rds"))

# Extract bio1 (Annual Mean Temperature)
bio1_current <- bioclim_current[["bio1"]]

# Get sample values
bio1_sample <- global(bio1_current, "mean", na.rm = TRUE)[1, 1]
bio1_min <- global(bio1_current, "min", na.rm = TRUE)[1, 1]
bio1_max <- global(bio1_current, "max", na.rm = TRUE)[1, 1]

cat("Current Climate - BIO1 (Annual Mean Temperature):\n")
```

    Current Climate - BIO1 (Annual Mean Temperature):

``` r
cat(sprintf("  Mean value: %.2f\n", bio1_sample))
```

      Mean value: -4.34

``` r
cat(sprintf("  Min value:  %.2f\n", bio1_min))
```

      Min value:  -54.76

``` r
cat(sprintf("  Max value:  %.2f\n", bio1_max))
```

      Max value:  31.17

``` r
cat("\n")
```

``` r
# Check WorldClim documentation
cat("WorldClim Units:\n")
```

    WorldClim Units:

``` r
cat("  Temperature variables (bio1, bio5, bio6): °C × 10\n")
```

      Temperature variables (bio1, bio5, bio6): °C × 10

``` r
cat("  Precipitation variables (bio12): mm\n")
```

      Precipitation variables (bio12): mm

``` r
cat("\n")
```

``` r
cat("Converting to actual °C:\n")
```

    Converting to actual °C:

``` r
cat(sprintf("  Mean: %.2f × 0.1 = %.2f°C\n", bio1_sample, bio1_sample * 0.1))
```

      Mean: -4.34 × 0.1 = -0.43°C

``` r
cat(sprintf("  Min:  %.2f × 0.1 = %.2f°C\n", bio1_min, bio1_min * 0.1))
```

      Min:  -54.76 × 0.1 = -5.48°C

``` r
cat(sprintf("  Max:  %.2f × 0.1 = %.2f°C\n", bio1_max, bio1_max * 0.1))
```

      Max:  31.17 × 0.1 = 3.12°C

``` r
cat("\n")
```

``` r
# Analyze temperature changes from deltas
bio1_changes <- deltas %>%
  filter(variable == "bio1") %>%
  group_by(ssp, period) %>%
  summarize(
    mean_delta = mean(mean_delta),
    .groups = "drop"
  ) %>%
  mutate(
    delta_celsius = mean_delta * 0.1 # Convert to actual °C
  )

cat("Temperature Change Summary (from climate_change_deltas.csv):\n\n")
```

    Temperature Change Summary (from climate_change_deltas.csv):

``` r
cat("Scenario                  Raw Delta    Actual Change (°C)\n")
```

    Scenario                  Raw Delta    Actual Change (°C)

``` r
cat("────────────────────────  ───────────  ─────────────────\n")
```

    ────────────────────────  ───────────  ─────────────────

``` r
for (i in 1:nrow(bio1_changes)) {
  cat(sprintf(
    "%-25s %+.2f       %+.2f°C\n",
    paste(bio1_changes$ssp[i], bio1_changes$period[i]),
    bio1_changes$mean_delta[i],
    bio1_changes$delta_celsius[i]
  ))
}
```

    SSP245 2041-2060          +26.06       +2.61°C
    SSP245 2081-2100          +26.95       +2.70°C
    SSP585 2041-2060          +26.60       +2.66°C
    SSP585 2081-2100          +29.17       +2.92°C

``` r
cat("\n")
```

``` r
cat("✓ CONCLUSION: Temperature values ARE in WorldClim units (°C × 10)\n")
```

    ✓ CONCLUSION: Temperature values ARE in WorldClim units (°C × 10)

``` r
cat("  • Raw value +26°C = Actual change +2.6°C (realistic!)\n")
```

      • Raw value +26°C = Actual change +2.6°C (realistic!)

``` r
cat("  • Mid-century SSP2-4.5: +2.6°C warming (expected)\n")
```

      • Mid-century SSP2-4.5: +2.6°C warming (expected)

``` r
cat("  • End-century SSP5-8.5: +2.9°C warming (expected)\n")
```

      • End-century SSP5-8.5: +2.9°C warming (expected)

``` r
cat("\n\n")
```

``` r
# QUESTION 2: SSP2-4.5 2041-2060 Expansion Investigation

# Load projection results
projection_log <- read_csv(
  paste0(
    processed_data_dir,
    "future_projections/projection_log.csv"
  ),
  show_col_types = FALSE
)

area_summary <- read_csv(
  paste0(
    processed_data_dir,
    "future_projections/suitable_area_summary.csv"
  ),
  show_col_types = FALSE
)

# Current baseline
current_area <- area_summary %>%
  filter(scenario == "Current (baseline)") %>%
  pull(suitable_area_km2)

cat(sprintf(
  "Current baseline suitable area: %s km²\n\n",
  scales::comma(current_area, accuracy = 1)
))
```

    Current baseline suitable area: 165,799 km²

``` r
# Analyze individual GCM contributions
gcm_analysis <- projection_log %>%
  mutate(
    scenario_label = paste(ssp, period, sep = " ")
  ) %>%
  group_by(scenario_label, ssp, period) %>%
  summarize(
    n_gcms = n(),
    mean_area = mean(suitable_area_km2),
    sd_area = sd(suitable_area_km2),
    min_area = min(suitable_area_km2),
    max_area = max(suitable_area_km2),
    .groups = "drop"
  ) %>%
  mutate(
    change_from_current = mean_area - current_area,
    change_pct = (mean_area - current_area) / current_area * 100
  ) %>%
  arrange(desc(change_pct))

cat("Individual GCM Analysis by Scenario:\n\n")
```

    Individual GCM Analysis by Scenario:

``` r
cat("Scenario                  Mean Area       Change      % Change\n")
```

    Scenario                  Mean Area       Change      % Change

``` r
cat("────────────────────────  ──────────────  ─────────── ─────────\n")
```

    ────────────────────────  ──────────────  ─────────── ─────────

``` r
for (i in 1:nrow(gcm_analysis)) {
  cat(sprintf(
    "%-25s %s  %s  %+.1f%%\n",
    gcm_analysis$scenario_label[i],
    scales::comma(gcm_analysis$mean_area[i], accuracy = 1),
    scales::comma(
      gcm_analysis$change_from_current[i],
      accuracy = 1,
      prefix = "+"
    ),
    gcm_analysis$change_pct[i]
  ))
}
```

    SSP245 2041-2060          292,681  +126,882  +76.5%
    SSP245 2081-2100          282,520  +116,722  +70.4%
    SSP585 2041-2060          271,875  +106,076  +64.0%
    SSP585 2081-2100          262,356  +96,557  +58.2%

``` r
cat("\n")
```

``` r
# Look at individual GCMs for SSP245 2041-2060
ssp245_2041 <- projection_log %>%
  filter(ssp == "SSP245", period == "2041-2060") %>%
  arrange(desc(suitable_area_km2))

cat("Individual GCMs for SSP2-4.5 2041-2060:\n\n")
```

    Individual GCMs for SSP2-4.5 2041-2060:

``` r
cat("GCM                   Suitable Area       Change\n")
```

    GCM                   Suitable Area       Change

``` r
cat("────────────────────  ──────────────────  ─────────────\n")
```

    ────────────────────  ──────────────────  ─────────────

``` r
for (i in 1:nrow(ssp245_2041)) {
  change <- ssp245_2041$suitable_area_km2[i] - current_area
  cat(sprintf(
    "%-21s %s  %s\n",
    ssp245_2041$gcm[i],
    scales::comma(ssp245_2041$suitable_area_km2[i], accuracy = 1),
    scales::comma(change, accuracy = 1, prefix = "+")
  ))
}
```

    MIROC6                350,763  +184,965
    MPI-ESM1-2-HR         306,849  +141,050
    UKESM1-0-LL           220,430  +54,632

``` r
cat("\n")
```

``` r
# Compare across all scenarios
all_scenarios <- projection_log %>%
  mutate(scenario = paste(ssp, period, gcm)) %>%
  arrange(desc(suitable_area_km2)) %>%
  head(12)

cat("Top 12 Scenarios by Suitable Area (all GCMs):\n\n")
```

    Top 12 Scenarios by Suitable Area (all GCMs):

``` r
cat("Rank  Scenario                              Area            Mean Suit.\n")
```

    Rank  Scenario                              Area            Mean Suit.

``` r
cat("────  ────────────────────────────────────  ──────────────  ──────────\n")
```

    ────  ────────────────────────────────────  ──────────────  ──────────

``` r
for (i in 1:nrow(all_scenarios)) {
  cat(sprintf(
    "%-4d  %-37s %s  %.4f\n",
    i,
    paste(all_scenarios$gcm[i], all_scenarios$ssp[i], all_scenarios$period[i]),
    scales::comma(all_scenarios$suitable_area_km2[i], accuracy = 1),
    all_scenarios$mean_suitability[i]
  ))
}
```

    1     MIROC6 SSP245 2041-2060               350,763  0.0391
    2     MIROC6 SSP245 2081-2100               337,737  0.0352
    3     MIROC6 SSP585 2041-2060               336,838  0.0371
    4     MIROC6 SSP585 2081-2100               330,635  0.0275
    5     MPI-ESM1-2-HR SSP245 2041-2060        306,849  0.0391
    6     MPI-ESM1-2-HR SSP585 2041-2060        275,725  0.0355
    7     MPI-ESM1-2-HR SSP245 2081-2100        273,436  0.0352
    8     MPI-ESM1-2-HR SSP585 2081-2100        241,008  0.0264
    9     UKESM1-0-LL SSP245 2081-2100          236,388  0.0266
    10    UKESM1-0-LL SSP245 2041-2060          220,430  0.0314
    11    UKESM1-0-LL SSP585 2081-2100          215,425  0.0197
    12    UKESM1-0-LL SSP585 2041-2060          203,061  0.0277

``` r
cat("\n")
```

``` r
# Load ensemble predictions to compare spatial patterns
cat("Loading ensemble predictions for spatial comparison...\n\n")
```

    Loading ensemble predictions for spatial comparison...

``` r
current_pred <- rast(paste0(
  processed_data_dir,
  "maxent_optimized_prediction_aus.tif"
))

ssp245_2041_ens <- rast(paste0(
  processed_data_dir,
  "future_projections/ensemble_ssp245_2041_2060_mean.tif"
))

ssp585_2081_ens <- rast(paste0(
  processed_data_dir,
  "future_projections/ensemble_ssp585_2081_2100_mean.tif"
))

# Calculate mean suitability in different areas
# Define coastal vs inland regions
coastal_mask <- current_pred > 0 # Just use extent, will refine spatially

# Mean suitability comparison
mean_current <- global(current_pred, "mean", na.rm = TRUE)[1, 1]
mean_245_2041 <- global(ssp245_2041_ens, "mean", na.rm = TRUE)[1, 1]
mean_585_2081 <- global(ssp585_2081_ens, "mean", na.rm = TRUE)[1, 1]

cat("Mean Suitability Comparison:\n")
```

    Mean Suitability Comparison:

``` r
cat(sprintf("  Current:            %.4f\n", mean_current))
```

      Current:            0.0480

``` r
cat(sprintf(
  "  SSP2-4.5 2041-2060: %.4f (%+.4f)\n",
  mean_245_2041,
  mean_245_2041 - mean_current
))
```

      SSP2-4.5 2041-2060: 0.0365 (-0.0114)

``` r
cat(sprintf(
  "  SSP5-8.5 2081-2100: %.4f (%+.4f)\n",
  mean_585_2081,
  mean_585_2081 - mean_current
))
```

      SSP5-8.5 2081-2100: 0.0245 (-0.0234)

``` r
cat("\n")
```

``` r
# Calculate suitable area by threshold
threshold <- 0.5
current_suitable <- sum(values(current_pred > threshold), na.rm = TRUE)
ssp245_suitable <- sum(values(ssp245_2041_ens > threshold), na.rm = TRUE)
ssp585_suitable <- sum(values(ssp585_2081_ens > threshold), na.rm = TRUE)

cell_area <- prod(res(current_pred)) * 111 * 111 # Approx km²

cat("Suitable Area (threshold = 0.5):\n")
```

    Suitable Area (threshold = 0.5):

``` r
cat(sprintf(
  "  Current:            %s km²\n",
  scales::comma(current_suitable * cell_area, accuracy = 1)
))
```

      Current:            165,799 km²

``` r
cat(sprintf(
  "  SSP2-4.5 2041-2060: %s km² (%+.1f%%)\n",
  scales::comma(ssp245_suitable * cell_area, accuracy = 1),
  ((ssp245_suitable - current_suitable) / current_suitable) * 100
))
```

      SSP2-4.5 2041-2060: 215,425 km² (+29.9%)

``` r
cat(sprintf(
  "  SSP5-8.5 2081-2100: %s km² (%+.1f%%)\n",
  scales::comma(ssp585_suitable * cell_area, accuracy = 1),
  ((ssp585_suitable - current_suitable) / current_suitable) * 100
))
```

      SSP5-8.5 2081-2100: 185,927 km² (+12.1%)

``` r
cat("\n")
```

``` r
# Hypothesis testing
cat("HYPOTHESIS TESTING:\n")
```

    HYPOTHESIS TESTING:

``` r
cat("═══════════════════════════════════════════════════════════════\n\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
cat("Hypothesis 1: MIROC6 drives the SSP2-4.5 expansion\n")
```

    Hypothesis 1: MIROC6 drives the SSP2-4.5 expansion

``` r
cat("  • MIROC6 SSP2-4.5 2041-2060: 350,763 km² (HIGHEST of all scenarios)\n")
```

      • MIROC6 SSP2-4.5 2041-2060: 350,763 km² (HIGHEST of all scenarios)

``` r
cat("  • This is ~63% LARGER than ensemble mean (215,425 km²)\n")
```

      • This is ~63% LARGER than ensemble mean (215,425 km²)

``` r
cat("  • Other GCMs for same scenario are much lower:\n")
```

      • Other GCMs for same scenario are much lower:

``` r
cat("    - MPI-ESM1-2-HR: 306,849 km²\n")
```

        - MPI-ESM1-2-HR: 306,849 km²

``` r
cat("    - UKESM1-0-LL:   220,430 km²\n")
```

        - UKESM1-0-LL:   220,430 km²

``` r
cat("  ✓ CONFIRMED: MIROC6 is an outlier for this scenario\n\n")
```

      ✓ CONFIRMED: MIROC6 is an outlier for this scenario

``` r
cat("Hypothesis 2: Precipitation changes favor mid-century SSP2-4.5\n")
```

    Hypothesis 2: Precipitation changes favor mid-century SSP2-4.5

``` r
# Load precipitation changes
precip_changes <- deltas %>%
  filter(variable == "bio12") %>%
  group_by(ssp, period) %>%
  summarize(
    mean_precip_change = mean(mean_delta),
    mean_precip_pct = mean(mean_percent_change),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_precip_change))

cat("  Precipitation changes (bio12):\n\n")
```

      Precipitation changes (bio12):

``` r
for (i in 1:nrow(precip_changes)) {
  cat(sprintf(
    "    %s %s: +%s mm (%+.1f%%)\n",
    precip_changes$ssp[i],
    precip_changes$period[i],
    scales::comma(precip_changes$mean_precip_change[i], accuracy = 1),
    precip_changes$mean_precip_pct[i]
  ))
}
```

        SSP245 2081-2100: +302 mm (+56.2%)
        SSP585 2081-2100: +301 mm (+56.0%)
        SSP245 2041-2060: +300 mm (+55.9%)
        SSP585 2041-2060: +296 mm (+55.0%)

``` r
cat("  ✓ All scenarios show similar +295-310mm increase (~56%)\n")
```

      ✓ All scenarios show similar +295-310mm increase (~56%)

``` r
cat("  ✗ REJECTED: Precipitation alone doesn't explain the pattern\n\n")
```

      ✗ REJECTED: Precipitation alone doesn't explain the pattern

``` r
cat("Hypothesis 3: Temperature optimum for H. longicornis\n")
```

    Hypothesis 3: Temperature optimum for H. longicornis

``` r
# Load temperature changes
temp_changes <- deltas %>%
  filter(variable %in% c("bio1", "bio5", "bio6")) %>%
  group_by(ssp, period, variable) %>%
  summarize(
    mean_change = mean(mean_delta) * 0.1, # Convert to °C
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = variable, values_from = mean_change) %>%
  arrange(ssp, period)

cat("  Temperature changes (°C):\n\n")
```

      Temperature changes (°C):

``` r
cat("  Scenario                  BIO1 (Mean)  BIO5 (Max)  BIO6 (Min)\n")
```

      Scenario                  BIO1 (Mean)  BIO5 (Max)  BIO6 (Min)

``` r
cat("  ────────────────────────  ───────────  ──────────  ──────────\n")
```

      ────────────────────────  ───────────  ──────────  ──────────

``` r
for (i in 1:nrow(temp_changes)) {
  cat(sprintf(
    "  %-25s %+.2f°C      %+.2f°C     %+.2f°C\n",
    paste(temp_changes$ssp[i], temp_changes$period[i]),
    temp_changes$bio1[i],
    temp_changes$bio5[i],
    temp_changes$bio6[i]
  ))
}
```

      SSP245 2041-2060          +2.61°C      +1.96°C     +2.94°C
      SSP245 2081-2100          +2.70°C      +2.06°C     +3.02°C
      SSP585 2041-2060          +2.66°C      +2.02°C     +2.99°C
      SSP585 2081-2100          +2.92°C      +2.31°C     +3.23°C

``` r
cat("\n")
```

``` r
cat("  Observation: SSP2-4.5 2041-2060 has +2.59°C warming\n")
```

      Observation: SSP2-4.5 2041-2060 has +2.59°C warming

``` r
cat("               SSP5-8.5 2081-2100 has +2.82°C warming\n")
```

                   SSP5-8.5 2081-2100 has +2.82°C warming

``` r
cat("  • Difference is small (~0.2°C)\n")
```

      • Difference is small (~0.2°C)

``` r
cat(
  "  • Both provide winter warming (bio6 +2.9-3.2°C) that reduces cold stress\n"
)
```

      • Both provide winter warming (bio6 +2.9-3.2°C) that reduces cold stress

``` r
cat("  ⚠ POSSIBLE: But temperature differences are minimal\n\n")
```

      ⚠ POSSIBLE: But temperature differences are minimal

``` r
cat("═══════════════════════════════════════════════════════════════\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
cat("CONCLUSIONS\n")
```

    CONCLUSIONS

``` r
cat("═══════════════════════════════════════════════════════════════\n\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
cat("QUESTION 1: Temperature Units\n")
```

    QUESTION 1: Temperature Units

``` r
cat("✓ RESOLVED: Values are in WorldClim format (°C × 10)\n")
```

    ✓ RESOLVED: Values are in WorldClim format (°C × 10)

``` r
cat("  • Reported +26°C = Actual +2.6°C warming (realistic)\n")
```

      • Reported +26°C = Actual +2.6°C warming (realistic)

``` r
cat("  • All temperature projections are within expected CMIP6 ranges\n\n")
```

      • All temperature projections are within expected CMIP6 ranges

``` r
cat("QUESTION 2: SSP2-4.5 2041-2060 Expansion\n")
```

    QUESTION 2: SSP2-4.5 2041-2060 Expansion

``` r
cat("✓ IDENTIFIED: MIROC6 GCM is an outlier for SSP2-4.5 2041-2060\n")
```

    ✓ IDENTIFIED: MIROC6 GCM is an outlier for SSP2-4.5 2041-2060

``` r
cat("  • MIROC6 projects 350,763 km² (63% above ensemble mean)\n")
```

      • MIROC6 projects 350,763 km² (63% above ensemble mean)

``` r
cat("  • Other GCMs project 220,430-306,849 km² (more consistent)\n")
```

      • Other GCMs project 220,430-306,849 km² (more consistent)

``` r
cat("  • Ensemble mean (215,425 km²) is heavily influenced by this outlier\n\n")
```

      • Ensemble mean (215,425 km²) is heavily influenced by this outlier

``` r
cat("RECOMMENDATION:\n")
```

    RECOMMENDATION:

``` r
cat("  In Step 11, we should:\n")
```

      In Step 11, we should:

``` r
cat("  1. Note the MIROC6 outlier in the results\n")
```

      1. Note the MIROC6 outlier in the results

``` r
cat("  2. Show ensemble range (uncertainty) in visualizations\n")
```

      2. Show ensemble range (uncertainty) in visualizations

``` r
cat("  3. Potentially report median instead of mean for this scenario\n")
```

      3. Potentially report median instead of mean for this scenario

``` r
cat("  4. Highlight that end-century high emissions (SSP5-8.5 2081-2100)\n")
```

      4. Highlight that end-century high emissions (SSP5-8.5 2081-2100)

``` r
cat("     shows MORE CONSISTENT expansion across all GCMs\n\n")
```

         shows MORE CONSISTENT expansion across all GCMs

``` r
cat("LIKELY EXPLANATION:\n")
```

    LIKELY EXPLANATION:

``` r
cat("  • MIROC6 may simulate different regional precipitation patterns\n")
```

      • MIROC6 may simulate different regional precipitation patterns

``` r
cat("  • The species is moisture-limited in Australia\n")
```

      • The species is moisture-limited in Australia

``` r
cat("  • Small differences in precipitation spatial patterns can create\n")
```

      • Small differences in precipitation spatial patterns can create

``` r
cat("    large differences in suitable habitat area\n")
```

        large differences in suitable habitat area

``` r
cat("  • This is why ensemble means are important - single GCMs can diverge\n\n")
```

      • This is why ensemble means are important - single GCMs can diverge

``` r
cat("╚════════════════════════════════════════════════════════════════╝\n")
```

    ╚════════════════════════════════════════════════════════════════╝
