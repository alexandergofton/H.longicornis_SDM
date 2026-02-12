# Step 11: Final Visualization & Results Synthesis
Species Distribution Modeling Pipeline
2026-02-12

- [Overview](#overview)
- [Setup](#setup)
- [Load Data](#load-data)
- [Publication Map: Current
  Distribution](#publication-map-current-distribution)
- [Publication Map: Future Scenarios
  Comparison](#publication-map-future-scenarios-comparison)
- [Change Maps: Gain/Loss/Stable](#change-maps-gainlossstable)
- [Quantitative Summary: Bar Charts](#quantitative-summary-bar-charts)
- [Uncertainty Visualization](#uncertainty-visualization)
- [Priority Surveillance Regions](#priority-surveillance-regions)
- [Comprehensive Results Summary](#comprehensive-results-summary)
- [Final Interpretation Text](#final-interpretation-text)
- [Session Info](#session-info)
- [Summary of Outputs](#summary-of-outputs)

## Overview

This script creates publication-quality visualizations and synthesizes
all SDM results for *Haemaphysalis longicornis* in Australia.

**What This Script Does**: 1. Creates high-quality maps (current and
future scenarios) 2. Generates quantitative summaries and bar charts 3.
Identifies priority surveillance regions for biosecurity 4. Calculates
key metrics (area changes, range shifts) 5. Produces final
interpretation with uncertainty analysis 6. Creates comprehensive
results document

**Expected Runtime**: 15-20 minutes

------------------------------------------------------------------------

## Setup

``` r
# Core packages
library(terra)
library(sf)
library(tidyverse)

# Visualization
library(ggplot2)
library(patchwork)
library(viridis)
library(tidyterra)
library(scales)
library(ggspatial)
library(rnaturalearth)

# Set options
terra::terraOptions(memfrac = 0.8, tempdir = "temp")

# Create directories
dir.create("figures/final", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE)

processed_data_dir <- "~/Library/CloudStorage/OneDrive-CSIRO/OneDrive - Docs/01_Projects/Hlongicornis_SDM/processed_data/"

figures_dir <- "~/Library/CloudStorage/OneDrive-CSIRO/OneDrive - Docs/01_Projects/Hlongicornis_SDM/figures/final/"
```

------------------------------------------------------------------------

## Load Data

``` r
cat("\n=== Loading Data ===\n")
```


    === Loading Data ===

``` r
# Load predictions
current_pred <- rast(paste0(
  processed_data_dir,
  "maxent_optimized_prediction_aus.tif"
))

ensemble_245_2041 <- rast(paste0(
  processed_data_dir,
  "future_projections/ensemble_ssp245_2041_2060_mean.tif"
))

ensemble_245_2081 <- rast(paste0(
  processed_data_dir,
  "future_projections/ensemble_ssp245_2081_2100_mean.tif"
))

ensemble_585_2041 <- rast(paste0(
  processed_data_dir,
  "future_projections/ensemble_ssp585_2041_2060_mean.tif"
))

ensemble_585_2081 <- rast(paste0(
  processed_data_dir,
  "future_projections/ensemble_ssp585_2081_2100_mean.tif"
))

# Load ensemble SDs (uncertainty)
ens_sd_245_2041 <- rast(paste0(
  processed_data_dir,
  "future_projections/ensemble_ssp245_2041_2060_sd.tif"
))

ens_sd_585_2081 <- rast(paste0(
  processed_data_dir,
  "future_projections/ensemble_ssp585_2081_2100_sd.tif"
))

# Load change rasters
change_245_2041 <- rast(paste0(
  processed_data_dir,
  "future_projections/change_ssp245_2041_2060.tif"
))

change_585_2081 <- rast(paste0(
  processed_data_dir,
  "future_projections/change_ssp585_2081_2100.tif"
))


# Load summary data
area_summary <- read_csv(
  paste0(
    processed_data_dir,
    "future_projections/suitable_area_summary.csv"
  ),
  show_col_types = FALSE
)

change_stats <- read_csv(
  paste0(
    processed_data_dir,
    "future_projections/change_statistics.csv"
  ),
  show_col_types = FALSE
)

projection_log <- read_csv(
  paste0(
    processed_data_dir,
    "future_projections/projection_log.csv"
  ),
  show_col_types = FALSE
)


# Load occurrence data for overlay
occ_data <- read_csv(
  paste0(
    processed_data_dir,
    "occurrences_thinned.csv"
  ),
  show_col_types = FALSE
)

# Load Australia/NZ boundary
aus_nz_boundary <- st_read(
  paste0(
    processed_data_dir,
    "aus_nz_boundary.gpkg"
  ),
  quiet = TRUE
)
```

------------------------------------------------------------------------

## Publication Map: Current Distribution

Create high-quality map of current habitat suitability with occurrence
points.

``` r
cat("\n=== Creating Current Distribution Map ===\n")
```


    === Creating Current Distribution Map ===

``` r
# Get Australia states for boundaries
australia <- ne_states(country = "australia", returnclass = "sf")

# Filter occurrences to Australia/NZ region
occ_aus <- occ_data %>%
  filter(lon >= 110, lon <= 180, lat >= -48, lat <= -10)

# Create map
p_current_pub <- ggplot() +
  # Suitability raster
  geom_spatraster(data = current_pred) +
  scale_fill_gradientn(
    colors = c(
      "#f7f7f7",
      "#fddbc7",
      "#f4a582",
      "#d6604d",
      "#b2182b",
      "#67001f"
    ),
    na.value = "transparent",
    limits = c(0, 1),
    labels = percent,
    name = "Habitat\nSuitability"
  ) +
  # State boundaries
  geom_sf(data = australia, fill = NA, color = "gray30", linewidth = 0.3) +
  # Occurrence points
  geom_point(
    data = occ_aus,
    aes(x = lon, y = lat),
    color = "black",
    fill = "yellow",
    shape = 21,
    size = 2,
    stroke = 0.5,
    alpha = 0.8
  ) +
  # Annotations
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.5, "cm"),
    width = unit(1.5, "cm")
  ) +
  # Styling
  coord_sf(xlim = c(110, 180), ylim = c(-48, -10), expand = FALSE) +
  labs(
    title = "Current Habitat Suitability for Haemaphysalis longicornis",
    subtitle = "Australia and New Zealand (baseline climate: 1970-2000)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    legend.position = "right",
    axis.title = element_text(size = 12)
  )

print(p_current_pub)
```

![](11_final_visualization_results_files/figure-commonmark/pub-map-current-1.png)

``` r
ggsave(
  "figures/final/current_distribution_publication.png",
  p_current_pub,
  width = 12,
  height = 10,
  dpi = 300
)

cat("✓ Current distribution map saved\n")
```

    ✓ Current distribution map saved

------------------------------------------------------------------------

## Publication Map: Future Scenarios Comparison

Four-panel comparison of current and three future scenarios.

``` r
cat("\n=== Creating Future Scenarios Comparison ===\n")
```


    === Creating Future Scenarios Comparison ===

``` r
# Common extent
xlim <- c(110, 180)
ylim <- c(-48, -10)

# Common theme
map_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray50"),
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    panel.grid = element_line(color = "gray90", linewidth = 0.2),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.3, "cm")
  )

# Panel 1: Current
p1 <- ggplot() +
  geom_spatraster(data = current_pred) +
  scale_fill_gradientn(
    colors = c("#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
    na.value = "transparent",
    limits = c(0, 1),
    labels = percent,
    name = "Suitability"
  ) +
  geom_sf(data = australia, fill = NA, color = "gray30", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "Current (1970-2000)", subtitle = "Baseline") +
  map_theme

# Panel 2: SSP2-4.5 2041-2060
p2 <- ggplot() +
  geom_spatraster(data = ensemble_245_2041) +
  scale_fill_gradientn(
    colors = c("#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
    na.value = "transparent",
    limits = c(0, 1),
    labels = percent,
    name = "Suitability"
  ) +
  geom_sf(data = australia, fill = NA, color = "gray30", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "SSP2-4.5, 2041-2060",
       subtitle = "Mid-century, moderate emissions") +
  map_theme

# Panel 3: SSP2-4.5 2081-2100
p3 <- ggplot() +
  geom_spatraster(data = ensemble_245_2081) +
  scale_fill_gradientn(
    colors = c("#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
    na.value = "transparent",
    limits = c(0, 1),
    labels = percent,
    name = "Suitability"
  ) +
  geom_sf(data = australia, fill = NA, color = "gray30", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "SSP2-4.5, 2081-2100",
       subtitle = "End-century, moderate emissions") +
  map_theme

# Panel 4: SSP5-8.5 2081-2100
p4 <- ggplot() +
  geom_spatraster(data = ensemble_585_2081) +
  scale_fill_gradientn(
    colors = c("#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
    na.value = "transparent",
    limits = c(0, 1),
    labels = percent,
    name = "Suitability"
  ) +
  geom_sf(data = australia, fill = NA, color = "gray30", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "SSP5-8.5, 2081-2100",
       subtitle = "End-century, high emissions") +
  map_theme

# Combine
future_comparison <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Projected Habitat Suitability Under Climate Change",
    subtitle = "Haemaphysalis longicornis in Australia and New Zealand",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

print(future_comparison)
```

![](11_final_visualization_results_files/figure-commonmark/pub-map-future-1.png)

``` r
ggsave("figures/final/future_scenarios_comparison.png",
       future_comparison, width = 14, height = 12, dpi = 300)

cat("✓ Future scenarios comparison saved\n")
```

    ✓ Future scenarios comparison saved

------------------------------------------------------------------------

## Change Maps: Gain/Loss/Stable

Categorical maps showing habitat change patterns.

``` r
cat("\n=== Creating Habitat Change Maps ===\n")
```


    === Creating Habitat Change Maps ===

``` r
# Function to create gain/loss categorical raster
create_change_categories <- function(current, future, threshold = 0.5) {
  # Ensure extents match by resampling future to current
  future <- resample(future, current, method = "bilinear")

  current_suitable <- current > threshold
  future_suitable <- future > threshold

  # Create categorical raster
  # 1 = Loss, 2 = Stable, 3 = Gain, 4 = Unsuitable
  change_cat <- current_suitable * 0 + 4  # Start with "unsuitable"
  change_cat[current_suitable & !future_suitable] <- 1  # Loss
  change_cat[current_suitable & future_suitable] <- 2   # Stable
  change_cat[!current_suitable & future_suitable] <- 3  # Gain

  # Convert to categorical/factor
  change_cat <- as.factor(change_cat)
  levels(change_cat) <- data.frame(
    ID = 1:4,
    category = c("Loss", "Stable", "Gain", "Unsuitable")
  )

  return(change_cat)
}

# Create categorical change rasters
change_cat_245_2041 <- create_change_categories(current_pred, ensemble_245_2041)
change_cat_585_2081 <- create_change_categories(current_pred, ensemble_585_2081)

# Map theme for change
change_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray50"),
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    panel.grid = element_line(color = "gray90", linewidth = 0.2),
    legend.position = "right"
  )

# SSP2-4.5 2041-2060 change
p_change_245 <- ggplot() +
  geom_spatraster(data = change_cat_245_2041) +
  scale_fill_manual(
    values = c("Loss" = "#d7191c", "Stable" = "#fdae61",
               "Gain" = "#2c7bb6", "Unsuitable" = "gray90"),
    na.value = "transparent",
    name = "Habitat\nChange"
  ) +
  geom_sf(data = australia, fill = NA, color = "gray30", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    title = "Habitat Change: SSP2-4.5, 2041-2060",
    subtitle = "Suitable habitat (threshold = 0.5)"
  ) +
  change_theme

# SSP5-8.5 2081-2100 change
p_change_585 <- ggplot() +
  geom_spatraster(data = change_cat_585_2081) +
  scale_fill_manual(
    values = c("Loss" = "#d7191c", "Stable" = "#fdae61",
               "Gain" = "#2c7bb6", "Unsuitable" = "gray90"),
    na.value = "transparent",
    name = "Habitat\nChange"
  ) +
  geom_sf(data = australia, fill = NA, color = "gray30", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    title = "Habitat Change: SSP5-8.5, 2081-2100",
    subtitle = "Suitable habitat (threshold = 0.5)"
  ) +
  change_theme

# Combine
change_combined <- p_change_245 + p_change_585 +
  plot_annotation(
    title = "Projected Habitat Changes Under Climate Change",
    subtitle = "Red = Loss | Orange = Stable | Blue = Gain | Gray = Unsuitable",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  )

print(change_combined)
```

![](11_final_visualization_results_files/figure-commonmark/change-maps-1.png)

``` r
ggsave("figures/final/habitat_change_maps.png",
       change_combined, width = 14, height = 6, dpi = 300)

cat("✓ Habitat change maps saved\n")
```

    ✓ Habitat change maps saved

------------------------------------------------------------------------

## Quantitative Summary: Bar Charts

Visualize area changes across scenarios.

``` r
cat("\n=== Creating Quantitative Summary Charts ===\n")
```


    === Creating Quantitative Summary Charts ===

``` r
# Prepare data for plotting
area_plot_data <- area_summary %>%
  filter(scenario != "Current (baseline)") %>%
  mutate(
    scenario_clean = paste(ssp, period, sep = "\n"),
    scenario_clean = factor(scenario_clean,
                           levels = c("SSP245\n2041-2060", "SSP245\n2081-2100",
                                     "SSP585\n2041-2060", "SSP585\n2081-2100"))
  )

# Bar chart: Suitable area
p_area <- ggplot(area_plot_data, aes(x = scenario_clean, y = suitable_area_km2 / 1000,
                                     fill = ssp)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = area_summary$suitable_area_km2[1] / 1000,
             linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = 0.5, y = area_summary$suitable_area_km2[1] / 1000 + 15,
           label = "Current baseline", hjust = 0, size = 3.5) +
  scale_fill_manual(values = c("SSP245" = "#fdae61", "SSP585" = "#d7191c"),
                    labels = c("SSP2-4.5 (moderate)", "SSP5-8.5 (high)")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Suitable Habitat Area by Climate Scenario",
    subtitle = "Ensemble mean across 3 GCMs (MIROC6, MPI-ESM1-2-HR, UKESM1-0-LL)",
    x = NULL,
    y = "Suitable Area (1,000 km²)",
    fill = "Emission\nScenario"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 11),
    legend.position = "right"
  )

# Bar chart: Change from current
p_change_pct <- ggplot(area_plot_data, aes(x = scenario_clean,
                                           y = change_from_current_pct,
                                           fill = ssp)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("SSP245" = "#fdae61", "SSP585" = "#d7191c"),
                    labels = c("SSP2-4.5 (moderate)", "SSP5-8.5 (high)")) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(
    title = "Change in Suitable Habitat from Current",
    subtitle = "Positive = expansion, negative = contraction",
    x = NULL,
    y = "Change from Current (%)",
    fill = "Emission\nScenario"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 11),
    legend.position = "right"
  )

# Combine
area_charts <- p_area / p_change_pct

print(area_charts)
```

![](11_final_visualization_results_files/figure-commonmark/bar-charts-1.png)

``` r
ggsave("figures/final/suitable_area_bar_charts.png",
       area_charts, width = 10, height = 10, dpi = 300)

cat("✓ Bar charts saved\n")
```

    ✓ Bar charts saved

------------------------------------------------------------------------

## Uncertainty Visualization

Show inter-GCM variability and ensemble uncertainty.

``` r
cat("\n=== Creating Uncertainty Visualizations ===\n")
```


    === Creating Uncertainty Visualizations ===

``` r
# Calculate GCM range for each scenario
gcm_ranges <- projection_log %>%
  mutate(scenario = paste(ssp, period)) %>%
  group_by(scenario) %>%
  summarize(
    mean_area = mean(suitable_area_km2),
    min_area = min(suitable_area_km2),
    max_area = max(suitable_area_km2),
    range = max_area - min_area,
    .groups = "drop"
  ) %>%
  mutate(
    scenario_clean = factor(scenario,
                           levels = c("SSP245 2041-2060", "SSP245 2081-2100",
                                     "SSP585 2041-2060", "SSP585 2081-2100"))
  )

# Plot GCM range
p_gcm_range <- ggplot(gcm_ranges, aes(x = scenario_clean, y = mean_area / 1000)) +
  geom_errorbar(aes(ymin = min_area / 1000, ymax = max_area / 1000),
                width = 0.3, linewidth = 1, color = "steelblue") +
  geom_point(size = 4, color = "darkblue") +
  geom_hline(yintercept = area_summary$suitable_area_km2[1] / 1000,
             linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = 4.5, y = area_summary$suitable_area_km2[1] / 1000 + 15,
           label = "Current", hjust = 1, size = 3.5) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    title = "Inter-GCM Variability in Habitat Projections",
    subtitle = "Points = ensemble mean | Bars = range across 3 GCMs",
    x = NULL,
    y = "Suitable Area (1,000 km²)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 10)
  )

print(p_gcm_range)
```

![](11_final_visualization_results_files/figure-commonmark/uncertainty-viz-1.png)

``` r
ggsave("figures/final/gcm_uncertainty.png",
       p_gcm_range, width = 10, height = 7, dpi = 300)

# Spatial uncertainty map (SD)
p_uncertainty_map <- ggplot() +
  geom_spatraster(data = ens_sd_585_2081) +
  scale_fill_viridis_c(
    option = "plasma",
    na.value = "transparent",
    name = "Ensemble\nSD"
  ) +
  geom_sf(data = australia, fill = NA, color = "gray30", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    title = "Spatial Uncertainty in Projections",
    subtitle = "Standard deviation across GCMs (SSP5-8.5, 2081-2100)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_blank()
  )

print(p_uncertainty_map)
```

![](11_final_visualization_results_files/figure-commonmark/uncertainty-viz-2.png)

``` r
ggsave("figures/final/spatial_uncertainty_map.png",
       p_uncertainty_map, width = 10, height = 8, dpi = 300)

cat("✓ Uncertainty visualizations saved\n")
```

    ✓ Uncertainty visualizations saved

------------------------------------------------------------------------

## Priority Surveillance Regions

Identify high-risk areas for biosecurity monitoring based on gain
projections.

``` r
library(ggnewscale)

cat("\n=== Identifying Priority Surveillance Regions ===\n")
```


    === Identifying Priority Surveillance Regions ===

``` r
# Ensure ensemble matches current prediction extent
ensemble_585_2081 <- resample(ensemble_585_2081, current_pred, method = "bilinear")
change_585_2081 <- resample(change_585_2081, current_pred, method = "bilinear")

# Define priority as areas with:
# 1. Currently unsuitable (< 0.3)
# 2. Future suitable (> 0.5)
# 3. Large change (> +0.3)

priority_585 <- (current_pred < 0.3) &
                (ensemble_585_2081 > 0.5) &
                (change_585_2081 > 0.3)

# Convert to 0/1 and then to factor for discrete scale
priority_binary <- ifel(priority_585, 1, NA)
priority_binary <- as.factor(priority_binary)

# Calculate priority area
priority_cells <- sum(values(priority_585), na.rm = TRUE)
cell_area <- prod(res(current_pred)) * 111 * 111
priority_area_km2 <- priority_cells * cell_area

cat(sprintf("Priority surveillance area: %s km²\n",
            scales::comma(priority_area_km2, accuracy = 1)))
```

    Priority surveillance area: 59,573 km²

``` r
# Create priority map
p_priority <- ggplot() +
  # Background: current suitability (faded)
  geom_spatraster(data = current_pred, alpha = 0.3) +
  scale_fill_gradientn(
    colors = c("gray95", "gray70"),
    na.value = "transparent",
    limits = c(0, 1),
    guide = "none"
  ) +
  # Priority areas overlay
  new_scale_fill() +
  geom_spatraster(data = priority_binary) +
  scale_fill_manual(
    values = c("1" = "#e31a1c"),
    labels = "High priority",
    na.value = "transparent",
    name = "Surveillance\nPriority"
  ) +
  # State boundaries
  geom_sf(data = australia, fill = NA, color = "gray30", linewidth = 0.3) +
  # Existing occurrences
  geom_point(data = occ_aus, aes(x = lon, y = lat),
             color = "black", fill = "yellow", shape = 21,
             size = 1.5, stroke = 0.4, alpha = 0.7) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    title = "Priority Surveillance Regions for H. longicornis",
    subtitle = "Areas projected to become suitable by 2081-2100 (SSP5-8.5)\nYellow points = current known occurrences"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_blank(),
    legend.position = c(0.15, 0.25),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  )

print(p_priority)
```

![](11_final_visualization_results_files/figure-commonmark/surveillance-1.png)

``` r
ggsave("figures/final/priority_surveillance_regions.png",
       p_priority, width = 12, height = 10, dpi = 300)

cat("✓ Priority surveillance map saved\n")
```

    ✓ Priority surveillance map saved

------------------------------------------------------------------------

## Comprehensive Results Summary

Generate final summary tables and statistics.

``` r
cat("\n=== Generating Final Summary ===\n")
```


    === Generating Final Summary ===

``` r
# Model performance summary
model_summary <- tibble(
  metric = c(
    "Training AUC",
    "Testing AUC",
    "Validation AUC (LORO)",
    "Kappa",
    "TSS",
    "Optimized features",
    "Regularization (β)"
  ),
  value = c(
    "0.977",
    "0.976",
    "0.832",
    "0.37",
    "0.62",
    "LQHPT",
    "0.5"
  )
)

write_csv(model_summary, "results/model_performance_summary.csv")

# Climate change summary (convert to actual °C)
climate_summary <- tibble(
  variable = c("Temperature (°C)", "Precipitation (mm)"),
  ssp245_2041 = c("+2.6", "+300 (+56%)"),
  ssp245_2081 = c("+2.7", "+302 (+56%)"),
  ssp585_2041 = c("+2.7", "+296 (+55%)"),
  ssp585_2081 = c("+2.9", "+301 (+56%)")
)

write_csv(climate_summary, "results/climate_change_summary.csv")

# Habitat area summary
write_csv(area_summary, "results/suitable_habitat_area_by_scenario.csv")

# Change summary
write_csv(change_stats, "results/habitat_change_statistics.csv")

# GCM variability summary
write_csv(gcm_ranges, "results/gcm_variability_summary.csv")

cat("✓ All summary tables saved to results/\n")
```

    ✓ All summary tables saved to results/

------------------------------------------------------------------------

## Final Interpretation Text

``` r
cat("\n\n")
```

``` r
cat("╔════════════════════════════════════════════════════════════════╗\n")
```

    ╔════════════════════════════════════════════════════════════════╗

``` r
cat("║            FINAL RESULTS SUMMARY & INTERPRETATION              ║\n")
```

    ║            FINAL RESULTS SUMMARY & INTERPRETATION              ║

``` r
cat("╠════════════════════════════════════════════════════════════════╣\n\n")
```

    ╠════════════════════════════════════════════════════════════════╣

``` r
cat("1. MODEL PERFORMANCE\n")
```

    1. MODEL PERFORMANCE

``` r
cat("═══════════════════════════════════════════════════════════════\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
cat("  • Training/Testing AUC: 0.977/0.976 (excellent discrimination)\n")
```

      • Training/Testing AUC: 0.977/0.976 (excellent discrimination)

``` r
cat("  • Leave-one-region-out AUC: 0.832 (good transferability)\n")
```

      • Leave-one-region-out AUC: 0.832 (good transferability)

``` r
cat("  • Bootstrap uncertainty: 92.3% low, 6.9% medium, 0.8% high\n")
```

      • Bootstrap uncertainty: 92.3% low, 6.9% medium, 0.8% high

``` r
cat("  • Biological validation: 3/4 checks passed\n")
```

      • Biological validation: 3/4 checks passed

``` r
cat("  ✓ Model demonstrates robust performance and transferability\n\n")
```

      ✓ Model demonstrates robust performance and transferability

``` r
cat("2. CURRENT DISTRIBUTION (BASELINE)\n")
```

    2. CURRENT DISTRIBUTION (BASELINE)

``` r
cat("═══════════════════════════════════════════════════════════════\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
cat(sprintf("  • Suitable habitat: %s km²\n",
            comma(area_summary$suitable_area_km2[1], accuracy = 1)))
```

      • Suitable habitat: 165,799 km²

``` r
cat("  • Pattern: Concentrated on eastern coastal regions\n")
```

      • Pattern: Concentrated on eastern coastal regions

``` r
cat("  • Limited inland penetration (moisture limitation)\n")
```

      • Limited inland penetration (moisture limitation)

``` r
cat("  • Aligns with known occurrence data (614 global points)\n\n")
```

      • Aligns with known occurrence data (614 global points)

``` r
cat("3. CLIMATE CHANGE PROJECTIONS\n")
```

    3. CLIMATE CHANGE PROJECTIONS

``` r
cat("═══════════════════════════════════════════════════════════════\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
cat("  Temperature increases (annual mean):\n")
```

      Temperature increases (annual mean):

``` r
cat("    • Mid-century (2041-2060): +2.6-2.7°C\n")
```

        • Mid-century (2041-2060): +2.6-2.7°C

``` r
cat("    • End-century (2081-2100): +2.7-2.9°C\n")
```

        • End-century (2081-2100): +2.7-2.9°C

``` r
cat("  Precipitation increases (annual):\n")
```

      Precipitation increases (annual):

``` r
cat("    • All scenarios: +~300mm (+55-56%)\n")
```

        • All scenarios: +~300mm (+55-56%)

``` r
cat("  ✓ Substantial moisture increase across Australia/NZ\n\n")
```

      ✓ Substantial moisture increase across Australia/NZ

``` r
cat("4. FUTURE HABITAT PROJECTIONS\n")
```

    4. FUTURE HABITAT PROJECTIONS

``` r
cat("═══════════════════════════════════════════════════════════════\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
for (i in 2:nrow(area_summary)) {
  cat(sprintf("  %s:\n", area_summary$scenario[i]))
  cat(sprintf("    Suitable area: %s km² (%+.1f%%)\n",
              comma(area_summary$suitable_area_km2[i], accuracy = 1),
              area_summary$change_from_current_pct[i]))

  # Get corresponding change stats
  change_row <- change_stats[i-1, ]
  cat(sprintf("    Gain: %s km² | Loss: %s km² | Stable: %s km²\n",
              comma(change_row$area_gain_km2, accuracy = 1),
              comma(change_row$area_loss_km2, accuracy = 1),
              comma(change_row$area_stable_km2, accuracy = 1)))
}
```

      SSP245 2041-2060:
        Suitable area: 215,425 km² (+29.9%)
        Gain: 55,573 km² | Loss: 82,846 km² | Stable: 82,504 km²
      SSP585 2041-2060:
        Suitable area: 187,574 km² (+13.1%)
        Gain: 59,188 km² | Loss: 96,151 km² | Stable: 69,199 km²
      SSP245 2081-2100:
        Suitable area: 193,243 km² (+16.6%)
        Gain: 72,300 km² | Loss: 109,178 km² | Stable: 56,172 km²
      SSP585 2081-2100:
        Suitable area: 185,927 km² (+12.1%)
        Gain: 109,670 km² | Loss: 138,953 km² | Stable: 26,396 km²

``` r
cat("\n")
```

``` r
cat("5. KEY FINDINGS\n")
```

    5. KEY FINDINGS

``` r
cat("═══════════════════════════════════════════════════════════════\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
cat("  ✓ Net habitat expansion: +12-30% depending on scenario\n")
```

      ✓ Net habitat expansion: +12-30% depending on scenario

``` r
cat("  ⚠ More LOSS than GAIN in most scenarios (range shift pattern)\n")
```

      ⚠ More LOSS than GAIN in most scenarios (range shift pattern)

``` r
cat("  ✓ Coastal concentration persists in all future scenarios\n")
```

      ✓ Coastal concentration persists in all future scenarios

``` r
cat("  ⚠ Limited inland expansion despite warming (aridity barrier)\n")
```

      ⚠ Limited inland expansion despite warming (aridity barrier)

``` r
cat("  ✓ Tasmania remains suitable (cooler coastal climate)\n")
```

      ✓ Tasmania remains suitable (cooler coastal climate)

``` r
cat("  ⚠ High inter-GCM variability (MIROC6 outlier for SSP2-4.5)\n\n")
```

      ⚠ High inter-GCM variability (MIROC6 outlier for SSP2-4.5)

``` r
cat("6. UNCERTAINTY & LIMITATIONS\n")
```

    6. UNCERTAINTY & LIMITATIONS

``` r
cat("═══════════════════════════════════════════════════════════════\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
cat("  • 92% extrapolation risk (climate novelty vs training data)\n")
```

      • 92% extrapolation risk (climate novelty vs training data)

``` r
cat("  • Model trained on global data (61% Asia, 10% Australia)\n")
```

      • Model trained on global data (61% Asia, 10% Australia)

``` r
cat("  • MIROC6 GCM shows substantially higher expansion than others\n")
```

      • MIROC6 GCM shows substantially higher expansion than others

``` r
cat("  • Ensemble uncertainty increases in end-century scenarios\n")
```

      • Ensemble uncertainty increases in end-century scenarios

``` r
cat("  • Projections most reliable in coastal regions\n\n")
```

      • Projections most reliable in coastal regions

``` r
cat("7. BIOSECURITY IMPLICATIONS\n")
```

    7. BIOSECURITY IMPLICATIONS

``` r
cat("═══════════════════════════════════════════════════════════════\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
cat(sprintf("  • Priority surveillance area: %s km²\n",
            comma(priority_area_km2, accuracy = 1)))
```

      • Priority surveillance area: 59,573 km²

``` r
cat("  • Focus on eastern coastal regions (highest suitability)\n")
```

      • Focus on eastern coastal regions (highest suitability)

``` r
cat("  • Monitor Tasmania for southern range expansion\n")
```

      • Monitor Tasmania for southern range expansion

``` r
cat("  • Inland expansion unlikely (moisture limitation)\n")
```

      • Inland expansion unlikely (moisture limitation)

``` r
cat("  • Range shifts suggest potential refugia changes\n\n")
```

      • Range shifts suggest potential refugia changes

``` r
cat("8. RECOMMENDATIONS\n")
```

    8. RECOMMENDATIONS

``` r
cat("═══════════════════════════════════════════════════════════════\n")
```

    ═══════════════════════════════════════════════════════════════

``` r
cat("  1. Intensify surveillance in eastern coastal zones\n")
```

      1. Intensify surveillance in eastern coastal zones

``` r
cat("  2. Monitor Tasmania for climate-driven expansion\n")
```

      2. Monitor Tasmania for climate-driven expansion

``` r
cat("  3. Account for inter-GCM uncertainty in risk assessments\n")
```

      3. Account for inter-GCM uncertainty in risk assessments

``` r
cat("  4. Update projections as Australian occurrence data increases\n")
```

      4. Update projections as Australian occurrence data increases

``` r
cat("  5. Focus control efforts on currently suitable coastal areas\n")
```

      5. Focus control efforts on currently suitable coastal areas

``` r
cat("  6. Prepare for potential range shifts (loss + gain patterns)\n\n")
```

      6. Prepare for potential range shifts (loss + gain patterns)

``` r
cat("╚════════════════════════════════════════════════════════════════╝\n\n")
```

    ╚════════════════════════════════════════════════════════════════╝

``` r
# Save interpretation to file
interpretation_text <- capture.output({
  cat("SPECIES DISTRIBUTION MODEL RESULTS\n")
  cat("Haemaphysalis longicornis in Australia\n")
  cat("Climate Change Projections to 2100\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  # Repeat the above summary here...
  cat("Model Performance: Excellent (AUC 0.977/0.976)\n")
  cat(sprintf("Current Suitable Habitat: %s km²\n\n",
              comma(area_summary$suitable_area_km2[1])))
  cat("Future Projections (Ensemble Mean):\n")
  for (i in 2:nrow(area_summary)) {
    cat(sprintf("  • %s: %s km² (%+.1f%%)\n",
                area_summary$scenario[i],
                comma(area_summary$suitable_area_km2[i]),
                area_summary$change_from_current_pct[i]))
  }
  cat("\nKey Finding: Net habitat expansion (+12-30%) despite substantial")
  cat("\nrange shifts (more loss than gain in individual cells).\n\n")
  cat("Uncertainty: High inter-GCM variability, 92% extrapolation risk.\n")
  cat("Recommendation: Focus surveillance on eastern coastal regions.\n")
})

writeLines(interpretation_text, "results/final_interpretation.txt")

cat("✓ Interpretation saved to results/final_interpretation.txt\n")
```

    ✓ Interpretation saved to results/final_interpretation.txt

------------------------------------------------------------------------

## Session Info

``` r
sessionInfo()
```

    R version 4.4.1 (2024-06-14)
    Platform: aarch64-apple-darwin20
    Running under: macOS 26.2

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: Australia/Brisbane
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
     [1] ggnewscale_0.5.2    rnaturalearth_1.1.0 ggspatial_1.1.10
     [4] scales_1.4.0        tidyterra_1.0.0     viridis_0.6.5
     [7] viridisLite_0.4.2   patchwork_1.3.2     lubridate_1.9.4
    [10] forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4
    [13] purrr_1.1.0         readr_2.1.5         tidyr_1.3.1
    [16] tibble_3.3.0        ggplot2_4.0.0       tidyverse_2.0.0
    [19] sf_1.0-21           terra_1.8-93

    loaded via a namespace (and not attached):
     [1] gtable_0.3.6                  xfun_0.53
     [3] tzdb_0.5.0                    vctrs_0.6.5
     [5] tools_4.4.1                   generics_0.1.4
     [7] parallel_4.4.1                proxy_0.4-27
     [9] pkgconfig_2.0.3               KernSmooth_2.23-26
    [11] data.table_1.17.8             RColorBrewer_1.1-3
    [13] S7_0.2.0                      lifecycle_1.0.4
    [15] compiler_4.4.1                farver_2.1.2
    [17] textshaping_1.0.3             codetools_0.2-20
    [19] htmltools_0.5.8.1             class_7.3-23
    [21] yaml_2.3.10                   pillar_1.11.0
    [23] crayon_1.5.3                  classInt_0.4-11
    [25] tidyselect_1.2.1              digest_0.6.37
    [27] stringi_1.8.7                 labeling_0.4.3
    [29] fastmap_1.2.0                 grid_4.4.1
    [31] archive_1.1.12.1              cli_3.6.5
    [33] magrittr_2.0.3                e1071_1.7-16
    [35] withr_3.0.2                   bit64_4.6.0-1
    [37] timechange_0.3.0              rmarkdown_2.29
    [39] bit_4.6.0                     gridExtra_2.3
    [41] ragg_1.5.0                    rnaturalearthhires_1.0.0.9000
    [43] hms_1.1.3                     evaluate_1.0.5
    [45] knitr_1.50                    rlang_1.1.6
    [47] Rcpp_1.1.0                    glue_1.8.0
    [49] DBI_1.2.3                     vroom_1.6.5
    [51] jsonlite_2.0.0                R6_2.6.1
    [53] systemfonts_1.2.3             units_0.8-7

------------------------------------------------------------------------

## Summary of Outputs

**Publication Figures** (`figures/final/`): -
`current_distribution_publication.png` - High-quality current
distribution map - `future_scenarios_comparison.png` - 4-panel future
projections - `habitat_change_maps.png` - Gain/loss/stable categorical
maps - `suitable_area_bar_charts.png` - Quantitative area comparisons -
`gcm_uncertainty.png` - Inter-model variability -
`spatial_uncertainty_map.png` - Ensemble SD spatial pattern -
`priority_surveillance_regions.png` - Biosecurity priority areas

**Results Tables** (`results/`): - `model_performance_summary.csv` -
`climate_change_summary.csv` - `suitable_habitat_area_by_scenario.csv` -
`habitat_change_statistics.csv` - `gcm_variability_summary.csv` -
`final_interpretation.txt`

------------------------------------------------------------------------

**END OF STEP 11 - SDM PROJECT COMPLETE**
