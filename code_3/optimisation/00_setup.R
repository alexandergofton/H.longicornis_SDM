# Shared setup for v3 optimisation scripts
options(java.parameters = "-Xmx12g")
suppressPackageStartupMessages({
  library(tidyverse); library(terra); library(sf); library(ENMeval)
})
base_dir     <- "~/Library/CloudStorage/OneDrive-CSIRO/OneDrive - Docs/01_Projects/SDM/Hlongicornis_SDM"
processed_v2 <- file.path(base_dir, "processed_data_2")
processed_v3 <- file.path(base_dir, "processed_data_3")
outputs_v3   <- file.path(base_dir, "outputs_3")
opt_dir      <- file.path(base_dir, "code_3", "optimisation")
terraOptions(memfrac = 0.6, tempdir = file.path(base_dir, "temp"))

env_all <- rast(file.path(processed_v2, "env_all_global.tif"))
all_vars <- c(paste0("bio", c(1:4, 5:7, 10:17)), "pet_idx", "aridity_idx", "cmi_idx", "vpd_annual")

build_data <- function() {
  occ <- read_csv(file.path(processed_v2, "occurrences_with_climate.csv"), show_col_types = FALSE) %>%
    mutate(region = ifelse(region == "Other", "East_Asia", region))   # Himalayan record -> native range
  bg <- read_csv(file.path(processed_v2, "background_points.csv"), show_col_types = FALSE) %>% select(lon, lat)
  bg_env <- terra::extract(env_all, as.matrix(bg[, c("lon", "lat")]))
  bg <- bind_cols(bg, bg_env) %>% drop_na(all_of(all_vars))
  # assign each background point to the region of its nearest occurrence (M = buffers around occurrences)
  occ_sf <- st_as_sf(occ, coords = c("lon", "lat"), crs = 4326)
  bg_sf  <- st_as_sf(bg,  coords = c("lon", "lat"), crs = 4326)
  bg$region <- occ$region[st_nearest_feature(bg_sf, occ_sf)]
  list(occ = occ, bg = bg)
}
