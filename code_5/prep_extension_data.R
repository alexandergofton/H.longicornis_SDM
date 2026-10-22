# Data preparation for the code_5 extensions (run once from code_5/: Rscript prep_extension_data.R).
# The same code is embedded in Hlongicornis_impact_extensions.qmd (chunk `prep-data`), where it is skipped
# when the outputs exist. Writes to data/ (raw inputs) and data/future_climate_cmip6_monthly/ (derived).
suppressPackageStartupMessages({ library(terra); library(geodata); library(parallel); library(dplyr); library(tidyr) })
base_dir <- normalizePath("..")
dat_dir  <- file.path(base_dir, "data")
mon_dir  <- file.path(dat_dir, "future_climate_cmip6_monthly")
pac_dir  <- file.path(dat_dir, "pacific")
for (d in c(mon_dir, pac_dir, file.path(dat_dir, "landcover_worldcover"), file.path(dat_dir, "population_gpw4"),
            file.path(dat_dir, "elevation_30s"))) dir.create(d, showWarnings = FALSE)
Sys.setenv(GDAL_HTTP_MAX_RETRY = "5", GDAL_HTTP_RETRY_DELAY = "10", GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
           CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif")
n_cores <- min(6, detectCores() - 1)

# ---- begin: shared with qmd chunk `prep-data` ----
cur_suit <- rast(file.path(base_dir, "processed_data_3", "suitability_current_aunz.tif"))   # Aus/NZ 2.5' land grid
cells    <- which(!is.na(values(cur_suit)))
env_all  <- rast(file.path(base_dir, "processed_data_2", "env_all_global.tif"))
wc_dir   <- file.path(dat_dir, "worldclim", "climate", "wc2.1_2.5m")

# (a) Static layers -------------------------------------------------------------------------------
# Current monthly water-vapour pressure (WorldClim 2.1; kPa) for saturation deficit
if (!file.exists(file.path(wc_dir, "wc2.1_2.5m_vapr_12.tif")))
  invisible(geodata::worldclim_global(var = "vapr", res = 2.5, path = file.path(dat_dir, "worldclim")))
# Human population density 2020 (GPW v4, 2.5')
pop_file <- file.path(dat_dir, "population_gpw4", "gpw_v4_population_density_2020_2.5m.tif")
if (!file.exists(pop_file))
  writeRaster(rast("/vsicurl/https://geodata.ucdavis.edu/geodata/pop/gpw_v4_population_density_rev11_2020_2.5m.tif"),
              pop_file, overwrite = TRUE)
# Population counts at 2.5', aggregated from 30" density x 30" cell area (handles coastlines correctly)
popc_aunz <- file.path(dat_dir, "population_gpw4", "gpw_v4_population_count_2020_2.5m_aunz.tif")
popc_pac  <- file.path(dat_dir, "population_gpw4", "gpw_v4_population_count_2020_2.5m_pacific.tif")
if (!file.exists(popc_aunz) || !file.exists(popc_pac)) {
  r30 <- rast("/vsicurl/https://geodata.ucdavis.edu/geodata/pop/gpw_v4_population_density_rev11_2020_30s.tif")
  cnt <- function(e, dx = 0) { a <- crop(r30, e) * 1; a <- a * cellSize(a, unit = "km"); if (dx != 0) a <- shift(a, dx = dx)
    aggregate(a, 5, sum, na.rm = TRUE) }
  writeRaster(cnt(ext(110, 180, -48, -8)), popc_aunz, overwrite = TRUE)
  writeRaster(merge(cnt(ext(140, 180, -30, 0)), cnt(ext(-180, -130, -30, 0), dx = 360)), popc_pac, overwrite = TRUE)
}
# ESA WorldCover fractional cover at 30" (trees, grassland, cropland), Australia/NZ window, stored as percent
for (v in c("trees", "grassland", "cropland")) {
  f <- file.path(dat_dir, "landcover_worldcover", sprintf("WorldCover_%s_30s_aunz.tif", v))
  if (!file.exists(f)) {
    r <- crop(rast(sprintf("/vsicurl/https://geodata.ucdavis.edu/geodata/landuse/WorldCover_%s_30s.tif", v)),
              ext(110, 180, -48, -8))
    writeRaster(round(r * 100), f, datatype = "INT1U", overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
  }
}
# 30" elevation (SRTM-based) for Australia and New Zealand
elev30_file <- file.path(dat_dir, "elevation_30s", "elevation_30s_aunz.tif")
if (!file.exists(elev30_file)) {
  e1 <- geodata::elevation_30s(country = "AUS", path = file.path(dat_dir, "elevation_30s"))
  e2 <- geodata::elevation_30s(country = "NZL", path = file.path(dat_dir, "elevation_30s"))
  writeRaster(merge(e1, e2), elev30_file, overwrite = TRUE)
}

# (b) Pacific/PNG grid: 140E-130W, 30S-0, in 0-360 longitudes --------------------------------------
rd_pac <- function(r) {
  a <- crop(r, ext(140, 180, -30, 0)); b <- shift(crop(r, ext(-180, -130, -30, 0)), dx = 360)
  merge(a, b)
}
pac_tmpl_file <- file.path(pac_dir, "pacific_env_current.tif")
if (!file.exists(pac_tmpl_file)) writeRaster(rd_pac(env_all), pac_tmpl_file, overwrite = TRUE)
pac_env  <- rast(pac_tmpl_file)
pac_cells <- which(!is.na(values(pac_env[["bio1"]])) & !is.na(values(pac_env[["aridity_idx"]])))

# (c) Future monthly climate --------------------------------------------------------------------------
gcm_core <- c("ACCESS-CM2", "ACCESS-ESM1-5", "CMCC-ESM2", "CNRM-ESM2-1", "EC-Earth3-Veg",
              "GISS-E2-1-G", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")
scen_m <- bind_rows(expand_grid(gcm = gcm_core, ssp = c("245", "370"), period = c("2041-2060", "2081-2100")),
                    expand_grid(gcm = gcm_core, ssp = "245", period = "2021-2040")) |>
  mutate(id = sprintf("%s_ssp%s_%s", gcm, ssp, period),
         f_aunz = file.path(mon_dir, paste0(id, "_tmin_tmax_aunz.tif")),
         f_pac  = file.path(mon_dir, paste0(id, "_monthly_pacific.tif")))
prep_month <- function(i) {
  s <- scen_m[i, ]
  if (file.exists(s$f_aunz) && file.exists(s$f_pac)) return("cached")
  url <- sprintf("/vsicurl/https://geodata.ucdavis.edu/cmip6/2.5m/%s/ssp%s/wc2.1_2.5m_%%s_%s_ssp%s_%s.tif",
                 s$gcm, s$ssp, s$gcm, s$ssp, s$period)
  tryCatch({
    au <- list(); pc <- list()
    for (v in c("tmin", "tmax", "prec")) {                      # one read of the 48S-0 band per variable
      band <- crop(rast(sprintf(url, v)), ext(-180, 180, -48, 0)) * 1
      if (v != "prec") au[[v]] <- mask(crop(band, cur_suit), cur_suit)
      pc[[v]] <- rd_pac(band); rm(band); gc()
    }
    writeRaster(round(rast(au) * 10), s$f_aunz, datatype = "INT2S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
    writeRaster(round(rast(pc) * 10), s$f_pac, datatype = "INT2S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
    "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
}
todo <- which(!(file.exists(scen_m$f_aunz) & file.exists(scen_m$f_pac)))
if (length(todo)) {
  cl <- makeCluster(min(4, n_cores))                                 # ~0.5 GB per worker while a band is in memory
  clusterEvalQ(cl, { suppressPackageStartupMessages(library(terra))
    Sys.setenv(GDAL_HTTP_MAX_RETRY = "5", GDAL_HTTP_RETRY_DELAY = "10", GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
               CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif"); NULL })
  cs_w <- wrap(cur_suit)
  clusterExport(cl, c("scen_m", "cs_w", "rd_pac"), envir = environment())
  clusterEvalQ(cl, { cur_suit <- unwrap(cs_w); NULL })
  prep_res <- parLapplyLB(cl, todo, prep_month)
  stopCluster(cl)
  print(table(unlist(prep_res)))
}
# ---- end: shared with qmd chunk ----
