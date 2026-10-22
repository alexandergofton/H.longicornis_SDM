# Prepare future-climate predictor stacks for Australia/New Zealand (CMIP6, WorldClim 2.1 downscaling).
# The same code is embedded in Hlongicornis_future_SDM.qmd (chunk `future-climate`); this script lets the
# slow download step (~30-60 min) be run once outside the render. Run from code_4/: Rscript prep_future_climate.R
suppressPackageStartupMessages({ library(terra); library(envirem); library(parallel); library(dplyr); library(tidyr) })
base_dir <- normalizePath("..")
fut_dir  <- file.path(base_dir, "data", "future_climate_cmip6")
dir.create(fut_dir, showWarnings = FALSE)
n_cores  <- min(6, detectCores() - 1)

gcm_core <- c("ACCESS-CM2", "ACCESS-ESM1-5", "CMCC-ESM2", "CNRM-ESM2-1", "EC-Earth3-Veg",
              "GISS-E2-1-G", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")
scen <- bind_rows(
  expand_grid(gcm = gcm_core, ssp = c("126", "245", "370", "585"),
              period = c("2021-2040", "2041-2060", "2061-2080", "2081-2100"), set = "core"),
  expand_grid(gcm = c("UKESM1-0-LL", "CNRM-CM6-1"), ssp = c("245", "585"),
              period = c("2041-2060", "2081-2100"), set = "v2_only")) |>
  mutate(id = sprintf("%s_ssp%s_%s", gcm, ssp, period), file = file.path(fut_dir, paste0(id, ".tif")))

cur_suit <- rast(file.path(base_dir, "processed_data_3", "suitability_current_aunz.tif"))   # Aus/NZ land grid
cells    <- which(!is.na(values(cur_suit)))
env_all  <- rast(file.path(base_dir, "processed_data_2", "env_all_global.tif"))

# ---- begin: shared with qmd chunk `future-climate` ----
# Bioclimatic variables (WorldClim/dismo definitions) from cells x 12 monthly matrices
bioclim_m <- function(pr, tn, tx) {
  p1 <- pr + 1
  q  <- sapply(1:12, function(m) rowSums(pr[, ((m - 1 + 0:2) %% 12) + 1]))        # rolling 3-month totals
  cbind(bio1  = rowMeans((tn + tx) / 2),
        bio2  = rowMeans(tx - tn),
        bio12 = rowSums(pr),
        bio13 = do.call(pmax, as.data.frame(pr)),
        bio15 = 100 * sqrt(rowSums((p1 - rowMeans(p1))^2) / 11) / rowMeans(p1),
        bio17 = do.call(pmin, as.data.frame(q)))
}
# Thornthwaite aridity index (envirem::aridityIndexThornthwaite) with Hargreaves PET (envirem::monthlyPET);
# extraterrestrial radiation by latitude and month as in ENVIREM (envirem:::calcSolRad, year = 40)
lat_rows <- yFromRow(cur_suit, 1:nrow(cur_suit))
RA <- sapply(1:12, function(m) sapply(lat_rows, function(y) envirem:::calcSolRad(year = 40, lat = y, month = m)))
RA_m <- RA[rowFromCell(cur_suit, cells), ]
aridity_m <- function(pr, tn, tx) {
  pet <- 0.0023 * (RA_m * 30) * ((tn + tx) / 2 + 17.8) * sqrt(tx - tn)
  pet[pet < 0] <- 0
  dry <- pr < pet
  n <- rowSums(pet * dry); d <- rowSums((pet - pr) * dry)
  ifelse(n == 0, 0, 100 * d / n)
}
read_monthly <- function(f) values(crop(rast(f), cur_suit))[cells, ]
wc_dir <- file.path(base_dir, "data", "worldclim", "climate", "wc2.1_2.5m")
cur_m  <- lapply(c(pr = "prec", tn = "tmin", tx = "tmax"), function(v)
  read_monthly(file.path(wc_dir, sprintf("wc2.1_2.5m_%s_%02d.tif", v, 1:12))))
ai_wc21_cur <- aridity_m(cur_m$pr, cur_m$tn, cur_m$tx)                   # current, same formulation
ai_envirem  <- values(crop(env_all[["aridity_idx"]], cur_suit))[cells]  # ENVIREM baseline used in training
fut_vars <- c("bio1", "bio2", "bio12", "bio13", "bio15", "bio17", "aridity_idx")
prep_one <- function(i) {
  s <- scen[i, ]
  if (file.exists(s$file)) return("cached")
  url <- sprintf("/vsicurl/https://geodata.ucdavis.edu/cmip6/2.5m/%s/ssp%s/wc2.1_2.5m_%%s_%s_ssp%s_%s.tif",
                 s$gcm, s$ssp, s$gcm, s$ssp, s$period)
  tryCatch({
    m <- lapply(c(pr = "prec", tn = "tmin", tx = "tmax"), function(v) read_monthly(sprintf(url, v)))
    ai_in <- pmin(pmax(ai_envirem + (aridity_m(m$pr, m$tn, m$tx) - ai_wc21_cur), 0), 100)   # anomaly on training baseline
    x <- cbind(bioclim_m(m$pr, m$tn, m$tx), aridity_idx = ai_in)
    r <- rast(cur_suit, nlyrs = length(fut_vars)); names(r) <- fut_vars
    v <- matrix(NA_real_, ncell(r), length(fut_vars)); v[cells, ] <- x[, fut_vars]; values(r) <- v
    writeRaster(r, s$file, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
    "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
}
# Socket (PSOCK) workers: forked processes cannot use libcurl safely on macOS
todo <- which(!file.exists(scen$file))
if (length(todo)) {
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, { suppressPackageStartupMessages(library(terra))
    Sys.setenv(GDAL_HTTP_MAX_RETRY = "5", GDAL_HTTP_RETRY_DELAY = "10", GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
               CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif"); NULL })
  cs_w <- wrap(cur_suit)
  clusterExport(cl, c("scen", "cells", "RA_m", "ai_wc21_cur", "ai_envirem", "fut_vars", "bioclim_m", "aridity_m", "cs_w"),
                envir = environment())
  clusterEvalQ(cl, { cur_suit <- unwrap(cs_w); read_monthly <- function(f) values(crop(rast(f), cur_suit))[cells, ]; NULL })
  prep_res <- parLapplyLB(cl, todo, prep_one)
  stopCluster(cl)
  print(table(unlist(prep_res)))
}
# ---- end: shared with qmd chunk ----
