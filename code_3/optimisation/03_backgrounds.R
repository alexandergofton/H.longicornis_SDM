source("00_setup.R")
d <- readRDS(file.path(processed_v3, "model_data_allvars.rds")); occ <- d$occ
envM <- rast(file.path(processed_v2, "env_selected_M.tif"))[[1]]
occ_sf <- st_as_sf(occ, coords = c("lon","lat"), crs = 4326)
regs <- c("East_Asia","North_America","Oceania")

# region raster for M (nearest-occurrence region), at 0.25 deg for efficiency, then disaggregated
m6 <- aggregate(!is.na(envM), 6, fun = "max", na.rm = TRUE); m6[m6 == 0] <- NA
xy6 <- as.data.frame(m6, xy = TRUE, na.rm = TRUE)[, 1:2]
reg6 <- match(occ$region[st_nearest_feature(st_as_sf(xy6, coords=c("x","y"), crs=4326), occ_sf)], regs)
regR6 <- rast(m6); regR6[cellFromXY(regR6, as.matrix(xy6))] <- reg6
regR <- resample(regR6, envM, method = "near"); regR <- mask(regR, envM); names(regR) <- "region"
writeRaster(regR, file.path(processed_v3, "M_region_raster.tif"), overwrite = TRUE)

cells <- which(!is.na(values(regR))); creg <- values(regR)[cells]
area_c <- values(cellSize(regR, unit = "km"))[cells]
reg_area <- tapply(area_c, creg, sum); share <- reg_area / sum(reg_area)
cat("Region area share in M:\n"); print(round(share, 3))

n_bg <- 10000
n_reg <- round(n_bg * share)

# ---- target-group bias surface, normalised WITHIN each region
ix <- read_csv(file.path(processed_v2, "ixodidae_points.csv"), show_col_types = FALSE)
z  <- read_csv(file.path(base_dir, "data", "Zhang_et_al_2019_ms_dataset_v1.csv"), show_col_types = FALSE) %>%
  filter(!is.na(lon), !is.na(lat)) %>% select(lon, lat)
tg <- bind_rows(ix, z) %>% distinct(round(lon, 1), round(lat, 1), .keep_all = TRUE) %>% select(lon, lat)
cnt <- rasterize(as.matrix(tg), m6, fun = "count", background = 0)       # 0.25 deg counts
gw  <- focalMat(cnt, 1, "Gauss")                                          # sigma = 1 deg (~100 km)
kde <- focal(cnt, gw, fun = "sum", na.rm = TRUE)
kde <- resample(kde, envM, method = "bilinear")
kv <- values(kde)[cells]; kv[is.na(kv)] <- 0
bias <- numeric(length(cells))
for (r in seq_along(regs)) {
  i <- creg == r
  b <- kv[i] / max(kv[i])
  bias[i] <- pmax(b, 0.02)              # floor so unsurveyed areas remain available
}
bias_r <- rast(regR); bias_r[cells] <- bias; names(bias_r) <- "bias"
writeRaster(bias_r, file.path(processed_v3, "bias_surface_regional.tif"), overwrite = TRUE)
cat("Target-group records used:", nrow(tg), "\n")

sample_bg <- function(prob, seed) {
  set.seed(seed)
  idx <- unlist(lapply(seq_along(regs), function(r) {
    i <- which(creg == r)
    sample(i, n_reg[r], prob = prob[i] * area_c[i])   # area-weight to correct lat/lon cell size
  }))
  ce <- cells[idx]; xy <- xyFromCell(regR, ce)
  out <- tibble(lon = xy[,1], lat = xy[,2], region = regs[creg[idx]])
  bind_cols(out, terra::extract(env_all, xy)) %>% drop_na(all_of(all_vars))
}
bg_uniform  <- sample_bg(rep(1, length(cells)), 101)
bg_regbias  <- sample_bg(bias, 102)
bg_eval     <- sample_bg(rep(1, length(cells)), 999)   # independent evaluation background (uniform)
bg_current  <- d$bg
for (nm in c("bg_uniform","bg_regbias","bg_eval","bg_current")) { cat(nm, ":\n"); print(table(get(nm)$region)) }

# ---- spatial-block folds: 5-degree cells randomly assigned to 5 folds
set.seed(2024)
blk_id <- function(lon, lat) paste(floor(lon/5), floor(lat/5))
all_blk <- unique(c(blk_id(occ$lon, occ$lat), blk_id(bg_eval$lon, bg_eval$lat), blk_id(bg_uniform$lon, bg_uniform$lat),
                    blk_id(bg_regbias$lon, bg_regbias$lat), blk_id(bg_current$lon, bg_current$lat)))
blk_fold <- setNames(sample(rep(1:5, length.out = length(all_blk))), all_blk)
addf <- function(x) x %>% mutate(block = blk_fold[blk_id(lon, lat)])
saveRDS(list(occ = addf(occ), bg_current = addf(bg_current), bg_uniform = addf(bg_uniform),
             bg_regbias = addf(bg_regbias), bg_eval = addf(bg_eval), regs = regs, n_reg = n_reg),
        file.path(processed_v3, "modelling_sets.rds"))
print(table(addf(occ)$block, occ$region))
