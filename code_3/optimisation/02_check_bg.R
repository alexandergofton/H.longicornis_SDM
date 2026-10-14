source("00_setup.R")
d <- readRDS(file.path(processed_v3, "model_data_allvars.rds")); bg <- d$bg; occ <- d$occ
envM <- rast(file.path(processed_v2, "env_selected_M.tif"))[[1]]
m4 <- aggregate(!is.na(envM), 4, fun = "max", na.rm = TRUE); m4[m4 == 0] <- NA
xy <- as.data.frame(m4, xy = TRUE, na.rm = TRUE)[, 1:2]
ca <- extract(cellSize(m4, unit = "km"), as.matrix(xy))[, 1]
occ_sf <- st_as_sf(occ, coords=c("lon","lat"), crs=4326)
reg <- occ$region[st_nearest_feature(st_as_sf(xy, coords=c("x","y"), crs=4326), occ_sf)]
print(tapply(ca, reg, sum))
ix <- read_csv(file.path(processed_v2, "ixodidae_points.csv"), show_col_types=FALSE)
ix$region <- occ$region[st_nearest_feature(st_as_sf(ix, coords=c("lon","lat"), crs=4326), occ_sf)]
cat("Ixodidae target-group points by region:\n"); print(table(ix$region))
z <- read_csv(file.path(base_dir,"data","Zhang_et_al_2019_ms_dataset_v1.csv"), show_col_types=FALSE)
cat("Zhang et al. tick records:", nrow(z), " species:", n_distinct(z$tick_sp), "\n")
print(summary(z[,c("lon","lat")]))
