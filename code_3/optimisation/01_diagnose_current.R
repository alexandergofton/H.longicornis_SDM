source("00_setup.R")
d <- build_data(); occ <- d$occ; bg <- d$bg
saveRDS(d, file.path(processed_v3, "model_data_allvars.rds"))
cat("occ:", nrow(occ), " bg:", nrow(bg), "\n"); print(table(occ$region)); print(table(bg$region))

# ---- collinearity on background (study area) vs occurrences
v0 <- c("bio2","bio3","bio13","aridity_idx","vpd_annual")
cat("\nCorrelation of current set on BACKGROUND:\n"); print(round(cor(bg[, v0]), 2))
cat("\nVIF (background):\n"); print(usdm::vif(as.data.frame(bg[, v0])))
cat("\nFull |r| on background for all candidates:\n")
cm <- cor(bg[, all_vars], method = "spearman"); print(round(cm, 2))

# ---- refit current model (LQHP, rm=2) with partitions
regs <- c("East_Asia","North_America","Oceania")
grp_o <- match(occ$region, regs); grp_b <- match(bg$region, regs)
occs <- occ %>% select(lon, lat, all_of(v0)) %>% as.data.frame()
bgs  <- bg  %>% select(lon, lat, all_of(v0)) %>% as.data.frame()

run <- function(part, vbg, ug = NULL) ENMevaluate(occs = occs, bg = bgs, algorithm = "maxent.jar",
  tune.args = list(fc = "LQHP", rm = 2), partitions = part, user.grp = ug, doClamp = TRUE,
  other.settings = list(validation.bg = vbg), quiet = TRUE)
e_block_full <- run("block", "full")
e_block_part <- run("block", "partition")
e_loro_full  <- run("user", "full", list(occs.grp = grp_o, bg.grp = grp_b))
e_loro_part  <- run("user", "partition", list(occs.grp = grp_o, bg.grp = grp_b))
show <- function(e, lab) { r <- e@results.partitions; r$setting <- lab; r }
res <- bind_rows(show(e_block_full,"block/full bg"), show(e_block_part,"block/partition bg"),
                 show(e_loro_full,"LORO/full bg"), show(e_loro_part,"LORO/partition bg"))
res$fold_name <- ifelse(grepl("LORO", res$setting), regs[res$fold], paste0("block", res$fold))
print(res %>% select(setting, fold_name, auc.val, cbi.val, or.10p, or.mtp))
write_csv(res, file.path(outputs_v3, "diag_current_model_partitions.csv"))

# ---- full model predictions at Oceanian records & tropical sanity checks
m <- e_block_full@models[[1]]
p_occ <- dismo::predict(m, occs[, v0]); p_bg <- dismo::predict(m, bgs[, v0])
thr <- 0.2138
cat("\nShare of records predicted suitable (>= MaxSS 0.214) by region:\n")
print(tibble(region = occ$region, p = p_occ) %>% group_by(region) %>% summarise(n = n(), mean_suit = mean(p), pct_suitable = mean(p >= thr)*100))
env_sel <- env_all[[v0]]
boxes <- list(PNG_lowlands = ext(141, 147, -9, -3), Java = ext(105, 115, -8.5, -6), FarNorthQLD_lowland = ext(142, 146, -18, -13),
              SE_QLD_NSW = ext(150, 154, -34, -27), Victoria = ext(143, 149, -39, -36.5), NZ_North = ext(172.5, 178.5, -41.5, -34.5))
for (nm in names(boxes)) {
  pr <- predict(crop(env_sel, boxes[[nm]]), m, type = "cloglog")
  cat(sprintf("%-22s mean suit = %.3f ; %% cells >= thr = %.1f\n", nm, global(pr, mean, na.rm=TRUE)[[1]], 100*global(pr >= thr, mean, na.rm=TRUE)[[1]]))
}
