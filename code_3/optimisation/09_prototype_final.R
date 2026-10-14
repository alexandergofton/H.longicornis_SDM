source("00_setup.R"); source("04_cv_functions.R")
vars <- c("bio1","bio2","bio13","bio17","aridity_idx")
m <- fit_mx(S$occ, S$bg_regbias, vars, "LQH", 8)
print(m@results[grep("contribution|permutation|Training.AUC|Regularized.training.gain", rownames(m@results)), , drop = FALSE])
lam <- m@lambdas; cat("n nonzero features:", sum(sapply(strsplit(lam[1:(length(lam)-4)], ","), function(z) as.numeric(z[2]) != 0)), "\n")
env_au <- crop(env_all[[vars]], ext(110, 180, -48, -8))
pr <- predict(env_au, m, args = "outputformat=cloglog", na.rm = TRUE)
po <- pred_mx(m, S$occ, vars); thr10 <- quantile(po, 0.1)
pe <- pred_mx(m, S$bg_regbias, vars)
ev <- dismo::evaluate(po, pe); thr_ss <- ev@t[which.max(ev@TPR + ev@TNR)]
cat("p10 thr:", thr10, " maxSSS thr:", thr_ss, "\n")
print(tibble(region = S$occ$region, p = po) %>% group_by(region) %>% summarise(mean = mean(p), pct_p10 = mean(p >= thr10)*100, pct_ss = mean(p >= thr_ss)*100))
boxes <- list(PNG_lowlands = ext(141, 147, -9, -3), Java = ext(105, 115, -8.5, -6), FarNorthQLD = ext(142, 146, -18, -13),
              SE_QLD_NSW = ext(150, 154, -34, -27), Victoria = ext(143, 149, -39, -36.5), NZ_North = ext(172.5, 178.5, -41.5, -34.5),
              Tasmania = ext(144.5, 148.5, -43.6, -40.6), SW_WA = ext(114.5, 118.5, -35.2, -31.5), CentralDesert = ext(125, 140, -30, -20))
for (nm in names(boxes)) { x <- crop(pr, boxes[[nm]])
  cat(sprintf("%-14s mean = %.3f ; %%>=p10 = %.1f\n", nm, global(x, mean, na.rm=TRUE)[[1]], 100*global(x >= thr10, mean, na.rm=TRUE)[[1]])) }
writeRaster(pr, "/tmp/hl_sdm/proto_pred.tif", overwrite = TRUE)
png("/tmp/hl_sdm/figs/proto.png", 1400, 1000, res = 130); plot(pr, col = hcl.colors(100, "viridis"), range = c(0,1)); points(S$occ$lon, S$occ$lat, pch = 16, cex = .4, col = "red"); dev.off()
