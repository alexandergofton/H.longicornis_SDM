source("00_setup.R"); source("04_cv_functions.R")
cfg <- list(H = list(v = c("bio1","bio12","bio15","aridity_idx"), fc = "LQ", rm = 3),
            C = list(v = c("bio1","bio2","bio12","bio15","aridity_idx"), fc = "LQ", rm = 2),
            G = list(v = c("bio1","bio2","bio13","bio17","aridity_idx"), fc = "LQH", rm = 4))
fit_one <- function(cf) function(p, a) { m <- fit_mx(p, a, cf$v, cf$fc, cf$rm); function(x) pred_mx(m, x, cf$v) }
ens <- function(p, a) { fs <- lapply(cfg, function(cf) fit_one(cf)(p, a)); function(x) rowMeans(sapply(fs, function(f) f(x))) }
B <- S$bg_regbias
out <- map_dfr(c(names(cfg), "consensus"), function(nm) {
  ff <- if (nm == "consensus") ens else fit_one(cfg[[nm]])
  cv <- bind_rows(run_cv(ff, B, "loro"), run_cv(ff, B, "block"))
  f <- ff(S$occ, B); po <- f(S$occ); thr <- quantile(po, .1)
  env_au <- crop(env_all, ext(110, 180, -48, -8))
  boxes <- list(SW_WA = ext(114.5, 118.5, -35.2, -31.5), Vic = ext(143, 149, -39, -36.5), Tas = ext(144.5, 148.5, -43.6, -40.6), FNQ = ext(142,146,-18,-13), SA_SE = ext(139.5,141,-38.1,-36))
  bx <- sapply(boxes, function(b) { e <- as.data.frame(crop(env_au, b), na.rm = TRUE); mean(f(e) >= thr) * 100 })
  wa <- S$occ %>% filter(lon < 125, region == "Oceania"); vic <- S$occ %>% filter(lon > 140, lon < 150, lat < -36.5)
  bind_cols(tibble(model = nm), summarise_cv(cv) %>% dplyr::select(loro_auc, loro_cbi, loro_or10, oce_auc, oce_cbi, blk_auc, blk_cbi, blk_or10),
            as_tibble(as.list(round(bx, 1))), tibble(WA_recs_suit = mean(f(wa) >= thr), Vic_recs_suit = mean(f(vic) >= thr), oce_recs_suit = mean(po[S$occ$region=="Oceania"] >= thr)))
})
print(as.data.frame(out %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
