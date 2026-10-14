source("00_setup.R"); source("04_cv_functions.R")
library(parallel)
sets <- list(
  S0_current   = c("bio2","bio3","bio13","aridity_idx","vpd_annual"),
  S1_limits    = c("bio5","bio6","bio13","bio17","aridity_idx"),
  S2_means     = c("bio1","bio2","bio12","bio15","aridity_idx"),
  S3_limits_dr = c("bio2","bio5","bio6","bio13","aridity_idx"),
  S4_vpd       = c("vpd_annual","bio4","bio13","bio17","aridity_idx"),
  S5_seasonal  = c("bio5","bio6","bio12","bio15","aridity_idx"))
B <- S$bg_uniform
coll <- imap_dfr(sets, function(v, nm) {
  cm <- abs(cor(B[, v], method = "spearman")); diag(cm) <- 0
  tibble(set = nm, vars = paste(v, collapse = ", "), max_abs_r = max(cm), max_vif = max(usdm::vif(as.data.frame(B[, v]))$VIF))
})
print(coll); write_csv(coll, file.path(outputs_v3, "opt_02_varset_collinearity.csv"))

grid <- expand_grid(set = names(sets), fc = c("LQ","LQH","LQHP"), rm = c(1, 2, 4), bg = c("bg_uniform","bg_regbias"))
cl <- makeCluster(7)
clusterEvalQ(cl, { options(java.parameters = "-Xmx1500m"); setwd(getwd()); source("00_setup.R"); source("04_cv_functions.R"); NULL })
clusterExport(cl, c("sets","grid"))
res <- parLapplyLB(cl, seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]; ff <- mx_fitter(sets[[g$set]], g$fc, g$rm); Bt <- S[[g$bg]]
  cv <- bind_rows(run_cv(ff, Bt, "loro"), run_cv(ff, Bt, "block"))
  bind_cols(g, summarise_cv(cv))
})
stopCluster(cl)
res <- bind_rows(res); write_csv(res, file.path(outputs_v3, "opt_03_varset_complexity_grid.csv"))
res %>% group_by(set, bg) %>% summarise(across(c(loro_auc, loro_cbi, oce_auc, oce_cbi, blk_auc, blk_cbi, loro_or10), mean), .groups = "drop") %>%
  arrange(desc(loro_auc)) %>% mutate(across(where(is.numeric), ~round(.x,3))) %>% print(n = 50)
res %>% arrange(desc(loro_auc + blk_auc)) %>% head(20) %>% mutate(across(where(is.numeric), ~round(.x,3))) %>% print(width = 200)
