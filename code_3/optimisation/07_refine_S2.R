source("00_setup.R"); source("04_cv_functions.R")
library(parallel)
sets <- list(
  S2_means   = c("bio1","bio2","bio12","bio15","aridity_idx"),
  S2a_noBio2 = c("bio1","bio12","bio15","aridity_idx"),
  S2b_bio11  = c("bio11","bio2","bio12","bio15","aridity_idx"),
  S2c_bio10  = c("bio10","bio2","bio12","bio15","aridity_idx"),
  S2d_bio13  = c("bio1","bio2","bio13","bio15","aridity_idx"),
  S2e_cmi    = c("bio1","bio2","bio12","bio15","cmi_idx"),
  S2f_bio17  = c("bio1","bio2","bio17","bio13","aridity_idx"))
B <- S$bg_uniform
coll <- imap_dfr(sets, function(v, nm) { cm <- abs(cor(B[, v], method = "spearman")); diag(cm) <- 0
  tibble(set = nm, max_abs_r = max(cm), max_vif = max(usdm::vif(as.data.frame(B[, v]))$VIF)) })
print(coll)
grid <- expand_grid(set = names(sets), fc = c("L","LQ","LQH","LQHP"), rm = c(1, 2, 4, 6, 8), bg = c("bg_uniform","bg_regbias"))
cl <- makeCluster(7)
clusterEvalQ(cl, { options(java.parameters = "-Xmx1500m"); source("00_setup.R"); source("04_cv_functions.R"); NULL })
clusterExport(cl, c("sets","grid"))
res <- bind_rows(parLapplyLB(cl, seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]; ff <- mx_fitter(sets[[g$set]], g$fc, g$rm); Bt <- S[[g$bg]]
  cv <- bind_rows(run_cv(ff, Bt, "loro"), run_cv(ff, Bt, "block"))
  bind_cols(g, summarise_cv(cv))
}))
stopCluster(cl)
write_csv(res, file.path(outputs_v3, "opt_04_S2_refinement_grid.csv"))
res %>% group_by(set, bg) %>% summarise(across(c(loro_auc, loro_cbi, loro_or10, oce_auc, oce_cbi, blk_auc, blk_cbi, blk_or10), mean), .groups = "drop") %>%
  arrange(desc(loro_auc)) %>% mutate(across(where(is.numeric), ~round(.x,3))) %>% print(n = 50)
res %>% mutate(score = (loro_auc + blk_auc)/2 + (loro_cbi + blk_cbi)/4 - loro_or10/2) %>% arrange(desc(score)) %>% head(25) %>%
  mutate(across(where(is.numeric), ~round(.x,3))) %>% select(set, fc, rm, bg, score, loro_auc, loro_cbi, loro_or10, oce_auc, oce_cbi, na_auc, ea_auc, blk_auc, blk_cbi, blk_or10) %>% print(n=25, width = 250)
