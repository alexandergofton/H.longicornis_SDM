source("00_setup.R"); source("04_cv_functions.R")
v0 <- c("bio2","bio3","bio13","aridity_idx","vpd_annual")
grid <- expand_grid(bg = c("bg_current","bg_uniform","bg_regbias"),
                    spec = list(c("LQHP","2"), c("LQ","2"), c("LQH","3")))
res <- pmap_dfr(grid, function(bg, spec) {
  ff <- mx_fitter(v0, spec[1], as.numeric(spec[2])); B <- S[[bg]]
  cv <- bind_rows(run_cv(ff, B, "loro"), run_cv(ff, B, "block"))
  # full-model regional discrimination (in-sample) & NA sensitivity
  f <- ff(S$occ, B); po <- f(S$occ); pe <- f(S$bg_eval); thr <- quantile(po, 0.1)
  insamp <- map_dbl(S$regs, ~ auc_fun(po[S$occ$region == .x], pe[S$bg_eval$region == .x]))
  out <- bind_cols(tibble(bg = bg, fc = spec[1], rm = spec[2]), summarise_cv(cv),
            tibble(full_auc_EA = insamp[1], full_auc_NA = insamp[2], full_auc_OC = insamp[3],
                   pct_NA_above_p10 = mean(po[S$occ$region=="North_America"] >= thr)*100))
  print(out %>% mutate(across(where(is.numeric), ~round(.x, 3))) %>% as.data.frame()); out
})
write_csv(res, file.path(outputs_v3, "opt_01_background_comparison.csv"))
