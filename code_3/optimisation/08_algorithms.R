source("00_setup.R"); source("04_cv_functions.R")
suppressPackageStartupMessages({ library(mgcv); library(gbm); library(ranger) })
vars <- c("bio1","bio2","bio13","bio17","aridity_idx")
prep <- function(p, a) { d <- bind_rows(mutate(p[, vars], y = 1), mutate(a[, vars], y = 0)) %>% as.data.frame()
  d$w <- ifelse(d$y == 1, 1, sum(d$y == 1) / sum(d$y == 0)); d }
fitters <- list(
  maxent = function(p, a) { m <- fit_mx(p, a, vars, "LQH", 8); function(x) pred_mx(m, x, vars) },
  glm = function(p, a) { d <- prep(p, a); sc <- scale(d[, vars]); ctr <- attr(sc, "scaled:center"); s <- attr(sc, "scaled:scale")
    dd <- as.data.frame(sc); dd$y <- d$y
    f <- as.formula(paste("y ~", paste(c(vars, paste0("I(", vars, "^2)")), collapse = " + ")))
    m <- suppressWarnings(glm(f, data = dd, family = binomial, weights = d$w))
    function(x) { nx <- as.data.frame(scale(as.data.frame(x[, vars]), ctr, s)); as.numeric(predict(m, nx, type = "response")) } },
  gam = function(p, a) { d <- prep(p, a)
    f <- as.formula(paste("y ~", paste0("s(", vars, ", k = 5)", collapse = " + ")))
    m <- gam(f, data = d, family = binomial, weights = w, method = "REML")
    function(x) as.numeric(predict(m, as.data.frame(x[, vars]), type = "response")) },
  brt = function(p, a) { d <- prep(p, a)
    m <- gbm(y ~ ., data = d[, c(vars, "y")], weights = d$w, distribution = "bernoulli", n.trees = 2000,
             interaction.depth = 2, shrinkage = 0.005, bag.fraction = 0.75, n.minobsinnode = 10, verbose = FALSE)
    function(x) as.numeric(predict(m, as.data.frame(x[, vars]), n.trees = 2000, type = "response")) },
  rf_down = function(p, a) { d <- prep(p, a); d$y <- factor(d$y); np <- sum(d$y == "1")
    m <- ranger(y ~ ., data = d[, c(vars, "y")], num.trees = 1000, probability = TRUE, replace = TRUE,
                sample.fraction = c(np, np) / nrow(d), min.node.size = 5)
    function(x) as.numeric(predict(m, as.data.frame(x[, vars]))$predictions[, "1"]) })
# ensemble: convert each model to percentile of its own training-background predictions, then average
ens_fitter <- function(members) function(p, a) {
  fs <- lapply(fitters[members], function(f) f(p, a))
  refs <- lapply(fs, function(f) ecdf(f(a)))
  function(x) rowMeans(sapply(seq_along(fs), function(i) refs[[i]](fs[[i]](x)))) }
B <- S$bg_regbias
res <- map_dfr(c(names(fitters), "ens_all", "ens_mx_gam_brt", "ens_mx_glm_gam"), function(nm) {
  ff <- switch(nm, ens_all = ens_fitter(names(fitters)), ens_mx_gam_brt = ens_fitter(c("maxent","gam","brt")),
               ens_mx_glm_gam = ens_fitter(c("maxent","glm","gam")), fitters[[nm]])
  t0 <- Sys.time(); cv <- bind_rows(run_cv(ff, B, "loro"), run_cv(ff, B, "block"))
  out <- bind_cols(tibble(model = nm), summarise_cv(cv)); print(out %>% mutate(across(where(is.numeric), ~round(.x, 3))) %>% as.data.frame())
  cat("  time:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  write_csv(cv %>% mutate(model = nm), file.path(outputs_v3, paste0("opt_05_cv_folds_", nm, ".csv"))); out })
write_csv(res, file.path(outputs_v3, "opt_05_algorithm_comparison.csv"))
