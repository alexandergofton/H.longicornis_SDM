# CV framework: common independent evaluation background for all candidate models
suppressPackageStartupMessages({ library(dismo); library(ecospat) })
S <- readRDS(file.path(processed_v3, "modelling_sets.rds"))

fc_args <- function(fc) {
  f <- strsplit(fc, "")[[1]]
  c(paste0("linear=", "L" %in% f), paste0("quadratic=", "Q" %in% f), paste0("hinge=", "H" %in% f),
    paste0("product=", "P" %in% f), "threshold=false", "autofeature=false")
}
fit_mx <- function(p, a, vars, fc, rm) {
  x <- rbind(as.data.frame(p[, vars]), as.data.frame(a[, vars]))
  dismo::maxent(x = x, p = c(rep(1, nrow(p)), rep(0, nrow(a))),
                args = c(fc_args(fc), paste0("betamultiplier=", rm), "doclamp=true", "addsamplestobackground=true"))
}
pred_mx <- function(m, x, vars) as.numeric(dismo::predict(m, as.data.frame(x[, vars]), args = "outputformat=cloglog"))

auc_fun <- function(pp, pa) { r <- rank(c(pp, pa)); (sum(r[seq_along(pp)]) - length(pp)*(length(pp)+1)/2) / (length(pp)*length(pa)) }
cbi_fun <- function(pp, pa) tryCatch(ecospat.boyce(fit = c(pa, pp), obs = pp, nclass = 0, PEplot = FALSE)$cor, error = function(e) NA)

# generic CV: fit_fun(p, a) returns a predict function f(x)
run_cv <- function(fit_fun, bg_train, scheme = c("loro", "block")) {
  scheme <- match.arg(scheme); occ <- S$occ; ev <- S$bg_eval
  folds <- if (scheme == "loro") S$regs else 1:5
  key_o <- if (scheme == "loro") occ$region else occ$block
  key_b <- if (scheme == "loro") bg_train$region else bg_train$block
  key_e <- if (scheme == "loro") ev$region else ev$block
  map_dfr(folds, function(k) {
    tr_p <- occ[key_o != k, ]; te_p <- occ[key_o == k, ]
    if (nrow(te_p) < 5) return(NULL)
    f <- fit_fun(tr_p, bg_train[key_b != k, ])
    pp <- f(te_p); pe <- f(ev[key_e == k, ]); ptr <- f(tr_p)
    thr10 <- quantile(ptr, 0.10)
    tibble(scheme = scheme, fold = as.character(k), n_test = nrow(te_p),
           auc = auc_fun(pp, pe), cbi = cbi_fun(pp, pe), or10 = mean(pp < thr10))
  })
}
mx_fitter <- function(vars, fc, rm) function(p, a) { m <- fit_mx(p, a, vars, fc, rm); function(x) pred_mx(m, x, vars) }

summarise_cv <- function(cv) {
  loro <- cv %>% filter(scheme == "loro"); blk <- cv %>% filter(scheme == "block")
  tibble(loro_auc = mean(loro$auc), loro_cbi = mean(loro$cbi, na.rm = TRUE), loro_or10 = mean(loro$or10),
         oce_auc = loro$auc[loro$fold == "Oceania"], oce_cbi = loro$cbi[loro$fold == "Oceania"],
         na_auc = loro$auc[loro$fold == "North_America"], ea_auc = loro$auc[loro$fold == "East_Asia"],
         blk_auc = mean(blk$auc), blk_cbi = mean(blk$cbi, na.rm = TRUE), blk_or10 = mean(blk$or10))
}
