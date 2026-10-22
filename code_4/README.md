# code_4 — Future-climate projections for *H. longicornis* (Australia & New Zealand)

* `Hlongicornis_future_SDM.qmd` / `.html` — manuscript-style report (abstract, introduction, methods,
  results, discussion, extensions for practical impact, supplementary material) with all code embedded.
  Includes a critical re-evaluation of the previous projections in `code_2/` (Steps 8–10; Table S1).
  Render from this folder with `quarto render Hlongicornis_future_SDM.qmd`.
* `prep_future_climate.R` — standalone version of the slow data-preparation step (the same code is
  embedded in the qmd and is skipped when the files exist). Reads the Australia/NZ window of the
  WorldClim 2.1 CMIP6 monthly tmin/tmax/prec directly from the server and writes one 7-layer predictor
  stack per scenario to `data/future_climate_cmip6/` (168 files, ~6 MB each).
* `references.bib`, `manuscript.css` — bibliography and styling (green = main text, grey = supplementary).

Inputs: the `code_3` consensus model (`processed_data_3/`, `outputs_3/key_results.rds`), occurrences and
global predictors from `processed_data_2/`, WorldClim current monthly climate in `data/worldclim/`.

Outputs: `outputs_4/` (scenario, ensemble and bootstrap summaries, variance partitioning, key_results.rds),
`figures_4/` (300-dpi PNGs), `processed_data_4/` (ensemble rasters per SSP × period: median suitability,
GCM agreement, GCM SD, fraction novel, change class; consensus and bootstrap lambdas).
`processed_data_4/cache/` holds the results of the slow steps (projections, bootstrap projections,
attribution); delete it to recompute (~1 h with 6 cores). A clean render without cache takes ~1 h.

Notes: run time and memory were tuned for an 8-core, 16 GB Mac. Parallel downloads use socket clusters
(forked processes cannot use libcurl on macOS); projections use forked workers.
