# code_7 — Future-climate projections for *H. longicornis* (Australia only)

Australia-only version of `code_4`. Uses the code_6 model and projection grid (continental Australia and
Tasmania); the model itself is the code_3 consensus, calibrated on all records including New Zealand and
the Pacific islands. GCM ensemble, scenarios, predictor derivation and analyses are as in code_4; summaries
are reported for Australia and its states/territories only.

* `Hlongicornis_future_SDM_Australia.qmd` / `.html` — render from this folder with
  `quarto render Hlongicornis_future_SDM_Australia.qmd` (~1 h with 6 cores from a clean cache).
* Future predictor stacks are shared with code_4 (`data/future_climate_cmip6/`, Australia/NZ grid) and
  cropped to the Australian grid on read. If they are missing, run `code_4/prep_future_climate.R` first.

Inputs: `processed_data_6/`, `outputs_6/key_results.rds`, `processed_data_2/`, WorldClim in `data/worldclim/`.
Outputs: `outputs_7/`, `figures_7/`, `processed_data_7/` (ensemble rasters, consensus/bootstrap lambdas);
`processed_data_7/cache/` holds projections, bootstrap projections and attribution (delete to recompute).
