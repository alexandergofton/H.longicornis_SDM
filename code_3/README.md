# code_3 — Current-climate SDM for *H. longicornis* (consensus model)

* `Hlongicornis_current_SDM.qmd` / `.html` — manuscript-style report (abstract, introduction, methods,
  results, discussion, supplementary material) with all code embedded. Render from this folder with
  `quarto render Hlongicornis_current_SDM.qmd` (~30–40 min from a clean cache; inputs come from
  `processed_data_2/` and `data/`).
* `references.bib`, `manuscript.css` — bibliography and styling (green boxes = primary figures/tables,
  grey boxes = supplementary).
* `optimisation/` — exploratory scripts used to diagnose the v2 model and test improvements
  (background design, predictor sets, complexity, algorithms, consensus). Not needed to render the report.

Outputs: `outputs_3/` (tables, CV results, key_results.rds), `figures_3/` (300-dpi PNGs of all figures),
`processed_data_3/` (backgrounds, bias surface, suitability / SD / MESS rasters, fitted models).
