# code_9 — Combined manuscript: *H. longicornis* in Australia (current, future, management)

* `Hlongicornis_Australia_manuscript.qmd` / `.html` — single journal-style paper combining the Australia-only
  analyses `code_6` (current climate), `code_7` (future climate) and `code_8` (management extensions), with
  Supplementary Methods (model development v1→v2→v3, AU+NZ vs Australia-only, future predictors, lambda
  prediction, re-evaluation of code_2 projections, process model, gravity/surveillance, warming levels),
  supplementary tables and figures.
* Fits no models: all numbers are read from `outputs_2–8/`, `processed_data_6–8/` (caches, lambdas) and
  figures from `figures_6–8/`. Renders in ~1 min: `quarto render Hlongicornis_Australia_manuscript.qmd`.
  Re-render after any upstream analysis is re-run.
* New figures (workflow is a Mermaid diagram in the qmd): `figures_9/Fig7_change_synthesis.png`,
  `figures_9/FigS1_model_development.png`.
* `references.bib`, `manuscript.css` — copied from code_8.
* ODMAP protocol (Table S24) and data products (Table S25). Rendering also writes `outputs_9/`:
  `rasters/` (cloud-optimised GeoTIFFs: current, 16 future SSP × period stacks, 0.25° establishment risk,
  plus `raster_band_descriptions.csv`), `Hlongicornis_interactive_map.html` (standalone Leaflet map with all
  layers; Esri satellite basemap tiles need internet; legends bottom left, layer control top right) and `ODMAP_Hlongicornis_Australia.csv`.
* Figures 4, 6, 10 and 11 each have an interactive map directly beneath the static figure.
