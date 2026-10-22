# code_6 — Current-climate SDM for *H. longicornis*, projected to Australia only

Australia-only version of `code_3`. Calibration is unchanged: all 389 thinned records (East Asia, North
America and Oceania, including New Zealand and the Pacific islands), the same region-balanced target-group
background, predictor sets, tuning, consensus members, seeds and threshold. Only the projection domain,
outputs and text differ. Predictions, maps and area summaries cover **continental Australia and Tasmania**:
the Australian states and territories (Natural Earth 1:10 m, including near-shore islands), excluding
Lord Howe Island, Macquarie Island and the external territories (112–154°E, 9–44°S).

* `Hlongicornis_current_SDM_Australia.qmd` / `.html` — manuscript-style report. Render from this folder
  with `quarto render Hlongicornis_current_SDM_Australia.qmd` (~40 min from a clean cache).
* New relative to code_3: transfer evaluation restricted to Australian records (Oceania-withheld
  consensus vs. Australian evaluation background), Australian in-sample metrics, areas by state/territory.
* `references.bib`, `manuscript.css` — as code_3. Exploratory optimisation scripts: see `code_3/optimisation/`.

Outputs: `outputs_6/` (tables, CV results, key_results.rds), `figures_6/` (300-dpi PNGs),
`processed_data_6/` (backgrounds, bias surface, fitted models, `suitability_current_aus.tif`,
`suitability_bootstrap_sd_aus.tif`, `mess_current_aus.tif`, `known_records_aus.csv`).
`processed_data_6/cache/` holds the slow calibration steps; delete it to recompute.
