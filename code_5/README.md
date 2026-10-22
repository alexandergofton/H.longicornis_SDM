# code_5 — Management-oriented extensions of the H. longicornis SDM

* `Hlongicornis_impact_extensions.qmd` / `.html` — manuscript-style report (same structure as code_3/code_4)
  covering extensions 1, 2, 3, 4, 5, 7 and 8 proposed in code_4 Section 5:
  livestock/people/economic exposure; gravity-model introduction pressure and establishment risk; a
  process-based (degree-day + saturation-deficit) life-cycle model; microclimatic refugia; global warming
  levels; risk-based surveillance design; Papua New Guinea/Pacific projection and One Health exposure.
  Render from this folder with `quarto render Hlongicornis_impact_extensions.qmd` (~10 min with cache).
* `prep_extension_data.R` — standalone data download/preparation (also embedded in the qmd; skipped when
  files exist): WorldClim vapour pressure, GPW v4 population, ESA WorldCover, 30" elevation, Pacific grid,
  and CMIP6 monthly tmin/tmax/prec for 50 GCM x scenario combinations (data/future_climate_cmip6_monthly/).
  GLW3 livestock (data/livestock_glw3/) and GMIA v5 irrigation (data/irrigation_gmia5/) were downloaded
  manually from Harvard Dataverse and FAO; the warming-level table is fetched automatically.

Inputs: code_3 model/rasters (processed_data_3/), code_4 projections and caches (processed_data_4/, outputs_4/).
Outputs: outputs_5/ (CSV summaries, key_results.rds), figures_5/ (300-dpi PNGs), processed_data_5/cache/
(slow steps; delete to recompute, ~30 min).
Only parthenogenetic H. longicornis occur in Australia/NZ/Pacific; no lineage analysis is included.
