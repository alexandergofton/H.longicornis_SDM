# code_8 — Management-oriented extensions of the H. longicornis SDM (Australia only)

Australia-only version of `code_5`, covering extensions 1, 2, 3, 4, 5 and 7 (livestock/people/economic
exposure; gravity-model introduction pressure; process-based life-cycle model; microclimatic refugia;
global warming levels; surveillance design). Extension 8 (Papua New Guinea and the Pacific) is omitted.
The gravity model now has two landmasses (mainland Australia and Tasmania, Bass Strait penalty 0.1).

* `Hlongicornis_impact_extensions_Australia.qmd` / `.html` — render from this folder with
  `quarto render Hlongicornis_impact_extensions_Australia.qmd`.
* Additional data are shared with code_5 (stored on the Australia/NZ window, cropped on read). If any are
  missing the report stops; run `code_5/prep_extension_data.R` (and see code_5/README.md for GLW3/GMIA).

Inputs: code_6 (`processed_data_6/`), code_7 (`processed_data_7/`, its cache, `outputs_7/`).
Outputs: `outputs_8/`, `figures_8/`, `processed_data_8/cache/` (slow steps; delete to recompute).
Only parthenogenetic H. longicornis occur in Australia; no lineage analysis is included.
