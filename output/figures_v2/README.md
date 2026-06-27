# Figure outputs — v2

A restyled set of every manuscript/exploratory figure. **The statistics and
fits are byte-for-byte the same as v1** — only the presentation changed. v1 is
conserved in `output/figures/`; nothing there was touched.

## v2 design brief
- **Large text** — bigger base sizes and per-element scaling so figures read
  without zooming.
- **Better use of empty space** — bottom (stacked) legends, trimmed margins,
  larger panels, and short axis titles. The time-matching detail
  (May = 0, July = 2 mai, Sep = 4 mai) now lives in the caption, not on the axis.
- **No main titles** — identification is the caption's job. (The confounding
  caveat on the pooled figures is kept as a small footnote caption, not a title.)
- **Accessibility** —
  - colour-blind-safe Okabe-Ito palette for the three plots
    (blue / vermillion / bluish-green);
  - **redundant encoding** so colour is never the only channel: fitted lines also
    differ by linetype (solid / dashed / dot-dash) and points by shape (sampling
    month: triangle / circle / square);
  - high-contrast annotation labels and a perceptually-uniform `magma` field on
    the map.

## How these are generated
The v2 scripts source the shared style module `R/v2_style.R` and write here:

| script                          | outputs |
|---------------------------------|---------|
| `R/03_response_figures_v2.R`    | 12 single-panel + 4 per-species 3-panel |
| `R/04_pooled_response_figures_v2.R` | 4 pooled 3-panel |
| `R/06_disease_influence_map_v2.R`   | `disease_influence_map.png` |
| `R/07_combined_panel_v2.R`      | `combined_panel.png` (manuscript Figure 5) |

Run from the repo root, e.g. `Rscript R/07_combined_panel_v2.R`. They read the
same inputs as v1 (`output/tables/ani_disease.csv`, `bioassay_data_v_final.xlsx`).

### Note on the map (06)
The v1 map used the `ggnewscale` package to place a second `fill` scale on the
tree markers. v2 removes that dependency: the disease-influence field stays on
`fill` (magma) and each tree's disease score moves to the `colour` aesthetic
(a light→dark blue ring on a white marker), which also reads better for
colour-vision-deficient viewers because the cool ring contrasts with the warm
field.
