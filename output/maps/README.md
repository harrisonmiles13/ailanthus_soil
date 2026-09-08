# Site maps

Satellite-view maps of the three Ailanthus stands, the 187 mapped stems, and the
15 soil-collection points.

Trees and soil samples were recorded in the field as **polar offsets** (compass
bearing + distance in metres) from a stand centre. The stand centres are the only
absolute coordinates in the study:

| Stand | Latitude | Longitude |
|---|---|---|
| Formerly attenuated strain (`vnaa140_2019`) | 37.31776 | -76.88626 |
| Virulent strain (`vnaa140_2023`) | 37.31766 | -76.88659 |
| No-fungus control (`vnaa140_control`) | 37.31725 | -76.88538 |

Both scripts below project the offsets the same way, using the bearing convention
already in `R/06_disease_influence_map.R` (`x = d·sin θ` east, `y = d·cos θ` north,
θ clockwise from north).

## `satellite_map.html` — interactive

Open it in a browser. Imagery tiles are fetched by the browser, so it needs a
network connection; nothing is baked into the file.

- Esri World Imagery basemap (USGS imagery and OpenStreetMap available in the
  layer switcher)
- Stems sized by DBH, colourable by September disease score, AUDPC at 4 MAI, or
  stand; click any stem for its full disease/AUDPC record
- Cyan diamonds = soil-collection points (stand centre + 10 m N/E/S/W, each
  sampled in May, July and September)
- Dashed white ring = the 10 m soil-sampling radius

Regenerate with:

```
python3 python/build_satellite_map.py
```

which also writes `output/spatial/site_points.csv` and `site_points.geojson`
(every point in WGS84, ready for QGIS/ArcGIS or `sf`).

## Static figures

**Automatic:** pushing a change to `R/09_satellite_map.R` or the data file
triggers `.github/workflows/satellite-map.yml`, which renders the figures on
a GitHub-hosted runner (open internet access to the tile servers, unlike the
sandboxes this repo is otherwise developed in) and commits the PNGs straight
to `output/figures/` on `main`. Trigger it by hand from the Actions tab
("Render satellite map figures" → Run workflow) whenever you want a refresh
without a code change.

**Locally**, from the **repository root** (every path in the script is
relative to it):

```
Rscript R/09_satellite_map.R
```

On Windows, if `Rscript` is not on `PATH`, call it by full path — e.g.

```
"C:\Users\harri\R\R-4.4.1\bin\Rscript.exe" R/09_satellite_map.R
```

adjusting for the R version actually installed. First run needs the mapping
packages:

```r
install.packages(c("sf", "maptiles", "tidyterra"))
```

It writes, at 300 dpi, two versions of every panel:

| File | Contents |
|---|---|
| `satellite_map_overview.png` | all three stands, on Esri World Imagery |
| `satellite_map_attenuated.png` | formerly attenuated strain, zoomed, on imagery |
| `satellite_map_virulent.png` | virulent strain, zoomed, on imagery |
| `satellite_map_no_fungus_control.png` | no-fungus control, zoomed, on imagery |
| `satellite_map_overview_blank.png` | all three stands, no basemap, transparent background |
| `satellite_map_<stand>_blank.png` | one zoomed panel per stand, no basemap, transparent background |

Each figure is sized to its own mapped aspect ratio, so the panel fills it
rather than leaving blank bands. Imagery is downloaded on first run and cached
under `output/spatial/tile_cache/` (git-ignored).

A convex hull of the mapped stems per stand is still written to
`site_points.geojson` for GIS use, but neither map draws it.

### The `_blank` figures

A cached basemap can be the wrong vintage for the vegetation being described —
Esri's mosaic for this site, for instance, does not necessarily show canopy
cover from the study period. The `_blank` figures carry only the geometry
(stems, soil points, sampling ring, stand centre, a 5 m scale bar) on a
transparent background, with **no basemap fetch at all** — so they can be
layered in image-editing or GIS software over an independent photo of the
site from the study dates, using the scale bar and the known 10 m sampling
rings to register scale and rotation by eye.

They need only `readxl`, `dplyr`, `ggplot2`, `tibble` and `sf` — not
`maptiles`/`tidyterra` — and no network access, so they render identically
anywhere R runs. Toggle `RENDER_BASEMAP_MAPS` / `RENDER_BLANK_MAPS` at the top
of the script to produce only one set.

## Bearing datum

Both scripts expose `DECLINATION_DEG`, applied to every bearing before
projection, and both default it to **0** — bearings are treated as true north,
exactly as the ANI analyses treat them. If the field compass was uncorrected
magnetic, set it to the local declination (≈ −11° for this site in 2023) to
rotate the whole point cloud onto true north. Stand centres are unaffected.
