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
- Coloured outline = convex hull of the mapped stems; dashed white ring = the
  10 m soil-sampling radius

Regenerate with:

```
python3 python/build_satellite_map.py
```

which also writes `output/spatial/site_points.csv` and `site_points.geojson`
(every point in WGS84, ready for QGIS/ArcGIS or `sf`).

## Static figures

```
Rscript R/09_satellite_map.R
```

writes `output/figures/satellite_map_overview.png` plus one zoomed panel per
stand. It needs `sf`, `maptiles` and `tidyterra`, and downloads Esri imagery on
first run (cached under `output/spatial/tile_cache/`).

## Bearing datum

Both scripts expose `DECLINATION_DEG`, applied to every bearing before
projection, and both default it to **0** — bearings are treated as true north,
exactly as the ANI analyses treat them. If the field compass was uncorrected
magnetic, set it to the local declination (≈ −11° for this site in 2023) to
rotate the whole point cloud onto true north. Stand centres are unaffected.
