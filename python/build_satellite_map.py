#!/usr/bin/env python3
"""
build_satellite_map.py
---------------------------------------------------------------------------
Georeference the study site and build a satellite-view map of it.

Every tree in `tree_data` and every soil-collection point in
`bioassay_primary` is recorded as a POLAR offset (compass bearing, distance
in metres) from its stand centre.  The three stand centres are the only
absolute positions we have:

    formerly attenuated strain  37.31776, -76.88626
    virulent strain             37.31766, -76.88659
    no-fungus control           37.31725, -76.88538

This script turns those polar offsets into WGS84 latitude/longitude, then
writes:

    output/spatial/site_points.csv        every point, flat table
    output/spatial/site_points.geojson    same, plus stand hulls + 10 m rings
    output/maps/satellite_map.html        interactive Esri satellite map

The HTML map pulls its imagery tiles from the browser, so open it locally
(or serve it) rather than expecting tiles to be baked into the file.

Conventions
-----------
* Bearing is measured clockwise from north, matching the polar_to_xy() used
  throughout R/ (x = d*sin(theta) east, y = d*cos(theta) north).
* DECLINATION_DEG rotates every bearing before projection.  It defaults to 0,
  i.e. bearings are treated as TRUE north, exactly as the ANI analyses treat
  them.  If the field compass was uncorrected magnetic, set this to the local
  declination (about -11 deg, i.e. 11 deg W, for Williamsburg VA in 2023) to
  rotate the whole point cloud onto true north.
"""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
DATA_FILE = ROOT / "bioassay_data_v_final.xlsx"
SPATIAL_DIR = ROOT / "output" / "spatial"
MAP_DIR = ROOT / "output" / "maps"

DECLINATION_DEG = 0.0          # see module docstring
SAMPLING_RADIUS_M = 10.0       # radius of the soil-sampling design

# stand centres, as supplied from the field GPS
STAND_CENTRES = {
    "attenuated":        (37.31776, -76.88626),
    "virulent":          (37.31766, -76.88659),
    "no_fungus_control": (37.31725, -76.88538),
}

# tree_data$Plot / bioassay_primary$treatment -> canonical stand key
TREATMENT_TO_STAND = {
    "vnaa140_2019":    "attenuated",
    "vnaa140_2023":    "virulent",
    "vnaa140_control": "no_fungus_control",
}

STAND_LABEL = {
    "no_fungus_control": "No-fungus control",
    "attenuated":        "Formerly attenuated strain",
    "virulent":          "Virulent strain",
}
STAND_ORDER = ["no_fungus_control", "attenuated", "virulent"]

# Okabe-Ito, matching R/v2_style.R PLOT_COLOURS
STAND_COLOUR = {
    "no_fungus_control": "#0072B2",
    "attenuated":        "#D55E00",
    "virulent":          "#009E73",
}


# --- geodesy ---------------------------------------------------------------
def metres_per_degree(lat_deg: float) -> tuple[float, float]:
    """Local metres per degree of latitude and longitude (WGS84 series)."""
    phi = math.radians(lat_deg)
    m_lat = (111132.92 - 559.82 * math.cos(2 * phi)
             + 1.175 * math.cos(4 * phi) - 0.0023 * math.cos(6 * phi))
    m_lon = (111412.84 * math.cos(phi) - 93.5 * math.cos(3 * phi)
             + 0.118 * math.cos(5 * phi))
    return m_lat, m_lon


def polar_to_xy(distance_m: float, bearing_deg: float) -> tuple[float, float]:
    """Compass polar offset -> local (east, north) metres. Mirrors R/polar_to_xy."""
    distance_m = 0.0 if distance_m is None else float(distance_m)
    bearing_deg = 0.0 if bearing_deg is None else float(bearing_deg)
    rad = math.radians(bearing_deg + DECLINATION_DEG)
    return distance_m * math.sin(rad), distance_m * math.cos(rad)


def offset_to_lonlat(centre: tuple[float, float], east_m: float,
                     north_m: float) -> tuple[float, float]:
    """(east, north) metres from a centre -> (lon, lat) in WGS84."""
    lat0, lon0 = centre
    m_lat, m_lon = metres_per_degree(lat0)
    return lon0 + east_m / m_lon, lat0 + north_m / m_lat


# --- workbook readers ------------------------------------------------------
def read_sheet(name: str) -> list[dict]:
    wb = load_workbook(DATA_FILE, data_only=True, read_only=True)
    ws = wb[name]
    rows = ws.iter_rows(values_only=True)
    header = [str(h) for h in next(rows)]
    out = [dict(zip(header, r)) for r in rows if any(v is not None for v in r)]
    wb.close()
    return out


def build_trees() -> list[dict]:
    trees = []
    for r in read_sheet("tree_data"):
        stand = r["Plot"]
        if stand not in STAND_CENTRES:
            raise KeyError(f"unmapped Plot value in tree_data: {stand!r}")
        east, north = polar_to_xy(r["Distance_m"], r["Bearing"])
        lon, lat = offset_to_lonlat(STAND_CENTRES[stand], east, north)
        trees.append({
            "kind": "tree",
            "stand": stand,
            "stand_label": STAND_LABEL[stand],
            "id": f"{stand}-T{r['Tree_num']}",
            "tree_num": r["Tree_num"],
            "bearing_deg": r["Bearing"],
            "distance_m": r["Distance_m"],
            "east_m": round(east, 3),
            "north_m": round(north, 3),
            "latitude": round(lat, 8),
            "longitude": round(lon, 8),
            "dbh_cm": r["dbh_cm"],
            "disease_initial": r["disease_initial"],
            "disease_2_mai": r["disease_2_mai"],
            "disease_4_mai": r["disease_4_mai"],
            "audpc_adj_2mai": r["audpc_adj_2mai"],
            "audpc_adj_4mai": r["audpc_adj_4mai"],
        })
    return trees


def build_soil() -> list[dict]:
    """One record per distinct sampling POSITION (the three months repeat it)."""
    seen: dict[tuple, dict] = {}
    for r in read_sheet("bioassay_primary"):
        treatment = r["treatment"]
        if treatment == "neg_control":
            continue                      # off-site reference soils, no coordinates
        stand = TREATMENT_TO_STAND[treatment]
        bearing, distance = r["bearing"], r["distance"]
        key = (stand, bearing, distance)
        if key in seen:
            seen[key]["months"].append(r["month"])
            seen[key]["pooled_samples"].append(r["pooled_sample"])
            continue
        east, north = polar_to_xy(distance, bearing)
        lon, lat = offset_to_lonlat(STAND_CENTRES[stand], east, north)
        position = "Stand centre" if not distance else f"{int(distance)} m @ {int(bearing):03d}°"
        seen[key] = {
            "kind": "soil",
            "stand": stand,
            "stand_label": STAND_LABEL[stand],
            "id": f"{stand}-S-{position}",
            "position": position,
            "bearing_deg": bearing,
            "distance_m": distance,
            "east_m": round(east, 3),
            "north_m": round(north, 3),
            "latitude": round(lat, 8),
            "longitude": round(lon, 8),
            "months": [r["month"]],
            "pooled_samples": [r["pooled_sample"]],
        }
    return list(seen.values())


def build_centres() -> list[dict]:
    return [{
        "kind": "centre",
        "stand": s,
        "stand_label": STAND_LABEL[s],
        "id": f"{s}-centre",
        "east_m": 0.0,
        "north_m": 0.0,
        "latitude": STAND_CENTRES[s][0],
        "longitude": STAND_CENTRES[s][1],
    } for s in STAND_ORDER]


# --- derived geometry ------------------------------------------------------
def convex_hull(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Monotone-chain hull; returns a closed ring."""
    pts = sorted(set(points))
    if len(pts) < 3:
        return pts
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    ring = lower[:-1] + upper[:-1]
    return ring + [ring[0]]


def stand_hulls(trees: list[dict]) -> dict[str, list[tuple[float, float]]]:
    """Mapped-stem convex hull per stand, as a closed (lon, lat) ring."""
    hulls = {}
    for stand in STAND_ORDER:
        local = [(t["east_m"], t["north_m"]) for t in trees if t["stand"] == stand]
        centre = STAND_CENTRES[stand]
        hulls[stand] = [offset_to_lonlat(centre, e, n) for e, n in convex_hull(local)]
    return hulls


def sampling_rings() -> dict[str, list[tuple[float, float]]]:
    """The 10 m soil-sampling ring per stand, as a closed (lon, lat) ring."""
    rings = {}
    for stand in STAND_ORDER:
        centre = STAND_CENTRES[stand]
        ring = []
        for i in range(181):
            theta = 2 * math.pi * i / 180
            ring.append(offset_to_lonlat(centre,
                                         SAMPLING_RADIUS_M * math.sin(theta),
                                         SAMPLING_RADIUS_M * math.cos(theta)))
        rings[stand] = ring
    return rings


# --- writers ---------------------------------------------------------------
CSV_COLUMNS = ["kind", "stand", "stand_label", "id", "tree_num", "position",
               "bearing_deg", "distance_m", "east_m", "north_m",
               "latitude", "longitude", "dbh_cm",
               "disease_initial", "disease_2_mai", "disease_4_mai",
               "audpc_adj_2mai", "audpc_adj_4mai", "months"]


def write_csv(records: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for rec in records:
            row = dict(rec)
            if "months" in row:
                row["months"] = ";".join(str(m) for m in row["months"])
            w.writerow(row)


def write_geojson(records: list[dict], hulls, rings, path: Path) -> None:
    features = []
    for rec in records:
        props = {k: v for k, v in rec.items()
                 if k not in ("latitude", "longitude", "pooled_samples")}
        features.append({
            "type": "Feature",
            "geometry": {"type": "Point",
                         "coordinates": [rec["longitude"], rec["latitude"]]},
            "properties": props,
        })
    for stand in STAND_ORDER:
        features.append({
            "type": "Feature",
            "geometry": {"type": "Polygon",
                         "coordinates": [[list(p) for p in hulls[stand]]]},
            "properties": {"kind": "stand_hull", "stand": stand,
                           "stand_label": STAND_LABEL[stand]},
        })
        features.append({
            "type": "Feature",
            "geometry": {"type": "Polygon",
                         "coordinates": [[list(p) for p in rings[stand]]]},
            "properties": {"kind": "sampling_ring", "stand": stand,
                           "stand_label": STAND_LABEL[stand],
                           "radius_m": SAMPLING_RADIUS_M},
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(
        {"type": "FeatureCollection",
         "crs": {"type": "name",
                 "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},
         "features": features}, indent=1))


def write_html(trees, soil, centres, hulls, rings, path: Path) -> None:
    template = (Path(__file__).parent / "satellite_map_template.html").read_text()
    payload = {
        "trees": trees,
        "soil": [{k: v for k, v in s.items() if k != "pooled_samples"} for s in soil],
        "centres": centres,
        "hulls": {k: [list(p) for p in v] for k, v in hulls.items()},
        "rings": {k: [list(p) for p in v] for k, v in rings.items()},
        "standOrder": STAND_ORDER,
        "standLabel": STAND_LABEL,
        "standColour": STAND_COLOUR,
        "samplingRadiusM": SAMPLING_RADIUS_M,
        "declinationDeg": DECLINATION_DEG,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(template.replace(
        "/*__SITE_DATA__*/null",
        json.dumps(payload, separators=(",", ":"))))


def main() -> None:
    trees, soil, centres = build_trees(), build_soil(), build_centres()
    hulls, rings = stand_hulls(trees), sampling_rings()
    records = centres + trees + soil

    write_csv(records, SPATIAL_DIR / "site_points.csv")
    write_geojson(records, hulls, rings, SPATIAL_DIR / "site_points.geojson")
    write_html(trees, soil, centres, hulls, rings, MAP_DIR / "satellite_map.html")

    lats = [r["latitude"] for r in records]
    lons = [r["longitude"] for r in records]
    print(f"stands {len(centres)} | trees {len(trees)} | soil points {len(soil)}")
    for stand in STAND_ORDER:
        n_t = sum(1 for t in trees if t["stand"] == stand)
        n_s = sum(1 for s in soil if s["stand"] == stand)
        print(f"  {STAND_LABEL[stand]:<28} {n_t:>3} trees  {n_s} soil points")
    print(f"extent  lat {min(lats):.6f}..{max(lats):.6f}  "
          f"lon {min(lons):.6f}..{max(lons):.6f}")
    for p in (SPATIAL_DIR / "site_points.csv",
              SPATIAL_DIR / "site_points.geojson",
              MAP_DIR / "satellite_map.html"):
        print("wrote", os.path.relpath(p, ROOT))


if __name__ == "__main__":
    main()
