#!/usr/bin/env python3
"""Browser GUI for selecting upper/lower Figure 6 height bands."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from urllib.parse import urlparse

os.environ.setdefault("MPLBACKEND", "Agg")

import h5py
import numpy as np

DEFAULT_STATISTICS_H5 = Path("figs/paper_refresh_20260721_current/paper_orbit_catalogue_statistics.h5")
DEFAULT_RADIANTS_H5 = Path("figs/paper_refresh_20260721_current/paper_radiant_results.h5")
DEFAULT_SELECTION_H5 = Path("figs/height_band_selection.h5")
DEFAULT_COMPARISON_PNG = Path("figs/height_band_comparison.png")


def polygon_to_path(polygon: np.ndarray) -> str:
    polygon = np.asarray(polygon, dtype=np.float64).reshape(-1, 2)
    if len(polygon) == 0:
        return ""
    pieces = [f"M {polygon[0, 0]:.6g},{polygon[0, 1]:.6g}"]
    pieces.extend(f"L {x:.6g},{y:.6g}" for x, y in polygon[1:])
    pieces.append("Z")
    return " ".join(pieces)


def path_to_polygon(path: str) -> np.ndarray:
    values = [float(x) for x in re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", path or "")]
    if len(values) < 6:
        return np.empty((0, 2), dtype=np.float64)
    if len(values) % 2:
        values = values[:-1]
    return np.asarray(values, dtype=np.float64).reshape(-1, 2)


def read_polygon(h5: h5py.File, name: str) -> np.ndarray:
    if f"polygons/{name}" not in h5:
        return np.empty((0, 2), dtype=np.float64)
    return np.asarray(h5[f"polygons/{name}"], dtype=np.float64).reshape(-1, 2)


def load_polygons(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        return {"upper": np.empty((0, 2)), "lower": np.empty((0, 2))}
    with h5py.File(path, "r") as h5:
        return {"upper": read_polygon(h5, "upper"), "lower": read_polygon(h5, "lower")}


def mask_from_polygon(x: np.ndarray, y: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    from matplotlib.path import Path as MplPath

    if len(polygon) < 3:
        return np.zeros(len(x), dtype=bool)
    points = np.column_stack((np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)))
    return MplPath(np.asarray(polygon, dtype=np.float64)).contains_points(points)


def load_statistics(statistics_h5: Path):
    with h5py.File(statistics_h5, "r") as h5:
        return (
            np.asarray(h5["height_velocity_count"], dtype=np.float64),
            np.asarray(h5["height_edges_km"], dtype=np.float64),
            np.asarray(h5["speed_edges_km_s"], dtype=np.float64),
        )


def load_rows_and_masks(radiants_h5: Path, polygons: dict[str, np.ndarray]):
    with h5py.File(radiants_h5, "r") as h5:
        rows = h5["radiants"][()]
    height = np.asarray(rows["first_alt_km"], dtype=np.float64)
    speed = np.asarray(rows["speed_km_s"], dtype=np.float64)
    good = np.isfinite(height) & np.isfinite(speed) & (height > 0.0) & (speed > 0.0)
    rows = rows[good]
    height = height[good]
    speed = speed[good]
    return (
        rows,
        mask_from_polygon(speed, height, polygons["upper"]),
        mask_from_polygon(speed, height, polygons["lower"]),
    )


def write_band_selection_group(h5: h5py.File, name: str, rows: np.ndarray, mask: np.ndarray) -> np.ndarray:
    selected = rows[mask]
    sample_idx = np.asarray(selected["sample_idx"], dtype=np.int64)
    group = h5.create_group(name)
    group.attrs["count"] = int(len(selected))
    group.attrs["event_id_definition"] = "event_id is the orbit/radiant event sample_idx at 1 MHz sample rate"
    group.create_dataset("event_id", data=sample_idx, compression="gzip", shuffle=True)
    group.create_dataset("sample_idx", data=sample_idx, compression="gzip", shuffle=True)
    group.create_dataset("radiants", data=selected, compression="gzip", shuffle=True)
    scalar_fields = {
        "initial_detection_height_km": "first_alt_km",
        "v_g_km_s": "speed_km_s",
        "radiant_ra_gcrs_deg": "radiant_ra_gcrs_deg",
        "radiant_dec_gcrs_deg": "radiant_dec_gcrs_deg",
        "radiant_lambda_ecliptic_deg": "radiant_lambda_ecliptic_deg",
        "radiant_beta_ecliptic_deg": "radiant_beta_ecliptic_deg",
        "lambda_minus_sun_deg": "lambda_minus_sun_deg",
        "plot_longitude_deg": "plot_longitude_deg",
        "solar_longitude_deg": "sun_lambda_ecliptic_deg",
        "combined_score": "combined_score",
    }
    for out_name, row_name in scalar_fields.items():
        if row_name in selected.dtype.names:
            group.create_dataset(out_name, data=np.asarray(selected[row_name]), compression="gzip", shuffle=True)
    if "hypothesis" in selected.dtype.names:
        group.create_dataset("hypothesis", data=np.asarray(selected["hypothesis"]), compression="gzip", shuffle=True)
    return sample_idx


def save_basic_selection(path: Path, polygons: dict[str, np.ndarray], radiants_h5: Path, statistics_h5: Path):
    rows, mask_upper, mask_lower = load_rows_and_masks(radiants_h5, polygons)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5:
        h5.attrs["script"] = Path(__file__).name
        h5.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        h5.attrs["radiants_h5"] = str(radiants_h5)
        h5.attrs["statistics_h5"] = str(statistics_h5)
        h5.attrs["polygon_columns"] = "speed_km_s, initial_detection_height_km"
        h5.attrs["schema"] = "height_band_selection_v2"
        h5.attrs["event_id_definition"] = "event_id is the orbit/radiant event sample_idx at 1 MHz sample rate"
        h5.attrs["selection_source"] = "closed polygons in geocentric speed versus initial detection height"
        h5.attrs["upper_count"] = int(np.sum(mask_upper))
        h5.attrs["lower_count"] = int(np.sum(mask_lower))
        grp = h5.create_group("polygons")
        for name, polygon in polygons.items():
            grp.create_dataset(name, data=np.asarray(polygon, dtype=np.float32))
        upper_sample_idx = write_band_selection_group(h5, "upper", rows, mask_upper)
        lower_sample_idx = write_band_selection_group(h5, "lower", rows, mask_lower)
        h5.create_dataset("upper_sample_idx", data=upper_sample_idx, compression="gzip", shuffle=True)
        h5.create_dataset("lower_sample_idx", data=lower_sample_idx, compression="gzip", shuffle=True)
        h5.create_dataset("upper_event_id", data=upper_sample_idx, compression="gzip", shuffle=True)
        h5.create_dataset("lower_event_id", data=lower_sample_idx, compression="gzip", shuffle=True)
        selected_event_id = np.concatenate((upper_sample_idx, lower_sample_idx))
        selected_band = np.concatenate(
            (
                np.full(len(upper_sample_idx), b"upper", dtype="S8"),
                np.full(len(lower_sample_idx), b"lower", dtype="S8"),
            )
        )
        order = np.argsort(selected_event_id)
        h5.create_dataset("selected_event_id", data=selected_event_id[order], compression="gzip", shuffle=True)
        h5.create_dataset("selected_sample_idx", data=selected_event_id[order], compression="gzip", shuffle=True)
        h5.create_dataset("selected_band", data=selected_band[order], compression="gzip", shuffle=True)
    return int(np.sum(mask_upper)), int(np.sum(mask_lower))


class HeightBandServer(HTTPServer):
    def __init__(self, address, handler, args):
        super().__init__(address, handler)
        self.args = args
        self.message = ""


class Handler(BaseHTTPRequestHandler):
    server: HeightBandServer

    def log_message(self, fmt, *args):  # noqa: D401
        print(f"{self.address_string()} - {fmt % args}")

    def send_json(self, payload: dict, status: int = 200) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        if path == "/data.json":
            self.send_json(self.data_payload())
            return
        if path != "/":
            self.send_error(404)
            return
        body = self.html().encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        if path != "/save":
            self.send_error(404)
            return
        try:
            n = int(self.headers.get("Content-Length", "0"))
            request = json.loads(self.rfile.read(n).decode("utf-8"))
            if "upperPoints" in request or "lowerPoints" in request:
                polygons = {
                    "upper": np.asarray(request.get("upperPoints", []), dtype=np.float64).reshape(-1, 2),
                    "lower": np.asarray(request.get("lowerPoints", []), dtype=np.float64).reshape(-1, 2),
                }
            else:
                polygons = {
                    "upper": path_to_polygon(request.get("upperPath", "")),
                    "lower": path_to_polygon(request.get("lowerPath", "")),
                }
            upper_count, lower_count = save_basic_selection(
                self.server.args.selection_h5,
                polygons,
                self.server.args.radiants_h5,
                self.server.args.statistics_h5,
            )
            import matplotlib

            matplotlib.use("Agg", force=True)
            from select_height_bands import (
                load_band_data,
                make_comparison_plot,
                selected_event_table,
                write_summary_text,
            )

            data = load_band_data(self.server.args.radiants_h5, polygons)
            upper_events = selected_event_table(self.server.args.orbit_metadata_dir, data.sample_idx[data.mask_upper])
            lower_events = selected_event_table(self.server.args.orbit_metadata_dir, data.sample_idx[data.mask_lower])
            make_comparison_plot(
                self.server.args.comparison_output,
                self.server.args.statistics_h5,
                polygons,
                data,
                upper_events,
                lower_events,
            )
            write_summary_text(self.server.args.comparison_output.with_suffix(".csv"), upper_events, lower_events, data)
            self.server.message = (
                f"Saved upper={upper_count} lower={lower_count}; "
                f"{self.server.args.selection_h5}; {self.server.args.comparison_output}"
            )
            self.send_json(
                {
                    "ok": True,
                    "upperCount": upper_count,
                    "lowerCount": lower_count,
                    "selection": str(self.server.args.selection_h5),
                    "comparison": str(self.server.args.comparison_output),
                    "summary": str(self.server.args.comparison_output.with_suffix(".csv")),
                }
            )
        except Exception as exc:
            self.send_json({"ok": False, "error": repr(exc)}, status=500)

    def data_payload(self) -> dict:
        hv, height_edges, speed_edges = load_statistics(self.server.args.statistics_h5)
        z = np.asarray(hv, dtype=np.float64)
        z[z <= 0.0] = np.nan
        z = np.log10(z)
        z_payload = [[None if not np.isfinite(v) else float(v) for v in row] for row in z.tolist()]
        polygons = load_polygons(self.server.args.selection_h5)
        return {
            "speedEdges": speed_edges.tolist(),
            "heightEdges": height_edges.tolist(),
            "logCount": z_payload,
            "upperPath": polygon_to_path(polygons["upper"]),
            "lowerPath": polygon_to_path(polygons["lower"]),
            "upperPoints": polygons["upper"].tolist(),
            "lowerPoints": polygons["lower"].tolist(),
            "message": self.server.message,
        }

    def html(self) -> str:
        return """<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>PANSY height-band selector</title>
  <style>
    body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; background: #111; color: #eee; }
    #bar { display: flex; gap: 10px; align-items: center; padding: 10px 14px; background: #1e1e1e; border-bottom: 1px solid #333; }
    button { font-size: 14px; padding: 7px 12px; border: 1px solid #666; border-radius: 4px; background: #2b2b2b; color: #eee; cursor: pointer; }
    button.active { border-color: #fff; box-shadow: 0 0 0 2px #fff4 inset; }
    #status { margin-left: 12px; font-size: 13px; color: #ccc; }
    #wrap { position: relative; width: 100vw; height: calc(100vh - 55px); overflow: hidden; }
    #heat, #overlay { position: absolute; inset: 0; width: 100%; height: 100%; }
    text { fill: #eee; font-size: 14px; }
    .axis { stroke: #eee; stroke-width: 1.2; }
    .grid { stroke: #444; stroke-width: 0.7; }
    .upper { stroke: #54a0ff; fill: rgba(84,160,255,0.13); }
    .lower { stroke: #ff9f43; fill: rgba(255,159,67,0.13); }
    .active-line { fill: none; stroke-width: 3; }
    .point { stroke: white; stroke-width: 1; }
  </style>
</head>
<body>
<div id="bar">
  <button id="upper" class="active">Upper band</button>
  <button id="lower">Lower band</button>
  <button id="clear">Clear active</button>
  <button id="save">Save + compare</button>
  <span id="status">Left-click points. Right-click or double-click closes the active polygon.</span>
</div>
<div id="wrap">
  <canvas id="heat"></canvas>
  <svg id="overlay"></svg>
</div>
<script>
let active = "upper";
let upperPoints = [];
let lowerPoints = [];
let current = [];
let payload = null;
const color = {upper: "#54a0ff", lower: "#ff9f43"};
const margin = {l: 72, r: 24, t: 24, b: 58};

function buttonState() {
  document.getElementById("upper").classList.toggle("active", active === "upper");
  document.getElementById("lower").classList.toggle("active", active === "lower");
}

function canvasSize() {
  const wrap = document.getElementById("wrap");
  const w = wrap.clientWidth, h = wrap.clientHeight;
  const dpr = window.devicePixelRatio || 1;
  const canvas = document.getElementById("heat");
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  canvas.style.width = w + "px";
  canvas.style.height = h + "px";
  const svg = document.getElementById("overlay");
  svg.setAttribute("viewBox", `0 0 ${w} ${h}`);
  return {w, h, dpr};
}

function plotBox() {
  const wrap = document.getElementById("wrap");
  const w = wrap.clientWidth, h = wrap.clientHeight;
  return {x0: margin.l, y0: margin.t, x1: w - margin.r, y1: h - margin.b, w, h};
}

function sx(x) {
  const b = plotBox();
  return b.x0 + (x - 0) / 80 * (b.x1 - b.x0);
}

function sy(y) {
  const b = plotBox();
  return b.y1 - (y - 60) / 100 * (b.y1 - b.y0);
}

function invx(px) {
  const b = plotBox();
  return (px - b.x0) / (b.x1 - b.x0) * 80;
}

function invy(py) {
  const b = plotBox();
  return 60 + (b.y1 - py) / (b.y1 - b.y0) * 100;
}

function magma(t) {
  t = Math.max(0, Math.min(1, t));
  const stops = [
    [0.00, [0, 0, 4]], [0.25, [80, 18, 123]], [0.50, [183, 55, 121]],
    [0.75, [251, 136, 97]], [1.00, [252, 253, 191]]
  ];
  for (let i = 1; i < stops.length; i++) {
    if (t <= stops[i][0]) {
      const [t0, c0] = stops[i - 1], [t1, c1] = stops[i];
      const u = (t - t0) / (t1 - t0);
      return c0.map((v, j) => Math.round(v + u * (c1[j] - v)));
    }
  }
  return stops[stops.length - 1][1];
}

function drawCanvas() {
  if (!payload) return;
  const {w, h, dpr} = canvasSize();
  const ctx = document.getElementById("heat").getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  ctx.fillStyle = "#111"; ctx.fillRect(0, 0, w, h);
  const values = payload.logCount.flat().filter(v => v !== null);
  const vmin = Math.min(...values), vmax = Math.max(...values);
  const he = payload.heightEdges, se = payload.speedEdges, z = payload.logCount;
  for (let iy = 0; iy < z.length; iy++) {
    for (let ix = 0; ix < z[iy].length; ix++) {
      const v = z[iy][ix];
      ctx.fillStyle = v === null ? "#000" : `rgb(${magma((v - vmin) / (vmax - vmin)).join(",")})`;
      const x0 = sx(se[ix]), x1 = sx(se[ix + 1]);
      const y0 = sy(he[iy]), y1 = sy(he[iy + 1]);
      ctx.fillRect(x0, y1, x1 - x0 + 1, y0 - y1 + 1);
    }
  }
  drawOverlay();
}

function svg(tag, attrs) {
  const e = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  return e;
}

function polyPoints(points) {
  return points.map(p => `${sx(p[0])},${sy(p[1])}`).join(" ");
}

function drawPoly(root, points, band, closed) {
  if (!points.length) return;
  root.appendChild(svg(closed ? "polygon" : "polyline", {
    points: polyPoints(points), class: closed ? band : "active-line",
    stroke: color[band], fill: closed ? "" : "none"
  }));
  for (const p of points) root.appendChild(svg("circle", {
    cx: sx(p[0]), cy: sy(p[1]), r: 4.5, fill: color[band], class: "point"
  }));
}

function drawOverlay() {
  const root = document.getElementById("overlay");
  const b = plotBox();
  root.replaceChildren();
  for (let x = 0; x <= 80; x += 10) {
    root.appendChild(svg("line", {x1: sx(x), y1: b.y0, x2: sx(x), y2: b.y1, class: "grid"}));
    const t = svg("text", {x: sx(x), y: b.y1 + 24, "text-anchor": "middle"}); t.textContent = x; root.appendChild(t);
  }
  for (let y = 60; y <= 160; y += 10) {
    root.appendChild(svg("line", {x1: b.x0, y1: sy(y), x2: b.x1, y2: sy(y), class: "grid"}));
    const t = svg("text", {x: b.x0 - 10, y: sy(y) + 4, "text-anchor": "end"}); t.textContent = y; root.appendChild(t);
  }
  root.appendChild(svg("rect", {x: b.x0, y: b.y0, width: b.x1 - b.x0, height: b.y1 - b.y0, fill: "none", class: "axis"}));
  const xl = svg("text", {x: (b.x0 + b.x1) / 2, y: b.h - 16, "text-anchor": "middle"}); xl.textContent = "Geocentric velocity, v_g (km/s)"; root.appendChild(xl);
  const yl = svg("text", {x: 18, y: (b.y0 + b.y1) / 2, transform: `rotate(-90 18 ${(b.y0 + b.y1) / 2})`, "text-anchor": "middle"}); yl.textContent = "Initial detection height (km)"; root.appendChild(yl);
  drawPoly(root, upperPoints, "upper", true);
  drawPoly(root, lowerPoints, "lower", true);
  drawPoly(root, current, active, false);
}

function finish() {
  if (current.length >= 3) {
    if (active === "upper") upperPoints = current.slice(); else lowerPoints = current.slice();
  }
  current = [];
  drawOverlay();
}

fetch("/data.json").then(r => r.json()).then(d => {
  payload = d;
  upperPoints = d.upperPoints || [];
  lowerPoints = d.lowerPoints || [];
  document.getElementById("status").textContent = d.message || "Left-click points. Right-click or double-click closes the active polygon.";
  drawCanvas();
});

document.getElementById("upper").onclick = () => {
  active = "upper"; buttonState();
  current = [];
  drawOverlay();
};
document.getElementById("lower").onclick = () => {
  active = "lower"; buttonState();
  current = [];
  drawOverlay();
};
document.getElementById("clear").onclick = () => {
  current = [];
  if (active === "upper") upperPoints = []; else lowerPoints = [];
  drawOverlay();
};
document.getElementById("save").onclick = () => {
  document.getElementById("status").textContent = "Saving and generating comparison plot...";
  fetch("/save", {
    method: "POST",
    headers: {"Content-Type": "application/json"},
    body: JSON.stringify({upperPoints, lowerPoints})
  }).then(r => r.json()).then(d => {
    if (!d.ok) throw new Error(d.error);
    document.getElementById("status").textContent =
      `Saved upper=${d.upperCount} lower=${d.lowerCount}; ${d.comparison}`;
  }).catch(err => {
    document.getElementById("status").textContent = "Save failed: " + err;
  });
};
document.getElementById("overlay").addEventListener("contextmenu", ev => ev.preventDefault());
document.getElementById("overlay").addEventListener("mousedown", ev => {
  const rect = document.getElementById("wrap").getBoundingClientRect();
  const px = ev.clientX - rect.left;
  const py = ev.clientY - rect.top;
  const b = plotBox();
  if (px < b.x0 || px > b.x1 || py < b.y0 || py > b.y1) return;
  if (ev.button === 2) {
    finish();
    return;
  }
  if (ev.button !== 0) return;
  const x = invx(px), y = invy(py);
  if (ev.detail >= 2) {
    finish();
    return;
  }
  current.push([x, y]);
  drawOverlay();
});
window.addEventListener("resize", drawCanvas);
</script>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--statistics-h5", type=Path, default=DEFAULT_STATISTICS_H5)
    parser.add_argument("--radiants-h5", type=Path, default=DEFAULT_RADIANTS_H5)
    parser.add_argument("--selection-h5", type=Path, default=DEFAULT_SELECTION_H5)
    parser.add_argument("--comparison-output", type=Path, default=DEFAULT_COMPARISON_PNG)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=None)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8787)
    parser.add_argument("--no-open", action="store_true")
    args = parser.parse_args()

    server = HeightBandServer((args.host, args.port), Handler, args)
    url = f"http://{args.host}:{server.server_port}/"
    print(url, flush=True)
    if not args.no_open:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
