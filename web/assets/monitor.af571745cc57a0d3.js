// web/lib/interactive_series.js?v=mobile-1
import init, { SeriesView, PlotViewport } from "../pkg/releases/7caedeeb810f57e4/rustplotlib_web.js";
var PlotControlGroup = class {
  element = document.createElement("div");
  members = /* @__PURE__ */ new Set();
  buttons = /* @__PURE__ */ new Map();
  mode = "scroll";
  constructor() {
    this.element.className = "rustplotlib-controls";
    this.element.setAttribute("role", "toolbar");
    this.element.setAttribute("aria-label", "Plot navigation");
    this.element.style.cssText = "display:flex;gap:5px;align-items:center;flex-wrap:wrap;margin:5px 0;font:12px system-ui";
    const button = (text, action) => {
      const b = document.createElement("button");
      b.type = "button";
      b.textContent = text;
      b.onclick = action;
      this.element.append(b);
      return b;
    };
    for (const [mode, label] of [["scroll", "Scroll"], ["zoom", "Zoom"], ["pan", "Pan"]]) this.buttons.set(mode, button(label, () => this.setMode(mode)));
    button("+", () => this.zoom(0.7)).setAttribute("aria-label", "Zoom in");
    button("\u2212", () => this.zoom(1 / 0.7)).setAttribute("aria-label", "Zoom out");
    button("Reset view", () => this.reset());
    const help = document.createElement("span");
    help.textContent = "Scroll page normally \xB7 select Zoom/Pan to drag";
    help.style.cssText = "color:var(--muted,#526b7d);font-size:11px";
    this.element.append(help);
    this.setMode("scroll");
  }
  add(member) {
    this.members.add(member);
    member.modeChanged(this.mode);
    return () => this.members.delete(member);
  }
  setMode(mode) {
    this.mode = mode;
    for (const [key, b] of this.buttons) b.setAttribute("aria-pressed", String(key === mode));
    this.members.forEach((m) => m.modeChanged(mode));
  }
  reset() {
    this.setMode("scroll");
    this.members.forEach((m) => {
      m.viewport.reset();
      m.changed();
    });
  }
  zoom(factor) {
    this.members.forEach((m) => {
      m.viewport.zoom(0.5, 0.5, factor);
      m.changed();
    });
  }
};
function bindPlotNavigation(options) {
  const { canvas, viewport, controls, point, changed, box } = options, abort = new AbortController(), signal = abort.signal;
  let drag;
  const cancel = () => {
    if (drag && canvas.hasPointerCapture(drag.id)) canvas.releasePointerCapture(drag.id);
    drag = void 0;
    box?.();
  };
  const remove = controls.add({ viewport, changed, modeChanged(mode) {
    cancel();
    canvas.style.touchAction = mode === "scroll" ? "pan-y pinch-zoom" : "none";
    canvas.style.cursor = mode === "pan" ? "grab" : mode === "zoom" ? "crosshair" : "auto";
  } });
  canvas.addEventListener("wheel", (e) => {
    if (controls.mode === "scroll") return;
    e.preventDefault();
    const [x, y] = point(e);
    viewport.zoom(x, y, Math.exp(Math.max(-1, Math.min(1, e.deltaY * 2e-3))));
    changed();
  }, { passive: false, signal });
  canvas.addEventListener("dblclick", (e) => {
    e.preventDefault();
    controls.reset();
  }, { signal });
  canvas.addEventListener("keydown", (e) => {
    if (e.key === "Escape" || e.key === "Home") {
      e.preventDefault();
      controls.reset();
    }
  }, { signal });
  canvas.addEventListener("pointerdown", (e) => {
    if (e.button !== 0 || controls.mode === "scroll") return;
    e.preventDefault();
    const [x, y] = point(e);
    drag = { id: e.pointerId, x, y, lastX: x, lastY: y, pan: controls.mode === "pan" || e.shiftKey };
    canvas.setPointerCapture(e.pointerId);
  }, { signal });
  canvas.addEventListener("pointermove", (e) => {
    const d = drag;
    if (!d || d.id !== e.pointerId) return;
    const [x, y] = point(e);
    if (d.pan) {
      viewport.pan(d.lastX - x, d.lastY - y);
      d.lastX = x;
      d.lastY = y;
      changed();
    } else box?.([d.x, d.y, x, y]);
  }, { signal });
  canvas.addEventListener("pointerup", (e) => {
    const d = drag;
    if (!d || d.id !== e.pointerId) return;
    const [x, y] = point(e);
    cancel();
    if (!d.pan) {
      if (Math.abs(x - d.x) + Math.abs(y - d.y) > 0.015) {
        viewport.zoom_box(d.x, d.y, x, y);
        changed();
      } else options.select?.(x, y);
    }
  }, { signal });
  canvas.addEventListener("pointercancel", cancel, { signal });
  canvas.addEventListener("lostpointercapture", cancel, { signal });
  return () => {
    cancel();
    remove();
    abort.abort();
  };
}
var ready;
var fontReady;
function loadPlotFont() {
  return fontReady ??= new FontFace("RustplotlibPagella", `url("${new URL("../fonts/texgyrepagella-regular.otf", import.meta.url)}")`).load().then((font) => {
    document.fonts.add(font);
  }).catch((error) => {
    fontReady = void 0;
    throw error;
  });
}
var InteractiveSeriesPlot = class _InteractiveSeriesPlot {
  constructor(root, spec, height, controls) {
    this.root = root;
    this.spec = spec;
    this.height = height;
    this.specJSON = JSON.stringify(spec);
    const group = controls || new PlotControlGroup();
    const stage = document.createElement("div");
    stage.style.cssText = "position:relative;background:white;";
    this.canvas.style.cssText = "display:block;width:100%";
    this.canvas.tabIndex = 0;
    this.box.style.cssText = "position:absolute;display:none;pointer-events:none;border:1px dashed #12669a;background:#12669a22";
    this.colorButton.type = "button";
    this.colorButton.setAttribute("aria-label", "Adjust color limits");
    this.colorButton.title = "Adjust color limits";
    this.colorButton.style.cssText = "position:absolute;z-index:2;border:0;background:transparent;padding:0;cursor:pointer";
    this.colorButton.onclick = () => this.editColors(stage);
    stage.append(this.canvas, this.box, this.colorButton);
    root.replaceChildren(...controls ? [stage] : [group.element, stage]);
    this.view = new SeriesView(this.canvas);
    this.observer = new ResizeObserver(() => this.schedule());
    this.observer.observe(root);
    const point = (e) => {
      const r = this.canvas.getBoundingClientRect();
      return [Math.max(0, Math.min(1, (e.clientX - r.left - this.rect[0]) / this.rect[2])), Math.max(0, Math.min(1, 1 - (e.clientY - r.top - this.rect[1]) / this.rect[3]))];
    };
    this.unbind = bindPlotNavigation({ canvas: this.canvas, viewport: this.viewport, controls: group, point, changed: () => this.draw(), box: (b) => {
      if (!b) {
        this.box.style.display = "none";
        return;
      }
      const [x0, y0, x, y] = b;
      Object.assign(this.box.style, { display: "block", left: this.rect[0] + Math.min(x, x0) * this.rect[2] + "px", top: this.rect[1] + (1 - Math.max(y, y0)) * this.rect[3] + "px", width: Math.abs(x - x0) * this.rect[2] + "px", height: Math.abs(y - y0) * this.rect[3] + "px" });
    } });
    this.draw();
  }
  viewport = new PlotViewport();
  canvas = document.createElement("canvas");
  view;
  observer;
  rect = [0, 0, 1, 1];
  box = document.createElement("div");
  unbind;
  pending = 0;
  specJSON;
  colorButton = document.createElement("button");
  editor;
  colorOverride;
  static async mount(root, spec, options = {}) {
    ready ??= init().catch((error) => {
      ready = void 0;
      throw error;
    });
    await Promise.all([ready, loadPlotFont()]);
    return new _InteractiveSeriesPlot(root, spec, options.height || 300, options.controls);
  }
  setHeight(height) {
    this.height = height;
    this.schedule();
  }
  setSpec(spec) {
    if (this.colorOverride) {
      spec = { ...spec, ...spec.field ? { field: { ...spec.field, limits: this.colorOverride } } : {}, ...spec.sky ? { sky: { ...spec.sky, limits: this.colorOverride } } : {} };
    }
    const next = JSON.stringify(spec);
    if (next === this.specJSON) return;
    this.spec = spec;
    this.specJSON = next;
    this.draw();
  }
  editColors(stage) {
    this.editor?.remove();
    const editor = document.createElement("div");
    this.editor = editor;
    editor.setAttribute("role", "dialog");
    editor.setAttribute("aria-label", "Color limits");
    editor.style.cssText = "position:absolute;right:10px;top:20px;z-index:5;display:grid;gap:10px;padding:14px;background:white;color:#172134;border:1px solid #9baab5;border-radius:5px;box-shadow:0 3px 15px #0003;font:14px system-ui";
    const limits = Array.from(this.view.color_limits()), inputs = [], outputs = [];
    for (const [i, name] of ["Minimum value", "Maximum value"].entries()) {
      const label = document.createElement("label"), input = document.createElement("input"), output = document.createElement("output");
      label.textContent = name + " ";
      input.type = "range";
      input.setAttribute("aria-label", name);
      const logarithmic = this.spec.sky?.paper_style;
      input.min = String(logarithmic ? limits[0] / 10 : Math.min(0, limits[0] * 2));
      input.max = String(logarithmic ? limits[1] * 2 : Math.max(30, limits[1] * 2));
      input.step = logarithmic ? String(limits[0] / 10) : "0.1";
      input.value = String(limits[i]);
      output.value = String(limits[i]);
      input.style.display = "block";
      label.append(input, output);
      editor.append(label);
      inputs.push(input);
      outputs.push(output);
    }
    const change = () => {
      const low = Number(inputs[0].value), high = Number(inputs[1].value);
      outputs[0].value = String(low);
      outputs[1].value = String(high);
      if (high > low) {
        this.colorOverride = [low, high];
        this.setSpec(this.spec);
      }
    };
    inputs.forEach((input) => input.oninput = change);
    const done = document.createElement("button");
    done.textContent = "Done";
    done.onclick = () => {
      editor.remove();
      this.editor = void 0;
    };
    editor.append(done);
    stage.append(editor);
    inputs[1].focus();
  }
  reset() {
    this.viewport.reset();
    this.draw();
  }
  schedule() {
    cancelAnimationFrame(this.pending);
    this.pending = requestAnimationFrame(() => this.draw());
  }
  draw() {
    const style = getComputedStyle(this.root);
    const width = Math.max(this.spec.compact_axes ? 150 : 260, Math.min(1600, this.root.clientWidth - parseFloat(style.paddingLeft) - parseFloat(style.paddingRight)));
    this.canvas.style.width = width + "px";
    this.canvas.style.height = this.height + "px";
    this.rect = Array.from(this.view.render_plot(this.specJSON, width, this.height, Math.min(devicePixelRatio || 1, 2), this.viewport.bounds()));
    this.colorButton.hidden = !this.spec.field && !this.spec.sky;
    Object.assign(this.colorButton.style, { left: this.rect[0] + this.rect[2] + 8 + "px", top: this.rect[1] + "px", width: "28px", height: this.rect[3] + "px" });
    this.canvas.setAttribute("aria-label", `${this.spec.ylabel} versus ${this.spec.xlabel}`);
    this.root.dataset.renderer = "rustplotlib-wasm-series";
    this.root.dataset.viewBounds = JSON.stringify(Array.from(this.viewport.bounds()));
    this.root.dataset.seriesCount = String(this.spec.series.length);
    this.root.dataset.sampleCount = String(this.spec.series.reduce((n, s) => n + s.points.length, 0));
  }
  dispose() {
    this.unbind();
    this.editor?.remove();
    this.observer.disconnect();
    cancelAnimationFrame(this.pending);
    this.view.free();
    this.viewport.free();
    this.root.replaceChildren();
  }
};

// web/health.mjs
function monitorHealth(ops, packet, now = Date.now() / 1e3) {
  const unknown = (label, detail) => ({ label, state: "unknown", detail });
  const opsTime = Date.parse(ops?.generated_utc) / 1e3, opsFresh = Number.isFinite(opsTime) && now - opsTime >= -60 && now - opsTime < 1800;
  let disk = unknown("Disks", "Awaiting a current disk report"), receiver = unknown("Receiver", "Awaiting a current service report"), phase = unknown("Phases", "Awaiting a current watchdog check");
  if (opsFresh) {
    const disks = ops.disks || [], critical = disks.filter((d) => d.critical || d.used_percent >= 95), warning = disks.filter((d) => d.warning || d.used_percent >= 80);
    const worst = [...disks].filter((d) => Number.isFinite(d.used_percent)).sort((a, b) => b.used_percent - a.used_percent)[0];
    disk = { label: "Disks", state: critical.length ? "bad" : disks.some((d) => d.error) || !disks.length ? "unknown" : warning.length ? "warn" : "ok", detail: worst ? `${worst.path.replace("/media/", "")}: ${worst.used_percent.toFixed(1)}% \xB7 ${(worst.free_bytes / 2 ** 30).toFixed(0)} GiB free` : "Disk measurements unavailable" };
    const services = ops.services || [], failed = services.filter((s) => s.active_state !== "active");
    receiver = { label: "Receiver", state: failed.length ? "bad" : services.length ? "ok" : "unknown", detail: failed.length ? failed.map((s) => s.name.replace("pansy-", "").replace(".service", "") + ": " + s.active_state).join(", ") : "Acquisition and processing services active" };
    const p = ops.phase_health;
    if (p && Number.isFinite(p.checked_unix_s) && now - p.checked_unix_s >= -60 && now - p.checked_unix_s < 900 && typeof p.ok === "boolean") phase = { label: "Phases", state: p.ok ? "ok" : "bad", detail: p.reason };
    else if (p) phase = unknown("Phases", p.ok == null ? p.reason : "Last watchdog phase result is over 15 minutes old");
  }
  const source = packet?.generated_unix_s, raw = packet?.freshness?.Raw;
  if (Number.isFinite(source) && now - source < 1800 && now - source >= -60 && Number.isFinite(raw) && source - raw > 600) receiver = { label: "Receiver", state: "bad", detail: `Raw recording stalled: ${Math.round((source - raw) / 60)} min behind at source` };
  else if (receiver.state === "ok" && (!Number.isFinite(source) || now - source >= 1800 || !Number.isFinite(raw))) receiver = unknown("Receiver", "Services active; current raw-recording health unconfirmed");
  const ages = Object.values(packet?.freshness || {}).filter(Number.isFinite).map((t) => Math.max(0, source - t));
  const processing = Number.isFinite(source) && now - source < 1800 && ages.length ? { label: "Processing", state: Math.max(...ages) >= 43200 ? "bad" : "ok", detail: `Largest source lag: ${(Math.max(...ages) / 3600).toFixed(1)} h \xB7 red at 12 h` } : unknown("Processing", "Awaiting current source-stage ages");
  const items = [disk, receiver, phase, processing];
  return { items, state: items.some((s) => s.state === "bad") ? "bad" : items.some((s) => s.state === "unknown") ? "unknown" : items.some((s) => s.state === "warn") ? "warn" : "ok" };
}

// web/monitor_data.mjs?v=restore-2
var colors = ["#0057d9", "#e52521", "#00872f", "#9b19be", "#ef8a00", "#009eaa", "#a44d00", "#202020"];
function chart(series, xlabel, ylabel, { origin, bounds } = {}) {
  const points = series.flatMap((s) => s.points).filter((p) => p.every((v) => v !== null && Number.isFinite(v)));
  if (!points.length && !bounds) throw Error(`No numeric observations for ${ylabel}`);
  if (!bounds) {
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const [x, y] of points) {
      x0 = Math.min(x0, x);
      x1 = Math.max(x1, x);
      y0 = Math.min(y0, y);
      y1 = Math.max(y1, y);
    }
    const dx = Math.max((x1 - x0) * 0.025, 1e-5), dy = Math.max((y1 - y0) * 0.1, Math.abs(y1) * 0.01, 1e-5);
    bounds = [x0 - dx, y0 - dy, x1 + dx, y1 + dy];
  }
  return { bounds, xlabel, ylabel, series, ...origin !== void 0 ? { utc_origin_seconds: origin } : {} };
}
function zenithTemperature(model) {
  const d = model.noise_model_data, fits = model.noise_zenith_model_status;
  if (!d?.times_h?.length || !fits || fits.length !== 7) throw Error("Waiting for the zenith-only temperature calibration");
  const series = [];
  fits.forEach((fit, ch) => {
    const gain = fit.power_per_kelvin;
    if (!(gain > 0) || !Number.isFinite(gain) || !Number.isFinite(fit.receiver_temp_k)) return;
    series.push({
      name: `RX ${ch}`,
      color: colors[ch],
      line: true,
      radius: 1,
      points: d.times_h.flatMap((t, i) => {
        if (model.noise_model_origin_unix_s + t * 3600 < model.noise_model_observation_stop_unix_s - 86400) return [];
        const p = d.noise_power[i]?.[0]?.[ch];
        return [[t, p > 0 && Number.isFinite(p) ? p / gain : null]];
      })
    });
  });
  return chart(series, "Time (UTC)", "Zenith system noise temperature (K)", { origin: model.noise_model_origin_unix_s });
}
function fitCharts(packet) {
  const fits = packet.fits || [], origin = packet.generated_unix_s - 24 * 3600;
  const bins = new Array(24).fill(0);
  for (const row of fits) {
    const i = Math.floor((row[0] - origin) / 3600);
    if (i >= 0 && i < 24) bins[i]++;
  }
  const counts = chart([{
    color: colors[0],
    name: "Fitted meteors / hour",
    line: true,
    radius: 2,
    points: bins.map((n, i) => [i + 0.5, n])
  }], "Time (UTC)", "Fitted meteors per hour", { origin });
  const speed = chart([{
    color: colors[0],
    line: false,
    radius: 1.4,
    points: fits.map((r) => [(r[0] - origin) / 3600, r[1]])
  }], "Time (UTC)", "Meteor speed (km/s)", { origin, ...!fits.length ? { bounds: [0, 0, 24, 80] } : {} });
  const radiant = {
    ...chart([], "Solar-centered ecliptic longitude (deg)", "Ecliptic latitude (deg)", { bounds: [-180, -90, 180, 90] }),
    sky: { lonlat: fits.filter((r) => Number.isFinite(r[3]) && Number.isFinite(r[4])).map((r) => [r[3], r[4]]), nside: 32, central_longitude: 270, west_positive: true, paper_style: true }
  };
  return { counts, speed, radiant };
}
function phaseChart(packet) {
  const rows = packet.phase || [], origin = packet.generated_unix_s - 86400;
  const spec = chart(
    Array.from({ length: 8 }, (_, ch) => ({
      name: ch === 7 ? "TX ref" : `RX ${ch}`,
      color: colors[ch],
      line: true,
      radius: 2.2,
      line_width: 2.2,
      points: rows.flatMap((row, i) => {
        const p = [(row[0] - origin) / 3600, row[ch + 1] ?? null], previous = rows[i - 1];
        return previous && Math.abs(row[ch + 1] - previous[ch + 1]) > 180 ? [[p[0], null], p] : [p];
      })
    })),
    "Time (UTC)",
    "Relative phase (deg)",
    { origin }
  );
  spec.bounds = [0, spec.bounds[1] - 2, 24, spec.bounds[3] + 2];
  return spec;
}
function rawChart(raw, ch) {
  const channel = raw.channels[ch];
  return { ...chart(
    ["real", "imag"].map((component, i) => ({
      color: colors[i],
      line: true,
      radius: 0,
      points: channel[component].map((v, n) => [n * raw.sample_step_us, v])
    })),
    "Time (\xB5s)",
    ch === 7 ? "TX reference \xB7 ADC" : `RX ${ch} \xB7 ADC`
  ), compact_axes: true, hide_x_labels: ch < 6, panel_label: ch === 7 ? "TX ref \xB7 ADC" : `RX ${ch} \xB7 ADC` };
}
function latencyChart(packet, now = packet.generated_unix_s) {
  const names = ["Raw", "TX", "Matched filter", "Detections", "Cuts", "Modes", "Fits"];
  const labels = ["Raw", "TX", "MF", "Events", "Cuts", "Modes", "Fits"];
  const series = names.flatMap((name, i) => {
    const t = packet.freshness?.[name];
    if (!Number.isFinite(t)) return [];
    const age = Math.max(0, now - t) / 60, healthy = age < 12 * 60;
    return [{ color: healthy ? "#138a38" : "#d32626", line: false, points: [[i, age]] }];
  });
  const max = Math.max(1, ...series.map((s) => s.points[0][1] * 1.15));
  return { ...chart(series, "Processing stage", "Latest sample age (min)", { bounds: [-0.6, 0, 6.6, max] }), bars: true, categories: labels };
}
function meteorCharts(m) {
  const series = m.seconds, mean = series.reduce((a, b) => a + b, 0) / series.length, out = {};
  ["east_km", "north_km", "up_km"].forEach((name, i) => {
    const points = series.map((t, j) => [t, m[name][j]]), err = m.residual_std_km?.[i];
    const traces = [{ color: colors[0], line: false, radius: 2.5, points, errors: Number.isFinite(err) ? points.map(() => 2 * err) : [] }];
    if (m.r0_km && m.v0_km_s) traces.push({ color: "#e32626", line: true, radius: 0, points: series.map((t) => [t, m.r0_km[i] + m.v0_km_s[i] * (t - mean)]) });
    const spec = chart(traces, "Time since first echo (s)", ["East (km)", "North (km)", "Up (km)"][i]);
    if (Number.isFinite(err)) {
      spec.bounds[1] -= 2 * err;
      spec.bounds[3] += 2 * err;
    }
    out[["east", "north", "up"][i]] = spec;
  });
  const height = m.up_km.reduce((a, b) => a + b, 0) / m.up_km.length, extent = height * Math.tan(20 * Math.PI / 180) * 1.04;
  const footprints = [[0, 0], [0, 10], [90, 10], [180, 10], [270, 10]].map(([az, za]) => {
    const a = az * Math.PI / 180, z = za * Math.PI / 180, b = 10 * Math.PI / 180;
    const d = [Math.sin(z) * Math.sin(a), Math.sin(z) * Math.cos(a), Math.cos(z)], e = [Math.cos(a), -Math.sin(a), 0], q = [Math.cos(z) * Math.sin(a), Math.cos(z) * Math.cos(a), -Math.sin(z)];
    return { color: "#8b969f", line: true, radius: 0, line_width: 0.6, fill: "#00000013", points: Array.from({ length: 97 }, (_, i) => {
      const p = i * 2 * Math.PI / 96, u = d.map((v, j) => Math.cos(b) * v + Math.sin(b) * (Math.cos(p) * e[j] + Math.sin(p) * q[j]));
      return [height * u[0] / u[2], height * u[1] / u[2]];
    }) };
  });
  out.ground = { ...chart([...footprints, { color: colors[0], line: false, radius: 2.5, points: series.map((t, i) => [m.east_km[i], m.north_km[i]]) }], "East (km)", "North (km)", { bounds: [-extent, -extent, extent, extent] }), equal_aspect: true };
  out.signal = chart([{ color: colors[0], line: false, radius: 2, points: series.map((t, i) => [t, m.snr_db[i]]) }], "Time since first echo (s)", "SNR (dB)");
  if (m.doppler_km_s) out.doppler = chart([{ color: "#148d28", line: false, radius: 3, points: series.map((t, i) => [t, m.doppler_km_s[i]]) }], "Time since first echo (s)", "Measured Doppler (km/s)");
  else if (m.match_function) out.doppler = chart([{ color: "#148d28", line: false, radius: 3, points: series.map((t, i) => [t, m.match_function[i]]) }], "Time since first echo (s)", "Match function");
  if (m.detection_trace) {
    const d = m.detection_trace;
    out.signal = chart([{ color: colors[0], line: false, radius: 2.5, points: d.seconds.map((t, i) => [t, d.range_km[i]]) }], "Event time (s)", "Detected range (km)");
    if (!m.doppler_km_s) out.doppler = chart([{ color: "#148d28", line: false, radius: 2.5, points: d.seconds.map((t, i) => [t, d.doppler_km_s[i]]) }], "Event time (s)", "Measured Doppler (km/s)");
  }
  if (m.range_time) {
    const r = m.range_time;
    out.signal = { ...chart([], "Time since first echo (s)", "Range (km)", { bounds: [Math.min(...r.seconds), Math.min(...r.range_km), Math.max(...r.seconds), Math.max(...r.range_km)] }), field: { x: axisEdges(r.seconds), y: axisEdges(r.range_km), values: r.snr_db.flat(), limits: [0, 20], label: "SNR (dB)", palette: "turbo" } };
  }
  return out;
}
var axisEdges = (values) => [values[0] - (values[1] - values[0]) / 2, ...values.slice(1).map((v, i) => (v + values[i]) / 2), values.at(-1) + (values.at(-1) - values.at(-2)) / 2];
function rtiChart(manifest2) {
  if (manifest2.schema === "pansy_zenith_display_v1") {
    const nx = manifest2.x.length - 1, values = new Array(nx * (manifest2.y.length - 1)).fill(null);
    for (const [c, r, v] of manifest2.cells) values[r * nx + c] = v;
    return { ...chart([], "Time (UTC)", "Zenith altitude (km)", { origin: manifest2.axis_start_us / 1e6, bounds: [0, 50, 24, 120] }), field: { x: manifest2.x, y: manifest2.y, values, limits: [0, 13], palette: "plasma", label: "SNR (dB)", missing_value: 0 } };
  }
  const data = rtiData(manifest2)[0];
  return { ...chart([], "Time (UTC)", "Zenith altitude (km)", { origin: manifest2.axis_start_us / 1e6, bounds: [0, 50, 24, 120] }), field: { x: Array.from(data.x), y: Array.from(data.y), values: Array.from(data.z, (v) => Number.isFinite(v) ? v : null), limits: [0, 13], palette: "plasma", label: "SNR (dB)", missing_value: 0 } };
}
function rtiData(manifest2) {
  const frames = [...manifest2.frame_packets || []].sort((a, b) => a.start_us - b.start_us);
  const geom = manifest2.geometry, first = manifest2.axis_start_us, nr = geom.slant_range_km.length;
  if (!frames.length) throw Error("No PMSE mode summaries have arrived for this UTC day");
  const columns = [], edges = [];
  for (const f of frames) {
    const left = Math.max(first, f.start_us), right = Math.min(manifest2.axis_stop_us, f.stop_us);
    if (right <= left) continue;
    if (edges.length && left < edges.at(-1)) throw Error("Overlapping PMSE intervals");
    if (edges.length && left > edges.at(-1)) {
      columns.push(null);
      edges.push(left);
    }
    if (!edges.length) edges.push(left);
    columns.push(f);
    edges.push(right);
  }
  const x = Float32Array.from(edges, (t) => (t - first) / 36e8), nx = columns.length;
  const axisEdges2 = (values) => Float32Array.from([
    values[0] - (values[1] - values[0]) / 2,
    ...values.slice(1).map((v, i) => (v + values[i]) / 2),
    values.at(-1) + (values.at(-1) - values.at(-2)) / 2
  ]);
  return Array.from({ length: 5 }, (_, beam) => {
    const z = new Float32Array(nx * nr).fill(NaN);
    columns.forEach((f, c) => {
      if (!f) return;
      f.index.forEach(([b, bin, r], i) => {
        if (b !== beam || r < 0 || r >= nr || !Number.isFinite(f.snr_db[i])) return;
        const cell = r * nx + c;
        z[cell] = Number.isNaN(z[cell]) ? f.snr_db[i] : Math.max(z[cell], f.snr_db[i]);
      });
    });
    return { x, y: axisEdges2(geom.altitude_km[beam]), z };
  });
}

// web/monitor.js
var $ = (id) => document.getElementById(id);
var plots = /* @__PURE__ */ new Map();
var pending = /* @__PURE__ */ new Map();
var telemetry;
var operations;
var rti;
var chosenDay = "";
var manifest;
var followLive = true;
var cardControls = /* @__PURE__ */ new Map();
var plotHeight = (id) => id.startsWith("raw-") ? Number(id.slice(4)) < 6 ? 105 : 132 : id === "meteorGround" ? 280 : id.startsWith("meteor") ? 220 : ["rti", "processing"].includes(id) ? 255 : id === "speed" ? 275 : 300;
var rtiRetry;
var checkingTelemetry = false;
var telemetryError = "";
var utc = (t) => new Date(t * 1e3).toISOString().replace("T", " ").slice(0, 19) + " UTC";
function message(id, text) {
  $(id).textContent = text;
}
function freshness() {
  if (!telemetry || !Number.isFinite(telemetry.generated_unix_s)) return;
  const age = Date.now() / 1e3 - telemetry.generated_unix_s, raw = telemetry.freshness?.Raw, rawAge = Date.now() / 1e3 - raw;
  message("freshness", `Source telemetry ${utc(telemetry.generated_unix_s)} \xB7 raw voltage ${Number.isFinite(rawAge) ? Math.max(0, rawAge).toFixed(0) + " s old" : "unavailable"}${age > 1800 ? " \xB7 telemetry delayed" : ""}${checkingTelemetry ? " \xB7 checking for updates\u2026" : ""}${telemetryError ? " \xB7 " + telemetryError : ""}`);
  $("freshness").classList.toggle("stale", !!telemetryError || age > 1800 || !Number.isFinite(raw) || rawAge > 600);
}
function snapshot(name) {
  try {
    return JSON.parse(sessionStorage.getItem("pansy-monitor-v2:" + name) || "null");
  } catch {
    return null;
  }
}
function saveSnapshot(name, data) {
  try {
    sessionStorage.setItem("pansy-monitor-v2:" + name, JSON.stringify(data));
  } catch {
  }
}
async function json(url) {
  const controller = new AbortController(), timer = setTimeout(() => controller.abort(), 3e4);
  try {
    const response = await fetch(url, { cache: "no-cache", signal: controller.signal });
    if (!response.ok) throw Object.assign(Error(`HTTP ${response.status}`), { status: response.status });
    return await response.json();
  } finally {
    clearTimeout(timer);
  }
}
async function draw(id, spec) {
  const height = plotHeight(id), controls = cardControls.get($(id).closest(".card"));
  if (!plots.has(id)) {
    if (!pending.has(id)) pending.set(id, InteractiveSeriesPlot.mount($(id), spec, { height, controls }).then((p) => {
      plots.set(id, p);
      return p;
    }).finally(() => pending.delete(id)));
    await pending.get(id);
  }
  plots.get(id).setSpec(spec);
}
async function attempt(id, fn) {
  try {
    await fn();
    $(id).classList.remove("error");
  } catch (error) {
    if (!plots.has(id) && !$(id).querySelector("canvas")) message(id, `Awaiting data: ${error.message}`);
    console.warn(id, error);
  }
}
function table(root, headers, rows) {
  const t = document.createElement("table"), head = document.createElement("tr");
  headers.forEach((s) => {
    const th = document.createElement("th");
    th.textContent = s;
    head.append(th);
  });
  t.append(head);
  for (const row of rows) {
    const tr = document.createElement("tr");
    row.forEach((s) => {
      const td = document.createElement("td");
      td.textContent = s;
      tr.append(td);
    });
    t.append(tr);
  }
  root.replaceChildren(t);
}
async function refreshNoise() {
  try {
    message("noiseStatus", "Loading receiver temperatures\u2026");
    let model;
    try {
      model = await json("/rvd/data/monitor-noise.json");
    } catch (error) {
      if (error.status !== 404) throw error;
      model = await json("/rvd/data/sky-noise/status.json");
    }
    await draw("noiseTemperature", zenithTemperature(model));
    const fits = model.noise_zenith_model_status;
    table($("noiseFits"), ["Receiver", "Zenith fitted Trec (K)", "Fit RMS", "Calibration status"], fits.map((f, ch) => [
      `RX ${ch}`,
      Number.isFinite(f.receiver_temp_k) ? f.receiver_temp_k.toFixed(0) : "Unavailable",
      Number.isFinite(f.fractional_rms) ? (100 * f.fractional_rms).toFixed(1) + "%" : "\u2014",
      f.status
    ]));
    const age = Date.now() / 1e3 - model.noise_model_observation_stop_unix_s;
    message("noiseStatus", `Zenith only \xB7 measured system temperature, not receiver offset \xB7 observations through ${utc(model.noise_model_observation_stop_unix_s)}${age > 3600 ? " \xB7 source data delayed" : ""}`);
  } catch (error) {
    message("noiseStatus", `Temperature refresh delayed: ${error.message}. Last successful plot is retained.`);
  }
}
async function refreshRTI() {
  try {
    if (!rti) {
      const cached = snapshot("overview");
      if (cached) {
        try {
          await draw("rti", rtiChart(cached));
          rti = plots.get("rti");
          message("pmseStatus", `${cached.selected_day} \xB7 source spectra through ${utc(cached.latest_mode_us / 1e6)} \xB7 checking for updates\u2026`);
        } catch {
        }
      }
    }
    let latest;
    try {
      latest = await json("/rvd/data/monitor-overview.json");
    } catch {
      latest = await json("/rvd/data/manifest.json");
    }
    const days = latest.available_days || [];
    if (followLive || !chosenDay || !days.includes(chosenDay)) chosenDay = days.at(-1) || latest.selected_day;
    $("day").replaceChildren(...days.map((day) => {
      const option = document.createElement("option");
      option.value = day;
      option.textContent = day;
      return option;
    }));
    $("day").value = chosenDay;
    const data = latest.selected_day === chosenDay ? latest : await json(`/rvd/data/days/${chosenDay}.json`);
    if (!(data.frame_count || data.frame_packets?.length)) throw Error("No completed mesosphere modes available");
    manifest = data;
    if (data.schema === "pansy_zenith_display_v1") saveSnapshot("overview", data);
    await draw("rti", rtiChart(data));
    rti = plots.get("rti");
    clearTimeout(rtiRetry);
    message("pmseStatus", `${chosenDay} \xB7 ${data.frame_count ?? data.frame_packets.length} completed modes \xB7 source spectra through ${utc(data.latest_mode_us / 1e6)}`);
  } catch (error) {
    message("pmseStatus", `PMSE refresh delayed: ${error.message}. Last successful plot is retained.`);
    if (!rti) {
      clearTimeout(rtiRetry);
      rtiRetry = setTimeout(() => void refreshRTI(), 1e4);
    }
  }
}
function renderRaw() {
  if (!telemetry?.raw) return;
  return Promise.all(telemetry.raw.channels.map((_, ch) => attempt(`raw-${ch}`, () => draw(`raw-${ch}`, rawChart(telemetry.raw, ch)))));
}
async function renderTelemetry(data, checking = false) {
  telemetry = data;
  checkingTelemetry = checking;
  freshness();
  const fits = fitCharts(data);
  for (const [id, select] of [["counts", "counts"], ["speed", "speed"], ["radiant", "radiant"]]) await attempt(id, () => draw(id, fits[select]));
  await attempt("phase", () => draw("phase", phaseChart(data)));
  if (data.meteor) {
    const m = data.meteor;
    message("meteorStatus", utc(m.time_unix_s));
    const views = meteorCharts(m), ids = { east: "meteorEast", north: "meteorNorth", up: "meteorUp", ground: "meteorGround", signal: "meteorSignal", doppler: "meteorDoppler" };
    for (const [key, spec] of Object.entries(views)) await attempt(ids[key], () => draw(ids[key], spec));
    const speed = m.v0_km_s ? Math.hypot(...m.v0_km_s).toFixed(1) + " km/s \xB7 " : "";
    message("meteorDetail", speed + (m.residual_std_km ? "Blue: measured points, error bars \xB12 residual RMS \xB7 red: saved straight-line fit. " : "") + (m.range_time ? "" : m.detection_trace ? "Range and Doppler show saved upstream detected peaks; the full range\u2013time field appears after the next fitted event." : "This legacy fit lacks range\u2013time and measured-Doppler arrays; showing the available SNR and match-function diagnostics."));
  }
  if (data.raw) {
    message("rawStatus", `${utc(data.raw.start_unix_s)} \xB7 ${data.raw.units}`);
    await renderRaw();
  } else message("rawStatus", `Raw trace unavailable: ${data.errors?.raw || "no sample yet"}`);
  await attempt("processing", () => draw("processing", latencyChart(data)));
  if (!document.body.dataset.primaryReadyMs) document.body.dataset.primaryReadyMs = String(Math.round(performance.now()));
}
async function refreshTelemetry() {
  checkingTelemetry = true;
  freshness();
  const incoming = json("./monitor-data.json").then((data) => ({ data })).catch((error) => ({ error }));
  try {
    if (!telemetry) {
      const cached = snapshot("telemetry");
      if (cached) {
        try {
          await renderTelemetry(cached, true);
        } catch {
        }
      }
    }
    const { data, error } = await incoming;
    if (error) throw error;
    telemetryError = "";
    await renderTelemetry(data);
    saveSnapshot("telemetry", data);
  } catch (error) {
    telemetryError = `Refresh delayed: ${error.message}`;
    checkingTelemetry = false;
    if (telemetry) freshness();
    else message("freshness", `Telemetry refresh delayed: ${error.message}.`);
    $("freshness").classList.add("stale");
  }
}
async function refreshOperations() {
  try {
    const ops = await json("./ops_status.json");
    operations = ops;
    renderHealth();
    message("operationsStatus", `Source service report ${ops.generated_utc}`);
    const services = ops.services.map((s) => {
      const node = document.createElement("span");
      node.className = "service" + (s.active_state === "active" ? "" : " bad");
      node.textContent = s.name.replace("pansy-", "").replace(".service", "") + ": " + s.active_state;
      return node;
    });
    $("serviceStatus").replaceChildren(...services);
    const inactive = ops.services.filter((s) => s.active_state !== "active").length;
    message("healthSummary", inactive ? `${inactive} services not active` : "all services active");
    $("healthSummary").className = inactive ? "bad" : "ok";
    table($("diskStatus"), ["Disk", "Used", "Free (GiB)", "Health"], ops.disks.map((d) => [d.path, d.used_percent?.toFixed(1) + "%", d.free_bytes == null ? "Unavailable" : (d.free_bytes / 2 ** 30).toFixed(1), d.error ? "Unknown" : d.critical ? "CRITICAL \u226595%" : d.warning ? "Warning \u226580%" : "OK"]));
  } catch (error) {
    message("operationsStatus", `Operations refresh delayed: ${error.message}`);
  }
  try {
    const backup = await json("./status.json");
    const description = document.createElement("p");
    description.textContent = `Backup report ${backup.generated_utc || "undated"}`;
    const remote = document.createElement("div"), mirror = document.createElement("div");
    $("backupStatus").replaceChildren(description, remote, mirror);
    const rows = (channels) => (channels || []).map((c) => [c.channel, c.latest_file_date || "Unavailable", c.rsync_active ? "Transferring" : String(c.last_rsync_exit_code ?? "Unknown")]);
    table(remote, ["Source metadata channel", "Latest UTC day", "Transfer exit/state"], rows(backup.channels));
    table(mirror, ["Mirror metadata channel", "Latest UTC day", "Transfer exit/state"], rows(backup.local_mirror?.channels));
  } catch (error) {
    message("backupStatus", `Backup report unavailable: ${error.message}`);
  }
}
for (let ch = 0; ch < 8; ch++) {
  const root = document.createElement("div");
  root.id = `raw-${ch}`;
  root.className = "plot";
  $("rawVoltage").append(root);
}
var lastHealth = "";
function renderHealth() {
  const health = monitorHealth(operations, telemetry);
  const key = JSON.stringify(health);
  if (key === lastHealth) return;
  lastHealth = key;
  $("systemHealth").dataset.state = health.state;
  message("overallHealth", { ok: "SYSTEM OK", bad: "ACTION REQUIRED", warn: "CHECK DISK HEADROOM", unknown: "HEALTH NOT FULLY CONFIRMED" }[health.state]);
  $("healthItems").replaceChildren(...health.items.map((item) => {
    const node = document.createElement("div");
    node.className = "health-item";
    node.dataset.state = item.state;
    const label = document.createElement("strong"), detail = document.createElement("span");
    label.textContent = `${item.label}: ${item.state === "bad" ? "FAULT" : item.state === "warn" ? "WARNING" : item.state === "unknown" ? "UNKNOWN" : "OK"}`;
    detail.textContent = item.detail;
    node.append(label, detail);
    return node;
  }));
}
var maximized;
function restoreCard() {
  if (!maximized) return;
  const { card, button, scroll, title } = maximized;
  card.classList.remove("maximized");
  button.textContent = "Maximize";
  button.setAttribute("aria-label", `Maximize ${title}`);
  button.setAttribute("aria-expanded", "false");
  document.body.classList.remove("has-maximized");
  for (const [id, p] of plots) if (card.contains($(id))) p.setHeight(plotHeight(id));
  maximized = void 0;
  window.scrollTo(0, scroll);
  button.focus({ preventScroll: true });
}
for (const card of document.querySelectorAll("section.card")) {
  const top = document.createElement("div");
  top.className = "card-top";
  const heading = card.querySelector(":scope > .heading,:scope > h2"), title = card.querySelector("h2")?.textContent || "panel";
  if (heading) top.append(heading);
  const actions = document.createElement("div");
  actions.className = "card-actions";
  if (card.querySelector(".plot,.raw-grid")) {
    const controls = new PlotControlGroup();
    cardControls.set(card, controls);
    actions.append(controls.element);
  }
  const button = document.createElement("button");
  button.type = "button";
  button.className = "maximize";
  button.textContent = "Maximize";
  button.setAttribute("aria-expanded", "false");
  button.setAttribute("aria-label", `Maximize ${title}`);
  button.onclick = () => {
    if (maximized?.card === card) {
      restoreCard();
      return;
    }
    restoreCard();
    maximized = { card, button, scroll: window.scrollY, title };
    card.classList.add("maximized");
    document.body.classList.add("has-maximized");
    button.textContent = "Restore";
    button.setAttribute("aria-label", `Restore ${title}`);
    button.setAttribute("aria-expanded", "true");
    for (const [id, p] of plots) if (card.contains($(id))) p.setHeight(id.startsWith("raw-") ? Math.max(plotHeight(id), Math.floor((innerHeight - 200) / 4)) : Math.max(plotHeight(id), Math.min(650, innerHeight - 220)));
    card.scrollTop = 0;
    button.focus({ preventScroll: true });
  };
  actions.append(button);
  top.append(actions);
  card.prepend(top);
}
document.addEventListener("keydown", (event) => {
  if (event.key === "Escape") restoreCard();
});
$("day").onchange = () => {
  chosenDay = $("day").value;
  followLive = chosenDay === $("day").options[$("day").options.length - 1]?.value;
  void refreshRTI();
};
function clock() {
  message("utc", utc(Date.now() / 1e3));
  freshness();
  renderHealth();
}
clock();
setInterval(clock, 1e3);
var refreshing = false;
async function refresh() {
  if (refreshing) return;
  refreshing = true;
  $("refresh").disabled = true;
  try {
    await Promise.allSettled([refreshRTI(), refreshNoise(), refreshTelemetry(), refreshOperations()]);
  } finally {
    refreshing = false;
    $("refresh").disabled = false;
    if (!document.body.dataset.initialReadyMs) document.body.dataset.initialReadyMs = String(Math.round(performance.now()));
  }
}
$("refresh").onclick = () => void refresh();
void refresh();
setInterval(() => void refresh(), 3e5);
