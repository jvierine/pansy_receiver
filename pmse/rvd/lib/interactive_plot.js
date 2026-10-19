// ../rustplotlib-pansy-rvd/web/src/interactive_plot.ts
import init, { FieldView, PlotViewport } from "../pkg/rustplotlib_web.js";
var initialization;
var clamp = (v) => Math.max(0, Math.min(1, v));
var svgURL = (svg) => URL.createObjectURL(new Blob([svg], { type: "image/svg+xml" }));
var InteractiveFieldPlot = class _InteractiveFieldPlot {
  constructor(root, options) {
    this.root = root;
    this.options = options;
    const s = options.spec;
    if (!s.xlabel?.trim() || !s.ylabel?.trim() || !s.colorbar_label?.trim()) throw Error("Interactive plots require x, y and colorbar labels");
    if (!s.panels.length || s.bounds.length !== 4 || !s.bounds.every(Number.isFinite)) throw Error("Invalid interactive figure specification");
    this.onSelect = options.onSelect;
    this.limits = options.colorLimits || [0, 1];
    const toolbar = document.createElement("div");
    toolbar.style.cssText = "display:flex;gap:12px;align-items:center;font:inherit;margin:6px 0";
    this.resetButton.textContent = "Reset view";
    this.resetButton.onclick = () => this.reset();
    const help = document.createElement("span");
    help.textContent = "Wheel: zoom \xB7 drag: box zoom \xB7 Shift-drag: pan \xB7 double-click: reset";
    help.style.fontSize = "0.875em";
    toolbar.append(this.resetButton, help);
    root.append(toolbar);
    root.style.overflowX = "auto";
    this.stage.style.cssText = "position:relative;background:white;touch-action:none";
    this.axes.style.cssText = "position:absolute;inset:0;width:100%;height:100%;pointer-events:none";
    this.axes.alt = `${s.ylabel} versus ${s.xlabel}; ${s.colorbar_label}`;
    this.stage.append(this.axes);
    root.append(this.stage);
    s.panels.forEach((name, index) => {
      const canvas = document.createElement("canvas");
      canvas.width = 10;
      canvas.height = 10;
      canvas.tabIndex = 0;
      canvas.setAttribute("aria-label", `${name}; ${s.ylabel} versus ${s.xlabel}`);
      canvas.style.cssText = "position:absolute;z-index:1;touch-action:none";
      const label = document.createElement("img");
      label.alt = name;
      label.style.cssText = "position:absolute;width:80px;height:24px;z-index:2;pointer-events:none";
      this.stage.append(canvas, label);
      const panel = new InteractiveFieldPanel(this, canvas, label, index);
      this.panels.push(panel);
      this.bind(panel);
    });
    this.selectionBox.style.cssText = "display:none;position:absolute;z-index:3;border:1px dashed #1d5577;background:#3e97bc22;pointer-events:none";
    this.stage.append(this.selectionBox);
    this.observer = new ResizeObserver(() => this.queueFrame());
  }
  viewport = new PlotViewport();
  panels = [];
  stage = document.createElement("div");
  axes = document.createElement("img");
  resetButton = document.createElement("button");
  observer;
  limits;
  revision = 0;
  timer;
  abort;
  frame;
  frameURL;
  labelURLs = [];
  drag;
  selectionBox = document.createElement("div");
  onSelect;
  static async mount(root, options) {
    initialization ??= init();
    await initialization;
    const plot = new _InteractiveFieldPlot(root, options);
    await plot.requestFrame();
    plot.observer.observe(root);
    return plot;
  }
  setColorLimits(low, high) {
    if (!Number.isFinite(low) || !Number.isFinite(high) || high <= low) throw Error("Invalid color limits");
    if (this.limits[0] === low && this.limits[1] === high) return;
    this.limits = [low, high];
    this.panels.forEach((p) => p.field.set_scale(low, high, 0, 0.1, 2));
    this.queueFrame();
  }
  setSpec(spec) {
    this.options.spec = spec;
    this.reset();
  }
  reset() {
    this.viewport.reset();
    this.changed();
  }
  physicalBounds() {
    const v = this.viewport.bounds(), b = this.options.spec.bounds;
    return [b[0] + v[0] * (b[2] - b[0]), b[1] + v[1] * (b[3] - b[1]), b[0] + v[2] * (b[2] - b[0]), b[1] + v[3] * (b[3] - b[1])];
  }
  applyView() {
    this.panels.forEach((p) => p.applyView());
  }
  changed() {
    this.applyView();
    this.axes.style.visibility = "hidden";
    this.queueFrame();
  }
  queueFrame() {
    clearTimeout(this.timer);
    this.timer = window.setTimeout(() => this.requestFrame().catch((e) => this.reportError(e)), 100);
  }
  reportError(e) {
    this.root.setAttribute("data-plot-error", String(e));
    this.root.dispatchEvent(new CustomEvent("rustplotlib-error", { detail: e }));
  }
  async requestFrame() {
    const request = ++this.revision;
    this.abort?.abort();
    this.abort = new AbortController();
    const width = Math.max(360, Math.min(1200, Math.round(this.root.clientWidth)));
    const endpoint = this.options.frameEndpoint || new URL("./api/frame", location.href).href;
    const response = await fetch(endpoint, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      signal: this.abort.signal,
      body: JSON.stringify({ ...this.options.frameParameters, id: this.options.id, width, viewport: Array.from(this.viewport.bounds()), color_limits: this.limits })
    });
    if (!response.ok) throw Error(`Native plot frame failed (${response.status})`);
    const frame = await response.json();
    if (request !== this.revision) return;
    this.frame = frame;
    this.stage.style.width = frame.width + "px";
    this.stage.style.height = frame.height + "px";
    if (this.frameURL) URL.revokeObjectURL(this.frameURL);
    this.frameURL = svgURL(frame.svg);
    this.axes.src = this.frameURL;
    this.labelURLs.forEach(URL.revokeObjectURL);
    this.labelURLs = [];
    frame.panels.forEach(([x, y, w, h], i) => {
      const p = this.panels[i];
      Object.assign(p.canvas.style, { left: x * 100 + "%", top: y * 100 + "%", width: w * 100 + "%", height: h * 100 + "%" });
      const dpr = Math.min(devicePixelRatio || 1, 2);
      p.canvas.width = Math.max(2, Math.round(w * frame.width * dpr));
      p.canvas.height = Math.max(2, Math.round(h * frame.height * dpr));
      const url = svgURL(frame.panel_labels[i]);
      this.labelURLs.push(url);
      p.label.src = url;
      p.label.style.left = `calc(${x * 100}% + 4px)`;
      p.label.style.top = `calc(${y * 100}% + 4px)`;
      p.field.set_scale(this.limits[0], this.limits[1], 0, 0.1, 2);
      p.applyView();
    });
    this.axes.style.visibility = "visible";
    this.root.removeAttribute("data-plot-error");
  }
  point(panel, event) {
    const r = panel.canvas.getBoundingClientRect();
    return [clamp((event.clientX - r.left) / r.width), clamp(1 - (event.clientY - r.top) / r.height)];
  }
  bind(panel) {
    const c = panel.canvas;
    c.addEventListener("wheel", (e) => {
      e.preventDefault();
      const [x, y] = this.point(panel, e);
      this.viewport.zoom(x, y, Math.exp(Math.max(-1, Math.min(1, e.deltaY * 2e-3))));
      this.changed();
    }, { passive: false });
    c.addEventListener("dblclick", (e) => {
      e.preventDefault();
      this.reset();
    });
    c.addEventListener("keydown", (e) => {
      if (e.key === "Escape" || e.key === "Home") {
        e.preventDefault();
        this.reset();
      }
    });
    c.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      e.preventDefault();
      const [x, y] = this.point(panel, e);
      this.drag = { panel: panel.index, x, y, lastX: x, lastY: y, pan: e.shiftKey };
      c.setPointerCapture(e.pointerId);
    });
    c.addEventListener("pointermove", (e) => {
      const d = this.drag;
      if (!d || d.panel !== panel.index) return;
      const [x, y] = this.point(panel, e);
      if (d.pan) {
        this.viewport.pan(d.lastX - x, d.lastY - y);
        d.lastX = x;
        d.lastY = y;
        this.changed();
      } else {
        const r = c.getBoundingClientRect(), s = this.stage.getBoundingClientRect();
        Object.assign(this.selectionBox.style, { display: "block", left: r.left - s.left + Math.min(x, d.x) * r.width + "px", top: r.top - s.top + (1 - Math.max(y, d.y)) * r.height + "px", width: Math.abs(x - d.x) * r.width + "px", height: Math.abs(y - d.y) * r.height + "px" });
      }
    });
    c.addEventListener("pointerup", (e) => {
      const d = this.drag;
      if (!d) return;
      this.drag = void 0;
      this.selectionBox.style.display = "none";
      const [x, y] = this.point(panel, e);
      if (!d.pan) {
        if (Math.abs(x - d.x) + Math.abs(y - d.y) > 0.015) {
          if (this.viewport.zoom_box(d.x, d.y, x, y)) this.changed();
        } else {
          const b = this.physicalBounds();
          this.onSelect?.({ panel: panel.index, x: b[0] + x * (b[2] - b[0]), y: b[1] + y * (b[3] - b[1]) });
        }
      }
    });
    c.addEventListener("pointercancel", () => {
      this.drag = void 0;
      this.selectionBox.style.display = "none";
    });
  }
  dispose() {
    this.observer.disconnect();
    clearTimeout(this.timer);
    this.abort?.abort();
    this.panels.forEach((p) => p.field.free());
    this.viewport.free();
    if (this.frameURL) URL.revokeObjectURL(this.frameURL);
    this.labelURLs.forEach(URL.revokeObjectURL);
    this.root.replaceChildren();
  }
};
var InteractiveFieldPanel = class {
  constructor(plot, canvas, label, index) {
    this.plot = plot;
    this.canvas = canvas;
    this.label = label;
    this.index = index;
    this.field = new FieldView(canvas);
    this.field.set_missing_color(0.18, 0.23, 0.3, 1);
  }
  field;
  dataBounds;
  set_rectilinear_data(x, y, z) {
    this.dataBounds = [x[0], y[0], x[x.length - 1], y[y.length - 1]];
    this.field.set_rectilinear_data(x, y, z);
    this.applyView();
  }
  set_scale(low, high, _mode = 0, _width = 0.1, _palette = 2) {
    this.plot.setColorLimits(low, high);
  }
  set_moment_data(values, velocity, width) {
    this.field.set_moment_data(values, velocity, width);
    this.applyView();
  }
  applyView() {
    if (!this.dataBounds) return;
    const b = this.plot.physicalBounds(), d = this.dataBounds;
    this.field.set_view((b[0] - d[0]) / (d[2] - d[0]), (b[1] - d[1]) / (d[3] - d[1]), (b[2] - d[0]) / (d[2] - d[0]), (b[3] - d[1]) / (d[3] - d[1]));
  }
};
export {
  InteractiveFieldPanel,
  InteractiveFieldPlot
};
