// pmse/rvd/viewer.js
import init3 from "/rvd/pkg/releases/e7774d827ba5b8e0/rustplotlib_web.js";

// pmse/rvd/lib/interactive_plot.js?v=loading-3
import init, { FieldView, PlotViewport, interactive_field_frame } from "/rvd/pkg/releases/e7774d827ba5b8e0/rustplotlib_web.js";

// pmse/rvd/palette.mjs
var palette_default = [[0.03137255, 0.15294118, 0.27843137], [0.12941176, 0.21176471, 0.42352941], [0.12156863, 0.3254902, 0.60392157], [0.07843137, 0.4745098, 0.68627451], [0.10196078, 0.62352941, 0.67058824], [0.41176471, 0.74901961, 0.70588235], [0.79607843, 0.82745098, 0.60784314], [0.95294118, 0.88627451, 0.6627451], [1, 0.97647059, 0.89019608]];

// pmse/rvd/lib/interactive_plot.js?v=loading-3
var initialization;
var clamp = (v) => Math.max(0, Math.min(1, v));
var svgURL = (svg) => URL.createObjectURL(new Blob([svg], { type: "image/svg+xml" }));
var InteractiveFieldPlot = class _InteractiveFieldPlot {
  constructor(root, options2) {
    this.root = root;
    this.options = options2;
    const s = options2.spec;
    if (!s.xlabel?.trim() || !s.ylabel?.trim() || !s.colorbar_label?.trim()) throw Error("Interactive plots require x, y and colorbar labels");
    if (!s.panels.length || s.bounds.length !== 4 || !s.bounds.every(Number.isFinite)) throw Error("Invalid interactive figure specification");
    this.onSelect = options2.onSelect;
    this.limits = [...options2.colorLimits || [0, 1]];
    this.initialLimits = [...this.limits];
    this.momentMode = options2.frameParameters?.hsv === 1;
    this.initialMomentLimits = [options2.frameParameters?.velocity_limit || 30, options2.frameParameters?.width_limit || 10];
    const toolbar = document.createElement("div");
    toolbar.style.cssText = "display:flex;flex-wrap:wrap;gap:10px;align-items:center;font:inherit;margin:0 6px 10px;color:var(--muted,#526b7d)";
    this.resetButton.textContent = "Reset view";
    this.resetButton.onclick = () => this.reset();
    const help = document.createElement("span");
    help.textContent = "Scroll to zoom \xB7 drag to zoom \xB7 Shift-drag to pan \xB7 click colorbar to adjust";
    help.style.fontSize = "0.75rem";
    toolbar.append(this.resetButton, help);
    root.append(toolbar);
    root.style.overflowX = "auto";
    this.stage.style.cssText = "position:relative;background:var(--paper,white);touch-action:none";
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
    this.selectionBox.style.cssText = "display:none;position:absolute;z-index:3;border:1px dashed var(--blue,#1d5577);background:#63a8c922;pointer-events:none";
    this.stage.append(this.selectionBox);
    this.observer = new ResizeObserver(() => {
      const width = Math.max(360, Math.min(1200, Math.round(this.root.clientWidth)));
      if (width !== this.frame?.width) this.queueFrame();
    });
  }
  viewport = new PlotViewport();
  panels = [];
  stage = document.createElement("div");
  axes = document.createElement("img");
  resetButton = document.createElement("button");
  observer;
  limits;
  initialLimits;
  initialMomentLimits;
  revision = 0;
  timer;
  abort;
  frame;
  frameURL;
  labelURLs = [];
  drag;
  selectionBox = document.createElement("div");
  colorButtons = [];
  editor;
  momentMode = false;
  momentLimits;
  onSelect;
  static async mount(root, options2) {
    initialization ??= init().catch((error) => {
      initialization = void 0;
      throw error;
    });
    await initialization;
    const plot2 = new _InteractiveFieldPlot(root, options2);
    try {
      await plot2.requestFrame();
      plot2.observer.observe(root);
      return plot2;
    } catch (error) {
      plot2.dispose();
      throw error;
    }
  }
  setColorLimits(low, high) {
    if (!Number.isFinite(low) || !Number.isFinite(high) || high <= low) throw Error("Invalid color limits");
    if (this.limits[0] === low && this.limits[1] === high) return;
    this.limits = [low, high];
    this.panels.forEach((p) => p.applyView());
    this.queueFrame();
  }
  getColorLimits() {
    return [...this.limits];
  }
  getColorStops() {
    return this.frame?.color_stops;
  }
  getDefaultHsvLimits() {
    return [...this.initialMomentLimits];
  }
  getHsvLimits() {
    return this.momentLimits ? [...this.momentLimits] : void 0;
  }
  setMomentMode(enabled) {
    if (this.momentMode === enabled) return;
    this.momentMode = enabled;
    this.queueFrame();
  }
  setHsvLimits(velocity, width) {
    if (!Number.isFinite(velocity) || velocity <= 0 || !Number.isFinite(width) || width <= 0) throw Error("Invalid HSV limits");
    if (this.momentLimits?.[0] === velocity && this.momentLimits?.[1] === width) return;
    this.momentLimits = [velocity, width];
    this.panels.forEach((p) => p.applyHsvLimits(velocity, width));
    this.queueFrame();
  }
  openColorEditor() {
    this.editor?.remove();
    const dialog = document.createElement("dialog");
    this.editor = dialog;
    dialog.setAttribute("aria-label", "Adjust color limits");
    dialog.style.cssText = "position:fixed;box-sizing:border-box;width:min(440px,calc(100vw - 32px));max-width:none;max-height:calc(100dvh - 32px);overflow:auto;padding:24px;border:1px solid var(--line,#8b99aa);border-radius:10px;background:var(--paper,white);color:var(--ink,#172134);font:15px system-ui;box-shadow:0 12px 40px #0003";
    const title = document.createElement("h2");
    title.textContent = "Color scale";
    title.style.cssText = "font:inherit;font-size:1.15em;font-weight:600;margin:0 0 6px";
    dialog.append(title);
    const description = document.createElement("p");
    description.textContent = this.options.spec.colorbar_label;
    description.style.cssText = "margin:0 0 18px;color:var(--muted,#526b7d)";
    dialog.append(description);
    const span = this.limits[1] - this.limits[0], base = [this.limits[0] - 4 * span, this.limits[1] + 4 * span];
    const row = (label, value, min, max, step, apply) => {
      const line = document.createElement("label");
      line.style.cssText = "display:grid;grid-template-columns:minmax(0,1fr) 90px;gap:8px;align-items:center;margin:16px 0";
      const text = document.createElement("span");
      text.textContent = label;
      text.style.gridColumn = "1 / -1";
      const slider = document.createElement("input");
      slider.type = "range";
      slider.min = String(min);
      slider.max = String(max);
      slider.step = String(step);
      slider.value = String(value);
      slider.setAttribute("aria-label", label);
      slider.style.cssText = "min-width:0;width:100%;margin:0";
      const number = document.createElement("input");
      number.type = "number";
      number.step = "any";
      number.value = String(value);
      number.setAttribute("aria-label", label + " value");
      number.style.cssText = "box-sizing:border-box;min-width:0;width:100%";
      const update = () => {
        const v = number.valueAsNumber, valid = Number.isFinite(v) && apply(v);
        number.setCustomValidity(valid ? "" : "Enter a finite value; maximum must exceed minimum and HSV limits must be positive.");
        number.setAttribute("aria-invalid", String(!valid));
        if (valid) {
          slider.min = String(Math.min(+slider.min, v));
          slider.max = String(Math.max(+slider.max, v));
          slider.value = String(v);
        }
      };
      slider.oninput = () => {
        number.value = slider.value;
        update();
      };
      number.oninput = update;
      line.append(text, slider, number);
      dialog.append(line);
    };
    row("Minimum", this.limits[0], base[0], base[1], span / 500, (v) => {
      if (v >= this.limits[1]) return false;
      this.setColorLimits(v, this.limits[1]);
      return true;
    });
    row("Maximum", this.limits[1], base[0], base[1], span / 500, (v) => {
      if (v <= this.limits[0]) return false;
      this.setColorLimits(this.limits[0], v);
      return true;
    });
    if (this.momentMode) {
      const limits = this.momentLimits || [this.options.frameParameters?.velocity_limit || 30, this.options.frameParameters?.width_limit || 10];
      row("Velocity +/- m/s", limits[0], 0.1, Math.max(limits[0] * 4, 200), 0.1, (v) => {
        if (v <= 0) return false;
        this.setHsvLimits(v, (this.momentLimits || limits)[1]);
        return true;
      });
      row("Width m/s", limits[1], 0.1, Math.max(limits[1] * 4, 100), 0.1, (v) => {
        if (v <= 0) return false;
        this.setHsvLimits((this.momentLimits || limits)[0], v);
        return true;
      });
    }
    const actions = document.createElement("div");
    actions.style.cssText = "display:flex;justify-content:space-between;gap:12px;margin-top:24px";
    const reset = document.createElement("button");
    reset.textContent = "Reset scale";
    reset.onclick = () => {
      this.setColorLimits(...this.initialLimits);
      if (this.momentMode) this.setHsvLimits(...this.initialMomentLimits);
      this.openColorEditor();
    };
    const close = document.createElement("button");
    close.textContent = "Done";
    close.onclick = () => dialog.close();
    actions.append(reset, close);
    dialog.append(actions);
    document.body.append(dialog);
    dialog.addEventListener("click", (event) => {
      if (event.target === dialog) {
        const r = dialog.getBoundingClientRect();
        const e = event;
        if (e.clientX < r.left || e.clientX > r.right || e.clientY < r.top || e.clientY > r.bottom) dialog.close();
      }
    });
    dialog.addEventListener("close", () => dialog.remove(), { once: true });
    dialog.showModal();
  }
  setSpec(spec) {
    if (JSON.stringify(spec) === JSON.stringify(this.options.spec)) return;
    this.options.spec = spec;
    this.reset();
  }
  setFrameParameters(parameters) {
    if (JSON.stringify(parameters) === JSON.stringify(this.options.frameParameters)) return;
    this.options.frameParameters = parameters;
    this.queueFrame();
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
    if (e instanceof DOMException && e.name === "AbortError") return;
    this.root.setAttribute("data-plot-error", String(e));
    this.root.dispatchEvent(new CustomEvent("rustplotlib-error", { detail: e }));
  }
  async requestFrame() {
    const request = ++this.revision;
    this.abort?.abort();
    this.abort = new AbortController();
    const width = Math.max(360, Math.min(1200, Math.round(this.root.clientWidth)));
    const parameters = { ...this.options.frameParameters, ...this.momentLimits ? { velocity_limit: this.momentLimits[0], width_limit: this.momentLimits[1] } : {} };
    let frame;
    if (this.options.frameEndpoint) {
      const response = await fetch(this.options.frameEndpoint, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        signal: this.abort.signal,
        body: JSON.stringify({ ...this.options.frameParameters, ...this.momentLimits ? { velocity_limit: this.momentLimits[0], width_limit: this.momentLimits[1] } : {}, hsv: this.momentMode ? 1 : 0, id: this.options.id, width, viewport: Array.from(this.viewport.bounds()), color_limits: this.limits })
      });
      if (!response.ok) throw Error(`Native plot frame failed (${response.status})`);
      frame = await response.json();
    } else {
      const spec = { ...this.options.spec, color_stops: this.options.spec.color_stops || palette_default, ...this.momentMode ? { hsv_key: { velocity_limit: parameters.velocity_limit || 30, width_limit: parameters.width_limit || 10 } } : {} };
      frame = JSON.parse(interactive_field_frame(JSON.stringify(spec), JSON.stringify({ width, viewport: Array.from(this.viewport.bounds()), color_limits: this.limits })));
    }
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
      p.label.style.display = this.options.spec.panels[i] ? "block" : "none";
      p.label.style.left = `calc(${x * 100}% + 4px)`;
      p.label.style.top = `calc(${y * 100}% + 4px)`;
      p.applyView();
    });
    this.axes.style.visibility = "visible";
    this.root.removeAttribute("data-plot-error");
    this.colorButtons.forEach((b) => b.remove());
    this.colorButtons = [];
    (frame.color_regions || []).forEach(([x, y, w, h]) => {
      const button = document.createElement("button");
      button.type = "button";
      button.setAttribute("aria-label", "Adjust color limits");
      button.setAttribute("aria-haspopup", "dialog");
      button.title = `Adjust ${this.options.spec.colorbar_label} color scale`;
      button.style.cssText = `position:absolute;left:${x * 100}%;top:${y * 100}%;width:${Math.max(24, w * frame.width)}px;height:${h * 100}%;z-index:4;background:transparent;border:0;padding:0;cursor:pointer`;
      button.onclick = () => this.openColorEditor();
      this.stage.append(button);
      this.colorButtons.push(button);
    });
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
    this.editor?.remove();
    this.observer.disconnect();
    clearTimeout(this.timer);
    this.abort?.abort();
    this.panels.forEach((p) => p.dispose());
    this.viewport.free();
    if (this.frameURL) URL.revokeObjectURL(this.frameURL);
    this.labelURLs.forEach(URL.revokeObjectURL);
    this.root.replaceChildren();
  }
};
var InteractiveFieldPanel = class {
  constructor(plot2, canvas, label, index) {
    this.plot = plot2;
    this.canvas = canvas;
    this.label = label;
    this.index = index;
    this.visibility = new IntersectionObserver(() => this.applyView(), { rootMargin: "40px" });
    this.visibility.observe(canvas);
    canvas.addEventListener("webglcontextlost", this.onLost);
    canvas.addEventListener("webglcontextrestored", this.onRestored);
  }
  gpu;
  dataBounds;
  data;
  dirty = true;
  lost = false;
  canRestore = false;
  restoring = false;
  disposed = false;
  extension = null;
  visibility;
  restoreTimer;
  moments;
  onLost = (event) => {
    event.preventDefault();
    this.lost = true;
    this.canRestore = true;
    this.gpu?.free();
    this.gpu = void 0;
    delete this.canvas.dataset.gpuReady;
    if (this.isVisible()) this.restoreTimer = setTimeout(() => this.applyView(), 50);
  };
  onRestored = () => {
    this.lost = false;
    this.canRestore = this.restoring = false;
    if (!this.disposed) this.applyView();
  };
  /** Low-level debug/pixel access may be temporarily unavailable while the GPU restores. */
  get field() {
    const field = this.ensureField();
    if (field && this.dataBounds) this.renderField(field);
    return field;
  }
  isVisible() {
    const r = this.canvas.getBoundingClientRect();
    return this.canvas.isConnected && r.width > 0 && r.height > 0 && r.bottom > -40 && r.top < innerHeight + 40 && r.right > 0 && r.left < innerWidth;
  }
  ensureField() {
    if (this.disposed) return;
    if (this.lost) {
      if (this.canRestore && !this.restoring) {
        this.restoring = true;
        this.extension?.restoreContext();
      }
      return;
    }
    if (!this.gpu) {
      this.gpu = new FieldView(this.canvas);
      this.extension = this.canvas.getContext("webgl2").getExtension("WEBGL_lose_context");
      this.gpu.set_missing_color(...this.plot.options.missingColor || [0.6196, 0.6745, 0.7098, 1]);
      this.dirty = true;
    }
    return this.gpu;
  }
  suspend() {
    if (!this.gpu || !this.extension) return;
    this.gpu.free();
    this.gpu = void 0;
    this.lost = true;
    this.canRestore = this.restoring = false;
    delete this.canvas.dataset.gpuReady;
    this.extension.loseContext();
  }
  set_rectilinear_data(x, y, z) {
    this.moments = void 0;
    this.plot.setMomentMode(false);
    this.dataBounds = [x[0], y[0], x[x.length - 1], y[y.length - 1]];
    this.data = { x, y, z };
    this.dirty = true;
    this.applyView();
  }
  set_scale(low, high, _mode = 0, _width = 0.1, _palette = 2) {
    this.plot.setColorLimits(low, high);
  }
  set_moment_data(values, velocity, width) {
    this.moments = values;
    this.plot.setMomentMode(true);
    this.dirty = true;
    this.plot.setHsvLimits(velocity, width);
    this.applyView();
  }
  applyHsvLimits(_velocity, _width) {
    if (this.moments) {
      this.dirty = true;
      this.applyView();
    }
  }
  applyView() {
    if (this.disposed) return;
    if (!this.isVisible()) {
      this.suspend();
      return;
    }
    if (!this.dataBounds) return;
    const field = this.ensureField();
    if (field) this.renderField(field);
  }
  renderField(field) {
    if (this.dirty && this.data) {
      field.set_rectilinear_data(this.data.x, this.data.y, this.data.z);
      if (this.moments) field.set_moment_data(this.moments, ...this.plot.getHsvLimits() || this.plot.getDefaultHsvLimits());
      this.dirty = false;
    }
    field.set_scale(...this.plot.getColorLimits(), 0, 0.1, 2);
    const stops = this.plot.getColorStops();
    if (stops) field.set_colormap(Float32Array.from(stops));
    const b = this.plot.physicalBounds(), d = this.dataBounds;
    field.set_view((b[0] - d[0]) / (d[2] - d[0]), (b[1] - d[1]) / (d[3] - d[1]), (b[2] - d[0]) / (d[2] - d[0]), (b[3] - d[1]) / (d[3] - d[1]));
    this.canvas.dataset.gpuReady = "true";
  }
  dispose() {
    this.disposed = true;
    clearTimeout(this.restoreTimer);
    this.visibility.disconnect();
    this.canvas.removeEventListener("webglcontextlost", this.onLost);
    this.canvas.removeEventListener("webglcontextrestored", this.onRestored);
    this.suspend();
    this.gpu?.free();
    this.gpu = void 0;
    this.data = this.moments = void 0;
  }
};

// pmse/rvd/imaging.js?v=loading-3
var $ = (id) => document.getElementById(id);
var worker;
var geometry;
var plot;
var plotMount;
var sequence = 0;
var previousRange;
var lastMoments;
var lastPacket;
var lastKey;
var echoFilter;
var busy = false;
var fetchAbort;
var timeout;
function setEchoFilter(resolver) {
  echoFilter = resolver;
  if (lastPacket) reconstruct();
}
var grid = 61;
function edges(extent) {
  const d = 2 * extent / (grid - 1);
  return Float32Array.from({ length: grid + 1 }, (_, i) => -extent - d / 2 + i * d);
}
function styleParameters(range) {
  return { range_km: range, hsv: $("display").value === "hsv" ? 1 : 0, velocity_limit: +$("imageVelocity").value, width_limit: +$("imageWidth").value };
}
function repaint() {
  if (plot && lastMoments) {
    const extent = previousRange * Math.sin(20 * Math.PI / 180);
    plot.panels[0].set_rectilinear_data(edges(extent), edges(extent), Float32Array.from({ length: grid * grid }, (_, i) => lastMoments[4 * i]));
    if ($("display").value === "hsv") plot.panels[0].set_moment_data(lastMoments, ...plot.getHsvLimits() || [+$("imageVelocity").value, +$("imageWidth").value]);
    plot.setFrameParameters(styleParameters(previousRange));
  }
}
function options() {
  return { method: $("sparseImaging").checked ? $("imageMethod").value : "tsvd", strength: $("sparseImaging").checked ? +$("imageStrength").value : 1, reweight: $("imageReweight").checked };
}
function methodNote() {
  const o = options();
  $("sparseOptions").hidden = o.method === "tsvd";
  $("reweightControl").hidden = o.method !== "sparse";
  $("imageMethodNote").textContent = o.method === "wavelet" ? "Positive source power with 2-D Haar-wavelet sparsity. Uncheck for the legacy TSVD comparison." : o.method === "sparse" ? "Positive pixel \u2113\u2081 for compact echoes; optional reweighting and selected-support power refit." : "Legacy 18-mode TSVD comparison: squared inverse magnitude, not sparse source-power reconstruction.";
}
function makeWorker() {
  const active = worker = new Worker(new URL("/rvd/imaging_worker.js?v=persistence-2", location.href), { type: "module" });
  worker.onmessage = async ({ data }) => {
    if (data.id !== sequence) return;
    busy = false;
    clearTimeout(timeout);
    $("imagePlot").setAttribute("aria-busy", "false");
    if (data.error) {
      $("imageStatus").textContent = data.error;
      return;
    }
    const { header, moments } = data, range = header.selected_range_km, extent = range * Math.sin(20 * Math.PI / 180);
    try {
      if (!plot) {
        if (!plotMount) {
          previousRange = range;
          plotMount = InteractiveFieldPlot.mount($("imagePlot"), {
            id: "image",
            frameParameters: styleParameters(range),
            spec: { bounds: [-extent, -extent, extent, extent], xlabel: "East (km)", ylabel: "North (km)", colorbar_label: "Relative power (dB)", panels: [""] },
            colorLimits: [-30, 0]
          });
        }
        try {
          plot = await plotMount;
        } finally {
          plotMount = void 0;
        }
      }
      if (data.id !== sequence) return;
      if (previousRange !== range) {
        previousRange = range;
        plot.setFrameParameters(styleParameters(range));
        plot.setSpec({ bounds: [-extent, -extent, extent, extent], xlabel: "East (km)", ylabel: "North (km)", colorbar_label: "Relative power (dB)", panels: [""] });
      }
      lastMoments = moments;
      plot.panels[0].set_rectilinear_data(edges(extent), edges(extent), Float32Array.from({ length: grid * grid }, (_, i) => moments[4 * i]));
      repaint();
      const signal = Array.from({ length: grid * grid }, (_, i) => moments[4 * i + 3]).some((x) => x > 0);
      const d = data.diagnostics, method = d.method === "wavelet" ? "PMSE wavelet \u2113\u2081" : d.method === "sparse" ? "Compact pixel \u2113\u2081" : "Legacy TSVD";
      $("imageStatus").textContent = `${new Date(header.start_us / 1e3).toISOString()} \u2014 ${new Date(header.stop_us / 1e3).toISOString().slice(11, 19)} UTC \xB7 slant range ${range.toFixed(2)} km \xB7 all five beams \xB7 ${header.looks} looks \xB7 ${header.index.length} retained spectral cells \xB7 Rust/WASM ${data.compute_ms.toFixed(0)} ms. ${method}. ` + (header.hidden_single_mode_cells != null ? `${header.hidden_single_mode_cells} single-mode spectral cells hidden. ` : "") + (signal ? "" : "No reconstructed power at this selection. ") + `Historical calibration: not independently validated for this observation.`;
      $("imageDiagnostics").textContent = d.cells ? `Covariance fit: mean relative residual ${(100 * d.relative_residual).toFixed(1)}% \xB7 weighted RMS ${d.residual_rms.toFixed(2)} \xB7 ${d.cells} cells \xB7 ${(d.iterations / d.cells).toFixed(0)} mean iterations. ` + (d.unconverged ? `${d.unconverged} cells reached the iteration limit; reconstruction is approximate. ` : "Iterate tolerance reached. ") + "Noise weights assume independent Gaussian looks; correlated noise and detection censoring are not modeled." : "Legacy TSVD: no sparse-fit residual available.";
      if (d.cells && d.relative_residual > 0.5) $("imageDiagnostics").textContent += " Large residual: check calibration, angular support and source model.";
      $("imageDiagnostics").classList.toggle("warn", d.unconverged > 0 || d.relative_residual > 0.5);
      $("imagePlot").dataset.method = d.method;
    } catch (e) {
      if (data.id === sequence) $("imageStatus").textContent = String(e);
    }
  };
  worker.onerror = (event) => {
    if (worker !== active) return;
    busy = false;
    clearTimeout(timeout);
    $("imagePlot").setAttribute("aria-busy", "false");
    $("imageStatus").textContent = `Imaging worker failed: ${event.message}`;
    active.terminate();
    worker = void 0;
  };
}
function reconstruct() {
  methodNote();
  if (!lastPacket) return;
  if ($("sparseImaging").checked && !$("imageStrength").checkValidity()) {
    $("imageStrength").reportValidity();
    return;
  }
  const id = ++sequence;
  if (busy) {
    worker.terminate();
    worker = void 0;
    busy = false;
  }
  clearTimeout(timeout);
  if (!worker) makeWorker();
  const packet = lastPacket.slice(0);
  busy = true;
  $("imagePlot").setAttribute("aria-busy", "true");
  $("imageDiagnostics").textContent = "";
  $("imageStatus").textContent = "Reconstructing all five beams in Rust/WebAssembly\u2026";
  worker.postMessage({ id, packet, geometry, grid, options: options(), allowedGates: echoFilter?.(lastKey) }, [packet]);
  timeout = setTimeout(() => {
    if (id === sequence && busy) {
      worker.terminate();
      worker = void 0;
      busy = false;
      $("imagePlot").setAttribute("aria-busy", "false");
      $("imageStatus").textContent = "Reconstruction exceeded 60 seconds. Select a different method or echo to retry.";
    }
  }, 6e4);
}
async function imageSelection(key, slantRange, offsetSeconds = 0, spectraReady = true) {
  lastKey = key;
  const id = ++sequence;
  fetchAbort?.abort();
  fetchAbort = new AbortController();
  const signal = fetchAbort.signal;
  lastPacket = void 0;
  if (busy) {
    worker.terminate();
    worker = void 0;
    busy = false;
  }
  clearTimeout(timeout);
  if (!spectraReady) {
    $("imagePlot").setAttribute("aria-busy", "false");
    $("imageStatus").textContent = "These observations are available in the RTI; their complex spectra are still transferring from Syowa. Retry after the next update.";
    $("imageDiagnostics").textContent = "";
    return;
  }
  $("imageStatus").textContent = "Fetching existing complex spectra\u2026";
  $("imageDiagnostics").textContent = "";
  try {
    geometry ??= await (await fetch("./api/geometry")).json();
    if (id !== sequence) return;
    const response = await fetch(`./api/spectra?key=${encodeURIComponent(key)}&range_km=${slantRange}&offset_s=${offsetSeconds}`, { signal });
    if (!response.ok) throw Error((await response.json()).error || "Spectra unavailable");
    const packet = await response.arrayBuffer();
    if (id !== sequence) return;
    if (packet.byteLength > 64 * 2 ** 20) throw Error("Selected spectra exceed the 64 MiB browser buffer limit");
    lastPacket = packet;
    reconstruct();
  } catch (e) {
    if (id === sequence && e.name !== "AbortError") $("imageStatus").textContent = String(e);
  }
}
for (const id of ["imageVelocity", "imageWidth", "display"]) $(id).addEventListener("change", () => {
  if (id !== "display") plot?.setHsvLimits(+$("imageVelocity").value, +$("imageWidth").value);
  repaint();
});
for (const id of ["sparseImaging", "imageMethod", "imageStrength", "imageReweight"]) $(id).addEventListener("change", reconstruct);
methodNote();

// pmse/rvd/echo_filter.mjs?v=1
var gateKey = ([beam, , range]) => `${beam}:${range}`;
function persistentEchoes(frames2, maxGapUs = 6e8) {
  const present = frames2.map((frame) => new Set(frame.index.filter((_, i) => Number.isFinite(frame.snr_db[i]) && frame.snr_db[i] > 0).map(gateKey)));
  let hiddenCells = 0, hiddenGates = 0;
  const allowed = /* @__PURE__ */ new Map();
  const filtered = frames2.map((frame, i) => {
    const previous = i > 0 && frame.start_us - frames2[i - 1].stop_us <= maxGapUs ? present[i - 1] : /* @__PURE__ */ new Set();
    const next = i + 1 < frames2.length && frames2[i + 1].start_us - frame.stop_us <= maxGapUs ? present[i + 1] : /* @__PURE__ */ new Set();
    const gates = new Set([...present[i]].filter((key) => previous.has(key) || next.has(key)));
    allowed.set(String(frame.key), gates);
    hiddenGates += present[i].size - gates.size;
    const keep = frame.index.map((cell, j) => gates.has(gateKey(cell)) ? j : -1).filter((j) => j >= 0);
    hiddenCells += frame.index.length - keep.length;
    const result = { ...frame, index: keep.map((j) => frame.index[j]), snr_db: keep.map((j) => frame.snr_db[j]) };
    for (const name of ["look_dominance", "look_persistence"]) if (frame[name]) result[name] = keep.map((j) => frame[name][j]);
    return result;
  });
  return { frames: filtered, allowed, hiddenCells, hiddenGates };
}

// pmse/rvd/lib/interactive_series.js?v=runtime-c70022f0
import init2, { SeriesView, PlotViewport as PlotViewport2 } from "/rvd/pkg/releases/e7774d827ba5b8e0/rustplotlib_web.js";
var ready;
var InteractiveSeriesPlot = class _InteractiveSeriesPlot {
  constructor(root, spec, height) {
    this.root = root;
    this.spec = spec;
    this.height = height;
    this.specJSON = JSON.stringify(spec);
    const toolbar = document.createElement("div"), reset = document.createElement("button"), help = document.createElement("span");
    toolbar.style.cssText = "display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin:0 6px 10px";
    reset.textContent = "Reset view";
    reset.onclick = () => this.reset();
    help.textContent = "Scroll to zoom \xB7 drag to zoom \xB7 Shift-drag to pan \xB7 double-click to reset";
    help.style.cssText = "font-size:.75rem;color:var(--muted,#526b7d)";
    toolbar.append(reset, help);
    const stage = document.createElement("div");
    stage.style.cssText = "position:relative;background:white;";
    this.canvas.style.cssText = "display:block;width:100%;touch-action:none";
    this.canvas.tabIndex = 0;
    this.box.style.cssText = "position:absolute;display:none;pointer-events:none;border:1px dashed #12669a;background:#12669a22";
    stage.append(this.canvas, this.box);
    root.replaceChildren(toolbar, stage);
    this.view = new SeriesView(this.canvas);
    this.observer = new ResizeObserver(() => this.schedule());
    this.observer.observe(root);
    const point = (e) => {
      const r = this.canvas.getBoundingClientRect();
      return [Math.max(0, Math.min(1, (e.clientX - r.left - this.rect[0]) / this.rect[2])), Math.max(0, Math.min(1, 1 - (e.clientY - r.top - this.rect[1]) / this.rect[3]))];
    };
    this.canvas.addEventListener("wheel", (e) => {
      e.preventDefault();
      const [x, y] = point(e);
      this.viewport.zoom(x, y, Math.exp(Math.max(-1, Math.min(1, e.deltaY * 2e-3))));
      this.draw();
    }, { passive: false });
    this.canvas.ondblclick = () => this.reset();
    this.canvas.onkeydown = (e) => {
      if (e.key === "Escape" || e.key === "Home") {
        e.preventDefault();
        this.reset();
      }
    };
    this.canvas.onpointerdown = (e) => {
      if (e.button !== 0) return;
      e.preventDefault();
      const [x, y] = point(e);
      this.drag = { x, y, lastX: x, lastY: y, pan: e.shiftKey };
      this.canvas.setPointerCapture(e.pointerId);
    };
    this.canvas.onpointermove = (e) => {
      const d = this.drag;
      if (!d) return;
      const [x, y] = point(e);
      if (d.pan) {
        this.viewport.pan(d.lastX - x, d.lastY - y);
        d.lastX = x;
        d.lastY = y;
        this.draw();
      } else Object.assign(this.box.style, { display: "block", left: this.rect[0] + Math.min(x, d.x) * this.rect[2] + "px", top: this.rect[1] + (1 - Math.max(y, d.y)) * this.rect[3] + "px", width: Math.abs(x - d.x) * this.rect[2] + "px", height: Math.abs(y - d.y) * this.rect[3] + "px" });
    };
    this.canvas.onpointerup = (e) => {
      const d = this.drag;
      if (!d) return;
      this.drag = void 0;
      this.box.style.display = "none";
      const [x, y] = point(e);
      if (!d.pan && Math.abs(x - d.x) + Math.abs(y - d.y) > 0.015) this.viewport.zoom_box(d.x, d.y, x, y);
      this.draw();
    };
    this.canvas.onpointercancel = () => {
      this.drag = void 0;
      this.box.style.display = "none";
    };
    this.draw();
  }
  viewport = new PlotViewport2();
  canvas = document.createElement("canvas");
  view;
  observer;
  rect = [0, 0, 1, 1];
  box = document.createElement("div");
  drag;
  pending = 0;
  specJSON;
  static async mount(root, spec, options2 = {}) {
    ready ??= init2();
    await ready;
    return new _InteractiveSeriesPlot(root, spec, options2.height || 300);
  }
  setSpec(spec) {
    this.spec = spec;
    this.specJSON = JSON.stringify(spec);
    this.draw();
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
    const width = Math.max(280, Math.min(1200, this.root.clientWidth - parseFloat(style.paddingLeft) - parseFloat(style.paddingRight)));
    this.canvas.style.width = width + "px";
    this.canvas.style.height = this.height + "px";
    this.rect = Array.from(this.view.render_plot(this.specJSON, width, this.height, Math.min(devicePixelRatio || 1, 2), this.viewport.bounds()));
    this.canvas.setAttribute("aria-label", `${this.spec.ylabel} versus ${this.spec.xlabel}`);
    this.root.dataset.renderer = "rustplotlib-wasm-series";
    this.root.dataset.viewBounds = JSON.stringify(Array.from(this.viewport.bounds()));
    this.root.dataset.seriesCount = String(this.spec.series.length);
    this.root.dataset.sampleCount = String(this.spec.series.reduce((n, s) => n + s.points.length, 0));
  }
  dispose() {
    this.observer.disconnect();
    cancelAnimationFrame(this.pending);
    this.view.free();
    this.viewport.free();
    this.root.replaceChildren();
  }
};

// pmse/rvd/diagnostic_data.mjs?v=1
var colors = ["#202020", "#2474a8", "#b65c18", "#22854a", "#b52e48"];
var names = ["Zenith", "North", "East", "South", "West"];
function chart(points, xlabel, ylabel, origin, { zero = false, series } = {}) {
  const finite = points.filter((p) => Number.isFinite(p[0]) && Number.isFinite(p[1]));
  if (!finite.length) throw Error(`No numeric data for ${ylabel}`);
  const xs = finite.map((p) => p[0]), ys = finite.map((p) => p[1]);
  let x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = zero ? 0 : Math.min(...ys), y1 = Math.max(...ys);
  if (x1 <= x0) x1 = x0 + 1;
  const pad = Math.max((y1 - y0) * 0.12, Math.abs(y1) * 1e-3, 1e-5);
  y0 = zero ? -pad * 0.15 : y0 - pad;
  y1 += pad;
  return { bounds: [x0, y0, x1, y1], xlabel, ylabel, utc_origin_seconds: origin, series: series || [{ color: "#215d86", points, line: true, radius: 2 }] };
}
function noiseSeries(model, module) {
  const data = model.noise_model_data;
  if (!data?.times_h?.length || !Number.isFinite(model.noise_model_origin_unix_s)) throw Error("Numeric sky-noise data have not arrived yet");
  return ["Receiver noise power (dB)", "Modeled sky temperature (K)"].map((ylabel, panel) => {
    const series = names.map((name, beam) => ({ name, color: colors[beam], line: true, radius: 1.7, points: data.times_h.map((t, i) => [t, panel === 0 ? 10 * Math.log10(data.noise_power[i][beam][module]) : data.sky_temperature_k[i][beam]]) }));
    return chart(series.flatMap((s) => s.points), "Time (UTC)", ylabel, model.noise_model_origin_unix_s, { series });
  });
}
function archiveSeries(archive) {
  if (!archive?.hourly_unix_s?.length) throw Error("Archive inventory data have not arrived yet");
  const origin = archive.hourly_unix_s[0], x = archive.hourly_unix_s.map((t) => (t - origin) / 3600);
  const hourly = archive.hourly_bytes.map((b) => b / 2 ** 20);
  let total = 0;
  const cumulative = [[x[0], 0], ...x.map((t, i) => [t + 1, total += hourly[i]])];
  const retention = x.map((t, i) => [t + 0.5, archive.hourly_dense_retained_fraction?.[i] == null ? null : 100 * archive.hourly_dense_retained_fraction[i]]);
  const specs = [
    chart(x.map((t, i) => [t + 0.5, hourly[i]]), "Time (UTC)", "Data per UTC hour (MiB)", origin, { zero: true }),
    chart(cumulative, "Time (UTC)", "Accumulated data (MiB)", origin, { zero: true })
  ];
  if (retention.some((p) => p[1] != null)) specs.push(chart(retention, "Time (UTC)", "Full spectra retained (%)", origin, { zero: true }));
  for (const spec of specs) {
    spec.bounds[0] = 0;
    spec.bounds[2] = Math.max(1, x.at(-1) + 1);
  }
  return specs;
}

// pmse/rvd/diagnostics.js?v=runtime-c70022f0
var plots = /* @__PURE__ */ new Map();
var mounting = /* @__PURE__ */ new Map();
async function draw(id, spec) {
  const root = document.getElementById(id);
  if (!root) return;
  if (!plots.has(id)) {
    if (!mounting.has(id)) mounting.set(id, InteractiveSeriesPlot.mount(root, spec, { height: 310 }).then((p) => {
      plots.set(id, p);
      return p;
    }).finally(() => mounting.delete(id)));
    await mounting.get(id);
  }
  plots.get(id).setSpec(spec);
}
async function renderNoiseDiagnostics(model, module) {
  const specs = noiseSeries(model, module);
  await Promise.all(specs.map((s, i) => draw(["noisePowerPlot", "skyTemperaturePlot"][i], s)));
}
async function renderArchiveDiagnostics(archive) {
  const specs = archiveSeries(archive);
  await Promise.all(specs.map((s, i) => draw(["hourlyDataPlot", "cumulativeDataPlot", "retentionPlot"][i], s)));
}

// pmse/rvd/loading.mjs?v=1
async function readJSON(response, onProgress = () => {
}) {
  if (!response.ok) throw Error(`Data unavailable (${response.status})`);
  if (!response.body?.getReader) return response.json();
  const size = Number(response.headers.get("Content-Length"));
  const total = !response.headers.get("Content-Encoding") && size > 0 ? size : null;
  const reader = response.body.getReader(), decoder = new TextDecoder(), parts = [];
  let received = 0;
  onProgress({ received, total });
  try {
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      received += value.byteLength;
      parts.push(decoder.decode(value, { stream: true }));
      onProgress({ received, total });
    }
    parts.push(decoder.decode());
    return JSON.parse(parts.join(""));
  } finally {
    reader.releaseLock();
  }
}

// pmse/rvd/viewer.js
var $2 = (id) => document.getElementById(id);
var names2 = ["Zenith", "North", "East", "South", "West"];
var fetchData = (url) => fetch(url, { cache: "no-cache", signal: AbortSignal.timeout(3e4) });
var initialResponse = fetchData("./data/manifest.json").catch((error) => ({ error }));
var manifest;
var frames = [];
var plots2 = [];
var rd;
var playing = null;
var selected = 0;
var revision = "";
var rtiPlot;
var rdPlot;
var tdPlot;
var chosenDay = "";
var refreshSequence = 0;
var originalFrames = [];
var summaryCache = /* @__PURE__ */ new WeakMap();
var catalogue;
var showLoading = true;
var snapshotKey = "pansy-rvd-day-view-v1";
function saveView(packet) {
  try {
    const saved = JSON.stringify(packet);
    if (saved.length < 25e5) sessionStorage.setItem(snapshotKey, saved);
  } catch {
  }
}
async function restoreView() {
  try {
    const packet = JSON.parse(sessionStorage.getItem(snapshotKey) || "null");
    if (!packet?.frame_packets?.length || !packet.geometry || !packet.selected_day) return;
    chosenDay = packet.selected_day;
    manifest = packet;
    originalFrames = packet.frame_packets;
    applyEchoFilter();
    const option = document.createElement("option");
    option.value = chosenDay;
    option.textContent = chosenDay;
    $2("day").replaceChildren(option);
    selected = frames.length - 1;
    $2("frame").max = String(selected);
    $2("frame").value = String(selected);
    await mountPlots();
    renderRTI();
    renderRD();
    revision = packet.revision;
    $2("status").textContent = `Saved view \xB7 ${chosenDay} UTC \xB7 observations through ${utc(packet.latest_mode_us)} \xB7 checking for updates\u2026`;
    document.body.dataset.plotsReadyMs = String(Math.round(performance.now()));
    loaded();
  } catch {
    revision = "";
    try {
      sessionStorage.removeItem(snapshotKey);
    } catch {
    }
  }
}
function loading(stage, progress) {
  if (!showLoading) return;
  $2("loadProgress").hidden = false;
  $2("loadStage").textContent = stage;
  const bar = $2("loadBar");
  if (progress?.total) {
    bar.max = progress.total;
    bar.value = progress.received;
  } else bar.removeAttribute("value");
  $2("loadDetail").textContent = progress?.received ? `${(progress.received / 2 ** 20).toFixed(2)} MiB received (decoded)` : "Please wait\u2026";
}
function loaded() {
  showLoading = false;
  $2("loadProgress").hidden = true;
}
function applyEchoFilter() {
  summaryCache = /* @__PURE__ */ new WeakMap();
  const enabled = $2("hideSingleMode").checked;
  const filtered = enabled ? persistentEchoes(originalFrames) : null;
  frames = filtered ? filtered.frames : originalFrames;
  setEchoFilter(filtered ? (key) => [...filtered.allowed.get(String(key)) || []] : void 0);
  $2("filterStatus").textContent = filtered ? `Single-mode filter: ${filtered.hiddenGates.toLocaleString()} beam/range detections hidden. Requires the same range gate in this beam in a neighboring mode (gap \u226410 min). New echoes wait for confirmation; short real echoes can also be hidden.` : "";
}
var skyNoise;
var noiseRequest;
var archiveInventory;
var secondaryMounts = /* @__PURE__ */ new Map();
var visiblePlots = new IntersectionObserver((entries) => {
  for (const entry of entries) if (entry.isIntersecting) {
    visiblePlots.unobserve(entry.target);
    secondaryMounts.get(entry.target)?.().then(() => {
      renderRTI();
      renderRD();
    }).catch(failure);
  }
}, { rootMargin: "200px" });
var axisStart = () => manifest.axis_start_us ?? manifest.first_mode_us;
var axisStop = () => manifest.axis_stop_us ?? manifest.latest_mode_us;
var dayNumber = () => Number(manifest.selected_day.replaceAll("-", ""));
var utc = (us) => new Date(us / 1e3).toISOString().replace("T", " ").slice(0, 19) + " UTC";
var edges2 = (centers) => Float32Array.from([...centers.map((x, i) => i ? (centers[i - 1] + x) / 2 : x - (centers[1] - x) / 2), centers.at(-1) + (centers.at(-1) - centers.at(-2)) / 2]);
async function mountPlots() {
  const hours = (axisStop() - axisStart()) / 36e8, origin = axisStart() / 1e6;
  const spec = { bounds: [0, 50, Math.max(hours, 1e-3), 120], xlabel: "Time (UTC)", ylabel: "Nominal altitude (km)", colorbar_label: "SNR (dB)", panels: names2, utc_origin_seconds: origin };
  if (!rtiPlot) {
    rtiPlot = await InteractiveFieldPlot.mount($2("rti"), { id: "rti", spec, frameParameters: { day: dayNumber() }, colorLimits: [+$2("low").value, +$2("high").value], onSelect: ({ panel, x, y }) => {
      const time = axisStart() + x * 36e8;
      selected = frames.reduce((best, f, i) => Math.abs((f.start_us + f.stop_us) / 2 - time) < Math.abs((frames[best].start_us + frames[best].stop_us) / 2 - time) ? i : best, 0);
      $2("rdbeam").value = String(panel);
      $2("frame").value = String(selected);
      renderRD();
      const range = y / Math.sin(manifest.geometry.beam_elevation_deg[panel] * Math.PI / 180);
      const offset = Math.max(0, Math.floor((time - frames[selected].start_us) / 12e7) * 120);
      imageSelection(frames[selected].key, range, Math.min(offset, Math.max(0, Math.floor((frames[selected].duration_s - 1) / 120) * 120)), frames[selected].spectra_ready !== false);
    } });
    plots2 = rtiPlot.panels.map((p) => ({ canvas: p.canvas, view: p, owner: rtiPlot }));
  } else rtiPlot.setSpec(spec);
  const mountRD = async () => {
    const [v0, v1] = manifest.velocity_limits_mps;
    if (!rdPlot) {
      rdPlot = await InteractiveFieldPlot.mount($2("rd"), { id: "rd", frameParameters: { day: dayNumber() }, spec: { bounds: [v0, 50, v1, 120], xlabel: "Radial velocity (m/s)", ylabel: "Nominal altitude (km)", colorbar_label: "SNR (dB)", panels: ["Spectrum"] }, colorLimits: [+$2("low").value, +$2("high").value], onSelect: ({ y }) => imageSelection(frames[selected].key, y / Math.sin(manifest.geometry.beam_elevation_deg[+$2("rdbeam").value] * Math.PI / 180), 0, frames[selected].spectra_ready !== false) });
      rd = { view: rdPlot.panels[0], owner: rdPlot };
    }
    rdPlot.setFrameParameters({ day: dayNumber() });
  };
  const mountTD = async () => {
    const [v0, v1] = manifest.velocity_limits_mps;
    const tdSpec = { bounds: [0, v0, Math.max(hours, 1e-3), v1], xlabel: "Time (UTC)", ylabel: "Radial velocity (m/s)", colorbar_label: "SNR (dB)", panels: names2, utc_origin_seconds: origin };
    if (!tdPlot) {
      tdPlot = await InteractiveFieldPlot.mount($2("timeDoppler"), { id: "td", frameParameters: { day: dayNumber() }, spec: tdSpec, colorLimits: [+$2("low").value, +$2("high").value], onSelect: ({ panel, x }) => {
        const time = axisStart() + x * 36e8;
        selected = frames.reduce((best, f, i) => Math.abs((f.start_us + f.stop_us) / 2 - time) < Math.abs((frames[best].start_us + frames[best].stop_us) / 2 - time) ? i : best, 0);
        $2("rdbeam").value = String(panel);
        $2("frame").value = String(selected);
        renderRD();
      } });
    } else tdPlot.setSpec(tdSpec);
    tdPlot.setFrameParameters({ day: dayNumber() });
  };
  for (const [root, mount, ready2] of [[$2("rd"), mountRD, rdPlot], [$2("timeDoppler"), mountTD, tdPlot]]) {
    secondaryMounts.set(root, mount);
    if (ready2) await mount();
    else visiblePlots.observe(root);
  }
}
function scale(p) {
  const [low, high] = p.owner.getColorLimits();
  p.view.set_scale(low, high, 0, 0.1, 2);
}
function summaries(frame, beam, nr) {
  const cached = summaryCache.get(frame);
  if (cached?.[beam]) return cached[beam];
  const peak = new Float32Array(nr), sum = new Float64Array(nr), vsum = new Float64Array(nr), v2sum = new Float64Array(nr), vel = manifest.geometry.velocity_mps;
  const nyquist = manifest.velocity_limits_mps[1];
  frame.index.forEach(([b, f, r], i) => {
    if (b !== beam) return;
    const w = 10 ** (frame.snr_db[i] / 10), angle = Math.PI * vel[f] / nyquist;
    peak[r] = Math.max(peak[r], frame.snr_db[i]);
    sum[r] += w;
    vsum[r] += w * Math.cos(angle);
    v2sum[r] += w * Math.sin(angle);
  });
  const moments = new Float32Array(nr * 4), variance = new Float64Array(nr);
  for (let r = 0; r < nr; r++) {
    moments[4 * r] = peak[r];
    if (sum[r]) moments[4 * r + 1] = Math.atan2(v2sum[r], vsum[r]) * nyquist / Math.PI;
  }
  frame.index.forEach(([b, f, r], i) => {
    if (b !== beam) return;
    const w = 10 ** (frame.snr_db[i] / 10), period = 2 * nyquist, dv = ((vel[f] - moments[4 * r + 1] + nyquist) % period + period) % period - nyquist;
    variance[r] += w * dv * dv;
  });
  for (let r = 0; r < nr; r++) if (sum[r]) moments[4 * r + 2] = Math.sqrt(variance[r] / sum[r]);
  const result = { peak, moments }, cache = cached || [];
  cache[beam] = result;
  summaryCache.set(frame, cache);
  return result;
}
function renderRTI() {
  if (!rtiPlot || !frames.length) return;
  const style = { day: dayNumber(), hsv: $2("display").value === "hsv" ? 1 : 0, velocity_limit: +$2("vlim").value, width_limit: +$2("wlim").value };
  rtiPlot.setFrameParameters(style);
  tdPlot?.setFrameParameters(style);
  const geom = manifest.geometry, nr = geom.slant_range_km.length, first = axisStart(), columns = [], xb = [];
  for (let i = 0; i < frames.length; i++) {
    const f = frames[i], previous = frames[i - 1], next = frames[i + 1];
    const left = previous && f.start_us - previous.stop_us <= 6e8 ? (previous.stop_us + f.start_us) / 2 : f.start_us;
    const right = next && next.start_us - f.stop_us <= 6e8 ? (f.stop_us + next.start_us) / 2 : f.stop_us;
    if (xb.length && left > xb.at(-1) + 1) {
      columns.push(null);
      xb.push(left);
    }
    if (!xb.length) xb.push(left);
    columns.push(f);
    xb.push(right);
  }
  const x = Float32Array.from(xb, (t) => (t - first) / 36e8), nx = columns.length;
  const nativeNf = geom.velocity_mps.length, nf = nativeNf + 1;
  tdPlot?.panels.forEach((panel, beam) => {
    const value = new Float32Array(nx * nf), moment = new Float32Array(nx * nf * 4);
    columns.forEach((frame, c) => {
      if (!frame) {
        for (let f = 0; f < nf; f++) {
          value[f * nx + c] = NaN;
          moment[4 * (f * nx + c)] = NaN;
        }
        return;
      }
      const summary = summaries(frame, beam, nr);
      frame.index.forEach(([b, f, r], i) => {
        if (b === beam && geom.altitude_km[beam][r] >= 50 && geom.altitude_km[beam][r] <= 120) {
          const cell = f * nx + c, snr = frame.snr_db[i];
          if (snr > value[cell]) {
            value[cell] = snr;
            moment[4 * cell] = snr;
            moment[4 * cell + 1] = geom.velocity_mps[f];
            moment[4 * cell + 2] = summary.moments[4 * r + 2];
            if (f === 0) {
              const last = nativeNf * nx + c;
              value[last] = snr;
              moment.set(moment.subarray(4 * cell, 4 * cell + 4), 4 * last);
              moment[4 * last + 1] = manifest.velocity_limits_mps[1];
            }
          }
        }
      });
    });
    panel.set_rectilinear_data(x, edges2([...geom.velocity_mps, manifest.velocity_limits_mps[1]]), value);
    panel.set_scale(...tdPlot.getColorLimits());
    if (style.hsv) panel.set_moment_data(moment, ...tdPlot.getHsvLimits() || [+$2("vlim").value, +$2("wlim").value]);
  });
  plots2.forEach((p, beam) => {
    const val = new Float32Array(nx * nr), moment = new Float32Array(nx * nr * 4);
    columns.forEach((frame, c) => {
      const s = frame ? summaries(frame, beam, nr) : null;
      for (let r = 0; r < nr; r++) {
        val[r * nx + c] = s ? s.peak[r] : NaN;
        for (let k = 0; k < 4; k++) moment[4 * (r * nx + c) + k] = s ? s.moments[4 * r + k] : k === 0 ? NaN : 0;
      }
    });
    p.view.set_rectilinear_data(x, edges2(geom.altitude_km[beam]), val);
    scale(p);
    if ($2("display").value === "hsv") p.view.set_moment_data(moment, ...rtiPlot.getHsvLimits() || [+$2("vlim").value, +$2("wlim").value]);
  });
  const hsv = $2("display").value === "hsv";
  $2("colorLabel").textContent = hsv ? "Hue: radial velocity; saturation: spectral width; brightness: dB SNR. Click the HSV key to adjust." : "Peak excess signal/noise (dB). Click a colorbar to adjust its limits.";
}
function renderRD() {
  if (!rd || !frames.length) return;
  const frame = frames[selected], beam = +$2("rdbeam").value, g = manifest.geometry, nr = g.slant_range_km.length;
  const bins = g.velocity_mps.map((v, i) => i), nx = bins.length + 1;
  const lookup = new Map(bins.map((f, i) => [f, i])), value = new Float32Array(nx * nr), moment = new Float32Array(nx * nr * 4), summary = summaries(frame, beam, nr);
  let count = 0, transients = 0;
  frame.index.forEach(([b, f, r], i) => {
    if (b === beam && lookup.has(f)) {
      const cell = r * nx + lookup.get(f);
      value[cell] = frame.snr_db[i];
      moment[4 * cell] = frame.snr_db[i];
      moment[4 * cell + 1] = g.velocity_mps[f];
      moment[4 * cell + 2] = summary.moments[4 * r + 2];
      if (f === 0) {
        const last = r * nx + nx - 1;
        value[last] = value[cell];
        moment.set(moment.subarray(4 * cell, 4 * cell + 4), 4 * last);
        moment[4 * last + 1] = manifest.velocity_limits_mps[1];
      }
      count++;
      if (frame.look_dominance?.[i] > 0.5) transients++;
    }
  });
  rd.view.set_rectilinear_data(edges2([...g.velocity_mps, manifest.velocity_limits_mps[1]]), edges2(g.altitude_km[beam]), value);
  scale(rd);
  const hsv = $2("display").value === "hsv";
  rdPlot.setFrameParameters({ day: dayNumber(), hsv: hsv ? 1 : 0, velocity_limit: +$2("vlim").value, width_limit: +$2("wlim").value });
  if (hsv) rd.view.set_moment_data(moment, ...rdPlot.getHsvLimits() || [+$2("vlim").value, +$2("wlim").value]);
  $2("frameStatus").textContent = `${utc(frame.start_us)} \xB7 ${names2[beam]} \xB7 ${frame.duration_s.toFixed(2)} s mode, ${frame.coverage_s.toFixed(2)} s valid support \xB7 ${frame.looks} looks \xB7 ${count} retained cells \xB7 16 ms complementary pairs \xB7 ${frame.quality}`;
}
async function refresh({ useIndex = false } = {}) {
  const request = ++refreshSequence;
  const update = (stage, progress) => {
    if (request === refreshSequence) loading(stage, progress);
  };
  update("Downloading observations\u2026");
  let m = catalogue;
  if (!useIndex || !m) {
    const response = initialResponse ? await initialResponse : await fetchData("./data/manifest.json");
    initialResponse = null;
    if (response.error) throw response.error;
    m = await readJSON(response, (p) => update("Downloading observations\u2026", p));
    catalogue = m;
  }
  if (request !== refreshSequence) return;
  renderProcessingStatus(m);
  const days = m.available_days || [];
  if (!chosenDay || !days.includes(chosenDay)) chosenDay = days.at(-1) || m.selected_day;
  const control = $2("day");
  control.replaceChildren(...days.map((day) => {
    const option = document.createElement("option");
    option.value = day;
    option.textContent = day;
    return option;
  }));
  control.value = chosenDay;
  if (skyNoise) renderNoise();
  else if (m.noise_model_status) renderNoise(m);
  if ($2("skyNoiseSection").open) void refreshSkyNoise();
  if (m.transfer) {
    const a = m.transfer.archive, t = m.transfer.syowa_to_revontuli;
    $2("transfer").textContent = `Delivered native archive: ${(a.bytes / 2 ** 20).toFixed(2)} MiB \xB7 ${a.records} native records \xB7 ${a.files} HDF5 files. ` + (a.dense_complex_sample_retained_fraction != null ? `Complex spectra retain ${(100 * a.dense_complex_sample_retained_fraction).toPrecision(3)}% of the full dense spectra (${a.retained_complex_samples.toLocaleString()} / ${a.dense_complex_samples.toLocaleString()} complex samples; all beams, Doppler bins, range gates, receiver channels and coherent looks, before altitude, velocity and detection masks). ` : "Complex-spectra retained fraction: awaiting complete dense-spectrum counts. ") + "Storage: complex64 spectra (8 bytes/sample); float32 noise power and uncertainty (4 bytes/value), with lossless compression. " + (a.first_mode_us ? `${utc(a.first_mode_us)} through ${utc(a.latest_mode_us)}. ` : "") + (a.projected_bytes_per_day === null ? "Daily growth: interval too short for extrapolation. " : `Provisional growth: ${(a.projected_bytes_per_day / 2 ** 20).toFixed(1)} MiB/day (${(a.projected_bytes_per_day / 86400 / 1024).toFixed(2)} KiB/s average). `) + (t ? `Latest Syowa replication: ${(t.literal_bytes / 2 ** 20).toFixed(2)} MiB in ${t.elapsed_seconds.toFixed(1)} s (${(t.payload_bytes_per_second / 1024).toFixed(2)} KiB/s payload; cap ${t.bandwidth_cap_kib_s} KiB/s). ` : "Transfer measurement pending. ") + "Web copy: seven-day retention; Revontuli archive retained.";
  }
  renderArchivePlot(m);
  if (m.selected_day !== chosenDay || !m.frame_packets) {
    update(`Downloading ${chosenDay}\u2026`);
    let daily = await fetchData(`./data/days/${encodeURIComponent(chosenDay)}.json`);
    if (daily.status === 404) daily = await fetchData(`./api/day?day=${encodeURIComponent(chosenDay)}`);
    m = await readJSON(daily, (p) => update(`Downloading ${chosenDay}\u2026`, p));
  }
  if (request !== refreshSequence) return;
  saveView(m);
  const age = (Date.now() - m.latest_mode_us / 1e3) / 6e4;
  $2("status").textContent = m.frame_count ? `${m.frame_count} observing modes \xB7 latest ${utc(m.latest_mode_us)} \xB7 ${Math.max(0, age).toFixed(0)} min behind now` : `${chosenDay} UTC \xB7 no observations delivered yet`;
  if (m.spectra_transfer_pending) $2("status").textContent += ` \xB7 complex spectra transferring for ${m.spectra_transfer_pending} modes`;
  const download = $2("downloadSpectra");
  download.href = `./api/download?day=${encodeURIComponent(chosenDay)}`;
  download.setAttribute("download", `pansy-sparse-spectra-${chosenDay}.tar`);
  download.setAttribute("aria-disabled", m.frame_count ? "false" : "true");
  download.title = `Native complex spectra and metadata for ${chosenDay} UTC (.tar)`;
  $2("status").className = age > 10 ? "status warn" : "status";
  $2("error").textContent = "";
  if (m.revision === revision) {
    loaded();
    return;
  }
  manifest = m;
  update("Preparing echo data\u2026");
  originalFrames = m.frame_packets;
  applyEchoFilter();
  revision = m.revision;
  selected = frames.length - 1;
  $2("frame").max = String(Math.max(0, selected));
  $2("frame").value = String(Math.max(0, selected));
  if (!frames.length) {
    $2("status").textContent = `${chosenDay} UTC \xB7 no complete 16 ms records yet`;
    loaded();
    return;
  }
  update("Drawing Rustplotlib axes\u2026");
  await new Promise((resolve) => requestAnimationFrame(resolve));
  await mountPlots();
  if (request !== refreshSequence) return;
  update("Drawing echoes\u2026");
  $2("sidecar").href = "./data/" + m.sidecar;
  renderRTI();
  renderRD();
  await rtiPlot.requestFrame();
  if (request !== refreshSequence) return;
  requestAnimationFrame(() => {
    if (request === refreshSequence) {
      const now = String(Math.round(performance.now()));
      if (!document.body.dataset.plotsReadyMs) document.body.dataset.plotsReadyMs = now;
      document.body.dataset.latestRefreshMs = now;
      loaded();
    }
  });
}
function renderArchivePlot(m) {
  archiveInventory = m.transfer?.archive;
  if ($2("archiveSection").open) void drawArchive();
}
async function drawArchive() {
  try {
    await renderArchiveDiagnostics(archiveInventory);
    $2("archivePlotStatus").textContent = "Numeric inventory rendered in this browser by Rustplotlib. Zoom and pan each plot; no plot images are downloaded.";
  } catch (error) {
    $2("archivePlotStatus").textContent = String(error);
  }
}
function renderProcessingStatus(m) {
  const s = m.source_processing, node = $2("processingStatus");
  node.hidden = !s;
  if (!s) return;
  const b = s.backfill || {}, live = s.live || {}, parts = [];
  if (live.last_written_stop_us) parts.push(`Source live processing (${live.state}): ${utc(live.last_written_stop_us)}.`);
  if (b.current_index && b.total_modes) {
    parts.push(`Rectangular backfill (${b.complete ? "complete" : b.state}): ${b.current_index} / ${b.total_modes} modes; working on ${utc(b.current_start_us)}.`);
    const frontier = new Date(b.current_start_us / 1e3).toISOString().slice(0, 10);
    const delivered = Object.entries(m.transfer?.archive?.days || {}).filter(([d, v]) => d <= frontier && v.latest_mode_us).sort(([a], [b2]) => a.localeCompare(b2)).at(-1);
    parts.push(delivered ? `Historical voltages delivered to viewer through ${utc(delivered[1].latest_mode_us)}.` : "Historical voltage delivery pending.");
  }
  const age = Math.max(0, (Date.now() / 1e3 - s.generated_unix_s) / 60);
  parts.push(`Processing status updated ${age.toFixed(0)} min ago; processing and transfer are separate.`);
  if (m.plot_generation_status === "delayed") parts.push("Optional plot refresh delayed; native observations continue updating.");
  node.textContent = parts.join(" ");
}
async function refreshSkyNoise() {
  if (noiseRequest) return noiseRequest;
  noiseRequest = (async () => {
    try {
      const response = await fetchData("./data/sky-noise/status.json");
      if (!response.ok) throw Error(`Sky-noise diagnostics unavailable (${response.status})`);
      const value = await response.json();
      if (value.noise_model_status?.length !== 7 || !value.noise_model_data?.times_h?.length) throw Error("Numeric sky-noise data have not arrived yet");
      skyNoise = value;
      renderNoise();
    } catch (error) {
      if (skyNoise) {
        skyNoise.noise_plot = { ...skyNoise.noise_plot, status: "delayed" };
        renderNoise();
      } else {
        $2("noiseStatus").textContent = "Numeric sky-noise data are not available yet; the independent analysis job will retry. Observations remain available.";
      }
    }
  })();
  try {
    await noiseRequest;
  } finally {
    noiseRequest = null;
  }
}
function renderNoise(model = skyNoise || manifest) {
  const module = +$2("noiseModule").value, fit = model?.noise_model_status?.[module];
  if (!fit) {
    void refreshSkyNoise();
    return;
  }
  $2("noiseStatus").textContent = `PyGDSM GSM2008 \xB7 one-way receive pattern \xB7 90% crossed-Yagi / 10% isotropic element. Module ${module}: ` + (fit.receiver_temp_k === null ? fit.status : `T_receiver = ${fit.receiver_temp_k.toFixed(0)} K; fractional RMS ${(fit.fractional_rms * 100).toFixed(1)}%; ${fit.status}`) + ". Top: measured noise power; bottom: modeled sky temperature. Fits update hourly and require \u226512 hours of data.";
  if (model.noise_model_observation_stop_unix_s) $2("noiseStatus").textContent += ` Measurements through ${utc(model.noise_model_observation_stop_unix_s * 1e6)}.`;
  if (model.noise_plot?.status === "delayed") $2("noiseStatus").textContent += " Refresh delayed; using cached numeric data.";
  if ($2("skyNoiseSection").open) void renderNoiseDiagnostics(model, module).catch((error) => {
    $2("noiseStatus").textContent += " " + String(error);
  });
}
function failure(e) {
  revision = "";
  loading("Loading failed \u2014 retrying\u2026");
  $2("error").textContent = e.message;
  console.error(e);
}
try {
  loading("Loading plotting engine\u2026");
  await init3();
  if (new URLSearchParams(location.search).has("embed")) document.body.classList.add("embedded");
  document.body.dataset.display = $2("display").value;
  for (const id of ["display", "low", "high", "vlim", "wlim"]) $2(id).onchange = () => {
    try {
      document.body.dataset.display = $2("display").value;
      if (id === "low" || id === "high") for (const plot2 of [rtiPlot, rdPlot, tdPlot]) plot2?.setColorLimits(+$2("low").value, +$2("high").value);
      if (id === "vlim" || id === "wlim") for (const plot2 of [rtiPlot, rdPlot, tdPlot]) plot2?.setHsvLimits(+$2("vlim").value, +$2("wlim").value);
      renderRTI();
      renderRD();
    } catch (e) {
      failure(e);
    }
  };
  $2("rdbeam").onchange = renderRD;
  $2("frame").oninput = () => {
    selected = +$2("frame").value;
    renderRD();
  };
  $2("hideSingleMode").onchange = () => {
    applyEchoFilter();
    renderRTI();
    renderRD();
  };
  $2("noiseModule").onchange = () => renderNoise();
  $2("skyNoiseSection").ontoggle = () => {
    if ($2("skyNoiseSection").open) void refreshSkyNoise();
  };
  $2("archiveSection").ontoggle = () => {
    if ($2("archiveSection").open) void drawArchive();
  };
  $2("day").onchange = () => {
    showLoading = true;
    loading("Loading selected day\u2026");
    chosenDay = $2("day").value;
    revision = "";
    $2("downloadSpectra").setAttribute("aria-disabled", "true");
    $2("downloadSpectra").removeAttribute("href");
    if (playing) {
      clearInterval(playing);
      playing = null;
      $2("play").textContent = "Play";
    }
    refresh({ useIndex: true }).catch(failure);
  };
  $2("play").onclick = () => {
    if (playing) {
      clearInterval(playing);
      playing = null;
      $2("play").textContent = "Play";
    } else {
      $2("play").textContent = "Pause";
      playing = setInterval(() => {
        selected = (selected + 1) % frames.length;
        $2("frame").value = selected;
        renderRD();
      }, 500);
    }
  };
  let resizeTimer;
  window.addEventListener("resize", () => {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(() => {
      if (frames.length) {
        renderRTI();
        renderRD();
      }
    }, 150);
  });
  async function pollObservations() {
    let delay = 6e4;
    try {
      await refresh();
    } catch (e) {
      failure(e);
      delay = 5e3;
    } finally {
      setTimeout(pollObservations, delay);
    }
  }
  await restoreView();
  await pollObservations();
} catch (e) {
  failure(e);
}
