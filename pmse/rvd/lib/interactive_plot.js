import init, { FieldView, PlotViewport } from "../pkg/releases/c70022f0996b1a74/rustplotlib_web.js";
let initialization;
const clamp = (v) => Math.max(0, Math.min(1, v));
const svgURL = (svg) => URL.createObjectURL(new Blob([svg], { type: "image/svg+xml" }));
class InteractiveFieldPlot {
  constructor(root, options) {
    this.root = root;
    this.options = options;
    const s = options.spec;
    if (!s.xlabel?.trim() || !s.ylabel?.trim() || !s.colorbar_label?.trim()) throw Error("Interactive plots require x, y and colorbar labels");
    if (!s.panels.length || s.bounds.length !== 4 || !s.bounds.every(Number.isFinite)) throw Error("Invalid interactive figure specification");
    this.onSelect = options.onSelect;
    this.limits = [...options.colorLimits || [0, 1]];
    this.initialLimits = [...this.limits];
    this.momentMode = options.frameParameters?.hsv === 1;
    this.initialMomentLimits = [options.frameParameters?.velocity_limit || 30, options.frameParameters?.width_limit || 10];
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
  static async mount(root, options) {
    initialization ??= init();
    await initialization;
    const plot = new InteractiveFieldPlot(root, options);
    await plot.requestFrame();
    plot.observer.observe(root);
    return plot;
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
    const endpoint = this.options.frameEndpoint || new URL("./api/frame", location.href).href;
    const response = await fetch(endpoint, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      signal: this.abort.signal,
      body: JSON.stringify({ ...this.options.frameParameters, ...this.momentLimits ? { velocity_limit: this.momentLimits[0], width_limit: this.momentLimits[1] } : {}, hsv: this.momentMode ? 1 : 0, id: this.options.id, width, viewport: Array.from(this.viewport.bounds()), color_limits: this.limits })
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
}
class InteractiveFieldPanel {
  constructor(plot, canvas, label, index) {
    this.plot = plot;
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
}
export {
  InteractiveFieldPanel,
  InteractiveFieldPlot
};
