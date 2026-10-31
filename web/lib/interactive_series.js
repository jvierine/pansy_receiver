// ../rustplotlib-monitor-mobile.QBv0hs/web/src/interactive_series.ts
import init, { SeriesView, PlotViewport } from "../pkg/releases/7caedeeb810f57e4/rustplotlib_web.js";

// ../rustplotlib-monitor-mobile.QBv0hs/web/src/plot_controls.ts
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

// ../rustplotlib-monitor-mobile.QBv0hs/web/src/interactive_series.ts
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
export {
  InteractiveSeriesPlot,
  PlotControlGroup
};
