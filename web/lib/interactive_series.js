import init, { SeriesView, PlotViewport } from "../pkg/rustplotlib_web.js";
let ready;
class InteractiveSeriesPlot {
  constructor(root, spec, height, compact = false) {
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
    root.replaceChildren(...compact ? [stage] : [toolbar, stage]);
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
  viewport = new PlotViewport();
  canvas = document.createElement("canvas");
  view;
  observer;
  rect = [0, 0, 1, 1];
  box = document.createElement("div");
  drag;
  pending = 0;
  specJSON;
  static async mount(root, spec, options = {}) {
    ready ??= init().catch((error) => {
      ready = void 0;
      throw error;
    });
    await ready;
    return new InteractiveSeriesPlot(root, spec, options.height || 300, options.compact);
  }
  setSpec(spec) {
    const next = JSON.stringify(spec);
    if (next === this.specJSON) return;
    this.spec = spec;
    this.specJSON = next;
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
}
export {
  InteractiveSeriesPlot
};
