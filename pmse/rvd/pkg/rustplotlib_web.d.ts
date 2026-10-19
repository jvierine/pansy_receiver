/* tslint:disable */
/* eslint-disable */

/**
 * Upload a rectangular field once, then change scales/views without re-upload.
 */
export class FieldView {
    free(): void;
    [Symbol.dispose](): void;
    draw(): void;
    constructor(canvas: HTMLCanvasElement);
    pixel(x: number, y: number): Uint8Array;
    renderer(): string;
    set_colormap(rgb: Float32Array): void;
    set_data(nx: number, ny: number, values: Float32Array): void;
    /**
     * Missing observations can be distinguished from values below the scale.
     */
    set_missing_color(red: number, green: number, blue: number, alpha: number): void;
    /**
     * Read a rendered pixel (bottom-left origin), useful for renderer checks.
     * Replace an already-sized field with [SNR dB, velocity, width, reserved]
     * cells. Existing rectilinear edges and view are preserved. HSV conversion
     * is performed by the WebGL shader, not by the host language.
     */
    set_moment_data(values: Float32Array, velocity_limit: number, width_limit: number): void;
    set_rectilinear_data(x_edges: Float32Array, y_edges: Float32Array, values: Float32Array): void;
    /**
     * mode: 0 linear, 1 logarithmic, 2 signed asinh; palette: 0 diverging, 1 sequential, 2 ember.
     */
    set_scale(low: number, high: number, mode: number, width: number, palette: number): void;
    set_view(x0: number, y0: number, x1: number, y1: number): void;
}

/**
 * Geographic RGBA layers on a rotatable globe. Input rows run south to north
 * and columns run west to east from -180 to 180 degrees.
 */
export class GlobeView {
    free(): void;
    [Symbol.dispose](): void;
    draw(): void;
    constructor(canvas: HTMLCanvasElement);
    set_camera(yaw: number, pitch: number): void;
    set_rgba(width: number, height: number, rgba: Uint8Array): void;
}

/**
 * Colored 3-D line geometry, e.g. traced vector fields. Coordinates stay in
 * physical units; the camera centre and extent use those same units.
 */
export class PathView {
    free(): void;
    [Symbol.dispose](): void;
    draw(): void;
    constructor(canvas: HTMLCanvasElement);
    set_camera(yaw: number, pitch: number, extent: number): void;
    set_paths(xyz_rgb: Float32Array): void;
}

export class PlotViewport {
    free(): void;
    [Symbol.dispose](): void;
    bounds(): Float64Array;
    constructor();
    /**
     * Pan by fractions of the currently visible extent.
     */
    pan(dx: number, dy: number): void;
    reset(): void;
    /**
     * Zoom to a dragged rectangle in current viewport coordinates (lower-left origin).
     * Reversed drags are accepted; empty or non-finite boxes leave the view unchanged.
     */
    zoom_box(x0: number, y0: number, x1: number, y1: number): boolean;
    /**
     * Zoom around pointer coordinates measured from the lower-left of the axes.
     */
    zoom(x: number, y: number, factor: number): void;
}

export class SeriesView {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * CSS dimensions and device pixel ratio are separate, preserving font size.
     */
    begin(width: number, height: number, dpr: number, x0: number, x1: number, y0: number, y1: number, left: string, right: string): void;
    constructor(canvas: HTMLCanvasElement);
    /**
     * Paired x,y samples; nonfinite samples break lines. Markers stay in CSS px.
     */
    series(xy: Float64Array, color: string, line: boolean, radius: number): void;
}

export class VolumeView {
    free(): void;
    [Symbol.dispose](): void;
    draw(): void;
    constructor(canvas: HTMLCanvasElement);
    renderer(): string;
    /**
     * Batch data, transfer and camera changes before one GPU draw.
     */
    set_auto_draw(enabled: boolean): void;
    /**
     * Preserve physical box proportions instead of stretching data to a cube.
     */
    set_box_extent(x: number, y: number, z: number): void;
    set_camera(yaw: number, pitch: number): void;
    set_colormap(rgb: Float32Array): void;
    set_data(nx: number, ny: number, nz: number, values: Float32Array): void;
    /**
     * Render a spherical atmosphere with a quadrant removed. Radii are in
     * the unit texture cube (maximum radius 0.5); data remain Cartesian.
     */
    set_spherical_shell(inner: number, outer: number, enabled: boolean): void;
    set_transfer(low: number, high: number, density: number, cut: number): void;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_fieldview_free: (a: number, b: number) => void;
    readonly __wbg_globeview_free: (a: number, b: number) => void;
    readonly __wbg_pathview_free: (a: number, b: number) => void;
    readonly __wbg_plotviewport_free: (a: number, b: number) => void;
    readonly __wbg_seriesview_free: (a: number, b: number) => void;
    readonly __wbg_volumeview_free: (a: number, b: number) => void;
    readonly fieldview_draw: (a: number) => void;
    readonly fieldview_new: (a: any) => [number, number, number];
    readonly fieldview_pixel: (a: number, b: number, c: number) => [number, number, number, number];
    readonly fieldview_renderer: (a: number) => [number, number];
    readonly fieldview_set_colormap: (a: number, b: number, c: number) => [number, number];
    readonly fieldview_set_data: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly fieldview_set_missing_color: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly fieldview_set_moment_data: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly fieldview_set_rectilinear_data: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number];
    readonly fieldview_set_scale: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly fieldview_set_view: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly globeview_draw: (a: number) => void;
    readonly globeview_new: (a: any) => [number, number, number];
    readonly globeview_set_camera: (a: number, b: number, c: number) => [number, number];
    readonly globeview_set_rgba: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly pathview_draw: (a: number) => void;
    readonly pathview_new: (a: any) => [number, number, number];
    readonly pathview_set_camera: (a: number, b: number, c: number, d: number) => [number, number];
    readonly pathview_set_paths: (a: number, b: number, c: number) => [number, number];
    readonly plotviewport_bounds: (a: number) => [number, number];
    readonly plotviewport_new: () => number;
    readonly plotviewport_pan: (a: number, b: number, c: number) => void;
    readonly plotviewport_reset: (a: number) => void;
    readonly plotviewport_zoom: (a: number, b: number, c: number, d: number) => void;
    readonly plotviewport_zoom_box: (a: number, b: number, c: number, d: number, e: number) => number;
    readonly seriesview_begin: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number) => [number, number];
    readonly seriesview_new: (a: any) => [number, number, number];
    readonly seriesview_series: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number];
    readonly volumeview_draw: (a: number) => void;
    readonly volumeview_new: (a: any) => [number, number, number];
    readonly volumeview_renderer: (a: number) => [number, number];
    readonly volumeview_set_auto_draw: (a: number, b: number) => void;
    readonly volumeview_set_box_extent: (a: number, b: number, c: number, d: number) => [number, number];
    readonly volumeview_set_camera: (a: number, b: number, c: number) => void;
    readonly volumeview_set_colormap: (a: number, b: number, c: number) => [number, number];
    readonly volumeview_set_data: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly volumeview_set_spherical_shell: (a: number, b: number, c: number, d: number) => [number, number];
    readonly volumeview_set_transfer: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly __wbindgen_exn_store: (a: number) => void;
    readonly __externref_table_alloc: () => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
