/* tslint:disable */
/* eslint-disable */

export class BrowserImager {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Spectra are complex64 [look,7,cell], transported without phase changes.
     * Output [pixel,4] = relative dB power, circular velocity, width, support.
     */
    image(header_json: string, geometry_json: string, values: Float32Array): Float32Array;
    /**
     * Default is positive Haar-wavelet analysis L1. Legacy image() retains
     * its historical normalization and TSVD reconstruction.
     */
    image_with_options(header_json: string, geometry_json: string, values: Float32Array, options_json: string): Float32Array;
    last_diagnostics(): string;
    constructor(geometry_json: string, height_km: number, grid: number, fixed_slant_range: boolean);
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_browserimager_free: (a: number, b: number) => void;
    readonly browserimager_image: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number, number, number];
    readonly browserimager_image_with_options: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number) => [number, number, number, number];
    readonly browserimager_last_diagnostics: (a: number) => [number, number, number, number];
    readonly browserimager_new: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly pansy_coherent: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => number;
    readonly pansy_complementary: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => number;
    readonly pansy_decode_selected: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number) => number;
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
