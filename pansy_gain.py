"""Reusable PANSY transmit and receive antenna gain models.

The PANSY transmitter is a sparse array of antenna modules.  Each module is a
19-element subarray, so the modeled transmit power gain is the product of the
single-module subarray power pattern and the sparse module-center array-factor
power pattern.  The interferometric receive channels use individual modules,
so their receive gain is the same single-module power pattern.
"""

from __future__ import annotations

import numpy as np

import pansy_config as pc
import pansy_modes as pmm

MIN_FIELD_AMPLITUDE = 1e-8
MIN_POWER_GAIN = MIN_FIELD_AMPLITUDE**2


def _is_ready(value) -> bool:
    try:
        return int(value) == 1
    except ValueError:
        return str(value).strip().lower() in {"true", "yes", "ready"}


def tx_beam_unit_vectors() -> np.ndarray:
    """Transmit beam pointing vectors in the interferometry u/v/w convention."""
    mode = pmm.get_m_mode()
    vecs = []
    for az_deg, za_deg in mode["beam_pos_az_za"]:
        el_deg = 90.0 - za_deg
        p_h = np.cos(np.deg2rad(el_deg))
        w = -np.sin(np.deg2rad(el_deg))
        v = p_h * np.cos(-np.deg2rad(az_deg))
        u = -p_h * np.sin(-np.deg2rad(az_deg))
        vecs.append([u, v, w])
    return np.asarray(vecs, dtype=np.float64)


def tx_module_positions() -> list[np.ndarray]:
    """Ready transmit antenna positions grouped by physical module."""
    groups = []
    for name in sorted(pc.modules):
        pos = []
        for ant in pc.antenna.values():
            if ant["name"] != name or ant["serial"] == "RFTX":
                continue
            if _is_ready(ant["ready"]):
                pos.append([ant["x"], ant["y"], ant["z"]])
        if pos:
            groups.append(np.asarray(pos, dtype=np.float64))
    if not groups:
        raise RuntimeError("no ready transmit antenna modules found")
    return groups


def tx_module_center_positions(modules: list[np.ndarray] | None = None) -> np.ndarray:
    """Ready transmit module centers with the array phase origin removed."""
    modules = tx_module_positions() if modules is None else modules
    centers = np.asarray([np.mean(pos, axis=0) for pos in modules], dtype=np.float64)
    return centers - np.mean(centers, axis=0, keepdims=True)


def rx_module_positions(channel: int | str | None = None) -> list[np.ndarray] | np.ndarray:
    """Receive module antenna positions relative to each module center.

    Parameters
    ----------
    channel
        If ``None``, return a list for all connected interferometric receive
        modules.  If an integer, use that receiver channel index from
        ``cfg/connections.txt``.  If a string, use that module name directly.
    """
    names = [name for name in pc.connections if name != "RFTX" and name in pc.modules]
    if channel is None:
        return [np.asarray(pc.modules[name], dtype=np.float64) - pc.module_center[name] for name in names]
    if isinstance(channel, int):
        name = pc.connections[channel]
    else:
        name = str(channel)
    if name not in pc.modules:
        raise KeyError(f"unknown PANSY receive module {name!r}")
    return np.asarray(pc.modules[name], dtype=np.float64) - pc.module_center[name]


def _as_direction_array(uvw: np.ndarray) -> tuple[np.ndarray, tuple[int, ...]]:
    uvw = np.asarray(uvw, dtype=np.float64)
    if uvw.shape[-1] != 3:
        raise ValueError("uvw must have shape (..., 3)")
    original_shape = uvw.shape[:-1]
    return uvw.reshape(-1, 3), original_shape


def module_field_pattern(
    uvw: np.ndarray,
    steer: np.ndarray | None = None,
    module_positions: np.ndarray | None = None,
) -> np.ndarray:
    """Complex single-module subarray voltage pattern.

    The returned field is normalized to one at ``steer``.  Use
    :func:`module_power_gain` for a power gain suitable for radar-equation/RCS
    work.
    """
    dirs, original_shape = _as_direction_array(uvw)
    steer = np.asarray([0.0, 0.0, -1.0] if steer is None else steer, dtype=np.float64)
    if module_positions is None:
        rx_modules = rx_module_positions()
        module_positions = rx_modules[0]
    rel_pos = np.asarray(module_positions, dtype=np.float64)
    rel_pos = rel_pos - np.mean(rel_pos, axis=0, keepdims=True)
    out = np.full(len(dirs), np.nan + 0j, dtype=np.complex128)
    good = np.all(np.isfinite(dirs), axis=1)
    if np.any(good):
        delta = dirs[good] - steer
        phase = (2.0 * np.pi / pc.wavelength) * (rel_pos @ delta.T)
        out[good] = np.mean(np.exp(1j * phase), axis=0)
    return out.reshape(original_shape)


def module_power_gain(
    uvw: np.ndarray,
    steer: np.ndarray | None = None,
    module_positions: np.ndarray | None = None,
) -> np.ndarray:
    """Single-module subarray receive/transmit power gain, linear units."""
    field = module_field_pattern(uvw, steer=steer, module_positions=module_positions)
    return np.maximum(np.abs(field) ** 2, MIN_POWER_GAIN)


def sparse_tx_array_field(
    uvw: np.ndarray,
    beam_id: int | np.ndarray,
    beam_vecs: np.ndarray | None = None,
    module_centers: np.ndarray | None = None,
) -> np.ndarray:
    """Sparse module-center transmit array voltage pattern."""
    dirs, original_shape = _as_direction_array(uvw)
    beam_vecs = tx_beam_unit_vectors() if beam_vecs is None else np.asarray(beam_vecs, dtype=np.float64)
    module_centers = tx_module_center_positions() if module_centers is None else np.asarray(module_centers, dtype=np.float64)
    beam_arr = np.asarray(beam_id, dtype=np.int64)
    if beam_arr.ndim == 0:
        beam_arr = np.full(len(dirs), int(beam_arr), dtype=np.int64)
    else:
        beam_arr = np.broadcast_to(beam_arr, original_shape).reshape(-1)
    out = np.full(len(dirs), np.nan + 0j, dtype=np.complex128)
    good = np.all(np.isfinite(dirs), axis=1) & (beam_arr >= 0) & (beam_arr < len(beam_vecs))
    if np.any(good):
        for beam in np.unique(beam_arr[good]):
            idx = good & (beam_arr == beam)
            delta = dirs[idx] - beam_vecs[beam]
            phase = (2.0 * np.pi / pc.wavelength) * (module_centers @ delta.T)
            out[idx] = np.mean(np.exp(1j * phase), axis=0)
    return out.reshape(original_shape)


def sparse_tx_array_power_gain(
    uvw: np.ndarray,
    beam_id: int | np.ndarray,
    beam_vecs: np.ndarray | None = None,
    module_centers: np.ndarray | None = None,
) -> np.ndarray:
    """Sparse module-center transmit array-factor power gain, linear units."""
    field = sparse_tx_array_field(uvw, beam_id, beam_vecs=beam_vecs, module_centers=module_centers)
    return np.maximum(np.abs(field) ** 2, MIN_POWER_GAIN)


def tx_power_gain(
    uvw: np.ndarray,
    beam_id: int | np.ndarray,
    beam_vecs: np.ndarray | None = None,
    modules: list[np.ndarray] | None = None,
) -> np.ndarray:
    """Full transmit power gain in linear units.

    ``tx_gain = single_module_gain * sparse_transmit_array_gain``.
    This assumes all transmit modules share the same subarray geometry; if small
    module-to-module differences are present, the first ready TX module is used
    as the representative single-module pattern.
    """
    modules = tx_module_positions() if modules is None else modules
    beam_vecs = tx_beam_unit_vectors() if beam_vecs is None else np.asarray(beam_vecs, dtype=np.float64)
    module_centers = tx_module_center_positions(modules)
    rel_module = modules[0] - np.mean(modules[0], axis=0, keepdims=True)
    dirs, original_shape = _as_direction_array(uvw)
    beam_arr = np.asarray(beam_id, dtype=np.int64)
    if beam_arr.ndim == 0:
        beam_arr = np.full(len(dirs), int(beam_arr), dtype=np.int64)
    else:
        beam_arr = np.broadcast_to(beam_arr, original_shape).reshape(-1)
    out = np.full(len(dirs), np.nan, dtype=np.float64)
    good = np.all(np.isfinite(dirs), axis=1) & (beam_arr >= 0) & (beam_arr < len(beam_vecs))
    if np.any(good):
        for beam in np.unique(beam_arr[good]):
            idx = good & (beam_arr == beam)
            steer = beam_vecs[beam]
            module_gain = module_power_gain(dirs[idx], steer=steer, module_positions=rel_module)
            sparse_gain = sparse_tx_array_power_gain(
                dirs[idx],
                int(beam),
                beam_vecs=beam_vecs,
                module_centers=module_centers,
            )
            out[idx] = module_gain * sparse_gain
    return out.reshape(original_shape)


def rx_power_gain(uvw: np.ndarray, channel: int | str | None = None, steer: np.ndarray | None = None) -> np.ndarray:
    """Receive power gain for one interferometric module, in linear units."""
    if channel is None:
        module_positions = None
    else:
        module_positions = rx_module_positions(channel)
    return module_power_gain(uvw, steer=steer, module_positions=module_positions)


def power_to_db(power: np.ndarray, floor: float = MIN_POWER_GAIN) -> np.ndarray:
    """Convert linear power gain to dB with a finite floor."""
    return 10.0 * np.log10(np.maximum(np.asarray(power, dtype=np.float64), floor))


def tx_gain_db(uvw: np.ndarray, beam_id: int | np.ndarray, beam_vecs: np.ndarray | None = None) -> np.ndarray:
    """Full transmit power gain in dB."""
    return power_to_db(tx_power_gain(uvw, beam_id, beam_vecs=beam_vecs))


def rx_gain_db(uvw: np.ndarray, channel: int | str | None = None, steer: np.ndarray | None = None) -> np.ndarray:
    """Receive module power gain in dB."""
    return power_to_db(rx_power_gain(uvw, channel=channel, steer=steer))


def two_way_power_gain(
    uvw: np.ndarray,
    beam_id: int | np.ndarray,
    beam_vecs: np.ndarray | None = None,
    rx_channel: int | str | None = None,
) -> np.ndarray:
    """Two-way meteor radar antenna gain in linear units.

    The two-way gain is the steered transmit power gain multiplied by the
    single-module receive power gain for the interferometric receiver.  Both
    transmit and receive module patterns are steered to the active beam.
    """
    dirs, original_shape = _as_direction_array(uvw)
    beam_vecs = tx_beam_unit_vectors() if beam_vecs is None else np.asarray(beam_vecs, dtype=np.float64)
    beam_arr = np.asarray(beam_id, dtype=np.int64)
    if beam_arr.ndim == 0:
        beam_arr = np.full(len(dirs), int(beam_arr), dtype=np.int64)
    else:
        beam_arr = np.broadcast_to(beam_arr, original_shape).reshape(-1)
    tx_gain = tx_power_gain(dirs, beam_arr, beam_vecs=beam_vecs).reshape(-1)
    rx_gain = np.full(len(dirs), np.nan, dtype=np.float64)
    module_positions = None if rx_channel is None else rx_module_positions(rx_channel)
    good = np.all(np.isfinite(dirs), axis=1) & (beam_arr >= 0) & (beam_arr < len(beam_vecs))
    if np.any(good):
        for beam in np.unique(beam_arr[good]):
            idx = good & (beam_arr == beam)
            rx_gain[idx] = module_power_gain(dirs[idx], steer=beam_vecs[beam], module_positions=module_positions)
    return (tx_gain * rx_gain).reshape(original_shape)


def two_way_gain_db(
    uvw: np.ndarray,
    beam_id: int | np.ndarray,
    beam_vecs: np.ndarray | None = None,
    rx_channel: int | str | None = None,
) -> np.ndarray:
    """Two-way meteor radar antenna gain in dB."""
    return power_to_db(two_way_power_gain(uvw, beam_id, beam_vecs=beam_vecs, rx_channel=rx_channel))


def precompute_tx_gain_maps(
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    valid: np.ndarray,
    beam_vecs: np.ndarray | None = None,
) -> np.ndarray:
    """Precompute full TX gain maps for all beams on a u/v/w sky grid."""
    beam_vecs = tx_beam_unit_vectors() if beam_vecs is None else np.asarray(beam_vecs, dtype=np.float64)
    uvw = np.column_stack([u[valid], v[valid], w[valid]])
    gain_maps = np.full((len(beam_vecs),) + u.shape, np.nan, dtype=np.float32)
    for beam_i in range(len(beam_vecs)):
        flat = np.full(u.shape, np.nan, dtype=np.float32)
        flat[valid] = tx_gain_db(uvw, beam_i, beam_vecs=beam_vecs).astype(np.float32)
        gain_maps[beam_i] = flat
    return gain_maps


def precompute_two_way_gain_maps(
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    valid: np.ndarray,
    beam_vecs: np.ndarray | None = None,
    rx_channel: int | str | None = None,
) -> np.ndarray:
    """Precompute two-way TX-RX gain maps for all beams on a u/v/w sky grid."""
    beam_vecs = tx_beam_unit_vectors() if beam_vecs is None else np.asarray(beam_vecs, dtype=np.float64)
    uvw = np.column_stack([u[valid], v[valid], w[valid]])
    gain_maps = np.full((len(beam_vecs),) + u.shape, np.nan, dtype=np.float32)
    for beam_i in range(len(beam_vecs)):
        flat = np.full(u.shape, np.nan, dtype=np.float32)
        flat[valid] = two_way_gain_db(uvw, beam_i, beam_vecs=beam_vecs, rx_channel=rx_channel).astype(np.float32)
        gain_maps[beam_i] = flat
    return gain_maps
