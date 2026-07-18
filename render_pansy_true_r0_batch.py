#!/usr/bin/env python3
from __future__ import annotations
import argparse, datetime, os, signal, traceback
from pathlib import Path
import numpy as np
import h5py
import digital_rf as drf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy.constants as c
import scipy.fftpack as fp
import scipy.optimize as opt

_default_repo = Path(__file__).resolve().parent
if not (_default_repo / "pansy_config.py").exists():
    _default_repo = Path("/home/j/src/pansy_receiver")
REPO = Path(os.environ.get("PANSY_RECEIVER_REPO", _default_repo))
os.chdir(REPO)
import sys
sys.path.insert(0, str(REPO))
import pansy_interferometry as pint
import pansy_modes as pmm
import pansy_config as pc
import fit_best_alias_physics_models as phys
import plot_interferometric_disambiguation as pid
from range_interpolation_study import PowerOnlyRangeDopplerSearch, DRG_KM, quadratic_peak_2d
from interferometer_alias_diagnostics import amp_scale

BASE = Path(os.environ.get("PANSY_BASE", "/mnt/data/juha/pansy"))
SRC_BATCH = BASE / "events/current_distribution_plots/example_event_batch_200"
OUT_BATCH = BASE / "events/current_distribution_plots/example_event_batch_200_true_r0_profile"
PLOTS = OUT_BATCH / "plots"
PROFILES = OUT_BATCH / "profiles"
HIGHRES = OUT_BATCH / "highres_fft_i2_p16"
PROFILE_INPUT_DIR = None
BEAT_H5_DIR = None
DEFAULT_SIGMA_POS = 0.5
DEFAULT_SIGMA_DOP = 1.0
RHO_METEOROID = 3000.0
DIAG_OVERRIDE_DIR = os.environ.get("PANSY_DIAG_OVERRIDE_DIR")

class Timeout(Exception):
    pass

def _timeout(signum, frame):
    raise Timeout()

class RDS:
    def __init__(self, txlen, echolen, nch):
        self.txlen = txlen
        self.fdec = 8
        self.fftlen = 512
        self.n_rg = echolen - txlen
        self.nch = nch
        self.fvec = np.fft.fftshift(np.fft.fftfreq(self.fftlen, d=self.fdec / 1e6))
        self.dopv = self.fvec * c.c / 2.0 / pc.freq
        idx = np.arange(self.txlen, dtype=np.int64)
        self.idx_mat = idx[None, :] + np.arange(self.n_rg, dtype=np.int64)[:, None]
    def decim(self, Z):
        nw = int(self.txlen / self.fdec)
        out = np.zeros((self.n_rg, nw), np.complex64)
        idx2 = np.arange(nw, dtype=np.int64) * self.fdec
        for ii in range(self.fdec):
            out += Z[:, idx2 + ii]
        return out
    def mf(self, tx, echo):
        MF = self.mf_matrix(tx, echo)
        noise = np.nanmedian(MF)
        return np.nanmax(MF, axis=1), noise
    def mf_matrix(self, tx, echo):
        MF = np.zeros((self.n_rg, self.fftlen), np.float32)
        for ch in range(self.nch):
            Z = echo[ch, self.idx_mat] * tx[None, :]
            ZF = np.fft.fftshift(fp.fft(self.decim(Z), self.fftlen, axis=1), axes=1)
            MF += ZF.real**2 + ZF.imag**2
        return MF

def highres_fft_measurements(sample_idx, cut, raw_idx, interp=2, fft_pad=16):
    HIGHRES.mkdir(parents=True, exist_ok=True)
    path = HIGHRES / f"highres_fft_i{interp}_p{fft_pad}_{sample_idx}.h5"
    if path.exists():
        with h5py.File(path, "r") as h:
            return {
                "range_km": h["range_km"][()],
                "doppler_km_s": h["doppler_mps"][()] / 1e3,
                "snr": h["snr"][()],
                "range_grid_m": float(h.attrs["range_grid_m"]),
                "doppler_step_mps": float(h.attrs["doppler_step_mps"]),
                "path": str(path),
            }
    z_rx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(cut["zrx_echoes_im"], dtype=np.complex64)
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.complex64)
    z_rx = z_rx * amp_scale()[None, :, None]
    delays = np.asarray(cut["delays"], dtype=np.int64)
    rds = PowerOnlyRangeDopplerSearch(
        txlen=z_tx.shape[1],
        echolen=z_rx.shape[2],
        n_channels=z_rx.shape[1],
        interp=int(interp),
        fft_pad=int(fft_pad),
    )
    range_km = np.full(len(raw_idx), np.nan, dtype=np.float64)
    doppler_mps = np.full(len(raw_idx), np.nan, dtype=np.float64)
    snr = np.full(len(raw_idx), np.nan, dtype=np.float64)
    for ii, pulse in enumerate(np.asarray(raw_idx, dtype=np.int64)):
        if pulse < 0 or pulse >= z_rx.shape[0]:
            continue
        mf = rds.mf(z_tx[pulse], z_rx[pulse])
        noise = float(np.nanmedian(mf))
        pprof = np.nanmax(mf, axis=1)
        ri = int(np.nanargmax(pprof))
        di = int(np.nanargmax(mf[ri]))
        log_mf = np.log(np.maximum(mf, max(noise, 1e-12)))
        dr, dd = quadratic_peak_2d(log_mf, ri, di)
        dop_step = float(rds.dopv[1] - rds.dopv[0])
        range_km[ii] = (delays[pulse] + (ri + dr) / float(interp)) * DRG_KM
        doppler_mps[ii] = rds.dopv[di] + dd * dop_step
        snr[ii] = (pprof[ri] - noise) / max(noise, 1e-12)
    with h5py.File(path, "w") as h:
        h["raw_idx"] = np.asarray(raw_idx, dtype=np.int64)
        h["range_km"] = range_km
        h["doppler_mps"] = doppler_mps
        h["snr"] = snr
        h.attrs["interp"] = int(interp)
        h.attrs["fft_pad"] = int(fft_pad)
        h.attrs["range_grid_m"] = float(DRG_KM * 1e3 / float(interp))
        h.attrs["doppler_step_mps"] = float(abs(rds.dopv[1] - rds.dopv[0]))
        h.attrs["source"] = "render_pansy_true_r0_batch.py highres_fft_measurements"
    return {
        "range_km": range_km,
        "doppler_km_s": doppler_mps / 1e3,
        "snr": snr,
        "range_grid_m": float(DRG_KM * 1e3 / float(interp)),
        "doppler_step_mps": float(abs(rds.dopv[1] - rds.dopv[0])),
        "path": str(path),
    }

def beam_pixmap():
    mm = pmm.get_m_mode()
    u, v, w = pint.uv_coverage(N=450, max_zenith_angle=20.0)
    bitmap = np.zeros(u.shape)
    for ii in range(5):
        az = mm["beam_pos_az_za"][ii][0]
        el = 90.0 - mm["beam_pos_az_za"][ii][1]
        ph = np.cos(np.pi * el / 180.0)
        pw = -np.sin(np.pi * el / 180.0)
        pv = ph * np.cos(-np.pi * az / 180.0)
        pu = -ph * np.sin(-np.pi * az / 180.0)
        bitmap[np.rad2deg(np.arccos(np.clip(u * pu + v * pv + pw * w, -1, 1))) < 10.0] += 1
    return u * 100, v * 100, bitmap

EWM, NSM, BMAP = beam_pixmap()

def r_um_to_mass_kg(r_um):
    r_um = np.asarray(r_um, dtype=float)
    return (4.0 / 3.0) * np.pi * RHO_METEOROID * (np.maximum(r_um, 1e-300) * 1e-6) ** 3

def mass_kg_to_r_um(m_kg):
    m_kg = np.asarray(m_kg, dtype=float)
    return ((3.0 * np.maximum(m_kg, 1e-300)) / (4.0 * np.pi * RHO_METEOROID)) ** (1.0 / 3.0) * 1e6

def log_grid_weights(radius_um, density_log_radius):
    radius_um = np.asarray(radius_um, dtype=float)
    density_log_radius = np.asarray(density_log_radius, dtype=float)
    x = np.log(radius_um)
    if len(x) == 1:
        w = np.array([1.0])
    else:
        edges = np.empty(len(x) + 1, dtype=float)
        edges[1:-1] = 0.5 * (x[:-1] + x[1:])
        edges[0] = x[0] - 0.5 * (x[1] - x[0])
        edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
        w = density_log_radius * np.diff(edges)
    w[~np.isfinite(w)] = 0.0
    w = np.maximum(w, 0.0)
    sw = np.sum(w)
    return w / sw if sw > 0 else w

def weighted_quantile(values, weights, qs):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    qs = np.asarray(qs, dtype=float)
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if np.count_nonzero(good) == 0:
        return np.full(qs.shape, np.nan, dtype=float)
    v = values[good]
    w = weights[good]
    order = np.argsort(v)
    v = v[order]
    w = w[order]
    cw = np.cumsum(w)
    cw /= cw[-1]
    return np.interp(qs, cw, v, left=v[0], right=v[-1])

def weighted_pointwise_interval(samples, weights, qs=(0.025, 0.975)):
    samples = np.asarray(samples, dtype=float)
    out = np.full((len(qs), samples.shape[1]), np.nan, dtype=float)
    for ii in range(samples.shape[1]):
        out[:, ii] = weighted_quantile(samples[:, ii], weights, qs)
    return out

def interval_text(lo, hi, unit="", sci=False):
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return "n/a"
    if sci:
        return f"{lo:.2g}-{hi:.2g}{unit}"
    if hi < 10:
        return f"{lo:.2f}-{hi:.2f}{unit}"
    return f"{lo:.0f}-{hi:.0f}{unit}"

def fmt_metric(value, precision=2):
    return f"{value:.{precision}f}" if np.isfinite(value) else "n/a"

def group_attr_float(group, name):
    try:
        return float(group.attrs[name])
    except Exception:
        return np.nan

def load_rows(limit=None, sample_index_h5=None):
    if sample_index_h5 is not None:
        with h5py.File(sample_index_h5, "r") as h:
            samples = np.asarray(h["sample_idx"], dtype=np.int64)
        if limit is not None:
            samples = samples[:limit]
        return [
            (
                rank,
                int(sample),
                datetime.datetime.fromtimestamp(
                    float(sample) / 1e6, tz=datetime.timezone.utc
                ).strftime("%Y-%m-%d"),
            )
            for rank, sample in enumerate(samples, start=1)
        ]

    rows = []
    with (SRC_BATCH / "candidate_manifest.tsv").open() as f:
        next(f)
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 4:
                continue
            rows.append((int(p[0]), int(p[2]), p[3]))
            if limit and len(rows) >= limit:
                break
    return rows

def diagnostic_path(sample_idx, day):
    if DIAG_OVERRIDE_DIR:
        override = Path(DIAG_OVERRIDE_DIR) / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"
        if override.exists():
            return override
    return BASE / "events" / day / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"

def coherent_doppler_path(sample_idx):
    return BASE / "events/current_distribution_plots/coherent_doppler_tests_no_ipp" / f"coherent_doppler_{sample_idx}.h5"

def load_coherent_doppler_overlay(sample_idx):
    path = coherent_doppler_path(sample_idx)
    if not path.exists():
        return None
    with h5py.File(path, "r") as h:
        g = h["estimates"]
        return {
            "t_rel_s": np.asarray(g["t_rel_s"], dtype=np.float64),
            "doppler_km_s": np.asarray(g["doppler_coherent_steered_mps"], dtype=np.float64) / 1e3,
            "fft_pad": int(h.attrs.get("fft_pad", 1)),
            "n_ipp": int(h.attrs.get("integrated_pulses_nominal", 1)),
        }

def lag_doppler_5pulse(cut, obs_abs_s, observed_range, coarse_doppler_km_s, lag_samples=64, half_width=2):
    zrx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(cut["zrx_echoes_im"], dtype=np.complex64)
    ztx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.complex64)
    scale = pid.amp_scale()
    for i in range(zrx.shape[1]):
        zrx[:, i, :] *= scale[i]
    raw_abs_s = np.asarray(cut["tx_idx"], dtype=np.float64) / 1e6
    delays = np.asarray(cut["delays"], dtype=np.int64)
    drg_km = c.c / 2.0 / 1e6 / 1e3
    raw_idx = np.asarray([int(np.nanargmin(np.abs(raw_abs_s - tt))) for tt in obs_abs_s], dtype=np.int64)
    range_bin = np.rint(np.asarray(observed_range, dtype=np.float64) / drg_km - delays[raw_idx]).astype(np.int64)
    corr = np.full(len(obs_abs_s), np.nan + 1j * np.nan, dtype=np.complex128)
    for ii, (pulse, rg) in enumerate(zip(raw_idx, range_bin)):
        if not np.isfinite(observed_range[ii]) or rg < 0:
            continue
        stop = int(rg) + ztx.shape[1]
        if stop > zrx.shape[2] or lag_samples <= 0 or lag_samples >= ztx.shape[1]:
            continue
        acc = 0.0j
        for ch in range(zrx.shape[1]):
            series = zrx[pulse, ch, int(rg):stop] * ztx[pulse]
            acc += np.sum(series[lag_samples:] * np.conj(series[:-lag_samples]))
        corr[ii] = acc
    lag_s = float(lag_samples) / 1e6
    period_mps = pc.wavelength / (2.0 * lag_s)
    out = np.full(len(obs_abs_s), np.nan, dtype=np.float64)
    for ii in range(len(corr)):
        lo = max(0, ii - int(half_width))
        hi = min(len(corr), ii + int(half_width) + 1)
        window = corr[lo:hi]
        window = window[np.isfinite(window.real) & np.isfinite(window.imag)]
        if len(window) == 0:
            continue
        csum = np.sum(window)
        if not np.isfinite(csum.real) or not np.isfinite(csum.imag) or abs(csum) == 0.0:
            continue
        wrapped = np.angle(csum) * pc.wavelength / (4.0 * np.pi * lag_s)
        coarse = np.nanmedian(np.asarray(coarse_doppler_km_s[lo:hi], dtype=np.float64)) * 1e3
        if np.isfinite(coarse):
            wrapped += np.rint((coarse - wrapped) / period_mps) * period_mps
        out[ii] = wrapped / 1e3
    return out

def true_radius_profile(sample_idx, day, grid_n=41, timeout_s=5):
    if PROFILE_INPUT_DIR is not None:
        profile_path = PROFILE_INPUT_DIR / f"mass_profile_{sample_idx}.h5"
        if not profile_path.exists():
            profile_path = PROFILE_INPUT_DIR / f"mass_profile_with_acceleration_{sample_idx}.h5"
        if not profile_path.exists():
            raise FileNotFoundError(f"missing completed mass profile: {profile_path}")
        with h5py.File(profile_path, "r") as h:
            radius_um = np.asarray(h["profile/radius_um"], dtype=float)
            params6 = np.asarray(h["profile/parameters6"], dtype=float)
            success = np.asarray(h["profile/success"], dtype=bool)
            probability_key = (
                "profile/relative_probability_density_log_radius"
                if "profile/relative_probability_density_log_radius" in h
                else "profile/probability_density_log_radius"
            )
            best_radius_key = (
                "result/free_best_radius_um"
                if "result/free_best_radius_um" in h
                else "result/best_radius_um"
            )
            return {
                "radius_um": radius_um,
                "prob": np.asarray(h[probability_key], dtype=float),
                "params6": params6,
                "success": success,
                "free_radius_um": float(h[best_radius_key][()]),
                "path": str(profile_path),
            }

    profile_path = PROFILES / f"radius_profile_{sample_idx}_true_1um_10000um_fit_rms_sigma.h5"
    if profile_path.exists():
        with h5py.File(profile_path, "r") as h:
            params6 = h["params6"][()] if "params6" in h else np.full((len(h["radius_um"]), 6), np.nan)
            success = h["success"][()] if "success" in h else np.isfinite(params6).all(axis=1)
            return {
                "radius_um": h["radius_um"][()],
                "prob": h["relative_probability_density_log_radius"][()],
                "params6": params6,
                "success": success,
                "free_radius_um": float(h.attrs.get("free_radius_um", np.nan)),
                "path": str(profile_path),
            }
    diag = diagnostic_path(sample_idx, day)
    obs = phys.load_best_alias(diag)
    rho, _ = phys.pbal.density_interpolator(obs["sample_epoch_unix"])
    t = obs["t_s"]
    points = obs["points_km"]
    dop = obs["doppler_km_s"]
    finite = np.isfinite(t) & np.all(np.isfinite(points), axis=1) & np.isfinite(dop)
    with h5py.File(diag, "r") as h:
        sel = h.attrs["selected_hypothesis"]
        sel = sel.decode() if isinstance(sel, bytes) else sel
        g = h["hypotheses"][sel]
        best_params = np.asarray(g["physics_ceplecha_params"][()], float)
        best_model = np.asarray(g["physics_ceplecha_model"][()], float)
        best_vel = np.asarray(g["physics_ceplecha_velocity_km_s"][()], float)
        best_pred = phys.predicted_doppler(best_model, best_vel)
        keep = np.asarray(g["physics_ceplecha_keep"][()], bool) if "physics_ceplecha_keep" in g else finite.copy()
    if np.count_nonzero(keep) < 20:
        keep = finite.copy()
    best_pos_res = points[keep] - best_model[keep]
    best_dop_res = dop[keep] - best_pred[keep]
    sigma_pos = float(np.sqrt(np.nanmean(best_pos_res**2)))
    sigma_dop = float(np.sqrt(np.nanmean(best_dop_res**2)))
    if not np.isfinite(sigma_pos) or sigma_pos <= 0:
        sigma_pos = DEFAULT_SIGMA_POS
    if not np.isfinite(sigma_dop) or sigma_dop <= 0:
        sigma_dop = DEFAULT_SIGMA_DOP
    best_raw = np.r_[((points[keep] - best_model[keep]) / sigma_pos).ravel(), (dop[keep] - best_pred[keep]) / sigma_dop]
    free_chi2 = float(np.sum(best_raw**2))
    free_radius_um = float(10 ** best_params[6] * 1e6)
    radius_um = np.geomspace(1.0, 10000.0, int(grid_n))
    chi2 = np.full(radius_um.shape, np.nan, float)
    ok = np.zeros(radius_um.shape, bool)
    params6 = np.full((len(radius_um), 6), np.nan, float)
    bounds = (np.array([-np.inf, -np.inf, 20e3, -90e3, -90e3, -90e3]), np.array([np.inf, np.inf, 220e3, 90e3, 90e3, 90e3]))
    x_scale = np.array([1e5, 1e5, 1e5, 7e4, 7e4, 7e4])
    def residual6(p6, lr):
        p = np.r_[p6, lr]
        pos, vel, rad, mass, success, msg = phys.propagate_shrinking_radius_model(p, t, rho)
        pred = phys.predicted_doppler(pos, vel)
        return np.r_[((points[keep] - pos[keep]) / sigma_pos).ravel(), (dop[keep] - pred[keep]) / sigma_dop]
    old_handler = signal.signal(signal.SIGALRM, _timeout)
    order = np.argsort(np.abs(np.log10(radius_um * 1e-6) - best_params[6]))
    last = best_params[:6].copy()
    try:
        for idx in order:
            lr = float(np.log10(radius_um[idx] * 1e-6))
            best = None
            try:
                signal.alarm(int(timeout_s))
                for st in (last, best_params[:6]):
                    res = opt.least_squares(lambda x: residual6(x, lr), st, bounds=bounds, x_scale=x_scale, loss="linear", max_nfev=35)
                    if best is None or res.cost < best.cost:
                        best = res
                signal.alarm(0)
                rr = residual6(best.x, lr)
                chi2[idx] = float(np.sum(rr**2))
                params6[idx] = best.x
                ok[idx] = bool(best.success)
                last = best.x.copy()
            except Exception:
                signal.alarm(0)
                ok[idx] = False
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
    finite_prof = np.isfinite(chi2)
    min_chi2 = float(np.nanmin(np.r_[chi2[finite_prof], free_chi2])) if np.any(finite_prof) else free_chi2
    dchi2 = chi2 - min_chi2
    prob = np.exp(-0.5 * np.clip(dchi2, 0, 700))
    prob[~np.isfinite(prob)] = 0.0
    integ = np.trapz(prob, np.log(radius_um))
    if integ > 0:
        prob /= integ
    PROFILES.mkdir(parents=True, exist_ok=True)
    with h5py.File(profile_path, "w") as h:
        h["radius_um"] = radius_um
        h["chi2"] = chi2
        h["delta_chi2"] = dchi2
        h["relative_probability_density_log_radius"] = prob
        h["success"] = ok
        h["params6"] = params6
        h.attrs["free_radius_um"] = free_radius_um
        h.attrs["free_chi2"] = free_chi2
        h.attrs["min_chi2"] = min_chi2
        h.attrs["n_keep"] = int(np.count_nonzero(keep))
        h.attrs["sigma_pos_km"] = float(sigma_pos)
        h.attrs["sigma_dop_km_s"] = float(sigma_dop)
        h.attrs["sigma_model"] = "fit residual RMS; one pooled position-component sigma and one Doppler sigma, fixed across r0 profile"
        h.attrs["radius_max_um"] = 10000.0
        h.attrs["note"] = "True fixed-r0 profile likelihood. Timed-out points have zero probability."
    return {
        "radius_um": radius_um,
        "prob": prob,
        "params6": params6,
        "success": ok,
        "free_radius_um": free_radius_um,
        "path": str(profile_path),
    }

def render(row, grid_n=41):
    rank, sample, day = row
    out = PLOTS / f"rank_{rank:03d}_{sample}.png"
    profile = true_radius_profile(sample, day, grid_n=grid_n)
    profile_radius_um = profile["radius_um"]
    profile_prob = profile["prob"]
    profile_params6 = profile["params6"]
    profile_success = profile["success"]
    free_radius_um = profile["free_radius_um"]
    profile_path = profile["path"]
    diag = diagnostic_path(sample, day)
    sample_epoch_unix = None
    observed_range = None
    stored_redchi = stored_pos_rms = stored_dop_rms = np.nan
    fit_n = np.nan
    fit_keep = None
    selected_range_time_component = -1
    with h5py.File(diag, "r") as h:
        sel = h.attrs["selected_hypothesis"]
        sel = sel.decode() if isinstance(sel, bytes) else sel
        sample_epoch_unix = float(h.attrs["sample_epoch_unix"])
        selected_range_time_component = int(h.attrs.get("selected_range_time_component", -1))
        ts = datetime.datetime.utcfromtimestamp(sample_epoch_unix).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"
        g = h["hypotheses"][sel]
        t = g["t_rel_s"][()]
        pos = g["position_enu_km"][()]
        model = g["physics_ceplecha_model"][()] if "physics_ceplecha_model" in g else g["selection_model"][()]
        snr = g["snr"][()]
        dop = g["doppler_mps"][()] / 1e3
        if "physics_ceplecha_velocity_km_s" in g:
            cepl_vel = g["physics_ceplecha_velocity_km_s"][()]
            pred = phys.predicted_doppler(model, cepl_vel)
        else:
            cepl_vel = None
            pred = g["selection_pred_doppler_km_s"][()] if "selection_pred_doppler_km_s" in g else np.full_like(t, np.nan)
        observed_range = g["range_km"][()] if "range_km" in g else None
        fit_keep = g["physics_ceplecha_keep"][()] if "physics_ceplecha_keep" in g else np.ones(len(t), dtype=bool)
        stored_redchi = group_attr_float(g, "physics_ceplecha_reduced_chi2")
        stored_pos_rms = group_attr_float(g, "physics_ceplecha_pos_rms_km")
        stored_dop_rms = group_attr_float(g, "physics_ceplecha_dop_rms_km_s")
        fit_n = group_attr_float(g, "physics_ceplecha_n")
    o = np.argsort(t)
    t = t[o]; pos = pos[o]; model = model[o]; snr = snr[o]; dop = dop[o]; pred = pred[o]
    if observed_range is not None:
        observed_range = observed_range[o]
    if fit_keep is not None:
        fit_keep = np.asarray(fit_keep, dtype=bool)[o]
    else:
        fit_keep = np.ones(len(t), dtype=bool)
    if not np.any(fit_keep):
        fit_keep = np.ones(len(t), dtype=bool)
    if np.any(fit_keep):
        t_fit_origin = float(t[fit_keep][0])
    else:
        t_fit_origin = float(t[0])
    t_plot = t - t_fit_origin
    line_mask = fit_keep & np.isfinite(t_plot)
    if cepl_vel is not None:
        cepl_vel = cepl_vel[o]
    snr_db = 10 * np.log10(np.maximum(snr, 1e-12))
    range_model = np.linalg.norm(model, axis=1)
    pos_residual = pos - model
    residual_mask = fit_keep & np.all(np.isfinite(pos_residual), axis=1)
    if not np.any(residual_mask):
        residual_mask = np.all(np.isfinite(pos_residual), axis=1)
    pos_rms = float(np.sqrt(np.nanmean(np.sum(pos_residual[residual_mask] ** 2, axis=1))))
    ew_rms = float(np.sqrt(np.nanmean(pos_residual[residual_mask, 0] ** 2)))
    ns_rms = float(np.sqrt(np.nanmean(pos_residual[residual_mask, 1] ** 2)))
    up_rms = float(np.sqrt(np.nanmean(pos_residual[residual_mask, 2] ** 2)))
    dop_rms = float(np.sqrt(np.nanmean((dop - pred) ** 2)))
    if observed_range is not None:
        range_rms = float(np.sqrt(np.nanmean((observed_range - range_model) ** 2)))
    else:
        range_rms = np.nan
    if not np.isfinite(stored_pos_rms):
        stored_pos_rms = pos_rms
    if not np.isfinite(stored_dop_rms):
        stored_dop_rms = dop_rms
    path_length_km = float(np.nansum(np.linalg.norm(np.diff(model, axis=0), axis=1)))
    profile_weights = log_grid_weights(profile_radius_um, profile_prob)
    valid_profile = np.isfinite(profile_weights) & (profile_weights > 0) & np.asarray(profile_success, bool) & np.all(np.isfinite(profile_params6), axis=1)
    up_interval = dop_interval = speed_interval = range_interval = None
    fixed_r0_dopplers = {}
    fixed_r0_speeds = {}
    best_profile_doppler = pred
    if np.count_nonzero(valid_profile) >= 2:
        rho, _ = phys.pbal.density_interpolator(sample_epoch_unix)
        pos_samples = []
        dop_samples = []
        speed_samples = []
        range_samples = []
        weights_used = []
        for radius_um, p6, w in zip(profile_radius_um[valid_profile], profile_params6[valid_profile], profile_weights[valid_profile]):
            p = np.r_[p6, np.log10(radius_um * 1e-6)]
            try:
                ppos, pvel, _rad, _mass, _success, _msg = phys.propagate_shrinking_radius_model(p, t, rho)
                ppred = phys.predicted_doppler(ppos, pvel)
            except Exception:
                continue
            if not (np.all(np.isfinite(ppos)) and np.all(np.isfinite(pvel)) and np.all(np.isfinite(ppred))):
                continue
            pos_samples.append(ppos)
            dop_samples.append(ppred)
            speed_samples.append(np.linalg.norm(pvel, axis=1))
            range_samples.append(np.linalg.norm(ppos, axis=1))
            weights_used.append(w)
        if len(weights_used) >= 2:
            weights_used = np.asarray(weights_used, dtype=float)
            weights_used /= np.sum(weights_used)
            pos_samples = np.asarray(pos_samples)
            dop_samples = np.asarray(dop_samples)
            speed_samples = np.asarray(speed_samples)
            range_samples = np.asarray(range_samples)
            up_interval = weighted_pointwise_interval(pos_samples[:, :, 2], weights_used)
            dop_interval = weighted_pointwise_interval(dop_samples, weights_used)
            speed_interval = weighted_pointwise_interval(speed_samples, weights_used)
            range_interval = weighted_pointwise_interval(range_samples, weights_used)
        for target_um in (10.0, 100.0, 1000.0):
            valid_idx = np.flatnonzero(valid_profile)
            if len(valid_idx) == 0:
                continue
            idx = valid_idx[np.argmin(np.abs(np.log(profile_radius_um[valid_idx]) - np.log(target_um)))]
            try:
                p = np.r_[profile_params6[idx], np.log10(profile_radius_um[idx] * 1e-6)]
                ppos, pvel, _rad, _mass, _success, _msg = phys.propagate_shrinking_radius_model(p, t, rho)
                ppred = phys.predicted_doppler(ppos, pvel)
            except Exception:
                continue
            pspeed = np.linalg.norm(pvel, axis=1)
            if np.all(np.isfinite(ppred)) and np.all(np.isfinite(pspeed)):
                fixed_r0_dopplers[target_um] = ppred
                fixed_r0_speeds[target_um] = pspeed
        valid_idx = np.flatnonzero(valid_profile)
        if len(valid_idx):
            idx = valid_idx[np.argmin(np.abs(np.log(profile_radius_um[valid_idx]) - np.log(free_radius_um)))]
            try:
                p = np.r_[profile_params6[idx], np.log10(profile_radius_um[idx] * 1e-6)]
                ppos, pvel, _rad, _mass, _success, _msg = phys.propagate_shrinking_radius_model(p, t, rho)
                candidate_doppler = phys.predicted_doppler(ppos, pvel)
                if np.all(np.isfinite(candidate_doppler)):
                    best_profile_doppler = candidate_doppler
            except Exception:
                pass
    r0_q = weighted_quantile(profile_radius_um, profile_weights, [0.025, 0.975])
    m0_q = r_um_to_mass_kg(r0_q)
    cut = drf.DigitalMetadataReader(str(BASE / "metadata/cut")).read(sample - 1, sample + 1)
    if sample not in cut:
        raise RuntimeError("missing cut")
    cut = cut[sample]
    obs_rti = pid.recompute_cut_observables(cut, interp=1)
    snr_gate = np.asarray(obs_rti["snr"], dtype=np.float64) > 7.0
    obs_for_clock = {
        key: val[snr_gate] if isinstance(val, np.ndarray) and len(val) == len(snr_gate) else val
        for key, val in obs_rti.items()
    }
    if selected_range_time_component >= 0:
        segments = pid.split_observations_by_range_time(obs_for_clock, min_points=3)
        if selected_range_time_component < len(segments):
            obs_for_clock = pid.subset_pulse_observations(obs_for_clock, segments[selected_range_time_component])
    if len(obs_for_clock.get("tx_idx", [])) == 0:
        raise RuntimeError("cannot reconstruct diagnostic measurement clock from cut observables")
    diagnostic_zero_abs_s = float(np.asarray(obs_for_clock["tx_idx"], dtype=np.float64)[0]) / 1e6
    obs_abs_s = diagnostic_zero_abs_s + np.asarray(t, dtype=np.float64)
    fit_origin_abs_s = diagnostic_zero_abs_s + float(t_fit_origin)
    t_plot = obs_abs_s - fit_origin_abs_s
    line_mask = fit_keep & np.isfinite(t_plot)
    raw_abs_s = np.asarray(cut["tx_idx"], dtype=np.float64) / 1e6
    raw_t = raw_abs_s - fit_origin_abs_s
    good_range = np.isfinite(obs_abs_s) & np.isfinite(range_model)
    if np.count_nonzero(good_range) < 2:
        raise RuntimeError("not enough finite fit samples for raw-pulse interpolation")
    raw_range_model = np.interp(
        raw_abs_s,
        obs_abs_s[good_range],
        range_model[good_range],
        left=range_model[good_range][0],
        right=range_model[good_range][-1],
    )
    rti_abs_s = np.asarray(obs_rti["tx_idx"], dtype=np.float64) / 1e6
    RTI_raw_db = 10 * np.log10(np.maximum(np.asarray(obs_rti["rti_snr"], dtype=np.float64), 1e-12))
    raw_for_measurement = np.asarray([int(np.nanargmin(np.abs(rti_abs_s - tt))) for tt in obs_abs_s], dtype=np.int64)
    RTI = RTI_raw_db[raw_for_measurement, :]
    rti_t = t_plot
    rvec = np.asarray(obs_rti["range_grid_km"], dtype=np.float64)
    highres_raw_idx = np.asarray([int(np.nanargmin(np.abs(raw_abs_s - tt))) for tt in obs_abs_s], dtype=np.int64)
    highres = highres_fft_measurements(sample, cut, highres_raw_idx, interp=2, fft_pad=16)
    plot_range = np.asarray(highres["range_km"], dtype=np.float64)
    plot_dop = np.asarray(highres["doppler_km_s"], dtype=np.float64)
    plot_range_rms = float(np.sqrt(np.nanmean((plot_range - range_model) ** 2)))
    plot_dop_rms = float(np.sqrt(np.nanmean((plot_dop - pred) ** 2)))
    both = np.r_[snr_db[np.isfinite(snr_db)], RTI[np.isfinite(RTI)]]
    snr_vmin = 0.0
    snr_vmax = max(18.0, float(np.nanpercentile(both, 99.7)))
    prob_rel = profile_prob / np.nanmax(profile_prob) if np.nanmax(profile_prob) > 0 else profile_prob
    fig, axs = plt.subplots(2, 4, figsize=(16, 8), dpi=120)
    ax = axs[0, 0]
    ax.pcolormesh(EWM, NSM, BMAP, cmap="gist_yarg", vmax=5, shading="auto")
    dropped_mask = (~fit_keep) & np.isfinite(t_plot)
    kept_mask = fit_keep & np.isfinite(t_plot)
    if np.any(dropped_mask):
        ax.scatter(pos[dropped_mask, 0], pos[dropped_mask, 1], color="0.75", s=5, edgecolors="none", zorder=2)
    sc = ax.scatter(pos[kept_mask, 0], pos[kept_mask, 1], c=snr_db[kept_mask], cmap="plasma", vmin=snr_vmin, vmax=snr_vmax, s=5, edgecolors="none", zorder=3)
    ax.plot(model[line_mask, 0], model[line_mask, 1], color="tab:blue", lw=1.0, zorder=4)
    ax.text(pos[0, 0], pos[0, 1], r"$t_0$", fontsize=7)
    ax.text(
        0.03,
        0.97,
        f"Path {path_length_km:.1f} km\nEW RMS {fmt_metric(ew_rms)} km\nNS RMS {fmt_metric(ns_rms)} km",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    cax = inset_axes(ax, width="4%", height="35%", loc="lower right", borderpad=2.0)
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("SNR (dB)", fontsize=7)
    cb.ax.yaxis.set_label_position("left")
    cb.ax.yaxis.set_ticks_position("right")
    cb.ax.tick_params(labelsize=7, length=2)
    ax.set_xlabel("EW (km)"); ax.set_ylabel("NS (km)"); ax.set_xlim(-35, 35); ax.set_ylim(-35, 35); ax.set_aspect("equal", adjustable="box"); ax.grid(alpha=0.15, lw=0.4)
    ax = axs[0, 1]
    if up_interval is not None:
        ax.fill_between(t_plot[line_mask], up_interval[0, line_mask], up_interval[1, line_mask], color="tab:blue", alpha=0.18, lw=0)
    if np.any(dropped_mask):
        ax.scatter(t_plot[dropped_mask], pos[dropped_mask, 2], color="0.75", s=5, edgecolors="none")
    ax.scatter(t_plot[kept_mask], pos[kept_mask, 2], color="black", s=5, edgecolors="none"); ax.plot(t_plot[line_mask], model[line_mask, 2], color="tab:blue", lw=1.0)
    ax.text(0.03, 0.97, f"Up RMS {fmt_metric(up_rms)} km", transform=ax.transAxes, ha="left", va="top", fontsize=8)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Up (km)"); ax.set_title(f"#{rank:03d} {ts}", fontsize=10); ax.grid(alpha=0.2, lw=0.4)
    ax = axs[0, 2]
    if dop_interval is not None:
        ax.fill_between(t_plot[line_mask], dop_interval[0, line_mask], dop_interval[1, line_mask], color="tab:blue", alpha=0.18, lw=0)
    valid = np.isfinite(plot_dop)
    dropped = valid & dropped_mask
    kept = valid & kept_mask
    if np.any(dropped):
        ax.scatter(t_plot[dropped], plot_dop[dropped], color="0.75", s=5, edgecolors="none")
    ax.scatter(t_plot[kept], plot_dop[kept], color="black", s=5, edgecolors="none"); ax.plot(t_plot[line_mask], pred[line_mask], color="tab:blue", lw=1.0, label="free fit")
    fixed_colors = {10.0: "tab:purple", 100.0: "tab:orange", 1000.0: "tab:red"}
    for target_um, fixed_pred in fixed_r0_dopplers.items():
        fixed_mask = line_mask & np.isfinite(fixed_pred)
        if np.any(fixed_mask):
            ax.plot(
                t_plot[fixed_mask],
                fixed_pred[fixed_mask],
                color=fixed_colors.get(target_um, "0.35"),
                lw=0.9,
                ls="--",
                label=rf"$r_0={target_um:g}\,\mu$m",
            )
    ax.text(0.03, 0.97, f"Dop RMS {fmt_metric(plot_dop_rms)} km/s\nI2 P16", transform=ax.transAxes, ha="left", va="top", fontsize=8)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Doppler (km/s)"); ax.grid(alpha=0.2, lw=0.4)
    if fixed_r0_dopplers:
        ax.legend(loc="lower right", fontsize=7, frameon=False)
    ax = axs[0, 3]
    if BEAT_H5_DIR is not None:
        from fit_inter_pulse_acceleration_mass import load_nonoverlapping_phase_acceleration

        beat_path = BEAT_H5_DIR / f"inter_pulse_phase_{sample}.h5"
        phase_data = load_nonoverlapping_phase_acceleration(beat_path)
        phase_samples = phase_data["samples"]
        phase_dt = 0.5 * (phase_samples["first_dt_s"] + phase_samples["second_dt_s"])
        phase_acceleration = (
            pc.wavelength
            * phase_samples["observed_delta_phase_rad"]
            / (4.0 * np.pi * phase_dt**2)
        )
        with h5py.File(profile_path, "r") as profile_handle:
            acceleration_key = "phase_acceleration/measured_radial_acceleration_mps2"
            if acceleration_key in profile_handle:
                fitted_acceleration = np.asarray(profile_handle[acceleration_key], dtype=float)
                if fitted_acceleration.shape == phase_acceleration.shape:
                    phase_acceleration = fitted_acceleration
        phase_acceleration_std = (
            pc.wavelength
            * phase_samples["formal_phase_std_rad"]
            / (4.0 * np.pi * phase_dt**2)
        )
        phase_time_abs = 0.5 * (phase_samples["first_time_s"] + phase_samples["second_time_s"])
        phase_time = phase_time_abs - fit_origin_abs_s
        first_model_velocity = np.interp(phase_samples["first_time_s"], obs_abs_s, best_profile_doppler * 1e3)
        second_model_velocity = np.interp(phase_samples["second_time_s"], obs_abs_s, best_profile_doppler * 1e3)
        model_acceleration = (
            (second_model_velocity - first_model_velocity)
            / (phase_samples["second_time_s"] - phase_samples["first_time_s"])
        )
        fit_good = (
            np.isfinite(phase_time)
            & np.isfinite(phase_acceleration)
            & np.isfinite(phase_acceleration_std)
        )
        fit_center = float(np.nanmedian(phase_time[fit_good]))
        fit_x = phase_time - fit_center
        fit_keep_phase = fit_good.copy()
        phase_line = np.asarray([0.0, np.nanmedian(phase_acceleration[fit_good])])
        for _ in range(5):
            phase_line = np.polyfit(
                fit_x[fit_keep_phase],
                phase_acceleration[fit_keep_phase],
                1,
                w=1.0 / np.maximum(phase_acceleration_std[fit_keep_phase], 1.0),
            )
            phase_residual = phase_acceleration - np.polyval(phase_line, fit_x)
            phase_median = np.nanmedian(phase_residual[fit_keep_phase])
            phase_sigma = 1.4826 * np.nanmedian(
                np.abs(phase_residual[fit_keep_phase] - phase_median)
            )
            if not np.isfinite(phase_sigma) or phase_sigma <= 0.0:
                break
            fit_keep_phase = fit_good & (np.abs(phase_residual - phase_median) < 4.0 * phase_sigma)
        phase_order = np.argsort(phase_time)
        ax.errorbar(
            phase_time,
            phase_acceleration / 1e3,
            yerr=phase_acceleration_std / 1e3,
            fmt=".",
            ms=3,
            color="black",
            ecolor="0.78",
            elinewidth=0.5,
            label="decoded beat phase",
        )
        ax.plot(
            phase_time[phase_order],
            model_acceleration[phase_order] / 1e3,
            color="tab:blue",
            lw=1.0,
            label="shrinking-radius model",
        )
        ax.plot(
            phase_time[phase_order],
            np.polyval(phase_line, fit_x[phase_order]) / 1e3,
            color="tab:red",
            lw=1.1,
            label="weighted measurement fit",
        )
        ambiguity = pc.wavelength / (4.0 * np.nanmedian(phase_dt) ** 2) / 1e3
        ax.axhline(ambiguity, color="0.5", ls="--", lw=0.7)
        ax.axhline(-ambiguity, color="0.5", ls="--", lw=0.7)
        ax.text(
            0.03,
            0.97,
            f"N {np.count_nonzero(fit_good)}\n"
            + f"fit {np.polyval(phase_line, 0.0) / 1e3:.2f} km s$^{{-2}}$\n"
            + f"model {np.nanmedian(model_acceleration) / 1e3:.2f} km s$^{{-2}}$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
        )
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"Radial acceleration (km s$^{-2}$)")
        ax.grid(alpha=0.2, lw=0.4)
        ax.legend(frameon=False, fontsize=6.5, loc="lower left")
    else:
        ax.axis("off")
    ax = axs[1, 0]
    ax.pcolormesh(rti_t, rvec, RTI.T, cmap="plasma", shading="auto", vmin=snr_vmin, vmax=snr_vmax)
    if range_interval is not None:
        ax.fill_between(t_plot[line_mask], range_interval[0, line_mask], range_interval[1, line_mask], color="tab:blue", alpha=0.18, lw=0)
    if plot_range is not None:
        range_valid = np.isfinite(plot_range) & np.isfinite(t_plot)
        range_dropped = range_valid & dropped_mask
        range_kept = range_valid & kept_mask
        if np.any(range_dropped):
            ax.scatter(t_plot[range_dropped], plot_range[range_dropped], color="0.65", s=5, edgecolors="none", alpha=0.75, zorder=5)
        if np.any(range_kept):
            ax.scatter(t_plot[range_kept], plot_range[range_kept], color="black", s=6, edgecolors="white", linewidths=0.15, alpha=0.9, zorder=6)
    ax.plot(t_plot[line_mask], range_model[line_mask], color="tab:blue", lw=1.0)
    ax.text(0.03, 0.97, f"Range RMS {fmt_metric(plot_range_rms)} km\nI2 P16", transform=ax.transAxes, ha="left", va="top", fontsize=8, color="white")
    x_min = float(np.nanmin(t_plot))
    x_max = float(np.nanmax(t_plot))
    ylim_values = range_model
    if plot_range is not None:
        ylim_values = np.r_[ylim_values, plot_range[np.isfinite(plot_range)]]
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Range (km)"); ax.set_xlim(x_min - 0.05, x_max + 0.05); ax.set_ylim(max(0, np.nanmin(ylim_values) - 4), np.nanmax(ylim_values) + 4)
    ax = axs[1, 1]
    ax.plot(profile_radius_um, prob_rel, color="black", lw=1.0)
    ax.fill_between(profile_radius_um, 0, prob_rel, color="0.75", alpha=0.7, lw=0)
    ax.axvline(free_radius_um, color="tab:blue", lw=1.0)
    ax.axvspan(r0_q[0], r0_q[1], color="tab:blue", alpha=0.12, lw=0)
    ax.text(
        0.04,
        0.94,
        "95% "
        + r"$r_0$ "
        + interval_text(r0_q[0], r0_q[1], r" $\mu$m")
        + "\n95% "
        + r"$m_0$ "
        + interval_text(m0_q[0], m0_q[1], " kg", sci=True),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    ax.set_xscale("log"); ax.set_xlim(1, 10000); ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"Initial radius $r_0$ ($\mu$m)"); ax.set_ylabel("Relative probability"); ax.grid(alpha=0.25, lw=0.4, which="both")
    sec = ax.secondary_xaxis("top", functions=(r_um_to_mass_kg, mass_kg_to_r_um))
    sec.set_xscale("log"); sec.set_xlabel(r"Initial mass $m_0$ (kg)"); sec.tick_params(labelsize=7, pad=1)
    ax = axs[1, 2]
    if cepl_vel is not None:
        if speed_interval is not None:
            ax.fill_between(t_plot[line_mask], speed_interval[0, line_mask], speed_interval[1, line_mask], color="tab:blue", alpha=0.18, lw=0)
        ax.plot(t_plot[line_mask], np.linalg.norm(cepl_vel[line_mask], axis=1), color="tab:blue", lw=1.0, label="free fit")
    for target_um, fixed_speed in fixed_r0_speeds.items():
        fixed_mask = line_mask & np.isfinite(fixed_speed)
        if np.any(fixed_mask):
            ax.plot(
                t_plot[fixed_mask],
                fixed_speed[fixed_mask],
                color=fixed_colors.get(target_um, "0.35"),
                lw=0.9,
                ls="--",
                label=rf"$r_0={target_um:g}\,\mu$m",
            )
    ax.text(
        0.03,
        0.97,
        r"$\chi^2_\nu$ "
        + fmt_metric(stored_redchi, 2)
        + f"\n3D RMS {fmt_metric(stored_pos_rms)} km\nN {int(fit_n) if np.isfinite(fit_n) else len(t)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Speed (km/s)"); ax.grid(alpha=0.2, lw=0.4)
    if fixed_r0_speeds:
        ax.legend(loc="lower left", fontsize=7, frameon=False)
    axs[1, 3].axis("off")
    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.08, top=0.90, wspace=0.32, hspace=0.40)
    PLOTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out, profile_path

def main():
    global BASE, SRC_BATCH, OUT_BATCH, PLOTS, PROFILES, HIGHRES, PROFILE_INPUT_DIR, BEAT_H5_DIR

    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=10)
    ap.add_argument("--start-rank", type=int, default=1)
    ap.add_argument("--grid-n", type=int, default=41)
    ap.add_argument("--sample-index-h5", type=Path)
    ap.add_argument("--base", type=Path)
    ap.add_argument("--output-dir", type=Path)
    ap.add_argument("--profile-input-dir", type=Path)
    ap.add_argument("--beat-h5-dir", type=Path)
    ap.add_argument("--worker-index", type=int, default=0)
    ap.add_argument("--worker-count", type=int, default=1)
    args = ap.parse_args()
    if args.worker_count < 1 or not 0 <= args.worker_index < args.worker_count:
        ap.error("worker-index must satisfy 0 <= worker-index < worker-count")
    if args.base is not None:
        BASE = args.base
        SRC_BATCH = BASE / "events/current_distribution_plots/example_event_batch_200"
    if args.output_dir is not None:
        OUT_BATCH = args.output_dir
        PLOTS = OUT_BATCH / "plots"
        PROFILES = OUT_BATCH / "profiles"
        HIGHRES = OUT_BATCH / "highres_fft_i2_p16"
    if args.profile_input_dir is not None:
        PROFILE_INPUT_DIR = args.profile_input_dir
    if args.beat_h5_dir is not None:
        BEAT_H5_DIR = args.beat_h5_dir
    OUT_BATCH.mkdir(parents=True, exist_ok=True); PLOTS.mkdir(parents=True, exist_ok=True); PROFILES.mkdir(parents=True, exist_ok=True)
    rows = [
        row for row in load_rows(limit=args.limit, sample_index_h5=args.sample_index_h5)
        if row[0] >= args.start_rank
    ]
    rows = rows[args.worker_index::args.worker_count]
    if args.sample_index_h5 is None:
        (OUT_BATCH / "candidate_manifest.tsv").write_text((SRC_BATCH / "candidate_manifest.tsv").read_text())
    summary_path = OUT_BATCH / f"render_summary_worker_{args.worker_index:03d}.tsv"
    with summary_path.open("a") as sf:
        for row in rows:
            rank, sample, day = row
            try:
                out, prof = render(row, grid_n=args.grid_n)
                status = "ok"
                print(rank, sample, status, out, prof, flush=True)
            except Exception as exc:
                status = f"ERR {type(exc).__name__}: {exc}"
                print(rank, sample, status, flush=True)
                traceback.print_exc()
            sf.write(f"{rank}\t{sample}\t{status}\n")
            sf.flush()
if __name__ == "__main__":
    main()
