#!/usr/bin/env python3
"""Inspect interferometer alias branches for one PANSY meteor cut."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import digital_rf as drf
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as c
import scipy.fftpack as fp
import scipy.ndimage as ndi
import h5py

import pansy_interferometry as pint
import pansy_config as pc
import pansy_modes as pmm


def amp_scale():
    with h5py.File("data/mesocal.h5", "r") as handle:
        pwr = handle["pwr"][()]
    return np.real(np.sqrt(pwr[0]) / np.sqrt(pwr[0:7]))


class RangeDopplerSearch:
    """Same matched-filter bank logic used by process_cut_meteor.py."""

    def __init__(self, txlen, echolen, n_channels, interp=1):
        self.txlen = txlen * interp
        self.fdec = 8 * interp
        self.fftlen = 512
        if int(self.txlen * interp / self.fdec) > self.fftlen:
            self.fftlen = int(self.txlen * interp / self.fdec)
        self.echolen = echolen * interp
        self.n_channels = n_channels
        self.interp = interp
        self.n_rg = self.echolen - self.txlen
        self.rg = np.arange(self.n_rg, dtype=np.int64)
        self.idx_mat = np.zeros((self.n_rg, self.txlen), dtype=np.int64)
        self.idx = np.arange(self.txlen, dtype=np.int64)
        for ri in range(self.n_rg):
            self.idx_mat[ri, :] = self.idx + self.rg[ri]
        self.z_rx = np.zeros((self.n_channels, self.echolen), dtype=np.complex64)
        self.rangev = self.rg * c.c / 2 / 1e6 / 1e3
        self.fvec = np.fft.fftshift(np.fft.fftfreq(self.fftlen, d=self.fdec / (interp * 1e6)))
        self.dopv = self.fvec * c.c / 2.0 / pc.freq
        self.ch_pairs = list(itertools.combinations(np.arange(n_channels), 2))
        self.n_pairs = len(self.ch_pairs)

    def decim(self, z):
        new_width = int(self.txlen / self.fdec)
        z2 = np.zeros((self.n_rg, new_width), dtype=np.complex64)
        idx = np.arange(new_width, dtype=np.int64)
        idx2 = np.arange(new_width, dtype=np.int64) * self.fdec
        for i in range(self.fdec):
            z2[:, idx] += z[:, idx2 + i]
        return z2

    def mf(self, tx, echo):
        txi = np.repeat(tx, self.interp)
        for chi in range(self.n_channels):
            self.z_rx[chi, :] = np.repeat(echo[chi, :], self.interp)
        mf = np.zeros((self.n_rg, self.fftlen), dtype=np.float32)
        xc = np.zeros((self.n_pairs, self.n_rg, self.fftlen), dtype=np.complex64)

        for chi in range(self.n_channels):
            z = self.z_rx[chi, self.idx_mat] * txi[None, :]
            zd = self.decim(z)
            zf = np.fft.fftshift(fp.fft(zd, self.fftlen, axis=1), axes=1)
            mf += zf.real**2.0 + zf.imag**2.0

        for chip, pair in enumerate(self.ch_pairs):
            z1 = self.z_rx[pair[0], self.idx_mat] * txi[None, :]
            z2 = self.z_rx[pair[1], self.idx_mat] * txi[None, :]
            zd1 = self.decim(z1)
            zd2 = self.decim(z2)
            zf1 = np.fft.fftshift(fp.fft(zd1, self.fftlen, axis=1), axes=1)
            zf2 = np.fft.fftshift(fp.fft(zd2, self.fftlen, axis=1), axes=1)
            xc[chip, :, :] = zf1 * np.conj(zf2)

        noise_floor = np.median(mf)
        pprof = np.max(mf, axis=1)
        dop_idx = np.argmax(mf, axis=1)
        rgmax = int(np.argmax(pprof))
        xc_peak = xc[:, rgmax, dop_idx[rgmax]]
        peak_dopv = self.dopv[np.argmax(mf, axis=1)]
        return mf, pprof, peak_dopv, noise_floor, xc_peak, rgmax


def load_cut(cut_dir: Path, sample_idx: int) -> dict:
    dm = drf.DigitalMetadataReader(str(cut_dir))
    records = dm.read(sample_idx - 1, sample_idx + 1)
    if sample_idx not in records:
        raise RuntimeError(f"cut sample {sample_idx} not found in {cut_dir}")
    return records[sample_idx]


def recompute_cut_observables(cut: dict, interp: int = 1):
    z_rx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(cut["zrx_echoes_im"], dtype=np.complex64)
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.complex64)
    scale = amp_scale()
    for i in range(z_rx.shape[1]):
        z_rx[:, i, :] *= scale[i]

    txlen = z_tx.shape[1]
    echolen = z_rx.shape[2]
    rds = RangeDopplerSearch(txlen=txlen, echolen=echolen, interp=interp, n_channels=z_rx.shape[1])

    delays = np.asarray(cut["delays"], dtype=np.int64)
    tx_idx = np.asarray(cut["tx_idx"], dtype=np.float64)
    beam_id = np.asarray(cut["beam_id"], dtype=np.int64)

    peak_rg = []
    peak_dop = []
    snr = []
    xct = []
    rti_snr = np.zeros((z_rx.shape[0], 1600), dtype=np.float32)
    dti_mps = np.zeros((z_rx.shape[0], 1600), dtype=np.float32)
    drg = c.c / 1e6 / 2 / 1e3
    ridx = np.arange(rds.n_rg, dtype=int)
    for ti in range(z_rx.shape[0]):
        _mf, pprof, peak_dopv, noise_floor, xc, _rgmax = rds.mf(z_tx[ti, :], z_rx[ti, :, :])
        max_rg = int(np.argmax(pprof))
        peak_rg.append((delays[ti] + max_rg / interp) * drg)
        peak_dop.append(peak_dopv[max_rg])
        snr.append((pprof[max_rg] - noise_floor) / noise_floor)
        xct.append(xc)
        rti_snr[ti, ridx + delays[ti]] = (pprof - noise_floor) / noise_floor
        dti_mps[ti, ridx + delays[ti]] = peak_dopv

    return {
        "tx_idx": tx_idx,
        "beam_id": beam_id,
        "range_km": np.asarray(peak_rg),
        "doppler_mps": np.asarray(peak_dop),
        "snr": np.asarray(snr),
        "xc": np.asarray(xct),
        "rti_snr": rti_snr,
        "dti_mps": dti_mps,
        "range_grid_km": np.arange(1600, dtype=np.float64) * drg,
        "radar_frequency_hz": float(pc.freq),
    }


def local_maxima_2d(values: np.ndarray, n_keep: int):
    neighborhood = ndi.maximum_filter(values, size=9, mode="nearest")
    mask = values == neighborhood
    idx = np.flatnonzero(mask.ravel())
    if len(idx) == 0:
        idx = np.arange(values.size)
    order = np.argsort(values.ravel()[idx])[::-1]
    return idx[order[:n_keep]]


def top_alias_peaks(xc_row, beam, phasecal, ch_pairs, dmat, uu, vv, ww, n_keep=8):
    z = np.exp(1j * (np.angle(xc_row) + phasecal[beam, ch_pairs[:, 0]] - phasecal[beam, ch_pairs[:, 1]]))
    match = np.abs(pint.mf(z, dmat, uu, vv, ww))
    flat_idx = local_maxima_2d(match, n_keep)
    return [
        {
            "u": float(uu.ravel()[idx]),
            "v": float(vv.ravel()[idx]),
            "w": float(ww.ravel()[idx]),
            "match": float(match.ravel()[idx]),
        }
        for idx in flat_idx
    ]


def tx_array_positions() -> np.ndarray:
    """Ready transmit antenna positions from the PANSY array description."""
    pos = []
    for serial, ant in pc.antenna.items():
        if serial == "RFTX":
            continue
        try:
            ready = int(ant["ready"]) == 1
        except ValueError:
            ready = str(ant["ready"]).strip().lower() in {"true", "yes", "ready"}
        if ready:
            pos.append([ant["x"], ant["y"], ant["z"]])
    if not pos:
        raise RuntimeError("no ready transmit antennas found in cfg/antpos.csv")
    return np.asarray(pos, dtype=np.float64)


def beam_unit_vectors() -> np.ndarray:
    """Beam pointing vectors in the same u/v/w convention used by imaging."""
    mode = pmm.get_m_mode()
    out = []
    for az_deg, za_deg in mode["beam_pos_az_za"]:
        el_deg = 90.0 - za_deg
        p_h = np.cos(np.deg2rad(el_deg))
        w = -np.sin(np.deg2rad(el_deg))
        v = p_h * np.cos(-np.deg2rad(az_deg))
        u = -p_h * np.sin(-np.deg2rad(az_deg))
        out.append([u, v, w])
    return np.asarray(out, dtype=np.float64)


def tx_array_gain_db(uvw: np.ndarray, beam_id: np.ndarray, tx_pos: np.ndarray, beam_vecs: np.ndarray) -> np.ndarray:
    """Relative steered transmit array-factor power gain in dB.

    The gain is normalized so that a target exactly at the commanded beam
    pointing direction has 0 dB gain. Sidelobes and grating lobes are retained
    because the full ready-antenna geometry is used.
    """
    uvw = np.asarray(uvw, dtype=np.float64)
    beam_id = np.asarray(beam_id, dtype=np.int64)
    k0 = 2.0 * np.pi / pc.wavelength
    gain = np.full(len(uvw), np.nan, dtype=np.float64)
    for i, direction in enumerate(uvw):
        if not np.all(np.isfinite(direction)):
            continue
        steer = beam_vecs[beam_id[i]]
        phase = k0 * (tx_pos @ (direction - steer))
        af = np.abs(np.mean(np.exp(1j * phase)))
        gain[i] = 20.0 * np.log10(max(af, 1e-8))
    return gain


def beam_snr_consistency(obs: dict, branch: dict, tx_pos: np.ndarray, beam_vecs: np.ndarray) -> dict:
    """Compare range-normalized matched-filter SNR to TX array-factor gain."""
    keep = branch["keep"]
    if np.count_nonzero(keep) < 8:
        return {"beam_rms_db": np.nan, "beam_corr": np.nan, "beam_gain_span_db": np.nan}

    gain_db = tx_array_gain_db(branch["uvw"][keep], obs["beam_id"][keep], tx_pos, beam_vecs)
    snr_db = 10.0 * np.log10(np.maximum(obs["snr"][keep], 1e-6))
    range_km = np.maximum(obs["range_km"][keep], 1e-6)
    snr_corr_db = snr_db + 40.0 * np.log10(range_km / np.median(range_km))

    good = np.isfinite(gain_db) & np.isfinite(snr_corr_db)
    if np.count_nonzero(good) < 8:
        return {"beam_rms_db": np.nan, "beam_corr": np.nan, "beam_gain_span_db": np.nan}

    gain_rel = gain_db[good] - np.median(gain_db[good])
    snr_rel = snr_corr_db[good] - np.median(snr_corr_db[good])
    resid = snr_rel - gain_rel
    corr = np.corrcoef(snr_rel, gain_rel)[0, 1] if np.std(gain_rel) > 1e-6 and np.std(snr_rel) > 1e-6 else np.nan
    out = {
        "tx_gain_db": np.full(len(keep), np.nan, dtype=np.float64),
        "range_corrected_snr_db": np.full(len(keep), np.nan, dtype=np.float64),
        "beam_rms_db": float(np.sqrt(np.mean(resid**2))),
        "beam_corr": float(corr),
        "beam_gain_span_db": float(np.nanpercentile(gain_db[good], 95) - np.nanpercentile(gain_db[good], 5)),
    }
    out["tx_gain_db"] = tx_array_gain_db(branch["uvw"], obs["beam_id"], tx_pos, beam_vecs)
    out["range_corrected_snr_db"][keep] = snr_db + 40.0 * np.log10(range_km / np.median(range_km))
    return out


def cluster_alias_centers(peaks, radius=0.035, min_weight=50.0):
    centers = []
    for peak in sorted(peaks, key=lambda p: p["weight"], reverse=True):
        uv = np.array([peak["u"], peak["v"]])
        assigned = False
        for center in centers:
            if np.linalg.norm(uv - center["uv"]) < radius:
                w0 = center["weight"]
                center["uv"] = (center["uv"] * w0 + uv * peak["weight"]) / (w0 + peak["weight"])
                center["weight"] += peak["weight"]
                center["n"] += 1
                assigned = True
                break
        if not assigned:
            centers.append({"uv": uv, "weight": peak["weight"], "n": 1})
    centers = [c0 for c0 in centers if c0["weight"] >= min_weight]
    centers.sort(key=lambda c0: c0["weight"], reverse=True)
    return centers


def fit_branch(obs, per_pulse_peaks, center, max_dist=0.06):
    t = obs["tx_idx"] / 1e6
    t_rel = t - t[0]
    points = []
    uvw = []
    keep = []
    matches = []
    for rg, peaks in zip(obs["range_km"], per_pulse_peaks):
        if not peaks:
            keep.append(False)
            points.append([np.nan, np.nan, np.nan])
            uvw.append([np.nan, np.nan, np.nan])
            matches.append(np.nan)
            continue
        uv0 = center["uv"]
        dist = np.asarray([np.linalg.norm(np.array([p["u"], p["v"]]) - uv0) for p in peaks])
        ii = int(np.argmin(dist))
        p = peaks[ii]
        good = bool(dist[ii] < max_dist)
        keep.append(good)
        los_up = np.sqrt(max(0.0, 1.0 - p["u"] ** 2 - p["v"] ** 2))
        points.append([rg * p["u"], rg * p["v"], rg * los_up])
        uvw.append([p["u"], p["v"], p["w"]])
        matches.append(p["match"])

    points = np.asarray(points, dtype=float)
    uvw = np.asarray(uvw, dtype=float)
    keep = np.asarray(keep, dtype=bool)
    matches = np.asarray(matches, dtype=float)
    if np.count_nonzero(keep) < 8:
        return None

    A = np.column_stack([np.ones(np.count_nonzero(keep)), t_rel[keep]])
    coeff = np.linalg.lstsq(A, points[keep], rcond=None)[0]
    model = np.column_stack([
        coeff[0, 0] + coeff[1, 0] * t_rel,
        coeff[0, 1] + coeff[1, 1] * t_rel,
        coeff[0, 2] + coeff[1, 2] * t_rel,
    ])
    v_km_s = coeff[1]
    rg_model = np.linalg.norm(model, axis=1)
    los_rate_km_s = np.sum(model * v_km_s[None, :], axis=1) / np.maximum(rg_model, 1e-6)
    dop_km_s = obs["doppler_mps"] / 1e3
    dop_resid = los_rate_km_s[keep] - dop_km_s[keep]
    pos_resid = points[keep] - model[keep]
    line_pts = points[keep]
    line_center = np.mean(line_pts, axis=0)
    _, _, vh = np.linalg.svd(line_pts - line_center, full_matrices=False)
    line_dir = vh[0]
    line_perp = (line_pts - line_center) - np.outer((line_pts - line_center) @ line_dir, line_dir)
    v_dir = v_km_s / np.maximum(np.linalg.norm(v_km_s), 1e-12)
    transverse_resid = pos_resid - np.outer(pos_resid @ v_dir, v_dir)
    return {
        "center": center,
        "keep": keep,
        "points": points,
        "uvw": uvw,
        "matches": matches,
        "model": model,
        "v_km_s": v_km_s,
        "speed_km_s": float(np.linalg.norm(v_km_s)),
        "los_rate_km_s": los_rate_km_s,
        "doppler_rms_km_s": float(np.sqrt(np.mean(dop_resid**2))),
        "pos_rms_km": float(np.sqrt(np.mean(pos_resid**2))),
        "line_rms_km": float(np.sqrt(np.mean(np.sum(line_perp**2, axis=1)))),
        "transverse_rms_km": float(np.sqrt(np.mean(np.sum(transverse_resid**2, axis=1)))),
        "n": int(np.count_nonzero(keep)),
    }


def plot_aliases(obs, branches, all_peak_cloud, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    t_rel = obs["tx_idx"] / 1e6 - obs["tx_idx"][0] / 1e6
    fig, axes = plt.subplots(4, 2, figsize=(13, 17), constrained_layout=True)

    ax = axes[0, 0]
    rti_db = 10.0 * np.log10(np.maximum(obs["rti_snr"], 1e-6))
    vmax = np.nanpercentile(rti_db[np.isfinite(rti_db)], 99.7)
    im = ax.pcolormesh(
        t_rel,
        obs["range_grid_km"],
        rti_db.T,
        shading="auto",
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
    )
    ax.scatter(t_rel, obs["range_km"], c=10.0 * np.log10(obs["snr"]), s=8, cmap="magma", edgecolors="none")
    ax.set_ylabel("Range (km)")
    ax.set_title("Range-Doppler matched-filter RTI")
    fig.colorbar(im, ax=ax, label="Matched-filter SNR (dB)")

    ax = axes[0, 1]
    ax.scatter(t_rel, obs["doppler_mps"] / 1e3, c=10.0 * np.log10(obs["snr"]), s=10, cmap="magma")
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Doppler/range rate (km/s)")
    ax.set_title(f"Peak Doppler from matched filter, f={obs['radar_frequency_hz']/1e6:.3f} MHz")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    sc = ax.scatter([p["u"] for p in all_peak_cloud], [p["v"] for p in all_peak_cloud], c=[p["weight"] for p in all_peak_cloud], s=4, cmap="viridis", alpha=0.35)
    for bi, branch in enumerate(branches):
        uv = branch["center"]["uv"]
        ax.scatter([uv[0]], [uv[1]], s=100, marker="x", label=f"branch {bi}")
    ax.set_xlabel("u")
    ax.set_ylabel("v")
    ax.set_title("Interferometer alias peak cloud")
    ax.set_aspect("equal")
    ax.legend(fontsize=8)
    fig.colorbar(sc, ax=ax, label="SNR-weighted phase match")

    ax = axes[1, 1]
    ax.scatter(t_rel, obs["doppler_mps"] / 1e3, s=10, color="black", label="measured Doppler")
    for bi, branch in enumerate(branches[:8]):
        ax.plot(t_rel, branch["los_rate_km_s"], lw=1.2, label=f"{bi}: {branch['doppler_rms_km_s']:.1f} km/s")
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Range rate (km/s)")
    ax.set_title("Alias branch Doppler consistency")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=2)

    ax = axes[2, 0]
    for bi, branch in enumerate(branches[:8]):
        pts = branch["points"]
        keep = branch["keep"]
        ax.plot(pts[keep, 0], pts[keep, 1], ".", ms=2, label=f"{bi}: |v|={branch['speed_km_s']:.1f}")
    ax.set_xlabel("East (km)")
    ax.set_ylabel("North (km)")
    ax.set_title("Candidate alias tracks")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)

    ax = axes[2, 1]
    labels = [str(i) for i in range(len(branches))]
    dop = [b["doppler_rms_km_s"] for b in branches]
    pos = [b["line_rms_km"] for b in branches]
    speed = [b["speed_km_s"] for b in branches]
    ax.scatter(dop, speed, c=pos, s=70, cmap="magma")
    for i, (x, y) in enumerate(zip(dop, speed)):
        ax.text(x, y, labels[i])
    ax.set_xlabel("Doppler RMS mismatch (km/s)")
    ax.set_ylabel("Branch speed (km/s)")
    ax.set_title("Doppler and straight-line score")
    ax.grid(True, alpha=0.3)
    fig.colorbar(ax.collections[0], ax=ax, label="Line RMS (km)")

    ax = axes[3, 0]
    for bi, branch in enumerate(branches[:8]):
        keep = branch["keep"]
        if "range_corrected_snr_db" not in branch:
            continue
        snr_db = branch["range_corrected_snr_db"][keep]
        gain_db = branch["tx_gain_db"][keep]
        if len(snr_db) == 0:
            continue
        snr_rel = snr_db - np.nanmedian(snr_db)
        gain_rel = gain_db - np.nanmedian(gain_db)
        ax.plot(t_rel[keep], snr_rel, ".", ms=2, alpha=0.5)
        ax.plot(t_rel[keep], gain_rel, lw=1.0, label=f"{bi}: {branch['beam_rms_db']:.1f} dB")
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Relative dB")
    ax.set_title("Range-normalized SNR vs TX beam gain")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=2)

    ax = axes[3, 1]
    beam = [b.get("beam_rms_db", np.nan) for b in branches]
    line = [b["line_rms_km"] for b in branches]
    sc = ax.scatter(line, beam, c=dop, s=70, cmap="viridis")
    for i, (x, y) in enumerate(zip(line, beam)):
        ax.text(x, y, labels[i])
    ax.set_xlabel("Straight-line RMS (km)")
    ax.set_ylabel("SNR/beam RMS mismatch (dB)")
    ax.set_title("Alias branch consistency checks")
    ax.grid(True, alpha=0.3)
    fig.colorbar(sc, ax=ax, label="Doppler RMS (km/s)")

    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot interferometer alias branches for one PANSY cut.")
    parser.add_argument("--sample-idx", type=int, default=1746494797212678)
    parser.add_argument("--cut-dir", type=Path, default=Path("data/metadata/cut"))
    parser.add_argument("--output", type=Path, default=Path("plots/orbit_determination_example/aliases/alias_branches_1746494797212678.png"))
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--grid-n", type=int, default=360)
    args = parser.parse_args()

    cut = load_cut(args.cut_dir, args.sample_idx)
    obs_all = recompute_cut_observables(cut)
    good = obs_all["snr"] > args.snr_threshold
    obs = {
        key: val[good] if isinstance(val, np.ndarray) and len(val) == len(good) else val
        for key, val in obs_all.items()
    }

    antpos = pint.get_antpos()
    ch_pairs = np.asarray(list(itertools.combinations(np.arange(7), 2)))
    dmat = pint.pair_mat(ch_pairs, antpos)
    phasecal = pint.get_phasecal()
    uu, vv, ww = pint.uv_coverage(N=args.grid_n, max_zenith_angle=35.0)

    per_pulse_peaks = []
    peak_cloud = []
    for xc, beam, snr in zip(obs["xc"], obs["beam_id"], obs["snr"]):
        peaks = top_alias_peaks(xc, int(beam), phasecal, ch_pairs, dmat, uu, vv, ww, n_keep=10)
        per_pulse_peaks.append(peaks)
        for peak in peaks:
            peak_cloud.append({"u": peak["u"], "v": peak["v"], "weight": peak["match"] * float(snr)})

    centers = cluster_alias_centers(peak_cloud, radius=0.035, min_weight=500.0)[:20]
    tx_pos = tx_array_positions()
    beam_vecs = beam_unit_vectors()
    branches = []
    for center in centers:
        branch = fit_branch(obs, per_pulse_peaks, center)
        if branch is not None:
            branch.update(beam_snr_consistency(obs, branch, tx_pos, beam_vecs))
            branches.append(branch)
    branches.sort(key=lambda b: (b["line_rms_km"], b["doppler_rms_km_s"], b.get("beam_rms_db", np.inf)))
    plot_aliases(obs, branches, peak_cloud, args.output)

    print(args.output)
    print("branch line_rms_km transverse_rms_km doppler_rms_km_s beam_rms_db beam_corr speed_km_s u v n")
    for i, b in enumerate(branches[:12]):
        print(
            i,
            f"{b['line_rms_km']:.3f}",
            f"{b['transverse_rms_km']:.3f}",
            f"{b['doppler_rms_km_s']:.3f}",
            f"{b.get('beam_rms_db', np.nan):.2f}",
            f"{b.get('beam_corr', np.nan):.2f}",
            f"{b['speed_km_s']:.3f}",
            f"{b['center']['uv'][0]:.4f}",
            f"{b['center']['uv'][1]:.4f}",
            b["n"],
        )


if __name__ == "__main__":
    main()
