#!/usr/bin/env python3
"""Plot the omega-Eridanid orbit ensemble and solar-longitude activity."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from radiant_visibility import centered_plot_longitude_deg, radiant_exposure_hours_points


DEFAULT_CATALOGUE = Path(__file__).parent / "figs" / "pansy_maarsy_keplerian_catalogue.h5"
DEFAULT_EXPOSURE = Path(__file__).parent / "figs" / "paper_refresh_20260721_current" / "paper_radiant_results.h5"
DEFAULT_OUTPUT = Path.home() / "src" / "pansy_paper" / "paper_omega_eridanids.png"

SOLAR_RANGE_DEG = (100.0, 120.0)
CORE_SOLAR_RANGE_DEG = (108.2, 112.2)
VG_RANGE_KM_S = (44.76, 55.19)
SEMIMAJOR_AXIS_RANGE_AU = (1.90, 6.52)
HEALPIX_NSIDE = 32
RADIANT_PIXELS = np.asarray(
    [
        10524, 10525, 10527,
        10641, 10642, 10643, 10644, 10645,
        10752, 10753, 10754, 10755, 10756, 10757,
        10860, 10862, 10863, 10864, 10865, 10866,
        10966, 10969,
    ],
    dtype=np.int64,
)
MEAN_RA_DEG = 54.33
MEAN_DEC_DEG = -30.63
MEAN_VG_KM_S = 48.65
MEAN_KEPLER = np.asarray([3.4039, 0.7089, 91.29, 290.14, 317.33, 42.24, 0.8999])
MEAN_SC_LON_DEG = 291.04
MEAN_BETA_DEG = -48.21
OBLIQUITY_DEG = 23.4392911


def wrap180(values):
    return (np.asarray(values, dtype=np.float64) + 180.0) % 360.0 - 180.0


def circular_mean_std_deg(values: np.ndarray) -> tuple[float, float]:
    """Return the circular mean and sample-like angular standard deviation."""
    angle = np.deg2rad(np.asarray(values, dtype=np.float64))
    finite = np.isfinite(angle)
    angle = angle[finite]
    if angle.size == 0:
        return np.nan, np.nan
    vector = np.mean(np.exp(1j * angle))
    mean = np.mod(np.rad2deg(np.angle(vector)), 360.0)
    if np.isclose(mean, 360.0):
        mean = 0.0
    residual = wrap180(np.rad2deg(angle) - mean)
    std = np.std(residual, ddof=1) if residual.size > 1 else 0.0
    return float(mean), float(std)


def ecliptic_to_equatorial_deg(lon_deg, lat_deg) -> tuple[np.ndarray, np.ndarray]:
    """Rotate mean-ecliptic radiant coordinates to equatorial coordinates."""
    lon = np.deg2rad(np.asarray(lon_deg, dtype=np.float64))
    lat = np.deg2rad(np.asarray(lat_deg, dtype=np.float64))
    eps = np.deg2rad(OBLIQUITY_DEG)
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    y_eq = np.cos(eps) * y - np.sin(eps) * z
    z_eq = np.sin(eps) * y + np.cos(eps) * z
    ra = np.mod(np.rad2deg(np.arctan2(y_eq, x)), 360.0)
    dec = np.rad2deg(np.arcsin(np.clip(z_eq, -1.0, 1.0)))
    return ra, dec


def ang2pix_ring(nside: int, lon_deg, lat_deg) -> np.ndarray:
    """HEALPix RING angular index, matching Radview's @hscmap/healpix implementation."""
    z = np.sin(np.deg2rad(np.asarray(lat_deg, dtype=np.float64)))
    za = np.abs(z)
    tt = np.mod(np.asarray(lon_deg, dtype=np.float64), 360.0) / 90.0
    z, tt = np.broadcast_arrays(z, tt)
    out = np.empty(z.shape, dtype=np.int64)
    equatorial = za <= (2.0 / 3.0)
    ncap = 2 * nside * (nside - 1)

    temp1 = nside * (0.5 + tt[equatorial])
    temp2 = nside * z[equatorial] * 0.75
    jp = np.floor(temp1 - temp2).astype(np.int64)
    jm = np.floor(temp1 + temp2).astype(np.int64)
    ir = nside + 1 + jp - jm
    kshift = 1 - (ir & 1)
    ip = np.floor((jp + jm - nside + kshift + 1) / 2).astype(np.int64) + 1
    ip = np.mod(ip - 1, 4 * nside) + 1
    out[equatorial] = ncap + (ir - 1) * 4 * nside + ip - 1

    polar = ~equatorial
    tp = tt[polar] - np.floor(tt[polar])
    tmp = nside * np.sqrt(3.0 * (1.0 - za[polar]))
    jp = np.floor(tp * tmp).astype(np.int64)
    jm = np.floor((1.0 - tp) * tmp).astype(np.int64)
    ir = jp + jm + 1
    ip = np.floor(tt[polar] * ir).astype(np.int64) + 1
    npix = 12 * nside * nside
    north = z[polar] > 0.0
    polar_pixels = np.empty(ir.shape, dtype=np.int64)
    polar_pixels[north] = 2 * ir[north] * (ir[north] - 1) + ip[north] - 1
    polar_pixels[~north] = npix - 2 * ir[~north] * (ir[~north] + 1) + ip[~north] - 1
    out[polar] = polar_pixels
    return out


def load_catalogue(path: Path) -> dict[str, np.ndarray]:
    with h5py.File(path, "r") as h5:
        solar = np.asarray(h5["solar_longitude_deg"], dtype=np.float64)
        source = np.asarray(h5["source_id"], dtype=np.int8)
        preselect = (source == 0) & (solar >= SOLAR_RANGE_DEG[0]) & (solar <= SOLAR_RANGE_DEG[1])
        index = np.flatnonzero(preselect)
        return {
            "index": index,
            "solar": solar[index],
            "vg": np.asarray(h5["vg_km_s"][index], dtype=np.float64),
            "lon": np.asarray(h5["sun_centered_lon_deg"][index], dtype=np.float64),
            "ecliptic_lon": np.asarray(h5["ecliptic_lon_deg"][index], dtype=np.float64),
            "beta": np.asarray(h5["ecliptic_lat_deg"][index], dtype=np.float64),
            "kepler": np.asarray(h5["kepler"][index], dtype=np.float64),
        }


def selection_masks(rows: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    pixel = ang2pix_ring(HEALPIX_NSIDE, rows["lon"], rows["beta"])
    velocity = (rows["vg"] >= VG_RANGE_KM_S[0]) & (rows["vg"] <= VG_RANGE_KM_S[1])
    semimajor_axis = (
        (rows["kepler"][:, 0] >= SEMIMAJOR_AXIS_RANGE_AU[0])
        & (rows["kepler"][:, 0] <= SEMIMAJOR_AXIS_RANGE_AU[1])
    )
    distribution_filter = velocity & semimajor_axis
    shower = distribution_filter & np.isin(pixel, RADIANT_PIXELS)
    core = shower & (rows["solar"] >= CORE_SOLAR_RANGE_DEG[0]) & (rows["solar"] <= CORE_SOLAR_RANGE_DEG[1])
    return shower, core


def interpolate_observation_solar_longitude(h5: h5py.File) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    catalogue = h5["radiants"]
    meteor_epoch = np.asarray(catalogue["epoch_unix"], dtype=np.float64)
    meteor_solar = np.asarray(catalogue["sun_lambda_ecliptic_deg"], dtype=np.float64)
    order = np.argsort(meteor_epoch)
    meteor_epoch = meteor_epoch[order]
    meteor_solar = np.rad2deg(np.unwrap(np.deg2rad(meteor_solar[order])))
    observation_epoch = np.asarray(h5["observation_sample_idx"], dtype=np.float64) / float(
        h5.attrs["observation_sample_rate_hz"]
    )
    observation_hours = np.asarray(h5["observation_hours"], dtype=np.float64)
    observation_solar = np.mod(np.interp(observation_epoch, meteor_epoch, meteor_solar), 360.0)
    return observation_epoch, observation_solar, observation_hours


def activity_profile(
    rows: dict[str, np.ndarray],
    shower_mask: np.ndarray,
    exposure_path: Path,
    bin_width_deg: float,
) -> dict[str, np.ndarray]:
    edges = np.arange(SOLAR_RANGE_DEG[0], SOLAR_RANGE_DEG[1] + bin_width_deg, bin_width_deg)
    centers = 0.5 * (edges[:-1] + edges[1:])
    shower_count, _ = np.histogram(rows["solar"][shower_mask], bins=edges)

    with h5py.File(exposure_path, "r") as h5:
        observation_epoch, observation_solar, observation_hours = interpolate_observation_solar_longitude(h5)
    exposure = np.zeros_like(centers)
    plot_lon = float(centered_plot_longitude_deg(MEAN_SC_LON_DEG))
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        keep = (observation_solar >= lo) & (observation_solar < hi)
        exposure[i] = float(
            radiant_exposure_hours_points(
                observation_epoch[keep],
                observation_solar[keep],
                observation_hours[keep],
                np.asarray([plot_lon]),
                np.asarray([MEAN_BETA_DEG]),
            )[0]
        )

    rate = np.divide(shower_count, exposure, out=np.full_like(exposure, np.nan), where=exposure > 0.0)
    uncertainty = np.divide(
        np.sqrt(shower_count),
        exposure,
        out=np.full_like(exposure, np.nan),
        where=exposure > 0.0,
    )
    return {
        "centers": centers,
        "rate": rate,
        "uncertainty": uncertainty,
        "shower_count": shower_count,
        "exposure": exposure,
    }


def rotation_matrix(node_deg: float, inc_deg: float, argp_deg: float) -> np.ndarray:
    node, inc, argp = np.deg2rad([node_deg, inc_deg, argp_deg])
    cn, sn = np.cos(node), np.sin(node)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(argp), np.sin(argp)
    return np.asarray(
        [
            [cn * cw - sn * sw * ci, -cn * sw - sn * cw * ci, sn * si],
            [sn * cw + cn * sw * ci, -sn * sw + cn * cw * ci, -cn * si],
            [sw * si, cw * si, ci],
        ]
    )


def orbit_xyz(kepler: np.ndarray, samples: int = 600) -> np.ndarray:
    a, e, inc, node, argp, _nu, q = map(float, kepler)
    if not np.all(np.isfinite([a, e, inc, node, argp, q])) or e < 0.0 or q <= 0.0:
        return np.empty((3, 0))
    if e < 1.0 and a > 0.0:
        anomaly = np.linspace(0.0, 2.0 * np.pi, samples)
        p = a * (1.0 - e * e)
    else:
        limit = np.arccos(np.clip(-1.0 / max(e, 1.0 + 1e-9), -1.0, 1.0)) - 1e-3
        anomaly = np.linspace(-limit, limit, samples)
        p = q * (1.0 + e)
    radius = p / (1.0 + e * np.cos(anomaly))
    perifocal = np.vstack((radius * np.cos(anomaly), radius * np.sin(anomaly), np.zeros_like(radius)))
    return rotation_matrix(node, inc, argp) @ perifocal


def draw_orbit_panel(ax, orbits: np.ndarray, coordinates: tuple[int, int], limit_au: float) -> None:
    theta = np.linspace(0.0, 2.0 * np.pi, 720)
    ix, iy = coordinates
    for radius, label, color in [(1.0, "Earth", "C0"), (5.203, "Jupiter", "0.45")]:
        circle = np.vstack((radius * np.cos(theta), radius * np.sin(theta), np.zeros_like(theta)))
        ax.plot(circle[ix], circle[iy], color=color, lw=1.0, label=f"{label} orbit")
    for kepler in orbits:
        xyz = orbit_xyz(kepler)
        if xyz.shape[1]:
            ax.plot(xyz[ix], xyz[iy], color="0.55", lw=0.45, alpha=0.16)
    mean_xyz = orbit_xyz(MEAN_KEPLER, samples=1200)
    ax.plot(mean_xyz[ix], mean_xyz[iy], color="black", lw=2.4, label=r"Mean $\omega$-Eridanid orbit")
    ax.scatter(0.0, 0.0, s=55, color="#f5b82e", edgecolor="black", linewidth=0.5, zorder=5)
    ax.set_xlim(-limit_au, limit_au)
    ax.set_ylim(-limit_au, limit_au)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(color="0.90", lw=0.6)


def draw_mean_perifocal_panel(ax, orbits: np.ndarray) -> None:
    """Project into the mean orbit plane with its ecliptic intersection horizontal."""
    theta = np.linspace(0.0, 2.0 * np.pi, 720)
    mean_rotation = rotation_matrix(MEAN_KEPLER[3], MEAN_KEPLER[2], MEAN_KEPLER[4])
    argp = np.deg2rad(MEAN_KEPLER[4])
    node_rotation = np.asarray(
        [
            [np.cos(argp), -np.sin(argp)],
            [np.sin(argp), np.cos(argp)],
        ]
    )

    def project(xyz: np.ndarray) -> np.ndarray:
        return node_rotation @ (mean_rotation.T @ xyz)[:2]

    for radius, label, color in ((1.0, "Earth", "C0"), (5.2044, "Jupiter", "0.45")):
        ecliptic = np.vstack((radius * np.cos(theta), radius * np.sin(theta), np.zeros_like(theta)))
        projected = project(ecliptic)
        ax.plot(projected[0], projected[1], color=color, lw=1.0, label=f"{label} orbit")
    for kepler in orbits:
        xyz = orbit_xyz(kepler)
        if xyz.shape[1]:
            projected = project(xyz)
            ax.plot(projected[0], projected[1], color="0.55", lw=0.45, alpha=0.16)
    mean_projected = project(orbit_xyz(MEAN_KEPLER, samples=1200))
    ax.plot(mean_projected[0], mean_projected[1], color="black", lw=2.4, label=r"Mean $\omega$-Eridanid orbit")
    ax.scatter(0.0, 0.0, s=55, color="#f5b82e", edgecolor="black", linewidth=0.5, zorder=5)
    ax.set_xlim(-10.0, 5.0)
    ax.set_ylim(-8.5, 6.5)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(color="0.90", lw=0.6)


def make_figure(rows: dict[str, np.ndarray], core: np.ndarray, profile: dict[str, np.ndarray], output: Path) -> None:
    if int(np.sum(core)) != 76:
        raise RuntimeError(f"omega-Eridanid selection changed: expected 76 meteors, found {int(np.sum(core))}")
    orbits = rows["kepler"][core]
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6), constrained_layout=True)
    draw_orbit_panel(axes[0], orbits, (0, 1), 11.5)
    draw_mean_perifocal_panel(axes[1], orbits)
    axes[0].set_xlabel("Ecliptic X (AU)")
    axes[0].set_ylabel("Ecliptic Y (AU)")
    axes[1].set_xlabel("Ecliptic line of nodes (AU)")
    axes[1].set_ylabel("In-plane perpendicular (AU)")
    axes[0].set_title("Viewed from the north ecliptic pole")
    axes[1].set_title("Mean orbital plane")
    axes[0].legend(loc="lower left", frameon=False, fontsize=7.5)

    ra, dec = ecliptic_to_equatorial_deg(rows["ecliptic_lon"][core], rows["beta"][core])
    ra_mean, ra_std = circular_mean_std_deg(ra)
    linear_values = np.column_stack((dec, rows["vg"][core], orbits[:, 0], orbits[:, 1], orbits[:, 6]))
    linear_mean = np.nanmean(linear_values, axis=0)
    linear_std = np.nanstd(linear_values, axis=0, ddof=1)
    angular_stats = [circular_mean_std_deg(orbits[:, column]) for column in (2, 3, 4, 5)]

    summary = (
        r"$N=76$" "\n"
        rf"$\alpha_g={MEAN_RA_DEG:.2f}\pm{ra_std:.2f}^\circ$, "
        rf"$\delta_g={MEAN_DEC_DEG:.2f}\pm{linear_std[0]:.2f}^\circ$" "\n"
        rf"$v_g={MEAN_VG_KM_S:.2f}\pm{linear_std[1]:.2f}$ km s$^{{-1}}$" "\n"
        rf"$a={MEAN_KEPLER[0]:.2f}\pm{linear_std[2]:.2f}$ AU, "
        rf"$e={MEAN_KEPLER[1]:.3f}\pm{linear_std[3]:.3f}$" "\n"
        rf"$q={MEAN_KEPLER[6]:.3f}\pm{linear_std[4]:.3f}$ AU" "\n"
        rf"$i={MEAN_KEPLER[2]:.2f}\pm{angular_stats[0][1]:.2f}^\circ$, "
        rf"$\Omega={MEAN_KEPLER[3]:.2f}\pm{angular_stats[1][1]:.2f}^\circ$" "\n"
        rf"$\omega={MEAN_KEPLER[4]:.2f}\pm{angular_stats[2][1]:.2f}^\circ$, "
        rf"$\nu={MEAN_KEPLER[5]:.2f}\pm{angular_stats[3][1]:.2f}^\circ$"
    )
    axes[1].text(0.02, 0.02, summary, transform=axes[1].transAxes, ha="left", va="bottom", fontsize=7.6)

    ax = axes[2]
    x = profile["centers"]
    y = profile["rate"]
    dy = profile["uncertainty"]
    ax.axhline(0.0, color="0.5", lw=0.7)
    ax.fill_between(x, y - dy, y + dy, color="C0", alpha=0.18, linewidth=0)
    ax.plot(x, y, color="C0", marker="o", ms=3.2, lw=1.4, label="Selected rate (left axis)")
    ax.axvline(110.13, color="black", ls=":", lw=1.0)
    ax.set_xlim(*SOLAR_RANGE_DEG)
    ax.set_xlabel(r"Solar longitude, $\lambda_\odot$ (deg)")
    ax.set_ylabel(r"Exposure-corrected detected rate (h$^{-1}$)")
    ax.set_title("Activity profile")
    ax.grid(color="0.90", lw=0.6)

    count_axis = ax.twinx()
    count_axis.plot(
        x,
        profile["shower_count"],
        color="0.25",
        marker="s",
        markerfacecolor="none",
        markeredgewidth=0.8,
        ms=3.2,
        lw=0.8,
        ls=":",
        alpha=0.75,
        label="Raw count (right axis)",
    )
    count_axis.set_ylim(bottom=0.0)
    count_axis.set_ylabel("Raw selected count per bin", color="0.25")
    count_axis.tick_params(axis="y", colors="0.25")
    count_axis.spines["right"].set_color("0.25")
    rate_handles, rate_labels = ax.get_legend_handles_labels()
    count_handles, count_labels = count_axis.get_legend_handles_labels()
    ax.legend(
        rate_handles + count_handles,
        rate_labels + count_labels,
        loc="upper right",
        frameon=False,
        fontsize=8,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=DEFAULT_CATALOGUE)
    parser.add_argument("--exposure", type=Path, default=DEFAULT_EXPOSURE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bin-width-deg", type=float, default=1.0)
    args = parser.parse_args()

    rows = load_catalogue(args.catalogue)
    shower, core = selection_masks(rows)
    profile = activity_profile(rows, shower, args.exposure, args.bin_width_deg)
    make_figure(rows, core, profile, args.output)
    print(f"selected core meteors: {int(np.sum(core))}")
    print(f"selected 100-120 deg meteors: {int(np.sum(shower))}")
    print(args.output)


if __name__ == "__main__":
    main()
