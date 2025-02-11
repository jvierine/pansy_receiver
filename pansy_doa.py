import argparse
import numpy as np
import scipy.optimize as sciopt
import matplotlib.pyplot as plt
from tqdm import tqdm
import pyant


def make_pansy_met_rec():
    beam = pyant.beam_of_radar("pansy", "array")

    groups = [
        "E1",
        "E2",
        "D6",
        "D8",
        "D3",
        "H2",
        "H4",
    ]

    keep = np.full((len(beam.antennas),), False, dtype=bool)
    for ind, grp in enumerate(beam.meta["groups"]):
        keep[ind] = str(grp) in groups

    beam.antennas = [grp for ind, grp in enumerate(beam.antennas) if keep[ind]]

    return beam


def plot_pansy_met_rec(args, beam):
    fig, ax = pyant.plotting.antenna_configuration(beam.antennas)
    ax.axis("equal")

    fig, ax = plt.subplots()
    pyant.plotting.gain_heatmap(
        beam,
        resolution=300,
        min_elevation=80.0,
        centered=False,
        ax=ax,
    )
    ax.set_title("PANSY radar")

    plt.show()


def plot_pansy_ambig_slice(args, beam: pyant.models.Array):
    beam.sph_point(0, 90, degrees=True)
    k_in = pyant.coordinates.sph_to_cart(
        np.array([0.0, 90.0, 1.0]),
        degrees=True,
    )

    resolution = 300

    S, K, k, inds, kx, ky = pyant.plotting.compute_k_grid(
        pointing=beam.pointing[:, 0],
        resolution=resolution,
        centered=False,
        cmin=1,
    )

    psi0 = beam.channel_signals(k_in)
    psi0 = psi0 / np.linalg.norm(psi0)

    psi_map = beam.channel_signals(k[:, inds])
    psi_map = psi_map / np.linalg.norm(psi_map, axis=0)[None, :]

    S[inds] = np.abs(np.sum(psi_map * np.conj(psi0)[:, None], axis=0))

    S = S.reshape(resolution, resolution)

    fig, ax = plt.subplots()
    pmesh = ax.pcolor(K[:, :, 0], K[:, :, 1], S, cmap="jet")
    fig.colorbar(pmesh, ax=ax)
    plt.show()


def calc_pansy_ambig_peaks(args, beam: pyant.models.Array):
    beam.sph_point(0, 90, degrees=True)
    k_in = pyant.coordinates.sph_to_cart(
        np.array([0.0, 85.0, 1.0]),
        degrees=True,
    )

    resolution = 500
    min_dist = 0.07

    S, K, k, inds, kx, ky = pyant.plotting.compute_k_grid(
        pointing=beam.pointing[:, 0],
        resolution=resolution,
        centered=False,
        cmin=1,
    )

    psi0 = beam.channel_signals(k_in)
    psi0 = psi0 / np.linalg.norm(psi0)

    psi_map = beam.channel_signals(k[:, inds])
    psi_map = psi_map / np.linalg.norm(psi_map, axis=0)[None, :]

    S[inds] = np.abs(np.sum(psi_map * np.conj(psi0)[:, None], axis=0))

    S_sort = np.argsort(S[inds])[::-1]
    p_dist = np.linalg.norm(k[:2, inds] - k_in[:2, None], axis=0)
    peak_ind = np.argmax(p_dist[S_sort] > min_dist)
    k_peak0 = k[:, inds][:, S_sort[peak_ind]]
    k_peak = k_peak0.copy()

    def min_func(k_flat):
        xy2 = k_flat[0] ** 2 + k_flat[1] ** 2
        if xy2 >= 1:
            return 0.0
        _k = np.array(
            [
                k_flat[0],
                k_flat[1],
                np.sqrt(1 - xy2),
            ]
        )
        psi_step = beam.channel_signals(_k)
        psi_step = psi_step / np.linalg.norm(psi_step)
        S = -np.abs(np.sum(psi_step * np.conj(psi0), axis=0))
        return S

    res = sciopt.minimize(
        min_func,
        x0=k_peak[:2],
        bounds=[(-1, 1), (-1, 1)],
    )
    k_peak[:2] = res.x
    k_peak[2] = np.sqrt(1 - k_peak[0]**2 - k_peak[1]**2)
    peak_S = -res.fun
    print(peak_S)

    S = S.reshape(resolution, resolution)

    fig, ax = plt.subplots()
    pmesh = ax.pcolor(K[:, :, 0], K[:, :, 1], S, cmap="bone")
    ax.plot(k_in[0], k_in[1], "xg", markersize=15)
    ax.plot(k_peak[0], k_peak[1], "or")
    fig.colorbar(pmesh, ax=ax)
    plt.show()


def sim_pansy_doa(args, beam: pyant.models.Array):
    beam.sph_point(0, 90, degrees=True)
    k_in = pyant.coordinates.sph_to_cart(
        np.array([0.0, 85.0, 1.0]),
        degrees=True,
    )

    SNRdBs = np.linspace(5, 30, num=50)
    k_err_vec = np.empty_like(SNRdBs)
    # k_errs_lim = 0.01
    for sni, SNRdB in tqdm(enumerate(SNRdBs), total=len(SNRdBs)):
        samples = 300
        # SNRdB = 100.0
        SNR = 10**(SNRdB/10.0)
        sigma_c = 1 / np.sqrt(SNR * 2 * beam.channels)
        rand_samps = np.random.randn(beam.channels, samples, 2)
        xi = sigma_c*(rand_samps[:, :, 0] + 1j * rand_samps[:, :, 1])

        resolution = 200
        S, K, k, inds, kx, ky = pyant.plotting.compute_k_grid(
            pointing=beam.pointing[:, 0],
            resolution=resolution,
            centered=False,
            cmin=np.cos(np.radians(60.0)),
        )

        psi0 = beam.channel_signals(k_in)
        psi0 = psi0 / np.linalg.norm(psi0)
        psi = psi0[:, None] + xi

        psi_map = beam.channel_signals(k[:, inds])
        psi_map = psi_map / np.linalg.norm(psi_map, axis=0)[None, :]

        k_finds = np.empty((3, samples), dtype=np.float64)
        k_errs = np.empty((samples, ), dtype=np.float64)
        # for ind in tqdm(range(samples), total=samples):
        for ind in range(samples):
            sig = np.conj(psi[:, ind])
            match = np.abs(np.sum(psi_map * sig[:, None], axis=0))
            max_ind = np.argmax(match)
            k_finds[:, ind] = k[:, inds][:, max_ind]
            k_errs[ind] = np.linalg.norm(k_in[:2] - k_finds[:2, ind])
        k_err_vec[sni] = np.mean(k_errs)

        # fig, axes = plt.subplots(1, 2)
        # axes[0].plot(k_in[0], k_in[1], ".r")
        # axes[0].plot(k_finds[0, :], k_finds[1, :], ".b")

        # axes[1].hist(k_errs)
        # plt.show()

    fig, ax = plt.subplots()
    ax.plot(SNRdBs, k_err_vec)
    plt.show()


task_map = {
    "plot-beam": plot_pansy_met_rec,
    "plot-ambig-slice": plot_pansy_ambig_slice,
    "calc-ambig-peaks": calc_pansy_ambig_peaks,
    "sim-doa-noise": sim_pansy_doa,
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("task", choices=list(task_map.keys()))
    args = parser.parse_args()

    beam = make_pansy_met_rec()

    task = task_map[args.task]
    task(args, beam)
