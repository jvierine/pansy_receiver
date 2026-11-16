"""Conditional joint range/Doppler/interferometric fits and radius profiles.

Profile intervals use a working Gaussian compressed-observable likelihood.
They exclude array orientation, alias order, atmosphere and calibration-systematic
uncertainty; they are not marginal posteriors or validated coverage statements.
"""

import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import solve_triangular
from .native import propagate

OTHER = np.array([0, 1, 2, 4, 5, 6, 7, 8, 9])
REF = 3


def radial(position, velocity):
    return np.sum(position * velocity, axis=1) / np.linalg.norm(position, axis=1)


class EventFit:
    def __init__(
        self, measurement, array, station, atmosphere, phasecal, calibration_std
    ):
        self.m = measurement
        self.array = array
        self.station = station
        self.alt, self.rho = atmosphere
        self.t = measurement["time_s"]
        self.n = len(self.t)
        self.a = measurement["amplitude"] * np.exp(-1j * phasecal)
        self.phase = np.angle(self.a[:, OTHER] * self.a[:, REF, None].conj())
        c = measurement["noise_covariance"]
        raw = measurement["amplitude"]
        phasecov = []
        diff = np.eye(10)[OTHER] - np.eye(10)[REF]
        for a in raw:
            rotation = np.exp(-1j * np.angle(a))
            cov = (c * rotation[:, None] * rotation.conj()[None, :]).real / (
                2 * np.maximum(abs(a[:, None]) * abs(a[None, :]), 1.0)
            )
            cov = diff @ cov @ diff.T
            val, vec = np.linalg.eigh(cov)
            cov = (vec * np.maximum(val, 1e-5)) @ vec.T
            phasecov.append(np.linalg.inv(np.linalg.cholesky(cov)))
        self.phasewhite = np.array(phasecov)
        self.rvwhite = np.array(
            [
                np.linalg.inv(np.linalg.cholesky(c))
                for c in measurement["range_velocity_covariance"]
            ]
        )
        self.calstd = np.maximum(
            np.sqrt(calibration_std[OTHER] ** 2 + calibration_std[REF] ** 2), 0.03
        )
        self.phase_scale = 1.0
        self.rv_scale = np.ones(2)
        self.keep = np.ones(self.n, bool)
        self.ppp = None
        self.bounds = (
            np.r_[[-15.0, -15.0, 60.0, -100.0, -100.0, -100.0], np.full(9, -2.0)],
            np.r_[[15.0, 15.0, 160.0, 100.0, 100.0, -0.01], np.full(9, 2.0)],
        )

    def model(self, p, radius=None):
        state = np.asarray(p[:6]) * 1000
        if radius is None:
            position = state[:3] + self.t[:, None] * state[3:]
            velocity = np.broadcast_to(state[3:], position.shape)
            r = np.full(self.n, np.nan)
        else:
            y = propagate(
                state, radius, self.t, self.alt, self.rho, self.station.altitude_m
            )
            position = y[:, :3]
            velocity = y[:, 3:6]
            r = y[:, 6]
        return position, velocity, r

    def residual(self, p, radius=None, use_ppp=True):
        count = (
            int(np.sum(self.keep)) * 11
            + 9
            + (0 if self.ppp is None or not use_ppp else int(np.sum(self.ppp["keep"])))
        )
        try:
            pos, vel, _ = self.model(p, radius)
            ranges = np.linalg.norm(pos, axis=1)
            uv = pos[:, :2] / ranges[:, None]
            response = self.array(uv)
            if not np.all(np.isfinite(response)) or np.any(abs(response) < 1e-8):
                raise ValueError("outside array-response domain")
            predicted = (
                np.angle(response[:, OTHER] * response[:, REF, None].conj())
                + p[6:][None, :]
            )
            phase = np.angle(np.exp(1j * (self.phase - predicted)))
            pr = np.einsum("nij,nj->ni", self.phasewhite, phase) / self.phase_scale
            rv = np.column_stack(
                [
                    ranges - self.m["range_m"],
                    radial(pos, vel) - self.m["radial_velocity_mps"],
                ]
            )
            rv = np.einsum("nij,nj->ni", self.rvwhite, rv) / self.rv_scale
            pieces = [rv[self.keep].ravel(), pr[self.keep].ravel(), p[6:] / self.calstd]
            if self.ppp is not None and use_ppp:
                pair = self.ppp
                i = pair["previous"]
                j = pair["current"]
                # Mean range rate over the actual pair, not instantaneous geometric acceleration.
                v = (ranges[j] - ranges[i]) / (self.t[j] - self.t[i])
                pieces.append(
                    ((v - pair["velocity_mps"]) / pair["sigma_mps"])[pair["keep"]]
                )
            result = np.concatenate(pieces)
            if not np.all(np.isfinite(result)):
                raise ValueError("nonfinite residual")
            return result
        except ValueError:
            return np.full(count, 1e6)

    def solve(self, p, radius=None, use_ppp=True, robust=False, max_nfev=100):
        return least_squares(
            lambda x: self.residual(x, radius, use_ppp),
            np.clip(
                p, np.array(self.bounds[0]) + 1e-7, np.array(self.bounds[1]) - 1e-7
            ),
            bounds=self.bounds,
            x_scale="jac",
            max_nfev=max_nfev,
            loss="soft_l1" if robust else "linear",
            ftol=2e-6,
            xtol=2e-6,
            gtol=2e-5,
            diff_step=1e-4,
        )

    def initialize(self, uv):
        points = (
            np.column_stack([uv, np.sqrt(np.maximum(1 - np.sum(uv**2, axis=1), 0))])
            * self.m["range_m"][:, None]
        )
        coef = np.linalg.lstsq(
            np.column_stack([np.ones(self.n), self.t]), points / 1000, rcond=None
        )[0]
        state = np.r_[coef[0], coef[1]]
        if state[5] > -0.1:
            state[5] = np.median(self.m["radial_velocity_mps"]) / 1000
        candidates = []
        for transverse in [state[3:5], np.zeros(2)]:
            seed = state.copy()
            seed[3:5] = transverse
            candidates.append(
                self.solve(np.r_[seed, np.zeros(9)], robust=True, max_nfev=150)
            )
        fit = min(candidates, key=lambda f: np.sum(np.minimum(f.fun**2, 100)))
        # Calibrate working scatter once; do not change it separately at each radius.
        residual = self.residual(fit.x, use_ppp=False)
        rv = residual[: self.n * 2].reshape(self.n, 2)
        self.rv_scale = np.maximum(
            1.0, 1.4826 * np.median(abs(rv - np.median(rv, axis=0)), axis=0)
        )
        ph = residual[self.n * 2 : self.n * 11]
        self.phase_scale = max(1.0, float(1.4826 * np.median(abs(ph - np.median(ph)))))
        fit = self.solve(fit.x, robust=True, max_nfev=100)
        residual = self.residual(fit.x, use_ppp=False)[: self.n * 2].reshape(self.n, 2)
        self.keep = np.max(abs(residual), axis=1) < 8
        if np.sum(self.keep) < 8:
            self.keep[:] = True
        return self.solve(fit.x, max_nfev=100)

    def profile(self, seed, radii, use_ppp=True):
        parameters = []
        cost = []
        success = []
        state = seed.copy()
        for radius in radii[::-1]:
            fit = self.solve(state, float(radius), use_ppp=use_ppp, max_nfev=70)
            valid = np.max(abs(fit.fun)) < 1e5 and np.all(np.isfinite(fit.x))
            parameters.append(fit.x)
            cost.append(np.sum(fit.fun**2) if valid else np.inf)
            success.append(valid and fit.success)
            if valid:
                state = fit.x
        return dict(
            radius_m=np.asarray(radii),
            parameters=np.array(parameters[::-1]),
            chi2=np.array(cost[::-1]),
            optimizer_success=np.array(success[::-1]),
        )


def profile_bounds(profile, threshold=3.84):
    radius = profile["radius_m"]
    delta = profile["chi2"] - np.min(profile["chi2"])
    k = int(np.argmin(delta))
    best = radius[k]
    # Connected component around the best fit; retain open endpoint flags.
    lo = k
    hi = k
    while lo > 0 and delta[lo - 1] <= threshold:
        lo -= 1
    while hi < len(radius) - 1 and delta[hi + 1] <= threshold:
        hi += 1
    lower = (
        radius[lo]
        if lo == 0
        else np.exp(
            np.interp(
                threshold,
                [delta[lo], delta[lo - 1]],
                np.log([radius[lo], radius[lo - 1]]),
            )
        )
    )
    upper = (
        np.inf
        if hi == len(radius) - 1
        else np.exp(
            np.interp(
                threshold,
                [delta[hi], delta[hi + 1]],
                np.log([radius[hi], radius[hi + 1]]),
            )
        )
    )
    return dict(
        best_radius_m=best,
        profile_lower_radius_m=lower,
        profile_upper_radius_m=upper,
        lower_open=lo == 0,
        upper_open=hi == len(radius) - 1,
        profile_threshold=threshold,
        status="conditional working-likelihood profile; not total uncertainty or a marginal posterior",
    )
