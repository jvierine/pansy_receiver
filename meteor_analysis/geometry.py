"""Element-sum response lookup, including the explicitly provisional active mask."""

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from .native import array_response


class ArrayModel:
    def __init__(self, path, wavelength, n=161, extent_deg=6.0):
        with h5py.File(path) as h:
            xy = h["xy_m"][()]
            self.membership = h["membership"][()]
            self.active = h["provisional_active_mask"][()].astype(np.uint8)
        self.xy = xy - (xy.min(axis=0) + xy.max(axis=0)) / 2
        self.wavelength = wavelength
        self.axis = np.linspace(
            -np.sin(np.radians(extent_deg)), np.sin(np.radians(extent_deg)), n
        )
        u, v = np.meshgrid(self.axis, self.axis, indexing="ij")
        response = array_response(
            self.xy,
            self.membership,
            self.active,
            np.column_stack([u.ravel(), v.ravel()]),
            wavelength,
        ).reshape(n, n, 10)
        self.interpolator = RegularGridInterpolator(
            (self.axis, self.axis),
            response,
            bounds_error=False,
            fill_value=complex(np.nan),
        )

    def __call__(self, uv):
        return self.interpolator(uv)

    def check(self, uv):
        exact = array_response(
            self.xy, self.membership, self.active, np.asarray(uv), self.wavelength
        )
        approximation = self(uv)
        use = abs(exact) > 0.05 * np.sum(self.membership * self.active[:, None], axis=0)
        return float(np.max(abs(np.angle(approximation[use] * exact[use].conj()))))
