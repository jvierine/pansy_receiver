"""Geodetic radar sites; array-axis orientation is a separate, explicit assumption."""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Station:
    name: str
    latitude_deg: float
    longitude_deg: float
    altitude_m: float
    coordinate_status: str = "surveyed"

    def enu_to_ecef(self):
        lat, lon = np.radians([self.latitude_deg, self.longitude_deg])
        return np.array(
            [
                [-np.sin(lon), -np.sin(lat) * np.cos(lon), np.cos(lat) * np.cos(lon)],
                [np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat) * np.sin(lon)],
                [0.0, np.cos(lat), np.sin(lat)],
            ]
        )


TIBET_MST = Station(
    "Tibet MST",
    30.1,
    90.5,
    4287.0,
    "approximate lat/lon from instrument presentation; altitude from antenna spreadsheet",
)
PANSY = Station(
    "PANSY",
    -69.010833,
    39.599722,
    0.0,
    "matches existing pansy_config and pansy_ballistic convention",
)
