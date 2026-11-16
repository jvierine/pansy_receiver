"""Explicit MSIS inputs and HDF5 caching; downloaded tabular data stays in memory."""

import csv
import datetime as dt
import io
import urllib.request
import hashlib
import numpy as np
import h5py
import pymsis

URL = "https://celestrak.org/SpaceData/SW-All.csv"


def build_profile(station, epoch, path):
    date = dt.datetime.fromtimestamp(epoch, dt.timezone.utc)
    payload = urllib.request.urlopen(URL, timeout=60).read()
    rows = {r["DATE"]: r for r in csv.DictReader(io.StringIO(payload.decode()))}
    today = rows[date.date().isoformat()]
    previous = rows[(date.date() - dt.timedelta(days=1)).isoformat()]
    f107 = float(previous["F10.7_OBS"])
    f107a = float(today["F10.7_OBS_CENTER81"])
    ap = float(today["AP_AVG"])
    altitude = np.r_[np.arange(50.0, 200.01, 0.25), np.arange(205.0, 600.01, 5.0)]
    time = np.datetime64(date.replace(tzinfo=None), "us")
    values = pymsis.calculate(
        [time],
        [station.longitude_deg],
        [station.latitude_deg],
        altitude,
        f107s=[f107],
        f107as=[f107a],
        aps=np.full((1, 7), ap),
        version=2.1,
        geomagnetic_activity=1,
    )
    rho = values[0, 0, 0, :, 0].astype(float)
    if np.any(~np.isfinite(rho)) or np.any(rho <= 0):
        raise ValueError("invalid MSIS densities")
    with h5py.File(path, "w") as h:
        h["altitude_m"] = altitude * 1000
        h["density_kg_m3"] = rho
        h["species_and_temperature"] = values[0, 0, 0]
        h.attrs.update(
            model="NRLMSIS 2.1",
            utc=date.isoformat(),
            latitude_deg=station.latitude_deg,
            longitude_deg=station.longitude_deg,
            f107_previous_day_sfu=f107,
            f107a_center81_sfu=f107a,
            daily_ap=ap,
            geomagnetic_activity=1,
            space_weather_source=URL,
            source_sha256=hashlib.sha256(payload).hexdigest(),
            pymsis_version=pymsis.__version__,
            weather_status="daily F10.7 quality retained below; centred 81-day average may include forecast values",
        )
        for label, row in [("observation_day", today), ("previous_day", previous)]:
            g = h.create_group("space_weather/" + label)
            for k, v in row.items():
                g.attrs[k] = v
    return altitude * 1000, rho
