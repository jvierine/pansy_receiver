# PANSY Zenodo staging data

This directory stages data products intended for the PANSY catalogue Zenodo
deposit. The shower-selection CSV files are the authoritative Radview exports
used to calculate the shower-averaged radiants and orbital elements in the
article.

| File | Shower | IAU number | Selected meteors |
|---|---|---:|---:|
| `-1_JUE_pansy01.csv` | July Eridanids (JUE), candidate | -1 | 83 |
| `001_CAP_pansy01.csv` | alpha Capricornids (CAP) | 1 | 283 |
| `115_DCS_pansy01.csv` | Daytime Capricornids-Sagittariids (DCS) | 115 | 249 |
| `1208_OES_pansy01.csv` | 62-Sagittariids (OES) | 1208 | 186 |

The `MetCod` column contains the PANSY event sample identifier used to join
each selection to the orbit metadata. DCS and OES share one boundary event
(`1769396179144494`); it is intentionally retained in both exported
selections.

## Orbit-metadata validation

Validated on 2026-07-27 against the selected-event records in
`/mnt/data/juha/pansy/metadata/orbit` on `revontuli.uit.no`.

- All 801 selection rows matched an orbit record by exact `sample_idx`.
- Every identifier is unique within its selection, and each identifier resolves
  to exactly one selected-event record in one metadata file. There are 800
  unique events across the four files because of the intentional DCS/OES
  overlap noted above.
- The CSV columns are rounded or binned export values. The table below gives
  the maximum absolute difference from the corresponding full-precision orbit
  metadata. Longitude differences use the shortest wrapped angle.

| Shower | Solar longitude (deg) | Sun-centered ecliptic longitude (deg) | Ecliptic latitude (deg) | Geocentric speed (km/s) |
|---|---:|---:|---:|---:|
| JUE | 0.0304 | 0.0397 | 0.0201 | 0.0191 |
| CAP | 0.0589 | 0.1552 | 0.0085 | 0.0117 |
| DCS | 0.1250 | 0.2494 | 0.0085 | 0.0121 |
| OES | 0.1250 | 0.2225 | 0.0084 | 0.0115 |

These bounded residuals are consistent with the coordinate and histogram-bin
precision of the exports. Full-precision orbit metadata, rather than the
rounded CSV display values, is used when calculating the article's shower
means.

## Activity-profile selection windows

The Figure 10 and Figure 11 activity profiles use the joint seven-parameter
windows below. Each bound is the full selected-event envelope in the current
Revontuli orbit metadata, rounded outward; consequently every CSV member is
retained. Right ascension is selected by the listed circular center and
half-width.

| Shower | \(v_g\) (km/s) | \(e\) | \(a\) (AU) | \(i\) (deg) | \(q\) (AU) | RA center +/- half-width (deg) | Dec (deg) |
|---|---|---|---|---|---|---|---|
| JUE | 45.5--52.5 | 0.625--0.830 | 2.2--4.8 | 84.0--99.0 | 0.79--0.96 | 54 +/- 9 | -34.5 to -25.5 |
| CAP | 19.5--27.0 | 0.690--0.830 | 1.6--3.0 | 4.0--12.0 | 0.40--0.69 | 301 +/- 17 | -15.5 to -5.0 |
| DCS | 22.0--28.5 | 0.725--0.845 | 1.7--3.0 | 4.0--10.5 | 0.40--0.56 | 312 +/- 11 | -30.5 to -22.0 |
| OES | 18.5--25.5 | 0.660--0.825 | 1.6--3.6 | 5.0--10.5 | 0.50--0.68 | 301 +/- 14 | -36.5 to -27.5 |

## SHA-256 checksums

```text
4436153b2f3f2018bc30a870e9a74f91b493fc036352bd795eb4527cac642a15  -1_JUE_pansy01.csv
a59e49f9f4c9a654da44d8828d33672d1be2a6dea09d95105954cf41c4225c1f  001_CAP_pansy01.csv
bd77b5844aab35ec1b9f53b700a625dd8ff31b36ea85667297eec0059879f50e  115_DCS_pansy01.csv
ffcf4ec9dbef9a60aba53b9814b74c36969368fe069ca7d44bf22e645137e7c5  1208_OES_pansy01.csv
```
