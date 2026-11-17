"""Compile a concise scientific overview and append the native-size event figures."""

import argparse
from pathlib import Path
import subprocess
import h5py
import numpy as np
from pypdf import PdfWriter


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("--include-channel-diagnostics", action="store_true",
                   help="Append the detailed receiver phase/amplitude pages")
    a = p.parse_args()
    directory = a.run / "report"
    directory.mkdir(exist_ok=True)
    rows = []
    with h5py.File(a.run / "catalogue.h5") as h:
        for name, g in h["events"].items():
            number = int(name.split("_")[-1]) + 1
            speed = float(
                g["preferred_geocentric_conditional/geocentric_speed_km_s"][()]
            )
            lo = float(g["profile_lower_radius_m"][()]) * 1e6
            hi = float(g["profile_upper_radius_m"][()]) * 1e6
            best = (
                "--"
                if bool(g["upper_open"][()])
                else f'{float(g["best_radius_m"][()])*1e6:.1f}'
            )
            bounds = f"{lo:.1f}--{hi:.1f}" if np.isfinite(hi) else f"$>{lo:.1f}$"
            rows.append(
                f'{number} & {speed:.1f} & {best} & {bounds} & {int(g["selected_ppp_pairs"][()])} & {int(g["selected_same_code_two_ipp_pairs"][()])} \\\\'
            )
    tex = (
        r"""\documentclass[11pt]{article}
% Scientific inputs and generator provenance are retained in HDF5 and source comments.
\usepackage[letterpaper,margin=0.72in]{geometry}
\usepackage[T1]{fontenc}\usepackage{newtxtext,newtxmath,booktabs,hyperref}
\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}
\begin{document}
\begin{center}\Large Tibet MST meteor trajectories\\[4pt]
\large Provisional radiant and shrinking-radius analysis\end{center}
\textbf{Observations.} Tibet MST, $30.1^\circ$ N, $90.5^\circ$ E, 4287 m;
50 MHz. Five files recorded on 14 November 2026, 05:51:56--05:52:35.864
Beijing time, interpreted as 13 November, 21:51:56--21:52:35.864 UTC.
Filename timezone and absolute subsecond timing remain unconfirmed.

\textbf{Important limitations.} The antenna axes and receive-phase sign are
not confirmed. The celestial radiant map assumes array $x$ east and $y$ north;
it is illustrative, not an absolute sky calibration. The lower map is in
array coordinates. All fits assume the nominal unaliased range. Even-IPP
range ambiguity remains unresolved. Receiver phases use the provisional
boresight-centred meteor ensemble calibration.

\textbf{Atmosphere and dynamics.} NRLMSIS 2.1 over Tibet at the observation
date and time: previous-day F10.7 = 115.5 sfu, centred 81-day mean = 114.2 sfu,
daily Ap = 3. The centred mean can include forecast values. A compact sphere
with bulk density $3000\ \mathrm{kg\,m^{-3}}$, free-molecular drag coefficient
$\Gamma=1$, and apparent ablation coefficient $\sigma=10^{-8}\ \mathrm{kg/J}$
is fitted jointly to range, within-pulse Doppler, receiver phases and retained
nonoverlapping pulse-phase velocity pairs. Wind and Coriolis terms are omitted.
Initial radius refers to the first retained observation, not the pre-atmospheric body.
Geocentric speed uses the backward-extrapolated 120 km state and the incoming
Earth-hyperbola asymptote, including the station's rotation through the inertial
coordinate transformation. All astronomical values remain conditional.

\begin{center}\begin{tabular}{rrrrrr}\toprule
Event & $v_g$ (km/s) & Best $r_0$ ($\mu$m) & Profile range ($\mu$m) & 1-IPP pairs & 2-IPP pairs\\\midrule
"""
        + "\n".join(rows)
        + r"""
\bottomrule\end{tabular}\end{center}
Profile ranges use $\Delta\chi^2=3.84$ in a working compressed-observable
likelihood; they are not marginal posteriors, calibrated-coverage intervals,
or total physical uncertainties. An open upper limit is not replaced by the
numerical grid endpoint. In this model $r_0$ scales directly with atmospheric
density, and inferred mass scales with its cube. Same-code 2-IPP pairs are
diagnostics, not additional independent likelihood terms.

\small Model references: Ceplecha et al. (1998),
\href{https://doi.org/10.1023/A:1005069928850}{Meteor Phenomena and Bodies};
Emmert et al. (2022), \href{https://doi.org/10.1029/2022JA030896}{NRLMSIS 2.1}.
Space weather and Earth orientation: \href{https://celestrak.org/SpaceData/}{CelesTrak}.
\end{document}
"""
    )
    path = directory / "overview.tex"
    path.write_text(tex)
    subprocess.run(
        [
            "latexmk",
            "-pdf",
            "-interaction=nonstopmode",
            "-halt-on-error",
            "overview.tex",
        ],
        cwd=directory,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    atlas = PdfWriter()
    atlas.append(directory / "overview.pdf")
    atlas.append(a.run / "figures/radiants.pdf")
    for path in sorted((a.run / "events").glob("event_*.h5")):
        name = path.stem
        atlas.append(a.run / "figures" / f"{name}.pdf")
        if a.include_channel_diagnostics:
            atlas.append(a.run / "figures" / f"{name}_channels.pdf")
    atlas.write(a.run / "tibet_meteor_event_atlas.pdf")
    print(a.run / "tibet_meteor_event_atlas.pdf")


if __name__ == "__main__":
    main()
