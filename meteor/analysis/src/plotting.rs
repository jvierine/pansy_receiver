//! Shared paper-style HDF5 event/radiant figures, with no generator labels.
use anyhow::Result;
use hdf5::{File, Group};
use rustplotlib::*;
use std::path::Path;

fn data(g: &Group, name: &str) -> Result<Vec<f64>> {
    Ok(g.dataset(name)?.read_raw::<f64>()?)
}
fn scalar(g: &Group, name: &str) -> Result<f64> {
    Ok(g.dataset(name)?.read_scalar::<f64>()?)
}
fn points(x: &[f64], y: &[f64]) -> Vec<[f64; 2]> {
    x.iter().zip(y).map(|(&a, &b)| [a, b]).collect()
}
fn limits(v: &[f64], padding: f64) -> [f64; 2] {
    let min = v
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::INFINITY, f64::min);
    let max = v
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    let d = (max - min).max(1e-4) * padding;
    [min - d, max + d]
}
fn scatter(p: &mut Panel, xy: &[[f64; 2]], color: [f64; 3], size: f64) -> Result<()> {
    p.marker_groups(&[MarkerGroup {
        label: "Measurement",
        points: xy,
        style: MarkerStyle {
            shape: MarkerShape::Circle,
            size_pt: size,
            color,
        },
    }])?;
    Ok(())
}
fn line(p: &mut Panel, xy: Vec<[f64; 2]>, color: &str, dash: bool) {
    let mut l = Line::new(xy, color);
    l.width_pt = 1.2;
    if dash {
        l.dash = Some("2,1".into());
    }
    p.push(Layer::Line(l));
}
fn axes(rect: Rect, x: [f64; 2], y: [f64; 2], xl: &str, yl: &str) -> Panel {
    let a = Axes::new(x, y, xl, yl);
    Panel::new(rect, a)
}

pub fn mst_figure(path: &Path, output: &Path) -> Result<()> {
    let h = File::open(path)?;
    let range = data(&h, "range_km")?;
    let freq = data(&h, "doppler_hz")?;
    let nr = range.len();
    let nf = freq.len();
    let edges = |x: &[f64]| {
        let d = x[1] - x[0];
        (0..=x.len())
            .map(|i| x[0] + (i as f64 - 0.5) * d)
            .collect::<Vec<_>>()
    };
    std::fs::create_dir_all(output)?;
    for overview in [true, false] {
        let rows = if overview { 2 } else { 5 };
        let mut f =
            ArticleLayout::new(179.324, 10.95)?.figure(1., if overview { 160. } else { 335. })?;
        let scale = f.scale(Scale::new(
            [-3., 20.],
            Normalization::Linear,
            ColorMap::sequential(),
        )?);
        for row in 0..rows {
            let rect = Rect::new(27., 12. + row as f64 * 64., 125., 46.);
            let (key, ch, label) = if overview {
                if row == 0 {
                    ("before_db", 2, "Zenith: before rejection")
                } else {
                    ("clean_db", 2, "Zenith: masked; coherent mean removed")
                }
            } else {
                (
                    "clean_db",
                    row,
                    [
                        "South receive",
                        "North receive",
                        "Zenith receive",
                        "East receive",
                        "West receive",
                    ][row],
                )
            };
            let raw = data(&h, key)?;
            let mut values = vec![0.; nr * nf];
            for r in 0..nr {
                for j in 0..nf {
                    values[r * nf + j] = raw[(ch * nr + r) * nf + j];
                }
            }
            let mut p = axes(
                rect,
                [-15., 15.],
                [9., 80.],
                "Doppler frequency (Hz)",
                "Range (km)",
            );
            p.push(Layer::Pcolormesh(Field::new(
                edges(&freq),
                edges(&range),
                values,
                scale,
            )?));
            if !overview || row == 1 {
                let mut band = Line::new(
                    vec![[-0.5, 9.], [0.5, 9.], [0.5, 80.], [-0.5, 80.], [-0.5, 9.]],
                    "#cccccc",
                );
                band.fill = Some("#cccccc".into());
                band.fill_opacity = 1.;
                band.width_pt = 0.;
                p.push(Layer::Line(band));
            }
            p.labels
                .push(Label::new(label, rect.x + rect.width / 2., rect.y - 4.));
            f.panels.push(p);
            f.colorbars.push(Colorbar::new(
                Rect::new(155., rect.y, 2., rect.height),
                scale,
                "Power/noise (dB)",
                vec![0., 10., 20.],
            ));
        }
        let name = if overview {
            "mst_zenith_comparison"
        } else {
            "mst_receive_beams"
        };
        f.save_for_article(output.join(format!("{name}.pdf")), &format!("{name}.pdf"))?;
        f.save(output.join(format!("{name}.png")))?;
    }
    Ok(())
}

pub fn event_figure(path: &Path, output: &Path, atmosphere: &Path) -> Result<()> {
    let h = File::open(path)?;
    let event = path.file_stem().unwrap().to_string_lossy();
    let t = data(&h, "measurement/time_s")?;
    let n = t.len();
    let tx = limits(&t, 0.03);
    let position = data(&h, "model/position_m")?;
    let velocity = data(&h, "model/velocity_mps")?;
    let uv = data(&h, "original_conditional_uv")?;
    let ranges = data(&h, "measurement/range_m")?;
    let obs_v = data(&h, "measurement/radial_velocity_mps")?
        .iter()
        .map(|v| v / 1000.)
        .collect::<Vec<_>>();
    let pred_v = data(&h, "model/radial_velocity_mps")?
        .iter()
        .map(|v| v / 1000.)
        .collect::<Vec<_>>();
    let covariance = data(&h, "measurement/range_velocity_covariance")?;
    let x = (0..n)
        .map(|i| uv[2 * i] * ranges[i] / 1000.)
        .collect::<Vec<_>>();
    let y = (0..n)
        .map(|i| uv[2 * i + 1] * ranges[i] / 1000.)
        .collect::<Vec<_>>();
    let px = position
        .as_chunks::<3>()
        .0
        .iter()
        .map(|p| p[0] / 1000.)
        .collect::<Vec<_>>();
    let py = position
        .as_chunks::<3>()
        .0
        .iter()
        .map(|p| p[1] / 1000.)
        .collect::<Vec<_>>();
    let station_alt = h.attr("altitude_m")?.read_scalar::<f64>()?;
    let observed_height = (0..n)
        .map(|i| {
            station_alt / 1000.
                + (1. - uv[2 * i].powi(2) - uv[2 * i + 1].powi(2))
                    .max(0.)
                    .sqrt()
                    * ranges[i]
                    / 1000.
        })
        .collect::<Vec<_>>();
    let height = position
        .as_chunks::<3>()
        .0
        .iter()
        .map(|p| (station_alt + p[2] + (p[0] * p[0] + p[1] * p[1]) / (2. * 6371000.)) / 1000.)
        .collect::<Vec<_>>();
    let mut f = ArticleLayout::new(179.324, 10.95)?.figure(1., 237.)?;
    f.labels.push(Label::new(
        format!(
            "{}: provisional array-frame trajectory",
            event.replace('_', " ")
        ),
        90.,
        3.,
    ));
    let rects = [
        Rect::new(21.5, 13., 48., 48.),
        Rect::new(112., 13., 61., 48.),
        Rect::new(23., 84., 61., 48.),
        Rect::new(112., 84., 46., 48.),
        Rect::new(23., 155., 61., 48.),
        Rect::new(112., 155., 61., 48.),
    ];
    let extent = x
        .iter()
        .chain(&y)
        .chain(&px)
        .chain(&py)
        .map(|v| v.abs())
        .fold(0., f64::max)
        .max(3.);
    let mut p = axes(
        rects[0],
        [-extent * 1.1, extent * 1.1],
        [-extent * 1.1, extent * 1.1],
        "Array x (km)",
        "Array y (km)",
    );
    let beam = ranges[0] / 1000. * 1.39f64.to_radians().tan();
    line(
        &mut p,
        (0..=120)
            .map(|i| {
                let a = i as f64 / 120. * std::f64::consts::TAU;
                [beam * a.cos(), beam * a.sin()]
            })
            .collect(),
        "#aaaaaa",
        true,
    );
    let snr = data(&h, "measurement/matched_filter_snr_db")?;
    let snr_scale = Scale::new([0., 30.], Normalization::Linear, ColorMap::sequential())?;
    p.scatter_colored(&points(&x, &y), &snr, &snr_scale, 2.8, 1., "#444444")?;
    line(&mut p, points(&px, &py), "#0072b2", false);
    p.labels.push(Label::new(
        format!(
            "RMS $x={:.0}$ m, $y={:.0}$ m",
            scalar(&h, "result/horizontal_x_rms_m")?,
            scalar(&h, "result/horizontal_y_rms_m")?
        ),
        rects[0].x + rects[0].width / 2.,
        rects[0].y - 4.,
    ));
    f.panels.push(p);
    let snr_scale_id = f.scale(snr_scale);
    f.colorbars.push(Colorbar::new(
        Rect::new(73., 13., 2., 48.),
        snr_scale_id,
        "SNR (dB)",
        vec![0., 15., 30.],
    ));
    let heights = observed_height
        .iter()
        .chain(&height)
        .copied()
        .collect::<Vec<_>>();
    let mut p = axes(
        rects[1],
        tx,
        limits(&heights, 0.08),
        "Time (s)",
        "Height MSL (km)",
    );
    scatter(
        &mut p,
        &points(&t, &observed_height),
        [0.35, 0.35, 0.35],
        2.,
    )?;
    line(&mut p, points(&t, &height), "#0072b2", false);
    p.labels.push(Label::new(
        format!("RMS {:.0} m", scalar(&h, "result/height_rms_m")?),
        rects[1].x + rects[1].width / 2.,
        rects[1].y - 4.,
    ));
    f.panels.push(p);
    let pair_t = data(&h, "pulse_phase/time_s")?;
    let pair_v = data(&h, "pulse_phase/velocity_mps")?
        .iter()
        .map(|v| v / 1000.)
        .collect::<Vec<_>>();
    let pair_sigma = data(&h, "pulse_phase/sigma_mps")?
        .iter()
        .map(|v| v / 1000.)
        .collect::<Vec<_>>();
    let keep = h.dataset("pulse_phase/keep")?.read_raw::<u8>()?;
    let all_v = obs_v
        .iter()
        .chain(&pred_v)
        .chain(&pair_v)
        .copied()
        .collect::<Vec<_>>();
    let mut p = axes(
        rects[2],
        tx,
        limits(&all_v, 0.1),
        "Time (s)",
        "Radial velocity (km/s)",
    );
    let bounds = (0..n)
        .map(|i| {
            let s = covariance[4 * i + 3].sqrt() / 1000.;
            [obs_v[i] - s, obs_v[i] + s]
        })
        .collect::<Vec<_>>();
    p.errorbars(
        &points(&t, &obs_v),
        None,
        Some(&bounds),
        "#999999",
        0.5,
        0.,
        0.6,
    )?;
    scatter(&mut p, &points(&t, &obs_v), [0.15, 0.15, 0.15], 2.)?;
    for accepted in [false, true] {
        let indices = (0..pair_t.len())
            .filter(|&i| (keep[i] != 0) == accepted)
            .collect::<Vec<_>>();
        let xy = indices
            .iter()
            .map(|&i| [pair_t[i], pair_v[i]])
            .collect::<Vec<_>>();
        let bounds = indices
            .iter()
            .map(|&i| [pair_v[i] - pair_sigma[i], pair_v[i] + pair_sigma[i]])
            .collect::<Vec<_>>();
        p.errorbars(
            &xy,
            None,
            Some(&bounds),
            if accepted { "#009e73" } else { "#cc6666" },
            0.5,
            0.,
            0.7,
        )?;
        scatter(
            &mut p,
            &xy,
            if accepted {
                [0., 0.55, 0.4]
            } else {
                [0.8, 0.3, 0.3]
            },
            3.,
        )?;
    }
    line(&mut p, points(&t, &pred_v), "#0072b2", false);
    let phase_rms = scalar(&h, "result/pulse_phase_velocity_rms_mps")?;
    let phase_text = if phase_rms.is_finite() {
        format!("{phase_rms:.0}")
    } else {
        "n/a".into()
    };
    p.labels.push(Label::new(
        format!(
            "RMS single {:.0}; phase {} m/s",
            scalar(&h, "result/radial_velocity_rms_mps")?,
            phase_text
        ),
        rects[2].x + rects[2].width / 2.,
        rects[2].y - 4.,
    ));
    f.panels.push(p);
    let te = data(&h, "plotting/time_edges_s")?;
    let re = data(&h, "plotting/range_edges_km")?;
    let scale = f.scale(Scale::new(
        [0., 25.],
        Normalization::Linear,
        ColorMap::sequential(),
    )?);
    let mut p = axes(
        rects[3],
        [te[0], *te.last().unwrap()],
        [re[0], *re.last().unwrap()],
        "Time (s)",
        "Range (km)",
    );
    p.push(Layer::Pcolormesh(Field::new(
        te,
        re,
        data(&h, "plotting/power_over_noise_db")?,
        scale,
    )?));
    let model_r = position
        .as_chunks::<3>()
        .0
        .iter()
        .map(|p| p.iter().map(|v| v * v).sum::<f64>().sqrt() / 1000.)
        .collect::<Vec<_>>();
    line(&mut p, points(&t, &model_r), "#ffffff", false);
    p.labels.push(Label::new(
        format!("RMS {:.0} m", scalar(&h, "result/range_rms_m")?),
        rects[3].x + rects[3].width / 2.,
        rects[3].y - 4.,
    ));
    f.panels.push(p);
    f.colorbars.push(Colorbar::new(
        Rect::new(160., 84., 2., 48.),
        scale,
        "Power/noise (dB)",
        vec![0., 10., 20.],
    ));
    let radius = data(&h, "profile_with_pulse_phase/radius_m")?
        .iter()
        .map(|v| v * 1e6)
        .collect::<Vec<_>>();
    let mut p = axes(
        rects[4],
        [1., 100000.],
        [0., 12.],
        r"Radius $r_0$ ($\mu$m)",
        r"Profile $\Delta\chi^2$",
    );
    p.axes.as_mut().unwrap().set_xscale(AxisScale::Log10)?;
    for (name, color) in [
        ("profile_single_pulse/chi2", "#888888"),
        ("profile_with_pulse_phase/chi2", "#0072b2"),
    ] {
        let cost = data(&h, name)?;
        let best = cost.iter().copied().fold(f64::INFINITY, f64::min);
        line(
            &mut p,
            // Keep off-panel profile values bounded before rasterization.
            points(
                &radius,
                &cost.iter().map(|v| (v - best).min(24.)).collect::<Vec<_>>(),
            ),
            color,
            false,
        );
    }
    line(&mut p, vec![[1., 3.84], [100000., 3.84]], "#444444", true);
    f.panels.push(p);
    let speed = velocity
        .as_chunks::<3>()
        .0
        .iter()
        .map(|v| v.iter().map(|x| x * x).sum::<f64>().sqrt() / 1000.)
        .collect::<Vec<_>>();
    let mut speed_curves = Vec::new();
    let atm = File::open(atmosphere)?;
    let alt = data(&atm, "altitude_m")?;
    let rho = data(&atm, "density_kg_m3")?;
    let parameters = data(&h, "profile_with_pulse_phase/parameters")?;
    for (target, color) in [(10., "#cc79a7"), (100., "#e69f00"), (1000., "#d55e00")] {
        let k = (0..radius.len())
            .min_by(|&i, &j| {
                (radius[i] / target)
                    .ln()
                    .abs()
                    .total_cmp(&(radius[j] / target).ln().abs())
            })
            .unwrap();
        let mut initial = [0.; 7];
        for j in 0..6 {
            initial[j] = parameters[k * 15 + j] * 1000.;
        }
        initial[6] = radius[k] * 1e-6;
        if let Some(states) = crate::propagate(
            initial,
            &t,
            &alt,
            &rho,
            station_alt,
            3000.,
            1e-8,
            0.0005,
            true,
        ) {
            let s = states
                .iter()
                .map(|v| v[3..6].iter().map(|x| x * x).sum::<f64>().sqrt() / 1000.)
                .collect::<Vec<_>>();
            speed_curves.push((s, color));
        }
    }
    let all_speeds = speed
        .iter()
        .copied()
        .chain(speed_curves.iter().flat_map(|(v, _)| v.iter().copied()))
        .collect::<Vec<_>>();
    let mut p = axes(
        rects[5],
        tx,
        limits(&all_speeds, 0.08),
        "Time (s)",
        "Speed (km/s)",
    );
    for (s, color) in speed_curves {
        line(&mut p, points(&t, &s), color, true);
    }
    line(&mut p, points(&t, &speed), "#0072b2", false);
    f.panels.push(p);
    f.labels.push(Label::new(
        "Grey: single pulse; green: pulse phase; red: ambiguous/rejected",
        90.,
        225.,
    ));
    f.labels.push(Label::new(
        r"Blue: fit; dashed: $r_0=10,100,1000\,\mu$m (magenta, amber, red)",
        90.,
        231.,
    ));
    std::fs::create_dir_all(output)?;
    f.save_for_article(output.join(format!("{event}.pdf")), &format!("{event}.pdf"))?;
    f.save(output.join(format!("{event}.png")))?;
    phase_figure(&h, &event, output)?;
    Ok(())
}

fn phase_figure(h: &File, event: &str, output: &Path) -> Result<()> {
    let t = data(h, "measurement/time_s")?;
    let measured = data(h, "plotting/channel_phase_rad")?;
    let model = data(h, "plotting/model_channel_phase_rad")?;
    let mut f = ArticleLayout::new(179.324, 10.95)?.figure(1., 302.)?;
    let names = [
        "S0",
        "S1",
        "S2",
        "L0 (reference)",
        "L1",
        "L2",
        "L3",
        "L4",
        "L5",
        "L6",
    ];
    for (ch, name) in names.iter().enumerate() {
        let rect = Rect::new(
            23. + (ch % 2) as f64 * 89.,
            8. + (ch / 2) as f64 * 59.,
            61.,
            39.,
        );
        let mut p = axes(
            rect,
            limits(&t, 0.03),
            [-3.3, 3.3],
            "Time (s)",
            "Phase (rad)",
        );
        let y = (0..t.len())
            .map(|i| measured[i * 10 + ch])
            .collect::<Vec<_>>();
        let m = (0..t.len()).map(|i| model[i * 10 + ch]).collect::<Vec<_>>();
        scatter(&mut p, &points(&t, &y), [0.35, 0.35, 0.35], 2.)?;
        // Leave phase wraps disconnected rather than drawing through the panel.
        let mut curve = Vec::new();
        for i in 0..t.len() {
            if i > 0 && (m[i] - m[i - 1]).abs() > 3. {
                curve.push([f64::NAN; 2]);
            }
            curve.push([t[i], m[i]]);
        }
        line(&mut p, curve, "#0072b2", false);
        p.labels
            .push(Label::new(*name, rect.x + rect.width / 2., rect.y + 4.));
        f.panels.push(p);
    }
    f.save_for_article(
        output.join(format!("{event}_channels.pdf")),
        &format!("{event}_channels.pdf"),
    )?;
    f.save(output.join(format!("{event}_channels.png")))?;
    Ok(())
}

fn hammer(lon: f64, lat: f64) -> [f64; 2] {
    let a = -((lon - 270. + 180.).rem_euclid(360.) - 180.).to_radians();
    let b = lat.to_radians();
    let den = (1. + b.cos() * (a / 2.).cos()).sqrt();
    [
        2. * 2f64.sqrt() * b.cos() * (a / 2.).sin() / den,
        2f64.sqrt() * b.sin() / den,
    ]
}

pub fn radiant_figure(run: &Path, output: &Path) -> Result<()> {
    let cat = File::open(run.join("catalogue.h5"))?;
    let mut f = ArticleLayout::new(179.324, 10.95)?.figure(1., 185.)?;
    let mut p = axes(
        Rect::new(30., 8., 135.3, 70.),
        [-2.9, 2.9],
        [-1.5, 1.5],
        "Sun-centred ecliptic longitude",
        "Ecliptic latitude",
    );
    p.axes.as_mut().unwrap().xticks = vec![-2. * 2f64.sqrt(), -2., 0., 2., 2. * 2f64.sqrt()];
    p.axes.as_mut().unwrap().xticklabels = Some(
        ["90", "0", "270", "180", "90"]
            .iter()
            .map(|s| format!("${s}^\\circ$"))
            .collect(),
    );
    p.axes.as_mut().unwrap().yticks = [-60f64, -30., 0., 30., 60.]
        .iter()
        .map(|v| 2f64.sqrt() * v.to_radians().sin())
        .collect();
    p.axes.as_mut().unwrap().yticklabels = Some(
        ["-60", "-30", "0", "30", "60"]
            .iter()
            .map(|s| format!("${s}^\\circ$"))
            .collect(),
    );
    line(
        &mut p,
        (0..=360)
            .map(|i| {
                let a = i as f64 * std::f64::consts::TAU / 360.;
                [2. * 2f64.sqrt() * a.cos(), 2f64.sqrt() * a.sin()]
            })
            .collect(),
        "#aaaaaa",
        false,
    );
    for lat in [-60., -30., 0., 30., 60.] {
        line(
            &mut p,
            (0..=360)
                .map(|i| hammer(90.00001 + i as f64 * (359.99998 / 360.), lat))
                .collect(),
            "#dddddd",
            false,
        );
    }
    for lon in [0., 90., 180., 270.] {
        line(
            &mut p,
            (-90..=90).map(|lat| hammer(lon, lat as f64)).collect(),
            "#dddddd",
            false,
        );
    }
    let mut q = axes(
        Rect::new(23., 115., 151., 48.),
        [0., 360.],
        [0., 90.],
        "Array azimuth (deg)",
        "Radiant zenith angle (deg)",
    );
    for name in cat.group("events")?.member_names()?.iter() {
        let g = cat.group(&format!("events/{name}"))?;
        if !g.link_exists("geocentric_first_detected_conditional") {
            continue;
        }
        let geo = g.group("preferred_geocentric_conditional")?;
        let lon = scalar(&geo, "sun_centered_ecliptic_longitude_deg")?;
        let lat = scalar(&geo, "ecliptic_latitude_deg")?;
        let speed = scalar(&geo, "geocentric_speed_km_s")?;
        let color = if speed < 30. {
            [0.8, 0.35, 0.]
        } else if speed < 50. {
            [0., 0.5, 0.4]
        } else {
            [0., 0.35, 0.7]
        };
        scatter(&mut p, &[hammer(lon, lat)], color, 5.)?;
        let az = scalar(&g, "array_radiant_azimuth_deg")?;
        let zen = scalar(&g, "array_radiant_zenith_angle_deg")?;
        scatter(&mut q, &[[az, zen]], color, 5.)?;
    }
    let mut legend = Legend::new(vec![
        LegendEntry::circle("below 30 km/s", "#cc5900"),
        LegendEntry::circle("30--50 km/s", "#008066"),
        LegendEntry::circle("above 50 km/s", "#0059b3"),
    ]);
    legend.allow_reused_marker_shapes = true;
    legend.position = LegendPosition::LowerLeft;
    q.legend = Some(legend);
    f.panels.push(p);
    f.panels.push(q);
    f.labels.push(Label::new(
        "Illustrative celestial map: x east, y north assumed",
        98.,
        97.,
    ));
    f.save_for_article(output.join("radiants.pdf"), "radiants.pdf")?;
    f.save(output.join("radiants.png"))?;
    Ok(())
}

pub fn plot_run(run: &Path, output: &Path) -> Result<()> {
    std::fs::create_dir_all(output)?;
    let mut paths = std::fs::read_dir(run.join("events"))?
        .filter_map(|p| p.ok())
        .map(|p| p.path())
        .filter(|p| p.extension().is_some_and(|v| v == "h5"))
        .collect::<Vec<_>>();
    paths.sort();
    for path in paths {
        event_figure(&path, output, &run.join("atmosphere.h5"))?;
        println!("Plotted {}", path.file_stem().unwrap().to_string_lossy());
    }
    radiant_figure(run, output)
}
