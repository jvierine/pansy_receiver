//! Full-band diagnostic plots; all plotting owned by Rustplotlib.
use anyhow::{Result, ensure};
use hdf5::{File, Group, H5Type};
use rustplotlib::*;
use std::path::Path;
fn read<T: H5Type>(g: &Group, n: &str) -> Result<Vec<T>> {
    Ok(g.dataset(n)?.read_raw()?)
}
fn edges(x: &[f64]) -> Vec<f64> {
    let mut e = vec![x[0] - (x[1] - x[0]) / 2.];
    e.extend(x.windows(2).map(|w| (w[0] + w[1]) / 2.));
    e.push(x[x.len() - 1] + (x[x.len() - 1] - x[x.len() - 2]) / 2.);
    e
}
fn time(h: f64) -> String {
    let s = (h * 3600.).round() as u32;
    format!("{:02}:{:02}:{:02}", s / 3600, (s / 60) % 60, s % 60)
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    ensure!(a.len() == 3, "usage: pansy-wing-context INPUT OUTPUT_DIR");
    let input = File::open(&a[1])?;
    let out = Path::new(&a[2]);
    std::fs::create_dir_all(out)?;
    for name in input
        .member_names()?
        .into_iter()
        .filter(|n| n.starts_with("episode"))
    {
        let g = input.group(&name)?;
        let vel = read::<f64>(&g, "velocity_mps")?;
        let range = read::<f64>(&g, "range_km")?;
        let selected = read::<f64>(&g, "selected_range_km")?[0];
        let r = read::<u32>(&g, "selected_range_index")?[0] as usize;
        let t = read::<f64>(&g, "selected_mode_hours")?[0];
        let start = read::<f64>(&g, "mode_start_hours")?;
        let stop = read::<f64>(&g, "mode_stop_hours")?;
        let window = read::<f64>(&g, "window_hours")?;
        let db = read::<f32>(&g, "time_doppler_db_mode_frequency")?;
        let rd = read::<f32>(&g, "snapshot_db_frequency_range")?;
        let nf = vel.len();
        let nr = range.len();
        let mut f = Figure::new(IntendedUse::ArticleWide);
        f.width_mm = 190.;
        f.height_mm = 278.;
        f.font_pt = 10.;
        f.raster_dpi = 400.;
        f.latex.cache = out.join(".latex-cache");
        f.labels.push(Label::new(
            format!("PANSY 2026-01-30: central-wing inspection, {} UTC", time(t)),
            100.,
            5.,
        ));
        f.scales.push(Scale::new(
            [0., 30.],
            Normalization::Linear,
            ColorMap::plasma(),
        )?);
        let centers: Vec<_> = start.iter().zip(&stop).map(|(a, b)| (a + b) / 2.).collect();
        let mut xe = vec![window[0]];
        let mut columns = Vec::new();
        for m in 0..start.len() {
            let l = if m > 0 && centers[m] - centers[m - 1] <= 1. / 6. {
                (centers[m - 1] + centers[m]) / 2.
            } else {
                start[m].max(window[0])
            };
            let rr = if m + 1 < start.len() && centers[m + 1] - centers[m] <= 1. / 6. {
                (centers[m] + centers[m + 1]) / 2.
            } else {
                stop[m].min(window[1])
            };
            if l > *xe.last().unwrap() + 1e-10 {
                xe.push(l);
                columns.push(None);
            }
            xe.push(rr);
            columns.push(Some(m));
        }
        if *xe.last().unwrap() < window[1] {
            xe.push(window[1]);
            columns.push(None);
        }
        let nx = columns.len();
        let mut pixels = vec![f64::NAN; nf * nx];
        for (x, m) in columns.iter().enumerate() {
            if let Some(m) = m {
                for v in 0..nf {
                    pixels[v * nx + x] = db[m * nf + v] as f64;
                }
            }
        }
        let mut ax = Axes::new(
            [window[0], window[1]],
            [-200., 200.],
            "Time (UTC)",
            r"Radial velocity (m s$^{-1}$)",
        );
        ax.yticks = vec![-200., -100., 0., 100., 200.];
        ax.xticks = vec![window[0], (window[0] + window[1]) / 2., window[1]];
        ax.xticklabels = Some(
            ax.xticks
                .iter()
                .map(|v| time(*v)[..5].to_string())
                .collect(),
        );
        let mut p = Panel::new(Rect::new(25., 18., 143., 44.), ax);
        p.push(Layer::Pcolormesh(Field::new(xe, edges(&vel), pixels, 0)?));
        let mut marker = Line::new(vec![[t, -200.], [t, 200.]], "#ffffff");
        marker.width_pt = 0.6;
        marker.dash = Some("1 1".into());
        p.push(Layer::Line(marker));
        f.panels.push(p);
        f.labels.push(Label::new(
            format!("Fixed slant range {selected:.2} km; dashed line = snapshot"),
            96.5,
            13.,
        ));
        let mut ax = Axes::new(
            [-200., 200.],
            [range[0], *range.last().unwrap()],
            r"Radial velocity (m s$^{-1}$)",
            "Slant range (km)",
        );
        ax.xticks = vec![-200., -100., 0., 100., 200.];
        let mut p = Panel::new(Rect::new(25., 87., 143., 53.), ax);
        let values = (0..nr)
            .flat_map(|r| {
                (0..nf).map({
                    let rd = &rd;
                    move |v| rd[v * nr + r] as f64
                })
            })
            .collect();
        p.push(Layer::Pcolormesh(Field::new(
            edges(&vel),
            edges(&range),
            values,
            0,
        )?));
        let mut line = Line::new(vec![[-200., selected], [200., selected]], "#ffffff");
        line.width_pt = 0.6;
        line.dash = Some("1 1".into());
        p.push(Layer::Line(line));
        f.panels.push(p);
        f.labels.push(Label::new(
            format!("Range--Doppler snapshot, {} UTC", time(t)),
            96.5,
            82.,
        ));
        let mut ax = Axes::new(
            [-200., 200.],
            [0., 40.],
            r"Radial velocity (m s$^{-1}$)",
            "Total/noise (dB)",
        );
        ax.xticks = vec![-200., -100., 0., 100., 200.];
        ax.yticks = vec![0., 10., 20., 30., 40.];
        let mut p = Panel::new(Rect::new(25., 165., 143., 32.), ax);
        p.push(Layer::Line(Line::new(
            vel.iter()
                .enumerate()
                .map(|(v, &x)| [x, rd[v * nr + r] as f64])
                .collect(),
            "#183e7c",
        )));
        f.panels.push(p);
        f.labels.push(Label::new(
            format!("Unmasked spectrum at {selected:.2} km"),
            96.5,
            160.,
        ));
        let fi = read::<u32>(&g, "frequency_index")?;
        let ri = read::<u32>(&g, "range_index")?;
        let theta = read::<f64>(&g, "capon_theta_deg")?;
        let flag = read::<u8>(&g, "decision")?;
        let mut ax = Axes::new(
            [-200., 200.],
            [0., 10.],
            r"Radial velocity (m s$^{-1}$)",
            "Zenith angle (deg)",
        );
        ax.xticks = vec![-200., -100., 0., 100., 200.];
        ax.yticks = vec![0., 2., 5., 10.];
        let mut p = Panel::new(Rect::new(25., 222., 143., 27.), ax);
        let mut line = Line::new(vec![[-200., 2.], [200., 2.]], "#888888");
        line.dash = Some("2 1".into());
        p.push(Layer::Line(line));
        // Separate line segments; do not join across invalid or untested bins.
        for (code, color) in [
            (0, "#16743f"),
            (1, "#ce6b16"),
            (3, "#ce6b16"),
            (4, "#a22e82"),
        ] {
            let mut points = vec![[f64::NAN, f64::NAN]; nf];
            for i in 0..fi.len() {
                if ri[i] as usize == r && flag[i] == code && theta[i].is_finite() {
                    points[fi[i] as usize] = [vel[fi[i] as usize], theta[i]];
                }
            }
            let mut line = Line::new(points, color);
            line.width_pt = 1.5;
            p.push(Layer::Line(line));
        }
        f.panels.push(p);
        f.labels.push(Label::new(
            "Capon: green retained; orange off-zenith/edge; purple code-ghost flag",
            96.5,
            216.,
        ));
        let mut cb = Colorbar::new(
            Rect::new(173., 18., 2., 122.),
            0,
            "Total/noise (dB)",
            vec![0., 10., 20., 30.],
        );
        cb.horizontal = false;
        f.colorbars.push(cb);
        f.labels.push(Label::new(
            "Angle gaps: untested or invalid covariance; central |v| < 15 was not screened",
            96.5,
            274.,
        ));
        for ext in ["png", "pdf"] {
            f.save(out.join(format!("2026-01-30-{name}-central-wing.{ext}")))?;
        }
        println!(
            "PLOT {name} {} UTC range={selected:.3} all_bins={nf}",
            time(t)
        );
    }
    Ok(())
}
