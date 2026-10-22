//! Native Rustplotlib figures from the small, binned HDF5 projection.
//! Quantities: peak total/noise, not noise-subtracted signal/noise. Missing
//! observations stay blank; zero or low velocities were not selected.
use anyhow::Result;
use hdf5::{File, Group, H5Type};
use rustplotlib::*;
use std::{fs, path::Path};
fn read<T: H5Type>(g: &Group, n: &str) -> Result<Vec<T>> {
    Ok(g.dataset(n)?.read_raw()?)
}
fn figure(out: &Path, height: f64) -> Figure {
    let mut f = Figure::new(IntendedUse::ArticleWide);
    f.width_mm = 190.;
    f.height_mm = height;
    f.font_pt = 10.;
    f.raster_dpi = 180.;
    f.latex.cache = out.join(".latex-cache");
    f.background = [1., 1., 1.];
    f
}
fn save(f: &Figure, out: &Path, name: &str) -> Result<()> {
    for ext in ["png", "pdf"] {
        f.save(out.join(format!("{name}.{ext}")))?;
    }
    Ok(())
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    anyhow::ensure!(
        args.len() == 3,
        "usage: pansy-legacy-event-plots SMALL_INPUT.h5 OUTPUT_DIR"
    );
    let input = File::open(&args[1])?;
    let out = Path::new(&args[2]);
    fs::create_dir_all(out)?;
    let date = String::from_utf8(read(&input, "date_utc")?)?;
    let o = input.group("overview")?;
    let t = read::<f64>(&o, "time_hours")?;
    let v = read::<f64>(&o, "velocity_mps")?;
    let s = read::<f64>(&o, "peak_total_noise_db")?;
    let d = read::<u8>(&o, "decision")?;
    let mut f = figure(out, 187.);
    f.labels.push(Label::new(
        format!("PANSY zenith events: {date} UTC"),
        98.,
        5.,
    ));
    let scale = Scale::new([5., 30.], Normalization::Linear, ColorMap::plasma())?;
    f.scales.push(scale.clone());
    for panel in 0..3 {
        let rect = Rect::new(24., 19. + panel as f64 * 53., 142., 38.);
        let ylim = if panel == 2 {
            [-60., 60.]
        } else {
            [-200., 200.]
        };
        let mut axes = Axes::new(
            [0., 24.],
            ylim,
            if panel == 2 { "Time (UTC)" } else { "" },
            r"Radial velocity (m s$^{-1}$)",
        );
        axes.xticks = vec![0., 4., 8., 12., 16., 20., 24.];
        axes.xticklabels = Some(
            ["00", "04", "08", "12", "16", "20", "24"]
                .iter()
                .map(|x| x.to_string())
                .collect(),
        );
        axes.yticks = if panel == 2 {
            vec![-60., -30., 0., 30., 60.]
        } else {
            vec![-200., -100., 0., 100., 200.]
        };
        axes.show_x_axis = panel == 2;
        let mut p = Panel::new(rect, axes);
        // One point per integration, 0.5 m/s bin and decision class; peak
        // over selected ranges. No interpolation over unobserved modes.
        let ids: Vec<_> = (0..t.len())
            .filter(|&i| (panel == 0 || d[i] == 0) && v[i] >= ylim[0] && v[i] <= ylim[1])
            .collect();
        let points: Vec<_> = ids.iter().map(|&i| [t[i], v[i]]).collect();
        let values: Vec<_> = ids.iter().map(|&i| s[i]).collect();
        p.scatter_colored(&points, &values, &scale, 1.5, 1., "none")?;
        f.panels.push(p);
        f.labels.push(Label::new(
            [
                "Before code-ghost exclusion",
                "After parent-compatible exclusion",
                "Retained candidates: central Doppler band",
            ][panel],
            95.,
            rect.y - 4.,
        ));
    }
    f.colorbars.push(Colorbar::new(
        Rect::new(172., 19., 3., 144.),
        0,
        "Peak total/noise (dB)",
        vec![5., 10., 15., 20., 25., 30.],
    ));
    save(&f, out, &format!("{date}-doppler-overview"))?;
    let mut f = figure(out, 203.);
    f.scales.push(scale.clone());
    f.labels.push(Label::new(
        format!("Longest retained candidate events: {date}"),
        98.,
        5.,
    ));
    for (k, name) in input.group("events")?.member_names()?.iter().enumerate() {
        let g = input.group(&format!("events/{name}"))?;
        let t = read::<f64>(&g, "time_hours")?;
        let v = read::<f64>(&g, "velocity_mps")?;
        let s = read::<f64>(&g, "peak_total_noise_db")?;
        if t.is_empty() {
            continue;
        }
        let col = k % 2;
        let row = k / 2;
        let rect = Rect::new(24. + col as f64 * 86., 19. + row as f64 * 56., 64., 37.);
        let x0 = t.iter().copied().fold(f64::INFINITY, f64::min);
        let x1 = t.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let dx = ((x1 - x0) * 0.07).max(0.02);
        let y0 = v.iter().copied().fold(f64::INFINITY, f64::min);
        let y1 = v.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let dy = ((y1 - y0) * 0.15).max(2.);
        let mut axes = Axes::new(
            [x0 - dx, x1 + dx],
            [y0 - dy, y1 + dy],
            "Time (UTC)",
            if col == 0 {
                r"Radial velocity (m s$^{-1}$)"
            } else {
                ""
            },
        );
        axes.xticks = vec![x0, (x0 + x1) / 2., x1];
        axes.xticklabels = Some(
            axes.xticks
                .iter()
                .map(|x| {
                    let minutes = (x * 60.).round() as i32;
                    format!("{:02}:{:02}", minutes / 60, minutes % 60)
                })
                .collect(),
        );
        axes.yticks = vec![
            (y0 - dy).ceil(),
            ((y0 + y1) / 2.).round(),
            (y1 + dy).floor(),
        ];
        let mut p = Panel::new(rect, axes);
        p.scatter_colored(
            &t.iter().zip(&v).map(|(&t, &v)| [t, v]).collect::<Vec<_>>(),
            &s,
            &scale,
            2.,
            1.,
            "none",
        )?;
        f.panels.push(p);
        let modes = read::<u64>(&g, "mode_count")?[0];
        f.labels.push(Label::new(
            format!("Event {name}: {modes} modes"),
            rect.x + rect.width / 2.,
            rect.y - 4.,
        ));
    }
    f.colorbars.push(Colorbar::new(
        Rect::new(25., 184., 148., 2.),
        0,
        "Peak total/noise (dB)",
        vec![5., 10., 15., 20., 25., 30.],
    ));
    f.colorbars[0].horizontal = true;
    save(&f, out, &format!("{date}-event-tracks"))?;
    Ok(())
}
