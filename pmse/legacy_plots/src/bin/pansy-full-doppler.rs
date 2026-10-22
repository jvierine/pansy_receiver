//! Full native spectral context. No threshold, candidate mask, velocity cut,
//! Doppler interpolation or cropping. Raster fields and all labels are native
//! Rustplotlib. One display column per mode, bounded by adjacent midpoints;
//! genuine long gaps (>10 minutes between integration centers) remain blank.
use anyhow::{Result, ensure};
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
    f.raster_dpi = 220.;
    f.latex.cache = out.join(".latex-cache");
    f
}
fn save(f: &Figure, out: &Path, name: &str) -> Result<()> {
    for ext in ["png", "pdf"] {
        f.save(out.join(format!("{name}.{ext}")))?;
    }
    Ok(())
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    ensure!(
        a.len() == 3 || (a.len() == 4 && a[3] == "--candidate-context"),
        "usage: pansy-full-doppler FULL_DOPPLER.h5 OUTPUT_DIR [--candidate-context]"
    );
    let input = File::open(&a[1])?;
    let out = Path::new(&a[2]);
    fs::create_dir_all(out)?;
    let date = String::from_utf8(read(&input, "date_utc")?)?;
    let velocity = read::<f64>(&input, "velocity_mps")?;
    let start = read::<f64>(&input, "mode_start_hours")?;
    let stop = read::<f64>(&input, "mode_stop_hours")?;
    let values = read::<f32>(&input, "peak_total_noise_db_mode_frequency")?;
    ensure!(
        values.len() == start.len() * velocity.len(),
        "full native spectrum shape mismatch"
    );
    let centers: Vec<_> = start.iter().zip(&stop).map(|(a, b)| (a + b) / 2.).collect();
    ensure!(
        centers.windows(2).all(|w| w[1] > w[0]),
        "integration times not increasing"
    );
    let mut xe = vec![0.];
    let mut columns = Vec::new();
    for m in 0..start.len() {
        let left = if m > 0 && centers[m] - centers[m - 1] <= 1. / 6. {
            (centers[m - 1] + centers[m]) / 2.
        } else {
            start[m]
        };
        let right = if m + 1 < start.len() && centers[m + 1] - centers[m] <= 1. / 6. {
            (centers[m] + centers[m + 1]) / 2.
        } else {
            stop[m]
        };
        if left > *xe.last().unwrap() + 1e-10 {
            xe.push(left);
            columns.push(None);
        }
        ensure!(
            right > *xe.last().unwrap(),
            "invalid integration display edges"
        );
        xe.push(right);
        columns.push(Some(m));
    }
    if *xe.last().unwrap() < 24. {
        xe.push(24.);
        columns.push(None);
    }
    let step = velocity[1] - velocity[0];
    let mut ye: Vec<_> = velocity.iter().map(|v| v - step / 2.).collect();
    ye.push(velocity.last().unwrap() + step / 2.);
    let nx = xe.len() - 1;
    let nf = velocity.len();
    let mut pixels = vec![f64::NAN; nx * nf];
    for x in 0..nx {
        if let Some(m) = columns[x] {
            for v in 0..nf {
                pixels[v * nx + x] = values[m * nf + v] as f64;
            }
        }
    }
    let scale = Scale::new([0., 30.], Normalization::Linear, ColorMap::plasma())?;
    let limit = (ye.iter().map(|v| v.abs()).fold(0., f64::max) / 10.).ceil() * 10.;
    let panel = |rect: Rect, xlim: [f64; 2]| -> Result<Panel> {
        let mut axes = Axes::new(
            xlim,
            [-limit, limit],
            "Time (UTC)",
            r"Radial velocity (m s$^{-1}$)",
        );
        axes.yticks = vec![-limit, -limit / 2., 0., limit / 2., limit];
        axes.xticks = if xlim == [0., 24.] {
            vec![0., 4., 8., 12., 16., 20., 24.]
        } else {
            vec![xlim[0], (xlim[0] + xlim[1]) / 2., xlim[1]]
        };
        axes.xticklabels = Some(
            axes.xticks
                .iter()
                .map(|h| {
                    let min = (h * 60.).round() as i32;
                    format!("{:02}:{:02}", min / 60, min % 60)
                })
                .collect(),
        );
        let mut p = Panel::new(rect, axes);
        p.push(Layer::Pcolormesh(Field::new(
            xe.clone(),
            ye.clone(),
            pixels.clone(),
            0,
        )?));
        Ok(p)
    };
    let mut f = figure(out, 112.);
    f.raster_dpi = 600.;
    f.scales.push(scale.clone());
    f.labels.push(Label::new(
        format!("PANSY {date}: full stored Doppler spectrum"),
        98.,
        5.,
    ));
    f.labels.push(Label::new(
        "Zenith, 80--92 km; one column per integration; all Dopplers",
        98.,
        11.,
    ));
    f.panels
        .push(panel(Rect::new(24., 20., 147., 62.), [0., 24.])?);
    let mut bar = Colorbar::new(
        Rect::new(24., 99., 147., 2.),
        0,
        "Peak total/noise (dB)",
        vec![0., 5., 10., 15., 20., 25., 30.],
    );
    bar.horizontal = true;
    f.colorbars.push(bar);
    // Horizontal colorbar labels need a real lower margin at article font size.
    f.height_mm = 119.;
    save(&f, out, &format!("{date}-FULL-doppler-overview"))?;
    // The overview is the default product. Event ranking is exploratory and
    // must never silently determine which spectra are presented to the user.
    if a.len() == 3 {
        println!(
            "Rendered {date}: ALL {nf} native Dopplers, {} modes, no candidate mask",
            start.len()
        );
        return Ok(());
    }
    let mut f = figure(out, 205.);
    f.scales.push(scale.clone());
    f.labels.push(Label::new(
        format!("Candidate intervals: full Doppler context, {date}"),
        98.,
        5.,
    ));
    for (k, name) in input.group("events")?.member_names()?.iter().enumerate() {
        let g = input.group(&format!("events/{name}"))?;
        let x0 = read::<f64>(&g, "start_hours")?[0];
        let x1 = read::<f64>(&g, "stop_hours")?[0];
        let rect = Rect::new(
            24. + (k % 2) as f64 * 86.,
            19. + (k / 2) as f64 * 56.,
            64.,
            37.,
        );
        f.panels.push(panel(rect, [x0, x1])?);
        f.labels.push(Label::new(
            format!("Candidate {name}"),
            rect.x + rect.width / 2.,
            rect.y - 4.,
        ));
        let mut individual = figure(out, 119.);
        individual.raster_dpi = 600.;
        individual.scales.push(scale.clone());
        individual.labels.push(Label::new(
            format!("Candidate {name}: full Doppler spectrum, {date}"),
            98.,
            5.,
        ));
        individual
            .panels
            .push(panel(Rect::new(24., 20., 147., 62.), [x0, x1])?);
        let mut cb = Colorbar::new(
            Rect::new(24., 99., 147., 2.),
            0,
            "Peak total/noise (dB)",
            vec![0., 5., 10., 15., 20., 25., 30.],
        );
        cb.horizontal = true;
        individual.colorbars.push(cb);
        save(
            &individual,
            out,
            &format!("{date}-candidate-{name}-FULL-doppler"),
        )?;
    }
    let mut bar = Colorbar::new(
        Rect::new(25., 184., 148., 2.),
        0,
        "Peak total/noise (dB)",
        vec![0., 5., 10., 15., 20., 25., 30.],
    );
    bar.horizontal = true;
    f.colorbars.push(bar);
    save(&f, out, &format!("{date}-candidate-FULL-doppler"))?;
    println!(
        "Rendered ALL {nf} Doppler bins for {} modes; min={:.3} max={:.3} m/s",
        start.len(),
        velocity[0],
        velocity[nf - 1]
    );
    Ok(())
}
