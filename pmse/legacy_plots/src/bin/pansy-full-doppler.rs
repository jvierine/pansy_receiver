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
    if a.len()==3 && a[1]=="--audit-koki" {
        let f=File::open(&a[2])?;
        let velocity=read::<f64>(&f,"velocity_mps")?;
        let start=read::<f64>(&f,"mode_start_hours")?;
        let stop=read::<f64>(&f,"mode_stop_hours")?;
        let cells=read::<u32>(&f,"accepted_extreme_cells_mode_frequency")?;
        ensure!(cells.len()==start.len()*velocity.len() && start.len()==stop.len(),"invalid binwise summary");
        for m in 0..start.len() {
            let row=&cells[m*velocity.len()..(m+1)*velocity.len()];
            if row.iter().any(|&x|x>0) {
                let bins:Vec<_>=row.iter().enumerate().filter(|(_,n)|**n>0).map(|(j,n)|format!("{:.2}({n})",velocity[j])).collect();
                println!("DETECTION mode={m} start_utc_hours={:.6} stop_utc_hours={:.6} velocity_mps(cells)={}",start[m],stop[m],bins.join(","));
            }
        }
        println!("TOTAL tested={} extreme_cells={} detected_modes={} input_errors={}",read::<u64>(&f,"tested_count")?[0],read::<u64>(&f,"extreme_cell_count")?[0],read::<u64>(&f,"detected_mode_count")?[0],read::<u64>(&f,"input_errors")?[0]);
        return Ok(());
    }
    if a.len() == 3 && a[1] == "--audit" {
        let f = File::open(&a[2])?;
        let start = read::<f64>(&f, "mode_start_hours")?;
        let stop = read::<f64>(&f, "mode_stop_hours")?;
        let errors = read::<u64>(&f, "input_errors")?[0];
        println!("AUDIT modes={} input_errors={errors}", start.len());
        let mut count = 0;
        for i in 1..start.len() {
            let gap = (start[i] + stop[i] - start[i - 1] - stop[i - 1]) * 30.;
            if gap > 10. {
                count += 1;
                if count <= 12 {
                    println!(
                        "WHITE_GAP previous_stop_utc_hours={:.6} next_start_utc_hours={:.6} uncovered_minutes={:.3} center_gap_minutes={gap:.3}",
                        stop[i - 1],
                        start[i],
                        (start[i] - stop[i - 1]) * 60.
                    );
                }
            }
        }
        println!("WHITE_GAPS count={count}");
        return Ok(());
    }
    if a.len() == 3 && a[1] == "--audit-index" {
        let f = File::open(&a[2])?;
        let days = read::<u32>(&f, "utc_day_yyyymmdd")?;
        let status = read::<u8>(&f, "status")?;
        for (d, s) in days.iter().zip(status) {
            if s == 2 {
                println!("FAILED_DAY {d}");
            }
        }
        return Ok(());
    }
    let capon_mode = a.len() == 5 && a[3] == "--capon";
    let koki_mode = a.len() == 5 && a[3] == "--koki";
    ensure!(
        a.len() == 3 || (a.len() == 4 && a[3] == "--candidate-context") || capon_mode || koki_mode,
        "usage: pansy-full-doppler FULL_DOPPLER.h5 OUTPUT_DIR [--candidate-context | --capon SIDECAR.h5]"
    );
    let input = File::open(&a[1])?;
    let out = Path::new(&a[2]);
    fs::create_dir_all(out)?;
    let date = String::from_utf8(read(&input, "date_utc")?)?;
    let velocity = read::<f64>(&input, "velocity_mps")?;
    let start = read::<f64>(&input, "mode_start_hours")?;
    let stop = read::<f64>(&input, "mode_stop_hours")?;
    let values = read::<f32>(&input, "peak_total_noise_db_mode_frequency")?;
    let peak_flags = if input.link_exists("extreme_peak_flag") {
        read::<u8>(&input, "extreme_peak_flag")?
    } else {
        vec![0; start.len()]
    };
    let peak_velocity = if input.link_exists("main_peak_velocity_mps") {
        read::<f64>(&input, "main_peak_velocity_mps")?
    } else {
        vec![f64::NAN; start.len()]
    };
    ensure!(
        peak_flags.len() == start.len() && peak_velocity.len() == start.len(),
        "flag shape mismatch"
    );
    let flagged = peak_flags.iter().filter(|&&x| x != 0).count();
    let decisions = if capon_mode {
        let c = File::open(&a[4])?;
        ensure!(
            String::from_utf8(read(&c, "date_utc")?)? == date,
            "Capon date mismatch"
        );
        read::<u8>(&c, "capon_decision")?
    } else {
        vec![0; start.len()]
    };
    ensure!(decisions.len() == start.len(), "Capon mode count mismatch");
    let detections = decisions.iter().filter(|&&x| x == 1).count();
    let koki = if koki_mode {
        Some(File::open(&a[4])?)
    } else {
        None
    };
    let cells = if let Some(k) = &koki {
        ensure!(
            String::from_utf8(read(k, "date_utc")?)? == date,
            "Koki date mismatch"
        );
        let v = read::<f64>(k, "velocity_mps")?;
        ensure!(
            v.len() == velocity.len() && v.iter().zip(&velocity).all(|(a, b)| (a - b).abs() < 1e-6),
            "Koki axis mismatch"
        );
        let c = read::<u32>(k, "accepted_extreme_cells_mode_frequency")?;
        ensure!(
            c.len() == start.len() * velocity.len(),
            "Koki shape mismatch"
        );
        c
    } else {
        vec![]
    };
    if input.link_exists("main_peak_range_km") {
        let ranges = read::<f64>(&input, "main_peak_range_km")?;
        let powers = read::<f32>(&input, "main_peak_total_noise_db")?;
        for m in 0..start.len() {
            if peak_flags[m] != 0 {
                println!(
                    "PEAK_FLAG date={date} center_utc_hours={:.6} velocity_mps={:.3} range_km={:.3} total_noise_db={:.2}",
                    (start[m] + stop[m]) / 2.,
                    peak_velocity[m],
                    ranges[m],
                    powers[m]
                );
            }
        }
    }
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
        for m in 0..start.len() {
            if koki_mode {
                for j in 0..velocity.len() {if cells[m*velocity.len()+j]>0 {
                    let x=centers[m];let y=velocity[j];let dx=(xlim[1]-xlim[0])*0.0017;let dy=limit*0.009;
                    let mut marker=Line::new(vec![[x-dx,y],[x,y+dy],[x+dx,y],[x,y-dy],[x-dx,y]],"#00e6a8");
                    marker.width_pt=0.7;p.push(Layer::Line(marker));
                }}
                continue;
            }
            if peak_flags[m] != 0 {
                let x = centers[m];
                let y = peak_velocity[m];
                let dx = (xlim[1] - xlim[0]) * 0.003;
                let dy = limit * 0.015;
                let points = if capon_mode && decisions[m] == 1 {
                    vec![
                        [x - dx, y],
                        [x, y + dy],
                        [x + dx, y],
                        [x, y - dy],
                        [x - dx, y],
                    ]
                } else {
                    vec![
                        [x - dx, y - dy],
                        [x + dx, y + dy],
                        [f64::NAN, f64::NAN],
                        [x - dx, y + dy],
                        [x + dx, y - dy],
                    ]
                };
                let mut marker = Line::new(
                    points,
                    if capon_mode && decisions[m] != 1 {
                        "#a0a0a0"
                    } else {
                        "#00e6a8"
                    },
                );
                marker.width_pt = if capon_mode && decisions[m] != 1 {
                    0.6
                } else {
                    1.5
                };
                p.push(Layer::Line(marker));
            }
        }
        Ok(p)
    };
    let mut f = figure(out, 112.);
    f.raster_dpi = if koki_mode {300.} else {600.};
    f.scales.push(scale.clone());
    f.labels.push(Label::new(
        format!("PANSY {date}: full stored Doppler spectrum"),
        98.,
        5.,
    ));
    f.labels.push(Label::new(
        if koki_mode {
            format!("XC2 zenith, 80--92 km; binwise leakage + Capon; |v| >= 15 m/s marked")
        } else if capon_mode {
            format!("XC2 zenith, 80--92 km; {detections} Capon detections; cyan diamonds, grey velocity flags")
        } else if input.link_exists("extreme_peak_flag") {
            format!("Zenith, 80--92 km; {flagged} flags: main peak within +/-90 has |v| >= 30 m/s")
        } else {
            "Zenith, 80--92 km; one column per integration; all Dopplers".to_string()
        },
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
    let plot_name = if koki_mode {
        format!("{date}-XC2-KOKI-overview-v1")
    } else if capon_mode {
        format!("{date}-XC2-CAPON-overview-v2")
    } else {
        format!("{date}-FULL-doppler-overview")
    };
    save(&f, out, &plot_name)?;
    if let Some(k) = &koki {
        let vx = read::<f64>(k, "histogram_velocity_mps")?;
        let series = [
            ("histogram_raw", "#737373", "Raw bins"),
            ("histogram_leakage_cleaned", "#df9518", "After leakage"),
            ("histogram_capon", "#087b9a", "After Capon"),
        ];
        let max = read::<u64>(k, "histogram_raw")?
            .into_iter()
            .max()
            .unwrap_or(1)
            .max(1) as f64;
        let mut h = figure(out, 142.);
        h.labels.push(Label::new(
            format!("PANSY {date}: zenith velocity distributions"),
            98.,
            5.,
        ));
        for row in 0..2 {
            let mut ax = Axes::new(
                [-90., 90.],
                [
                    0.,
                    if row == 0 {
                        max * 1.08
                    } else {
                        max.log10().max(1.) * 1.08
                    },
                ],
                "Radial velocity (m/s)",
                if row == 0 {
                    "Bin count"
                } else {
                    "log10(bin count)"
                },
            );
            ax.xticks = vec![-90., -60., -30., 0., 30., 60., 90.];
            let mut p = Panel::new(Rect::new(25., 15. + row as f64 * 60., 145., 45.), ax);
            for (name, color, _) in &series {
                let c = read::<u64>(k, name)?;
                let points = vx
                    .iter()
                    .zip(&c)
                    .map(|(&v, &n)| {
                        [
                            v,
                            if row == 0 {
                                n as f64
                            } else if n > 0 {
                                (n as f64).log10()
                            } else {
                                f64::NAN
                            },
                        ]
                    })
                    .collect();
                let mut line = Line::new(points, *color);
                line.width_pt = 0.9;
                p.push(Layer::Line(line));
            }
            h.panels.push(p);
        }
        for (i, (_, color, label)) in series.iter().enumerate() {
            let mut l = Label::new(*label, 48. + i as f64 * 48., 136.);
            l.color = (*color).to_owned();
            h.labels.push(l);
        }
        save(&h, out, &format!("{date}-XC2-KOKI-histogram-v1"))?;
        let mut html=format!("<!doctype html><meta charset=utf-8><title>PANSY {date} event spectra</title><style>body{{font:17px system-ui;max-width:1100px;margin:2em auto;padding:0 1em;background:white;color:#172432}}img{{width:100%;height:auto}}article{{margin:2em 0}}</style><h1>PANSY {date}: detected integrations</h1><p>One plot per flagged integration, with all accepted extreme Doppler bins marked. These are not independent-event classifications. Full native Doppler range; range-max power over 80–92 km.</p>");
        for m in 0..start.len() {
            let row=&cells[m*nf..(m+1)*nf];
            if !row.iter().any(|&n|n>0) {continue;}
            let seconds=(centers[m]*3600.).round() as usize;
            let time=format!("{:02}:{:02}:{:02}",seconds/3600,(seconds/60)%60,seconds%60);
            let stem=format!("{date}-XC2-KOKI-event-{m:06}-v2");
            if !out.join(format!("{stem}.png")).exists() {
                let used=fs::read_dir(out)?.filter_map(|e|e.ok()).filter(|e|e.file_name().to_string_lossy().contains("-XC2-KOKI-")).filter_map(|e|e.metadata().ok()).map(|m|m.len()).sum::<u64>();
                ensure!(used+1_000_000<=190_000_000,"event figure budget reached; no further rendering");
                let mut e=figure(out,151.);e.raster_dpi=180.;
                let context_peak=centers.iter().enumerate().filter(|(_,t)|(**t-centers[m]).abs()<=0.5).flat_map(|(i,_)|values[i*nf..(i+1)*nf].iter().copied()).fold(15f32,f32::max) as f64;
                let upper=(context_peak/5.).ceil()*5.;
                e.scales.push(Scale::new([0.,upper],Normalization::Linear,ColorMap::plasma())?);
                e.labels.push(Label::new(format!("PANSY {date}, {time} UTC"),98.,5.));
                e.labels.push(Label::new("Zenith, 80--92 km; accepted Doppler bins marked",98.,11.));
                e.panels.push(panel(Rect::new(25.,20.,145.,50.),[(centers[m]-0.5).max(0.),(centers[m]+0.5).min(24.)])?);
                let step=if upper<=40. {5.} else {10.};
                let mut cb=Colorbar::new(Rect::new(25.,84.,145.,2.),0,"Peak total/noise (dB)",(0..=(upper/step) as usize).map(|i|i as f64*step).collect());cb.horizontal=true;e.colorbars.push(cb);
                let peak=values[m*nf..(m+1)*nf].iter().copied().fold(10f32,f32::max) as f64;
                let mut ax=Axes::new([-limit,limit],[0.,(peak+2.).max(15.)],"Radial velocity (m/s)","Peak total/noise (dB)");
                ax.xticks=vec![-limit,-limit/2.,0.,limit/2.,limit];
                let mut p=Panel::new(Rect::new(25.,108.,145.,28.),ax);
                p.push(Layer::Line(Line::new(velocity.iter().enumerate().map(|(j,&v)|[v,values[m*nf+j] as f64]).collect(),"#31205c")));
                for j in 0..nf {if row[j]>0 {
                    let v=velocity[j];let y=values[m*nf+j] as f64;
                    let mut marker=Line::new(vec![[v-0.8,y],[v,y+0.4],[v+0.8,y],[v,y-0.4],[v-0.8,y]],"#008b73");marker.width_pt=1.;p.push(Layer::Line(marker));
                }}
                e.panels.push(p);e.save(out.join(format!("{stem}.png")))?;
            }
            html.push_str(&format!("<article><h2>{time} UTC: {} accepted extreme cells</h2><a href='{stem}.png'><img loading='lazy' src='{stem}.png' alt='{date} {time} full Doppler event context'></a></article>",row.iter().sum::<u32>()));
        }
        fs::write(out.join(format!("{date}-XC2-KOKI-events-v2.html")),html)?;
    }
    // The overview is the default product. Event ranking is exploratory and
    // must never silently determine which spectra are presented to the user.
    if a.len() == 3 || capon_mode || koki_mode {
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
