//! All axes, ticks, labels and colorbars are rustplotlib Figure objects.
//! The live canvas overlays the exact panel rectangle exported here.
use rustplotlib::*;
use anyhow::Result;
use serde_json::{json,Value};
use std::{fs,path::Path};

const WAVE_PAPER:[f64;3]=[250./255.,247./255.,237./255.];
// Scientific-notation tick labels (e.g. 1.00 x 10^4 K) exceed the former
// 20 mm gutter. Keep geometry checks enabled and reserve actual label space.
const AUX_LEFT_MM:f64=32.;

fn auxiliary_ylabel(figure:&mut Figure,rect:Rect,text:&str) {
    // Axes' fixed default title offset assumes short ticks. Use a native Label
    // in a dedicated outer gutter, separate from scientific-notation ticks.
    let mut label=Label::new(text,4.,rect.y+rect.height/2.);
    label.rotation=-90.;figure.labels.push(label);
}

fn wave_figure(width:f64,height:f64)->Result<Figure> {
    let mut figure=WebLayout::new(width,16.)?.figure(height)?;
    figure.background=WAVE_PAPER;
    Ok(figure)
}

fn save(path: &Path, text: &str) -> Result<()> {
    let tmp=path.with_extension("partial"); fs::write(&tmp,text)?; fs::rename(tmp,path)?; Ok(())
}
fn wave_palette()->Result<ColorMap> {
    ColorMap::new(serde_json::from_str(include_str!("../../rvd/palette.json"))?)
}
fn main()->Result<()> {
    let args:Vec<_>=std::env::args().collect();
    if args.get(1).map(String::as_str)==Some("--frame") {
        let mut value:Value=serde_json::from_reader(std::io::stdin())?;
        value["spec"]["color_stops"]=serde_json::from_str(include_str!("../../rvd/palette.json"))?;
        value["spec"]["background"]=json!(WAVE_PAPER);
        let spec=serde_json::from_value(value["spec"].clone())?;
        let request=serde_json::from_value(value["request"].clone())?;
        let cache=std::env::var("PANSY_FRAME_CACHE").unwrap_or_else(|_|"/tmp/pansy-native-frame-cache".into());
        let frame=rustplotlib::interactive::field_frame(&spec,&request,Path::new(&cache))?;
        println!("{}",serde_json::to_string(&frame)?);return Ok(());
    }
    let archive_only=args.get(1).map(String::as_str)==Some("--archive-only");
    let noise_only=args.get(1).map(String::as_str)==Some("--noise-only");
    let offset=usize::from(archive_only||noise_only);
    anyhow::ensure!(args.len()==3+offset,"Usage: pansy-pmse-plots [--archive-only|--noise-only] manifest.json output_directory");
    let spec:Value=serde_json::from_slice(&fs::read(&args[1+offset])?)?;
    let out=Path::new(&args[2+offset]);fs::create_dir_all(out)?;
    // The metadata-only web publisher needs storage diagnostics, not dozens
    // of static RTI frames or any radar/complex-voltage processing.
    if !archive_only&&!noise_only {
    let hours=(spec["latest_mode_us"].as_f64().unwrap()-spec["first_mode_us"].as_f64().unwrap())/3.6e9;
    let labels:Vec<String>=serde_json::from_value(spec["time_tick_labels"].clone())?;
    let mut layouts=serde_json::Map::new();
    for width in [360.,640.,960.] {
      for (low,high) in [(3.,6.),(3.,10.),(3.,15.),(3.,20.),(3.,30.),(3.,40.)] {
       for kind in ["rti","rd"] {
        let mut f=wave_figure(width,300.)?;
        f.latex.cache=out.join(".latex-cache");
        let rect=Rect::new(20.,7.,f.width_mm-50.,f.height_mm-24.);
        let mut axes=Axes::new(if kind=="rti"{[0.,hours.max(0.001)]}else{[-90.321448,90.321448]},
            [50.,120.],if kind=="rti"{"Time (UTC)"}else{r"Radial velocity (m s$^{-1}$)"},"Nominal altitude (km)");
        axes.yticks=vec![50.,70.,90.,110.];
        if kind=="rti" {
            axes.xticks=vec![0.,hours/2.,hours];axes.xticklabels=Some(labels.clone());
        } else { axes.xticks=vec![-90.,0.,90.]; }
        f.panels.push(Panel::new(rect,axes));
        f.scales.push(Scale::new([low,high],Normalization::Linear,
            wave_palette()?)?);
        let ticks=if high==6. {vec![3.,4.,5.,6.]} else {
            let mut ticks=vec![3.];
            for v in [5.,10.,20.,30.,40.] { if v<=high {ticks.push(v);} }
            if !ticks.contains(&high){ticks.push(high);} ticks
        };
        f.colorbars.push(Colorbar::new(Rect::new(rect.x+rect.width+3.,rect.y,3.,rect.height),0,"SNR (dB)",ticks));
        f.layout_colorbars()?;
        let key=format!("{kind}-{}-{low}-{high}",width as i32);
        let file=format!("{key}.svg");save(&out.join(&file),&f.svg()?)?;
       layouts.insert(key,json!({"file":file,"width":width,"height":300.,
            "panel":[rect.x/f.width_mm,rect.y/f.height_mm,rect.width/f.width_mm,rect.height/f.height_mm]}));
       }
       let mut f=wave_figure(width,700.)?;f.latex.cache=out.join(".latex-cache");
       f.scales.push(Scale::new([low,high],Normalization::Linear,
           wave_palette()?)?);
       let mut panels=vec![];let mut beam_labels=vec![];
       for beam in 0..5 {
           let rect=Rect::new(20.,4.+beam as f64*32.,f.width_mm-50.,30.);
           let mut axes=Axes::new([0.,hours.max(0.001)],[50.,120.],"Time (UTC)",
               if beam==2{"Nominal altitude (km)"}else{""});
           axes.yticks=vec![60.,90.,120.];axes.show_x_axis=beam==4;
           axes.xticks=vec![0.,hours/2.,hours];axes.xticklabels=Some(labels.clone());
           f.panels.push(Panel::new(rect,axes));
           panels.push([rect.x/f.width_mm,rect.y/f.height_mm,rect.width/f.width_mm,rect.height/f.height_mm]);
           let filename=format!("beam-{beam}.svg");
           let mut label=wave_figure(80.,24.)?;label.latex.cache=f.latex.cache.clone();
           label.labels.push(Label::new(["Zenith","North","East","South","West"][beam],label.width_mm/2.,label.height_mm/2.));
           save(&out.join(&filename),&label.svg()?)?;beam_labels.push(filename);
       }
       let r=f.panels[0].rect;
       let ticks=if high==6.{vec![3.,4.,5.,6.]}else{vec![3.,high]};
       f.colorbars.push(Colorbar::new(Rect::new(r.x+r.width+3.,r.y,3.,158.),0,"SNR (dB)",ticks));
       f.layout_colorbars()?;
       let key=format!("rti_stack-{}-{low}-{high}",width as i32);let file=format!("{key}.svg");
       save(&out.join(&file),&f.svg()?)?;
       layouts.insert(key,json!({"file":file,"width":width,"height":700.,"panels":panels,"beam_labels":beam_labels}));
      }
    }
    save(&out.join("layouts.json"),&serde_json::to_string(&layouts)?)?;
    }
    if !archive_only&&spec.get("noise_model_data").is_some() {
        let data=&spec["noise_model_data"];
        let times:Vec<f64>=serde_json::from_value(data["times_h"].clone())?;
        let sky:Vec<Vec<f64>>=serde_json::from_value(data["sky_temperature_k"].clone())?;
        let power:Vec<Vec<Vec<f64>>>=serde_json::from_value(data["noise_power"].clone())?;
        let labels:Vec<String>=serde_json::from_value(data["time_tick_labels"].clone())?;
        let colors=["#202020","#2474a8","#b65c18","#22854a","#b52e48"];
        let names=["Zenith","North","East","South","West"];
        for width in [640.,960.] {for channel in 0..7 {
            let mut f=wave_figure(width,550.)?;f.latex.cache=out.join(".latex-cache");
            let xmax=times.last().copied().unwrap_or(0.).max(0.01);
            for panel in 0..2 {
                let mut values=vec![];
                for t in 0..times.len(){for b in 0..5{values.push(if panel==0{10.*power[t][b][channel].log10()}else{sky[t][b]});}}
                let lo=values.iter().copied().fold(f64::INFINITY,f64::min);
                let hi=values.iter().copied().fold(f64::NEG_INFINITY,f64::max);
                let pad=((hi-lo)*0.15).max(if panel==0{0.1}else{10.});
                let mut axes=Axes::new([0.,xmax],[lo-pad,hi+pad],"Time (UTC)","");
                axes.xticks=vec![0.,xmax/2.,xmax];axes.xticklabels=Some(labels.clone());
                // Multi-day UTC tick labels include the date. Reserve their
                // half-width at the right edge without rescaling typography.
                let rect=Rect::new(AUX_LEFT_MM,10.+panel as f64*68.,f.width_mm-AUX_LEFT_MM-16.,43.);
                auxiliary_ylabel(&mut f,rect,if panel==0{"Receiver noise power (dB)"}else{r"Modeled $T_{\rm sky}$ (K)"});
                let mut p=Panel::new(rect,axes);
                for beam in 0..5 {
                    let points:Vec<_>=(0..times.len()).map(|t|[times[t],if panel==0{10.*power[t][beam][channel].log10()}else{sky[t][beam]}]).collect();
                    p.push(Layer::Line(Line::new(points.clone(),colors[beam])));
                    let rgb=match beam {0=>[0.12,0.12,0.12],1=>[0.14,0.45,0.66],2=>[0.71,0.36,0.09],3=>[0.13,0.52,0.29],_=>[0.71,0.18,0.28]};
                    p.markers(&mut f,&points,&MarkerStyle{shape:MarkerShape::Circle,size_pt:2.8,color:rgb})?;
                }
                if panel==0 {
                    let mut legend=Legend::new(names.iter().zip(colors).map(|(n,c)|LegendEntry::line(*n,c)).collect());
                    legend.columns=5;legend.position=LegendPosition::AboveCenter;p.legend=Some(legend);
                }
                f.panels.push(p);
            }
            save(&out.join(format!("noise-module-{channel}-{}.svg",width as i32)),&f.svg()?)?;
        }}
    }
    if !noise_only&&spec["transfer"]["archive"]["hourly_unix_s"].is_array() {
        let a=&spec["transfer"]["archive"];
        let hours:Vec<f64>=serde_json::from_value(a["hourly_unix_s"].clone())?;
        let bytes:Vec<f64>=serde_json::from_value(a["hourly_bytes"].clone())?;
        let labels:Vec<String>=serde_json::from_value(a["hour_tick_labels"].clone())?;
        let retention:Vec<Option<f64>>=serde_json::from_value(a["hourly_dense_retained_fraction"].clone())?;
        if !hours.is_empty(){for width in [640.,960.] {
            let mut f=wave_figure(width,810.)?;f.latex.cache=out.join(".latex-cache");
            let xmax=((a["latest_mode_us"].as_f64().unwrap()/1e6-hours[0])/3600.).max(0.01);
            for panel in 0..3 {
                let ymax=match panel {
                    0=>bytes.iter().copied().fold(0_f64,f64::max)/1048576.,
                    1=>bytes.iter().sum::<f64>()/1048576.,
                    _=>retention.iter().flatten().copied().fold(0_f64,f64::max)*100.,
                };
                let ymin=if panel==2 {-(ymax*0.04).max(0.00004)} else {0.};
                let mut axes=Axes::new([0.,xmax],[ymin,(ymax*1.15).max(0.001)],"Time (UTC)","");
                axes.xticks=vec![0.,xmax/2.,xmax];axes.xticklabels=Some(labels.clone());
                let rect=Rect::new(AUX_LEFT_MM,8.+panel as f64*68.,f.width_mm-AUX_LEFT_MM-16.,43.);
                auxiliary_ylabel(&mut f,rect,match panel {0=>"Data per UTC hour (MiB)",1=>"Accumulated data (MiB)",_=>"Full spectra retained (\\%)"});
                let mut p=Panel::new(rect,axes);
                let mut sum=0.;let mut cumulative=vec![[0.,0.]];
                for i in 0..hours.len(){let x=(hours[i]-hours[0])/3600.;let right=(x+1.).min(xmax);let y=bytes[i]/1048576.;
                    if panel==0{let mut bar=Line::new(vec![[x,0.],[right,0.],[right,y],[x,y],[x,0.]],"#215d86");bar.fill=Some("#5a9ec4".into());p.push(Layer::Line(bar));}
                    if panel==2 {if let Some(Some(fraction))=retention.get(i) {
                        let percent=fraction*100.;
                        let mut bar=Line::new(vec![[x,0.],[right,0.],[right,percent],[x,percent],[x,0.]],"#215d86");
                        bar.fill=Some("#5a9ec4".into());p.push(Layer::Line(bar));
                        if percent==0. {p.markers(&mut f,&[[(x+right)/2.,0.]],&MarkerStyle{shape:MarkerShape::Circle,size_pt:2.8,color:[0.13,0.36,0.53]})?;}
                    }}
                    sum+=y;cumulative.push([right,sum]);
                }
                if panel==1{p.push(Layer::Line(Line::new(cumulative,"#215d86")));}
                f.panels.push(p);
            }
            save(&out.join(format!("archive-growth-{}.svg",width as i32)),&f.svg()?)?;
        }}
    }
    Ok(())
}
