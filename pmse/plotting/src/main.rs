//! All axes, ticks, labels and colorbars are rustplotlib Figure objects.
//! The live canvas overlays the exact panel rectangle exported here.
use rustplotlib::*;
use anyhow::Result;
use serde_json::{json,Value};
use std::{fs,path::Path};

fn save(path: &Path, text: &str) -> Result<()> {
    let tmp=path.with_extension("partial"); fs::write(&tmp,text)?; fs::rename(tmp,path)?; Ok(())
}
fn main()->Result<()> {
    let args:Vec<_>=std::env::args().collect();
    anyhow::ensure!(args.len()==3,"Usage: pansy-pmse-plots manifest.json output_directory");
    let spec:Value=serde_json::from_slice(&fs::read(&args[1])?)?;
    let out=Path::new(&args[2]);fs::create_dir_all(out)?;
    let hours=(spec["latest_mode_us"].as_f64().unwrap()-spec["first_mode_us"].as_f64().unwrap())/3.6e9;
    let labels:Vec<String>=serde_json::from_value(spec["time_tick_labels"].clone())?;
    let mut layouts=serde_json::Map::new();
    for width in [360.,640.,960.] {
      for (low,high) in [(3.,6.),(3.,10.),(3.,15.),(3.,20.),(3.,30.),(3.,40.)] {
       for kind in ["rti","rd"] {
        let mut f=WebLayout::new(width,16.)?.figure(300.)?;
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
            ColorMap::new(vec![[0.08,0.08,0.16],[0.35,0.13,0.31],[0.71,0.24,0.25],[0.97,0.61,0.28],[1.,0.94,0.73]])?)?);
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
      }
    }
    save(&out.join("layouts.json"),&serde_json::to_string(&layouts)?)?;
    if spec.get("noise_model_data").is_some() {
        let data=&spec["noise_model_data"];
        let times:Vec<f64>=serde_json::from_value(data["times_h"].clone())?;
        let sky:Vec<Vec<f64>>=serde_json::from_value(data["sky_temperature_k"].clone())?;
        let power:Vec<Vec<Vec<f64>>>=serde_json::from_value(data["noise_power"].clone())?;
        let labels:Vec<String>=serde_json::from_value(data["time_tick_labels"].clone())?;
        let colors=["#202020","#2474a8","#b65c18","#22854a","#b52e48"];
        let names=["Zenith","North","East","South","West"];
        for width in [640.,960.] {for channel in 0..7 {
            let mut f=WebLayout::new(width,16.)?.figure(550.)?;f.latex.cache=out.join(".latex-cache");
            let xmax=times.last().copied().unwrap_or(0.).max(0.01);
            for panel in 0..2 {
                let mut values=vec![];
                for t in 0..times.len(){for b in 0..5{values.push(if panel==0{10.*power[t][b][channel].log10()}else{sky[t][b]});}}
                let lo=values.iter().copied().fold(f64::INFINITY,f64::min);
                let hi=values.iter().copied().fold(f64::NEG_INFINITY,f64::max);
                let pad=((hi-lo)*0.15).max(if panel==0{0.1}else{10.});
                let mut axes=Axes::new([0.,xmax],[lo-pad,hi+pad],"Time (UTC)",
                    if panel==0{"Receiver noise power (dB)"}else{r"Modeled $T_{\rm sky}$ (K)"});
                axes.xticks=vec![0.,xmax/2.,xmax];axes.xticklabels=Some(labels.clone());
                let rect=Rect::new(20.,10.+panel as f64*68.,f.width_mm-28.,43.);
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
    if spec["transfer"]["archive"]["hourly_unix_s"].is_array() {
        let a=&spec["transfer"]["archive"];
        let hours:Vec<f64>=serde_json::from_value(a["hourly_unix_s"].clone())?;
        let bytes:Vec<f64>=serde_json::from_value(a["hourly_bytes"].clone())?;
        let labels:Vec<String>=serde_json::from_value(a["hour_tick_labels"].clone())?;
        if !hours.is_empty(){for width in [640.,960.] {
            let mut f=WebLayout::new(width,16.)?.figure(550.)?;f.latex.cache=out.join(".latex-cache");
            let xmax=((a["latest_mode_us"].as_f64().unwrap()/1e6-hours[0])/3600.).max(0.01);
            for panel in 0..2 {
                let ymax=if panel==0{bytes.iter().copied().fold(0_f64,f64::max)/1048576.}else{bytes.iter().sum::<f64>()/1048576.};
                let mut axes=Axes::new([0.,xmax],[0.,(ymax*1.15).max(0.001)],"Time (UTC)",
                    if panel==0{"Data per UTC hour (MiB)"}else{"Accumulated data (MiB)"});
                axes.xticks=vec![0.,xmax/2.,xmax];axes.xticklabels=Some(labels.clone());
                let mut p=Panel::new(Rect::new(20.,8.+panel as f64*68.,f.width_mm-28.,43.),axes);
                let mut sum=0.;let mut cumulative=vec![[0.,0.]];
                for i in 0..hours.len(){let x=(hours[i]-hours[0])/3600.;let right=(x+1.).min(xmax);let y=bytes[i]/1048576.;
                    if panel==0{let mut bar=Line::new(vec![[x,0.],[right,0.],[right,y],[x,y],[x,0.]],"#215d86");bar.fill=Some("#5a9ec4".into());p.push(Layer::Line(bar));}
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
