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
        let ticks=(0..4).map(|i|low+(high-low)*i as f64/3.).collect();
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
    Ok(())
}
