//! Browser-only scientific imaging. Transport and UI contain no numerical inversion.
use pansy_pmse_dsp::imaging::{altitude_pixel,calibrated_covariance,channel_pairs,TsvdImager,doppler_moments};
use num_complex::Complex64;
use serde::Deserialize;
use wasm_bindgen::prelude::*;

#[derive(Deserialize)]
struct Geometry {
    positions_m: [[f64;3];7], phasecal_rad: [[f64;7];5],
    beam_azimuth_deg: [f64;5], beam_elevation_deg:[f64;5],
    frequency_hz:f64,
}
#[derive(Deserialize)]
struct Packet {
    index:Vec<[usize;3]>, ranges_km:Vec<f64>, velocity_mps:Vec<f64>,
    noise_power:Vec<Vec<Vec<f64>>>, looks:usize, nyquist_mps:f64,
}
struct Beam { pixels:Vec<usize>, imager:TsvdImager, weights:Vec<f64> }

#[wasm_bindgen]
pub struct BrowserImager { beams:Vec<Beam>, ranges:Vec<f64>, pixels:usize, step:f64, nyquist:f64 }

fn err(s:impl ToString)->JsValue { JsValue::from_str(&s.to_string()) }

#[wasm_bindgen]
impl BrowserImager {
    #[wasm_bindgen(constructor)]
    pub fn new(geometry_json:&str, height_km:f64, grid:usize, fixed_slant_range:bool)->Result<BrowserImager,JsValue> {
        if !(21..=101).contains(&grid) || !height_km.is_finite() || !(50. ..=120.).contains(&height_km) {
            return Err(err("Invalid image grid/altitude"));
        }
        let g:Geometry=serde_json::from_str(geometry_json).map_err(err)?;
        let extent=height_km*if fixed_slant_range {20f64.to_radians().sin()} else {20f64.to_radians().tan()};let step=2.*extent/(grid-1) as f64;
        let mut ranges=Vec::new();let mut directions=Vec::new();
        for y in 0..grid {for x in 0..grid {
            let east=-extent+x as f64*step;let north=-extent+y as f64*step;
            let h=if fixed_slant_range {(height_km*height_km-east*east-north*north).sqrt()}else{height_km};
            let (r,d)=altitude_pixel(east,north,h).map_err(err)?;
            ranges.push(r);directions.push(d);
        }}
        let mut beams=Vec::new();
        for b in 0..5 {
            let az=g.beam_azimuth_deg[b].to_radians();let el=g.beam_elevation_deg[b].to_radians();
            let center=[el.cos()*az.sin(),el.cos()*az.cos(),-el.sin()];
            let radius=9f64.to_radians();let mut pixels=Vec::new();let mut weights=Vec::new();
            for (p,d) in directions.iter().enumerate() {
                let angle=(0..3).map(|j|center[j]*d[j]).sum::<f64>().clamp(-1.,1.).acos();
                if angle<radius {pixels.push(p);weights.push((std::f64::consts::PI*angle/(2.*radius)).cos().powi(2));}
            }
            let dirs:Vec<_>=pixels.iter().map(|p|directions[*p]).collect();
            let imager=TsvdImager::new(&g.positions_m,&dirs,g.frequency_hz,18,1.).map_err(err)?;
            beams.push(Beam {pixels,imager,weights});
        }
        Ok(Self{beams,ranges,pixels:grid*grid,step,nyquist:0.})
    }

    /// Spectra are complex64 [look,7,cell], transported without phase changes.
    /// Output [pixel,4] = relative dB power, circular velocity, width, support.
    pub fn image(&mut self, header_json:&str, geometry_json:&str, values:&[f32])->Result<Vec<f32>,JsValue> {
        let p:Packet=serde_json::from_str(header_json).map_err(err)?;
        let g:Geometry=serde_json::from_str(geometry_json).map_err(err)?;
        let nc=p.index.len();let nf=p.velocity_mps.len();let nr=p.ranges_km.len();
        if nf<2 || nr<2 || p.looks==0 || p.looks>256 || nc>200_000 || values.len()!=2*p.looks*7*nc || !p.nyquist_mps.is_finite() || p.nyquist_mps<=0. {
            return Err(err("Invalid sparse image packet"));
        }
        if p.noise_power.len()!=7 || p.noise_power.iter().any(|a|a.len()!=5 || a.iter().any(|b|b.len()!=nr || b.iter().any(|n| !n.is_finite() || *n<=0.)))
            || p.ranges_km.windows(2).any(|r| !r[0].is_finite() || r[1]<=r[0]) {
            return Err(err("Invalid dense noise/range metadata"));
        }
        self.nyquist=p.nyquist_mps;
        let pairs=channel_pairs();let mut spectra=vec![0.;self.pixels*nf];let mut support=vec![0.;self.pixels];
        let dr=p.ranges_km[1]-p.ranges_km[0];
        for (b,beam) in self.beams.iter().enumerate() {
            for (j,&pixel) in beam.pixels.iter().enumerate() {support[pixel]+=beam.weights[j];}
            let mut channel=Vec::with_capacity(p.looks*7);
            for (cell,&[beam_id,f,r]) in p.index.iter().enumerate() {
                if beam_id>=5 || f>=nf || r>=nr {return Err(err("Sparse index outside axes"));}
                if beam_id!=b || !beam.pixels.iter().any(|&pixel|(self.ranges[pixel]-p.ranges_km[r]).abs()<dr) {continue;}
                channel.clear();
                for look in 0..p.looks {for ch in 0..7 {
                    let offset=2*((look*7+ch)*nc+cell);
                    channel.push(Complex64::new(values[offset] as f64,values[offset+1] as f64));
                }}
                let mut cov=calibrated_covariance(&channel,&g.phasecal_rad[b],1).map_err(err)?;
                let auto:Vec<f64>=cov[..7].iter().map(|z|z.re).collect();
                let noise:Vec<f64>=(0..7).map(|ch|p.noise_power[ch][b][r]).collect();
                let excess=(auto.iter().sum::<f64>()/noise.iter().sum::<f64>()-1.).max(0.);
                let signal:Vec<_>=(0..7).map(|ch|(auto[ch]-noise[ch]).max(noise[ch])).collect();
                for (z,&(a,c)) in cov.iter_mut().zip(&pairs) {*z/=(signal[a]*signal[c]).sqrt();}
                let power=beam.imager.power(&cov).map_err(err)?;
                let norm=power.iter().sum::<f64>()*self.step*self.step;
                if norm<=0. {continue;}
                for (j,&pixel) in beam.pixels.iter().enumerate() {
                    let fraction=(1.-(self.ranges[pixel]-p.ranges_km[r]).abs()/dr).max(0.);
                    spectra[pixel*nf+f]+=power[j]/norm*excess*fraction*beam.weights[j];
                }
            }
        }
        let mut total=vec![0.;self.pixels];let mut moment=vec![(0.,0.);self.pixels];
        for pixel in 0..self.pixels {
            if support[pixel]>0. {
                for power in &mut spectra[pixel*nf..(pixel+1)*nf] {*power/=support[pixel];}
                if let Some((sum,v,w))=doppler_moments(&spectra[pixel*nf..(pixel+1)*nf],&p.velocity_mps,p.nyquist_mps).map_err(err)? {
                    total[pixel]=sum;moment[pixel]=(v,w);
                }
            }
        }
        let peak=total.iter().copied().fold(0.,f64::max);let mut out=vec![0f32;self.pixels*4];
        for pixel in 0..self.pixels {
            out[4*pixel]=if support[pixel]==0. {f32::NAN} else if total[pixel]==0. {-120.} else {(10.*(total[pixel]/peak).log10()) as f32};
            out[4*pixel+1]=moment[pixel].0 as f32;out[4*pixel+2]=moment[pixel].1 as f32;
            out[4*pixel+3]=if total[pixel]>0. {1.}else{0.};
        }
        Ok(out)
    }
}
