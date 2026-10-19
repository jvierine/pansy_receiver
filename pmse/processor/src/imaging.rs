//! Native imaging kernels. I/O adapters must preserve missing-cell masks;
//! absence of a retained spectral bin is censoring, not measured zero power.
//! Geometry matches the historical positive steering/calibration convention.
use nalgebra::{Complex, DMatrix, DVector};
use std::f64::consts::PI;

type C = Complex<f64>;

/// Consecutive, half-open integration windows in sample coordinates.
/// Splitting never crosses a mesomode boundary and never drops the tail.
pub fn integration_windows(start: u64, stop: u64, cap_samples: u64)
    -> Result<Vec<(u64, u64)>, &'static str>
{
    if stop <= start || cap_samples == 0 { return Err("invalid integration bounds"); }
    let mut out = Vec::new();
    let mut first = start;
    while first < stop {
        let last = first.saturating_add(cap_samples).min(stop);
        out.push((first, last));
        first = last;
    }
    Ok(out)
}

/// Common horizontal plane above the radar, east/north/up coordinates.
/// The historical forward model uses negative vertical direction cosine.
pub fn altitude_pixel(east_km: f64, north_km: f64, height_km: f64)
    -> Result<(f64, [f64; 3]), &'static str>
{
    if ![east_km, north_km, height_km].iter().all(|v| v.is_finite()) || height_km <= 0. {
        return Err("invalid altitude-plane geometry");
    }
    let r = (east_km*east_km + north_km*north_km + height_km*height_km).sqrt();
    Ok((r, [east_km/r, north_km/r, -height_km/r]))
}

/// One covariance vector per retained range/Doppler cell; all seven channels
/// and all incoherent looks are required. Layout is [look, channel].
/// Phase sign is explicit: historical five-beam calibration uses +1.
pub fn calibrated_covariance(spectra: &[C], phases: &[f64; 7], phase_sign: i8)
    -> Result<Vec<C>, &'static str>
{
    if spectra.is_empty() || spectra.len()%7 != 0 || ![-1,1].contains(&phase_sign)
        || !phases.iter().all(|p| p.is_finite())
        || !spectra.iter().all(|s| s.re.is_finite() && s.im.is_finite()) {
        return Err("invalid spectra or calibration");
    }
    let pairs = channel_pairs();
    let mut result = vec![C::new(0.,0.); pairs.len()];
    for look in spectra.chunks_exact(7) {
        for (value, &(a,b)) in result.iter_mut().zip(&pairs) {
            *value += look[a]*look[b].conj();
        }
    }
    for (value, &(a,b)) in result.iter_mut().zip(&pairs) {
        *value *= C::from_polar(1./(spectra.len()/7) as f64,
            phase_sign as f64*(phases[a]-phases[b]));
    }
    Ok(result)
}

pub fn channel_pairs() -> Vec<(usize, usize)> {
    (0..7).map(|i| (i,i)).chain((0..7).flat_map(|i| ((i+1)..7).map(move |j| (i,j)))).collect()
}

/// Cache this operator per geometry; never run an SVD for every data frame.
pub struct TsvdImager {
    inverse: DMatrix<C>,
    pub singular_values: Vec<f64>,
    pub retained_modes: usize,
}

impl TsvdImager {
    /// Historical regularization is 1/(sigma+1). A zero ridge gives pure TSVD.
    pub fn new(positions_m: &[[f64; 3]; 7], directions: &[[f64; 3]],
        frequency_hz: f64, modes: usize, ridge: f64) -> Result<Self, &'static str>
    {
        if directions.is_empty() || modes == 0 || !frequency_hz.is_finite() || frequency_hz <= 0.
            || !ridge.is_finite() || ridge < 0.
            || !positions_m.iter().flatten().chain(directions.iter().flatten()).all(|x| x.is_finite())
            || directions.iter().any(|d| (d.iter().map(|v|v*v).sum::<f64>()-1.).abs()>1e-6) {
            return Err("invalid imaging geometry");
        }
        let pairs = channel_pairs();
        let k = 2.*PI*frequency_hz/299792458.;
        let forward = DMatrix::from_fn(pairs.len(), directions.len(), |row,col| {
            let (a,b)=pairs[row];
            let phase: f64 = (0..3).map(|j| (positions_m[b][j]-positions_m[a][j])*directions[col][j]).sum();
            C::from_polar(1., k*phase)
        });
        let svd = forward.svd(true,true);
        let s = svd.singular_values;
        let u = svd.u.ok_or("missing left singular vectors")?;
        let vh = svd.v_t.ok_or("missing right singular vectors")?;
        let keep = modes.min(s.len());
        let cutoff = s[0]*1e-10;
        let mut inverse = DMatrix::zeros(directions.len(), pairs.len());
        let mut actual = 0;
        for mode in 0..keep {
            if s[mode] <= cutoff { continue; }
            actual += 1;
            for p in 0..directions.len() {
                for row in 0..pairs.len() {
                    inverse[(p,row)] += vh[(mode,p)].conj()*u[(row,mode)].conj()/(s[mode]+ridge);
                }
            }
        }
        Ok(Self { inverse, singular_values: s.as_slice().to_vec(), retained_modes: actual })
    }

    /// Returns a relative reconstruction proxy, not calibrated brightness/SNR.
    pub fn power(&self, covariance: &[C]) -> Result<Vec<f64>, &'static str> {
        if covariance.len()!=self.inverse.ncols() ||
            !covariance.iter().all(|z| z.re.is_finite() && z.im.is_finite()) {
            return Err("invalid covariance");
        }
        Ok((&self.inverse*DVector::from_column_slice(covariance)).iter().map(|z|z.norm_sqr()).collect())
    }
}

/// Circular velocity moments across the full periodic Doppler interval.
/// Avoids a false zero-velocity centroid when power straddles Nyquist.
pub fn doppler_moments(power: &[f64], velocity: &[f64], nyquist: f64)
    -> Result<Option<(f64,f64,f64)>, &'static str>
{
    if power.len()!=velocity.len() || power.is_empty() || !nyquist.is_finite() || nyquist<=0.
        || power.iter().any(|p| !p.is_finite() || *p<0.)
        || velocity.iter().any(|v| !v.is_finite() || *v < -nyquist || *v >= nyquist) {
        return Err("invalid periodic Doppler spectrum");
    }
    let sum: f64 = power.iter().sum();
    if sum == 0. { return Ok(None); }
    let phasor: C = power.iter().zip(velocity).map(|(p,v)| C::from_polar(*p,PI*v/nyquist)).sum();
    let mean=phasor.arg()*nyquist/PI;
    let variance: f64 = power.iter().zip(velocity).map(|(p,v)| {
        let dv=(v-mean+nyquist).rem_euclid(2.*nyquist)-nyquist;
        p*dv*dv
    }).sum::<f64>()/sum;
    Ok(Some((sum,mean,variance.sqrt())))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cap_preserves_mode_and_tail() {
        assert_eq!(integration_windows(10,96,120).unwrap(), vec![(10,96)]);
        assert_eq!(integration_windows(10,311,120).unwrap(), vec![(10,130),(130,250),(250,311)]);
        assert!(integration_windows(10,10,120).is_err());
    }
    #[test]
    fn inclined_plane_range() {
        let (r,d)=altitude_pixel(15.,0.,85.290954301).unwrap();
        assert!(r>85.290954301);
        assert!((r*(-d[2])-85.290954301).abs()<1e-12);
    }
    #[test]
    fn covariance_calibration_sign() {
        let phase=[0.,0.2,-0.4,0.8,0.1,-0.3,0.6];
        let volt:Vec<C>=phase.iter().map(|&p|C::from_polar(1.,-p)).collect();
        for c in calibrated_covariance(&volt,&phase,1).unwrap() {
            assert!((c-C::new(1.,0.)).norm()<1e-12);
        }
    }
    #[test]
    fn tsvd_single_pixel_and_native_covariance() {
        let pos=std::array::from_fn(|i|[i as f64*4.,(i*i) as f64,0.]);
        let direction=[0.08,-0.06,-(1.-0.08f64.powi(2)-0.06f64.powi(2)).sqrt()];
        let k=2.*PI*47e6/299792458.;
        let volt: Vec<C>=pos.iter().map(|p| C::from_polar(1.,-k*(0..3).map(|j|p[j]*direction[j]).sum::<f64>())).collect();
        let cov=calibrated_covariance(&volt,&[0.;7],1).unwrap();
        let im=TsvdImager::new(&pos,&[direction],47e6,18,0.).unwrap();
        assert!((im.power(&cov).unwrap()[0]-1.).abs()<1e-10);
    }
    #[test]
    fn resolves_injected_direction_among_multiple_pixels() {
        let pos=std::array::from_fn(|i|[i as f64*9.,(i*i) as f64*3.,0.]);
        let directions:Vec<[f64;3]>=(-2..=2).map(|i| {
            let u=i as f64*0.07;
            [u,0.04,-(1.-u*u-0.04f64.powi(2)).sqrt()]
        }).collect();
        let k=2.*PI*47e6/299792458.;
        let target=3;
        let volt:Vec<C>=pos.iter().map(|p|C::from_polar(1.,
            -k*(0..3).map(|j|p[j]*directions[target][j]).sum::<f64>())).collect();
        let cov=calibrated_covariance(&volt,&[0.;7],1).unwrap();
        let im=TsvdImager::new(&pos,&directions,47e6,18,0.).unwrap();
        let power=im.power(&cov).unwrap();
        for (p,value) in power.iter().enumerate() {
            assert!((value-if p==target {1.} else {0.}).abs()<1e-9);
        }
    }
    #[test]
    fn wrap_moments_not_zero_velocity() {
        let (_,mean,width)=doppler_moments(&[1.,1.],&[-99.,99.],100.).unwrap().unwrap();
        assert!((mean.abs()-100.).abs()<1e-10);
        assert!((width-1.).abs()<1e-10);
        assert!(doppler_moments(&[0.],&[0.],100.).unwrap().is_none());
    }
}
