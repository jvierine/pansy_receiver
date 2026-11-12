"""Continuous single-pulse Doppler, conditional on the selected range gate.

The 512-bin detector remains unchanged. This fit evaluates the original complex
RX samples against conjugated measured TX, without boxcar decimation or extra
zero padding. Each receiver gets its own nuisance complex amplitude; frequency
is shared. It does not pretend that finer sampling increases coherent dwell.
"""
import numpy as np
from scipy.optimize import minimize_scalar


def noise_variance(samples):
    samples = np.asarray(samples, dtype=np.complex128)
    center = np.median(samples.real, axis=-1) + 1j*np.median(samples.imag, axis=-1)
    return np.maximum(np.median(abs(samples-center[..., None])**2, axis=-1)/np.log(2), 1e-20)


def refine_doppler(tx, rx, gate, coarse_mps, *, sample_rate_hz=1e6,
                   radar_frequency_hz=47e6, fft_bin_hz=125000/512):
    """Profile Gaussian likelihood with independent complex RX amplitudes.

    Sigma is a local, fixed-range, single-tone model standard error, inflated
    for excess residual power and a common noisy-TX reference floor. It is NOT
    a calibrated total velocity error (multipath, range errors and acceleration
    remain possible). Invalid/boundary fits return NaNs and valid=False; callers
    retain the coarse observation separately, rather than disguising fallback.
    """
    tx = np.asarray(tx, dtype=np.complex128)
    rx = np.asarray(rx, dtype=np.complex128)
    bad = dict(velocity_mps=np.nan, sigma_mps=np.nan, valid=False,
               coherence=np.nan, reduced_chi2=np.nan, reference_noise_variance=np.nan)
    if tx.ndim != 1 or rx.ndim != 2 or len(tx)<8 or gate<0 or gate+len(tx)>rx.shape[1]:
        return bad
    if not np.isfinite(tx).all() or not np.isfinite(rx).all() or not np.isfinite(coarse_mps):
        return bad
    power = abs(tx)**2
    energy = power.sum()
    if energy<=0 or sample_rate_hz<=0 or radar_frequency_hz<=0 or fft_bin_hz<=0:
        return bad
    quiet = np.ones(rx.shape[1], dtype=bool)
    quiet[max(0,gate-4):min(rx.shape[1],gate+len(tx)+4)] = False
    if np.count_nonzero(quiet)<16:
        return bad
    variance = noise_variance(rx[:, quiet])
    y = rx[:, gate:gate+len(tx)]
    t = np.arange(len(tx))/sample_rate_hz
    t -= np.sum(power*t)/energy
    moment = np.sum(power*t*t)
    if moment<=0:
        return bad
    decoded = y*tx.conj()[None,:]
    def amplitudes(f):
        return (decoded @ np.exp(-2j*np.pi*f*t))/energy
    def score(f):
        return float(np.sum(abs(amplitudes(f))**2/variance)*energy)
    conversion = 299792458.0/(2*radar_frequency_hz)
    center = coarse_mps/conversion
    half_width = max(3*fft_bin_hz, 750.0)
    grid = center+np.linspace(-half_width, half_width, 17)
    values = np.array([score(f) for f in grid])
    k = int(np.argmax(values))
    if k==0 or k==len(grid)-1:
        return bad
    opt = minimize_scalar(lambda offset:-score(grid[k]+offset),
                          bounds=(grid[k-1]-grid[k],grid[k+1]-grid[k]),
                          method='bounded',options={'xatol':.01,'maxiter':32})
    frequency = grid[k]+opt.x
    amplitude = amplitudes(frequency)
    information = np.sum(abs(amplitude)**2/variance)
    if not opt.success or not np.isfinite(information) or information<=0:
        return bad
    model = amplitude[:,None]*tx[None,:]*np.exp(2j*np.pi*frequency*t)[None,:]
    dof = max(1, 2*rx.shape[0]*(len(tx)-1)-1)
    reduced_chi2 = float(2*np.sum(abs(y-model)**2/variance[:,None])/dof)
    reference_quiet = power < .001*power.max()
    ref_noise = float(noise_variance(tx[reference_quiet])) if np.count_nonzero(reference_quiet)>=8 else 0.0
    variance_hz = (max(1.0,reduced_chi2)/information+ref_noise)/(8*np.pi**2*moment)
    coherence = float(np.sum(abs(amplitude)**2)*energy/max(float(np.sum(abs(y)**2)),1e-30))
    return dict(velocity_mps=float(frequency*conversion), sigma_mps=float(np.sqrt(variance_hz)*conversion),
                valid=True, coherence=coherence, reduced_chi2=reduced_chi2,
                reference_noise_variance=ref_noise)


def doppler_outliers(t, velocity, sigma):
    """Independent smooth Doppler QC; do not force it to the position solution.

    An echo's Doppler mask is separate from the trajectory mask. Quadratic time
    dependence allows ordinary changing radial projection/deceleration. Tiny
    formal errors cannot make the clipping tolerance arbitrarily small.
    """
    import itertools
    t,velocity,sigma = map(lambda x:np.asarray(x,dtype=float),(t,velocity,sigma))
    valid = np.isfinite(t)&np.isfinite(velocity)&np.isfinite(sigma)&(sigma>0)
    if np.count_nonzero(valid)<3:
        return np.zeros(len(t),dtype=bool),np.full(len(t),np.nan)
    indices=np.flatnonzero(valid);span=np.ptp(t[valid])
    if span<=0: return np.zeros(len(t),dtype=bool),np.full(len(t),np.nan)
    x=(t-np.mean(t[valid]))/span
    degree=2 if len(indices)>=8 else 1
    design=np.vander(x,N=degree+1,increasing=True)
    noise=np.maximum(sigma,50.0)
    # Fixed measurement scales: broad residual MAD must not turn an incoherent
    # two-echo mixture into an apparently good track by inflating its own noise.
    if len(indices)<=12:
        candidates=list(itertools.combinations(indices,degree+1))
    else:
        rng=np.random.default_rng(0)
        candidates=[rng.choice(indices,degree+1,replace=False) for _ in range(256)]
    best=None;best_score=(-1,-np.inf)
    for subset in candidates:
        subset=np.asarray(subset)
        if np.ptp(x[subset])<.25: continue
        coefficient,_,rank,_=np.linalg.lstsq(design[subset],velocity[subset],rcond=None)
        if rank<degree+1: continue
        model=design@coefficient
        residual=np.where(valid,(velocity-model)/noise,np.inf)
        keep=valid&(abs(residual)<=4.5)
        score=(int(keep.sum()),-float(np.sum(np.minimum(residual[valid]**2,4.5**2))))
        if score>best_score:best_score=score;best=keep
    minimum=max(3,degree+2 if degree==2 else 3,int(np.ceil(.6*len(indices))))
    if best is None or best.sum()<minimum:
        return np.zeros(len(t),dtype=bool),np.full(len(t),np.nan)
    keep=best
    for _ in range(8):
        coefficient,_,rank,_=np.linalg.lstsq(design[keep]/noise[keep,None],velocity[keep]/noise[keep],rcond=None)
        if rank<degree+1:break
        model=design@coefficient
        new=valid&(abs((velocity-model)/noise)<=4.5)
        if new.sum()<minimum:break
        stable=np.array_equal(new,keep);keep=new
        if stable:return keep,model
    # Nonconvergence or inadequate consensus is a failed QC result, not a
    # justification for expanding the acceptance threshold.
    return np.zeros(len(t),dtype=bool),np.full(len(t),np.nan)
