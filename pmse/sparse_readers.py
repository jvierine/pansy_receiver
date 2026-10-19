"""Explicit adapters, never a pretend dense XC array with zero-filled bins."""
import itertools
import numpy as np
from .sparse_record import channel_covariance


def cross_products(record,phasecal=None,amplitude=None,phase_sign=None):
    """Return selected-cell correlations in the historical pair order.

    C_cd=mean(S_c*conj(S_d)). The current five-beam imaging convention applies
    exp(+i*(phi_c-phi_d)); callers must explicitly request sign=+1 or -1 when
    supplying a phase calibration. This is not a claim of independently checked
    physical arrival-angle/vertical-velocity signs. No phasecal is applied by
    default. Missing sparse cells are NOT represented as measured zero power.
    """
    pairs=np.array([(i,i) for i in range(7)]+list(itertools.combinations(range(7),2)),dtype=np.int8)
    covariance=channel_covariance(record)
    values=covariance[:,pairs[:,0],pairs[:,1]]
    if amplitude is not None:
        amplitude=np.asarray(amplitude,dtype=float)
        if amplitude.shape!=(7,) or not np.isfinite(amplitude).all() or np.any(amplitude<=0):
            raise ValueError('Expected seven positive voltage gain factors')
        values*=amplitude[pairs[:,0]]*amplitude[pairs[:,1]]
    if phasecal is not None:
        if phase_sign not in (-1,1):
            raise ValueError('Explicitly specify phase_sign=+1 or -1 for this calibration')
        phasecal=np.asarray(phasecal,dtype=float)
        if phasecal.shape!=(5,7) or not np.isfinite(phasecal).all():
            raise ValueError('Expected finite beam-specific phase calibration [5,7]')
        beam=np.asarray(record['sparse_index'])[:,0]
        phi=phasecal[beam]
        values*=np.exp(1j*phase_sign*(phi[:,pairs[:,0]]-phi[:,pairs[:,1]]))
    return dict(sparse_index=np.asarray(record['sparse_index']).copy(),ch_pairs=pairs,
                cross_spectra=values,convention='mean(S_c*conj(S_d))',
                phasecal_applied=phasecal is not None,phase_sign=phase_sign,
                censored=True)
