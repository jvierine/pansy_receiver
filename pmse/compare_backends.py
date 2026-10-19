"""Validate independently produced native records; persist numerical evidence."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from .sparse_record import read_complete,channel_covariance


def compare(reference,candidate):
    exact=['sparse_index','range_sample','altitude_valid','look_bounds_sample',
           'pulses_per_beam_per_look','tx_cycle_start_sample','tx_cycle_verified','rx_cycle_used']
    for name in exact: np.testing.assert_array_equal(reference[name],candidate[name],err_msg=name)
    for name in ['slant_range_km','nominal_altitude_km','radial_velocity_mps','doppler_frequency_hz',
                 'coverage_seconds','mode_duration_seconds']:
        np.testing.assert_allclose(reference[name],candidate[name],rtol=0,atol=1e-10,err_msg=name)
    errors={}
    for name in ['noise_power','noise_power_standard_error','mean_power','power_variance',
                 'power_standard_error','detector_snr_db','spectra']:
        a=np.asarray(reference[name]); b=np.asarray(candidate[name]); scale=max(float(np.max(abs(a),initial=0)),1e-20)
        np.testing.assert_allclose(a,b,rtol=2e-4,atol=scale*1e-5,err_msg=name)
        errors[name+'_relative_l2']=float(np.linalg.norm((b-a).ravel())/max(np.linalg.norm(a.ravel()),1e-20))
    a=channel_covariance(reference); b=channel_covariance(candidate)
    errors['covariance_relative_l2']=float(np.linalg.norm((b-a).ravel())/max(np.linalg.norm(a.ravel()),1e-20))
    errors.update(passed=True,cells=len(candidate['sparse_index']),looks=int(candidate['n_coherent_looks']),
        coverage_seconds=float(candidate['coverage_seconds']),
        reference_wall_seconds=float(reference['processing_wall_seconds']),
        rust_wall_seconds=float(candidate['processing_wall_seconds']))
    errors['speedup']=errors['reference_wall_seconds']/errors['rust_wall_seconds']
    return errors


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True); p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--key',type=int,required=True); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); reference=read_complete(a.reference,a.key); candidate=read_complete(a.candidate,a.key)
    evidence=compare(reference,candidate); a.output.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(a.output,'w') as h:
        h.attrs.update(reference=str(a.reference),candidate=str(a.candidate),key=a.key,
            reference_commit=reference['source_commit'],candidate_commit=candidate['source_commit'],
            generator='pmse/compare_backends.py',tolerance='rtol 2e-4, atol 1e-5 global field peak')
        for k,v in evidence.items(): h[k]=v
    print(json.dumps(evidence,sort_keys=True))
