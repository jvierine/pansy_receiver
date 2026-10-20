"""Refresh sky-noise diagnostics on the analysis host, independently of the GUI.

Read dense noise summaries, never complex spectra. Reuse the precomputed sky
convolution; fit hourly, render native SVGs, publish their status last.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import time

import h5py
import numpy as np

from .publish_rvd import atomic_json
from .sky_noise import update


def fit_specification(path):
    with h5py.File(path,'r') as h:
        times=h['times_s'][()]
        sky=h['sky_temperature_k'][()]
        noise=h['noise_power'][()]
        fits=json.loads(h.attrs['fits_json'])
        generated=float(h.attrs['generated_unix_s'])
        fitted=float(h.attrs['fit_unix_s'])
    if not len(times) or sky.shape!=(len(times),5) or noise.shape!=(len(times),5,7) or len(fits)!=7:
        raise ValueError('Invalid cached sky-noise axes or module fits')
    if not np.isfinite(noise).all() or np.any(noise<=0) or not np.isfinite(sky).all():
        raise ValueError('Invalid sky-noise measurements')
    data=dict(times_h=((times-times[0])/3600).tolist(),
              time_tick_labels=[datetime.fromtimestamp(t,timezone.utc).strftime('%m-%d %H:%M')
                                for t in (times[0],(times[0]+times[-1])/2,times[-1])],
              sky_temperature_k=sky.tolist(),noise_power=noise.tolist())
    return dict(noise_model_data=data),dict(schema='pansy_sky_noise_v1',noise_model_status=fits,
        noise_model_updated_unix_s=generated,noise_model_fit_unix_s=fitted,
        noise_model_observation_stop_unix_s=float(times[-1]))


def publish(root,cache,fit,output,destination=None,renderer=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    update(root,cache,fit)
    spec,status=fit_specification(fit)
    revision=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:16]
    specification=output/'plot-spec.json';atomic_json(specification,spec)
    renderer=Path(renderer) if renderer is not None else Path(__file__).parent/'plotting/target/release/pansy-pmse-plots'
    subprocess.run([str(renderer),'--noise-only',str(specification),str(output)],
                   check=True,capture_output=True,text=True,timeout=180)
    files={str(ch):{str(w):f'noise-module-{ch}-{w}.svg' for w in (640,960)} for ch in range(7)}
    if not all((output/name).is_file() for widths in files.values() for name in widths.values()):
        raise OSError('Native renderer did not produce all sky-noise plots')
    status['noise_plot']=dict(status='ready',revision=revision,files=files,
        generated_utc=datetime.now(timezone.utc).isoformat())
    atomic_json(output/'status.json',status)
    if destination:
        subprocess.run(['rsync','-rt','--delay-updates','--exclude=.*','--exclude=*.json','--exclude=*.partial',
                        str(output)+'/',destination],check=True,timeout=180)
        subprocess.run(['rsync','-rt','--delay-updates',str(output/'status.json'),destination],check=True,timeout=30)
    print(f'Sky noise ready: seven modules, observations through {status["noise_model_observation_stop_unix_s"]}',flush=True)
    return status


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','cache','fit','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--destination');p.add_argument('--loop',action='store_true');a=p.parse_args()
    while True:
        try:publish(a.root,a.cache,a.fit,a.output,a.destination)
        except Exception:
            import traceback;traceback.print_exc()
            if not a.loop:raise
        if not a.loop:break
        time.sleep(600)


if __name__=='__main__':main()
