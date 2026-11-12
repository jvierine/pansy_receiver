#!/usr/bin/env python3
"""Bounded offline cut replay: never write production metadata or telemetry.

    python receiver/audit_simple_fit_precision.py --output /tmp/check.h5 CUT.h5...

Radiant calculation is stubbed because this checks measurement extraction,
trajectory/Doppler QC, and the actual Digital Metadata round trip, not orbit
astrometry. Only scientific measurement fields are copied to the audit HDF5.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import h5py
import numpy as np
import process_cut_meteor as processor
from meteor_metadata import StrictDigitalMetadataReader
from meteor_fit import simple_fit_meteor


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cuts',nargs='+',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--max-events',type=int,default=8)
    args=parser.parse_args()
    if not 1<=args.max_events<=32: parser.error('max-events must be 1..32')
    if args.output.exists(): parser.error('refusing to overwrite an existing audit')
    with h5py.File(args.output,'x') as report,tempfile.TemporaryDirectory(prefix='pansy-fit-roundtrip-') as tmp:
        report.attrs['purpose']='offline Doppler/trajectory audit, not a radiant catalogue'
        report.attrs['astrometry_stubbed']=True
        root=Path(tmp)/'metadata';root.mkdir()
        underlying=processor.simple_fit_writer(str(root));count=0
        class RecordingWriter:
            def write(self,samples,payload):
                self.last=int(np.atleast_1d(samples)[0])
                return underlying.write(samples,payload)
        writer=RecordingWriter()
        for path in args.cuts:
            before=hashlib.sha256(path.read_bytes()).hexdigest()
            with h5py.File(path,'r') as stream:
                for key in sorted(stream,key=int):
                    if count>=args.max_events:break
                    data={k:stream[key][k][()] for k in stream[key]}
                    if np.max(data['c_snr'])<20: continue
                    start=time.monotonic()
                    with patch.object(processor.pc,'simple_fit_metadata_dir',str(root)), \
                         patch.object(processor,'SIMPLE_FIT_WRITE_LOCK',str(Path(tmp)/'writer.lock')), \
                         patch.object(processor.sr,'get_radiant',return_value={'lat':0.,'lon':0.}), \
                         patch('monitor.event.save_event'):
                        accepted=processor.process_cut(data,writer,write_dm=True)
                    metrics={'input_key':int(key),'accepted':bool(accepted),'wall_seconds':time.monotonic()-start}
                    if accepted:
                        reader=StrictDigitalMetadataReader(str(root));lo,hi=reader.get_bounds()
                        records=reader.read(lo,hi)
                        latest=records[writer.last]
                        g=report.create_group(str(int(key)))
                        for field,value in latest.items():
                            if field not in ('eclat','eclon'):g[field]=value
                        assert 'doppler_mps' in reader.get_fields(),'Doppler missing from field discovery'
                        for field in ('doppler_mps','doppler_inlier','fit_inlier','doppler_std_mps'):
                            assert len(latest[field])==len(latest['txidx'])
                        t=np.asarray(latest['txidx'])/1e6;t-=t.mean()
                        keep=np.asarray(latest['doppler_inlier'],dtype=bool)
                        def rms(y):
                            if keep.sum()<4:return None
                            residual=y[keep]-np.polyval(np.polyfit(t[keep],y[keep],2),t[keep])
                            return float(np.sqrt(np.mean(residual**2)))
                        old=simple_fit_meteor(t,latest['ew'],latest['ns'],latest['up'])
                        metrics.update(pulses=len(t),position_rejected=int(np.count_nonzero(~np.asarray(latest['fit_inlier']))),
                            doppler_kept=int(keep.sum()),fft_bin_mps=float(latest['doppler_fft_bin_width_mps']),
                            coarse_smooth_rms_mps=rms(np.asarray(latest['doppler_coarse_mps'])),
                            refined_smooth_rms_mps=rms(np.asarray(latest['doppler_mps'])),
                            median_formal_sigma_mps=float(np.nanmedian(latest['doppler_std_mps'])),
                            old_position_rms_m=np.asarray(old[3:]).tolist(),new_position_rms_m=(1000*np.asarray(latest['std'])).tolist(),
                            velocity_change_mps=float(1000*np.linalg.norm(latest['v0']-old[1])))
                        metrics['old_position_rms_m']=(1000*np.asarray(metrics['old_position_rms_m'])).tolist()
                        g.attrs['input_file']=str(path);g.attrs['input_sha256']=before
                        g.attrs['metrics_json']=json.dumps(metrics)
                    print(json.dumps(metrics),flush=True);count+=1
            assert hashlib.sha256(path.read_bytes()).hexdigest()==before,'Input cut was modified'
        report.attrs['events_evaluated']=count


if __name__=='__main__':main()
