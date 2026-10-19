"""Read-only, ephemeral browser transport and native Rustplotlib axis service.

No radar imaging, covariance calculation, or new scientific products on server.
"""
import argparse
import hashlib
import json
import struct
import subprocess
import threading
from datetime import datetime,timezone,timedelta
from collections import OrderedDict
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse, parse_qs

import h5py
import numpy as np

HEIGHT_KM=85.290954301
_render_lock=threading.Lock()
_frames=OrderedDict()


def day_packet(root,day,manifest):
    """One UTC day from existing native metadata; never read voltage spectra."""
    first=datetime.strptime(day,'%Y-%m-%d').replace(tzinfo=timezone.utc)
    if first.strftime('%Y-%m-%d')!=day:raise ValueError('Invalid UTC day')
    start=round(first.timestamp()*1e6);stop=start+86_400_000_000
    packets=[]
    dirs=list(root.glob(day+'T*'))+list(root.glob((first-timedelta(days=1)).strftime('%Y-%m-%d')+'T23*'))
    for directory in sorted(dirs):
        for path in sorted(directory.glob('pmse@*.h5')):
            with h5py.File(path,'r') as h:
                for key in sorted(h):
                    if not key.isdigit():continue
                    g=h[key]
                    if 'complete' not in g or int(g['complete'][()])!=1:continue
                    begin=int(g['mode_start_sample'][()]);end=int(g['mode_stop_sample_exclusive'][()])
                    if end<=start or begin>=stop:continue
                    schema=g['schema_version'][()].decode()
                    if schema!='pansy_pmse_sparse16_v2':raise ValueError('Refusing mixed/legacy sparse product')
                    packets.append(dict(key=key,start_us=begin,stop_us=end,
                        duration_s=float(g['mode_duration_seconds'][()]),coverage_s=float(g['coverage_seconds'][()]),
                        looks=int(g['n_coherent_looks'][()]),quality=g['quality'][()].decode(),
                        product_schema=schema,index=g['sparse_index'][()].tolist(),snr_db=g['detector_snr_db'][()].tolist(),
                        look_dominance=None,look_persistence=None))
    packets.sort(key=lambda p:p['start_us'])
    revision=hashlib.sha256(('16ms:'+day+','.join(p['key'] for p in packets)).encode()).hexdigest()[:16]
    return dict(manifest,selected_day=day,revision=revision,axis_start_us=start,axis_stop_us=stop,
                frame_packets=packets,frame_count=len(packets),
                frames=[{k:p[k] for k in ('key','start_us','stop_us','duration_s','coverage_s','looks','quality')} for p in packets])


def geometry(reference):
    groups={}
    for line in (reference/'cfg/antpos.csv').read_text().splitlines()[1:]:
        fields=line.split(',')
        groups.setdefault(fields[1],[]).append([float(v) for v in fields[4:7]])
    connections=[line.split(',')[0] for line in (reference/'cfg/connections.txt').read_text().splitlines()]
    positions=[np.mean(groups[name],axis=0).tolist() for name in connections if name!='RFTX']
    phases=[];hashes=[]
    for beam in range(5):
        path=reference/f'data/phases{beam}.h5'
        with h5py.File(path) as h:phases.append(h['phasecal'][()].tolist())
        hashes.append(hashlib.sha256(path.read_bytes()).hexdigest())
    return dict(positions_m=positions,phasecal_rad=phases,phasecal_sha256=hashes,
                beam_azimuth_deg=[0,0,90,180,270],beam_elevation_deg=[90,80,80,80,80],
                frequency_hz=47e6,phase_sign=1,calibration_status='historical; contemporary calibration not independently validated')


def existing_packet(root, key, range_km, offset_s=0):
    """Select bytes from existing data only. Preserve all looks/channels/bins.

    Cap incoherent look selection to 120 s, aligned to the mesomode start.
    Original coherent looks are indivisible and never exceed the window.
    """
    import datetime
    hour=datetime.datetime.fromtimestamp(key/1e6,datetime.timezone.utc).strftime('%Y-%m-%dT%H-00-00')
    path=root/hour/f'pmse@{key//1_000_000//600*600}.h5'
    with h5py.File(path,'r') as h:
        g=h[str(key)]
        if int(g['complete'][()])!=1:raise ValueError('Record incomplete')
        index=g['sparse_index'][()];ranges=g['slant_range_km'][()]
        gate=int(np.argmin(abs(ranges-range_km)))
        cells=np.flatnonzero(index[:,2]==gate)
        bounds=g['look_bounds_sample'][()] if 'look_bounds_sample' in g else None
        starts=bounds[:,0] if bounds is not None else None
        if starts is None and float(g['mode_duration_seconds'][()])>120:
            raise ValueError('Long record lacks look times; refusing unbounded integration')
        start=int(g['mode_start_sample'][()]);stop=int(g['mode_stop_sample_exclusive'][()])
        begin=start+int(offset_s*1e6);end=min(stop,begin+120_000_000)
        lo=0;hi=g['spectra'].shape[0]
        if starts is not None:
            use=np.flatnonzero((starts>=begin)&(bounds[:,1]<=end))
            if not len(use):raise ValueError('No coherent looks in selected integration')
            lo=int(use[0]);hi=int(use[-1])+1
        z=g['spectra'][lo:hi,:,cells]
        if z.shape[0]==0:raise ValueError('No coherent looks in selected integration')
        velocity=g['radial_velocity_mps'][()]
        config=json.loads(g['config_json'][()].decode())
        native_ipp=5*config['ipp_samples']/config['sample_rate_hz']*(2 if config.get('complementary_pairs') else 1)
        nyquist=299792458./(4*config['radar_frequency_hz']*native_ipp)
        header=dict(key=str(key),start_us=begin,stop_us=end,range_gate=gate,selected_range_km=float(ranges[gate]),
                    index=index[cells].tolist(),ranges_km=ranges.tolist(),velocity_mps=velocity.tolist(),
                    noise_power=g['noise_power'][()].tolist(),looks=z.shape[0],nyquist_mps=nyquist,
                    spectrum_layout='look,channel,cell; little-endian complex64',
                    schema=g['schema_version'][()].decode(),censored=True)
        raw=np.ascontiguousarray(z,dtype='<c8').tobytes()
    encoded=json.dumps(header,separators=(',',':'),allow_nan=False).encode()
    padding=(-len(encoded))%4
    return struct.pack('<I',len(encoded))+encoded+b' '*padding+raw


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/mnt/data/juha/pansy/metadata/pmse_sparse16'))
    p.add_argument('--manifest',type=Path,default=Path('/mnt/data/juha/pansy/rvd-web/data/manifest.json'))
    p.add_argument('--port',type=int,default=18769)
    a=p.parse_args();reference=Path(__file__).parent/'reference';geo=geometry(reference)
    renderer=Path(__file__).parent/'plotting/target/release/pansy-pmse-plots'
    class Handler(BaseHTTPRequestHandler):
        def send(self,data,kind='application/json',status=200):
            self.send_response(status);self.send_header('Content-Type',kind);self.send_header('Content-Length',str(len(data)))
            self.send_header('Cache-Control','no-store');self.end_headers();self.wfile.write(data)
        def do_GET(self):
            try:
                url=urlparse(self.path)
                if url.path.endswith('/geometry'):
                    return self.send(json.dumps(geo).encode())
                if url.path.endswith('/day'):
                    q=parse_qs(url.query);m=json.loads(a.manifest.read_text())
                    return self.send(json.dumps(day_packet(a.root,q['day'][0],m),separators=(',',':'),allow_nan=False).encode())
                if url.path.endswith('/spectra'):
                    q=parse_qs(url.query);key=int(q['key'][0]);r=float(q['range_km'][0]);offset=float(q.get('offset_s',['0'])[0])
                    if not 50<=r<=123 or not 0<=offset<=86400 or key<=0:raise ValueError('Invalid selection')
                    return self.send(existing_packet(a.root,key,r,offset),'application/octet-stream')
                self.send(b'{}',status=404)
            except (ValueError,KeyError,OSError) as e:self.send(json.dumps(dict(error=str(e))).encode(),status=400)
        def do_POST(self):
            try:
                if not self.path.endswith('/frame'):return self.send(b'{}',status=404)
                length=int(self.headers.get('Content-Length','0'))
                if not 0<length<=4096:raise ValueError('Invalid request size')
                body=json.loads(self.rfile.read(length));m=json.loads(a.manifest.read_text());ident=body['id']
                day=datetime.strptime(str(int(body.get('day',int(m.get('selected_day','2026-10-19').replace('-',''))))),'%Y%m%d').replace(tzinfo=timezone.utc)
                origin=day.timestamp();duration=24.
                spec=dict(bounds=[0,50,max(duration,.001),120],xlabel='Time (UTC)',ylabel='Nominal altitude (km)',
                          colorbar_label='SNR (dB)',panels=['Zenith','North','East','South','West'],utc_origin_seconds=origin)
                if ident in ('rd','td'):
                    v=m['geometry']['velocity_mps'];v0,v1=m['velocity_limits_mps']
                    if ident=='rd':
                        spec.update(bounds=[v0,50,v1,120],xlabel=r'Radial velocity (m s$^{-1}$)',panels=['Spectrum'],utc_origin_seconds=None)
                    else:
                        spec.update(bounds=[0,v0,max(duration,.001),v1],ylabel=r'Radial velocity (m s$^{-1}$)')
                elif ident=='image':
                    radius=float(body.get('range_km',HEIGHT_KM))
                    if not 50<=radius<=123:raise ValueError('Invalid range')
                    extent=radius*np.sin(np.deg2rad(20))
                    spec.update(bounds=[-extent,-extent,extent,extent],xlabel='East (km)',ylabel='North (km)',
                                colorbar_label='Relative power (dB)',panels=[''],utc_origin_seconds=None,equal_aspect=True)
                elif ident!='rti':raise ValueError('Unknown plot')
                request={k:body[k] for k in ('width','viewport','color_limits')}
                if ident in ('image','rti') and body.get('hsv')==1:
                    velocity=float(body.get('velocity_limit',30));width=float(body.get('width_limit',10))
                    if not 0<velocity<=1000 or not 0<width<=1000:raise ValueError('Invalid HSV limits')
                    spec['hsv_key']=dict(velocity_limit=velocity,width_limit=width,
                        power_limits=request['color_limits'],power_label='Power (dB)' if ident=='image' else 'SNR (dB)')
                cachekey=json.dumps([spec,request],sort_keys=True)
                with _render_lock:
                    if cachekey not in _frames:
                        result=subprocess.run([str(renderer),'--frame'],input=json.dumps(dict(spec=spec,request=request)),
                                              text=True,capture_output=True,check=True,timeout=30)
                        _frames[cachekey]=result.stdout.encode()
                        while len(_frames)>128:_frames.popitem(last=False)
                    data=_frames[cachekey]
                self.send(data)
            except (ValueError,KeyError,OSError,subprocess.SubprocessError) as e:
                self.send(json.dumps(dict(error=str(e))).encode(),status=400)
    ThreadingHTTPServer(('127.0.0.1',a.port),Handler).serve_forever()


if __name__=='__main__':main()
