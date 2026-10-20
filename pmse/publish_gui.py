"""Prepare lightweight GUI metadata locally; never load or process complex spectra.

Syowa owns radar processing. Browser workers own imaging and spectral moments.
The web host indexes only changed native files into atomic local HTTP packets.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import gzip
import json
from pathlib import Path
import time

import h5py

from .publish_rvd import atomic_json
from .serve_browser import day_packet
from .archive_inventory import inventory
from .web_retention import sha256
from .archive_plot import refresh_archive_plot


class MetadataCache:
    """Persistent per-file metadata; replacements invalidate even at equal size."""
    def __init__(self, directory, source=None):
        self.directory=Path(directory)
        self.directory.mkdir(parents=True,exist_ok=True)
        self.pending={}
        self.source=Path(source) if source else None
        self.hashes={}

    def source_directories(self,root,patterns):
        return [root/d.name for pattern in patterns for d in self.source.glob(pattern)] if self.source else []

    def source_files(self,directory):
        return [directory/p.name[:-8] for p in (self.source/directory.name).glob('pmse@*.h5.json.gz')] if self.source else []

    def location(self, path):
        return self.directory/(hashlib.sha256(str(path).encode()).hexdigest()+'.json')

    @staticmethod
    def signature(path):
        st=path.stat()
        return [st.st_ino,st.st_mtime_ns,st.st_size]

    def get(self, path):
        if self.source:
            source=self.source/path.parent.name/(path.name+'.json.gz')
            if source.exists():
                hint=json.loads(gzip.decompress(source.read_bytes()))
                ready=False
                if path.is_file() and path.stat().st_size==hint['native_bytes']:
                    state=self.signature(path)
                    previous=self.hashes.get(str(path))
                    if previous is None or previous[0]!=state:
                        previous=(state,sha256(path));self.hashes[str(path)]=previous
                    ready=previous[1]==hint['native_sha256']
                return [dict(p,spectra_ready=ready) for p in hint['packets']]
        state=self.signature(path);self.pending[str(path)]=state
        try:
            cached=json.loads(self.location(path).read_text())
            if cached['signature']==state:return cached['packets']
        except (OSError,ValueError,KeyError):pass
        return None

    def put(self, path, packets):
        state=self.pending[str(path)]
        if state==self.signature(path):
            atomic_json(self.location(path),dict(signature=state,packets=packets))


def publish(root, output, cache):
    """Refresh changed metadata in the background, latest day first."""
    root=Path(root);output=Path(output)
    manifest_path=output/'manifest.json'
    base=json.loads(manifest_path.read_text())
    # Daily files contain packets themselves; do not nest yesterday's packets.
    base.pop('frame_packets',None)
    days={p.name[:10] for p in root.glob('20*T*') if any(p.glob('pmse@*.h5'))}
    if cache.source:days.update(p.name[:10] for p in cache.source.glob('20*T*') if any(p.glob('pmse@*.h5.json.gz')))
    days=sorted(days)
    if not days:return
    (output/'days').mkdir(exist_ok=True)
    base.update(available_days=days,generated_utc=datetime.now(timezone.utc).isoformat(),
                plot_generation_status='unavailable',gui_processing='metadata only; no complex spectra processing')
    status=root/'processing-status.json'
    if status.exists():base['source_processing']=json.loads(status.read_text())
    for day in reversed(days):
        try:
            packet=day_packet(root,day,base,cache)
        except (OSError,KeyError,ValueError) as error:
            print(f'GUI metadata not ready for {day}: {error}',flush=True)
            continue
        sidecar=output/f'gui-data-{day}.h5'
        packet['sidecar']=sidecar.name
        ready=False
        if sidecar.exists():
            with h5py.File(sidecar,'r') as h:ready=h.attrs.get('revision')==packet['revision']
        if not ready:
            temporary=sidecar.with_suffix('.h5.partial')
            with h5py.File(temporary,'w') as h:
                h.attrs.update(schema='pansy_rvd_plot_v2',revision=packet['revision'],
                               source=str(root),generator='pmse/publish_gui.py',
                               snr_definition='excess signal/noise',moment_censoring_db=3)
                for key,value in packet['geometry'].items():h[key]=value
                groups=h.create_group('frames')
                for item in packet['frame_packets']:
                    group=groups.create_group(item['key'])
                    for key,value in item.items():
                        if key=='key' or value is None:continue
                        if isinstance(value,str):group.attrs[key]=value
                        else:group[key]=value
            temporary.replace(sidecar)
        atomic_json(output/'days'/f'{day}.json',packet)
        if day==days[-1] and packet['frame_count']:
            # One compressed static fetch contains everything for the first plots.
            atomic_json(manifest_path,packet)
            print(f'GUI ready: {day}, {packet["frame_count"]} modes',flush=True)
    # Inventory reads dimensions and counts only; its per-file cache skips
    # unchanged chunks. It never loads the complex-valued spectra dataset.
    measured=inventory(root)
    current=json.loads(manifest_path.read_text())
    current.setdefault('transfer',{})['archive']=measured
    current['transfer']['archive'].update(destination=str(root),web_retention_days=7)
    refresh_archive_plot(current,output)
    atomic_json(manifest_path,current)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cache',type=Path,default=Path.home()/'.cache/pansy-rvd-metadata')
    parser.add_argument('--loop',action='store_true')
    args=parser.parse_args();cache=MetadataCache(args.cache,args.output/'source')
    while True:
        try:publish(args.root,args.output,cache)
        except Exception:
            import traceback;traceback.print_exc()
            if not args.loop:raise
        if not args.loop:break
        time.sleep(60)


if __name__=='__main__':main()
