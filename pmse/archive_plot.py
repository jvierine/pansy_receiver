"""Native storage diagnostics from inventory metadata, never voltage spectra."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
from urllib.request import urlopen

from .publish_rvd import atomic_json


def refresh_archive_plot(manifest, output, renderer=None):
    output=Path(output)
    archive=manifest.get('transfer',{}).get('archive',{})
    previous=manifest.get('archive_plot',{})
    files={str(width):f'axes/archive-growth-{width}.svg' for width in (640,960)}
    existing=all((output/name).is_file() for name in files.values())
    if not archive.get('hourly_unix_s'):
        manifest['archive_plot']=dict(status='unavailable')
        return
    revision=hashlib.sha256(json.dumps(archive,sort_keys=True).encode()).hexdigest()[:16]
    if existing and previous.get('status')=='ready' and previous.get('revision')==revision:
        return
    renderer=Path(renderer) if renderer is not None else Path(__file__).parent/'plotting/target/release/pansy-pmse-plots'
    try:
        specification=output/'archive-plot-spec.json'
        atomic_json(specification,dict(transfer=dict(archive=archive)))
        subprocess.run([str(renderer),'--archive-only',str(specification),str(output/'axes')],
                       check=True,capture_output=True,text=True,timeout=120)
        if not all((output/name).is_file() for name in files.values()):
            raise OSError('Native renderer did not produce both archive plots')
        manifest['archive_plot']=dict(status='ready',revision=revision,files=files,
            generated_utc=datetime.now(timezone.utc).isoformat())
    except (OSError,subprocess.SubprocessError) as error:
        # A diagnostic failure must not freeze current GUI metadata or hide
        # the last successful plots. Report their age separately from RTI axes.
        print(f'Archive plot refresh delayed: {error}',flush=True)
        fallback=dict(status='unavailable')
        if existing:
            generated=previous.get('generated_utc') or datetime.fromtimestamp(
                min((output/name).stat().st_mtime for name in files.values()),timezone.utc).isoformat()
            fallback=dict(status='delayed',revision=previous.get('revision','previous'),
                          generated_utc=generated,files=files)
        manifest['archive_plot']=fallback


def attach_archive_plot(manifest,output):
    """Web-host metadata only: use analysis-produced plots, never build them."""
    output=Path(output)
    try:
        status=json.loads((output/'archive-plot-status.json').read_text())
        expected={str(w):f'axes/archive-growth-{w}.svg' for w in (640,960)}
        if status.get('files')!=expected or not all((output/p).is_file() for p in expected.values()):
            raise ValueError('Archive diagnostic files have not arrived yet')
        manifest['archive_plot']=status
    except (OSError,ValueError):
        if 'archive_plot' in manifest:manifest['archive_plot']['status']='delayed'
        else:manifest['archive_plot']=dict(status='unavailable')


def publish_archive_plot(manifest_url,output,destination):
    """Render the delivered web inventory on the analysis host, then push it."""
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    with urlopen(manifest_url,timeout=30) as response:manifest=json.load(response)
    refresh_archive_plot(manifest,output)
    status=manifest['archive_plot']
    if status['status']!='ready':raise RuntimeError('Archive plot rendering is not ready')
    atomic_json(output/'archive-plot-status.json',status)
    subprocess.run(['rsync','-rt','--delay-updates','--exclude=.*','--exclude=*.json','--exclude=*.partial',
                    str(output)+'/',destination],check=True,timeout=120)
    subprocess.run(['rsync','-rt','--delay-updates',str(output/'archive-plot-status.json'),destination],check=True,timeout=30)
