"""Recolor PNG/GIF from an imaging HDF5 sidecar; no XC2 or inversion needed."""
import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
from PIL import Image


def colors(mean,width,db,cover,*,doppler_max,width_min,width_max,saturation_min,db_min,db_max):
    if not (doppler_max>0 and width_max>width_min>=0 and db_max>db_min
            and 0<=saturation_min<=1):
        raise ValueError('Invalid Doppler, width, saturation, or dB color limits')
    hsv=np.stack([.75*(1-np.clip((mean/doppler_max+1)/2,0,1)),
        1-(1-saturation_min)*np.clip((width-width_min)/(width_max-width_min),0,1),
        np.clip((db-db_min)/(db_max-db_min),0,1)],axis=-1)
    rgb=hsv_to_rgb(hsv)
    rgb[:,~cover]=0
    return rgb


def render(source,output,*,doppler_max=9.568,width_min=0,width_max=6.379,
           saturation_min=.2,db_min=-30,db_max=0,frame_ms=250,dpi=300):
    source,output=Path(source),Path(output)
    with h5py.File(source) as h:
        mean=h['centroid_mps'][()]; width=h['width_mps'][()]
        db=h['relative_power_db'][()]; cover=h['coverage_mask'][()]
        east=h['east_km'][()]; north=h['north_km'][()]
        hgt=float(h['height_km'][()]); times=h['time_unix_us'][()]
        az=np.deg2rad(h['beam_az_deg'][()]); el=np.deg2rad(h['beam_el_deg'][()])
        geometry=h.attrs['geometry']
    if geometry!='constant_altitude_plane':
        raise ValueError(f'Unsupported geometry: {geometry}; do not silently relabel')
    settings=dict(doppler_max=doppler_max,width_min=width_min,width_max=width_max,
                  saturation_min=saturation_min,db_min=db_min,db_max=db_max)
    rgb=colors(mean,width,db,cover,**settings)
    plt.rcParams.update({'font.size':11,'axes.titlesize':11,'axes.labelsize':11,
                         'xtick.labelsize':11,'ytick.labelsize':11})
    fig=plt.figure(figsize=(166/25.4,7.5),layout='constrained')
    gs=fig.add_gridspec(4,1,height_ratios=[12,1,1,1])
    ax=fig.add_subplot(gs[0])
    dx=(east[1]-east[0])/2; dy=(north[1]-north[0])/2
    im=ax.imshow(rgb[0],origin='lower',
        extent=[east[0]-dx,east[-1]+dx,north[0]-dy,north[-1]+dy],interpolation='nearest')
    ax.set(xlabel='East of radar (km)',ylabel='North of radar (km)')
    for b in range(5):
        x=hgt/np.tan(el[b])*np.sin(az[b]); y=hgt/np.tan(el[b])*np.cos(az[b])
        ax.plot(x,y,'+',color='white',ms=6)
        ax.text(x+1,y+1,f'TX {b}',color='white',fontsize=11)
    title=ax.set_title('')
    for i,(label,low,high,kind) in enumerate([
        ('Doppler centroid (m/s)',-doppler_max,doppler_max,'h'),
        ('Doppler RMS width (m/s)',width_min,width_max,'s'),
        ('Relative reconstructed power (dB)',db_min,db_max,'v')]):
        bar=fig.add_subplot(gs[i+1]); ramp=np.linspace(0,1,256)
        hc=np.zeros((1,256,3)); hc[:,:,0]=.08; hc[:,:,1:]=1
        if kind=='h': hc[:,:,0]=.75*(1-ramp)
        if kind=='s': hc[:,:,1]=1-(1-saturation_min)*ramp
        if kind=='v': hc[:,:,2]=ramp
        bar.imshow(hsv_to_rgb(hc),aspect='auto',extent=[low,high,0,1])
        bar.set_yticks([]); bar.set_xlabel(label)
    output.parent.mkdir(parents=True,exist_ok=True)
    frames=[]
    first_index=int(np.argmin(abs(times-1738210528737371)))
    for i,stamp in enumerate(times):
        im.set_data(rgb[i])
        utc=datetime.fromtimestamp(int(stamp)/1e6,timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
        title.set_text(f'{utc}\nHeight {hgt:.2f} km above radar | five-beam HSV')
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()))
        if i==first_index:
            fig.savefig(output.with_suffix('.png'),dpi=dpi)
    frames[0].save(output.with_suffix('.gif'),save_all=True,append_images=frames[1:],
                   duration=frame_ms,loop=0)
    plt.close(fig)
    with h5py.File(output.with_suffix('.style.h5'),'w') as h:
        h.attrs.update(settings,frame_ms=frame_ms,dpi=dpi,source=str(source.resolve()),
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            generator='pmse/render_five_beam_hsv.py',png_frame_index=first_index,
            hue_formula='0.75*(1-clip((centroid/doppler_max+1)/2,0,1))')
    print('Rendered',output,'from HDF5 only',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path)
    p.add_argument('--output',type=Path,required=True)
    for name,default in [('doppler-max',9.568),('width-min',0),('width-max',6.379),
                         ('saturation-min',.2),('db-min',-30),('db-max',0)]:
        p.add_argument('--'+name,type=float,default=default)
    p.add_argument('--frame-ms',type=int,default=250)
    p.add_argument('--dpi',type=int,default=300)
    render(**vars(p.parse_args()))


if __name__=='__main__':
    main()
