"""Reconstruct angular images from the archived first PANSY PMSE block.

Default convention matches image_points1/u_phases. --legacy-signs reproduces
the negative calibration and forward-matrix signs of image_fivebeam.py.
The original uncropped convenience file is unavailable, so this is a
reconstruction from the surviving cropped spectra, not pixel-identical output.
"""
import argparse
from pathlib import Path
import hashlib
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C = 299792458.0
F_RADAR = 47e6


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('xc2'))
    p.add_argument('--file', type=Path)
    p.add_argument('--group', default='1738210528737371')
    p.add_argument('--reference', type=Path, default=Path(__file__).parent/'reference')
    p.add_argument('--range-km', type=float, default=85.5)
    p.add_argument('--velocities', type=float, nargs='+', default=[-6, 0, 6])
    p.add_argument('--modes', type=int, default=18)
    p.add_argument('--legacy-signs', action='store_true')
    p.add_argument('--output', type=Path, default=Path('figures/first_image.png'))
    a = p.parse_args()
    source = a.file or a.root/'2025-01-30T04-00-00'/'xc@1738210500.h5'
    with h5py.File(source) as h:
        g = h[a.group]
        xc = g['xc_arr'][:, 0].astype(np.complex128)
        pairs = g['ch_pairs'][()]
        r = np.arange(g['r0'][()],g['r1'][()],g['rdec'][()])*C/2e9
        f = np.fft.fftshift(np.fft.fftfreq(int(g['n_fft'][()]),
                 d=5*int(g['ipp'][()])/1e6))[int(g['f0'][()]):int(g['f1'][()])]
    velocity = f*C/(2*F_RADAR)
    # The source CSV is an existing instrument geometry input, not an output.
    groups = {}
    for line in (a.reference/'cfg/antpos.csv').read_text().splitlines()[1:]:
        fields = line.split(',')
        groups.setdefault(fields[1], []).append([float(x) for x in fields[4:7]])
    connections = [line.split(',')[0] for line in
                   (a.reference/'cfg/connections.txt').read_text().splitlines()]
    positions = np.array([np.mean(groups[x],axis=0) for x in connections if x!='RFTX'])
    with h5py.File(a.reference/'data/phases0.h5') as h:
        phase = h['phasecal'][()]
    baseline = positions[pairs[:,1]]-positions[pairs[:,0]]
    # Match the original 100x100 candidate grid and nine-degree zenith mask.
    axis = np.linspace(-np.arcsin(np.pi*30/180),np.arcsin(np.pi*30/180),100)
    uu,vv = np.meshgrid(axis,axis)
    mask = uu**2+vv**2 < np.sin(np.deg2rad(9))**2
    directions = np.column_stack([uu[mask],vv[mask],-np.sqrt(1-uu[mask]**2-vv[mask]**2)])
    sign = -1 if a.legacy_signs else 1
    A = np.exp(sign*1j*2*np.pi*F_RADAR/C*(baseline@directions.T))
    U,s,Vh = np.linalg.svd(A,full_matrices=False)
    keep = min(a.modes,len(s))
    inverse = (Vh[:keep].conj().T/(s[:keep]+1))@U[:,:keep].conj().T
    # Synthetic point-source check of baseline/steering direction and sign.
    synthetic = directions[np.argmin(np.sum(directions[:,:2]**2,axis=1))]
    steering = np.exp(sign*1j*2*np.pi*F_RADAR/C*(baseline@synthetic))
    recovered = np.argmax(np.abs(A.conj().T@steering))
    assert np.allclose(directions[recovered], synthetic), 'Synthetic direction check failed'
    power = np.abs(xc[:7])
    outer = np.abs(f) >= .75*np.max(np.abs(f))
    noise = np.median(power[:,outer,:],axis=(1,2))
    ri = int(np.argmin(np.abs(r-a.range_km)))
    fi = [int(np.argmin(np.abs(velocity-v))) for v in a.velocities]
    phase_factor = np.exp(sign*1j*(phase[pairs[:,0]]-phase[pairs[:,1]]))
    images=[]
    for index in fi:
        # Use a noise-relative floor: an absolute tiny floor makes sub-noise
        # channels dominate the inversion with arbitrarily large coherence.
        signal = np.maximum(power[:,index,ri]-noise,noise)
        coherence = xc[:,index,ri]*phase_factor/np.sqrt(signal[pairs[:,0]]*signal[pairs[:,1]])
        result = np.abs(inverse@coherence)
        assert np.all(np.isfinite(result)) and np.max(result)>0
        images.append(result)
    plt.rcParams.update({'font.size':11,'axes.titlesize':11,'axes.labelsize':11,
                         'xtick.labelsize':11,'ytick.labelsize':11})
    nrows = int(np.ceil((len(fi)+1)/2))
    fig,axes=plt.subplots(nrows,2,figsize=(166/25.4,3.3*nrows),layout='constrained')
    axes=axes.flatten()
    for unused in axes[len(fi)+1:]:
        unused.set_visible(False)
    summed=power.sum(axis=0)
    noise_sum=np.median(summed[outer])
    m=axes[0].pcolormesh(velocity,r,10*np.log10(summed/noise_sum).T,
                         shading='auto',vmin=0,vmax=35,cmap='plasma')
    axes[0].axhline(r[ri],color='cyan',ls='--')
    axes[0].set(xlabel='Radial velocity (m/s)',ylabel='Range (km)',title='Zenith spectrum')
    fig.colorbar(m,ax=axes[0],label='Power / noise (dB)')
    for ax,index,result in zip(axes[1:],fi,images):
        # Render contiguous grid cells, avoiding white spaces between markers.
        angular_image = np.full(uu.shape, np.nan)
        angular_image[mask] = result
        cmap = plt.get_cmap('inferno').copy()
        cmap.set_bad(cmap(0.0))
        m=ax.pcolormesh(uu,vv,np.ma.masked_invalid(angular_image),
                        shading='nearest',cmap=cmap,vmin=0,
                        edgecolors='none',antialiased=False,rasterized=True)
        half_cell = (axis[1]-axis[0])/2
        extent = np.max(np.abs(directions[:,:2]))+half_cell
        ax.set_xlim(-extent,extent)
        ax.set_ylim(-extent,extent)
        ax.set(xlabel='u',ylabel='v',title=f'{velocity[index]:.2f} m/s',aspect='equal')
        fig.colorbar(m,ax=ax,label='Magnitude (a.u.)')
    convention='legacy signs' if a.legacy_signs else 'meteor calibration convention'
    fig.suptitle(f'2025-01-30 04:15:28 UTC | {r[ri]:.2f} km\n{keep} SVD modes | {convention}',fontsize=11)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(a.output,dpi=300)
    with h5py.File(a.output.with_suffix('.h5'),'w') as h:
        h['directions']=directions
        h['magnitude']=images
        h['u_grid']=uu
        h['v_grid']=vv
        h['imaged_mask']=mask
        h['velocity_mps']=velocity[fi]
        h['range_km']=r[ri]
        h['singular_values']=s
        h['antenna_positions_m']=positions
        h['phasecal_rad']=phase
        h.attrs['source']=str(source)
        h.attrs['group']=a.group
        h.attrs['source_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
        h.attrs['convention']=convention
        h.attrs['generator']='pmse/image_first_block.py'
    print(a.output, 'range',r[ri], 'velocities',velocity[fi], 'modes',keep)


if __name__=='__main__':
    main()
