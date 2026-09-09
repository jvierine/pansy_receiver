"""Compare all five PMSE angular reconstructions in one new figure.

Uses native XC2, per-beam phasecal, and the positive meteor-calibration sign.
Outputs the reconstructed numerical products and provenance to HDF5.
"""
import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C = 299792458.0
F_RADAR = 47e6


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--file', type=Path, required=True)
    p.add_argument('--group', default='1738210528737371')
    p.add_argument('--reference', type=Path, default=Path(__file__).parent/'reference')
    p.add_argument('--range-km', type=float, default=85.5)
    p.add_argument('--velocity-mps', type=float, default=0)
    p.add_argument('--radius-deg', type=float, default=20)
    p.add_argument('--grid-size', type=int, default=151)
    p.add_argument('--modes', type=int, default=18)
    p.add_argument('--output', type=Path, default=Path('figures/five_beam_image.png'))
    a = p.parse_args()
    if not 0 < a.radius_deg < 60 or a.grid_size < 5 or a.modes < 1:
        p.error('Require 0 < radius < 60 degrees, grid-size >= 5, modes >= 1')
    with h5py.File(a.file) as h:
        g = h[a.group]
        xc = g['xc_arr'][()].astype(np.complex128)
        pairs = g['ch_pairs'][()]
        ranges = np.arange(g['r0'][()],g['r1'][()],g['rdec'][()])*C/2e9
        f = np.fft.fftshift(np.fft.fftfreq(int(g['n_fft'][()]),
            d=5*int(g['ipp'][()])/1e6))[int(g['f0'][()]):int(g['f1'][()])]
    assert xc.shape[1] == 5
    assert np.array_equal(pairs[:7], np.column_stack([np.arange(7)]*2))
    velocity = f*C/(2*F_RADAR)
    if not ranges.min() <= a.range_km <= ranges.max():
        p.error('Requested range is outside the archived data')
    if not velocity.min() <= a.velocity_mps <= velocity.max():
        p.error('Requested velocity is outside the archived data')
    ri = int(np.argmin(abs(ranges-a.range_km)))
    fi = int(np.argmin(abs(velocity-a.velocity_mps)))
    groups = {}
    for line in (a.reference/'cfg/antpos.csv').read_text().splitlines()[1:]:
        fields = line.split(',')
        groups.setdefault(fields[1], []).append([float(x) for x in fields[4:7]])
    connections = [line.split(',')[0] for line in
                   (a.reference/'cfg/connections.txt').read_text().splitlines()]
    positions = np.array([np.mean(groups[x],axis=0) for x in connections if x!='RFTX'])
    phases = []
    calibration_hashes = []
    for beam in range(5):
        path = a.reference/f'data/phases{beam}.h5'
        with h5py.File(path) as h:
            phases.append(h['phasecal'][()])
        calibration_hashes.append(hashlib.sha256(path.read_bytes()).hexdigest())
    phases = np.asarray(phases)
    assert phases.shape == (5,7) and np.isfinite(phases).all()
    baseline = positions[pairs[:,1]]-positions[pairs[:,0]]
    axis = np.linspace(-np.sin(np.deg2rad(a.radius_deg)),
                       np.sin(np.deg2rad(a.radius_deg)),a.grid_size)
    uu,vv = np.meshgrid(axis,axis)
    mask = uu**2+vv**2 <= np.sin(np.deg2rad(a.radius_deg))**2
    directions = np.column_stack([uu[mask],vv[mask],-np.sqrt(1-uu[mask]**2-vv[mask]**2)])
    A = np.exp(1j*2*np.pi*F_RADAR/C*(baseline@directions.T))
    U,s,Vh = np.linalg.svd(A,full_matrices=False)
    keep = min(a.modes,len(s))
    inverse = (Vh[:keep].conj().T/(s[:keep]+1))@U[:,:keep].conj().T
    # Check phase correction and steering against a known source, separately
    # for each beam, rather than treating agreement between beams as validation.
    test_index = np.argmin(np.sum((directions[:,:2]-[0.08,-0.06])**2,axis=1))
    truth = A[:,test_index]
    power = np.abs(xc[:7])
    outer = abs(f) >= .75*np.max(abs(f))
    noise = np.median(power[:,:,outer,:],axis=(2,3))
    assert np.all(noise > 0)
    images, coherences, snr_db = [], [], []
    for beam in range(5):
        factor = np.exp(1j*(phases[beam,pairs[:,0]]-phases[beam,pairs[:,1]]))
        corrected_test = (truth/factor)*factor
        assert np.allclose(corrected_test,truth)
        assert np.argmax(abs(A.conj().T@corrected_test)) == test_index
        signal = np.maximum(power[:,beam,fi,ri]-noise[:,beam],noise[:,beam])
        coherence = xc[:,beam,fi,ri]*factor/np.sqrt(signal[pairs[:,0]]*signal[pairs[:,1]])
        images.append(abs(inverse@coherence))
        coherences.append(coherence)
        snr_db.append(10*np.log10(power[:,beam,fi,ri].sum()/noise[:,beam].sum()))
    images = np.asarray(images)
    assert np.isfinite(images).all() and np.all(images.max(axis=1)>0)
    plt.rcParams.update({'font.size':11,'axes.titlesize':11,'axes.labelsize':11,
                         'xtick.labelsize':11,'ytick.labelsize':11})
    fig,axes = plt.subplots(3,2,figsize=(166/25.4,8.2),layout='constrained')
    cmap = plt.get_cmap('inferno').copy()
    cmap.set_bad(cmap(0))
    vmax = images.max()
    for beam,ax in enumerate(axes.flat):
        if beam == 5:
            ax.axis('off')
            continue
        raster = np.full(uu.shape,np.nan)
        raster[mask] = images[beam]
        m = ax.pcolormesh(uu,vv,np.ma.masked_invalid(raster),shading='nearest',
                         cmap=cmap,vmin=0,vmax=vmax,edgecolors='none',antialiased=False)
        name = 'Zenith / TX 0' if beam==0 else f'TX {beam}'
        ax.set(xlabel='u',ylabel='v',title=f'{name} | {snr_db[beam]:.1f} dB',aspect='equal',
               xlim=(axis[0],axis[-1]),ylim=(axis[0],axis[-1]))
        ax.set_xticks([-0.3,0,0.3] if a.radius_deg >= 18 else [-0.1,0,0.1])
        ax.set_yticks(ax.get_xticks())
    fig.colorbar(m,ax=axes.ravel().tolist(),location='bottom',fraction=.045,
                 pad=.02,label='Reconstruction magnitude (a.u.; common scale)')
    stamp = datetime.fromtimestamp(int(a.group)/1e6,timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
    fig.suptitle(f'{stamp}\n{ranges[ri]:.2f} km | {velocity[fi]:.2f} m/s | {keep} SVD modes',fontsize=11)
    axes[2,1].text(0,.95,'Five TX beams\nSame range and Doppler bin\n\nBeam-specific calibration\nPositive meteor convention\n\nTitles: power / noise (dB)\nNot absolute echo intensity\n\nPreliminary calibration',
                   transform=axes[2,1].transAxes,va='top',fontsize=11)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(a.output,dpi=300)
    with h5py.File(a.output.with_suffix('.h5'),'w') as h:
        for key,value in dict(magnitude=images,coherence=coherences,beam=np.arange(5),
             directions=directions,u_grid=uu,v_grid=vv,imaged_mask=mask,
             velocity_mps=velocity[fi],range_km=ranges[ri],singular_values=s,
             phasecal_rad=phases,antenna_positions_m=positions,power_noise_db=snr_db,
             noise_power=noise,ch_pairs=pairs).items():
            h[key]=value
        h.attrs.update(source=str(a.file),group=a.group,
            source_sha256=hashlib.sha256(a.file.read_bytes()).hexdigest(),
            calibration_sha256=calibration_hashes,retained_modes=keep,
            radius_deg=a.radius_deg,convention='positive calibration and forward sign',
            generator='pmse/image_five_beams.py',shared_color_max=vmax)
    print(a.output, 'range', ranges[ri], 'velocity', velocity[fi], 'power/noise dB',snr_db,
          'image shape',images.shape,flush=True)


if __name__ == '__main__':
    main()
