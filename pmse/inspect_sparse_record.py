"""Read and validate sparse Digital Metadata; optionally render its detected cells."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
from pmse.sparse_record import read_complete,channel_covariance


def inspect(root,key,output=None):
    record=read_complete(root,key)
    indices=np.asarray(record['sparse_index'])
    z=np.asarray(record['spectra'])
    assert z.ndim==3 and z.shape[1]==7 and z.shape[2]==len(indices)
    assert np.isfinite(z).all()
    expected=(abs(z.astype('c16'))**2).mean(axis=0)
    np.testing.assert_allclose(record['mean_power'],expected,rtol=2e-6)
    covariance=channel_covariance(record)
    np.testing.assert_allclose(covariance,covariance.conj().transpose(0,2,1),atol=1e-8)
    look_count=int(record['n_coherent_looks'])
    report=dict(mode_start_utc=datetime.fromtimestamp(key/1e6,timezone.utc).isoformat(),
        mode_seconds=float(record.get('mode_duration_seconds',0)),
        covered_seconds=float(record.get('coverage_seconds',0)),
        looks=look_count,selected_cells=len(indices),tested_cells=int(record['n_tested_cells']),
        high_velocity_cells=int(np.count_nonzero(record['high_velocity_candidate'])),
        sparse_voltage_bytes=int(z.nbytes),
        null_relative_std=float(record['null_power_relative_standard_deviation']),
        processing_seconds=float(record.get('processing_wall_seconds',0)))
    print(json.dumps(report,indent=2),flush=True)
    if output is not None:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size':11,'axes.titlesize':11,'axes.labelsize':11,
                             'xtick.labelsize':11,'ytick.labelsize':11})
        fig,axes=plt.subplots(3,2,figsize=(166/25.4,8.7),layout='constrained')
        velocity=np.asarray(record['radial_velocity_mps'])
        height=np.asarray(record['nominal_altitude_km'])
        cmap=plt.get_cmap('plasma').copy(); cmap.set_bad(cmap(0))
        for beam,ax in enumerate(axes.flat):
            if beam==5:
                ax.axis('off'); continue
            power=np.full((len(velocity),height.shape[1]),np.nan)
            chosen=indices[:,0]==beam
            power[indices[chosen,1],indices[chosen,2]]=record['detector_snr_db'][chosen]
            mesh=ax.pcolormesh(velocity,height[beam],power.T,shading='nearest',cmap=cmap,
                              vmin=3,vmax=30,edgecolors='none',antialiased=False)
            ax.set(xlabel='Radial velocity (m/s)',ylabel='Nominal altitude (km)',
                   title='Zenith / TX 0' if beam==0 else f'TX {beam}',xlim=(-92,92),ylim=(50,120))
            ax.set_xticks([-90,0,90])
        fig.colorbar(mesh,ax=axes.ravel().tolist(),location='bottom',fraction=.03,
                     label='Detected ensemble signal / noise (dB)')
        stamp=datetime.fromtimestamp(key/1e6,timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
        fig.suptitle(f'{stamp}\n{report["mode_seconds"]:.2f} s | {look_count} coherent looks | sparse PMSE/PMWE',fontsize=11)
        output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(output,dpi=300); plt.close(fig)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path); p.add_argument('--key',type=int)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.key is None:
        import digital_rf as drf
        a.key=int(drf.DigitalMetadataReader(str(a.root)).get_bounds()[1])
    inspect(a.root,a.key,a.output)


if __name__=='__main__':
    main()
