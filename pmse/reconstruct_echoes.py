"""Inverse transform retained, windowed decoded complementary echoes."""
import argparse
from pathlib import Path
import h5py
import numpy as np
from .sparse_record import read_complete


def reconstruct(record,beam):
    index=np.asarray(record['sparse_index']);z=np.asarray(record['spectra'],dtype=np.complex64)
    nfft=len(record['doppler_frequency_hz']);nr=len(record['slant_range_km'])
    full=np.zeros((len(z),7,nfft,nr),dtype=np.complex64)
    select=index[:,0]==beam
    full[:,:,index[select,1],index[select,2]]=z[:,:,select]
    echo=np.fft.ifft(np.fft.ifftshift(full,axes=2),axis=2).astype(np.complex64)
    return echo


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    p.add_argument('--key',type=int,required=True);p.add_argument('--beam',type=int,choices=range(5),required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();record=read_complete(a.root,a.key)
    echo=reconstruct(record,a.beam)
    with h5py.File(a.output,'w') as h:
        h.create_dataset('decoded_echo',data=echo,compression='gzip',shuffle=True)
        for key in ['look_window','look_pair_time0_sample_by_beam','slant_range_km']:
            if key in record:h[key]=record[key]
        h.attrs.update(generator='pmse/reconstruct_echoes.py',beam=a.beam,key=a.key,
            axis_order='coherent_look,rx_channel,pair_time,range',
            description='Retained signal only; windowed complementary echo series; omitted bins zero-filled; no original ADC recovery')
