"""Inverse transform retained, windowed decoded complementary echoes."""
import argparse
from pathlib import Path
import h5py
import numpy as np
from .sparse_record import read_complete


def reconstruct(record,beam,restore_pair_scale=False):
    index=np.asarray(record['sparse_index']);z=np.asarray(record['spectra'],dtype=np.complex64)
    nfft=len(record['doppler_frequency_hz']);nr=len(record['slant_range_km'])
    full=np.zeros((len(z),7,nfft,nr),dtype=np.complex64)
    select=index[:,0]==beam
    full[:,:,index[select,1],index[select,2]]=z[:,:,select]
    echo=np.fft.ifft(np.fft.ifftshift(full,axes=2),axis=2).astype(np.complex64)
    if restore_pair_scale:
        weights=np.asarray(record['look_window'],dtype=np.float64)
        if weights.shape!=(len(z),nfft) or not np.isfinite(weights).all() or np.any(weights<0):
            raise ValueError('Invalid saved pair weights')
        scale=weights.max(axis=1)
        if np.any(scale<=0) or not np.all((weights==0)|(weights==scale[:,None])):
            raise ValueError('Pair-voltage scaling requires rectangular weights; cannot undo historical Hann')
        # A constant gain commutes with spectral projection. Do not divide by
        # zero at rejected/padded times or pretend their projected values are
        # recovered measurements.
        echo=(echo/scale[:,None,None,None]).astype(np.complex64)
    return echo


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    p.add_argument('--key',type=int,required=True);p.add_argument('--beam',type=int,choices=range(5),required=True)
    p.add_argument('--pair-voltage',action='store_true',help='Undo constant rectangular noise normalization; refuses historical Hann records')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();record=read_complete(a.root,a.key)
    echo=reconstruct(record,a.beam,restore_pair_scale=a.pair_voltage)
    with h5py.File(a.output,'w') as h:
        h.create_dataset('decoded_echo',data=echo,compression='gzip',shuffle=True)
        for key in ['look_window','look_pair_time0_sample_by_beam','slant_range_km']:
            if key in record:h[key]=record[key]
        h.attrs.update(generator='pmse/reconstruct_echoes.py',beam=a.beam,key=a.key,
            axis_order='coherent_look,rx_channel,pair_time,range',
            amplitude_domain='decoded complementary pair voltage' if a.pair_voltage else 'stored FFT normalization',
            description='Retained signal only; complementary echo series; omitted bins zero-filled; invalid/padded times are not measurements; no original ADC recovery')
