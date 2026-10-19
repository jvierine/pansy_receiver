"""Typed NumPy boundary to Rust DSP; no alternate normalization or I/O schema."""
import ctypes
import os
from pathlib import Path
import numpy as np


class RustDSP:
    def __init__(self, library=None):
        import sys
        suffix='dylib' if sys.platform=='darwin' else 'so'
        path=Path(library or os.environ.get('PANSY_PMSE_DSP_LIBRARY',
            Path(__file__).parent/f'processor/target/release/libpansy_pmse_dsp.{suffix}'))
        self.path=path.resolve()
        self.lib=ctypes.CDLL(str(self.path))
        ptr=ctypes.c_void_p; size=ctypes.c_size_t
        self.lib.pansy_decode_selected.argtypes=[ptr,ptr,size,size,size,ptr,size,ptr,ptr]
        self.lib.pansy_decode_selected.restype=ctypes.c_int
        self.lib.pansy_coherent.argtypes=[ptr,ptr,size,size,size,size,ptr]
        self.lib.pansy_coherent.restype=ctypes.c_int

    @staticmethod
    def checked(code):
        if code:
            raise ValueError(f'Rust DSP rejected input or failed: status {code}')

    def decode_selected(self,rx,tx,cfg,ranges):
        rx=np.ascontiguousarray(rx,dtype=np.complex64)
        tx=np.ascontiguousarray(tx,dtype=np.complex64)
        ranges=np.ascontiguousarray(ranges,dtype=np.uint32)
        if rx.ndim!=2 or rx.shape!=tx.shape or rx.shape[1]!=cfg.ipp_samples or not len(rx):
            raise ValueError('Expected matching nonempty [pulse,IPP] RX and measured TX')
        if ranges.ndim!=1 or not len(ranges) or np.any(ranges.astype(np.uint64)+cfg.tx_samples>cfg.ipp_samples):
            raise ValueError('Range gate extends past linear decoded pulse')
        out=np.empty((len(rx),len(ranges)),dtype=np.complex64)
        power=np.empty(len(rx),dtype=np.float64)
        self.checked(self.lib.pansy_decode_selected(rx.ctypes.data,tx.ctypes.data,len(rx),cfg.ipp_samples,
            cfg.tx_samples,ranges.ctypes.data,len(ranges),out.ctypes.data,power.ctypes.data))
        return out,power

    def coherent_spectrum(self,decoded,cfg,pulse_valid=None):
        from .sparse_record import coherent_window
        x=np.ascontiguousarray(decoded,dtype=np.complex64)
        if x.ndim!=4 or x.shape[:2]!=(7,5) or not 4<=x.shape[2]<=cfg.nfft:
            raise ValueError('Expected [7,5,pulse,range]')
        count=x.shape[2]
        window=coherent_window(count,pulse_valid)
        if window.shape!=(count,):
            raise ValueError('Invalid pulse mask length')
        out=np.empty((7,5,cfg.nfft,x.shape[3]),dtype=np.complex64)
        self.checked(self.lib.pansy_coherent(x.ctypes.data,window.ctypes.data,35,count,
            x.shape[3],cfg.nfft,out.ctypes.data))
        return out
