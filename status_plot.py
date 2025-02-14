import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc
import traceback

# how many receiver restarts in the log 
# sandra-rx events

# plot latest 24 hours of pmse

# plot latest 24 hours of meteors

d0=drf.DigitalMetadataReader(pc.mf_metadata_dir)
d1=drf.DigitalMetadataReader(pc.tx_metadata_dir)
d2=drf.DigitalMetadataReader(pc.detections_metadata_dir)
d3=drf.DigitalMetadataReader(pc.cut_metadata_dir)
d4=drf.DigitalMetadataReader(pc.mesomode_metadata_dir)
d5=drf.DigitalMetadataReader(pc.xc_metadata_dir)#meteor_cal_metadata_dir="/media/analysis/metadata/meteor_cal"
d6=drf.DigitalMetadataReader(pc.simple_fit_metadata_dir)
dr=drf.DigitalRFReader(pc.raw_voltage_dir)

tnow=time.time()*1e6
mfb=d0.get_bounds()
latest_mf=stuffr.unix2datestr(mfb[1]/1e6)

txb=d1.get_bounds()
latest_tx=stuffr.unix2datestr(txb[1]/1e6)

detb=d2.get_bounds()
latest_det=stuffr.unix2datestr(detb[1]/1e6)

cutb=d3.get_bounds()
latest_cut=stuffr.unix2datestr(cutb[1]/1e6)

modeb=d4.get_bounds()
latest_mode=stuffr.unix2datestr(modeb[1]/1e6)

xcb=d5.get_bounds()
latest_xc=stuffr.unix2datestr(xcb[1]/1e6)

fitb=d6.get_bounds()
latest_fit=stuffr.unix2datestr(fitb[1]/1e6)

b=dr.get_bounds("ch000")
latest_raw=stuffr.unix2datestr(b[1]/1e6)

print("latest raw voltage %s (%1.0f s behind)"%(latest_raw,(tnow-b[1])/1e6))
print("latest tx %s (%1.0f s behind)"%(latest_tx,(tnow-txb[1])/1e6))
print("latest mf %s (%1.0f s behind)"%(latest_mf,(tnow-mfb[1])/1e6))
print("latest det %s (%1.0f s behind)"%(latest_det,(tnow-detb[1])/1e6))
print("latest cut %s (%1.0f s behind)"%(latest_cut,(tnow-cutb[1]/1e6)))
print("latest mode %s (%1.0f s behind)"%(latest_mode,(tnow-modeb[1])/1e6))
print("latest xc %s (%1.0f s behind)"%(latest_xc,(tnow-xcb[1])/1e6))
print("latest fit %s (%1.0f s behind)"%(latest_fit,(tnow-fitb[1])/1e6))