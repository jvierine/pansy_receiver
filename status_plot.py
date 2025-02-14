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
d1=drf.DigitalMetadataReader(tx_metadata_dir)
d2=drf.DigitalMetadataReader(detections_metadata_dir)
d3=drf.DigitalMetadataReader(cut_metadata_dir)
d4=drf.DigitalMetadataReader(mesomode_metadata_dir)
d5=drf.DigitalMetadataReader(xc_metadata_dir)#meteor_cal_metadata_dir="/media/analysis/metadata/meteor_cal"
d6=drf.DigitalMetadataReader(simple_fit_metadata_dir)
dr=drf.DigitalRFReader(raw_voltage_dir)


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

print("latest tx %s"%(latest_tx))
print("latest mf %s"%(latest_mf))
print("latest det %s"%(latest_det))
print("latest cut %s"%(latest_cut))
print("latest mode %s"%(latest_mode))
print("latest xc %s"%(latest_xc))
print("latest fit %s"%(latest_fit))