import numpy as n
import pansy_detect as pd
import matplotlib.pyplot as plt
import digital_rf as drf
import os
import stuffr

metadata_dir = "/media/archive/metadata/tx"
os.system("mkdir -p %s"%(metadata_dir))

subdirectory_cadence_seconds = 3600
file_cadence_seconds = 60
samples_per_second_numerator = 1000000
samples_per_second_denominator = 1
file_name = "m_tx"

dmw = drf.DigitalMetadataWriter(
    metadata_dir,
    subdirectory_cadence_seconds,
    file_cadence_seconds,
    samples_per_second_numerator,
    samples_per_second_denominator,
    file_name,
)

# this is where the data is
d=drf.DigitalRFReader("/media/archive/")
# tx channel bounds
b=d.get_bounds("ch007")

dt=60000000
n_windows = int(n.floor((b[1]-b[0])/dt))
i0=b[0]
for i in range(n_windows):
    start_idx=pd.find_m_mode_start(d,
                                   i0=i0+i*dt,
                                   i1=i0+i*dt+dt,
                                   ch="ch007", # channel 007 is the transmit sample
                                   debug=False)
    print("%s found %d pulses"%(stuffr.unix2datestr((i0+i*dt)/1e6),len(start_idx)))

    if len(start_idx)>0:
        print("writing metadata")
        data_dict={}
        # let's use 1 as id of standard M-mode
        mode_id=n.array(n.repeat(1,len(start_idx)),dtype=n.int32)
        data_dict["id"]=mode_id
        dmw.write(start_idx,data_dict)
    
