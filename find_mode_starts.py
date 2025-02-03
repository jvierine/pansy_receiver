import numpy as n
import pansy_detect as pd
import matplotlib.pyplot as plt
import digital_rf as drf
import os
import stuffr
import time


def update_tx_pulses():
    """
    Find transmit pulses for the mesosphere mode.
    Start off where the current metadata ends.
    """
    metadata_dir = "/media/archive/metadata/tx"
    db=[-1,-1]
    if os.path.exists(metadata_dir):
        print("metadata directory exists. searching for last timestamp")
        try:
            dmr = drf.DigitalMetadataReader(metadata_dir)
            db = dmr.get_bounds()
            print(db)
        except:
            print("couldn't read metadata")
    else:
        os.system("mkdir -p %s"%(metadata_dir))

    # setup the directory and file cadence.
    # use 1 MHz, as this is the sample-rate and thus a
    # natural resolution for timing.
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


    i0=b[0]
    if db[1] != -1:
        # start where we left off, instead of the start
        i0=db[1]+10*1600
    print("starting at %s"%(stuffr.unix2datestr(i0/1e6)))

    dt=10000000
    n_windows = int(n.floor(((b[1]-1000000)-i0)/dt))
    

    for i in range(n_windows):
        idx0=i0+i*dt
        idx1=i0+i*dt+dt
        # search for the start of a continuous 20 IPP sequence
        start_idx=pd.find_m_mode_start(d,
                                       i0=idx0,
                                       i1=idx1,
                                       ch="ch007", # channel 007 is the transmit sample
                                       debug=False)
        print("%s found %d pulses"%(stuffr.unix2datestr((i0+i*dt)/1e6),20*len(start_idx)))

        if len(start_idx)>0:
            print("writing metadata")
            data_dict={}
            # let's use 1 as id of standard M-mode
            mode_id=n.array(n.repeat(1,len(start_idx)),dtype=n.uint8)
            gidx=n.where( (start_idx >= idx0) & (start_idx < idx1) )[0]
            if len(gidx)>0:
                print("%d in range"%(20*len(gidx)))
                data_dict["id"]=mode_id[gidx]
                dmw.write(start_idx[gidx],data_dict)
    

if __name__ == "__main__":
    while True:
        print("looking for new tx pulses")
        update_tx_pulses()
        time.sleep(10)
