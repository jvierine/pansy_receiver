import numpy as n
import pansy_detect as pd
import matplotlib.pyplot as plt
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc

def update_tx_pulses():
    """
    Find transmit pulses for the mesosphere mode.
    Start off where the current metadata ends.
    """
    metadata_dir=pc.tx_metadata_dir
    #metadata_dir = "/media/archive/metadata/tx"
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


        start_idx=pd.find_b13_mode_start(d,
                                         i0=idx0,
                                         i1=idx1)
        
        print("%s found %d pulses b13 mode"%(stuffr.unix2datestr((i0+i*dt)/1e6),len(start_idx)))

        if len(start_idx) > 0:
            data_dict={}
            # let's use 3 as id of standard barker13 4 ms ipp-mode
            mode_id=n.repeat(3,len(start_idx))
            gidx=n.where( (start_idx >= idx0) & (start_idx < idx1) )[0]
            mean_ipp=n.mean(n.diff(start_idx[gidx]))
            print("mean_ipp=%d"%(mean_ipp))
            if len(gidx)>0:
                try:
                    data_dict["id"]=mode_id[gidx]
                    dmw.write(start_idx[gidx],data_dict)
                except:
                    import traceback
                    traceback.print_exc()

        # find 7-bit Barker codes associated with ISR mdoe
        start_idx=pd.find_isr_mode_start(d,
                                         i0=idx0,
                                         i1=idx1)
        
        print("%s found %d pulses isr mode"%(stuffr.unix2datestr((i0+i*dt)/1e6),20*len(start_idx)))

        if len(start_idx) > 0:
            data_dict={}
            # let's use 2 as id of standard isr-mode
            mode_id=n.repeat(2,len(start_idx))
            gidx=n.where( (start_idx >= idx0) & (start_idx < idx1) )[0]
            
            if len(gidx)>0:
                try:
                    data_dict["id"]=mode_id[gidx]
                    dmw.write(start_idx[gidx],data_dict)
                except:
                    import traceback
                    traceback.print_exc()


        # search for the start of a continuous 20 IPP sequence
        start_idx=pd.find_m_mode_start(d,
                                       i0=idx0,
                                       i1=idx1,
                                       debug=False)
        
        print("%s found %d pulses mesosphere mode"%(stuffr.unix2datestr((i0+i*dt)/1e6),20*len(start_idx)))

        if len(start_idx)>0:
            data_dict={}
            # let's use 1 as id of standard M-mode
            mode_id=n.array(n.repeat(1,len(start_idx)),dtype=n.uint8)
            gidx=n.where( (start_idx >= idx0) & (start_idx < idx1) )[0]
            if len(gidx)>0:
                print("%d in range"%(20*len(gidx)))
                try:
                    data_dict["id"]=mode_id[gidx]
                    dmw.write(start_idx[gidx],data_dict)
                except:
                    import traceback
                    traceback.print_exc()
    

if __name__ == "__main__":
    while True:
        print("looking for new tx pulses")
        update_tx_pulses()
        time.sleep(10)
