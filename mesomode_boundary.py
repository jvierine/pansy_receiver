import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt
import pansy_modes as pm
import scipy.signal.windows as sw
import scipy.constants as c
import stuffr
import time
import pansy_config as pc
import os
import traceback

def find_blocks():
    """
    find contiguous blocks of mesosphere mode
    """

    dmr = drf.DigitalMetadataReader(pc.tx_metadata_dir)
    db = dmr.get_bounds()
    start_idx=db[0]
    # if we already have something. continue from end
    try:
        dmm = drf.DigitalMetadataReader(pc.mesomode_metadata_dir)
        mmb=dmm.get_bounds()
        start_idx=mmb[1]
    except Exception:
        print("no mm metadata; starting new mesomode metadata")

    subdirectory_cadence_seconds = 3600
    file_cadence_seconds = 600
    samples_per_second_numerator = 1000000
    samples_per_second_denominator = 1
    file_name = "mesomode"
    os.system("mkdir -p %s"%(pc.mesomode_metadata_dir))

    dmw = drf.DigitalMetadataWriter(
        pc.mesomode_metadata_dir,
        subdirectory_cadence_seconds,
        file_cadence_seconds,
        samples_per_second_numerator,
        samples_per_second_denominator,
        file_name,
    )

    max_gap = 20*1600+1600
    block=10*60*1000000
    processing_lag=2*block
    process_end=db[1]-processing_lag
    print("tx bounds %s - %s; processing through %s"%(
        stuffr.unix2datestr(db[0]/1e6),
        stuffr.unix2datestr(db[1]/1e6),
        stuffr.unix2datestr(process_end/1e6)))
    if start_idx >= process_end:
        print("not enough tx metadata yet for mesomode boundary processing")
        return

    i0=start_idx
#    meso_blocks=[]
    meso_start=-1
    meso_prev=-1
    last_written_end=-1
    while i0<process_end:
        i1=min(i0+block, process_end)
        data_dict = dmr.read(i0, i1, "id")
        if len(data_dict.keys()) == 0:
            print("no meso-mode")
        else:
            for k in sorted(data_dict.keys()):
                if data_dict[k] != 1:
                    continue
                if meso_start == -1:
                    meso_start = k
                    meso_prev=k
                if ((k-meso_prev) > 0) and ((k-meso_prev) < max_gap):
                    meso_prev=k
                if ((k-meso_prev) > 0) and ((k-meso_prev) > max_gap):
                    meso_end=meso_prev
                    #meso_blocks.append({"start":meso_start,"end":meso_end})
                    print("%s found meso mode %1.2f (s)"%(stuffr.unix2datestr(meso_start/1e6), (meso_end-meso_start)/1e6))
                    odata_dict={}
                    odata_dict["start"]=[meso_start]
                    odata_dict["end"]=[meso_end]
                    try:
                        dmw.write([meso_end],odata_dict)
                        last_written_end=meso_end
                    except Exception:
                        traceback.print_exc()
                    # start new
                    meso_start=k
                    meso_prev=k
        i0+=block
    if meso_start != -1 and meso_prev > last_written_end:
        print("%s writing open meso mode through %s (%1.2f s)"%(
            stuffr.unix2datestr(meso_start/1e6),
            stuffr.unix2datestr(meso_prev/1e6),
            (meso_prev-meso_start)/1e6))
        odata_dict={}
        odata_dict["start"]=[meso_start]
        odata_dict["end"]=[meso_prev]
        try:
            dmw.write([meso_prev],odata_dict)
        except Exception:
            traceback.print_exc()


if __name__ == "__main__":
    while True:
        find_blocks()
        time.sleep(300)
