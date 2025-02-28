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
    start_idx=dmr.get_bounds()[0]
    # if we already have something. continue from end
    try:
        dmm = drf.DigitalMetadataReader(pc.mesomode_metadata_dir)
        mmb=dmm.get_bounds()
        start_idx=mmb[1]
    except:
        print("no mm metadata")
        exit(0)

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

    db = dmr.get_bounds()
    max_gap = 20*1600+1600
    block=10*60*1000000
    i0=start_idx
#    meso_blocks=[]
    meso_start=-1
    meso_prev=-1
    while i0<(db[1]-2*block):
        i1=i0+block
        data_dict = dmr.read(i0, i1, "id")
        if len(data_dict.keys()) == 0:
            print("no meso-mode")
        else:
            for k in data_dict.keys():
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
                        dmw.write([meso_end],pdata_dict)
                    except:
                        traceback.print_exc()
                    # start new
                    meso_start=k
                    meso_prev=k
        i0+=block


if __name__ == "__main__":
    while True:
        find_blocks()
        time.sleep(3600)

