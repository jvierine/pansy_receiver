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


MESOMODE_MAX_GAP_SAMPLES = int(
    float(os.environ.get("PANSY_MESOMODE_MAX_GAP_SECONDS", "0.5")) * 1e6
)


def closed_mode_blocks(start_indices, process_end, max_gap=MESOMODE_MAX_GAP_SAMPLES):
    """Group detected M-mode cycle starts into blocks known to be closed."""
    starts = n.asarray(start_indices, dtype=n.int64)
    if starts.size == 0:
        return []

    starts = n.unique(starts)
    blocks = []
    block_start = int(starts[0])
    previous = block_start
    for current in starts[1:]:
        current = int(current)
        if current - previous <= max_gap:
            previous = current
            continue
        blocks.append((block_start, previous))
        block_start = current
        previous = current

    # process_end is already behind live TX metadata. If no cycle start has
    # appeared within max_gap, the trailing block is known to have ended.
    if int(process_end) - previous > max_gap:
        blocks.append((block_start, previous))
    return blocks


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
        start_idx=mmb[1]+1
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

    data_dict = dmr.read(start_idx, process_end, "id")
    meso_starts = [
        k for k in sorted(data_dict) if int(n.asarray(data_dict[k]).reshape(-1)[0]) == 1
    ]
    blocks = closed_mode_blocks(meso_starts, process_end)
    if not blocks:
        print("no closed meso-mode blocks")
        return

    for meso_start, meso_end in blocks:
        print("%s found meso mode %1.2f (s)"%(
            stuffr.unix2datestr(meso_start/1e6),
            (meso_end-meso_start)/1e6))
        odata_dict={}
        odata_dict["start"]=[meso_start]
        odata_dict["end"]=[meso_end]
        try:
            dmw.write([meso_end],odata_dict)
        except Exception:
            traceback.print_exc()


if __name__ == "__main__":
    while True:
        find_blocks()
        time.sleep(300)
