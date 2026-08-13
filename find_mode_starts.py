import numpy as n
import pansy_detect as pd
import matplotlib.pyplot as plt
import digital_rf as drf
import os
import json
from pathlib import Path
import stuffr
import time
import pansy_config as pc
import traceback

ENABLE_ISR_MODE = os.environ.get("PANSY_ENABLE_ISR_MODE", "0") == "1"
PROGRESS_PATH = Path(
    os.environ.get(
        "PANSY_FIND_MODE_PROGRESS_PATH",
        Path.home() / ".local/state/pansy-receiver/find_mode_starts.json",
    )
)
WINDOW_SAMPLES = int(os.environ.get("PANSY_FIND_MODE_WINDOW_SAMPLES", "10000000"))
RAW_LAG_SAMPLES = int(float(os.environ.get("PANSY_FIND_MODE_RAW_LAG_SECONDS", "5")) * 1e6)
MAX_WINDOWS_PER_PASS = int(os.environ.get("PANSY_FIND_MODE_MAX_WINDOWS_PER_PASS", "360"))
IDLE_SLEEP_SECONDS = float(os.environ.get("PANSY_FIND_MODE_IDLE_SLEEP_SECONDS", "10"))


def load_progress(path=PROGRESS_PATH):
    try:
        return int(json.loads(path.read_text())["next_sample"])
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def save_progress(next_sample, path=PROGRESS_PATH):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps({"next_sample": int(next_sample)}, indent=2) + "\n")
    temporary.replace(path)


def write_mode_metadata(writer, start_indices, idx0, idx1, mode_id):
    in_window = n.asarray(start_indices)[
        (n.asarray(start_indices) >= idx0) & (n.asarray(start_indices) < idx1)
    ]
    if in_window.size == 0:
        return 0
    writer.write(
        in_window,
        {"id": n.full(in_window.size, mode_id, dtype=n.uint8)},
    )
    return int(in_window.size)

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
        except Exception:
            traceback.print_exc()
            print("couldn't read metadata; starting new tx metadata")
            db=[-1,-1]
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


    i0 = load_progress()
    if i0 is None or i0 < b[0] or i0 > b[1]:
        i0=b[0]
        if db[1] != -1:
            # Bootstrap old installations from the latest output metadata.
            i0=db[1]+10*1600
    print("starting at %s"%(stuffr.unix2datestr(i0/1e6)))

    dt=WINDOW_SAMPLES
    safe_end = b[1] - RAW_LAG_SAMPLES
    n_windows = max(0, int(n.floor((safe_end-i0)/dt)))
    n_windows = min(n_windows, MAX_WINDOWS_PER_PASS)
    

    for i in range(n_windows):
        idx0=i0+i*dt
        idx1=i0+i*dt+dt


        # search for the start of a continuous 20 IPP sequence
        mesosphere_starts=pd.find_m_mode_start(d,
                                               i0=idx0,
                                               i1=idx1,
                                               debug=False)
        
        print("%s found %d pulses mesosphere mode"%(stuffr.unix2datestr(idx0/1e6),20*len(mesosphere_starts)))
        mesosphere_count = write_mode_metadata(dmw, mesosphere_starts, idx0, idx1, 1)
        if mesosphere_count:
            print("%d in range"%(20*mesosphere_count))
        else:
            # Barker-13 and mesosphere modes are mutually exclusive on this
            # time scale. Avoid its expensive FFT scan when M-mode is present.
            b13_starts=pd.find_b13_mode_start(d, i0=idx0, i1=idx1)
            print("%s found %d pulses b13 mode"%(stuffr.unix2datestr(idx0/1e6),len(b13_starts)))
            write_mode_metadata(dmw, b13_starts, idx0, idx1, 3)

            if ENABLE_ISR_MODE and len(b13_starts) == 0:
                isr_starts=pd.find_isr_mode_start(d, i0=idx0, i1=idx1)
                print("%s found %d pulses isr mode"%(stuffr.unix2datestr(idx0/1e6),20*len(isr_starts)))
                write_mode_metadata(dmw, isr_starts, idx0, idx1, 2)

        # Persist each completed window so a service restart never recreates a
        # long scan backlog.
        save_progress(idx1)
    next_sample = i0 + n_windows * dt
    if n_windows > 0:
        print(
            "mode scan progress %s; lag %.1f s"
            % (stuffr.unix2datestr(next_sample/1e6), max(0, b[1]-next_sample)/1e6)
        )
    return n_windows

if __name__ == "__main__":
    while True:
        print("looking for new tx pulses")
        n_processed = update_tx_pulses()
        if n_processed == 0:
            time.sleep(IDLE_SLEEP_SECONDS)
