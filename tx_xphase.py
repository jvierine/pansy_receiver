import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc
import cluster_mf as cmf
import traceback
import scipy.fftpack as fp
import itertools
import json
from pathlib import Path
import pansy_modes as pmm


PROGRESS_PATH = Path(
    os.environ.get(
        "PANSY_TX_XPHASE_PROGRESS_PATH",
        Path.home() / ".local/state/pansy-receiver/tx_xphase.json",
    )
)
TX_METADATA_LAG_SAMPLES = int(
    float(os.environ.get("PANSY_TX_XPHASE_METADATA_LAG_SECONDS", "5")) * 1e6
)
MAX_WINDOWS_PER_PASS = int(os.environ.get("PANSY_TX_XPHASE_MAX_WINDOWS_PER_PASS", "240"))


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


def mode_id(value):
    values = n.asarray(value).reshape(-1)
    if values.size == 0:
        return None
    try:
        return int(values[0])
    except (TypeError, ValueError):
        return None


def tx_xphase():
    mddir=pc.tx_metadata_dir
    dm = drf.DigitalMetadataReader(mddir)
    rdir=pc.raw_voltage_dir
    d = drf.DigitalRFReader(rdir)


    # setup the directory and file cadence.
    # use 1 MHz, as this is the sample-rate and thus a
    # natural resolution for timing.
    subdirectory_cadence_seconds = 3600
    file_cadence_seconds = 60
    samples_per_second_numerator = 1000000
    samples_per_second_denominator = 1
    file_name = "txphase"
    os.system("mkdir -p %s"%(pc.phase_metadata_dir))

    dmw = drf.DigitalMetadataWriter(
        pc.phase_metadata_dir,
        subdirectory_cadence_seconds,
        file_cadence_seconds,
        samples_per_second_numerator,
        samples_per_second_denominator,
        file_name,
    )


    channels=["ch000","ch001","ch002","ch003","ch004","ch005","ch006","ch007"]
    raw_bounds = [d.get_bounds(channel) for channel in channels]
    b = (max(bound[0] for bound in raw_bounds), min(bound[1] for bound in raw_bounds))
    tx_bounds = dm.get_bounds()
    dt=60000000

    start_idx=load_progress()
    if start_idx is None or start_idx < b[0] or start_idx > b[1]:
        start_idx=b[0]

        try:
            dmp = drf.DigitalMetadataReader(pc.phase_metadata_dir)
            dmp_b=dmp.get_bounds()
            start_idx=dmp_b[1]
        except Exception:
            traceback.print_exc()

    process_end = min(b[1] - 120, tx_bounds[1] - TX_METADATA_LAG_SAMPLES)
    n_block=max(0, int(n.floor((process_end-start_idx)/dt)))
    n_block=min(n_block,MAX_WINDOWS_PER_PASS)

    xphase=n.zeros(8,dtype=n.complex64)
    for bi in range(n_block):
        print(stuffr.unix2datestr((start_idx+bi*dt)/1e6))

        data=dm.read(start_idx+bi*dt,start_idx+(bi+1)*dt-1,"id")
        mesosphere_keys=sorted(k for k, value in data.items() if mode_id(value) == 1)
        if mesosphere_keys:
            k=mesosphere_keys[0]
            if b[0] <= k and k + 120 <= b[1]:
                try:
                    z0=d.read_vector_c81d(k,120,"ch000")
                    for j in range(8):                    
                        z1=d.read_vector_c81d(k,120,channels[j])
                        xphase[j]=n.mean(z0*n.conj(z1))
                        if False:
                            #plt.subplot(121)
                            plt.plot((z0*n.conj(z1)).real)
                            plt.plot((z0*n.conj(z1)).imag)

                            #plt.plot(z0.imag)
                            #plt.subplot(122)
                            #plt.plot(z1.real)
                            #plt.plot(z1.imag)
                            plt.show()
                    if n.mean(n.abs(xphase))>1e7:
                        print("%s %1.1f %1.1f %1.1f %1.1f %1.1f %1.1f %1.1f %1.1f"%(stuffr.unix2datestr(k/1e6),n.angle(xphase[0]),n.angle(xphase[1]),n.angle(xphase[2]),n.angle(xphase[3]),n.angle(xphase[4]),n.angle(xphase[5]),n.angle(xphase[6]),n.angle(xphase[7])))
                        dout={"xphase":xphase}
                        dmw.write(k,dout)
                    else:
                        print("%s %1.1g"%(stuffr.unix2datestr(k/1e6),n.mean(n.abs(xphase))))

                except:
                    pass
#                    print()
#                    import traceback
#                    traceback.print_exc()
    if n_block > 0:
        save_progress(start_idx+n_block*dt)
    return n_block
            


if __name__ == "__main__":
    sleep_seconds = float(os.environ.get("PANSY_TX_XPHASE_SLEEP_SECONDS", "60"))
    while True:
        n_processed = tx_xphase()
        if n_processed < MAX_WINDOWS_PER_PASS:
            time.sleep(sleep_seconds)
