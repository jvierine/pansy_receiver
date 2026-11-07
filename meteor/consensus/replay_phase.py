"""Bounded, read-only comparison of legacy and phase-aware RANSAC grouping."""
import argparse
from pathlib import Path
import sys
import time

import digital_rf
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import cluster_mf
from meteor_phase_grouping import PhaseGrouper


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata")
    parser.add_argument("raw")
    parser.add_argument("start",type=int)
    parser.add_argument("end",type=int)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    if not 0 < args.end-args.start <= 10_000_000 or args.output.exists():
        parser.error("require a new output and a window <=10 seconds")
    reader = digital_rf.DigitalMetadataReader(args.metadata)
    txp,t,r,v,s,b = cluster_mf.read_mf_output(reader,args.start,args.end)
    keep = (t >= args.start) & (t < args.end)
    txp,t,r,v,s,b = [x[keep] for x in (txp,t,r,v,s,b)]
    class Capture:
        def __init__(self): self.records=[]
        def write(self,keys,data): self.records.append(data)
    captured = Capture()
    cluster_mf.find_clusters(txp,t,r,v,s,b,captured)
    legacy = np.full(len(t),-1,dtype=np.int32)
    for i,record in enumerate(captured.records):
        legacy[np.isin(t,record["tx_idx"])] = i
    grouper = PhaseGrouper(args.raw)
    began = time.monotonic()
    try:
        labels, signatures = grouper.group(t,r,v,s,b)
    finally:
        grouper.close()
    with h5py.File(args.output,"x") as out:
        for key,value in dict(sample_idx=t.astype(np.int64),range_km=r.astype(np.float32),
                              doppler_m_s=v.astype(np.float32),snr_db=(10*np.log10(s)).astype(np.float32),
                              beam_id=(b%5).astype(np.int8),legacy_track_id=legacy,
                              track_id=labels,receiver_amplitudes=signatures).items():
            out[key]=value
        out.attrs["elapsed_s"]=np.float32(time.monotonic()-began)
        out.attrs["phase_coherence_min"]=np.float32(.65)
        out.attrs["max_phase_link_lag_s"]=np.float32(.017)
    print(f"Legacy tracks={len(captured.records)} points={np.count_nonzero(legacy>=0)}; "
          f"new tracks={len(np.unique(labels[labels>=0]))} points={np.count_nonzero(labels>=0)}",flush=True)


if __name__ == "__main__":
    main()
