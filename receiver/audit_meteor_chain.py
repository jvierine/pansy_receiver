"""Bounded read-only audit of completed meteor metadata and retained raw data."""
import argparse
import fcntl
import os
from pathlib import Path
import sys
import time
import h5py
import numpy as np
import digital_rf

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pansy_config as pc
from meteor_metadata import StrictDigitalMetadataReader as Reader
from meteor_progress import read_progress,search_checkpoint,inclusive_chunks


def gaps(blocks,start,end):
    cursor=start
    result=[]
    for first,length in sorted(blocks.items()):
        stop=min(end+1,int(first)+int(length))
        if stop<=cursor: continue
        if first>cursor: result.append((cursor,min(int(first),end+1)))
        cursor=max(cursor,stop)
    if cursor<=end: result.append((cursor,end+1))
    return result


def read(reader,a,b,columns=None,lock=None):
    for attempt in range(5):
        try:
            if lock:
                with open(lock,"a") as stream:
                    fcntl.flock(stream,fcntl.LOCK_SH)
                    return reader.read(a,b,columns)
            return reader.read(a,b,columns)
        except OSError:
            if attempt==4: raise
            time.sleep(.2)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--seconds",type=int,default=600)
    parser.add_argument("--output",required=True)
    parser.add_argument("--end-sample",type=int,help="inclusive historical audit end")
    args=parser.parse_args()
    if not 1<=args.seconds<=3600: parser.error("audit at most one hour per invocation")
    state=Path.home()/".local/state/pansy_receiver"
    checkpoints={name:read_progress(state/f"{name}_complete_v1.h5")
                 for name in ("cut","fit")}
    checkpoints["group"]=read_progress(state/"cluster_phase_v1.h5")
    checkpoints["search"]=read_progress(search_checkpoint())+1
    end=min(checkpoints.values())-1
    if args.end_sample is not None:
        if args.end_sample>end: parser.error("audit end is not completed by every stage")
        end=args.end_sample
    start=end-args.seconds*1_000_000+1
    tx=Reader(pc.tx_metadata_dir); mf=Reader(pc.mf_metadata_dir)
    det=Reader(pc.detections_metadata_dir); cut=Reader(pc.cut_metadata_dir)
    fit=Reader(pc.simple_fit_metadata_dir)
    raw=digital_rf.DigitalRFReader(pc.raw_voltage_dir)
    raw_bounds=[raw.get_bounds(f"ch{i:03d}") for i in range(8)]
    missing_mf=[]; bad_mf=[]; missing_cut=[]; bad_cut=[]; raw_gaps=[]
    counts=dict(tx_batches=0,mf_batches=0,detections=0,cuts=0,fits=0,cut_payloads_checked=0,low_snr_cuts=0)
    fields=["beam_pos_idx","tx_idxs","tx_pwr","max_snr","max_range","max_dopvel","noise_floor","tx_std"]
    for a,b in inclusive_chunks(start,end,60_000_000):
        expected={int(k) for k,v in read(tx,a,b,"id").items() if int(np.asarray(v).ravel()[0])==1}
        measured=read(mf,a,b,fields,"/tmp/pansy_mf_metadata.lock")
        missing_mf.extend(sorted(expected-set(measured)))
        counts["tx_batches"]+=len(expected); counts["mf_batches"]+=len(measured)
        for key,row in measured.items():
            if any(np.asarray(row[f]).size!=20 for f in fields): bad_mf.append(key)
            elif not np.array_equal(row["tx_idxs"],int(key)+np.arange(20)*1600): bad_mf.append(key)
        events=read(det,a,b,["tx_idx","beam"])
        cuts=read(cut,a,b,["c_tx_idx","c_snr"])
        for key,row in cuts.items():
            snr=np.asarray(row["c_snr"])
            if not snr.size: bad_cut.append(key)
            elif np.max(snr)<20: counts["low_snr_cuts"]+=1
        counts["detections"]+=len(events); counts["cuts"]+=len(cuts)
        counts["fits"]+=len(read(fit,a,b,"txidx"))
        for key,row in events.items():
            pulse=np.asarray(row["tx_idx"])
            eligible=max(np.bincount(np.asarray(row["beam"],dtype=int)%5,minlength=5))>=5
            retained=pulse.min()-97600>=max(v[0] for v in raw_bounds)
            if eligible and retained and key not in cuts: missing_cut.append(key)
        # Sample one full voltage cut per minute without retaining raw arrays.
        if cuts:
            key=next(iter(cuts))
            payload=read(cut,key,key)[key]
            n=len(payload["tx_idx"])
            if not n or any(len(payload[f])!=n for f in
                    ["beam_id","delays","ztx_pulses_re","ztx_pulses_im","zrx_echoes_re","zrx_echoes_im"]):
                bad_cut.append(key)
            counts["cut_payloads_checked"]+=1
    for ch,bounds in enumerate(raw_bounds):
        a=max(start,bounds[0]); b=min(end,bounds[1]-1)
        if a<=b:
            raw_gaps.extend((ch,lo,hi) for lo,hi in gaps(raw.get_continuous_blocks(a,b,f"ch{ch:03d}"),a,b))
    with h5py.File(args.output,"x") as out:
        out.attrs["created_unix_s"]=np.int64(time.time())
        out.attrs["start_sample"]=np.int64(start); out.attrs["end_sample_inclusive"]=np.int64(end)
        for name,value in counts.items(): out.attrs[name]=np.int32(value)
        for name,value in checkpoints.items(): out.attrs[f"{name}_next_sample"]=np.int64(value)
        for name,value in dict(missing_mf=missing_mf,bad_mf=bad_mf,missing_cut=missing_cut,bad_cut=bad_cut,raw_gaps=raw_gaps).items():
            out[name]=np.asarray(value,dtype=np.int64)
    print(counts,flush=True)
    print("Missing MF",len(missing_mf),"invalid MF",len(bad_mf),"missing cuts",len(missing_cut),
          "invalid sampled cuts",len(bad_cut),"raw gaps",len(raw_gaps),flush=True)
    if missing_mf or bad_mf or missing_cut or bad_cut or raw_gaps: sys.exit(1)


if __name__=="__main__": main()
