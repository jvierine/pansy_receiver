"""Publish numbers, not Matplotlib images, from the receiver computer.

This replaces only the status publisher. It never controls acquisition, writes
science metadata, or reads the sparse complex-voltage archive.
"""
import argparse
from datetime import datetime, timedelta, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import h5py
import numpy as np
from monitor.health import phase_health

SERVICES = ["pansy-uhd-rx", "pansy-uhd-watchdog", "pansy-delete-old-rf",
            "pansy-find-mode-starts", "pansy-status-plot", "pansy-mesomode-boundary",
            "pansy-cluster-mf", "pansy-quick-search-meteor", "pansy-pmse-sparse",
            "pansy-tx-xphase", "pansy-cut-meteors", "pansy-process-cut-meteor"]
CHANNELS = {"TX": "tx", "Matched filter": "mf", "Detections": "detections",
            "Cuts": "cut", "Modes": "mesomode", "Fits": "simple_meteor_fit"}


def files_between(root, start_s, stop_s):
    day = datetime.fromtimestamp(start_s, timezone.utc).date()
    end = datetime.fromtimestamp(stop_s, timezone.utc).date()
    while day <= end:
        yield from sorted(Path(root).glob(day.isoformat() + "T*/*.h5"))
        day += timedelta(days=1)


def metadata_rows(root, start_s, stop_s, fields):
    """Read only selected existing datasets; skip a concurrently written row."""
    for path in files_between(root, start_s, stop_s):
        try:
            with h5py.File(path, "r") as handle:
                for key in sorted((k for k in handle if k.isdecimal()), key=int):
                    if not start_s * 1e6 <= int(key) <= stop_s * 1e6:
                        continue
                    group = handle[key]
                    if not all(field in group for field in fields):
                        continue
                    yield int(key), {field: group[field][()] for field in fields}
        except OSError as exc:
            print("Metadata temporarily unreadable:", path.name, exc, flush=True)


def latest_sample(root):
    # The small metadata keys are enough; never load event waveforms here.
    for directory in sorted(Path(root).glob("20*T*"), reverse=True):
        for path in sorted(directory.glob("*.h5"), reverse=True):
            try:
                with h5py.File(path, "r") as handle:
                    keys = [int(k) for k in handle if k.isdecimal()]
                    if keys:
                        return max(keys) / 1e6
            except OSError:
                continue
    return None


def finite(value, digits=5):
    array = np.asarray(value)
    if array.ndim:
        return [finite(v, digits) for v in array]
    number = float(array)
    return round(number, digits) if np.isfinite(number) else None


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(json.dumps(value, separators=(",", ":"), allow_nan=False))
    temporary.replace(path)


def operations():
    result = subprocess.run(["systemctl", "--user", "show", *[s + ".service" for s in SERVICES],
                             "-p", "Id", "-p", "ActiveState", "-p", "SubState", "-p", "MainPID"],
                            capture_output=True, text=True, timeout=20, check=True)
    services = []
    for block in result.stdout.strip().split("\n\n"):
        row = dict(line.split("=", 1) for line in block.splitlines() if "=" in line)
        services.append(dict(name=row.get("Id", "unknown"), active=row.get("ActiveState") == "active",
                             active_state=row.get("ActiveState", "unknown"), sub_state=row.get("SubState"),
                             main_pid=int(row.get("MainPID", "0"))))
    disks = []
    for path in ("/media/analysis", "/media/old_analysis", "/media/archive"):
        try:
            usage = shutil.disk_usage(path)
            used = 100 * usage.used / usage.total
            disks.append(dict(path=path, total_bytes=usage.total, used_bytes=usage.used,
                              free_bytes=usage.free, used_percent=used, warning=used >= 80, critical=used >= 95))
        except OSError as exc:
            disks.append(dict(path=path, error=str(exc), warning=True))
    return dict(generated_utc=datetime.now(timezone.utc).isoformat(), services=services, disks=disks,
                phase_health=phase_health(Path(__file__).resolve().parents[1] / 'logs/pansy_uhd_watchdog.log'))


def packet(metadata, raw, now=None, reader_factory=None):
    now = time.time() if now is None else now
    metadata = Path(metadata)
    output = dict(schema="pansy_receiver_monitor_v1", generated_unix_s=now,
                  freshness={name: latest_sample(metadata / channel) for name, channel in CHANNELS.items()},
                  fits=[], phase=[], raw=None, errors={})
    fields = ("v0", "r0", "eclon", "eclat")
    latest = None
    for key, row in metadata_rows(metadata / "simple_meteor_fit", now - 24 * 3600, now, fields):
        output["fits"].append([key / 1e6, finite(np.linalg.norm(row["v0"])),
                               finite(np.linalg.norm(row["r0"])), finite(row["eclon"]), finite(row["eclat"])])
        if latest is None or key > latest:
            latest = key
    if latest is not None:
        for key, row in metadata_rows(metadata / "simple_meteor_fit", latest / 1e6, latest / 1e6,
                                       ("txidx", "ew", "ns", "up", "snr")):
            output["meteor"] = dict(time_unix_s=key / 1e6,
                seconds=finite((row["txidx"] - row["txidx"][0]) / 1e6, 7),
                east_km=finite(row["ew"]), north_km=finite(row["ns"]), up_km=finite(row["up"]),
                snr_db=finite(10 * np.log10(np.maximum(row["snr"], np.finfo(float).tiny))))
        # Read saved fit parameters/residuals, without recomputing the fit.
        for key,row in metadata_rows(metadata/'simple_meteor_fit',latest/1e6,latest/1e6,
                                     ('r0','v0','std','mfs')):
            if 'meteor' in output:
                output['meteor'].update(r0_km=finite(row['r0']),v0_km_s=finite(row['v0']),
                                         residual_std_km=finite(row['std']),match_function=finite(row['mfs']/21))
        # Associate an existing upstream detection trace by exact TX indices.
        # This reads saved scalars only, never cut waveforms or matched filters.
        for key,row in metadata_rows(metadata/'cut',latest/1e6-10,latest/1e6+1,
                                     ('c_tx_idx','c_range_km','c_snr','c_doppler')):
            if 'meteor' not in output:break
            times=np.asarray(row['c_tx_idx']).ravel()
            fit_times=np.rint(np.asarray(output['meteor']['seconds'])*1e6).astype(np.int64)+latest
            if not np.intersect1d(times,fit_times).size:continue
            output['meteor']['detection_trace']=dict(seconds=finite((times-latest)/1e6,7),
                range_km=finite(np.asarray(row['c_range_km']).ravel()),
                snr_db=finite(10*np.log10(np.maximum(np.asarray(row['c_snr']).ravel(),np.finfo(float).tiny))),
                doppler_km_s=finite(np.asarray(row['c_doppler']).ravel()/1000))
            break
        event=Path.home()/'.local/state/pansy-monitor/latest-event.h5'
        if event.exists() and 'meteor' in output:
            try:
                with h5py.File(event,'r') as h:
                    if int(h['key'][()])==latest:
                        nr,nt=h['power_snr'].shape
                        rs=max(1,(nr+127)//128);ts=max(1,(nt+255)//256)
                        output['meteor']['range_time']=dict(
                            seconds=finite(h['time_unix_s'][::ts]-latest/1e6,7),
                            range_km=finite(h['range_km'][::rs]),
                            display_time_stride=ts,display_range_stride=rs,
                            snr_db=finite(10*np.log10(np.maximum(h['power_snr'][::rs,::ts],np.finfo(float).tiny)),2))
                        output['meteor']['doppler_km_s']=finite(h['doppler_mps'][()]/1000)
            except OSError as exc:output['errors']['meteor_diagnostics']=str(exc)
    previous_bucket = None
    for key, row in metadata_rows(metadata / "phase", now - 24 * 3600, now, ("xphase",)):
        bucket = key // 600_000_000
        if bucket == previous_bucket:
            continue
        output["phase"].append([key / 1e6, *finite(np.angle(row["xphase"][:8], deg=True))])
        previous_bucket = bucket
    try:
        if reader_factory is None:
            from digital_rf import DigitalRFReader
            reader_factory = DigitalRFReader
        reader = reader_factory(str(raw))
        bounds = reader.get_bounds("ch000")
        raw_read_time = time.time()
        output["freshness"]["Raw"] = bounds[1] / 1e6
        # Display-only stride, no averaging or acquisition/normalization changes.
        begin = bounds[1] - 1_000_000
        channels = []
        for channel in range(8):
            voltage = reader.read_vector_1d(begin, 3200, "ch%03d" % channel)[::8]
            channels.append(dict(real=finite(voltage.real, 2), imag=finite(voltage.imag, 2)))
        output["raw"] = dict(start_unix_s=begin / 1e6, read_unix_s=raw_read_time, sample_step_us=8, channels=channels,
                             units="raw ADC units; display stride 8, no receiver equalization")
    except Exception as exc:
        output["errors"]["raw"] = str(exc)
    return output


def publish(metadata, raw, output, destination=None):
    output = Path(output)
    ops = output / "ops_status.json"
    atomic_json(ops, operations())
    def transfer(path):
        subprocess.run(["rsync", "-az", "--timeout=120", "--bwlimit",
                        os.environ.get("PANSY_STATUS_RSYNC_BWLIMIT", "1"), str(path), destination],
                       check=True, timeout=600)
    if destination:
        transfer(ops)
    data = output / "monitor-data.json"
    telemetry = packet(metadata, raw)
    atomic_json(data, telemetry)
    # Preserve the legacy publisher's operational clock-offset measurement.
    # This is telemetry only, not a write to science/acquisition products.
    if telemetry.get("raw") is not None:
        try:
            from digital_rf import DigitalMetadataWriter
            clock = Path(metadata) / "clock"
            clock.mkdir(parents=True, exist_ok=True)
            writer = DigitalMetadataWriter(str(clock), 3600, 60, 1_000_000, 1, "clock")
            now = telemetry["raw"]["read_unix_s"]
            writer.write(int(now * 1e6), dict(pc_tnow=now, raw_b1=telemetry["freshness"]["Raw"]))
        except Exception as exc:
            print("Clock telemetry delayed:", exc, flush=True)
    if destination:
        # Small service status first; keep the established satellite-link cap.
        transfer(data)
    print("Numeric receiver telemetry published", data.stat().st_size, "bytes", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, default=Path("/media/analysis/metadata"))
    parser.add_argument("--raw", type=Path, default=Path("/media/archive"))
    parser.add_argument("--output", type=Path, default=Path("/home/radar/.local/state/pansy-monitor"))
    parser.add_argument("--destination", default="j@juha.no:/var/www/html/pansy/")
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    while True:
        try:
            publish(args.metadata, args.raw, args.output, args.destination)
        except Exception:
            import traceback
            traceback.print_exc()
            if args.once:
                raise
        if args.once:
            break
        time.sleep(300)


if __name__ == "__main__":
    main()
