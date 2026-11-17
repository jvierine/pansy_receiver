"""Exercise independent Octave workers; this is test orchestration, not a runtime dependency."""

import argparse
from pathlib import Path
import subprocess
import time
import os
import h5py
from scipy.io import loadmat


def quote(value):
    return str(value).replace("'", "''")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("raw", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--workers", type=int, default=6)
    args = p.parse_args()
    assert not args.output.exists()
    args.output.mkdir(parents=True)
    logs = args.output.parent / (args.output.name + "_worker_logs")
    logs.mkdir()
    code = Path(__file__).resolve().parents[1]
    jobs = []
    started = time.monotonic()
    failures = []
    environment = dict(
        os.environ,
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
    )
    for worker in range(1, args.workers + 1):
        command = (
            f"addpath('{quote(code)}'); s=tibet_watch('{quote(args.raw)}','{quote(args.output)}',"
            f"struct('RunOnce',true,'StableSeconds',0,'Format','mat','WorkerCount',{args.workers},"
            f"'WorkerIndex',{worker},'DeleteRawAfterSuccess',false)); "
            f"save('-mat7-binary','{quote(logs/f'worker_{worker}.mat')}','s'); assert(s.failed==0);"
        )
        handle = (logs / f"worker_{worker}.log").open("w")
        process = subprocess.Popen(
            ["octave", "--no-gui", "--quiet", "--eval", command],
            stdout=handle,
            stderr=subprocess.STDOUT,
            env=environment,
        )
        jobs.append((worker, process, handle))
    results = []
    while jobs:
        for worker, process, handle in jobs[:]:
            if process.poll() is not None:
                handle.close()
                if process.returncode != 0:
                    failures.append(worker)
                    jobs.remove((worker, process, handle))
                    continue
                value = loadmat(logs / f"worker_{worker}.mat", simplify_cells=True)["s"]
                results.append(
                    [
                        worker,
                        value["processed"],
                        value["skipped_complete"],
                        value["failed"],
                        value["elapsed_seconds"],
                    ]
                )
                print("worker", worker, value, flush=True)
                jobs.remove((worker, process, handle))
        if jobs:
            time.sleep(1)
    wall = time.monotonic() - started
    assert not failures, (failures, logs)
    assert sum(row[1] for row in results) == len(list(args.raw.glob("20*.raw")))
    with h5py.File(logs / "benchmark.h5", "w") as h:
        h["workers"] = sorted(results)
        h["wall_seconds"] = wall
        h.attrs["columns"] = "worker index, processed, skipped, failed, worker seconds"
        h.attrs["scope"] = (
            "Cold-start batch throughput on development hardware; not sustained site real-time validation. No original raw files deleted."
        )
    subprocess.run(
        [
            "octave",
            "--no-gui",
            "--quiet",
            "--eval",
            f"addpath('{quote(code)}'); tibet_refresh_overviews('{quote(args.output)}');",
        ],
        check=True,
    )
    print("PASS", args.workers, "workers; wall seconds", wall, flush=True)


if __name__ == "__main__":
    main()
