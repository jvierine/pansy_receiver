# Tibet MST: Instructions for Yufang Tian

## First Run in Windows MATLAB

1. Extract the package into a permanent local folder. Keep every supplied `.m`
   file and `array_geometry.mat` together. Make this MATLAB's Current Folder.
2. Open `run_tibet_meteors.m`. Set `InputFolder` to the acquisition folder and
   `OutputFolder` to the retained-event folder, preferably on a different fast
   volume. Examples are `D:\TibetMST\raw` and `E:\TibetMST\meteors`.
3. Leave `DeleteRawAfterSuccess = false` during commissioning. Set
   `RunOnce = true` for a batch of complete pilot files and press Run.
4. Check that the terminal reports completed files, the event files reopen,
   and the overview figures show sensible detections. Compare selected events
   with the native reference analysis before committing a new experiment mode.
5. For live folder watching, set `RunOnce = false`. A source must remain
   unchanged for 20 seconds and all declared pulse markers/counters must pass.
   Polling is every 10 seconds. Acquisition should use unique filenames and
   should publish a file only after writing it; do not reuse old filenames.

No Python, Rust executable, or specialist MATLAB toolbox is needed at runtime.
Stop a worker with Ctrl+C. Acquisition is independent and is not stopped by this
program. Never disconnect its input/output drives during processing.

## Output Dates and Plots

Output is organized as:

```text
OutputFolder/
  YYYY-MM-DD/                    UTC day
    overview.png                 rolling overview for the day
    overview.pdf
    original_raw_filename_stem/
      catalogue.mat              completion record, tracks and compact power map
      overview.png               this source file
      overview.pdf
      event_00000.mat or .h5      ten-channel voltage cut and metadata
      retention.mat              present when raw cleanup was requested
```

The supplied recordings start on 14 November 2026 at 05:51:56 Beijing time, so
their output day is **2026-11-13 UTC**. This conversion assumes that filenames
use Beijing time; the binary calendar is invalid. UTC folders do not establish
absolute clock calibration. Files crossing midnight belong to their source's
start day.

The upper plot is matched-filter peak power divided by median power, taking
the maximum over Doppler. It is not calibrated SNR. A profile is saved every
16 IPPs; every IPP is still searched by the detector. The lower panel contains
retained candidate tracks, coloured by track. File gaps remain visible.
Per-source plots are permanent. Worker 1 refreshes a rolling day overview every
60 seconds, using the latest 120 completed files from each of the latest two
output days. The entire day is not loaded into RAM.

## Throughput and Independent Workers

Watching a folder is not proof that processing keeps up. This mode produces
roughly one source file every 8 seconds. The **combined** workers must finish
at least that quickly, with margin. Monitor file completion delay and the
number/age of waiting raw files over a sustained site pilot before enabling
cleanup.

The initial implementation needed 176 seconds in one worker and 108 seconds
in six workers for 39.3 seconds of recordings. Profiling found that Octave
fully sorted large arrays to obtain medians. Exact order-statistic selection
removed that bottleneck without changing the five-channel detector. The
optimized five-worker test needed 24.1 seconds, including per-file plots and
verified cuts, before the final shared overview refresh. All 16 events,
1,063 points and 12,056,232 saved integer components still match the reference.
This is faster than acquisition in a replay, not a sustained Windows site
acceptance test. A longer replay repeated the five recordings four times:
20 files, 157.3 seconds of recorded IPPs, processed in about one minute by
five workers, with compressed verified cuts and overview plots. Repeated
recordings benefit from filesystem caching; this is not a test of simultaneous
acquisition writes or a Windows MATLAB benchmark. More workers alone may not help when storage or CPU is
saturated. Use fast local acquisition/staging storage, sufficient
physical cores, and benchmark the actual Tibet computer with MATLAB. The decoder
reads about 60 MB/s of declared voltage records at acquisition rate, in addition
to acquisition writes and output I/O.

To try five independent workers, launch five dedicated MATLAB or Octave sessions with the
same input and output settings. In each session's Command Window, use a distinct
index and then run the launcher:

```matlab
WorkerCount = 5;
WorkerIndex = 1;  % 1, 2, 3, 4 or 5; one distinct index per session
DetectionMode = 'all';
run_tibet_meteors
```

All workers must use the same `WorkerCount`. The filename frame counter selects
one worker deterministically; fallback names use a stable name hash. Claim
directories protect against duplicate processing. One computation thread per
worker avoids thread oversubscription. Check MATLAB licensing and available RAM
before starting multiple sessions. No Parallel Computing Toolbox is required.
Do not run other analysis in a dedicated worker session.

## Optional Transmit-Matched Detection

`DetectionMode = 'transmit'` searches only the full-array receiver whose raw
steering codes match that pulse's transmit steering. For the present zenith
files this is MATLAB channel 3 (zero-based channel 2). An absent or ambiguous
match is an error, not a guess. Cuts still retain all ten interferometer
channels and the sampled transmit reference. Frame-by-frame scanning is read
from each header; only within-frame switching needs a per-pulse pointing file.

The default five-channel mode remains the reference detector. The one-channel
mode changes the detection statistic and uses a 20 dB peak/median-minus-one
threshold instead of 12 dB: the same numerical threshold would admit far more
single-channel noise peaks. Neither threshold is a calibrated false-alarm
probability. Compare event completeness and false alarms before adopting the
one-channel mode; do not assume its output is identical to the reference.
In the five-file test it recovered all 16 reference events, but did not retain
exactly the same pulse samples. Change `OutputFolder` when comparing modes;
the watcher refuses to reuse completed results from a different mode.
The full five-channel optimized mode already exceeds acquisition rate in the
five-file replay, so dropping channels is not required for that result.

## Retain Only Detected Voltages

After a successful site validation and a throughput pilot, set:

```matlab
DeleteRawAfterSuccess = true;
```

This is deliberately **off by default**. Before removing a source, the worker
requires full successful processing, atomically published output, saved-event
read-back verification, unchanged saved-event integrity checksums, unchanged
source size/modification time/header, and the requested per-source overview.
A retention record is written before and after removal. Failed files, changed
files, missing/changed cuts, unsupported modes and incomplete outputs retain the
raw source. A valid zero-detection source can be removed after its completion
record and diagnostics have been saved.

Only detected ten-channel echo voltages and sampled transmit history remain;
compact catalogues, diagnostic power profiles and plots are retained too. For
the five-file test, the average retained output is about **5.1 MiB per 900 MiB
source**, including diagnostics: about **175-fold smaller**, or **99.4% less
disk space**. Events alone average about **4.2 MiB per source**. These are sample
averages, not guaranteed limits; meteor rate, duration and compression change
the ratio. Roughly half of the source file capacity is undeclared padding,
which is not processed as observations.

Deletion prevents later searches for missed/weaker events. Preserve selected
calibration/audit raw files outside the watched folder. The validation and
retention tests did not delete any original Tibet recordings.
Read-back checks and atomic folder publication are not backups or a guarantee
against power-loss/write-cache failure. Protect the retained archive accordingly;
this code does not implement a cross-platform durable-storage transaction.

## Restarting and Errors

- Completed unchanged sources are skipped on restart; cleanup may finish on a
  restart if it was requested and all integrity checks still pass.
- A `.partial` source output is not complete. Keep its raw source. Stop the
  responsible worker, inspect the failure, and move the incomplete output aside
  before retrying. Never remove a claim directory while its worker is running.
- A stale `.lock` after a crash needs the same check: confirm no worker owns it
  before removing it. The raw source is never removed for an incomplete result.
- A graphics failure leaves valid cuts but prevents requested cleanup. Fix the
  graphics backend and regenerate the plots, rather than deleting the raw file
  manually. Daily plots can be rebuilt with
  `tibet_refresh_overviews(OutputFolder)`.
- A five-direction cycle of raw-file frames works automatically with the same
  coded pulse and IPP. `BeamTotal` is the number of frames in the cycle;
  `CurFrameNum` and transmit alpha/beta describe the current frame. No external
  schedule is needed for ordinary frame-by-frame scanning. Every cut records
  the frame beam ID, cycle count, steering codes and selected common sampled
  transmit reference. No geographic azimuth is guessed from unconfirmed axes.
- Only direction switching *within* one raw-file frame needs an explicit
  per-pulse schedule. Set `PointingFolder`; files are
  named `<raw stem>_pointing.mat` and contain the `pointing` fields described in
  the README. Unknown pointing and unvalidated waveform modes fail closed.
  The raw pulse records do not themselves contain a per-pulse direction code.

## What Octave Does

GNU Octave is a free MATLAB-like numerical interpreter, available on Windows,
macOS and Linux. It executes these same `.m` detector, cutting and plotting
functions directly. It is not a converter and does not call MATLAB, Python or
Rust to find the events. FFTs use its native FFT library, limited to one thread
per worker by default.

MATLAB's default output is HDF5 with shuffle + deflate level 4. The tested
Octave path uses compressed MATLAB v7 event files because it does not provide
the same high-level HDF5-writing API. Both store int16 interleaved I/Q and
sparse overflow words, recovering the original integers exactly. Octave needs
a working graphics backend for PNG/PDF output. This package was tested with
Octave 11.3 on macOS; Windows MATLAB and Windows Octave have not been executed
in the development tests. Test the actual site installation before cleanup.
