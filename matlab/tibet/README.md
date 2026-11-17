# Tibet Meteor Processing in MATLAB

Standalone MATLAB/Octave detector and ten-channel voltage cutter. No Rust,
Python, external executable, or MATLAB toolbox is needed to run it.

Open `run_tibet_meteors.m`, set `InputFolder` and `OutputFolder`, and run it
from this folder. Windows paths are supported. Completed, unchanged files are
skipped after restarting. Files must be stable before processing. Failed or
partial outputs are preserved and never marked complete.

See [Instructions for Yufang Tian](OPERATIONS.md) for Windows deployment,
parallel workers, live-throughput checks, and opt-in verified raw-file cleanup.
New outputs are grouped by UTC `YYYY-MM-DD`, then by source-file identifier.
Each source has a two-panel overview; worker 1 maintains a rolling day overview.
The five-file example averages about 5.1 MiB retained per 900 MiB input,
including cuts, metadata and plots. The current development setup has **not**
demonstrated sustained real-time throughput.

The detector searches all five full-array receivers with measured-code
range-Doppler matched filters. Cuts contain only S0-S2 and L0-L6, a common
sampled transmitter reference, two earlier references for ambiguity tests,
original pulse indices, receiver configuration, and transmit pointing.
No phase calibration is applied to the stored voltages.

MATLAB writes HDF5 with byte shuffle and deflate level 4. Octave, which lacks
the same high-level HDF5 writer, writes compressed MAT v7 files. Set `Format`
to `mat` to use that format in MATLAB too. All main voltage arrays are
interleaved int16. Rare larger values are retained exactly in sparse overflow
tables; `tibet_unpack_iq` reconstructs the original int32 values. There is no
quantization or loss from clipping. In MATLAB the echo dimensions are
`[I/Q, range sample, receiver, pulse]`.

The input is the validated 500 kHz, 1280 us IPP, 64 us complementary-code mode.
Other modes are rejected. A nontrivial multi-beam mode requires a MAT file
containing a `pointing` structure with one value per source pulse: `alpha`,
`beta`, `beam_id`, `zenith_deg`, and `azimuth_deg`. Set `PointingFolder` to its
folder; each schedule is named `<raw filename stem>_pointing.mat`. Unknown
directions are not guessed. The included geometry retains all 931 positions;
active/dummy membership and geographic phase orientation remain provisional.

Binary calendar fields in the current campaign are invalid. Filename dates
are interpreted as Beijing time, explicitly flagged as unconfirmed. The
search threshold is a detector statistic, not a calibrated SNR. Even-IPP
range aliases remain unresolved despite the earlier-code comparison.

Validation: the standalone Octave detector processed all five original files
(30,720 IPPs), reproducing all 16 native events and 1,063 detection points.
Pulse indices, range/Doppler bins, all 3,594 candidates and 1,228 code-alias
vetoes match. Interpolated Doppler differs by less than 0.01 Hz. All
12,056,232 saved I/Q components are bit-for-bit identical to the original
native cuts, and Rust replay checks all 10,630 channel amplitudes. A repeat
folder scan skips all five completed files. Windows MATLAB execution and
its native HDF5-writing branch have not been tested on a Windows machine.
