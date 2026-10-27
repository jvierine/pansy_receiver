# PANSY receiver user services

This directory contains the user-space systemd services for running the PANSY
receiver and processing chain on `syowa-meteor`.

The service set includes:

- `pansy-uhd-rx.service`: runs `pansy_uhd_rx`
- `pansy-uhd-watchdog.service`: restarts the receiver if Digital RF output or
  phase metadata stalls or jumps
- `pansy-delete-old-rf.service`: deletes old `rf*.h5` files under
  `/media/archive/ch*`
- `pansy-find-mode-starts.service`: runs `find_mode_starts.py`
- `pansy-status-plot.service`: runs `status_plot.py`
- `pansy-mesomode-boundary.service`: runs `mesomode_boundary.py`
- `pansy-cluster-mf.service`: runs `cluster_mf.py`
- `pansy-quick-search-meteor.service`: runs
  `mpirun -np 6 python3 quick_search_meteor.py`
- `pansy-meso-xc.service`: retired XC2 producer; disabled by default. Historical
  XC2 data and its readers remain available. An intentional installer rollback
  requires `PANSY_ENABLE_LEGACY_XC=1`; routine installation must not re-enable it.
  The live spectra producer is the separately managed `pansy-pmse-sparse.service`.
- `pansy-tx-xphase.service`: runs `tx_xphase.py`
- `pansy-cut-meteors.service`: runs `cut_meteors.py`
- `pansy-process-cut-meteor.service`: runs
  `mpirun -np 2 python3 process_cut_meteor.py`

The status monitor reads completed compressed-spectra metadata, including empty
noise-only modes, instead of using XC2 for its live PMSE panel or processing lag.
Its archive path can be overridden with `PANSY_PMSE_SPARSE_METADATA`. The panel
shows only retained zenith detections (peak SNR in dB); blank cells are not zero
power. Missing cells and NaNs use the colorbar's zero-level color for display
only; no zero-valued measurements are inserted. It does not read complex-voltage
arrays or change acquisition. Do not run
the full installer just to update this panel: restart only the status service.

Install or update on the receiver computer with:

```bash
cd /home/radar/src/git/pansy_receiver
bash receiver/install_user_service.sh
```

The installer writes the default configuration to:

```bash
~/.config/pansy-receiver/pansy-receiver.env
```

The installer also writes log rotation config to:

```bash
/etc/logrotate.d/pansy-receiver
```

The watchdog checks `ch000` through `ch007` under `/media/archive` every
10 seconds. One `pansy_uhd_rx` process owns all four USRPs and all eight
channels. If any channel is stale, the watchdog stops that combined receiver,
verifies that it has exited, waits 15 seconds for every UHD device session to
close, and then starts all receivers together. The receiver binary also treats
a stall in any individual stream as a failure of the combined receiver,
allowing systemd to reopen every
USRP and channel together instead of restarting one channel with an unrelated
phase. It measures transmit-pulse
cross-phase only at sample indices identified as mesosphere mode (`id=1`) in
the TX metadata. A phase jump larger than the configured threshold triggers
the same full receiver restart. Every
watchdog restart is recorded in
`~/.local/state/pansy-receiver/receiver_restart.json`; phase validation waits
90 seconds and only accepts a mesosphere-mode transmit pulse whose TX metadata
and raw voltage were both captured after that timestamp.
If the initial three phase-recovery attempts fail, the watchdog continues with
one full receiver restart every 10 minutes until a fresh TX pulse matches the
known-good phase pattern. It no longer leaves a bad phase state running
indefinitely after exhausting the initial retry burst.

The mode finder keeps its independent raw-voltage scan cursor in
`~/.local/state/pansy-receiver/find_mode_starts.json`. The cursor advances for
every completed window, including windows with no recognized transmit mode,
so transmitter-off intervals do not create an ever-growing rescan backlog.
The cross-phase monitor similarly tracks progress in
`~/.local/state/pansy-receiver/tx_xphase.json`, stays behind the TX metadata
writer, and computes one cross-phase sample only from a confirmed mesosphere
mode (`id=1`) pulse.

Useful commands:

```bash
systemctl --user --no-pager --type=service --state=running | grep pansy
systemctl --user status pansy-uhd-rx.service
systemctl --user status pansy-uhd-watchdog.service
journalctl --user -u pansy-quick-search-meteor.service -f
journalctl --user -u pansy-uhd-rx.service -f
journalctl --user -u pansy-uhd-watchdog.service -f
ls -lh /home/radar/src/git/pansy_receiver/logs
tail -f /home/radar/src/git/pansy_receiver/logs/pansy_uhd_rx.service.log
tail -f /home/radar/src/git/pansy_receiver/logs/pansy_uhd_watchdog.log
tail -f /home/radar/src/git/pansy_receiver/logs/pansy_quick_search_meteor.log
```
