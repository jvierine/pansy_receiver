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
- `pansy-meso-xc.service`: runs `meso_xc.py`
- `pansy-tx-xphase.service`: runs `tx_xphase.py`
- `pansy-cut-meteors.service`: runs `cut_meteors.py`
- `pansy-process-cut-meteor.service`: runs
  `mpirun -np 2 python3 process_cut_meteor.py`

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
10 seconds. If any channel is stale for more than 60 seconds, it stops
`pansy-uhd-rx.service`, waits for `pansy_uhd_rx` to exit, and starts the
service again. After a data-stall restart it also waits for the next phase
metadata sample in `/media/analysis/metadata/phase` and compares it with the
recent circular median. A phase jump larger than the configured threshold
causes another clean receiver restart.

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
