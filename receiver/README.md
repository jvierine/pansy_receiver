# PANSY UHD receiver user service

This directory contains the user-space systemd service for running
`pansy_uhd_rx` on `syowa-meteor`, plus a watchdog that restarts the receiver
when one or more Digital RF channels stop producing new HDF5 files.

Install or update on the receiver computer with:

```bash
cd /home/radar/src/git/pansy_receiver
bash receiver/install_user_service.sh
```

The installer writes the default configuration to:

```bash
~/.config/pansy-receiver/pansy-receiver.env
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
systemctl --user status pansy-uhd-rx.service
systemctl --user status pansy-uhd-watchdog.service
journalctl --user -u pansy-uhd-rx.service -f
journalctl --user -u pansy-uhd-watchdog.service -f
tail -f /home/radar/src/git/pansy_receiver/logs/pansy_uhd_rx.service.log
tail -f /home/radar/src/git/pansy_receiver/logs/pansy_uhd_watchdog.log
```
