# PANSY metadata backup on revontuli

This directory contains user-space systemd services for mirroring PANSY metadata
from the Antarctica receiver to `revontuli.uit.no`.

The backup is pull-only. It never uses `rsync --delete` and never removes files
from the Antarctica receiver.

## Install on revontuli

Clone or update this repository on `j@revontuli.uit.no`, then run:

```bash
cd ~/src/pansy_receiver
bash backup/install_user_services.sh
```

If the service should run while user `j` is logged out, enable lingering once:

```bash
sudo loginctl enable-linger j
```

The installer copies the default configuration to:

```bash
~/.config/pansy-backup/pansy-backup.env
```

Edit that file if paths or ports need to change.

Transfers share the configured `RSYNC_BWLIMIT_KB` budget (128 KiB/s by
default); `RSYNC_TOTAL_BWLIMIT_KB` overrides it. Cutouts and TX timing receive
larger shares than the smaller channels. Each channel retries independently
every `RSYNC_SLEEP_SECONDS`, so a cutout backlog cannot prevent a failed phase
or detection transfer from retrying. Per-channel results update
`backup_state/rsync_status.tsv` immediately. All phase history is synchronized,
including gaps older than 72 hours. Transfers include channel properties files.

The budget is shared, not applied independently to every channel. With the seven
default channels, cutouts receive 64 KiB/s, TX receives 26 KiB/s and each other
channel receives 7 KiB/s (125 KiB/s total after rounding). SSH/protocol overhead
is additional; this is not a host-wide network shaper and does not cap RVD.
Each channel runs ordinary `rsync -avz` until completion, without a wall-clock
cutoff or recent/history scheduling. The 120-second rsync timeout detects a
stalled connection; it does not interrupt a transfer that is making progress.
Interrupted files are retained in `.rsync-partial` for reuse on retry. While
catching up, newer files may wait behind older missing files in that channel.
To configure direct pulls with the default cap (configuration is
backed up; receiver services are never touched):

```bash
python3 backup/scripts/configure_direct_rsync.py \
  --config ~/.config/pansy-backup/pansy-backup.env \
  --unit-dir ~/.config/systemd/user --bwlimit-kb 128
systemctl --user daemon-reload
systemctl --user restart pansy-backup-rsync.service
```

## Transfer safety and interpreting status

- There is one worker set per state directory, enforced with a nonblocking lock;
  each channel also has its own lock. Detached SSH masters do not inherit them.
- Neither the source pull nor the local mirror deletes source or destination
  science files. Both keep interrupted files in hidden `.rsync-partial`
  directories; the local mirror excludes source partial directories.
- Set `LOCAL_METADATA_MOUNT` and `LOCAL_MIRROR_METADATA_MOUNT` to the actual
  storage mountpoints. Missing mounts stop the affected service rather than
  silently filling the root filesystem. Source and mirror must not overlap.
- A successful recent-only legacy pass does not verify the historical archive.
  New status rows explicitly identify a full-channel pass. While such a pass is
  running, archive completeness is not inferred from a recent live file.
- Freshness comes from actual HDF5 observation timestamps, not hour-directory
  mtimes or empty directory names. The check is bounded to the latest 48 hour
  directories and up to three files per directory. Without h5py, filename
  interval starts are used and labeled as such. Full archive counts/byte totals
  are explicitly unknown, not fabricated from directory sizes.
- Status JSON is replaced atomically and includes mount availability and free
  space. A successful local mirror only verifies the Revontuli-to-mirror copy,
  not the completeness of the original Syowa transfer.

To apply guards and point only the three transfer/status services at a checked
checkout, pass the following options to `configure_direct_rsync.py`:

```bash
--metadata-mount /mnt/data --mirror-mount /urdr/data --service-repo /path/to/checkout
```

The configurator backs up the environment file and never restarts services.
After reviewing its changes, reload systemd and restart just those services.
RVD sparse-spectra pulls have a separate bandwidth cap; the 128 KiB/s metadata
budget is not a cap on every network service on either host.

## Services

```bash
systemctl --user status pansy-backup-tunnel.service
systemctl --user status pansy-backup-rsync.service
systemctl --user status pansy-local-metadata-mirror.service
systemctl --user status pansy-detection-history.service
systemctl --user status pansy-backup-web.service
```

The services do the following:

- `pansy-backup-tunnel.service`: keeps the SSH tunnel open through
  `j@4.235.86.214` to `radar@localhost -p 3131`.
- `pansy-backup-rsync.service`: pulls metadata channels through
  `localhost:2222` into `/mnt/data/juha/pansy/metadata`.
- `pansy-local-metadata-mirror.service`: mirrors
  `/mnt/data/juha/pansy/metadata/` to `/urdr/data/juha/pansy/metadata/`
  without `--delete`.
- `pansy-detection-history.service`: counts the `detections` metadata channel
  by UTC day and writes `detections_daily.png`.
- `pansy-backup-web.service`: writes `status.json`, then
  deploys the web page to `j@4.235.86.214:/var/www/html/pansy`.

## Manual commands

Start or restart everything:

```bash
systemctl --user restart pansy-backup-tunnel.service
systemctl --user restart pansy-backup-rsync.service
systemctl --user restart pansy-local-metadata-mirror.service
systemctl --user restart pansy-detection-history.service
systemctl --user restart pansy-backup-web.service
```

Inspect logs:

```bash
journalctl --user -u pansy-backup-rsync.service -f
tail -f /mnt/data/juha/pansy/backup_state/rsync.log
tail -f /mnt/data/juha/pansy/backup_state/local_mirror.log
```

Deploy dynamic web status JSON and plot images once:

```bash
PANSY_RECEIVER_REPO=~/src/pansy_receiver \
  ~/src/pansy_receiver/backup/scripts/pansy_deploy_web.sh
```

Deploy the static web page and logo/favicon assets only after editing them:

```bash
PANSY_RECEIVER_REPO=~/src/pansy_receiver \
  ~/src/pansy_receiver/backup/scripts/pansy_deploy_web.sh --include-static
```
