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

Deploy the web page once:

```bash
PANSY_RECEIVER_REPO=~/src/pansy_receiver \
  ~/src/pansy_receiver/backup/scripts/pansy_deploy_web.sh
```
