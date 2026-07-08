# PANSY Receiver Agent Notes

Date: 2026-06-03

Local repo: `/Users/jvi019/src/pansy_receiver`
Remote repo: `/home/radar/src/git/pansy_receiver` on `syowa-meteor`
Remote SSH: `ssh -J j@4.235.86.214 -p 3131 radar@localhost`

Analysis local configuration:
- Use `/Users/jvi019/src/pansy_receiver` with `conda run -n base python ...` locally.
- Runtime analysis overrides can be placed in `~/.config/pansy-receiver/pansy-analysis.env`, or another file pointed to by `PANSY_ANALYSIS_CONFIG`.
- `DASST_KERNEL` or `DASST_KERNEL_PATH` in that config controls the local `de430.bsp` path. The process environment variable `DASST_KERNEL` takes precedence.
- If the configured DASST kernel file is missing, the orbit code downloads `de430.bsp` from the NASA NAIF generic kernels site into the configured path. Do not manually copy this file between machines.

Latest deployed commit: `1a4d0be Publish process and disk status`

What changed:
- Added user systemd services for the receiver processing chain under `receiver/systemd/`.
- Added shared service launchers in `receiver/scripts/`.
- Added log rotation template at `receiver/logrotate/pansy-receiver`.
- Installer is `receiver/install_user_service.sh` and supports noninteractive deploy with `PANSY_SUDO_PASSWORD`.
- Logs append under `/home/radar/src/git/pansy_receiver/logs`.

Remote mount state:
- `/dev/sdd1` UUID `ba610c5a-cfe8-4b7c-a9d6-ea1007eb3907` is mounted at `/media/analysis`.
- `/dev/sdb1` UUID `52c5062b-1110-44dd-b1d5-1c0cba7d2c1e` is mounted at `/media/old_analysis`.
- `/dev/sda1` UUID `4f808d0c-a30e-4491-b11b-252bb46e2508` remains mounted at `/media/archive`.
- `/etc/fstab` was updated; backup made as `/etc/fstab.pansy-20260603T080816Z.bak`.

Operational notes:
- Deploy code to remote machines only through git. Do not copy source files with
  `scp`, `rsync`, heredocs, or manual paste as a deployment mechanism. Commit and
  push locally, then `git pull`/`git clone` on the remote and install packages
  from the git checkout with the standard Python environment.
- For scientific/intermediate data products, use HDF5 (`.h5`) rather than NPZ.
  Do not create NPZ files for project data unless explicitly requested.
- Plot, annotate, and discuss SNR in dB unless the user explicitly asks for
  linear SNR. Label axes and stored derived fields clearly when converting.
- Standard local Python environment for this repo is the base conda environment:
  use `conda run -n base python ...` for local project scripts unless explicitly
  told otherwise.
- Standard analysis-server Python environment on `revontuli` is `/usr/bin/python3`
  with the user site `/home/j/.local/lib/python3.10/site-packages` enabled. Install
  project-adjacent Python packages there with `/usr/bin/python3 -m pip install
  --user ...`.
- `iau_meteor_showers` is a local git package, not a PyPI package. Source lives
  locally at `/Users/jvi019/src/iau_meteor_showers`; on `revontuli`, clone/pull it
  as a git checkout under `/home/j/src/iau_meteor_showers`, then install with
  `/usr/bin/python3 -m pip install --user /home/j/src/iau_meteor_showers`. After
  installation, plain `/usr/bin/python3 -c "import iau_meteor_showers"` should
  work without `PYTHONPATH`.

System noise calibration notes:
- The current local system-noise calibration sidecar is
  `/Users/j/src/pansy_receiver/figs/cut_noise_measurements_2025-05-10_guard25_mean_per_module.h5`.
  It was generated from Revontuli cut metadata for 2025-05-10 using 25-sample
  guards and mean raw-voltage power, preserving per-receiver-module powers.
  The matching Revontuli path is
  `/home/j/src/pansy_receiver/figs/cut_noise_measurements_2025-05-10_guard25_mean_per_module.h5`.
- The current paper system-noise model fit uses PyGDSM/GSM2008 survey maps
  convolved with the one-way single-module receive gain pattern, not two-way
  gain and not the full transmit array. The adopted receive element pattern is
  a blend of 90% canonical crossed-Yagi module element pattern and 10%
  isotropic element response (`--element-pattern-blend 0.1`).
- Adopted module-0 calibration from the 2025-05-10 fit:
  `T_rec = 535.6 K`, bootstrap 16/50/84 percentiles about
  `497.9 / 539.3 / 574.8 K`, and fractional system-noise model RMS `5.7%`.
  Fit minimizes relative power residuals and enforces `T_rec >= 273 K`.
- Per-module raw-power equalization factors from that fit, normalized to module
  0 power: module 0 `1.000`, 1 `2.501`, 2 `1.380`, 3 `1.192`, 4 `1.681`,
  5 `2.732`, 6 `9.067`. Voltage factors are the square roots of these.
- Useful local scripts are `save_cut_noise_measurements.py`,
  `fit_cut_noise_pygdsm.py`, `paper_plot_noise.py`, and `pansy_gain.py`.
  Regenerate the per-module measurement sidecar with a command of the form:
  `python save_cut_noise_measurements.py --day 2025-05-10 --metadata /mnt/data/juha/pansy/metadata/cut --output figs/cut_noise_measurements_2025-05-10_guard25_mean_per_module.h5 --workers 10 --chunk-seconds 60 --guard-samples 25 --per-module`.
  Regenerate the fit/plot with a command of the form:
  `python fit_cut_noise_pygdsm.py --day 2025-05-10 --measurements figs/cut_noise_measurements_2025-05-10_guard25_mean_per_module.h5 --output figs/cut_noise_pygdsm_fit_2025-05-10_guard25_module0_blend0p1_relative_paper.png --element-pattern-blend 0.1 --relative-fit --module 0 --n-bootstrap 200`.
- Paper/memo outputs from the current fit are in `/Users/j/src/pansy_paper`:
  `system_noise_model.png` for the paper and
  `memos/figures/system_noise_fit_2025-05-10_module0_blend0p1_relative.png`
  for the memo. The detailed writeup is
  `/Users/j/src/pansy_paper/memos/cut_padding_noise_debug.tex`.
- All PANSY user services were active after deployment.
- Fresh tx metadata was initialized at `/media/analysis/metadata/tx/dmd_properties.h5`.
- Downstream services use `pansy_wait_for_paths.sh` so they wait for required metadata instead of crash-looping on a fresh analysis disk.
- Fresh-disk runtime fixes were added for `find_mode_starts.py`, `quick_search_meteor.py`, `mesomode_boundary.py`, `tx_xphase.py`, and `meso_xc.py`.
- Remote repo has many pre-existing untracked operational files; tracked branch was clean at `d05828b` after deploy.

Cleanup notes:
- `pansy-delete-old-rf.service` runs `receiver/scripts/pansy_delete_old_rf.sh`, deleting `rf*.h5` under `/media/archive/ch*` with default `PANSY_DELETE_RF_MTIME=+10` and `PANSY_DELETE_RF_INTERVAL_SECONDS=60`.
- On 2026-06-03, `pansy_delete_old_rf.log` showed many `find: ... No such file or directory` messages for old `rf*.h5` files. This was a normal live-tree readdir/delete race: `find` saw a file, then it vanished before `find -delete` reached it.
- Commit `0b86652` added GNU `find -ignore_readdir_race` support, but initially put the option before the path list, which caused `find: paths must precede expression`.
- Commit `3adf230` fixed the option order. The live command should look like `find /media/archive/ch000 ... /media/archive/ch007 -ignore_readdir_race -type f -name rf*.h5 -mtime +10 -delete`.
- After deploying `3adf230`, `pansy-delete-old-rf.service` restarted active. Old race and option-order errors remain in the append-only log, so distinguish them from fresh lines after the restart.

Mesomode/xc notes:
- `pansy_find_mode_starts` was producing tx metadata correctly, but `mesomode_boundary.py` originally only wrote closed mesomode blocks after seeing a gap.
- On the fresh analysis disk, tx metadata initially had long mesosphere-mode stretches and no mesomode samples were written, leaving downstream readers with only `dmd_properties.h5`.
- Commit `0cbe323` makes `mesomode_boundary.py` write safe open mesomode blocks through a 20-minute lag horizon and retry every 5 minutes.
- Commit `b55be43` lets `meso_xc.py` start from the first mesomode block when `xc2` metadata is empty.
- Commit `2c8bc83` makes `meso_xc.py` retry every 5 minutes instead of hourly.
- After deploying `0cbe323`, `/media/analysis/metadata/mesomode` became readable and had blocks through about `2026-05-26T09:03:05Z`.
- At that moment, the only >60 s mesomode blocks were before the raw `ch000` ringbuffer start (`2026-05-26T08:53:11Z`), so `xc2` had no files yet; this should advance once a long enough mesomode block falls inside the raw ringbuffer.

Quick-search notes:
- Before commit `e8732c4`, `quick_search_meteor.py` used `db_mf[1]` even when `mf` metadata was unreadable. On a fresh analysis disk this made the search start near Unix time `-60s` and grind through empty minute windows.
- Commit `c5d5f83` replaced repeated per-rank `no meso mode found` logs with a rank-0 summary.
- Commit `e8732c4` changed quick search to use `/media/analysis/metadata/mesomode` blocks as the work queue and to start from the mesomode/raw/tx overlap when `mf` metadata is empty.
- Commit `f22824b` suppresses the repeated Digital RF `read_vector_c81d` deprecation warning in quick-search logs.
- After deploying `f22824b`, `pansy_quick_search_meteor.log` reported work queues such as `2026-05-26T09:22:40Z` to `2026-05-26T09:41:46Z` with 29 mesomode blocks, and no longer emitted the old `no meso mode found` spam.

Cluster notes:
- `pansy_cluster_mf.log` previously failed at `dm_mf.get_bounds()` when `mf` metadata existed but was not yet readable.
- Commit `d05828b` makes `cluster_mf.py` wait cleanly for readable `mf`/`detections` metadata, caps stale `/tmp/meteor_mf_*.h5` progress to actual `mf` bounds, and retries every 5 minutes instead of hourly.
- After deploying `d05828b`, `cluster_mf.py` processed 137 ten-second windows from `2026-05-26T09:14:00Z` to `2026-05-26T09:36:50Z`.
- Verification after deployment: `mf` metadata was readable with 17,943 keys through `2026-05-26T09:41:46Z`; `detections` metadata was readable with 74 keys through `2026-05-26T09:36:44Z`.

Cut/process notes:
- Commit `702cafb` makes `cut_meteors.py` wait for readable detections/cut metadata, skip cuts outside the raw ringbuffer, suppress the Digital RF warning, and retry every 5 minutes; it also makes `process_cut_meteor.py` wait cleanly for readable cut/simple_fit metadata instead of aborting when cut metadata exists but has no readable samples.
- After deploying `702cafb`, `cut` metadata was readable from `2026-05-26T09:32:50.32Z` to `2026-05-26T09:55:55.64Z`, and `simple_fit` initially had 6 keys through `2026-05-26T09:36:40.75Z`.
- `pansy_process_cut_meteor.log` then exposed a separate MPI race: both ranks could hit Digital RF's first `simple_fit` field creation at the same time, producing `ValueError: Unable to synchronously create dataset (name already exists)` and causing `mpirun` to abort.
- Commit `b68762e` serializes `simple_fit` metadata writes with `/tmp/pansy_simple_fit_metadata.lock` and reopens the Digital RF writer if a rank sees the first-schema `name already exists` race.
- After deploying `b68762e` and restarting `pansy-process-cut-meteor.service`, `simple_fit` advanced to 84 keys and spans `2026-05-26T09:32:50.33Z` to `2026-05-26T09:55:55.73Z`; the service remained active under `mpirun -np 2`.
- Commit `d751ed0` reduces `process_cut_meteor.py` latency: work chunks are now 30 seconds by default instead of 120, idle polling is 15 seconds instead of 60, and the loop immediately repolls after actual `simple_fit` writes. These are configurable with `PANSY_PROCESS_CUT_BLOCK_SECONDS`, `PANSY_PROCESS_CUT_SLEEP_SECONDS`, and `PANSY_PROCESS_CUT_IDLE_LOG_SECONDS`.
- `d751ed0` also throttles idle waiting logs to once every 5 minutes by default and has `process_cut()` return whether it actually wrote metadata, so low-SNR skipped cuts do not cause a tight loop.
- After deploying `d751ed0`, `pansy-process-cut-meteor.service` was active at remote HEAD `d751ed0`; `cut` remained at 107 keys through `2026-05-26T09:55:55.64Z`, and `simple_fit` remained readable at 84 keys through `2026-05-26T09:55:55.73Z`.
- The process log still contains some `pansy_interferometry.py:100 RuntimeWarning: invalid value encountered in arccos` messages during fitting. These warnings did not abort the MPI job, but they may be worth cleaning up later by clipping the arccos argument into `[-1, 1]`.
- Commit `cd8de58` fixes a serious `cut_meteors.py` backlog bug: it computed minute-window count as `(bm[1]-start_idx)/60` instead of `(bm[1]-start_idx)/(60*1e6)`, causing absurd logs such as `23357736 minute windows ... 2070`. After deploy, it logged a sane `21 minute windows` range.
- Commit `78f9845` skips cut records with no qualifying beams and makes `process_cut_meteor.py` skip malformed/empty cut metadata instead of throwing `IndexError: tuple index out of range`.
- Health sweep after `78f9845`: all PANSY services were active; raw voltage was current within about 1 second; metadata products were still backlogged around May 26 and catching up after the disk incident.
- Commit `89b5c3d` fixes a `process_cut_meteor.py` duplicate-output stall where a cut metadata key could be slightly after the `simple_fit` output sample. The processor now starts one second after the latest fit bound and logs `simple_fit already exists; skipping` instead of traceback spam for duplicate writes.
- After deploying `89b5c3d`, `pansy-process-cut-meteor.service` restarted active and the fresh log showed `process_cut_meteor waiting: start 2026-05-26T11:17:56.18Z is not before cut end 2026-05-26T11:17:55.23Z`, confirming it advanced past the previous duplicate-write loop.

Useful checks:
- `systemctl --user --no-pager --type=service --all | grep pansy`
- `tail -f /home/radar/src/git/pansy_receiver/logs/pansy_uhd_rx.service.log`
- `tail -f /home/radar/src/git/pansy_receiver/logs/pansy_find_mode_starts.log`
- `df -h /media/analysis /media/old_analysis /media/archive`

Web monitor notes:
- Public monitor URL: `https://juha.no/pansy/`.
- Web root on `juha.no`: `/var/www/html/pansy`.
- Main page: `/var/www/html/pansy/index.html`.
- Backup made before the 2026-06-03 web edit: `/var/www/html/pansy/index.html.pre-pansy-service-20260603T0846Z.bak`.
- Backup made before the 2026-06-03 disk-warning web edit: `/var/www/html/pansy/index.html.pre-disk-warning-20260603T1005Z.bak`.
- The page now has a "Receiver Operations" card documenting the systemd service deployment and analysis disk swap.
- The page now has a prominent "Disk Space Warning" card near the top. It says disk pressure triggered the recent incident, tells operators to check `df -h /media/analysis /media/old_analysis /media/archive`, and lists the 2026-06-03 UTC snapshot: `/media/analysis` 1% used, `/media/old_analysis` 85% used, `/media/archive` 72% used.
- The page now has "Syowa Service Status" and "Syowa Disk Usage" cards populated from `/var/www/html/pansy/ops_status.json`.
- Commit `1a4d0be` makes `status_plot.py` write `ops_status.json` with 12 user systemd service states and disk usage for `/media/analysis`, `/media/old_analysis`, and `/media/archive`; it rsyncs the JSON with the existing monitor artifacts using conservative `PANSY_STATUS_RSYNC_BWLIMIT=1`.
- Web backup made before adding process/disk cards: `/var/www/html/pansy/index.html.pre-ops-status-20260603T1103Z.bak`.
- Initial published `ops_status.json` on 2026-06-03T11:02Z showed 12 services, `/media/analysis` 0.0% used, `/media/old_analysis` 80.5% used, and `/media/archive` 68.2% used.
- The page now computes the age of `status.json` from `generated_utc` and warns when backup/mirror status is stale.
- As of the web edit, `status.json` was generated on 2026-05-28 and plot PNGs were last modified on 2026-06-01, so the browser page warns that backup status may be stale.
- `status_plot.py` on the Syowa receiver rsyncs plot/image artifacts to `/var/www/html/pansy`, but does not generate `status.json`.
- Commit `eff8adb` makes `status_plot.py` tolerate missing optional metadata such as `/media/analysis/metadata/mf_isr` and `/media/analysis/metadata/xc2`. Missing metadata is logged as unavailable and rendered as missing/NaN instead of crashing the status loop.
- Commit `fc9e416` adds a normal `if __name__ == "__main__"` guard around the status plot loop, so `import status_plot; status_plot.plot_status()` is safe for one-shot tests.
- Commit `fcdb0f7` briefly changed the status plot rsync bandwidth limit to default `1000`, but this was too aggressive for the limited Syowa internet link.
- Commit `e277155` restores the conservative default `PANSY_STATUS_RSYNC_BWLIMIT=1`. Keep this low unless the user explicitly asks to spend more bandwidth.
- After deploying the status plot fixes, a one-shot `status_plot.plot_status()` completed, rsynced `fit_data.h5`, `latest_hist.png`, `latest_meteor.png`, `latest_radiants.png`, `processing.png`, `raw.png`, and `status.png`; `pansy-status-plot.service` was restarted active after `e277155`.
- Commit `0ffcc0b` disables ISR mode metadata recording for now. `find_mode_starts.py` no longer searches for or writes ISR mode ID `2` unless `PANSY_ENABLE_ISR_MODE=1` is explicitly set.
- `0ffcc0b` also removes `mf_isr` from `status_plot.py`, so the status plot no longer reads or logs `/media/analysis/metadata/mf_isr`.
- After deploying `0ffcc0b`, `pansy-find-mode-starts.service` and `pansy-status-plot.service` were restarted active. Fresh `find_mode_starts` logs after restart showed only `b13` and `mesosphere mode` lines, not `isr mode` lines. Existing historical tx metadata was not deleted.
