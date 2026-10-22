# PANSY RVD data flow

Only the work laptop may originate a reverse SSH tunnel. The RVD data path
uses ordinary outbound SSH sessions and no port forwarding:

1. Syowa creates stable, lossless native HDF5 snapshots and publishes them to
   `juha.no:/mnt/shovel/rvd/pmse_sparse16` at the existing 10 KiB/s link cap.
2. Revontuli initiates incremental size/mtime pulls every 120 seconds into its
   retained local copy. No files are deleted from Revontuli by this operation.
   Full content verification is an explicit one-off `pull --checksum-audit`
   without `--loop`, so ordinary polling never rereads every cloud HDF5 file.
3. The public API, data packets, and Rustplotlib axes are served on juha.no
   against its local archive. Apache proxies only to a process on juha.no's
   own loopback interface. Public requests never connect into Revontuli.

Install source code through Git. Service templates in `pmse/systemd` have
explicit host roles: `pansy-rvd-source-publish` on Syowa,
`pansy-pmse-sparse-backup` on Revontuli, and `pansy-rvd-browser-api` plus
`pansy-rvd-web` on juha.no. Do not install the public API on Revontuli or
recreate the removed browser tunnel service.

Sky-noise analysis on Revontuli publishes numeric arrays and fitted metadata
in `data/sky-noise/status.json`; it no longer renders diagnostic SVGs. Archive
arrays already arrive in the normal GUI inventory. Rustplotlib's shared
interactive series API draws both sets of diagnostics, including axes and
legends, entirely in the browser. No chart images or public requests to
Revontuli are involved, and juha.no does not perform scientific reductions.

Publication receipts record source versions only after completed transfers.
Partial transfers remain pending. Revontuli receipts contain SHA256 hashes of
its retained files, not hashes asserted from an unverified remote copy.
Same-size reprocessed snapshots are detected through the fresh archive mtime
assigned by atomic publication; receipts hash changed files on Revontuli.
The pull operation and source publisher do not expire or delete archive files.

The seven legacy metadata backup channels share `RSYNC_TOTAL_BWLIMIT_KB=10`,
partitioned into 1 KiB/s per channel by default, including cutouts. Combined
with the separate 10 KiB/s RVD publisher and 1 KiB/s status publisher, this
budgets 18 KiB/s of file payload from Syowa. SSH/TCP overhead, file listings,
interactive sessions, and unrelated station traffic are additional; this is
not a router-level cap on the entire station internet connection.
