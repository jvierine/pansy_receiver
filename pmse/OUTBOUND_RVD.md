# PANSY RVD data flow

Only the work laptop may originate a reverse SSH tunnel. The RVD data path
uses ordinary outbound SSH sessions and no port forwarding:

1. Syowa creates stable, lossless native HDF5 snapshots and publishes them to
   `juha.no:/mnt/shovel/rvd/pmse_sparse16` at the existing 10 KiB/s link cap.
2. Revontuli initiates checksum-verified pulls of that archive into its retained
   local copy. No files are deleted from Revontuli by this operation.
3. The public API, data packets, and Rustplotlib axes are served on juha.no
   against its local archive. Apache proxies only to a process on juha.no's
   own loopback interface. Public requests never connect into Revontuli.

Install source code through Git. Service templates in `pmse/systemd` have
explicit host roles: `pansy-rvd-source-publish` on Syowa,
`pansy-pmse-sparse-backup` on Revontuli, and `pansy-rvd-browser-api` plus
`pansy-rvd-web` on juha.no. Do not install the public API on Revontuli or
recreate the removed browser tunnel service.

Publication receipts record source versions only after completed transfers.
Partial transfers remain pending. Revontuli receipts contain SHA256 hashes of
its retained files, not hashes asserted from an unverified remote copy.
Same-size reprocessed snapshots are checked by content on the pull leg.
The pull operation and source publisher do not expire or delete archive files.
