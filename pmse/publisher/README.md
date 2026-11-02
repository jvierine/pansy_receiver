# Incremental RVD publisher

`pansy-rvd-publisher` reads the small, source-produced `.h5.json.gz` indexes
under the local web directory. Hour-directory timestamps identify changes.
Only affected UTC days are rebuilt. Unchanged cycles preserve output mtimes.
Cached native checksums survive service restarts. Native files are opened only
for pending deliveries, with a configurable byte budget per 60-second cycle;
the archive is never recursively scanned or opened for inventory calculations.
Existing inventory and analysis-produced archive plots are retained.

```sh
cargo test --manifest-path pmse/publisher/Cargo.toml
cargo build --release --locked --manifest-path pmse/publisher/Cargo.toml
pmse/publisher/target/release/pansy-rvd-publisher \
  --native-root /mnt/shovel/rvd/pmse_sparse16 \
  --output /var/www/rvd-live/data \
  --python /home/j/venv/pansy-rvd/bin/python \
  --sidecar-script /path/to/checkout/pmse/write_gui_sidecar.py
```

Use `--once` for an isolated benchmark. `--state` and `--source` override the
default local cache and source-index locations. `--verify-bytes` defaults to
4194304 per cycle. Do not run two publishers against the same output directory.
The Python helper writes the downloadable HDF5 from the already prepared JSON
only when a day changes; it never reads native spectra or walks the archive.
Monitor SNR and noise projections are generated directly in Rust.

For the existing web-host service, substitute this checkout's absolute path
for `__PANSY_RECEIVER_REPO__` in `pmse/systemd/pansy-rvd-web-rust.conf` and
install the result as
`~/.config/systemd/user/pansy-rvd-web.service.d/30-rust-publisher.conf`.
Reload user units and restart `pansy-rvd-web.service`. This replaces the Python
publisher in that same service, avoiding competing writers. It retains the
2% CPU ceiling. Build with a separate transient unit so build limits do not
change the publisher's runtime limits.

On juha.no, `bash pmse/install_rust_publisher.sh` initializes the persistent
index once in a transient unit capped at 20% CPU, then starts the regular
publisher with its 2% ceiling. Initialization failure restarts the previous
publisher. The script uses passwordless sudo for that temporary build-free
initialization unit and records its CPU time and peak resident memory.
