#!/usr/bin/env python3
"""End-to-end PANSY meteor orbit-determination example.

This example is intentionally verbose and heavily commented. It is meant to be
read by someone who is new to the PANSY meteor processing chain and wants a
single script that starts from the public Zenodo example dataset and ends with
metadata products that can be consumed by later analysis scripts.

Dataset citation:

Vierinen, J. (2026). PANSY Meteor Head Echo Example Dataset (1.0.0)
[Data set]. Zenodo. https://doi.org/10.5281/zenodo.20694794

The workflow is:

1. Download the Zenodo tarball and SHA-256 checksum into ./data/zenodo/.
2. Verify the checksum.
3. Extract the archive into ./data/, producing ./data/metadata/simple_meteor_fit
   and the other supporting metadata channels.
4. Select one "nice" meteor from simple_meteor_fit.
5. Fit a robust MSIS-drag ballistic trajectory to the per-pulse ENU positions.
6. Reverse-propagate the fitted trajectory to an above-atmosphere state.
7. Convert that state to the GCRS reference frame.
8. Estimate heliocentric orbital parameters and uncertainties.
9. Write both derived products as Digital Metadata channels:
   ./data/metadata/ballistic_fit
   ./data/metadata/orbital_parameters

The derived products are metadata, not loose JSON sidecars. This is deliberate:
downstream processing can use the same DigitalMetadataReader pattern for raw
fits, ballistic fits, and orbital parameters.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import urllib.request
from pathlib import Path

import pansy_ballistic as pbal
import pansy_orbit as porb


# Zenodo record and file names. The archive uses zstd compression because it is
# smaller and faster than gzip for this 2.5 GB metadata package.
DOI = "10.5281/zenodo.20694794"
RECORD_URL = "https://zenodo.org/records/20694794"
ARCHIVE_NAME = (
    "pansy_meteor_head_echo_example_dataset_2025_eta_aquariids_"
    "metadata_20250506_no_tx.tar.zst"
)
ARCHIVE_URL = (
    "https://zenodo.org/api/records/20694794/files/"
    f"{ARCHIVE_NAME}/content"
)
SHA256_URL = (
    "https://zenodo.org/api/records/20694794/files/"
    f"{ARCHIVE_NAME}.sha256/content"
)


# The Zenodo example is exactly this UTC day, chosen around the 2025
# eta-Aquariids peak. Keeping the interval explicit avoids local-time ambiguity.
START = "2025-05-06T00:00:00+00:00"
END = "2025-05-07T00:00:00+00:00"


def progress_download(url: str, path: Path) -> None:
    """Download a URL to disk with coarse progress reports.

    The standard library is used here so the download step has no dependency on
    requests. This matters for a fresh user who has not yet installed the full
    science stack.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")

    with urllib.request.urlopen(url) as response, tmp.open("wb") as output:
        total_header = response.headers.get("Content-Length")
        total = int(total_header) if total_header else None
        received = 0
        next_report = 0

        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            output.write(chunk)
            received += len(chunk)

            # Print roughly every 250 MB. This keeps long downloads reassuring
            # without flooding the terminal.
            if received >= next_report:
                if total:
                    pct = 100.0 * received / total
                    print(
                        f"downloaded {received / 1e9:5.2f}/{total / 1e9:5.2f} GB "
                        f"({pct:4.1f}%)",
                        flush=True,
                    )
                else:
                    print(f"downloaded {received / 1e9:5.2f} GB", flush=True)
                next_report = received + 250 * 1024 * 1024

    tmp.replace(path)


def download_dataset(data_dir: Path, force: bool = False) -> tuple[Path, Path]:
    """Download the Zenodo archive and checksum if needed."""
    download_dir = data_dir / "zenodo"
    archive_path = download_dir / ARCHIVE_NAME
    sha_path = download_dir / f"{ARCHIVE_NAME}.sha256"

    if force or not sha_path.exists():
        print(f"Downloading checksum from {SHA256_URL}")
        progress_download(SHA256_URL, sha_path)
    else:
        print(f"Using existing checksum: {sha_path}")

    if force or not archive_path.exists():
        print(f"Downloading dataset from {ARCHIVE_URL}")
        progress_download(ARCHIVE_URL, archive_path)
    else:
        print(f"Using existing archive: {archive_path}")

    return archive_path, sha_path


def expected_sha256(sha_path: Path) -> str:
    """Read the first token from the sha256sum sidecar file."""
    first_field = sha_path.read_text().split()[0]
    if len(first_field) != 64:
        raise ValueError(f"Could not parse SHA-256 checksum from {sha_path}")
    return first_field.lower()


def verify_sha256(archive_path: Path, sha_path: Path) -> None:
    """Verify the downloaded archive before extraction."""
    expected = expected_sha256(sha_path)
    digest = hashlib.sha256()
    with archive_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    actual = digest.hexdigest()
    if actual != expected:
        raise RuntimeError(
            f"Checksum mismatch for {archive_path}\n"
            f"expected {expected}\n"
            f"actual   {actual}"
        )
    print(f"SHA-256 OK: {archive_path.name}")


def supports_tar_zstd() -> bool:
    """Return True if the system tar can directly read .tar.zst files."""
    result = subprocess.run(
        ["tar", "--zstd", "--help"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


def extract_archive(archive_path: Path, data_dir: Path, force: bool = False) -> None:
    """Extract the Zenodo archive into ./data.

    Extraction creates ./data/metadata/... while leaving the small calibration
    files that already live in ./data untouched.
    """
    metadata_dir = data_dir / "metadata" / "simple_meteor_fit"
    if metadata_dir.exists() and not force:
        print(f"Using existing extracted metadata: {data_dir / 'metadata'}")
        return

    data_dir.mkdir(parents=True, exist_ok=True)
    print(f"Extracting {archive_path.name} into {data_dir}")
    if shutil.which("tar") and supports_tar_zstd():
        subprocess.run(["tar", "--zstd", "-xf", str(archive_path), "-C", str(data_dir)], check=True)
    elif shutil.which("zstd") and shutil.which("tar"):
        # Portable fallback for systems where tar lacks --zstd support.
        with subprocess.Popen(["zstd", "-dc", str(archive_path)], stdout=subprocess.PIPE) as zstd:
            subprocess.run(["tar", "-xf", "-", "-C", str(data_dir)], stdin=zstd.stdout, check=True)
            if zstd.stdout:
                zstd.stdout.close()
            if zstd.wait() != 0:
                raise RuntimeError("zstd decompression failed")
    else:
        raise RuntimeError("Need either `tar --zstd` support or both `zstd` and `tar`.")


def extracted_dataset_exists(data_dir: Path) -> bool:
    """Return True if the extracted Zenodo metadata needed by the example exists."""
    return (data_dir / "metadata" / "simple_meteor_fit" / "dmd_properties.h5").exists()


def run_one_event(args: argparse.Namespace) -> None:
    """Fit one event and write both metadata products.

    Reanalysis behavior:
    Digital Metadata writers do not overwrite an existing sample key. When
    --reanalyze is used, this script deletes the existing top-level HDF5 group
    for the selected event in both derived metadata channels before writing the
    new products.
    """
    simple_dir = args.data_dir / "metadata" / "simple_meteor_fit"
    ballistic_dir = args.data_dir / "metadata" / "ballistic_fit"
    orbit_dir = args.data_dir / "metadata" / "orbital_parameters"

    # Select a reproducible event from simple_meteor_fit. The selection function
    # prefers events with many points, good SNR, and plausible geocentric speed.
    sample_idx, simple_record = pbal.select_nice_event(
        simple_dir,
        pbal.unix_us(args.start),
        pbal.unix_us(args.end),
    )
    print(f"Selected event sample index: {sample_idx}")

    # Reanalysis means "replace the derived products for this event." We delete
    # the matching metadata key before fitting so a failed later write does not
    # leave stale duplicate state.
    if args.reanalyze:
        deleted_ballistic = pbal.delete_metadata_key(
            ballistic_dir,
            sample_idx,
            pbal.BALLISTIC_WRITER_ARGS["file_name"],
        )
        deleted_orbit = pbal.delete_metadata_key(
            orbit_dir,
            sample_idx,
            porb.ORBIT_WRITER_ARGS["file_name"],
        )
        print(f"Reanalysis delete: ballistic_fit={deleted_ballistic}, orbital_parameters={deleted_orbit}")

    # The ballistic fit is robust in two passes:
    #   1. soft_l1 least-squares to reduce the leverage of wild points;
    #   2. per-pulse sigma clipping followed by a linear least-squares refit.
    # The fit returns both local ENU products and GCRS state vectors.
    ballistic = pbal.fit_ballistic_event(simple_record, sample_idx)
    pbal.write_fit(pbal.metadata_writer(ballistic_dir), ballistic)
    ballistic_plot = pbal.plot_fit(ballistic, args.plot_dir / "ballistic_fit")
    print(f"Wrote ballistic_fit metadata: {ballistic_dir}")
    print(f"Ballistic diagnostic plot: {ballistic_plot}")
    print(
        "Above-atmosphere extrapolation speed change: "
        f"{ballistic['speed_increase_to_above_atmosphere_mps'] / 1e3:.3f} km/s"
    )

    # Orbit determination consumes the above-atmosphere GCRS state. The
    # uncertainty is estimated by drawing ballistic parameter samples from the
    # fitted covariance, reverse-propagating each to above atmosphere, converting
    # each sample to GCRS, and recomputing orbital elements.
    orbit = porb.orbit_from_ballistic_fit(ballistic, n_samples=args.n_samples, seed=args.seed)
    porb.write_orbit(porb.metadata_writer(orbit_dir), orbit)
    orbit_plot = porb.plot_orbit(orbit, args.plot_dir / "orbital_parameters")
    print(f"Wrote orbital_parameters metadata: {orbit_dir}")
    print(f"Orbit diagnostic plot: {orbit_plot}")
    print("Kepler elements [a_au, e, i_deg, raan_deg, argp_deg, nu_deg, q_au]:")
    print(orbit["kepler"])
    print("1-sigma element uncertainties:")
    print(orbit["kepler_std"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one PANSY meteor orbit-determination example.")
    parser.add_argument("--data-dir", type=Path, default=Path("data"), help="Download and extraction directory.")
    parser.add_argument("--plot-dir", type=Path, default=Path("plots/orbit_determination_example"), help="Diagnostic plot directory.")
    parser.add_argument("--start", default=START, help="UTC start time for event selection.")
    parser.add_argument("--end", default=END, help="UTC end time for event selection.")
    parser.add_argument("--force-download", action="store_true", help="Download Zenodo files even if local copies exist.")
    parser.add_argument("--force-extract", action="store_true", help="Extract even if ./data/metadata already exists.")
    parser.add_argument("--skip-download", action="store_true", help="Require an already extracted dataset; do not contact Zenodo.")
    parser.add_argument("--download-only", action="store_true", help="Download and verify, but do not fit anything.")
    parser.add_argument("--reanalyze", action="store_true", help="Delete existing ballistic_fit and orbital_parameters keys before writing.")
    parser.add_argument("--n-samples", type=int, default=200, help="Monte Carlo samples for orbital-parameter uncertainty.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed for orbital uncertainty samples.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.skip_download:
        if not extracted_dataset_exists(args.data_dir):
            raise FileNotFoundError(
                f"Could not find extracted metadata in {args.data_dir / 'metadata' / 'simple_meteor_fit'}"
            )
    elif (
        extracted_dataset_exists(args.data_dir)
        and not args.force_download
        and not args.force_extract
        and not args.download_only
    ):
        print(f"Using existing extracted metadata: {args.data_dir / 'metadata'}")
    else:
        archive_path, sha_path = download_dataset(args.data_dir, force=args.force_download)
        verify_sha256(archive_path, sha_path)
        if args.download_only:
            return 0
        extract_archive(archive_path, args.data_dir, force=args.force_extract)

    run_one_event(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
