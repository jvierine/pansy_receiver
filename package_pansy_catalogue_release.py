#!/usr/bin/env python3
"""Create and verify the three PANSY catalogue Zenodo ZIP archives."""

import argparse
import hashlib
import shutil
import zipfile
from pathlib import Path


ROOT_NAME = "pansy_catalogue_v1"
COMMON_FILES = ("README.md", "release_summary.json", "SHA256SUMS", "zenodo_metadata.json")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def release_checksums(release):
    files = sorted(
        path
        for path in release.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS"
    )
    text = "\n".join(f"{sha256(path)}  {path.relative_to(release)}" for path in files)
    (release / "SHA256SUMS").write_text(text + "\n")


def write_archive(output, release, members):
    with zipfile.ZipFile(
        output,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=9,
        allowZip64=True,
    ) as archive:
        for path in members:
            archive.write(path, Path(ROOT_NAME) / path.relative_to(release))
    with zipfile.ZipFile(output, "r") as archive:
        bad = archive.testzip()
        if bad is not None:
            raise RuntimeError(f"CRC failure in {output}: {bad}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path(__file__).parent / "zenodo" / "catalogue_v1" / "zenodo_metadata.json",
    )
    args = parser.parse_args()
    release = args.release_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    shutil.copy2(args.metadata, release / "zenodo_metadata.json")
    release_checksums(release)
    common = [release / name for name in COMMON_FILES]
    products = {
        "level1": common
        + sorted(release.glob("pansy_level1_example_*.h5")),
        "level2": common
        + [release / "example_level2_radiant.py", release / "verify_level2_level3.py"]
        + sorted((release / "level2").glob("*.h5")),
        "level3": common
        + [release / "verify_level2_level3.py"]
        + sorted((release / "level3").glob("*.h5")),
    }

    archives = []
    for level, members in products.items():
        missing = [path for path in members if not path.is_file()]
        if missing:
            raise FileNotFoundError(", ".join(str(path) for path in missing))
        archive = output / f"pansy_catalogue_v1_{level}.zip"
        write_archive(archive, release, members)
        archives.append(archive)
        print(f"{archive.name}: {archive.stat().st_size} bytes")

    checksum_path = output / "pansy_catalogue_v1_ZIP_SHA256SUMS"
    checksum_path.write_text(
        "\n".join(f"{sha256(path)}  {path.name}" for path in archives) + "\n"
    )
    print(checksum_path)


if __name__ == "__main__":
    main()
