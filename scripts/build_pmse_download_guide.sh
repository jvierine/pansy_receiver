#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
memo_style_dir="${UIT_MEMO_STYLE_DIR:-$PWD/../uit-memo-style}"
memo_style_revision=06dd4b15d32ddeaf06f6c9df583079b84f494c87
if [[ ! -f "$memo_style_dir/uitmemo.sty" || ! -f "$memo_style_dir/assets/uit-circular-logo.pdf" ]]; then
  echo "Set UIT_MEMO_STYLE_DIR to the jvierine/uit-memo-style checkout." >&2
  exit 1
fi
memo_style_dir="$(cd "$memo_style_dir" && pwd)"
if [[ "$(git -C "$memo_style_dir" rev-parse HEAD)" != "$memo_style_revision" ]]; then
  echo "Expected UiT memo style revision $memo_style_revision." >&2
  exit 1
fi
git -C "$memo_style_dir" diff --quiet HEAD -- uitmemo.sty assets/uit-circular-logo.pdf
mkdir -p output/pdf tmp/pdfs/uit-download-guide
export TEXINPUTS="$memo_style_dir:${TEXINPUTS:-}:"
for pass in 1 2 3; do
  if ! pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error \
    -output-directory=tmp/pdfs/uit-download-guide -jobname=pansy-sparse-spectra-guide \
    "\def\uitstyledir{$memo_style_dir}\input{memos/013_sparse_spectra_download.tex}" \
    > tmp/pdfs/uit-download-guide/build.log; then
    tail -60 tmp/pdfs/uit-download-guide/build.log >&2
    exit 1
  fi
done
cp tmp/pdfs/uit-download-guide/pansy-sparse-spectra-guide.pdf output/pdf/pansy-sparse-spectra-guide.pdf
cp output/pdf/pansy-sparse-spectra-guide.pdf pmse/rvd/pansy-sparse-spectra-guide.pdf
printf '%s\n' "$PWD/output/pdf/pansy-sparse-spectra-guide.pdf"
