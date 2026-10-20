# Memos 

- [2024-09-27 Memo 001: Testing of the receiver at Shigaraki](001_shigaraki_testing.pdf)
- [2024-09-27 Memo 002: Observation log of PANSY tests](notes20240927.txt)
- [2024-09-27 Memo 003: Detecting PANSY modes](003_tx_detection0.pdf)
- [2024-10-02 Memo 004: Simple analysis of ST and M modes](004-range-doppler-test0.pdf)
- [2024-10-02 Memo 005: UPS failure analysis](005_ups_issues)
- [Memo 013: Sparse complex spectra downloads - file contents and data layout](../pmse/rvd/pansy-sparse-spectra-guide.pdf) ([LaTeX source](013_sparse_spectra_download.tex))

Memo 013 uses the unmodified `jvierine/uit-memo-style` package and UiT seal at
revision `06dd4b15d32ddeaf06f6c9df583079b84f494c87`. The private style is an
external build dependency; its sources are not vendored into this repository.
With that checkout next to this repository (or `UIT_MEMO_STYLE_DIR` pointing at
it), run `bash scripts/build_pmse_download_guide.sh`. This uses pdfLaTeX without
shell escape, checks the pinned style, and builds both the local PDF under
`output/pdf/` and the public static asset beside the download Help link.
The memo includes executable retained-spectrum inverse-FFT and noise-model
examples; `tests/test_pmse_download_memo.py` runs the printed Python examples,
including the FFT round trip and preservation of observed bins during synthetic
noise filling. `tests/test_pmse_reconstruction.py` tests the reversibility
boundaries explicitly. The guide now covers rectangular constant-gain reversal,
spectral resolution, and historical Hann compatibility. The noise simulation
does not change measured coefficients; the separately requested oldest-first
backfill reprocesses available raw recordings and replaces historical Hann
records with recoverable, atomic publication.
