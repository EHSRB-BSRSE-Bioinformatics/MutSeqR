# MutSeqR 0.99.0 (2025-06-19)

Initial public version.

### Major changes

- Added `filter_mut()` to workflow: germline identification via `vaf_cutoff`, region filtering, and depth correction now occur here instead of the import functions.
- `calculate_mut_freq()` is renamed to `calculate_mf()`.
- `calculate_mf()` no longer requires depth; users may:
  1. calculate depth from mutation data,
  2. supply a separate depth table, or
  3. omit depth entirely (only mutation counts returned).
- `correct_depth` option moved to `calculate_mf()`.
- `plot_spectra()`, `plot_trinucleotide()`, and `spectra_comparison()` now use `mf_data` instead of raw mutations.
- Output options added: VCF, FASTA, SigProfiler-compatible format, Excel workbook.
- Example dataset (~44MB) added.

### New features

- `render_report()` added for standardized summary reporting.

### Other

- Removed `custom_regions` parameter; replaced by generalized `regions` argument.
- Public release 🎉
- See the [vignette](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/articles/MutSeqR_introduction.html) for details.