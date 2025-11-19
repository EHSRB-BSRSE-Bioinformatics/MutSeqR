# Changelog

## MutSeqR 0.99.2 (2025-09-15)

Preparing for Bioconductor release. This change removes bmd_toxicR from
the package. ToxicR dependency is not supported by Bioconductor. See
ToxicR_archive branch for bmd_toxicr function. Modifies
signature_fitting() to use SigProfilerMatrixGenerator python dependency
rather than R dependency.

## MutSeqR 0.99.1 (2025-08-13)

Preparing for Bioconductor release. This change adds MutSeqRData to
suggests, and alters how the examples are run (now depends on the
ExperimentHub accessions).

## MutSeqR 0.99.0 (2025-06-19)

Initial public version.

#### Major changes

- Added
  [`filter_mut()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/filter_mut.md)
  to workflow: germline identification via `vaf_cutoff`, region
  filtering, and depth correction now occur here instead of the import
  functions.
- `calculate_mut_freq()` is renamed to
  [`calculate_mf()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/calculate_mf.md).
- [`calculate_mf()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/calculate_mf.md)
  no longer requires depth; users may:
  1.  calculate depth from mutation data,
  2.  supply a separate depth table, or
  3.  omit depth entirely (only mutation counts returned).
- `correct_depth` option moved to
  [`calculate_mf()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/calculate_mf.md).
- [`plot_spectra()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/plot_spectra.md),
  [`plot_trinucleotide()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/plot_trinucleotide.md),
  and
  [`spectra_comparison()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/spectra_comparison.md)
  now use `mf_data` instead of raw mutations.
- Output options added: VCF, FASTA, SigProfiler-compatible format, Excel
  workbook.
- Example dataset (~44MB) added.

#### New features

- [`render_report()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/render_report.md)
  added for standardized summary reporting.

#### Other

- Removed `custom_regions` parameter; replaced by generalized `regions`
  argument.
- Public release 🎉
- See the
  [vignette](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/articles/MutSeqR_introduction.html)
  for details.
