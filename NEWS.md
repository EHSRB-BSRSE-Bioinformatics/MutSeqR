# MutSeqR 0.99.4 (2025-12-09)

Preparing for Bioconductor release. We address comments from Bioconductor reviewers. Major changes include parameter validation and vectoriation of functions. filter_mut() "snv_in_germ_mnv" parameter was fixed such that it now only filters out overlapping snvs if their variation matches that of the germline mnv (previously it was blanket removing all snvs that overlap with germline mnvs). plot_lollipop() can now colour the plot by different subtype resolutions (previously just base_6). Fixed plot_spectra() axes labels after clustering (previously not showing labels).

# MutSeqR 0.99.3 (2025-10-10)

Fix dependencies for Bioconductor release.

# MutSeqR 0.99.2 (2025-09-15)

Preparing for Bioconductor release. This change removes bmd_toxicR from the package. ToxicR dependency is not supported by Bioconductor. See ToxicR_archive branch for bmd_toxicr function. Modifies signature_fitting() to use SigProfilerMatrixGenerator python dependency rather than R dependency.

# MutSeqR 0.99.1 (2025-08-13)

Preparing for Bioconductor release. This change adds MutSeqRData to suggests, and alters how the examples are run (now depends on the ExperimentHub accessions).

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
