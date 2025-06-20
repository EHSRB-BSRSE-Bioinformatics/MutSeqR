

# MutSeqR <a href="https://https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/"><img src="man/figures/temp-hex.png" align="right" height="138" /></a>

<!-- badges: start -->
  [![R-CMD-check](https://github.com/EHSRB-BSRSE-Bioinformatics/duplex-sequencing/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EHSRB-BSRSE-Bioinformatics/duplex-sequencing/actions/workflows/R-CMD-check.yaml)
  [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
  <!-- badges: end -->

## Overview
MutSeqR is an open-source R package to analyze error-corrected Next-Generation
Sequencing (ECS) data, empowering users with flexibility during exploratory analyses while ensuring compatibility across technologies.

<img src="man/figures/MutSeqR overview.png" align=center alt="A Flowchart showing MutSeqR's function utility and workflow: Data Import, Data Processing, Statistical Analyses, Visualization, Output. Includes a visual of a woman working at a computer.">

<p style="font-size: 0.5em;">
Figure transcript: 1. Data Import: Imports mutation data into the R environment. Binds data from multiple libraries into a single object. Joins sample and target region metadata to the mutation data. Retrieves trinucleotide context. Functions: import_mut_data() & import_vcf_data(). Inputs: Raw Mutation data (VCF or tabular), sample metadata, sequencing panel ranges. 2. Data Processing: Calculates mutation frequencies for gruops of interest. Calculates frequencies andproportions of mutation subtypes. Opinoal Variant filtering: eliminates putative germline variants, removes variants outside of specified regions, quality assurance filtering. Functions: calculate_mf() & filter_mut(). Inputs: mutation data. Outputs: mf data. 3. Statistical Analyses: Generalized linear modeling. Benchmark Dose Modeling. COSMIC signature analysis. Spectra comparison between groups. Unsupervised clustering based on mutation spectra. Functions: model_mf(), bmd_proast(), bmd_toxicr(), signature_fitting(), spectra_comparison(), cluster_spectra(). Inputs: mf data. 4. Visualization: Create figures to display mutation frequencies and the proportions of mutation subtypes. VIsualise statistical results. Visualise mutation distribution across genomic loci. View clonal expansion of mutations. Functions: plot_mf(), plot_mean_mf(), plot_model_mf(), plot_ci(), plot_spectra(), plot_trinucleotide(), plot_trinucleotide_heatmap(), plot_bubbles(), plot_radar(). 5. Output: Summary report RMarkdown file will faciliatte the generation of results. Output mutation data as VCF. Output sequences in FASTA format. Output spectra data in SigProfiler format. Export results to Excel workbook. Functions: render_report(), write_vcf_from_mut(), write_reference_fasta(), get_sigprofiler_output(), write_excel().
</p>

## Vignette

See the [vignette](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/articles/MutSeqR_introduction.html#introduction) for details on function utility.

## Recent Changes:
Changes: 2025-06-19
- We're Public!!
- See the Vignette for function details.
- Moved the correct_depth feature to calculate_mf

Changes: 2025-06-05
- Summary Report: Users may now run a standardized analysis workflow using render_report().

Changes: 2025-03-27
- Removed custom_regions parameter. Utility is now incorporated by regions parameter.

Major changes on 2025-03-24
- *filter_mut()* function added to workflow. This function filters the mutation_data: germline identification via vaf_cutoff, depth correction, and filtering variants based on regions have all been moved from the import functions to filter_mut(). calculate_mf(), plot_bubbles(), and signature_fitting() filter out variants using the filter_mut column instead of the is_germline column.
- calculate_mut_freq() is renamed to calculate_mf()
- total_depth and depth are no longer required for import
- calculate_mf() no longer requires depth. Users may choose to 1) calculate depth from mutation data, 2) supply precalculated depth values in a separate table, 3) No depth, mf is not calculated, only mutation sums.
- plot_spectra, plot_trinucleotide, and spectra_comparison are now supplied with mf_data instead of the mutation data.
- Example data has been added: Currently 44Mb

## Installation

Install the package from github:

```{r}
# install.packages("devtools")
devtools::install_github("EHSRB-BSRSE-Bioinformatics/MutSeqR", auth_token = "your personal_access_token from GitHub")
```

Load the package
```{r}
library(MutSeqR)
```

## Getting Help

If you encounter a clear bug, please file an issue with a minimal reproducible example on [Github](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/issues).

## Citation

To cite this package in publications use:

Meier M, Dodge A, Williams A, Esina E (2025). MutSeqR: Analysis of Error-Corrected Sequencing Data for Mutation Detection. R package version 0.99.0, https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/.