<!-- badges: start -->
  [![R-CMD-check](https://github.com/EHSRB-BSRSE-Bioinformatics/duplex-sequencing/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EHSRB-BSRSE-Bioinformatics/duplex-sequencing/actions/workflows/R-CMD-check.yaml)
  [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
  <!-- badges: end -->
  
# Change Report:
Changes: 2025-06-05
- Summary Report: Users may now run a standardized analysis workflow using render_report().
- Users may consult the MutSeqR vignette:

```{r}
  rmarkdown::render(
    input = system.file("vignette", "vignette.rmd", package="MutSeqR"),
    output_dir = NULL,
    output_file = "./Vignette_Rmd",
    envir = new.env()
  )
```

Change: 2025-06-19
- Moved the correct_depth feature to calculate_mf

Changes: 2025-03-27
- Removed custom_regions parameter. Utility is now incorporated by regions parameter.

Major changes on 2025-03-24
- *filter_mut()* function added to workflow. This function filters the mutation_data: germline identification via vaf_cutoff, depth correction, and filtering variants based on regions have all been moved from the import functions to filter_mut(). calculate_mf(), plot_bubbles(), and signature_fitting() filter out variants using the filter_mut column instead of the is_germline column.
- calculate_mut_freq() is renamed to calculate_mf()
- total_depth and depth are no longer required for import
- calculate_mf() no longer requires depth. Users may choose to 1) calculate depth from mutation data, 2) supply precalculated depth values in a separate table, 3) No depth, mf is not calculated, only mutation sums.
- plot_spectra, plot_trinucleotide, and spectra_comparison are now supplied with mf_data instead of the mutation data.
- Example data has been added: Currently 44Mb

  
For full details on function utility, see below.

# MutSeqR: Error-corrected Next-Generation Sequencing (ECS) Analysis For Mutagenicity Assessment

Error-corrected next-generation sequencing (ECS) uses various methods to combine multiple independent raw sequence reads derived from an original starting molecule, thereby subtracting out artifacts introduced during sequencing or library preparation. This results in a highly accurate representation of the original molecule. ECS is particularly useful for detecting rare somatic mutations (or mutations induced in germ cells), such as those that arise from mutagen exposure or other sources of DNA damage. ECS is a powerful tool for assessing the mutagenicity of chemicals, drugs, or other agents, and can be used to identify the mutational signatures of these agents. ECS can also be used to detect rare mutations in cancer or other diseases, and to track the clonal evolution of these diseases over time.

For more background on how ECS works and its context in regulatory toxicology testing and genetic toxicology, see the following articles:
- [Menon and Brash, 2023](10.1016/j.mrrev.2023.108471)
- [Marchetti et al., 2023a](https://doi.org/10.1038/d41573-023-00014-y)
- [Marchetti et al., 2023b](https://doi.org/10.1016/j.mrrev.2023.108466)
- [Kennedy et al., 2014](https://doi.org/10.1038/nprot.2014.170)

# Features
 - Import tabular or VCF mutation data
 - Filter variants
 - Summarise mutation frequency across samples, experimental groups, and mutation subtypes
 - Perform statistical analysis of mutation data
 - Visualize data

# Installation

Install the package from github:

```{r}
# install.packages("devtools")
devtools::install_github("EHSRB-BSRSE-Bioinformatics/MutSeqR", auth_token = "your personal_access_token from GitHub")
```

Load the package
```{r}
library(MutSeqR)
```

# Citation
To cite this package in publications use:

Meier M, Dodge A, Williams A, Esina E (2025). MutSeqR: Analysis of Error-Corrected Sequencing Data for Mutation Detection. R package version 0.99.0, https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/.