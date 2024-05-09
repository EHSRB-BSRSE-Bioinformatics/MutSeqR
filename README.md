<!-- badges: start -->
  [![R-CMD-check](https://github.com/EHSRB-BSRSE-Bioinformatics/duplex-sequencing/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EHSRB-BSRSE-Bioinformatics/duplex-sequencing/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->
  
# MutSeqR: Error-corrected Next-Generation Sequencing (ecNGS) Analysis For Mutagenicity Assessment

## What is ecNGS?

Error-corrected next-generation sequencing (ecNGS) uses various methods to combine multiple independent raw sequence reads derived from an original starting molecule, thereby subtracting out artifacts introduced during sequencing or library preparation. This results in a highly accurate representation of the original molecule. ecNGS is particularly useful for detecting rare somatic mutations (or mutations in germ cells), such as those that arise from mutagen exposure or other sources of DNA damage. ecNGS is a powerful tool for assessing the mutagenicity of chemicals, drugs, or other agents, and can be used to identify the mutational signatures of these agents. ecNGS can also be used to detect rare mutations in cancer or other diseases, and to track the clonal evolution of these diseases over time.

For more background on how ecNGS works and its context in regulatory toxicology testing and genetic toxicology, see the following articles:
- [Brash et al., 2023](10.1016/j.mrrev.2023.108471)
- [Marchetti et al., 2023a](https://doi.org/10.1038/d41573-023-00014-y)
- [Marchetti et al., 2023b](https://doi.org/10.1016/j.mrrev.2023.108466)
- [Kennedy et al., 2014](https://doi.org/10.1038/nprot.2014.170)

This R package is meant to facilitate the import, cleaning, and analysis of ecNGS data, beginning with a table of variant calls or a variant call file (VCF). The package is designed to be flexible and enable users to perform common statistical analyses and visualisations. Currently, it  has been tested primarily with TwinStrand Biosciences' Duplex Sequencing data, but it should be adaptable to other ecNGS methods as well.

## Installation

Install from github with:

```{r}
# install.packages("devtools")
devtools::install_github("EHSRB-BSRSE-Bioinformatics/MutSeqR", auth_token = "your personal_access_token from GitHub")
```

## Data import

The main goal of this package is to generate summary statistics,
visualization, exploratory analysis, and other post-processing tasks
such as mutational signature analysis or generalized linear modeling.
The main piece of information you want to import is the variant file.
Your variant file can be imported as either a `.mut` file using the
function `import_mut_data` or as a `.vcf` file using the function
`read_vcf`.

### Importing .mut files

General usage: Indicate the file path to your .mut file using the
mut_file parameter. If you have sample metadata, then you can indicate
the file path to your sample data file using the sample_data_file
parameter. Set the vaf_cutoff to flag ostensibly germline mutations that
have a variant allele fraction greater than this parameter. Finally,
load in the metadata for an interval list of genomic target regions 
using the regions parameter. 

```{r}
library(MutSeqR)
# mut_data <- "file path to .mut file"
# sample_data <- "file path to sample meta data"
mutation_data <-
  import_mut_data(mut_file = mut_data,
                  sample_data_file = sample_data,
                  vaf_cutoff = 0.1,
                  regions = "TSpanel_mouse" # Twinstand's Mouse Mutagenesis Panel
                  )
```

If you are not using one of TwinStrand's DuplexSeq™ Mutagenesis Panels, then
you will add your target regions' metadata using a custom_regions_file.
Use parameters to indicate your file's file path, delimiter, and whether
the region ranges coordinates are 0-based or 1-based. Mutation data and
region coordinates will be converted to 1-based.

```{r}
# mut <- "file path to .mut file"
# sample_data <- "file path to sample meta data"
mut_data <- import_mut_data(
              mut_file = mut,
              sample_data_file = sample_data,
              vaf_cutoff = 0.1,
              regions = "custom_interval",
              custom_regions_file = "file path to your regions file",
              rg_sep = "\t", # tab-delimited
              is_0_based = FALSE # Ranges are 1-based 
              )
```

Required columns for your `.mut` file are listed in the table below. We
recognize that column names may differ. Therefore, we have implemented
some default column name synonyms. If your column name matches one of
our listed synonyms, it will automatically be changed to match our set
values. For example, your `contig` column may be named `chr` or
`chromosome`. After importing your data, this synonymous column name
will be changed to `contig`. Column names are case-insensitive. A list
of column name synonyms are listed alongside the column definitions
below.

**Table 1.** **Required columns for .mut file import.**
| Column Name      | Definition                               | Synonyms            |
|-------------|------------------------------------------|---------------------|
| contig      | The reference sequence name.              | chr; chromosome; seqnames|
| start       | The 0-based start position of the feature | position           |
| end         | The half-open end position of the feature  |                   |
| sample      | The sample name.                          | sample_name; sample_id |
| ref         | The reference allele at this position.    |                     |
| alt         | The left-aligned, alternate allele at this position. | alt.value|
| alt_depth   | The read depth supporting the alternate allele.| var_depth      |
| depth_col   | The total read depth at this position. This column can be total_depth (excluding N-calls) or depth (including N-calls; if total_depth is not available).| informative_somatic_depth = total_depth|
| variation_type | The category to which this variant is assigned.| type; mut_type; variant_type|
| context     | The local reference trinucleotide context at this position (e.g. ATC - not necessarily the transcript codon).| sequence_context; flanking_sequence   |

### Importing .vcf files

Required fields for your `.vcf` file are listed in the table below.

|              | **Field Name** | **Definition** |
|--------------|-----------------|-----------------|
| **FIXED FIELDS** | `CHROM` | The reference sequence name. |
|              | `POS` | The 0-based start position of the feature in contig. |
|              | `REF` | The reference allele at this position. |
|              | `ALT` | The left-aligned, alternate allele at this position. |
| **FORMAT FIELDS** | `AD` | The allelic depths for the reference and alternate alleles in the order listed. |
|              | `DP` | The total read depth at this position (including N-calls). Equivalent to `depth`. |
|              | `VD` | Variant Depth. Equivalent to `alt_depth`. |
| **INFO FIELDS** | `TYPE` | The category to which this variant is assigned. Equivalent to `variation_type`. |
|              | `END` | The half-open end position of the feature in contig. |
| *SUGGESTED INFO FIELDS* | `sample` | An identifying field for your samples; either in the INFO field or as the header to the FORMAT field. |
|              | `SVTYPE` | Structural variant types; INV DUP DEL INS FUS. |
|              | `SVLEN` | Length of the structural variant in base pairs. |


The column variation_type/TYPE may contain these values:
| `variation_type` | Definition                                          |
|------------------|-----------------------------------------------------|
| no_variant       | No variation, the null-case.                        |
| snv              | Single nucleotide variant.                          |
| mnv              | Multiple nucleotide variant.                        |
| insertion        | Insertion, length of REF = 1bp.                     |
| deletion         | Deletion, length of ALT = 1bp.                      |
| complex          | Length of REF and ALT differ and are both > than 1 bp. |
| symbolic         | Structural variant or IUPAC ambiguity code.         |


It is critical that the variant file gets annotated with information
such as the genomic region from which the mutation originates, some
information about the sample from which the mutation originates, and
that each row (mutation) receives a calculated group depth and group
frequency (of interest for later analyses).

We do this by first importing the variant file as a data frame, and then
joining it with the metadata (e.g., sample_data_file, regions file), and
creating some columns that will be helpful for calculating frequencies
in later analyses. Finally, our functions provide the option to convert
the resulting data frame into a `granges` object. This facilitates use
in other packages and makes doing 'genome math' on the ranges
significantly easier.

Columns that are added to the resulting data frame are listed below.

| Column Name        | Definition                                        |
|--------------------|---------------------------------------------------|
| `nchar_ref`        | The length (in bp) of the reference allele.       |
| `nchar_alt`        | The length (in bp) of the alternate allele.       |
| `varlen`           | The length (in bp) of the variant.                |
| `total_depth`      | The total read depth at this position, excluding N-calls. |
| `vaf`              | The variant allele fraction. Calculated as `alt_depth`/`depth_col` where `depth_col` can be `total_depth` or `depth`. |
| `is_germline`      | TRUE or FALSE. Flags ostensible germline mutations (`vaf` > `vaf_cutoff`). |
| `ref_depth`        | The total read depth at the position calling for the reference allele. Calculated as `depth_col` - `alt_depth` where `depth_col` can be `total_depth` or `depth`. |
| `subtype`          | The substitution type for the snv variant (12-base spectrum; e.g., A>C). |
| `short_ref`        | The reference base at this position.              |
| `normalized_subtype` | The C/T-based substitution type for the snv variant (6-base spectrum; e.g., A>C -> T>G). |
| `normalized_ref`   | The reference base in C/T-base notation for this position (e.g., A -> T). |
| `context_with_mutation` | The substitution type for the snv variant including the two flanking nucleotides (192-trinucleotide spectrum; e.g., T[A>C]G). |
| `normalized_context_with_mutation` | The C/T-based substitution type for the snv variant including the two flanking nucleotides (96-base spectrum e.g., T[A>C]G -> C[T>G]A). |
| `normalized_context` | The trinucleotide context in C/T base notation for this position (e.g., TAG -> CTA). |
| `gc_content`       | % GC of the trinucleotide context at this position. |


### Metadata: an important consideration

The other important component of importing your data for proper use is
to assign each mutation to a biological sample, and also make sure that
some additional information about each sample is present (e.g., a
chemical treatment, a dose, etc.). This is done by providing a sample
data file (tab delimited, comma delimited, etc.; the choice is up to the
user, but the delimiter of the file must be specified as a parameter in
the function). Importantly, this is a file that would be analogous to
"colData", or "column data", a term often used in the `DESeq2` package.
Hence, it must contain some information about an existing column in your
variant file, which is typically going to be sample. So the first column
in your sample data file should indeed be `sample`. Then, additional
columns such as `dose` or `tissue` or `treatment` can be added, and
these columns will be joined with your variant file to capture that
information and associate it with each mutation.

### Other Notes

When first importing a `.mut` file, it is preferred to keep non-variant
rows. This allows the calculuation of mutation frequencies. The data set
can be pared down later to include only mutations of interest (SNVs,
indels, SVs, or any combination).

## Calculating Mutation Frequencies
The function calculate_mut_freq() summarises the mutation counts
across arbitrary groupings within the mutation data. Mutations
can be summarised across samples, experimental groups, and mutation
subtypes for later statistical analyses. Mutation frequency is calculated
by dividing the number of mutations by the total number of sequenced bases
in each group. The units for mutation frequency are mutations/bp.

### Grouping Mutations
Mutation counts and total sequenced bases are summed within
groups that can be designated using the cols_to_group parameter.
This parameter can be set to one or more columns in the mutation data
that represents experimental variables of interest. 

The following will return mutation counts and frequencies summed
across samples.
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "none"
            )
```

Alternatively, you can sum mutations by experimental groups
such as 'dose' or 'tissue', or both at the same time. Counts
and frequencies will be returned for every level of the designated
groups.
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = c("dose", "tissue"),
            subtype_resolution = "none"
            )
```

### Mutation Subtypes
Mutations can also be grouped by mutation subtype at varying
degrees of resolution using the 'subtype_resolution' parameter.
Mutations and total sequenced bases will be summed across groups
for each mutation subtype. The total number of sequenced bases
is calculated based on the sequence context in which a mutation
subtype occurs. For instance, C>T mutations will only occur
at positions with a C reference. Therefore,
the mutation frequency for C>T mutations is calculated as
the total number of C>T mutations divided by the total number
bases sequenced at a C reference position for a particular sample/group.

The function will also calculate the the proportion
of mutations for each subtype. The proportion
of mutations is calculated by dividing the number of mutations for
each subtype by the total number of mutations within a sample or group.
Proportions are then normalized to the sequencing depth. First proportions
are divided by the context-dependent number of sequenced bases for that sample/group.
Values are then divided by the sum of all values for that sample/group.

Subtype resolutions: 
1. type: the variation type. "snv", "mnv", "insertion", "deletion", "complex", and "symbolic" variants.
2.  base_6: the simple snv spectrum. The snv subtypes are normalized to their pyrimidine context. C>A, C>G, C>T, T>A, T>C, T>G.
3. base_12: The non-normalized snv subtypes. A>C, A>G, A>T, C>A, C>G, C>T, G>A, G>C, G>T, T>A, T>C, T>G.
4. base_96: the trinculeotide spectrum. The snv subtypes are reported in their pyrimidine context alongside their two flanking nucleotides. Ex. A[C>T]A.
5. base_192: The non-normalized snv subtypes are reported alongside their two flanking nucleotides. Ex. A[G>T]A.

Ex. The following code will return the simple mutation spectra for all samples. 
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "base_6"
            )
```

This function will also calculate mutation frequencies and proportions
on a specific subset of variation types, which can be set using
the 'variant_types' parameter. The 'variant_types' parameter can be
set to a character string of "types" values that the user wants included
in the mutation counts. By default the function will calculate
summary values based on all mutation types.

Ex. The following code will calculate mutation frequencies per sample
for only insertion and deletion mutations. 
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            variant_types = c("insertion", "deletion")
            )
```
Alternatively, the following code will return the mutation frequencies 
and proportions for single-nucleotide variants (snv) only, differentiated
into their trinucleotide spectra.
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "base_96",
            variant_types = "snv"
            )
```
### Mutation Counting Methods
 Mutations are counted based on two opposing assumptions. Frequencies and
 proportions are returned for both options.

The Minimum Independent (min) Mutation Counting Method counts each mutation
once, regardless of the number of reads that contain the non-reference
allele. This method assumes that multiple instances of the same 
mutation within a sample/library are the result of clonal expansion of a
single mutational event. This is likely an undercount of mutations 
because we expect some mutations to recur in mutation hotspot regions.

The Maximum Independent (max) Mutation Counting Method counts multiple
identical mutations at the same position within a sample/library as
independent mutation events. This is likely an overcount of mutations
since we do expect some recurrent mutations to arise through clonal
expansion. 

The Minimum Independent counting method is generally recommended
for characterising rare somatic mutations because the Maximum
Independent method tends to increase the sample variance of mutation
frequencies by a significant degree. 

### Variant Filtering
By default, germline variants will not be included the sumarised mutation counts, frequencies, or 
proportions. Germline mutations are identified in the mutation data using the
is_germline column. Users may choose to include them in their mutation counts
by setting the 'filter_germ' parameter to FALSE.

### Summary Table
The function will output the resulting mf_data as a data frame with the mutation frequency and proportion
calculated. If the 'summary' parameter is set to TRUE, the data frame will be a
summary table with the mutation frequency calculated for each group. If summary
is set to FALSE, the mutation frequency will be appended to each row of the original
mutation_data.

The summary table will include:
- cols_to__group: all columns used to group the data.
- _sum_: the min/max mutation counts for the group.
-_MF_: the min/max mutation frequency for the group.
_proportion_: the min/max proportion for the group.

Additional columns from the orginal mutation data can be retained using the 
'retain_metadata_cols' parameter. Retaining higher-order experimental groups
may be useful for later statistical analyses or plotting.

Ex. The following code will calculate the mutation frequencies for each sample
and retain the dose column for each sample to use in later analyes.
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            retain_metadata_cols = "dose"
            )
```

## Linear Modelling
The 'model_mf' function will fit a linear model to analyse the effect(s) of given factor(s) 
on mutation frequency and perform specified pairwise comparisons. Mutation data should first
be summarised by sample using the calculate_mut_freq function. The mf_data should be output as
a summary table. Be sure to retain the columns for experimental variables of interest using the
'retain_metadata_cols' parameter.

Users may specify factors and covariates for their model using the
'fixed_effects' and 'random_effects' parameters respectively. If more
than one fixed_effect is supplied, then users may specify whether they wish
to test the interaction between their fixed_effects using the
test_interaction parameter. 

Users must specify the columns in their mf_data that contain
the mutation counts and the total sequenced bases per sample using
the 'muts' and 'total_counts' parameters respectively. 

By default, the function will fit a generalized linear model with
a quasibinomial distribution. If a random effect is provided then the model
will fit a general linear mixed model with a binomial distribution. The 
dispersion family for the model can be customized using the 'family'
parameter.

Ex. The following code will fit a generalized linear model to study the
effect of dose on mutation frequency.
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            summary = TRUE,
            retain_metadata_cols = "dose"
            )

model_by_dose <- model_mf(mf_data = mf_data,
                          fixed_effects = "dose",
                          muts = "sample_sum_min",
                          total_count = "sample_group_depth"
                          )
```

Additional arguments can be passed to the model to further customize it to
the user's needs. Details on the arguments for the generalized linear
model can be found here \link[stats]{glm} and for the general linear mixed
model here \link[lme4]{glmer}. 

Ex. We can study the effects of dose on mutation frequency for individual
genomic loci from a panel of targets.

First, we summarise mutations by sample and by genomic target.
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = c("sample", "target")
            subtype_resolution = "none",
            summary = TRUE,
            retain_metadata_cols = "dose"
            )
```
We will then fit a general linear mixed model of mutation frequencies
using dose and target as fixed effects and sample as the random effect.
For more complicated models, we can increase the ____ to improve convergence. # Ask Andrew to explain better. 
```{r}
model_by_target <- model_mf(mf_data = mf_data,
  fixed_effects = c("dose", "target"),
  test_interaction = TRUE,
  random_effects = "sample",
  muts = "sample_target_sum_min",
  total_count = "sample_target_group_depth",
  control = lme4::glmerControl(check.conv.grad = lme4::.makeCC("warning",
                                                               tol = 3e-3,
                                                               relTol = NULL))
  )
```