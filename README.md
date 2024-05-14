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
The main piece of information you want to import is the genome variant file.
Your variant file can be imported as either a `.mut` file using the
function `import_mut_data` or as a `.vcf` file using the function
`read_vcf`.

When first importing a variant file, it is preferred to keep non-variant
rows. This allows the calculuation of mutation frequencies. The data set
can be pared down later to include only mutations of interest (SNVs,
indels, SVs, or any combination). Genome mut and genome vcf files will
provide a row for every position in the interval range, regardless of
whether or not a mutation call was made.

### Variant Filtering
#### Germline Variants
Set the vaf_cutoff to flag ostensibly germline mutations that
have a variant allele fraction greater than this parameter.
The variant allele fraction (VAF) is the fraction of haploid genomes in the
original sample that harbor a specific mutation at a specific base-pair
coordinate of the reference genome. Specifically, is it calculated
by dividing the number of variant reads by the total sequencing depth
at a specific base pair coordinate. The VAF is a good indicator of the
zygosity of a variant. In a typical diploid cell, a homozygous germline
variant will appear on both alleles, in every cell. As such, we expect this
variant to occur on every read - giving us a VAF = 1. A heterozygous
germline variant occurs on one of the two alleles in every cell, as such
we expect this variant to occur on about half of the reads, giving a VAF = 0.5.
Rare somatic variants occur in only a small portion of the cells, thus we 
expect them to appear in only a small percentage of the reads. Typical
VAF values for somatic variants will be less than 0.01 - 0.1. Setting the
vaf_cutoff parameter to 0.01 or 0.1 will flag all variants that have a VAF
greater than this value as germline within the is.germline column. Germline
variants are not included in the mutation counts when calculating mutation
frequencies.

#### Variants within target regions
Supply a regions interval list of genomic ranges of interest and filter
out mutations occuring outside of these regions. If users are targetting or are
interested in a known range of genomic regions, these regions may be specified using
the 'regions' parameter. Any variant that occurs outside of the specified
range will be filtered out of the variant file and returned to the users in
a seperate data frame. This includes variants that partially extend outside of
the regions such as large insertions/deletions or structural variants. Users may
choose to retain some or all of these variants using the 'range_buffer' parameter.
Setting this parameter to an integer will extend the range of the genomic regions
in which a variant can occur by the specified number of base-pairs.

The 'regions' parameter can be set to one of TwinStrand's DuplexSeq™
Mutagenesis Panels; "TSpanel_mouse", "TSpanel_human", or "TSpanel_rat". If you
are using an alternative panel then you may set the 'regions' parameter to 
"custom_interval" and  you will add your target regions' metadata using a custom_regions_file.
Use parameters to indicate your file's file path, delimiter, and whether
the region ranges coordinates are 0-based or 1-based. Mutation data and
region coordinates will be converted to 1-based. If you do not wish to specify
a regions list, then set the 'regions' parameter to "none".

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
                  regions = "TSpanel_mouse", # Twinstand's Mouse Mutagenesis Panel
                  range_buffer = 500 # Retain variants that extend 500 bp outside of the Mouse Mutageneis Panel's target ranges.
                  )
```

If you are using a custom target panel, provide the filepath
to the interval list of region ranges. This can be saved as any file type,
but be sure to specify the proper delimiter using the rg_sep parameter.
Required columns for a custom_regions_file are "contig", "start", and "end".


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
below. If your data contains a column that is synonymous to one of the
required columns, but the name is not included in our synonyms list,
your column name may be substituted using the 'custom_column_names'
parameter. Provide this parameter with a list of names to specify the meaning
of column headers.
```{r}
library(MutSeqR)
# mut_data <- "file path to .mut file"
# sample_data <- "file path to sample meta data"
mutation_data <-
  import_mut_data(mut_file = mut_data,
                  sample_data_file = sample_data,
                  vaf_cutoff = 0.1,
                  regions = "TSpanel_mouse", # Twinstand's Mouse Mutagenesis Panel
                  custom_column_names = list(total_depth = "my_custom_depth_name",
                                             sample = "my_custom_sample_column_name")
                  )
```


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
General usage is the same as for import_mut_file:
Indicate the file path to your .vcf file using the
mut_file parameter. If you have sample metadata, then you can indicate
the file path to your sample data file using the sample_data_file
parameter. Set the vaf_cutoff to flag ostensibly germline mutations that
have a variant allele fraction greater than this parameter. Finally,
load in the metadata for an interval list of genomic target regions 
using the regions parameter. 

```{r}
library(MutSeqR)
# mut_data <- "file path to .vcf file"
# sample_data <- "file path to sample meta data"
mutation_data <-
  import_vcf_data(vcf_file = mut_data,
                  sample_data_file = sample_data,
                  vaf_cutoff = 0.1,
                  regions = "TSpanel_mouse", # Twinstand's Mouse Mutagenesis Panel
                  range_buffer = 500 # Retain variants that extend 500 bp outside of the Mouse Mutageneis Panel's target ranges.
                  )
```

When importing data using a vcf file, the function will retrieve the sequence context information for
each position called in the variant file. The sequence context includes the reference base at the
specified base pair coordinate alongside its two flanking bases. Ex. ACT. If a user specifies a regions
interval file, then the function will retrieve the sequences of the specified genomic intervals from the
USCS database. When using one of TwinStrand's DuplexSeq™ Mutagenesis Panels, the reference genomes are
pre-set to human: GRCh38, mouse: mm10, rat" rn6. If you are supplying a custom_regions_file, then
you must supply the reference genome for your target regions using the 'genome' paramater. This
ensures that the function retrieves the proper sequences to populate the context column.
```{r}
library(MutSeqR)
# mut_data <- "file path to .vcf file"
# sample_data <- "file path to sample meta data"
mutation_data <-
  import_vcf_data(vcf_file = mut_data,
                  sample_data_file = sample_data,
                  vaf_cutoff = 0.1,
                  regions = "custom_interval", 
                  custom_regions_file = "file path to your regions file",
                  rg_sep = "\t", # tab-delimited
                  is_0_based = FALSE, # Ranges are 1-based 
                  genome = "mm10" # Will download target sequences from the mm10 reference genome
                  )
```

If you choose not to supply an interval list of target regions, then you must supply
both the species and the genome assembly version for your reference genome using the
'species' and 'genome' parameters respectively. The function will browse BSgenome
\link[BSgenome]{available.genomes} for the appropriate reference genome and install
the corresponding package. Context information will be extracted from the installed
BSgenome object. BSgenome offers genomes with masked sequences. If you wish to use
the masked version of the genome, set 'masked_BS_genome' to TRUE.

```{r}
library(MutSeqR)
# mut_data <- "file path to .vcf file"
# sample_data <- "file path to sample meta data"
mutation_data <-
  import_vcf_data(vcf_file = mut_data,
                  sample_data_file = sample_data,
                  vaf_cutoff = 0.1,
                  regions = "none", 
                  species = "mouse" # Common name or scientific name is acceptable
                  genome = "mm10", # assembly version
                  masked_BS_genome = FALSE # download the non-masked version of the genome.
                  )
# the function will install BSgenome.Mmusculus.UCSC.mm10 from BioConductor.
```

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

### mutation_data output
The functions will import the variant file(s) as a dataframe, join
it with the metadata, and create some columns that will be helpful
for calculating frequencies in later analyses. A list of the new columns
and their definitions can be found below. The functions will also
make some adjustments to the variation_type column. The following table
displays the categories for the different variation types. Some adjustments
may include, changing "indel" to "insertion" or "deletion".

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


Finally, our functions provide the option to convert
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

Similarly, if the user is using a target panel, they may supply
additional metadata columns in their custom_regions_file that
will be appended to the variant file. Metadata for the
TwinStrand's DuplexSeq™ Mutagenesis Panels include:
genic context, region chromatin state, region GC content, and
the regions' genes.

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

## Generalized Linear Modelling
An important component of analysing mutagencity data is
how mutation frequency changes based on experimental variables.

The 'model_mf' function will fit a generalized linear model to
analyse the effect(s) of given factor(s) on mutation frequency
and perform specified pairwise comparisons between levels of your
factors. Mutation data should first be summarised by sample using
the calculate_mut_freq function. The mf_data should be output as
a summary table. Be sure to retain the columns for experimental
variables of interest using the 'retain_metadata_cols' parameter.

Users may specify factors and covariates for their model using the
'fixed_effects' and 'random_effects' parameters respectively. If more
than one fixed_effect is supplied, then users may specify whether they wish
to test the interaction between their fixed_effects using the
'test_interaction' parameter. 

Users must specify the columns in their mf_data that contain
the mutation counts and the total sequenced bases per sample using
the 'muts' and 'total_counts' parameters respectively. 

By default, the function will fit a generalized linear model with
a quasibinomial distribution. If a random effect is provided than the model
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
We set 'test_interaction' to TRUE to study how the dose response might
change between the different targets.
For more complicated models such as this, we can increase the ____ to
improve convergence by supplying extra arguments directly to the
lme4::glmer function. # Ask Andrew to explain better. 
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

The function will provide model estimates for all levels of the fixed_effects

### Goodness of Fit
The model_mf function will output the model residuals appended to the mf_data.
Additionally, model residuals will be plotted as a histogram and a QQ-plot
so users can ensure a good model fit. We assume that residuals will follow a
normal distribution with a mean of 0.

### Pairwise Comparisons
The model_mf() function will also run specified pairwise comparisons
between the levels of the fixed_effects. The user must supply a constrast
table using the 'contrasts' parameter. This can either be a data frame 
or a file path to a text file. The table must consist of two columns,
each containing groups within the fixed_effects. The group in the first
column will be compared to the group in the second column. Users should
also provide the reference level for each fixed effect using the reference
level parameter. If the user specifies multiple pariwise comparisons, then
the p-values will be corrected using the Sidak method. 

Ex. Going back to our example in which we model the effect
of dose on mutation frequency; let's assume that we have four
dose groups: D1, D2, D3, and a vehicle control D0. The 'reference_level'
will be D0. Using the contrast table, we can specify pairwise comparisons
between each of the doses and the vehicle control (D1 vs. D0, D2 vs. D0, D3 vs. D0).
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            summary = TRUE,
            retain_metadata_cols = "dose"
            )
contrasts_table <- data.frame(col1 = c(D1, D2, D3),
                              col2 = c(D0, D0, D0))
model_by_dose <- model_mf(mf_data = mf_data,
                          fixed_effects = "dose",
                          muts = "sample_sum_min",
                          total_count = "sample_group_depth",
                          reference_level = "D0",
                          contrasts = contrasts_table
                          )
```

For multiple fixed effects, the user must include levels for all fixed_effects
in each value of the contrasts table. Within each value, the levels of the different
fixed_effects should be seperated by a colon.

Ex. Let's go back to our example modelling the effect of dose across multiple
genomic targets. We will define the levels of dose as D0, D1, D2, and D3, with D0
as the reference level. The genomic target factor will have levels chr1 and chr2, 
representing two genomic targets. We will arbitrarily set the reference level as
chr1 for this factor. We will create a contrasts table that compares each dose group to
the control dose D0 for both of the genomic targets. The order in which values occur
for both the reference level and the contrasts should match the order in which
the fixed_effects are listed. In this example "dose" levels will always preceed
"target" levels. 
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = c("sample", "target")
            subtype_resolution = "none",
            summary = TRUE,
            retain_metadata_cols = "dose"
            )
contrast_table <- data.frame(col1 = c("D1:chr1", D2:chr1, "D3:chr1", "D1:chr2", "D2:chr2", "D3:chr2"),
                             col2 = c("D0:chr1", "D0:chr1", "D0:chr1", "D0:chr2", "D0:chr2", "D0:chr2"))
model_by_target <- model_mf(mf_data = mf_data,
  fixed_effects = c("dose", "target"),
  test_interaction = TRUE,
  random_effects = "sample",
  muts = "sample_target_sum_min",
  total_count = "sample_target_group_depth",
  reference_level = c("D0", "chr1"),
  contrasts = contrast_table,
  control = lme4::glmerControl(check.conv.grad = lme4::.makeCC("warning",
                                                               tol = 3e-3,
                                                               relTol = NULL))
  )
```

### Output
The function will output a list of results.
- model_data: the supplied mf_data with added column for model residuals.
- summary: the summary of the model.
- anova: the analysis of variance for models with two or more effects. \link[car]{Anova}`(model) `
- residuals_histogram: the model residuals plotted as a histogram. This is
used to check whether the variance is normally distributed. A symmetric
bell-shaped histogram, evenly distributed around zero indicates that the
normality assumption is likely to be true.
- residuals_qq_plot: the model residuals plotted in a quantile-quantile plot.
 For a normal distribution, we expect points to roughly follow the y=x line.  
- point_estimates_matrix: the contrast matrix used to generate point-estimates for the fixed effects. 
- point_estimates: the point estimates for the fixed effects.
- pairwise_comparisons_matrix: the contrast matrix used to conduct the pairwise comparisons specified in the `contrasts`.
- pairwise_comparisons: the results of pairwise comparisons specified in the `contrasts`.

## Benchmark Dose Modelling
A benchmark dose (BMD) is a dose or concentration that produces a predetermined
change in the response rate of an adverse effect. This predetermined change in
response is called the benchmark response (BMR). In chemical risk assessment
the BMD can be used as point of departure (POD) to derive human health-based guidance 
value such as reference dose (RfD) or derived no-effect level (DNEL) or acceptable 
daily intake (ADI).

The BMD is also a useful for potency comparisons between substances...

The BMD is estimated by applying various mathmatical models to fit the dose-response data.
Some requirements that must be met before modelling the BMD. There must be a clear
dose-response trend in the mutaiton frequency data. We suggest using the 'model-mf'
function to test for significant increases in MF with dose prior to running a BMD analysis.
In general, studies with more dose groups and a graded monotonic response with dose will be
more useful for BMD analysis. A minimum of three dose groups + 1 control group is suggested.
Datasets in which a response is only observed at the high dose are usually not suitable for BMD modeling.
However, if the one elevated response is near the BMR, adequate BMD computation may result. For a better
estimate of the BMD, it is preferable to have studies with one or more doses near the level of the BMR.


Individual vs Summary data
# https://www.epa.gov/sites/default/files/2015-01/documents/benchmark_dose_guidance.pdf
It is preferable to provide information on individual subjects however, it is also possible to 
use summary information (mean + SD) concerning the measured effect, especially for continuous response variables such as mutation frequency.


The BMD is reported alonside its upper and lower confidence intervals; the BMDU and BMDL.
The BMDL is typically used to derive human health-based guidance values.

