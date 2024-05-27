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

This R package is meant to facilitate the import, cleaning, and analysis of ecNGS data, beginning with a table of variant calls or a variant call file (VCF). The package is designed to be flexible and enable users to perform common statistical analyses and visualisations. Currently, it has been tested primarily with TwinStrand Biosciences' Duplex Sequencing data, but it should be adaptable to other ecNGS methods as well.

## Installation

Install from github with:

```{r}
# install.packages("devtools")
devtools::install_github("EHSRB-BSRSE-Bioinformatics/MutSeqR", auth_token = "your personal_access_token from GitHub")
```

## Data import

The main goal of this package is to generate summary statistics, visualization, exploratory analysis, and other post-processing tasks such as mutational signature analysis or generalized linear modeling. The main piece of information you want to import is the **genome variant file**. Your variant file can be imported as either a `.mut` file using the function `import_mut_data` or as a `.vcf` file using the function `read_vcf`.

When first importing a variant file, it is preferred to keep non-variant rows. This allows the calculuation of mutation frequencies. The data set can be pared down later to include only mutations of interest (SNVs, indels, SVs, or any combination). Genome mut and genome vcf files will provide a row for every position in the interval range, regardless of whether or not a mutation call was made.

### Variant Filtering
#### Germline Variants
Set the `vaf_cutoff` to flag ostensibly germline mutations that have a **variant allele fraction** greater than this parameter. The variant allele fraction (VAF) is the fraction of haploid genomes in the original sample that harbor a specific mutation at a specific base-pair coordinate of the reference genome. Specifically, is it calculated by dividing the number of variant reads by the total sequencing depth at a specific base pair coordinate. The VAF is a good indicator of the zygosity of a variant. In a typical diploid cell, a homozygous germline variant will appear on both alleles, in every cell. As such, we expect this variant to occur on every read - giving us a VAF = 1. A heterozygous germline variant occurs on one of the two alleles in every cell, as such we expect this variant to occur on about half of the reads, giving a VAF = 0.5. Rare somatic variants occur in only a small portion of the cells, thus we expect them to appear in only a small percentage of the reads. Typical VAF values for somatic variants will be less than 0.01 - 0.1. Setting the `vaf_cutoff` parameter to 0.01 or 0.1 will flag all variants that have a VAF **greater** than this value as germline within the `is.germline` column. Germline variants are not included in the mutation counts when calculating mutation frequencies.

#### Variants within target regions
Supply a regions interval list of genomic ranges of interest and filter out mutations occuring outside of these regions. If you are targetting or are interested in a known range of genomic regions, these regions may be specified using the `regions` parameter. Any variant that occurs outside of the specified ranges will be filtered out of the variant file and returned in a seperate data frame. This includes variants that partially extend outside of the regions such as large insertions/deletions or structural variants. You may choose to retain some or all of these variants using the `range_buffer` parameter. Setting this parameter to an integer will extend the range of the genomic regions in which a variant can occur by the specified number of base-pairs.

The `regions` parameter can be set to one of TwinStrand's DuplexSeq™ Mutagenesis Panels; *TSpanel_mouse*, *TSpanel_human*, or *TSpanel_rat*. If you are using an alternative panel then you may set the `regions` parameter to  "custom_interval" and  you will add your target regions' metadata using a `custom_regions_file`. Use parameters to indicate your file's file path, delimiter, and whether the region ranges coordinates are 0-based or 1-based. Mutation data and region coordinates will be converted to 1-based. If you do not wish to specify a regions list, then set the `regions` parameter to *none*.

### Importing .mut files

`import_mut_data` General usage: Indicate the file path to your **.mut** file using the `mut_file` parameter. You may provide sample metadata uing the `sample_data_file` parameter. This parameter can take a R data frame, or it can read in a file if provided with a filepath. If using a filepath, specify the proper delimiter using the `sd_sep` parameter. Set the `vaf_cutoff` to flag ostensibly germline mutations that have a variant allele fraction greater than this parameter. Finally, load in the metadata for an interval list of genomic target regions  using the `regions` parameter. 

```{r}
library(MutSeqR)
# mut_data <- "file path to .mut file"
sample_data <- data.frame(sample = c("sample_1", "sample_2", "sample_3", "sample_4"),
                          dose = c("control", "control", "treated", "treated"))
mutation_data <-
  import_mut_data(mut_file = mut_data,
                  sample_data_file = sample_data,
                  vaf_cutoff = 0.1,
                  regions = "TSpanel_mouse", # Twinstand's Mouse Mutagenesis Panel
                  range_buffer = 500 # Retain variants that extend 500 bp outside of the Mouse Mutageneis Panel's target ranges.
                  )
```

If you are using a custom target panel, provide the the interval list of region ranges in `custom_regions_file`. This can be either a filepath or a R dataframe. If using a filepath, this can be saved as any file type, but be sure to specify the proper delimiter using the `rg_sep` parameter. Required columns for a `custom_regions_file` are *contig*, *start*, and *end*. Additional columns will be appended to the mutation data. It is recommended to include a column that provides a unique identifier for each genomic target. In this way, mutation counts can easily be summarised across targets for region-based analysis.

```{r}
# mut <- "file path to .mut file"
# sample_data <- "file path to sample meta data"
my_regions <- data.frame(contig = c("chr1", "chr1", "chr2", "chr3"),
                         start = c(69304218, 155235939, 50833176, 109633161),
                         end = c(69306617, 155238338, 50835575, 109635560),
                         label = c("target1", "target2", "target3", "target4"), #unique identifier for targets
                         genome = c("mm10", "mm10", "mm10", "mm10"),
                         genic_context = c("intergenic", "genic", "intergenic", "genic"))
mut_data <- import_mut_data(
              mut_file = mut,
              sample_data_file = sample_data,
              vaf_cutoff = 0.1,
              regions = "custom_interval",
              custom_regions_file = my_regions,
              is_0_based = FALSE # Ranges are 1-based 
              )
```

Required columns for your `.mut` file are listed in the table below. We recognize that column names may differ. Therefore, we have implemented some default column name synonyms. If your column name matches one of our listed synonyms, it will automatically be changed to match our set values. For example, your `contig` column may be named `chr` or `chromosome`. After importing your data, this synonymous column name will be changed to `contig`. Column names are case-insensitive. A list of column name synonyms are listed alongside the column definitions below. If your data contains a column that is synonymous to one of the required columns, but the name is not included in our synonyms list, your column name may be substituted using the `custom_column_names` parameter. Provide this parameter with a list of names to specify the meaning of column headers.
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
`import_vcf_data` General usage is the same as for import_mut_file: Indicate the file path to your **.vcf** file using the `vcf_file` parameter. If you have sample metadata, then you can indicate the file path to your sample data file using the `sample_data_file` parameter. Set the `vaf_cutoff` to flag ostensibly germline mutations that have a variant allele fraction greater than this parameter. Finally, load in the metadata for an interval list of genomic target regions using the `regions` parameter. 

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

When importing data using a vcf file, the function will retrieve the sequence context information for each position called in the variant file. The sequence context includes the reference base at the specified base pair coordinate alongside its two flanking bases. *Ex. ACT*. If a user specifies a regions interval file, then the function will retrieve the sequences of the specified genomic intervals from the USCS database. When using one of TwinStrand's DuplexSeq™ Mutagenesis Panels, the reference genomes are pre-set to human: GRCh38, mouse: mm10, rat: rn6. If you are supplying a `custom_regions_file`, then you must supply the reference genome for your target regions using the `genome` parameter. This ensures that the function retrieves the proper sequences to populate the context column.

```{r}
library(MutSeqR)
# mut_data <- "file path to .vcf file"
# sample_data <- "file path to sample meta data"
# my_regions <- "file path to your regions file"

mutation_data <-
  import_vcf_data(vcf_file = mut_data,
                  sample_data_file = sample_data,
                  vaf_cutoff = 0.1,
                  regions = "custom_interval", 
                  custom_regions_file = my_regions,
                  rg_sep = "\t", # tab-delimited
                  is_0_based = FALSE, # Ranges are 1-based 
                  genome = "mm10" # Will download target sequences from the mm10 reference genome
                  )
```

If you choose not to supply an interval list of target regions, then you must supply both the species and the genome assembly version for your reference genome using the `species` and `genome` parameters respectively. The function will browse [BSgenome::available.genomes](https://www.rdocumentation.org/packages/BSgenome/versions/1.40.1/topics/available.genomes) for the appropriate reference genome and install the corresponding package. Context information will be extracted from the installed BSgenome object. BSgenome offers genomes with masked sequences. If you wish to use the masked version of the genome, set `masked_BS_genome` to `TRUE`.

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

### Mutation_data Output
Both import functions will import the variant file(s) as a dataframe, join it with the metadata, and create some columns that will be helpful for calculating frequencies in later analyses. A list of the new columns and their definitions can be found below. The functions will also make some adjustments to the `variation_type` column. The following table displays the categories for the different variation types. Some adjustments may include, changing "indel" to "insertion" or "deletion".

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


Finally, our functions provide the option to convert the resulting data frame into a `granges` object. This facilitates use in other packages and makes doing "genome math" on the ranges significantly easier.

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
An important component of importing your data for proper use is to assign each mutation to a biological sample, and also make sure that some additional information about each sample is present (e.g., a chemical treatment, a dose, etc.). This is done by providing a sample data file. Importantly, it must contain some information about an existing column in your variant file, which is typically going to be sample. So the first column in your sample data file should indeed be `sample`. Then, additional columns such as `dose` or `tissue` or `treatment` can be added, and these columns will be joined with your variant file to capture that information and associate it with each mutation.

Similarly, if you are using a target panel, they may supply additional metadata columns in their `custom_regions_file` that will be appended to the variant file. Metadata for the TwinStrand's DuplexSeq™ Mutagenesis Panels include: genic context, region chromatin state, region GC content, and the regions' genes.

## Calculating Mutation Frequencies
The function `calculate_mut_freq` summarises the mutation counts across arbitrary groupings within the mutation data. Mutations can be summarised across samples, experimental groups, and mutation subtypes for later statistical analyses. Mutation frequency is calculated by dividing the number of mutations by the total number of sequenced bases in each group. The units for mutation frequency are mutations/bp.

### Grouping Mutations
Mutation counts and total sequenced bases are summed within groups that can be designated using the `cols_to_group` parameter. This parameter can be set to one or more columns in the mutation data that represent experimental variables of interest. 

The following will return mutation counts and frequencies summed across samples.
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "none"
            )
```

Alternatively, you can sum mutations by experimental groups such as `dose` or `tissue`, or both at the same time. Counts and frequencies will be returned for every level of the designated groups.
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = c("dose", "tissue"),
            subtype_resolution = "none"
            )
```

### Mutation Subtypes
Mutations can also be grouped by mutation subtype at varying degrees of resolution using the `subtype_resolution` parameter. Mutations and total sequenced bases will be summed across groups for each mutation subtype. The total number of sequenced bases is calculated based on the sequence context in which a mutation subtype occurs. For instance, C>T mutations will only occur at positions with a C reference. Therefore, the mutation frequency for C>T mutations is calculated as the total number of C>T mutations divided by the total number bases sequenced at a C reference position for a particular sample/group.

The function will also calculate the the proportion of mutations for each subtype. The proportion of mutations is calculated by dividing the number of mutations for each subtype by the total number of mutations within a sample or group. Proportions are then normalized to the sequencing depth. First proportions are divided by the context-dependent number of sequenced bases for that sample/group. Values are then divided by the sum of all values for that sample/group.

**Subtype resolutions:** 

* type: the variation type.
    + "snv", "mnv", "insertion", "deletion", "complex", and "symbolic" variants.
*  base_6: the simple snv spectrum. The snv subtypes are normalized to their pyrimidine context.
    + C>A, C>G, C>T, T>A, T>C, T>G.
* base_12: The non-normalized snv subtypes.
    + A>C, A>G, A>T, C>A, C>G, C>T, G>A, G>C, G>T, T>A, T>C, T>G.
* base_96: the trinculeotide spectrum. The snv subtypes are reported in their pyrimidine context alongside their two flanking nucleotides.
    + Ex. A[C>T]A.
* base_192: The non-normalized snv subtypes are reported alongside their two flanking nucleotides.
    + Ex. A[G>T]A.

*Ex. The following code will return the simple mutation spectra for all samples.*
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "base_6"
            )
```

This function will also calculate mutation frequencies and proportions on a specific subset of variation types, which can be set using the `variant_types` parameter. The `variant_types` parameter can be set to a character string of "types" values that the user wants included in the mutation counts. By default the function will calculate summary values based on all mutation types.

*Ex. The following code will calculate mutation frequencies per sample for only insertion and deletion mutations.*
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            variant_types = c("insertion", "deletion")
            )
```

*Ex. Alternatively, the following code will return the mutation frequencies  and proportions for single-nucleotide variants (snv) only, differentiated into their trinucleotide spectra.*
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "base_96", # trinucleotide resolution
            variant_types = "snv" # include on snv mutations
            )
```

### Mutation Counting Methods
 Mutations are counted based on two opposing assumptions. Frequencies and proportions are returned for both options.

The **Minimum Independent Mutation Counting Method** (min) counts each mutation once, regardless of the number of reads that contain the non-reference allele. This method assumes that multiple instances of the same  mutation within a sample/library are the result of *clonal expansion of a single mutational event*. This is likely an undercount of mutations because we expect some mutations to recur in mutation hotspot regions.

The **Maximum Independent Mutation Counting Method** (max) counts multiple identical mutations at the same position within a sample/library as
*independent mutation events*. This is likely an overcount of mutations since we do expect some recurrent mutations to arise through clonal expansion. 

The Minimum Independent counting method is generally recommended for characterising rare somatic mutations because the Maximum Independent method tends to increase the sample variance of mutation frequencies by a significant degree. 

### Variant Filtering
By default, germline variants will not be included in the sumarised mutation counts, frequencies, or  proportions. Germline mutations are identified in the mutation data using the `is_germline` column. Users may choose to include them in their mutation counts by setting the `filter_germ` parameter to FALSE.

### Summary Table
The function will output the resulting `mf_data` as a data frame with the mutation frequency and proportion calculated. If the `summary` parameter is set to `TRUE`, the data frame will be a summary table with the mutation frequency calculated for each group. If `summary` is set to `FALSE`, the mutation frequency will be appended to each row of the original `mutation_data`.

The summary table will include:
* `cols_to_group`: all columns used to group the data.
* `_sum_`: the min/max mutation counts for the group.
* `_MF_`: the min/max mutation frequency for the group.
* `_proportion_`: the min/max proportion for the group.

Additional columns from the orginal mutation data can be retained using the `retain_metadata_cols` parameter. Retaining higher-order experimental groups may be useful for later statistical analyses or plotting.

*Ex. The following code will calculate the mutation frequencies for each sample and retain the dose column for each sample to use in later analyes.*
```{r}
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            retain_metadata_cols = "dose"
            )
```

## Generalized Linear Modelling
An important component of analysing mutagencity data is how mutation frequency changes based on experimental variables.

The `model_mf` function will fit a generalized linear model to analyse the effect(s) of given factor(s) on mutation frequency and perform specified pairwise comparisons between levels of your factors. Mutation data should first be summarised by sample using the `calculate_mut_freq` function. The `mf_data` should be output as a summary table. Be sure to retain the columns for experimental
variables of interest using the `retain_metadata_cols` parameter.

You may specify factors and covariates for your model using the `fixed_effects` and `random_effects` parameters respectively. If more
than one `fixed_effect` is supplied, then you may specify whether you wish to test the interaction between your fixed_effects using the `test_interaction` parameter. 

You must specify the columns in your `mf_data` that contain the mutation counts and the total sequenced bases per sample using the `muts` and `total_counts` parameters respectively. 

By default, the function will fit a generalized linear model with a quasibinomial distribution. If a random effect is provided than the model will fit a general linear mixed model with a binomial distribution. The  dispersion family for the model can be customized using the `family`
parameter.

*Ex. The following code will fit a generalized linear model to study the effect of dose on mutation frequency.*
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

Additional arguments can be passed to the model to further customize it to your needs. Details on the arguments for the generalized linear model can be found here [stats::glm](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm) and for the general linear mixed model here [lme4::glmer](https://www.rdocumentation.org/packages/lme4/versions/1.1-35.3/topics/glmer). 

*Ex. We can study the effects of dose on mutation frequency for individual genomic loci from a panel of targets.*
```{r}
# Summarise mutations by sample and by genomic target.
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = c("sample", "target")
            subtype_resolution = "none",
            summary = TRUE,
            retain_metadata_cols = "dose"
            )
```
*We will then fit a general linear mixed model of mutation frequencies using dose and target as fixed effects and sample as the random effect. We set `test_interaction` to `TRUE` to study how the dose response might change between the different targets. For more complicated models such as this, we can increase the ____ to improve convergence by supplying extra arguments directly to the lme4::glmer function. # Ask Andrew to explain better.*
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
The `model_mf` function will output the model residuals appended to the `mf_data`. Additionally, model residuals will be plotted as a histogram and a QQ-plot so you can ensure a good model fit. We assume that residuals will follow a normal distribution with a mean of 0.

### Pairwise Comparisons
The `model_mf` function will also run specified pairwise comparisons between the levels of the `fixed_effects`. You must supply a constrast table using the `contrasts` parameter. This can either be a data frame  or a file path to a text file. The table must consist of two columns, each containing groups within the `fixed_effects`. The group in the first column will be compared to the group in the second column. You should also provide the reference level for each fixed effect using the `reference_level` parameter. If you specify multiple pairwise comparisons, then the p-values will be corrected for multiple comparisons using the Sidak method. 

*Ex. Let's go back to our example in which we model the effect of dose on mutation frequency. Let's assume that we have four dose groups: D1, D2, D3, and a vehicle control D0. The `reference_level` will be D0. Using the contrast table, we can specify pairwise comparisons between each of the doses and the vehicle control (D1 vs. D0, D2 vs. D0, D3 vs. D0).*
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

For multiple fixed effects, the user must include levels for all `fixed_effects` in each value of the contrasts table. Within each value, the levels of the different `fixed_effects` should be seperated by a colon.

*Ex. Let's go back to our example modelling the effect of dose across multiple genomic targets. We will define the levels of dose as D0, D1, D2, and D3, with D0 as the reference level. The genomic target factor will have levels chr1 and chr2, representing two genomic targets. We will arbitrarily set the reference level as chr1 for this factor. We will create a contrasts table that compares each dose group to the control dose D0 for both of the genomic targets. The order in which values occur for both the reference level and the contrasts should match the order in which the fixed_effects are listed. In this example "dose" levels will always preceed "target" levels.* 
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
model_by_target <- model_mf(
  mf_data = mf_data,
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

* model_data: the supplied mf_data with added column for model residuals.
* summary: the summary of the model.
* anova: the analysis of variance for models with two or more effects. See [car::Anova.](https://www.rdocumentation.org/packages/car/versions/1.0-9/topics/Anova)
* residuals_histogram: the model residuals plotted as a histogram. This is used to check whether the variance is normally distributed. A symmetric bell-shaped histogram, evenly distributed around zero indicates that the normality assumption is likely to be true.
* residuals_qq_plot: the model residuals plotted in a quantile-quantile plot. For a normal distribution, we expect points to roughly follow the y=x line.  
* point_estimates_matrix: the contrast matrix used to generate point-estimates for the fixed effects. 
* point_estimates: the point estimates for the fixed effects.
* pairwise_comparisons_matrix: the contrast matrix used to conduct the pairwise comparisons specified in the `contrasts`.
* pairwise_comparisons: the results of pairwise comparisons specified in the `contrasts`.

## Benchmark Dose Modelling
A **benchmark dose** (BMD) is a dose or concentration that produces a predetermined change in the response rate of an adverse effect. This predetermined change in response is called the **benchmark response** (BMR). In chemical risk assessment the BMD can be used as a point of departure (POD) to derive human health-based guidance  values such as the reference dose (RfD), the derived no-effect level (DNEL) or the acceptable daily intake (ADI).

The BMD is estimated by applying various mathmatical models to fit the dose-response data. Some requirements must be met before modelling the BMD. There must be a clear dose-response trend in the mutation frequency data. We suggest using the `model_mf` function to test for significant increases in MF with dose prior to running a BMD analysis. In general, studies with more dose groups and a graded monotonic response with dose will be more useful for BMD analysis. A minimum of three dose groups + 1 control group is suggested. Datasets in which a response is only observed at the high dose are usually not suitable for BMD modeling. However, if the one elevated response is near the BMR, adequate BMD computation may result. For a better estimate of the BMD, it is preferable to have studies with one or more doses near the level of the BMR.

MutSeqR can perform benchmark dose modeling of mutation frequencies using the [ToxicR](https://github.com/NIEHS/ToxicR) package, available on Github.

To install this package, use the following code:
``` {r}
library(devtools)
install_github("NIEHS/ToxicR")
```
*See the ToxicR repository on github for more information on installing the package. Differences apply for mac users*

We have two functions available to users:

* `mf_bmd` will fit a single continuous BMD model to the mutation frequency data.
* `bmd_ma` will fit a model average continuous BMD to the mutation frequency data. 

Protection and safety authorities recommend the use of model averaging to determine the benchmark dose. Model averaging incorporates information across multiple models to acount for model uncertainty. In most cases, this allows the BMD to be more accurately estimated.

### Choosing your BMR
One of the most important considerations for BMD modeling is choosing the appropriate benchmark response (BMR). The BMD will be estimated as the dose at which the BMR occurs.  Selecting a BMR involves making judgements about the statistical and biological characteristics of the dataset and about the applications for which the resuling BMDs will be used. There are several different definitions of the BMR. Our functions offer several options  that are commonly used for continuous data:

* Relative deviation (*rel*): the BMD represents the dose that changes the mean mutation frequency a certain percentage from the background dose. 
* Standard deviation (*sd*): the BMD represents the dose associated with the mean mutation frequency changing a specified number of standard deviations from the background mean. 
* Absolute deviation (*abs*): the  BMD represents the dose associated with a specified absolute deviation from the background mean. 
* Hybrid deviation (*hybrid*): the  BMD represents the dose that changes the probability of an adverse event by a specified amount. 

One of these options can be specified using the `bmr_type` parameter. The `bmr` parameter is set to a numeric value specifying the benchmark response, defined in relation to the calculation requested in `bmr_type`.

```{r}
# summarise mutation frequencies by sample
# retain the dose column in the summary table
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            summary = TRUE,
            retain_metadata_cols = "dose"
            )
# Fit a model-averaged BMD to the mutation frequency data.
bmd <- bmd_ma(mf_data,
  # specify the name of your dose column
              dose_col = "dose", 
  # specify the name of your response column(s)
              response_cols = c("sample_MF_min", "sample_MF_max"),
  # BMR is specified as a 50% relative increase from the background mean
              bmr_type = "rel",
              bmr = 0.5,
              ...)
```

Ideally, the BMR would be based on a consensus scientific definition of what  minimal level of change in mutation frequency is biologically significant. Currently, the default provided by this package calculates the BMD at a 50% relative increase in mutation frequency from the background. This BMR was selected based on previous recommendations for genotoxicity assessment by White et al., 2020.

### Models
Model averaging highly depends on the set of candidate models used. A sufficiently large set of models is needed to ensure that a well-fitting model is included in the averaging. The `bmd_ma` function uses the default EFSA models to average. These models are (normal then lognormal for each model): `exp-aerts`, `invexp-aerts`, `hill-aerts`, `lognormal-aerts`, `gamma-efsa`, `LMS`, `probit-aerts`, and `logistic-aerts`.

When using the `mf_bmd` function, the `model_type` parameter specifies the model that will be fit to the data. All EFSA models can be specified. Additionally, legacy continuous models based upon US EPA BMDS software can be specified: `hill`, `exp-3`, `exp-5`, `power`, `polynomial`. See R documentation ?ToxicR::single_continuous_fit for more details.

### Data Type
For both functions, dose-response data can be provided for individual subjects, or as a summary across dose groups. It is preferable to provide information on individual subjects however, in the case where this information is not available, summary data may be used.

Ex. Individual data
```{r}
# summarise mutation frequencies by sample
mf_data <- calculate_mut_freq(
            mutation_data = mut_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            summary = TRUE,
            retain_metadata_cols = "dose"
            )
# Fit a model-averaged BMD to the mutation frequency data.
bmd <- bmd_ma(mf_data,
              data_type = "individual"
              dose_col = "dose", 
              response_cols = c("sample_MF_min", "sample_MF_max"),
              bmr_type = "rel",
              bmr = 0.5,
              ...)
```

Ex. Summary data. `mf_data` should be a data frame containing the mean mutation frequency per dose. For each dose group, you must also provide the standard deviation and the sample size. Indicate the names of these four columns (mean response, dose, standard deviation, and sample size) using the `response_cols`, `dose_col`, `sd_col`, and `n_col` parameters respectively. 
```{r}
mf_data <- data.frame(
  dose = c("control", "D1", "D2", "D3"),
  mean_mf = c(1.3E-07, 3.3E-07, 6.8E-07, 1.0E-06),
  standard_deviation = c(1.0E-08, 1.6E-08, 2.3E-08, 2.8E-08),
  sample_size = c(6, 6, 6, 6))

bmd <- bmd_ma(mf_data,
              data_type = "summary"
              dose_col = "dose", 
              response_cols = "mean_mf",
              sd_col = "standard_deviation",
              n_col = "sample_size",
              bmr_type = "rel",
              bmr = 0.5,
              ...)
```

### Output
The BMD is reported alonside its upper and lower confidence intervals; the BMDU and BMDL. The BMDL is typically used to derive human health-based guidance values. The functions will also output several plots to visualise the results.

## Mutation Spectra Analysis
The mutation spectra is the pattern of mutation subtypes within a sample or group. The mutation spectra can inform on the mechanisms involved in mutagenesis.

### Comparison of Mutation Spectra Between Groups
We can compare the mutation spectra between experimental groups using the `spectra_comparison` function. This function will compare the proportion of mutation subtypes at any resolution between specified groups using a modified contingency table approach (Piegorsch and Bailer, 1994). 

This approach is applied to the mutation counts for each mutation subtype in a given group. The contingency table is represented as $R * T$ where R is the number of subtypes involved in the analysis, and T is the number of groups. The `spectra_comparison` function performs comparisons between T = 2 specified groups. The statistical hypothesis of homogeneity is that the proportion (count/group total) of each mutation subtype equals that of the other group. To test the significance of the homogeneity hypothesis, the $G^{2}$ likelihood ratio statistic is used: 

$$G^{2} = 2\  \sum_{i=1}^{R}\  \sum_{j=1}^{T}\  Y_{ij}\  log(\frac{Y_{ij}}{E_{ij}})$$

$Y_{ij}$ represents the mutation counts and $E_{ij}$ are the *expected* counts under the null hypothesis. The $G^{2}$ statistic possesses approximately a $\chi^{2}$ distribution in large sample sizes under the null hypothesis of no spectral differences. Thus, as the column totals become large, $G^{2}$ may be referred to a $\chi^{2}$ distribution with $(R -  1)(T - 1)$ degrees of freedom. It is important to note that the $G^{2}$ statistic may exhibit high false positive rates in small sample sizes when referred to a $\chi^{2}$ distribution. In such cases, we instead switch to an F-distribution. This has the effect of reducing the rate at which $G^{2}$ rejects each null hypothesis, providing greater stability in terms of false positive error rates. Thus when N/(R-1) < 20, where N  is the total mutation counts across both groups, the function will use a F-distribution, otherwise it will use a $\chi^{2}$-distribution.

This comparison assumes independance among the observations. Each tabled observation represents a sum of independent contributions to the total mutant count. We assume independance is valid for mutants derived from a  mixed population, however, mutants that are derived clonally from a single progenitor cell would violate this assumption. As such, it is recommended to use the **MFmin method** of mutation counting for spectral analyses  to ensure that all mutation counts are independant. In those cases where the independence may be invalid, and where additional, extra-multinomial sources of variability are present, more complex, hierarchical statistical models are required. This is currently outside the scope of this package.

The `spectra_comparison` function takes the imported mutation data. It will use the `calculate_mut_freq` function to sum each of the mutation subtypes across specified groups. Use the `subtype_resolution` and the `variant_types` parameters to specify the mutation subtypes that you wish to include in the analysis. Comparisons between groups are made based on an inputted contrasts table. The contrasts table will consist of two columns, each specifying a group to be contrasted against the other. 

*Ex. Consider a study in which we are studying the effect of a mutagenic chemical on the mutation spectra. Our samples were exposed to three doses of a mutagenic chemical (D1, D2, D3), or to the vehicle control (D0). We will compare the simple snv subtypes, alongside non-snv variants, of each of the three chemical dose groups to the control. In this way we can investigate if exposure to the mutagenic chemical leads to differences in the mutation spectrum. The function will output the $G^{2}$ statistic and p-value for each of the three comparisons listed in the `constrasts_table`. P-values are adjusted for multiple comparison using the Sidak method.*
```{r}
contrasts_table <- data.frame(col1 = c(D1, D2, D3),
                              col2 = c(D0, D0, D0))
simple_spectra <- spectra_comparison(mutation_data,
                                     cols_to_group = "dose",
                                     subtype_resolution = "base_6",
                                     variant_types =  c("snv",
                                                        "deletion",
                                                        "insertion",
                                                        "complex",
                                                        "mnv",
                                                        "symbolic"),
                                     mf_type = "min",
                                     contrasts = contrasts_table)
```

### Mutational Signatures Analysis
Mutational processes generate characteristic patterns of mutations, known as a mutational signatures. Distinct mutational signatures have been extracted from various cancer types and normal somatic tissues and deposited in the  Catalogue of Somatic Mutations in Cancer, or [COSMIC database](https://cancer.sanger.ac.uk/signatures/). These include signatures of single base substitutions (SBSs), doublet base substitutions (DBSs), small insertions and deletions (IDs) and copy number alterations (CNs). It is possible to assign mutational signatures to individual samples or groups using the `signature_fitting` function. This analysis provides the opportunity to identify the mutational processes involved in somatic mutagenesis within specific samples/groups. In the long run, it can provide evidence supporting the contribution of environmental mutagens to the mutation spectrum observed in human cancers. 

The `signature_fitting` function utilizes the [SigProfiler](https://github.com/AlexandrovLab) suite of tools developped by the Alexandrov lab. This function will create a virtual environment using reticulate to run python, which is required for the SigProfiler tools. It will also install several python dependencies using a conda virtual environment on first use, as well as the FASTA files for all chromosomes for your specified reference genome. As a result ~3Gb of storage must be available for the downloads of each genome. 

Somatic mutations in their 96-base trinucleotide context are summed across samples or experimental groups of interest to create a mutation count matrix. Mutational signatures are then assigned to each sample/group using refitting methods (Díaz-Gay et al., 2023). Signature refitting (also known as signature fitting or signature assignment) quantifies the contribution of a set of signatures to the mutational profile of a sample/group. The process is a numerical optimization approach that finds the combination of mutational signatures that most closely reconstructs the mutation count matrix. To quantify the number of mutations imprinted by each signature, the tool uses a custom implementation of the forward stagewise algorithm and it applies nonnegative least squares, based on the Lawson-Hanson method. 

Currently, `signature_fitting` offers fitting of COSMIC version 3.3 SBS signatures to the SBS96 matrix of any sample/group. For advanced use, including using a custom set of reference signatures, or fitting the DBS, ID, or CN signatures, it is suggested to use the SigProfiler python tools directly in python as described in their respective documentation [here](https://github.com/AlexandrovLab).

`signature_fitting` requires the installation of reticulate and [SigProfilerMatricGeneratorR](https://github.com/AlexandrovLab/SigProfilerMatrixGeneratorR) as well as a version of python 3.8 or newer.
```{r}
# Install reticulate
install.packages("reticulate")

# Install python
reticulate::install_python()

# Install SigProfilerMatrixGeneratorR from github using devtools.
library(devtools)
install_github("AlexandrovLab/SigProfilerMatrixGeneratorR")

```

`signature_fitting` will take the imported `mutation_data` and create mutational matrices for the somatic SNV mutations. Mutations are summed across levels of the `group` parameter. This can be set to individual samples or to an experimental group. 

The `project_genome` will be referenced for the creation of the mutational matrices.

The virtual environment can be specified with the `env_name` parameter. If no such environmnent exists, then the function will create one in which to store the dependencies and run the signature refitting. Specify your version of python using the `python_version` parameter (must be 3.8 or higher).

```{r}
signature_fitting(mutation_data,
                  project_name = "Example",
                  project_genome = "GRCh38", # human
                  env_name = "MutSeqR_example", # virtual environment
                  group = "sample",
                  output_path = NULL, # file path for output directory
                  python_version = "3.11" # your python version
                  )
```

#### Mutation Signatures Output
Results from the `signature_fitting` will be stored in an output folder. A filepath to a specific output directory can be designated using the `output_path` parameter.  If null, the output will be stored within your working directory. Results will be organized into subfolders based on the `group` parameter. The output structure is divided into three folders: input, output, and logs. 

The input folder contains copies of the `mutation_data`, following processing steps to make it compatible with the SigProfiler tools. It consists of a list of all the snv variants in each `group` alongside their genomic positions. This data serves as input for matrix generation.

The log folder contains the error and the log files for SigProfilerMatrixGeneration.

The output folder contains the results from matrix generation and signature refitting desribed in detail below. 

##### Matrix Generation
`signature_fitting` uses [SigProfilerMatrixGenerator](https://osf.io/s93d5/wiki/home/) to create mutational matrices (Bergstrom et al., 2019). Mutation matrices are created for single-base substitutions (SBS) and doublet-base substitutions (DBS), including matrices with extended sequence context and transcriptional strand bias. SBS and DBS matrices are stored in their respective folders in the output directory. **Only the SBS96 matrix is used for refitting**. 

**Single-Base Substitution Matrices:**
A single base substitution (SBS) is a mutation in which a single DNA base-pair is substituted with another single DNA base-pair. The most basic classfication catalogues SBSs into six distinct categories; C:G>A:T, C:G>G:C, C:G>T:A, T:A>A:T, T:A>C:G, T:A>G:C. However, it is common practice to report these mutations using the pyrimidine base of the Watson-Crick base-pair: C>A, C>G, C>T, T>A, T>C, T>G. These 6 mutation types can be further extended to include sequence context and strand information. The SBS-96, which includes the 6 mutations alongside their flanking nucleotides, is particularly useful for analysis of sequencing data. This classification is both simple enough to allow for visual inspection of mutational patterns and yet sufficiently complicated for seperating different sources of the same type of an SBS. The context can be further extended to include the flanking dinucleotides giving 1536 possible mutation classifications. However, use of the SBS-1536 requires a large number of somatic mutations, such as from the whole-sequencing of cancer samples with high mutational burden (> 2 mutations/Mb).

Generated matrices are described below. *Matrices are stored as `.all` files which can be viewed in a text-editor like notepad.*
| File                    |                                                                                                       |
|-------------------------|-------------------------------------------------------------------------------------------------------|
|*project_name.SBS6.all*  | 6-base. The 6 pyrimidine single-nucleotide variants. C > {A, G, or T} and T > {A, G, or C} = 6        |
|*project_name.SBS18.all* | 6-base. The 6 pyrimidine single-nucleotide variants within 3 transcriptional bias categories: Untranscribed (U), Transcribed (T), Non-Transcribed Region (N).|
|*project_name.SBS24.all* | 6-base. The 6 pyrimidine single-nucleotide variants within 4 transcriptional bias categories: Untranscribed (U), Transcribed (T), Bidirectional (B), Non-Transcribed Region (N).|
|**project_name.SBS96.all**  | 96-base. The 6 pyrimidine single-nucleotide variants alongside their flanking nucleotides (4 x 4 = 16 combinations). *Ex. A[C>G]T*|
|*project_name.SBS288.all*  |96-base. The 96-base single-nucleotide variants within 3 transcriptional bias categories (U, T, N).  |
|*project_name.SBS384.all* |96-base. The 96-base single-nucleotide variants within 4 transcriptional bias categories (U, T, N, B).|
|*project_name.SBS1536.all* | 1536-base. The 6 pyrimidine single-nucleotide variants alongside their flanking dinucleotides (16 x 16 = 256 combinations). *Ex. AA[C>G]TT*|
|*project_name.SBS4608.all* | 1536-base. The 1536-base single-nucleotide variants within 3 transcriptional bias categories (U, T, N).|
|*project_name.SBS6144.all* | 1536-base. The 1536-base single-nucleotide variants within 4 transcriptional bias categories (U, T, N, B).|

DBS Matrices
A doublet-base substitution (DBS) is a somatic mutation in which a set of two adjacent DNA base-pairs is *simultaneously* substituted with another set of two adjacent DNA base-pairs. **We do not recommend using the DBS matrices generated using `signature_fitting` for further analysis.** The `signature_fitting` function that we provide is designed to handle only the SBS mutations. All true multi-nucleotide variants, including doublets, are filtered out of the `mutation_data` prior to MatrixGeneration. However, the tool will still attempt to identify DBSs and will occasionally find two independent SBSs occuring next to each other simply by chance. If you wish to use DBS mutations in your signature analysis, please refer directly to the SigProfiler tools.

Other Usages
SigProfilerMatrix Generator also supports small indels, structural variants, and copy-number variants. `signature_fitting` will not generate these matrices. If you wish to utilise these features, please refer directly to the SigProfiler tools.

Barplots of the mutation matrices for all individual samples/groups can be found in the Plots folder. The number of mutations are plotted for each individual sample/group at the various subtype resolutions
 
##### Transcription Strand Bias (TSB)
SBS mutations will be tested for [transcription strand bias](https://osf.io/s93d5/wiki/5.%20Output%20-%20TSB/). These results will be stored in the `TSB` folder.

It is not possible to distinguish which of the two DNA strands a mutation originated on. However, one expects that mutations from the same type will be equally distributed across the two DNA strands. In other words, we expect mutations to occur at the same proportion as their reverse complement. For example, given a mutational process that causes purely C>G:T:A mutations, and a long repetitive sequence 5'- CGCGCGCGCGCGCGCGCGCGCG-3' on the reference genome, one would expect to see an equal number of C>T and G>A mutations. However, in many cases, an asymteric number of mutations are observed due to one of the strands having a higher propensity for being either damaged or repaired. A common example of this is transcription strand bias where the transcribed strand is subjected to higher rates of DNA repair as part of the transcriptional process compared to the untranscribed strand.

SigProfilerMatrixGenerator evaluates the transcriptional strand classificaiton of mutations within well-annotated protein coding genes of a reference genome. Mutations that occur outside of coding regions are classified as Non-transcribed (N). Mutations that occur within coding regions are classified as one of; transcribed (T), un-transcribed (U), bi-directional (B), or unknown. In order to classify mutations within coding regions, mutations are oriented based on the reference strand and their pyrimidine context. 

For example, consider that our reference sequence contains the coding sequence of a gene i.e. it is the coding/UN-transcribed DNA strand. A T>C mutation called on this reference sequence would be referred to as an untranscribed T>C (**U:T>C**). However, if instead an A>G mutation is called on this reference sequence, it would be referred to as a transcribed T>C (**T:T>C**). Purine-based mutations called from the reference sequence are converted to their pyrimidine context and this includes swapping to the complementary DNA strand. In this case, the reverse  complement of untranscribed A>G is transcribed T>C.

In rare cases, both strands of a genomic region code for a gene. Such mutations are annotated as bidirectional based on their pyrimidine context. For example, both T:A>C:G and A:T>G:C mutations in regions of bidirectional transcription will be annotated as bidirectional T>C (**B:T>C**) mutations. 

All SBS mutations will be classified within the four transcriptional bias categories:
* Transcribed (T)
    + The variant is on the transcribed (template) strand.
* Untranscribed (U)
    + The variant is on the untranscribed (coding) strand.
* Bidirectional (B)
    + The variant is on both strands and is transcribed either way.
* Nontranscribed (N)
    + The variant is in a non-coding region and is untranslated. 

The tool will then perform a transcription strand bias test which compares the number of transcribed and untranscribed mutations for each mutation type. For example, it will compare the number of transcribed T>C to untranscribed T>C mutations. Should there be a significant difference, it would indicate that T:A>C:G mutations are occuring at a higher rate on one of the strands compared to the other. Transcription strand bias tests will be included for the 6-base, 96-base and 288-base SBS mutation contexts. The output files contain the following information:
* the `group`
* the mutation type
* the enrichment value (# Transcribed / # untranscribed)
* the p-value, corrected for multiple comparisons using the false discrover rate method
* the false discovery rate q-value

Files include:
* *strandBiasTes_24.txt*: stats of the SBS6 variants
* *strandBiasTes_384.txt*: stats of the SBS96 variants
* *strandBiasTes_6144.txt*: stats of the SBS288 variants
* *significantResults_strandBiasTest.txt*: returns significant results from the three files above.

##### Signature Refitting Results

Results from the signature refitting perfomed by [SigProfilerAssignment](https://osf.io/mz79v/wiki/home/) will be stored within in the `Assignment_Solution` folder. `Assignment_Solution` consists of 3 subdirectories;  `Activities`, `Signatures`, and `Solution_Stats`. 

* Activities
    + *Assignment_Solution_Activities.txt*: the activity matrix for the selected signatures. The first column lists all of the samples/groups. All following columns list the calculated activity value for the respective signatures. The number of columns is the number of signatures identified. Signature activities correspond to the specific numbers of mutations from the original catalog caused by a particular mutational process.
    + *Assignment_Solution_Activity_Plots.pdf*: stacked barplot showing the number of mutations in each signature on the y-axis and the sample name on the x-axis.
    + *Assignment_Solution_TMB_plot.pdf*: tumor mutational burden plot. The y-axis is the somatic mutations per megabase and the x-axis is the number of samples/groups plotted over the total number of samples/groups included. The column names are the mutational signatures and the plot is ordered by the median somatic mutations per megabase.
    + *Decomposed_Mutation_Probabilities.txt*: the probabilities of each of the 96 mutation types in each sample/group. The probabilities refer to the probability of each  mutation type being caused by a specific signature. The first column lists all the samples/groups, the second column lists all the mutation types, and the following columns list the calculated probability value for the respective signatures.
    + *SampleReconstruction*: 
* Signatures
    + *Assignment_Solution_Signatures.txt*: the distribution of mutation types in the input mutational signatures. The first column lists all 96 of the mutation types. The following columns are the signatures.
    + *SBS_96_plots_Assignment_Solution.pdf*:  barplots for each signature identified that depicts the proportion of the mutation types for that signature. The top right corner also lists the total number of mutations and the percentage of total mutations assigned to the mutational signature.
* Solution_Stats
    + *Assignment_Solution_Samples_Stats.txt*: the accuracy metrics for the reconstruction. statistics for each sample including the total number of mutations, cosine similarity, L1 norm (calculated as the sum of the absolute values of the vector), L1 norm percentage, L2 norm (calculated as the square root of the sum of the squared vector values), and L2 norm percentage, along with the Kullback-Leibler divergence.
    + *Assignment_Solution_Signature_Assignment_log.txt*:  the events that occur when known signatures are assigned to an input sample. The information includes the L2 error and cosine similarity between the reconstructed and original sample within different composition steps.

# References
Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG, Alexandrov LB. SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. BMC Genomics. 2019 Aug 30;20(1):685. doi: 10.1186/s12864-019-6041-2. PMID: 31470794; PMCID: PMC6717374.

Díaz-Gay M, Vangara R, Barnes M, Wang X, Islam SMA, Vermes I, Narasimman NB, Yang T, Jiang Z, Moody S, Senkin S, Brennan P, Stratton MR, Alexandrov LB. Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment. bioRxiv [Preprint]. 2023 Jul 11:2023.07.10.548264. doi: 10.1101/2023.07.10.548264. Update in: Bioinformatics. 2023 Dec 1;39(12): PMID: 37502962; PMCID: PMC10369904.

Piegorsch WW, Bailer AJ. Statistical approaches for analyzing mutational spectra: some recommendations for categorical data. Genetics. 1994 Jan;136(1):403-16. doi: 10.1093/genetics/136.1.403. PMID: 8138174; PMCID: PMC1205789.

White PA, Long AS, Johnson GE. Quantitative Interpretation of Genetic Toxicity Dose-Response Data for Risk Assessment and Regulatory Decision-Making: Current Status and Emerging Priorities. Environ Mol Mutagen. 2020 Jan;61(1):66-83. doi: 10.1002/em.22351. Epub 2019 Dec 19. PMID: 31794061.







