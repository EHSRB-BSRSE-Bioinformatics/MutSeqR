<!-- badges: start -->
  [![R-CMD-check](https://github.com/EHSRB-BSRSE-Bioinformatics/duplex-sequencing/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EHSRB-BSRSE-Bioinformatics/duplex-sequencing/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->
  
# Change Report:
Major changes on 2025-03-24
- *filter_mut()* function added to workflow. This function filters the mutation_data: germline identification via vaf_cutoff, depth correction, and filtering variants based on regions have all been moved from the import functions to filter_mut(). calculate_mf(), plot_bubbles(), and signature_fitting() filter out variants using the filter_mut column instead of the is_germline column.
- calculate_mut_freq() is renamed to calculate_mf()
- total_depth and depth are no longer required for import
- calculate_mf() no longer requires depth. Users may choose to 1) calculate depth from mutation data, 2) supply precalculated depth values in a separate table, 3) No depth, mf is not calculated, only mutation sums.
- plot_spectra, plot_trinucleotide, and spectra_comparison are now supplied with mf_data instead of the mutation data.
- Example data has been added: Currently 44Mb

  
For full details on function utility, see below.

# MutSeqR: Error-corrected Next-Generation Sequencing (ECS) Analysis For Mutagenicity Assessment

## What is ECS?

Error-corrected next-generation sequencing (ECS) uses various methods to combine multiple independent raw sequence reads derived from an original starting molecule, thereby subtracting out artifacts introduced during sequencing or library preparation. This results in a highly accurate representation of the original molecule. ECS is particularly useful for detecting rare somatic mutations (or mutations induced in germ cells), such as those that arise from mutagen exposure or other sources of DNA damage. ECS is a powerful tool for assessing the mutagenicity of chemicals, drugs, or other agents, and can be used to identify the mutational signatures of these agents. ECS can also be used to detect rare mutations in cancer or other diseases, and to track the clonal evolution of these diseases over time.

For more background on how ECS works and its context in regulatory toxicology testing and genetic toxicology, see the following articles:
- [Menon and Brash, 2023](10.1016/j.mrrev.2023.108471)
- [Marchetti et al., 2023a](https://doi.org/10.1038/d41573-023-00014-y)
- [Marchetti et al., 2023b](https://doi.org/10.1016/j.mrrev.2023.108466)
- [Kennedy et al., 2014](https://doi.org/10.1038/nprot.2014.170)

This R package is meant to facilitate the import, cleaning, and analysis of ECS data, beginning with a table of variant calls or a variant call file (VCF). The package is designed to be flexible and enable users to perform common statistical analyses and visualisations.

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

## Data import

The main goal of MutSeqR is to generate summary statistics, visualizations, exploratory analyses, and other post-processing tasks such as mutational signature analysis or generalized linear modeling. Mutation data should be supplied as a table of variants with their genomic positions. Mutation data can be imported as either VCF files or as tabular data using the functions `import_vcf_data` and `import_mut_data`, respectively. It is reccomended that files include a record for every sequenced position, regardless of whether a variant was called or not, along with the `total_depth` for each record. This enables site-specific depth calculations that are required for the calculation of mutation subtype frequencies ad other site-specific frequemcies. The data set can be pared down later to include only mutations of interest (SNVs, indels, SVs, or any combination). 

Required columns for mutation data import:
| **Column** | **VCF Specification** | **Definition** |
|------------|-----------------------|----------------|
| contig | CHROM | The name of the reference sequence. |
| start | POS | The start position of the feature. 0-based coordinates are accepted but will be changed to 1-based during import. |
| end | INFO(END) | The half-open end position of the feature in contig. |
| ref | REF | The reference allele at this position. |
| alt | ALT | The left-aligned, normalized, alternate allele at this position. |
| sample | INFO(sample) or Genotype Header | A unique identifier for the sample library. For VCF files, this field may be provided in either the INFO field, or as the header to the GENOTYPE field. |
| *SUGGESTED FIELDS* | | |
| alt_depth | VD | The read depth supporting the alternate allele. If not included, the function will assume an alt_depth of 1 at variant sites. |
| total_depth or depth | AD or DP | The total read depth at this position. This column can be “total_depth” which excludes N-calls, or “depth”, which includes N-calls, if “total_depth” is not available. For VCF files, the total_depth is calculated as the sum of AD. DP is equivalent to "depth". |

VCF files should follow the VCF specification (version 4.5; Danecek et al. 2011). VCF files may be bg/g-zipped. Multiple sample VCF files are not supported. Multiple alt alleles called for a single position should be represented as sseparate rows in the data. All extra columns, INFO fields, and FORMAT fields will be retained upon import.

Upon import, records are categorized within the `variation_type` column based on their REF and ALT. Categories are listed below.
| variation_type | Definition |
|------------------|------------|
| no_variant | No variation, the null-case. |
| snv | Single nucleotide variant. |
| mnv | Multiple nucleotide variant. |
| insertion | Insertion. |
| deletion | Deletion. |
| complex | REF and ALT are of different lengths and nucleotide compositions. |
| symbolic | Structural variant |
| ambiguous | ALT contains IUPAC ambiguity codes. |
| uncategorized | The record does not fall into any of the preceding categories. |

Additional columns are created to further characterise variants.

| Column Name        | Definition                                        |
|--------------------|---------------------------------------------------|
|  short_ref | The reference base at the start position. |
| normalized_ref | The short_ref in C/T (pyrimidine) notation for this position. Ex. `A` -> `T`, `G` -> `C` |
| context | The trinucleotide context at this position. Consists of the reference base and the two flanking bases. Sequences are retrieved from the appropriate BS genome. Ex. `TAC` |
| normalized_context | The trinucleotide context in C/T (pyrimidine) notation for this position (Ex. `TAG` -> `CTA`) |
| variation_type | The type of variant (no_variant, snv, mnv, insertion, deletion, complex, sv, ambiguous, uncategorized) |
| subtype | The substitution type of the snv variant (12-base spectrum; Ex. `A>C`) |
| normalized_subtype | The snv subtype in C/T (pyrimidine) notation (6-base spectrum; Ex. `A>C` -> `T>G`) |
| context_with_mutation | The snv subtype including the two flanking nucleotides (192-base spectrum; Ex. `T[A>C]G`) |
| normalized_context_with_mutation | The snv subtype in C/T (pyrimidine) notation including the two flanking nucleotides (96-base spectrum; Ex. `T[A>C]G` -> `C[T>G]A`) |
| nchar_ref | The length (in bp) of the reference allele. |
| nchar_alt | The length (in bp) of the alternate allele. |
| varlen | The length (in bp) of the variant. |
| ref_depth | The depth of the reference allele. Calculated as `total_depth` - `alt_depth`, if applicable. |
| vaf | The variant allele fraction. Calculated as `alt_depth`/`total_depth` |
| gc_content | % GC of the trinucleotide context at this position. |
| is_known | A logical value indicating if the record is a known variant; i.e. ID field is not NULL. |
| row_has_duplicate | A logical value that  flags rows whose position is the same as that of at least one other row for the same sample. |

### General Usage:  `import_vcf_data` & `import_mut_data`
Indicate the file path to your mutation data file(s) using the `vcf_file`/`mut_file` parameter. This can be either a single file or a directory containing multiple files. Should you provide a directory, then all files within will be bound into a single data frame.

An important component of importing your data for proper use is to assign each mutation to a biological sample, and also make sure that some additional information about each sample is present (ex., a chemical treatment, a dose, etc.). This is done by providing `sample_data`. This parameter can take a data frame, or it can read in a file if provided with a filepath. If using a filepath, specify the proper delimiter using the `sd_sep` parameter. Sample metadata will be joined with the mutation data using the "sample" column to capture that information and associate it with each mutation.

Specify the appropriate BS genome with which to populate the context column by supplying the species, genome, and masked_BS_genome parameters. The function will browse [BSgenome::available.genomes](https://www.rdocumentation.org/packages/BSgenome/versions/1.40.1/topics/available.genomes) for the appropriate reference genome and install the corresponding package. Context information will be extracted from the installed BSgenome object. BSgenome offers genomes with masked sequences. If you wish to use the masked version of the genome, set `masked_BS_genome` to `TRUE`.

**Example Data**

*We provide an example data set taken from Leblanc et al., 2022. This data consists of 24 mouse bone marrow samples sequenced with Duplex Sequencing using Twinstrand's Mouse Mutagenesis Panel of twenty 2.4kb targeted genomic loci. Mice were exposed to three doses of benzo[a]pyrene (BaP) alongside vehicle controls, n = 6.*

*Example 1.1. Import the example .vcf.bgz file. Provided is the genomic vcf.gz file for sample dna00996.1. It is comprised of a record for all 48K positions sequenced for the Mouse Mutagenesis Panel with the alt_depth and the tota_depth values for each record.*
```{r}
example_file <- system.file("extdata", "example_import_vcf_data_cleaned.vcf.bgz", package = "MutSeqR")
sample_metadata <- data.frame(sample = "dna00996.1",
                          dose = "50",
                          dose_group = "High")
# Import the data
imported_example_data <- import_vcf_data(vcf_file = example_file,
                                         sample_data = sample_metadata,
                                         genome = "mm10",
                                         species = "mouse",
                                         masked_BS_genome = FALSE)

```

*Example 1.2. Import the example tabular data. This is the equivalent file to the example vcf file. It is stored as an .rds file. We will load the data frame and supply it the `import_mut_data`. The mut_file parameter can accept file paths or data frames as input.*
```{r}
example_file <- system.file("extdata", "example_import_mut_data.rds", package = "MutSeqR")
example_data <- readRDS(example_file)
sample_metadata <- data.frame(sample = "dna00996.1",
                              dose = "50",
                              dose_group = "High")
# Import the data
imported_example_data <- import_mut_data(mut_file = example_data,
                                         sample_data = sample_metadata,
                                         genome = "mm10",
                                         species = "mouse",
                                         masked_BS_genome = FALSE,
                                         is_0_based_mut = TRUE)
# is_0_based_mut indicates that the genomic coordinates are 0-based. These will
# be changed to 1-based upon import. 
```

#### Variants within target regions
Similar to sample metadata, you may supply a file containing the metadata of genomic regions to the `regions` & `custom_regions` parameters. Region metadata will be joined with mutation data by checking for overlap between the target region ranges and the position of the record.
The `regions` parameter can be set to one of TwinStrand's DuplexSeq™ Mutagenesis Panels; *TSpanel_mouse*, *TSpanel_human*, or *TSpanel_rat*. If you are using an alternative panel then you may set the `regions` parameter to  "custom" and  you will add your target regions' metadata using a `custom_regions`. You may supply your custom_regions file as either a data frame or a file path, which will be read in. Required columns are `contig`, `start`, and `end`. Use parameters to indicate your file's delimiter and whether the region coordinates are 0-based or 1-based. Mutation data and region coordinates will be converted to 1-based. If you do not wish to specify regions, then set the `regions` parameter to *none*.

*Example 1.3. Add the metadata for TwinStrand's Mouse Mutagenesis panel to our example vcf file.*
```{r}
imported_example_data <- import_vcf_data(vcf_file = example_file,
                                         sample_data = sample_metadata,
                                         genome = "mm10",
                                         species = "mouse",
                                         masked_BS_genome = FALSE,
                                         regions = "TSpanel_mouse")
```                                         
*To see an example of the region files, you can load the TSpanels:*
```{r}
region_example <- load_regions_file("TSpanel_mouse")
```

#### Column name synonyms
We recognize that column names may differ from those sepcified in the Required Columns table. Therefore, we have implemented some default column name synonyms. If your column name matches one of our listed synonyms, it will automatically be changed to match our set values. For example, your `contig` column may be named `chr` or `chromosome`. After importing your data, this synonymous column name will be changed to `contig`. Column names are case-insensitive. A list of column name synonyms are listed below. If your data contains a column that is synonymous to one of the required columns, but the name is not included in our synonyms list, your column name may be substituted using the `custom_column_names` parameter. Provide this parameter with a list of names to specify the meaning of column headers.
```{r}
library(MutSeqR)
# mut_data <- "file path to .mut file"
mutation_data <-
  import_mut_data(mut_file = mut_data,
                  custom_column_names = list(contig = "custom_contig_name",
                                             sample = "custom_sample_name"))
```
| Column | Synonyms |
|--------|----------|
| alt | alternate |
| alt_depth | alt_read_depth, var_depth, variant_depth, vd |
| context | sequence_context, trinucleotide_context, flanking_sequence |
| contig | chr, chromosome, seqnames |
| depth | dp |
| end | end_position |
| no_calls | n_calls, n_depth, no_depth |
| ref | reference, ref_value |
| sample | sample_name, sample_id |
| start | pos, position |
| total_depth | informative_somatic_depth |

#### Output
Mutation data can be output as either a data frame or a *GRanges* object for downstream analysis. Use the `output_granges` parameter to specify the output. *GRanges* may faciliate use in other packages and makes doing genomic based analyses on the ranges significantly easier. Most downstream analyses provided by MutSeqR will use a data frame.

## Variant Filtering
Following data import, the mutation data may be filtered using the `filter_mut()` function. This function will flag variants based on various parameters in the `filter_mut` column. Variants flagged as TRUE in the `filter_mut` column will be automatically excluded from downstream analyses such as `calculate_mf()`, `plot_bubbles`, and `signature_fitting`. When specified, this function may also remove records from the mutation data. Finally, the function offers other data cleaning options.

Flagging the variants in the `filter_mut` column will not remove them from the mutation data, however, they will be excluded from mutation counts in all downstream analyses. Filtered variants are retained in the data so that their `total_depth` values may still be used in frequency calculations.

By default, all filtering parameters are disabled. Users should be mindful of the filters that they use, ensuring first that they are applicable to their data.

### Correcting the Depth
There may be cases when the same genomic position is represented multiple times within the mutation data of a single sample. For instance, we require that multiple different alternate alleles for a single position be reported as seperate rows. Another common example of this that we've observed are deletion calls. Often, for each deletion called, a no_variant call at the same start position will also be present. For data that includes the `total_depth` of each record, we want to prevent double-counting the depth at these positions during downstream analyses. When set to TRUE, `correct_depth` parameter will correct the depth at these positions. For positions with the same sample, contig, and start values, the `total_depth` will be retained for only one row. All other rows in the group will have their `total_depth` set to 0. The import functions will automatically check for duplicated rows and return a warning advising to correct the depth should it find any in the mutation data. This is recommended for all users whose data contain `total_depth` values.

*NOTE: The case of deletions and no_variants: When identifying a deletion, some variant callers will re-align the sequences. As such a second total_depth value will be calculated, specifically for the deletion. The variant caller will then provide both the no_variant and the deletion call, the former with the initial depth, and the latter with the re-aligned depth value. Generally, when correcting the depth, the function will retain the total_depth of the first row in the group, and set the rest to 0. However, in the case of deletions, the function will prioritize retaining the re-aligned total_depth over the total_depth assigned to the no_variant. This prioritization can be disabled with the `correct_depth_to_indel` parameter.*

### Filtering Germline Variants
Germline variants may be identified and flagged for filtering by setting the `vaf_cutoff` parameter. The variant allele fraction (vaf) is the fraction of haploid genomes in the original sample that harbor a specific mutation at a specific base-pair coordinate of the reference genome. Specifically, is it calculated by dividing the number of variant reads by the total sequencing depth at a specific base pair coordinate (`alt_depth` / `total_depth`). In a typical diploid cell, a homozygous germline variant will appear on both alleles, in every cell. As such, we expect this variant to occur on every read giving us a vaf = 1. A heterozygous germline variant occurs on one of the two alleles in every cell, as such we expect this variant to occur on about half of the reads, giving a vaf = 0.5. Somatic variants occur in only a small portion of the cells, thus we expect them to appear in only a small percentage of the reads. Typical vaf values for somatic variants are less than 0.01 - 0.1. Setting the `vaf_cutoff` parameter will flag all variants that have a vaf **greater than this value as germline** within the `is_germline` column. It will also flag these variants in the `filter_mut` so as to exclude them from downstream analyses.

vaf germline filtering is only applicable to users whose data was sequenced to a sufficient depth. High coverage-low depth sequencing cannot be used to calculate the vaf, thus it is reccomended to filter germline mutations by contrasting against a database of known polymorphisms, or by using conventional whole-genome sequencing to identify germline variants for each sample. 

### Quality Assurance Filtering
The `filter_mut()` function offers some filtering options to ensure the quality of the mutation data.

`snv_in_germ_mnv` = TRUE will flag all snvs that overlap with germline mnvs. Germline mnvs are defined as having a vaf > vaf_cutoff. These snvs are often artifacts from variant calling. No-calls in reads supporting the germline mnv will create false minor-haplotypes from the original mnv that can appear as sub-clonal snvs, thus such variants are excluded from downstream analyses.

`rm_abnormal_vaf` = TRUE This parameter identifies rows with abnormal vaf values and **removes** them from the mutation data. Abnormal vaf values are defined as being between 0.05-0.45 or between 0.55-0.95. In a typical diploid organism, we expect variants to have a vaf ~0, 0.5, or 1, reflecting rare somatic variants, heterozygous variants, and homozygous variants respectively. Users should be aware of the ploidy of their model system when using this filter. Non-diploid organisms may exhibit different vafs.

`rm_filtered_mut_from_depth = TRUE` Variants flagged in the `filter_mut` column will have their `alt_depth` subtracted from the `total_depth`. When TRUE, this parameter treats flagged variants as No-calls. This does not apply to variants that were idenfied as germline variants.

### Custom Filtering
Users may use the `filter_mut()` function to flag or remove variants based on their own custom columns. Any record that contains the `custom_filter_val`value within the `custom_filter_col` column of the mutation data will be either flagged in the `filter_mut` column or, if specified by the `custom_filter_rm` parameter, removed from the mutation data.

### Filtering by Regions
Users may remove rows that are either within or outside of specified genomic regions. Provide the region ranges to the `regions` parameter. This may be provided as either a file path or a data frame. `regions` must contain `contig`, `start`, and `end`.  The function will check whether each record falls within the given regions. Users can define how this filter should be used with `regions_filter`. `region_filter = "remove_within"` will remove all rows whose positions overlaps with the provided regions. `region_filter = "keep_within"` will remove all rows whose positions are outside of the provided regions. By default, records that are > 1bp must start and end within the regions to count as being within the region. `allow_half_overlap = TRUE` allow records that only start or end within the regions but extend outside of them to be counted as being within the region.

### Retain your filtered rows
`return_filtered_rows = TRUE` The function will return both the filtered mutation data and the rows that were removed/flagged in a seperate data frame. The two dataframes will be returned inside of a list, with names `mutation_data` and `filtered_rows`. Default is FALSE.

### Example
*Example 2. The example file "example_mutation_data.rds is the output of import_mut_data() run on all 24 mouse libraries from (LeBlanc et al., 2022).*

*Filters used:*
- *Depth correction*
- *Filter germline variants: vaf < 0.01*
- *Filter snvs overlapping with germline variants and have their `alt_depth` removed from their `total_depth`.*
- *Remove records outside of the TwinStrand Mouse Mutagenesis Panel. Here we are using `load_regions_file() `to grab the TSpanel_mouse regions from the package files*
- *Filter variants that contain "EndRepairFillInArtifact" in the "filter" column. Their `alt_depth` will be removed from their `total_depth`.*
```{r}
# load the example data
example_file <- system.file("extdata",
                            "example_mutation_data.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Filter
filtered_example_mutation_data <- filter_mut(
  mutation_data = example_data,
  correct_depth = TRUE,
  vaf_cutoff = 0.01,
  regions = load_regions_file("TSpanel_mouse"),
  regions_filter = "keep_within",
  custom_filter_col = "filter",
  custom_filter_val = "EndRepairFillInArtifact",
  custom_filter_rm = FALSE,
  snv_in_germ_mnv = TRUE,
  rm_filtered_mut_from_depth = TRUE,
  return_filtered_rows = FALSE
)
```

## Calculating Mutation Frequencies
Mutation Frequency (MF) is calculated by dividing the sum of mutations by the sum of the `total_depth` for a given group (units mutations/bp). The function `calculate_mf()` sums mutation counts and `total_depth` across user-defined groupings within the mutation data to calculate the MF. Mutations can be summarised across samples, experimental groups, and mutation subtypes for later statistical analyses. Variants flagged as `TRUE` in the `filter_mut` column will be excluded from the mutation sums; however, the *total_depth* of these variants will be counted in the *total_depth* sum.

If the total_depth is not available, the function will sum the mutations across groups, without calculating the mutation frequencies.

### Mutation Counting Methods
 Mutations are counted based on two opposing assumptions (Dodge et al., 2023):

The **Minimum Independent Mutation Counting Method** (min) Each mutation is counted once, regardless of the number of reads that contain the non-reference allele. This method assumes that multiple instances of the same  mutation within a sample are the result of clonal expansion of a single mutational event. When summing mutations across groups using the Min method (`sum_min`), the `alt_depth` of each variant is set to 1. Ex. For 3 variants with `alt_depth` values of 1, 2, and 10, the `sum_min `= 3.

The **Maximum Independent Mutation Counting Method** (max) Multiple identical mutations at the same position within a sample are counted as independent mutation events. When summing mutations across groups using the Max method (`sum_max`), the `alt_depth` of each variant is summed unchanged. Ex. for 3 variants with `alt_depth` values of 1, 2, and 10, the `sum_max` = 13.

The Min and Max mutation counting methods undercount and overcount the mutations, respectively. We expect some recurrent mutations to be a result of clonal expansion. We also expect some recurrent mutations to arise independently of each other. As we expect the true number of independent mutations to be somewhere in the middle of the two counting methods, we calculate frequencies for both methods. However, the Min mutation counting method is generally recommended for statistical analyses because the Max method tends to increase the sample variance by a significant degree. 

### Grouping Mutations
Mutation counts and `total_depth` are summed across groups that can be designated using the `cols_to_group` parameter. This parameter can be set to one or more columns in the mutation data that represent experimental variables of interest. 

*Example 3.1. Calculate mutation sums and frequencies per sample. The file example_mutation_data_filtered.rds is the output of filter_mut() from Example 2*
```{r}
# load example data:
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

mf_data <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "none"
)
```

Alternatively, you can sum mutations by experimental groups such as `label` (the genomic target). Counts and frequencies will be returned for every level of the designated groups.

*Example 3.2. Calculate mutation sums and frequencies per sample and genomic target.*
```{r}
mf_data <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = c("sample", "label"),
  subtype_resolution = "none"
)
```

**calculate_mf() does not calculate the mean MF for any given group** If you want to calculate the mean MF for a given experimental variable, you may group by "sample" and retain the experimental variable in the summary table for averaging.

*Example 3.3. Calculate the mean MF per dose* 
```{r}
mf_data <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "none",
  retain_metadata_cols = "dose"
)
mean <- mf_data %>%
  dplyr::group_by(dose) %>%
  dplyr::summarise(mean_mf_min = mean(mf_min), SE = sd(mf_min)/sqrt(n()))
```
### Mutation Subtypes
Mutations can also be grouped by mutation subtype at varying degrees of resolution using the `subtype_resolution` parameter.
|Subtype resolutions| Definition | Example |
|-------------------|------------|---------|
| `type`            | The variation type | snv, mnv, insertion, deletion, complex, and symbolic variants |
| `base_6`          | The simple spectrum; snv subtypes in their pyrimidine context | C>A, C>G, C>T, T>A, T>C, T>G |
| `base_12`         | The snv subtypes | A>C, A>G, A>T, C>A, C>G, C>T, G>A, G>C, G>T, T>A, T>C, T>G |
| `base_96`         | The trinculeotide spectrum; the snv subtypes in their pyrimidine context alongside their two flanking nucleotides | A[C>T]A |
| `base_192`        | The snv subtypes reported alongside their two flanking nucleotides | A[G>T]A |

Mutations and `total_depth` will be summed across groups for each mutation subtype to calculate frequencies. For SNV subtypes, the `total_depth` is summed based on the sequence context in which the SNV subtype occurs (`subtype_depth`). In the simplest example, for the `base_6` SNV subtypes, the two possible reference bases are C or T; hence, the `total_depth` is summed separately for C:G positions and T:A positions. Thus, the MF for C>T mutations is calculated as the total number of C>T mutations divided by total_depth for C:G positions within the group: `sum` / `subtype_depth`. Non-snv mutation subtypes, such as mnvs, insertions, deletions, complex variants, and structural variants, will be calculated as their `sum` / `group_depth`, since they can occur in the context of any nucleotide.

Upon import of mutation data, columns are created that facilitate the grouping of SNV subtypes and their associated sequence context for the various resolutions. Below the columns associated with each subtype_resolution are defined:
|Subtype resolution | Subtype Column | Context Column |
|-------------------|----------------|----------------|
| `type` | `variation_type`| NA |
| `base_6` | `normalized_subtype` | `normalized_ref` |
| `base_12` | `subtype` | `short_ref` |
| `base_96` | `normalized_context_with_mutation` | `normalized_context` |
| `base_192` | `context_with_mutation` | `context` |

The function will also calculate the the proportion of mutations for each subtype, normalized to the `total_depth`:
$$P_s = \frac{\left(\frac{M_s}{D_s}\right)}{\sum_s \left(\frac{M_s}{D_s}\right)}$$
Where $P_s$ is the normalized mutation proportion for subtype $s$. $M_s$ is the group mutation sum for subtype $s$. $D_s$ is the group sum of the `subtype_depth` for subtype $s$.

If total_depth is not available for the mutation data, `calculate_mf()` will return the subtype mutation counts per group. It will also calculate subtype proportion, without normalizing to the total_depth:
$$P'_s = \frac{M_s}{M_{total}}$$
Where, $P'_s$ is the non-normalized mutation proportion of subtype $s$. $M_s$ is the group mutation sum for subtype $s$. $M_{total}$ is the total mutation sum for the group.


*Example 3.4. The following code will return the base_6 mutation spectra for all samples with mutation proportions normalized to depth.*
```{r}
mf_data <- calculate_mf(
            mutation_data = example_data,
            cols_to_group = "sample",
            subtype_resolution = "base_6"
            )
```

#### Selecting Variation Types
`calculate_mf()` can be used on a user-defined subset of `variation_type` values. The `variant_types` parameter can be set to a character string of `variation_type` values that the user wants to include in the mutation counts. When calculating group mutation sums, only variants of the specified `variation_type`s will be counted. The `total_depth` for records of excluded `variation_types` will still be included in the `group_depth` and the `subtype_depth`, if applicable.

By default the function will calculate summary values based on all mutation types.

*Example 3.5. Calculate global mutation frequencies per sample including only insertion and deletion mutations in the count.*
```{r}
mf_data <- calculate_mf(
            mutation_data = example_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            variant_types = c("insertion", "deletion")
            )
```

Users may also supply a list of variation_types to exclude to the `variant_types` parameter, so long as the value is preceeded by a "-".

*Example 3.6. Calculate mutation frequencies at the "type" subtype resolution, excluding ambiguous and uncategorized mutations.*
```{r}
mf_data <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "type",
  variant_types = c("-ambiguous", "-uncategorized")
)
```
*Example 3.7. Include only snv mutations at the base_96 resolution*
```{r}
mf_data <- calculate_mf(
            mutation_data = example_data,
            cols_to_group = "sample",
            subtype_resolution = "base_96",
            variant_types = "snv"
            )
```

### Precalculated Depth
For mutation data that does not include a `total_depth` value for each sequenced site, precalculated depths for the specified groups can be supplied separately to `calculate_mf()` through the `precalc_depth_data` parameter. Depth should be provided at the correct subtype resolution for each level of the specified grouping variable. This parameter will accept either a data frame or a file path to import the depth data. Required columns are 
1) The user-defined cols_to_group variable(s)
2) `group_depth`: the `total_depth` summed across the cols_to_group.
3) The context column for the specified `subtype_resolution`. Only applicable if using the SNV resolutions (base_6, base_12, base_96, base_192). Column names are listed in the table above.
4) `subtype_depth`: the `total_depth` summed across the cols_to_group for each context. Only applicable if using the SNV resolutions.

*Example 3.8. Use precalculated depth values to calculate the global per sample MF.*
```{r}
sample_depth <- data.frame(sample = unique(example_data$sample),
                           group_depth =c(565395266, 755574283, 639909215,
                                          675090988, 598104021, 611295330,
                                          648531765, 713240735, 669734626,
                                          684951248, 716913381, 692323218,
                                          297661400, 172863681, 672259724,
                                          740901132, 558051386, 733727643,
                                          703349287, 884821671, 743311822,
                                          799605045, 677693752, 701163532))
mf_data <- calculate_mf(mutation_data = example_data,
                        cols_to_group = "sample",
                        subtype_resolution = "none",
                        calculate_depth = FALSE,
                        precalc_depth_data = sample_depth)
```
*Example 3.9. Use precalculated depth values to calculate the base_6 per sample MF.*
```{r}
mf_data <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "base_6",
  calculate_depth = FALSE,
  precalc_depth = system.file("extdata",
                              "precalc_depth_base_6_example.text",
                              package = "MutSeqR")
)
# To view the example precalc depth file:
depth <- read.table(file = system.file("extdata", "precalc_depth_base_6_example.text", package = "MutSeqR"), header = TRUE)
View(depth)
```

### Summary Table
The function will output the resulting `mf_data` as a data frame with the MF and proportion calculated. If the `summary` parameter is set to `TRUE`, the data frame will be a summary table with the MF calculated for each group. If `summary` is set to `FALSE`, the MF will be appended to each row of the original `mutation_data`.

The summary table will include:
* `cols_to_group`: all columns used to group the data.
* Subtype column for the given resolution, if applicable.
* Context column for the given resolution, if applicable.
* `sum_min` & `sum_max`: the min/max mutation counts for the group/subtype.
* `group_depth`: the total_depth summed across the group.
* `subtype_depth`: the total_depth summed across the group, for the given context, if applicable.
* `MF_min` & `MF_max`: the min/max MF for the group/subtype.
* `proportion_min` & `proportion_max`: the min/max subtype proportion, if applicable.

Additional columns from the orginal mutation data can be retained using the `retain_metadata_cols` parameter. Retaining higher-order experimental groups may be useful for later statistical analyses or plotting. See above for example *calculating mean MF per dose*.


**NOTE** The summary table uses a pre-defined list of possible subtypes for each resolution. If a particular subtype within a given group is not recorded in the mutation data, the summary table will have no frame of reference for populating the `metadata_cols`. Thus, for subtypes that do not occur in the mutation data for a given group, the corresponding `metadata_col` will be NA.

## Plotting Mutation Frequency
You can visualize the results of `calculate_mf` using the `plot_mf()` and `plot_mean_mf`() functions. These functions offer variable aesthetics for easy visualization. The output is a ggplot object that can be modified using ggplot2.

*Example 3.10. Plot the Min and Max MF per sample, coloured and ordered by dose group.*
```{r}
# load example data:
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

mf_data <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "none",
  retain_metadata_cols = "dose_group"
)
# Define the order for dose groups
mf_data$dose_group <- factor(mf_data$dose_group,
                             levels = c("Control",
                                        "Low",
                                        "Medium",
                                        "High"))
plot <- plot_mf(mf_data = mf_data,
                group_col = "sample",
                plot_type = "bar",
                mf_type = "both",
                fill_col = "dose_group",
                group_order = "arranged",
                group_order_input = "dose_group")

```
![plot_mf](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot3.10.png)


Calculate and plot the mean MF for a user-defined group using `plot_mean_mf()`.

*Example 3.11. Plot the mean MF min per dose, including SEM and individual values coloured by dose.*
```{r}
plot_mean <- plot_mean_mf(mf_data = mf_data,
                          group_col = "dose_group",
                          mf_type = "min",
                          plot_type = "line",
                          fill_col = "dose_group",
                          plot_error_bars = TRUE,
                          plot_indiv_vals = TRUE,
                          add_labels = "none")
```
![plot_mean_mf](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot3.11.png)

## Generalized Linear Modelling
An important component of analysing mutagencity data is how MF changes based on experimental variables.

The `model_mf()` function can be used to quantitatively evaluate how MF changes based on experimental variables. The function models the proportion of reads that carry a mutation (i.e. the MF) across user-supplied fixed and/or random effects and interaction parameters. Depending on the supplied effects, `model_mf()` will automatically choose to fit either a generalized linear model (GLM) using the `glm()` function from the stats library or a generalized linear mixed model (GLMM) using the `glmer()` function from the lme4 library.

### Model Formula
- **Fixed effects**: The experimental variable/factor of interest.
- **Random effects**: Variable used to account for variability within a larger group. Covariate.
- **Interaction**: The product of two or more interacting fixed_effects. An interaction implies that the effect of one fixed effect changes depending on the levels of another fixed effect, indicating a non-additive reationship.
- **Response**: the function will model the MF, taking the mutation count (`sum`) and the `group_depth` as input.

The model formula is built as:

`cbind(sum, group_depth) ~ fixed_effect1 + fixed_effect2 + ... + (1|random_effect)`

If `test_interaction` is set to TRUE, the model will use the product of the fixed_effects:

`cbind(sum, group_depth) ~ fixed_effect1 * fixed_effect2 * ... + (1|random_effect)
`

### Family
The occurence of mutations is assumed to follow a binomial distribution as:
1) There is a finite number of sequenced bases.
2) A mutation at any given base is equally probable.
3) Mutations occur independently of other mutations.

To account for over-dispersion, the function will fit a GLM with a quasibinomial distribution.  If the dispersion parameter of the model is less than 1, the model will be refit using a binomial distribution. The GLMM will always use a binomial distribution.

The `model_mf` function will fit a generalized linear model to analyse the effect(s) of given fixed_effects(s) on MF and perform specified pairwise comparisons between levels of your factors. 

Additional arguments can be passed to the model to further customize it to your needs. Details on the arguments for the generalized linear model can be found at [stats::glm](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm) and for the general linear mixed model at [lme4::glmer](https://www.rdocumentation.org/packages/lme4/versions/1.1-35.3/topics/glmer). 

### Estimates and Comparisons
`model_mf()` provides estimates of the mean for each level of the fixed effects. Furthermore, pairwise comparisons can be performed based on a user-supplied table. Mean estimates and comparisons are conducted using the doBy R library (Halekoh and Højsgaard 2024). Estimates of the back-transformed standard errors are approximated using the delta method. The p-values are adjusted for multiple comparisons using the Sidak method. 

### Goodness of Fit
The `model_mf` function will output the Pearson residuals appended to the `mf_data`. Additionally, Pearson residuals will be plotted as a histogram and a QQ-plot to check for deviances from model assumptions. We assume that residuals will follow a normal distribution with a mean of 0. Thus, we expect the histogram to follow a bell curve and the QQ-plot to be plotted along the y=x line.

 ### General Usage: model_mf
Mutation data should first be summarised by sample using `calculate_mf()`. The `mf_data` should be output as a summary table. Be sure to retain the columns for experimental variables of interest using the `retain_metadata_cols` parameter.

You may specify factors and covariates for your model using the `fixed_effects` and `random_effects` parameters respectively. If more than one fixed_effect is supplied, then you may specify whether you wish to test the interaction between your fixed_effects using the `test_interaction` parameter. 

You must specify the columns in your `mf_data` that contain the mutation counts and the total_depth per sample using the `muts` and `total_counts` parameters respectively. 

To perform pairwise comparisons between levels of your fixed effects, supply a constrast table using the `contrasts` parameter. This can either be a data frame  or a file path to a text file that will be loaded into R. The table must consist of two columns, each containing levels within the fixed_effects. The level in the first column will be compared to the level in the second column for each row. You should also provide the reference level for each fixed effect using the `reference_level` parameter. If you specify multiple pairwise comparisons, then the p-values will automatically be corrected for multiple comparisons using the Sidak method. For multiple fixed effects, the user must include levels for all fixed_effects in each value of the contrasts table. Within each value, the levels of the different fixed_effects should be seperated by a colon.

### Output
The function will output a list of results.

* model_data: the supplied `mf_data` with added column for Pearson residuals.
* summary: the summary of the model.
* anova: the analysis of variance for models with two or more fixed_effects. See [car::Anova.](https://www.rdocumentation.org/packages/car/versions/1.0-9/topics/Anova)
* residuals_histogram: the Pearson residuals plotted as a histogram. This is used to check whether the variance is normally distributed. A symmetric bell-shaped histogram, evenly distributed around zero indicates that the normality assumption is likely to be true.
* residuals_qq_plot: the Pearson residuals plotted in a quantile-quantile plot. For a normal distribution, we expect points to roughly follow the y=x line.  
* point_estimates_matrix: the contrast matrix used to generate point-estimates for the fixed effects. 
* point_estimates: the point estimates for the fixed effects.
* pairwise_comparisons_matrix: the contrast matrix used to conduct the pairwise comparisons specified in the `contrasts`.
* pairwise_comparisons: the results of pairwise comparisons specified in the `contrasts`.
### Examples

*Example 4.1. Model the effect of dose on MF. Our example data consists of 24 mouse samples, exposed to 3 doses of BaP or a vehicle control. Dose Groups are : Control, Low, Medium, and High. We will determine if the MFmin of each BaP dose group is significantly increased from that of the Control group.*
```{r}
# load example data:
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)
# Calculate the MF per sample, retaining the dose_group in the summary table.
mf_data <- calculate_mf(
            mutation_data = example_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            retain_metadata_cols = "dose_group"
            )
# Create a contrasts table for pairwise comparisons
contrasts_table <- data.frame(col1 = c("Low", "Medium", "High"),
                              col2 = c("Control", "Control", "Control"))
# Run the model
model_by_dose <- model_mf(mf_data = mf_data,
                          fixed_effects = "dose_group",
                          muts = "sum_min",
                          total_count = "group_depth",
                          reference_level = "Control",
                          contrasts = contrasts_table
                          )
# View the results
model_by_dose$summary
model_by_dose$point_estimates
model_by_dose$pairwise_comparisons
```
![model_mf diagnostic plots](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot4.1.png)

*Example 4.2. Model the effects of dose and genomic locus on MF. Seqencing for the example data was done on a panel of 20 genomic targets. We will determine if the MF of each BaP dose group is significantly different from the Control individually for all 20 targets. In this model, dose group and target label will be our fixed effects. We include the interaction between the two fixed effects. Because sample will be a repeated measure, we will use it as a random effect.*
```{r}
# load example data:
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Summarise mutations by sample and by genomic target.
mf_data <- calculate_mf(
            mutation_data = example_data,
            cols_to_group = c("sample", "label"),
            subtype_resolution = "none",
            summary = TRUE,
            retain_metadata_cols = "dose_group"
            )

# Create a contrasts table for the pairwise comparisons.
combinations <- expand.grid(dose_group = unique(mf_data$dose_group),
                            label = unique(mf_data$label))
combinations <- combinations[combinations$dose_group != "Control", ]
combinations$col1 <- with(combinations, paste(dose_group, label, sep=":"))
combinations$col2 <- with(combinations, paste("Control", label, sep=":"))
contrasts2 <- combinations[, c("col1", "col2")]

# Run the model
# To improve convergence of the model, we will supply the control argument
# to the  glmer function
model_by_target <- model_mf(mf_data = mf_data,
  fixed_effects = c("dose_group", "label"),
  test_interaction = TRUE,
  random_effects = "sample",
  muts = "sum_min",
  total_count = "group_depth",
  contrasts = contrasts2,
  reference_level = c("Control", "chr1"),
  control = lme4::glmerControl(check.conv.grad = lme4::.makeCC("warning",
                                                               tol = 3e-3,
                                                               relTol = NULL))
  )
# View the results
model_by_target$summary
model_by_target$anova
model_by_target$point_estimates
model_by_target$pairwise_comparisons
```
![model_mf diagnostic plots](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot4.2.png)

### Visualize Model Results
Plot the results of `model_mf()` using `plot_model_mf()`. This function will create a bar or line plot of the point estimates. Users can include the estimated standard error as error bars with `plot_error_bars` and add significance labels based on the pairwise comparisons with `plot_signif`. The function can plot model results with up to two fixed effects. Users must specify which fixed effect should be represented on the x-axis with `x_effect`.

When adding significance labels, the function will generate a unique symbol for each unique level or combination of levels in column 2 of the contrasts table. If using two fixed effects, but your comparisons only contrast levels of one fixed effect, you may specify this using the `ref_effect` parameter. Supply it with the fixed effect that is being contrasted to reduce the number of unique symbols generated.

The other fixed effect will be represented by colour. The output is a ggplot object that can be modified with ggplot2.

*Example 4.1. model by dose*
```{r}
plot <- plot_model_mf(model_by_dose,
                      plot_type = "bar",
                      x_effect = "dose",
                      plot_error_bars = TRUE,
                      plot_signif = TRUE,
                      x_order = c("Control", "Low", "Medium", "High"),
                      x_label = "Dose Group",
                      y_label = "Estimated Mean Mutation Frequency (mutations/bp)")
```
![plot_model_mf Dose Model](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot4.1.2.png)

*Example 4.2. model by dose and genomic locus. In this example, we only made comparisons between dose groups. For each contrast, we held the label (target) constant. Thus, we will set the ref_effect to dose_group so that significance labels are generated to indicate differences in dose, not label.*
```{r}
# Define the order of the genomic targets for the x-axis: 
# We will order them from lowest to highest MF at the High dose.
label_order <- model_by_target$point_estimates %>%
 dplyr::filter(dose_group == "High") %>%
 dplyr::arrange(Estimate) %>%
 dplyr::pull(label)

# Define the order of the doses for the fill
dose_order <- c("Control", "Low", "Medium", "High")

plot <- plot_model_mf(model = model_by_target,
                      plot_type = "bar",
                      x_effect = "label",
                      plot_error_bars = TRUE,
                      plot_signif = TRUE,
                      ref_effect = "dose_group",
                      x_order = label_order,
                      fill_order = dose_order,
                      x_label = "Target",
                      y_label = "Mutation Frequency (mutations/bp)",
                      fill_label = "Dose",
                      plot_title = "",
                      custom_palette = c("#ef476f",
                                         "#ffd166",
                                         "#06d6a0",
                                         "#118ab2"))
plot <- plot + ggplot2::theme(axis.text.x = element_text(angle = 90))
```
![plot_model_mf Dose and Target Model](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot4.2.2.png)

## Benchmark Dose Modelling
Dose-response models are essential for quatitative risk assessment of mutagenicity, as they provide a framework to evaluate the levels at which exposure to a substance might cause an adverse effect. The **benchmark dose** (BMD) is a dose that produces a predetermined change in the measured response, defined as the **benchmark response** (BMR). The BMD is used as a point of departure to derive human health-based guidance values to inform regulatory risk assessment such as the reference dose (RfD), the derived no-effect level (DNEL) or the acceptable daily intake (ADI).

The BMD is estimated by applying various mathmatical models to fit the dose-response data. Some requirements must be met before modelling the BMD. There must be a clear dose-response trend in the MF data. We suggest using `model_mf()` to test for significant increases in MF with dose prior to running a BMD analysis. In general, studies with more dose groups and a graded monotonic response with dose will be more useful for BMD analysis. A minimum of three dose groups + 1 control group is suggested. Datasets in which a response is only observed at the high dose are usually not suitable for BMD modeling. However, if the one elevated response is near the BMR, adequate BMD computation may result. For a better estimate of the BMD, it is preferable to have studies with one or more doses near the level of the BMR.

Protection and safety authorities recommend the use of model averaging to determine the BMD and its confidence intervals. Model averaging incorporates information across multiple models to acount for model uncertainty, allowing the BMD to be more accurately estimated.

Ideally, the BMR would be based on a consensus scientific definition of what  minimal level of change in MF is biologically significant. Currently, the default provided by this package calculates the BMD at a 50% relative increase in MF from the background. This BMR was selected based on previous recommendations for genotoxicity assessment by White et al., 2020.

MutSeqR provides two functions for BMD modelling, each employing widely-used software designed to be consistent with methods used by regulatory authorities.
1) `bmd_proast()` runs a modified version of the [proast71.1](www.rivm.nl/en/proast) R library that is parametirized instead of menu-based.
2) `bmd_toxicr` uses the [ToxicR](https://github.com/NIEHS/ToxicR) library, available on Github.

### PROAST
`bmd_proast` will analyze the continuous, individual MF data following a log transformation. PROAST uses four families of nested models: exponential, Hill, inverse exponential, and log-normal. The Akaike information criterion (AIC) is used to select best fit. BMD confidence intervals are assessed by the Maximum Likelihood Profile method, or by model averaging via bootstrapping (recommended). BMR values are defined as a user-defined relative increase in MF from the control. Alternatively, users may set the BMR as one standard-deviation from the control.
#### General Usage: bmd_proast
Supply `bmd_proast` with per-sample mf data calculated using `calculate_mf()`, with the dose column retained in the summary table.

Specify the column that contains the numeric dose values in `dose_col`. The function can model more than one response variable at once. Supply all response variables to `response_col`. If you wish to include a covariate in the analysis, supply the covariate variable to `covariate_col`. PROAST will assess if the BMD values differ significantly between levels of the covariate and give a BMD estimate for each.

It is highly recommended to use model averaging when calculating the BMD confidence intervals. Specify the number of bootstraps to run with `num_bootstraps`. The recommended value is 200, but be aware that this may take some time to run. PROAST model averaging will return the upper and lower 90% BMD confidence intervals. MutSeqR calculates the model-averaged BMD value as the median BMD of all bootstrap runs.

Users may choose to generate model plots with `plot_results`. If TRUE, plots will be automatically saved to an output directory specified in `output_path`.

The function will output a data frame of final results, including a BMD estimate for each model family and the model averaging results, if applicable. Users may access the raw, unparsed PROAST results by setting `summary = FALSE`.

*Example 5.1. Calculate the BMD with model averaging for a 50% relative increase in MF from control. This will be calculated for both MFmin and MFmax*.
```{r}
# load example data:
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Calculate the MF per sample, retaining dose in the summary table.
mf_data <- calculate_mf(
            mutation_data = example_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            retain_metadata_cols = "dose"
            )

# Run PROAST            
proast_results <- bmd_proast(mf_data = mf_data,
                             dose_col = "dose",
                             response_col = c("mf_min", "mf_max"),
                             bmr = 0.5,
                             model_averaging = TRUE,
                             num_bootstraps = 200,
                             plot_results = TRUE,
                             output_path = "file.path.to.plot.output")
```

### ToxicR
`bmd_toxicr` will analyze continuous, individual or continuous, summary MF data, assuming either a normal (default) or log-normal distribution. This function employs Bayesian estimation using either the Laplace Maximum *a posteriori* approach (default) (Gelman et al. 1995) or Markov chain Monte Carlo (MCMC) simulation (Brooks et al. 2011). The default model parameter’s prior distributions are specified in (Wheeler et al. 2020), but they are user-modifiable. The default models are described in the European Food Safety Authority’s 2022 Guidance on the use of benchmark dose approach in risk assessment (EFSA Scientific Committee et al. 2022). Model averaging may be applied using described methodologies (Wheeler et al. 2020; 2022).

ToxicR offers several options for the BMR:
* **Relative deviation (*rel*)**: the BMD represents the dose that changes the mean MF a certain percentage from the background dose. 
* **Standard deviation (*sd*)**: the BMD represents the dose associated with the mean MF changing a specified number of standard deviations from the background mean. 
* **Absolute deviation (*abs*)**: the  BMD represents the dose associated with a specified absolute deviation from the background mean. 
* **Hybrid deviation (*hybrid*)**: the  BMD represents the dose that changes the probability of an adverse event by a specified amount. 

One of these options can be specified using the `bmr_type` parameter. The `bmr` parameter is set to a numeric value specifying the BMR, defined in relation to the calculation requested in `bmr_type`.

Selecting an appropriate BMR involves making judgements about the statistical and biological characteristics of the dataset and about the applications for which the resuling BMDs will be used.

#### General Usage: bmd_toxicr
To install ToxicR from Github
``` {r}
library(devtools)
install_github("NIEHS/ToxicR")
```
*See the ToxicR repository on github for more information on installing the package. Differences apply for mac users*

Supply `bmd_toxicr` with the per-sample mf data from `calculate_mf()`. Retain the dose column in the summary table. Specify the column that contains the numeric dose values in `dose_col`. The function can model more than one response variable at once. Supply all response variables to `response_col`.

Calculate the BMD with model averaging by setting `model_averaging = TRUE`. The confidence level for the upper and lower confidence intervals can be defined with the `alpha` parameter (default 90% CI).

Users may choose to generate model plots with `plot_results`. If TRUE, plots will be automatically saved to an output directory specified in `output_path`.

The function will return the BMD with its upper and lower confidence intervals for each response variable. If model averaging, a breakdown of the model averaging process can be returned alongside the results `ma_summary = TRUE`. This will return the estimate for each model with its associated posterior probability.

*Example 5.2. Calculate the BMD with model averaging for a 50% relative increase in MF from control. This will be calculated for both MFmin and MFmax*.
```{r}
# load example data:
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Calculate the MF per sample, retaining the dose in the summary table.
mf_data <- calculate_mf(
            mutation_data = example_data,
            cols_to_group = "sample",
            subtype_resolution = "none",
            retain_metadata_cols = "dose"
            )

# Run ToxicR
toxicr_results <- bmd_toxicr(mf_data = mf_data,
                             dose_col = "dose",
                             response_col = c("mf_min", "mf_max"),
                             bmr_type = "rel",
                             bmr = 0.5,
                             model_averaging = TRUE,
                             ma_summary = TRUE,
                             plot_results = TRUE,
                             output_path = "file.path.to.plot.output")
```

### Visualizing BMD Confidence Intervals
`plot_ci()` creates a confidence interval plot of BMD results for easy comparison of BMDs between response variables.

*Example 5.3. Compare the estimated BMD for MFmin (50% relative increase) by PROAST versus ToxicR*
```{r}
plot_results <- data.frame(Response = c("PROAST", "ToxicR"),
                           BMD = c(9.111, 9.641894),
                           BMDL = c(7.38, 8.032936),
                           BMDU = c(10.9, 10.97636))
plot <- plot_ci(data = plot_results,
                order = "asc",
                x_lab = "Dose (mg/kg-bw/d)",
                y_lab = "BMD Method")
```
![BMD CI plot](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot5.3.png)

## Mutation Spectra Analysis
The mutation spectra is the pattern of mutation subtypes within a sample or group. The mutation spectra can inform on the mechanisms involved in mutagenesis.

### Comparison of Mutation Spectra Between Groups
Compare the mutation spectra between experimental groups using the`spectra_comparison()` function. This function will compare the proportion of mutation subtypes at any resolution between user-defined groups using a modified contingency table approach (Piegorsch and Bailer, 1994). 

This approach is applied to the mutation counts of each subtype in a given group. The contingency table is represented as $R * T$ where $R$ is the number of subtypes, and $T$ is the number of groups. `spectra_comparison()` performs comparisons between $T = 2$ specified groups. The statistical hypothesis of homogeneity is that the proportion (count/group total) of each mutation subtype equals that of the other group. To test the significance of the homogeneity hypothesis, the $G^{2}$ likelihood ratio statistic is used: 

$$G^{2} = 2\  \sum_{i=1}^{R}\  \sum_{j=1}^{T}\  Y_{ij}\  log(\frac{Y_{ij}}{E_{ij}})$$

$Y_{ij}$ represents the mutation counts and $E_{ij}$ are the *expected* counts under the null hypothesis. The $G^{2}$ statistic possesses approximately a $\chi^{2}$ distribution in large sample sizes under the null hypothesis of no spectral differences. Thus, as the column totals become large, $G^{2}$ may be referred to a $\chi^{2}$ distribution with $(R -  1)(T - 1)$ degrees of freedom. It is important to note that the $G^{2}$ statistic may exhibit high false positive rates in small sample sizes when referred to a $\chi^{2}$ distribution. In such cases, we instead switch to an F-distribution. This has the effect of reducing the rate at which $G^{2}$ rejects each null hypothesis, providing greater stability in terms of false positive error rates. Thus when $N/(R-1) < 20$, where $N$  is the total mutation counts across both groups, the function will use a F-distribution, otherwise it will use a $\chi^{2}$-distribution.

This comparison assumes independance among the observations. Each tabled observation represents a sum of independent contributions to the total mutant count. We assume independance is valid for mutants derived from a  mixed population, however, mutants that are derived clonally from a single progenitor cell would violate this assumption. As such, it is recommended to use the **MFmin method** of mutation counting for spectral analyses  to ensure that all mutation counts are independant. In those cases where the independence may be invalid, and where additional, extra-multinomial sources of variability are present, more complex, hierarchical statistical models are required. This is currently outside the scope of this package.

#### General Usage: spectra_comparison()
The first step is to calculate the per-group mf data at the desired subtype_resolution using `calculate_mf()`. This mf_data is then supplied to `spectra_comparison()` along with a contrasts table to specify the comparisons. The contrasts table will consist of two columns, each specifying a group to be contrasted against the other. 

The function will output the $G^{2}$ statistic and p-value for each specified comparison listed in `constrasts`. P-values are adjusted for multiple comparison using the Sidak method.

*Example 6.1. In our example data, we are studying the mutagenic effect of BaP. Our samples were exposed to three doses of a BaP (Low, Medium, High), or to the vehicle control (Control). We will compare the base_6 snv subtypes, alongside non-snv variants, of each of the three dose groups to the control. In this way we can investigate if exposure to BaP leads to significant spectral differences.*
```{r}
# load example data:
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Calculate the MF per dose at the base_6 resolution.
mf_data <- calculate_mf(
            mutation_data = example_data,
            cols_to_group = "dose_group",
            subtype_resolution = "base_6"
            )

# Create the contrast table
contrasts_table <- data.frame(col1 = c("Low", "Medium", "High"),
                              col2 = c("Control", "Control", "Control"))
# Run the analysis
simple_spectra <- spectra_comparison(mf_data,
                                     cols_to_group = "dose_group",
                                     subtype_resolution = "base_6",
                                     mf_type = "min",
                                     contrasts = contrasts_table)
simple_spectra
# The base_6 spectra of all BaP dose groups is significantly different
# from Control.                                    
```

### Mutational Signatures Analysis
Mutational processes generate characteristic patterns of mutations, known as mutational signatures. Distinct mutational signatures have been extracted from various cancer types and normal tissues using data from the Catalogue of Somatic Mutations in Cancer, ([COSMIC](https://cancer.sanger.ac.uk/signatures/)) database. These include signatures of single base substitutions (SBSs), doublet base substitutions (DBSs), small insertions and deletions (IDs) and copy number alterations (CNs). It is possible to assign mutational signatures to individual samples or groups using the `signature_fitting` function. Linking ECS mutational profiles of specific mutagens to existing mutational signatures provides empirical evidence for the contribution of environmental mutagens to the mutations found in human cancers and informs on mutagenic mechanisms.

The `signature_fitting` function utilizes the [SigProfiler](https://github.com/AlexandrovLab) suite of tools (Díaz-Gay et al. 2023; Khandekar et al. 2023) to assign SBS signatures from the COSMIC database to the 96-base SNV subtypes of a given group by creating a virtual environment to run python using `reticulate`. `signature_fitting()` facilitates interoperability between these tools for users less familiar with python and assists users by coercing the mutation data to the necessary structure for the SigProfiler tools. 

This function will  install several python dependencies using a conda virtual environment on first use, as well as the FASTA files for all chromosomes for your specified reference genome. As a result ~3Gb of storage must be available for the downloads of each genome. 

SNVs in their 96-base trinucleotide context are summed across groups to create a mutation count matrix by `SigProfilerMatrixGeneratorR()` (SigProfilerMatrixGenerator; Khandekar et al. 2023). `Analyze.cosmic_fit ` (SigProfilerAssignment; Díaz-Gay et al. 2023) is then run to assign mutational signatures to each group using refitting methods, which quantifies the contribution of a set of signatures to the mutational profile of the group. The process is a numerical optimization approach that finds the combination of mutational signatures that most closely reconstructs the mutation count matrix. To quantify the number of mutations imprinted by each signature, the tool uses a custom implementation of the forward stagewise algorithm and it applies nonnegative least squares, based on the Lawson-Hanson method. Cosine similarity values, and other solution statistics, are generated to compare the reconstructed mutational profile to the original mutational profile of the group, with cosine values > 0.9 indicating a good reconstruction.

Currently, `signature_fitting` offers fitting of COSMIC version 3.4 SBS signatures to the SBS96 matrix of any sample/group. For advanced use, including using a custom set of reference signatures, or fitting the DBS, ID, or CN signatures, it is suggested to use the SigProfiler python tools directly in as described in their respective [documentation](https://github.com/AlexandrovLab).

#### General Usage: signature_fitting()
`signature_fitting` requires the installation of reticulate and [SigProfilerMatrixGeneratorR](https://github.com/AlexandrovLab/SigProfilerMatrixGeneratorR) as well as a version of python 3.8 or newer.
```{r}
# Install reticulate
install.packages("reticulate")

# Install python
reticulate::install_python()

# Install SigProfilerMatrixGeneratorR from github using devtools.
library(devtools)
install_github("AlexandrovLab/SigProfilerMatrixGeneratorR")

```

`signature_fitting` will take the imported (and filtered, if applicable) `mutation_data` and create mutational matrices for the somatic SNV mutations. Mutations are summed across levels of the `group` parameter. This can be set to individual samples or to an experimental group. 

The `project_genome` will be referenced for the creation of the mutational matrices. The reference genome will be installed if not already.

The virtual environment can be specified with the `env_name` parameter. If no such environmnent exists, then the function will create one in which to store the dependencies and run the signature refitting. Specify your version of python using the `python_version` parameter (must be 3.8 or higher).

*Example 6.2. Determine the COSMIC SBS signatures associated with each BaP dose group.*
```{r}
# Load the example data.
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Run Analysis
signature_fitting(mutation_data = example_data,
                  project_name = "Example",
                  project_genome = "mm10",
                  env_name = "MutSeqR",
                  group = "dose_group",
                  python_version = "3.11",
                  output_path = NULL)
```

#### Signature_fitting() Output
Results from the `signature_fitting` will be stored in an output folder. A filepath to a specific output directory can be designated using the `output_path` parameter.  If null, the output will be stored within your working directory. Results will be organized into subfolders based on the `group` parameter. The output structure is divided into three folders: input, output, and logs. 

The input folder contains "mutations.txt", a text file with the mutation_data coerced into the required format for SigProfilerMatrixGenerator. It consists of a list of all the snv variants in each `group` alongside their genomic positions. This data serves as input for matrix generation.

The log folder contains the error and the log files for SigProfilerMatrixGeneration.

The output folder contains the results from matrix generation and signature refitting desribed in detail below. 

**Matrix Generation Output**

`signature_fitting` uses [SigProfilerMatrixGenerator](https://osf.io/s93d5/wiki/home/) to create mutational matrices (Bergstrom et al., 2019). Mutation matrices are created for single-base substitutions (SBS) and doublet-base substitutions (DBS), including matrices with extended sequence context and transcriptional strand bias. SBS and DBS matrices are stored in their respective folders in the output directory. **Only the SBS96 matrix is used for refitting**. 

**Single-Base Substitution Matrices (SBS):**
Generated matrices are described below. *Matrices are stored as `.all` files which can be viewed in a text-editor like notepad.*
| File|                                                               |
|-----|---------------------------------------------------------------|
|*project_name.SBS6.all*  | 6-base. The 6 pyrimidine single-nucleotide variants. *C>A, C>G, C>T, T>A, T>C, or T>G*|
|*project_name.SBS18.all*| 6-base. The 6 pyrimidine single-nucleotide variants within 3 transcriptional bias categories: Untranscribed (U), Transcribed (T), Non-Transcribed Region (N).|
|*project_name.SBS24.all*| 6-base. The 6 pyrimidine single-nucleotide variants within 4 transcriptional bias categories: Untranscribed (U), Transcribed (T), Bidirectional (B), Non-Transcribed Region (N).|
|**project_name.SBS96.all**| 96-base. The 6 pyrimidine single-nucleotide variants alongside their flanking nucleotides (4 x 4 = 16 combinations). *Ex. A[C>G]T*|
|*project_name.SBS288.all*|96-base. The 96-base single-nucleotide variants within 3 transcriptional bias categories (U, T, N).|
|*project_name.SBS384.all* |96-base. The 96-base single-nucleotide variants within 4 transcriptional bias categories (U, T, N, B).|
|*project_name.SBS1536.all* | 1536-base. The 6 pyrimidine single-nucleotide variants alongside their flanking dinucleotides (16 x 16 = 256 combinations). *Ex. AA[C>G]TT*|
|*project_name.SBS4608.all* | 1536-base. The 1536-base single-nucleotide variants within 3 transcriptional bias categories (U, T, N).|
|*project_name.SBS6144.all* | 1536-base. The 1536-base single-nucleotide variants within 4 transcriptional bias categories (U, T, N, B).|

**Doublet-base Matrices (DBS)**: DBS are somatic mutations in which a set of two adjacent DNA base-pairs are simultaneously substituted with another set of two adjacent DNA base-pairs. **We do not recommend using the DBS matrices generated using `signature_fitting` for further analysis.** The `signature_fitting` function is designed to handle only SBS mutations. All true MNVs, including doublets, are filtered out of the `mutation_data` prior to MatrixGeneration. However, the tool will still attempt to identify DBSs and will occasionally find two independent SBSs occuring next to each other simply by chance. If you wish to use DBS mutations in your signature analysis, please refer directly to the SigProfiler tools.

*SigProfilerMatrix Generator also supports small indels, structural variants, and copy-number variants. `signature_fitting` will not generate these matrices. If you wish to utilise these features, please refer directly to the SigProfiler tools.*

Barplots of the mutation matrices for all groups can be found in the "plots" folder. The number of mutations are plotted for each group at the various subtype resolutions

**vcf_files**: This output folder provides text-based files containing the original mutations and their SigProfilerMatrixGenerator classification for each chromosome. The files are separated into dinucleotides (DBS), multinucleotide substitutions (MNS), and single nucleotide variants (SNV) folders containing the appropriate files. The headers are:
1) The group
2) the chromosome
3) the position
4) the SigProfilerMatrixGenerator classification
5) the strand {1, 0, -1}.

The headers for each file are the same with the exception of the MNS files which don't contain a matrix classification or a strand classification. As noted above the DBS and MNS matrices do no reflect the true mutation counts for these variant types. Only SBS/SNV mutations are included in the matrix generation.

**Transcription Strand Bias (TSB)**: SBS mutations will be tested for [transcription strand bias](https://osf.io/s93d5/wiki/5.%20Output%20-%20TSB/). These results will be stored in the `TSB` folder.

It is not possible to distinguish which of the two DNA strands a mutation originated on. However, one expects that mutations from the same type will be equally distributed across the two DNA strands. In other words, we expect mutations to occur at the same proportion as their reverse complement.

*Example, given a mutational process that causes purely **C>G:T:A** mutations, and a long repetitive sequence on the reference genome:*

**5'- CGCGCGCGCGCGCGCGCGCGCG-3'**

*One would expect to see an equal number of **C>T** and **G>A** mutations.*

However, in many cases, an asymteric number of mutations are observed due to one of the strands having a higher propensity for being either damaged or repaired. A common example of this is transcription strand bias where the transcribed strand is subjected to higher rates of DNA repair as part of the transcriptional process compared to the untranscribed strand.

`SigProfilerMatrixGenerator` evaluates the transcriptional strand classification of mutations within well-annotated protein coding genes of a reference genome. Mutations that occur outside of coding regions are classified as *Non-transcribed* (N). Mutations that occur within coding regions are classified as one of; *transcribed* (T), *un-transcribed* (U), *bi-directional* (B), or unknown. In order to classify mutations within coding regions, mutations are oriented based on the reference strand and their pyrimidine context. 

*For example, consider that our reference sequence contains the coding sequence of a gene i.e. it is the coding/UN-transcribed DNA strand. A **T>C** mutation called on this reference sequence would be referred to as an untranscribed **T>C** (**U:T>C**). However, if instead an **A>G** mutation is called on this reference sequence, it would be referred to as a transcribed **T>C** (**T:T>C**). Purine-based mutations called from the reference sequence are converted to their pyrimidine context and this includes swapping to the complementary DNA strand. In this case, the reverse  complement of untranscribed **A>G** is transcribed **T>C**.*

In rare cases, both strands of a genomic region code for a gene. Such mutations are annotated as bidirectional based on their pyrimidine context. *For example, both **T:A>C:G** and **A:T>G:C** mutations in regions of bidirectional transcription will be annotated as bidirectional **T>C** (**B:T>C**) mutations.* 

All SBS mutations will be classified within the four transcriptional bias categories:

|Category        | Description                                           |
|----------------|-------------------------------------------------------|
|Transcribed (T) |The variant is on the transcribed (template) strand.   |
|Untranscribed (U)|The variant is on the untranscribed (coding) strand.  |
|Bidirectional (B) |The variant is on both strands and is transcribed either way. |
|Nontranscribed (N) |The variant is in a non-coding region and is untranslated. |

The tool will then perform a transcription strand bias test which compares the number of transcribed and untranscribed mutations for each mutation type. *For example, it will compare the number of transcribed T>C to untranscribed T>C mutations. Should there be a significant difference, it would indicate that T:A>C:G mutations are occuring at a higher rate on one of the strands compared to the other.* Transcription strand bias tests will be included for the 6-base, 96-base and 1536-base SBS mutation contexts. 

The output files contain the following information:
* the `group`
* the mutation type
* the enrichment value (# Transcribed / # untranscribed)
* the p-value, corrected for multiple comparisons using the false discrover rate method
* the false discovery rate q-value


Files include:
| File | Description |
|------|-------------|
|*strandBiasTes_24.txt*| stats of the SBS6 variants|
|*strandBiasTes_384.txt*| stats of the SBS96 variants|
|*strandBiasTes_6144.txt*|stats of the SBS1536 variants|
|*significantResults_strandBiasTest.txt*|returns significant results from the three files above.|

**Signature Refitting Results**

Results from the signature refitting perfomed by [SigProfilerAssignment](https://osf.io/mz79v/wiki/home/) will be stored within the "Assignment_Solution" folder. "Assignment_Solution" consists of 3 subdirectories;  "Activities", "Signatures", and "Solution_Stats". 

**Activities**

| File | Description |
|------|-------------|
|*Assignment_Solution_Activities.txt* | This file contains the activity matrix for the selected signatures. The first column lists all of the samples/groups. All of the following columns list the calculated activity value for the respective signatures. Signature activities correspond to the specific numbers of mutations from the sample's original mutation matrix caused by a particular mutational process. |
|*Assignment_Solution_Activity_Plots.pdf* | This file contains a stacked barplot showing the number of mutations in each signature on the y-axis and the samples/groups on the x-axis.|
| *Assignment_Solution_TMB_plot.pdf* | This file contains a tumor mutational burden plot. The y-axis is the somatic mutations per megabase and the x-axis is the number of samples/groups plotted over the total number of samples/groups included. The column names are the mutational signatures and the plot is ordered by the median somatic mutations per megabase. |
| *Decomposed_Mutation_Probabilities.txt* | This file contains the probabilities of each of the 96 mutation types in each sample/group. The probabilities refer to the probability of each  mutation type being caused by a specific signature. The first column lists all the samples/groups, the second column lists all the mutation types, and the following columns list the calculated probability value for the respective signatures. |
| *SampleReconstruction* | This folder contains generated plots for each sample/group summarizing the assignment results. Each plot consists of three panels. (i) Original: a bar plot of the inputted 96SBS mutation matrix for the sample/group. (ii) Reconstructed: a bar plot of the reconstruction of the original mutation matrix. (iii) The mutational profiles for each of the known mutational signatures assigned to that sample/group, including the activities for each signature. Accuracy metrics for the reconstruction are displayed at the bottom of the figure.|

**Signatures**
| Files | Description |
|-------|-------------|
| *Assignment_Solution_Signatures.txt* | The distribution of mutation types in the input mutational signatures. The first column lists all 96 of the mutation types. The following columns are the signatures. |
| *SBS_96_plots_Assignment_Solution.pdf* | Barplots for each signature identified that depicts the proportion of the mutation types for that signature. The top right corner also lists the total number of mutations and the percentage of total mutations assigned to the mutational signature. |

**Solution_Stats**
| Files | Description |
|-------|-------------|
| *Assignment_Solution_Samples_Stats.txt* | The accuracy metrics for the reconstruction. statistics for each sample including the total number of mutations, cosine similarity, L1 norm (calculated as the sum of the absolute values of the vector), L1 norm percentage, L2 norm (calculated as the square root of the sum of the squared vector values), and L2 norm percentage, along with the Kullback-Leibler divergence. |
| *Assignment_Solution_Signature_Assignment_log.txt* | The events that occur when known signatures are assigned to an input sample. The information includes the L2 error and cosine similarity between the reconstructed and original sample within different composition steps. |

**Other Files**

*JOB_METADATA_SPA.txt*: This file contains the metadata about system and runtime.

#### Get the Input Files for the SigProfiler Webtool
Users may choose to use the [SigProfiler Webtool](https://cancer.sanger.ac.uk/signatures/assignment/) instead of using the `signature_fitting() `function. MutSeqR offers functions to coerce mutation data into the proper format for input files.

**Mutation Calling File**

`write_mutation_calling_file()` creates a simple text file from mutation data that can be used for mutation signatures analysis using the SigProfiler Assignment web application as a "mutation calling file". Signature analyses are done at the sample level when using mutation calling files. The file will be saved to your output directory, specified in `output_path`.

*Example 6.3. Analyze the COSMIC SBS signatures contributing to each of the 24 samples using the SigProfiler Web Tool. Output a mutation calling file that can be uploaded to the webtool.*
```{r}
# Load example data
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Get input file
write_mutation_calling_file(mutation_data = example_data,
                            project_name = "Example",
                            project_genome = "mm10",
                            output_path = NULL)
```

**Mutational Matrix**

 `write_mutational_matrix()` will sum mutations across user-defined groups before coercing the data into the proper format for input as a "mutational matrix". SNV subtypes can be resolved to either the base_6 or base_96 resolution. The file is saved to the specified output directory.

 *Example 6.4.  Analyze the COSMIC SBS signatures contributing to each dose group using the SigProfiler Web Tool. Output a mutational matrix that can be uploaded to the webtool.*
```{r}
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)
write_mutational_matrix(mutation_data = example_data,
                        group = "dose_group",
                        subtype_resolution = "base_96",
                        mf_type = "min",
                        output_path = NULL)
```
### Visualizing Mutation Spectra
The mutation spectra can be visualized with  `plot_spectra` which will create a stacked bar plot for user-defined groups at the desired subtype resolution. Mutation subtypes are represented by colour. The value can represent subtype count (`sum`), frequency (`mf`), or `proportion`. 

**Hierarchical Clustering**

`plot_spectra()` integrates `cluster_spectra()` which performs unsupervised hierarchical clustering of samples based on the mutation spectra. `cluster_spectra()` uses `dist()` from the stats library to compute the sample-to-sample distances using a user-defined distance measure (default Euclidean). The resulting distance matrix is passed to `hclust()` to cluster samples using the specified linkage method (default Ward). The function will output a dendrogram visually representing the clusters' relationships and hierarchy. The dendrogram will be overlaid on the `plot_spectra()` bar plot and samples will be ordered accordingly.

*Example 6.5. Plot the base_6 proportions for each dose group.*
```{r}
# load the example data
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Calculate the mf data at the 6-base resolution for each dose
# We will exclude ambiguous or uncategorized variants
mf_data <- calculate_mf(mutation_data = example_data,
                        cols_to_group = "dose_group",
                        subtype_resolution = "base_6",
                        variant_types = c("-ambiguous", "-uncategorized"))
# Set the desired order for the dose group:
mf_data$dose_group <- factor(mf_data$dose_group,
                             levels = c("Control",
                                        "Low",
                                        "Medium",
                                        "High"))
# Plot
plot <- plot_spectra(mf_data = mf_data,
                     group_col = "dose_group",
                     subtype_resolution = "base_6",
                     response = "proportion",
                     group_order = "arranged",
                     group_order_input = "dose_group",
                     x_lab = "Dose Group",
                     y_lab = "Subtype Proportion")                       
```
![plot_spectra](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot6.5.png)

*Example 6.6. Plot the base_6 mutation spectra per sample, with hierarchical clustering. For this example we have created a new sample column with more intuitive sample names: new_sample_id. These names correspond to their associated dose groups. We will see that samples largly cluster within their dose groups.*
```{r}
# load the example data
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Calculate the mf data at the 6-base resolution for each sample
mf_data <- calculate_mf(mutation_data = example_data,
                        cols_to_group = "new_sample_id",
                        subtype_resolution = "base_6",
                        variant_types = c("-ambiguous", "-uncategorized"))
# Plot
plot <- plot_spectra(mf_data = mf_data,
                     group_col = "new_sample_id",
                     subtype_resolution = "base_6",
                     response = "proportion",
                     group_order = "clustered",
                     x_lab = "Sample",
                     y_lab = "Subtype Proportion")                          
```                        
![plot_spectra with Clustering](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot6.6.png)

The 96-base SNV mutation subtypes can be vizualised using `plot_trinucleotide()`. This function creates a bar plot of the 96-base SNV spectrum for all levels of a user-defined group. Data can represent subtype mutation count (`sum`), frequency (`mf`), or `proportion`. Aesthetics are consistent with COSMIC trinucleotide plots. Plots are automatically saved to the specified output directory.

*Example 6.7. plot the base_96 mutation spectra proportions for each dose group.*
```{r}
# load the example data
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Calculate the mf data at the 96-base resolution for each dose
mf_data <- calculate_mf(mutation_data = example_data,
                        cols_to_group = "dose_group",
                        subtype_resolution = "base_96",
                        variant_types = "snv")
# Plot
plot_trinucleotide(mf_96 = mf_data,
                   group_col = "dose_group",
                   response = "proportion",
                   mf_type = "min",
                   output_path = "file.path.to.output.folder")                       
```
![plot_trinucleotide High Dose](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot6.7.png)


Another option for vizualizing the base-96 mutation spectra is `plot_trinucleotide_heatmap()`. This function creates a heatmap of the 96-base SNV proportions. Plots can be facetted by additional grouping variables. Heatmaps are useful for making comparisons between experimental variables when information density becomes too high to represent using traditional plots.

*Example 6.8. Plot the 96-base SNV spectrum for each sample, facetted by dose group.*
```{r}
# load the example data
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

# Calculate the mf data at the 96-base resolution for each sample
mf_data <- calculate_mf(mutation_data = example_data,
                        cols_to_group = "sample",
                        subtype_resolution = "base_96",
                        variant_types = "snv",
                        retain_metadata_cols = "dose_group")
mf_data$dose_group <- factor(mf_data$dose_group,
                             levels = c("Control",
                                        "Low",
                                        "Medium",
                                        "High"))                        
plot <- plot_trinucleotide_heatmap(mf_data = mf_data,
                                   group_col = "sample",
                                   facet_col = "dose_group")            
```
![plot_trinucleotide_heatmap](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot6.8.png)

## Visualize Recurrent Mutations
`plot_bubbles` is used to visually represent the distribution and density of recurrent mutations. Each mutation is in a given group is represented by a bubble whose size is scaled on either the `alt_depth` or the `vaf`. Thus a highly reccurent mutation is represented by a large bubble. These plots make it easy to determine if MFmax is driven by a few highly recurrent mutations versus serveral moderately recurrent mutations.
Plots can be facetted by user-defined groups, and bubbles can be coloured by any variable of interest to help discern patterns in mutation recurrence.

*Example 7. Plot mutations per dose group, bubbles coloured by base-6 subtype*
```{r}
# load the example data
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)
plot <- plot_bubbles(mutation_data = example_data,
                     size_by = "alt_depth",
                     facet_col = "dose_group",
                     color_by = "normalized_subtype")
```
![plot_bubbles](https://github.com/EHSRB-BSRSE-Bioinformatics/MutSeqR/blob/smm_seq/inst/extdata/Example_files/plot7.png)

## Retrieve Sequences of genomic target regions
`get_seq()` will retrive raw nucleotide sequences for specified genomic intervals. This function will install an appropriate BS genome library to retrieve sequences based on species, genome, and masked parameter.

TwinStrand's Mutagenesis Panels are stored in package files and can easily be retrieved. 

Sequences are returned within a *GRanges* object.

*Example 8.1. Retrieve the sequences for our example's target panel, TwinStrand's Mouse Mutagenesis Panel*
```{r}
regions_seq <- get_seq(regions = "TSpanel_mouse")
```

*Example 8.2. Retrieve sequences for a custom interval of regions. We will use the Human Mutagenesis Panel as an example.*
```{r}
# We will load the TSpanel_human regions file as an example
human <- load_regions_file("TSpanel_human")
regions_seq <- get_seq(regions = "custom",
                       custom_regions = human,
                       is_0_based = TRUE,
                       species = "human",
                       genome = "hg38",
                       masked = FALSE,
                       padding = 0)
```

Sequences can be exported as FASTA files with `write_reference_fasta()`. Supply this function with the GRanges object with the sequences of the regions. Each one will be written to a single FASTA file.

```{r}
write_reference_fasta(regions_seq, output_path = NULL)
```
## Exporting Results

Users can easily output data frames to an Excel workbook with `write_excel()`. This function can write single data frames or it can take a list of dataframes and write each one to a separate Excel sheet in a workbook.

In addition to data frames, `write_excel()` will also extract the mf_data, point_estimates, and pairwise_comparisons from `model_mf()` output to write to an excel workbook. Set `model_results` to TRUE if supplying the function with the output to model_mf().

*Example 9.1. Write MF data to excel workbook.*
```{r}
# Load the example data
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)

mf1 <- calculate_mf(example_data,
                    cols_to_group = "sample",
                    subtype_resolution = "none",
                    retain_metadata_cols = "dose")
mf2 <- calculate_mf(example_data,
                    cols_to_group = "sample",
                    subtype_resolution = "base_6",
                    variant_types = c("-ambiguous", "-uncategorized"))
mf3 <- calculate_mf(example_data,
                    cols_to_group = "dose",
                    subtype_resolution = "base_96",
                    variant_types = "snv")

# save a single data frame to an Excel file
write_excel(mf1, output_path, workbook_name = "test_single")

# Write multiple data frames to a list to export all at once.                    
list <- list(mf1, mf2, mf3)
names(list) <- c("mf1", "mf2", "mf3")

#save a list of data frames to an Excel file
write_excel(list, output_path, workbook_name = "test_list")

```

*Example 9.2. Export model results*
```{r}
# Run the model
model  <- model_mf(mf1,
                   fixed_effects = "dose",
                   reference_level = 0,
                   contrasts = data.frame(col1 = c(12.5, 25, 50),
                                          col2 = rep(0,3)))
write_excel(model,
            workbook_name = "Example_model",
            model_results = TRUE)
```

Mutation data can be written to a VCF file for downstream applications with `write_vcf_from_mut()`.

```{r}
example_file <- system.file("extdata",
                            "example_mutation_data_filtered.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)
 
write_vcf_from_mut(example_data)
```

# References
Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG, Alexandrov LB. SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. BMC Genomics. 2019 Aug 30;20(1):685. doi: 10.1186/s12864-019-6041-2. PMID: 31470794; PMCID: PMC6717374.

Brooks, Steve, Andrew Gelman, Galin Jones, and Xiao-Li Meng. 2011. Handbook of Markov Chain Monte Carlo. CRC Press.

Danecek, Petr, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, Eric Banks, Mark A. DePristo, Robert E. Handsaker, et al. 2011. “The Variant Call Format and VCFtools.” Bioinformatics 27 (15): 2156–58. https://doi.org/10.1093/bioinformatics/btr330.

Díaz-Gay M, Vangara R, Barnes M, Wang X, Islam SMA, Vermes I, Duke S, Narasimman NB, Yang T, Jiang Z, Moody S, Senkin S, Brennan P, Stratton MR, Alexandrov LB. Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment. Bioinformatics. 2023 Dec 1;39(12):btad756. doi: 10.1093/bioinformatics/btad756. PMID: 38096571; PMCID: PMC10746860.

Dodge, Annette E., Danielle P. M. LeBlanc, Gu Zhou, Andrew Williams, Matthew J. Meier, Phu Van, Fang Yin Lo, et al. 2023. “Duplex Sequencing Provides Detailed Characterization of Mutation Frequencies and Spectra in the Bone Marrow of MutaMouse Males Exposed to Procarbazine Hydrochloride.” Archives of Toxicology 97 (8): 2245–59. https://doi.org/10.1007/s00204-023-03527-y.

EFSA Scientific Committee, Simon John More, Vasileios Bampidis, Diane Benford, Claude Bragard, Thorhallur Ingi Halldorsson, Antonio F Hernández-Jerez, et al. 2022. “Guidance on the Use of the Benchmark Dose Approach in Risk Assessment.” EFSA Journal 20 (10): e07584. https://doi.org/10.2903/j.efsa.2022.7584.

Gelman, Andrew, John B. Carlin, Hal S. Stern, and Donald B. Rubin. 1995. Bayesian Data Analysis. New York: Chapman and Hall/CRC. https://doi.org/10.1201/9780429258411.

Halekoh, Ulrich, and Søren Højsgaard. 2024. “doBy: Groupwise Statistics, LSmeans, Linear Estimates, Utilities.” https://cran.r-project.org/web/packages/doBy/index.html.

Kennedy, Scott R., Michael W. Schmitt, Edward J. Fox, Brendan F. Kohrn, Jesse J. Salk, Eun Hyun Ahn, Marc J. Prindle, et al. 2014. “Detecting Ultralow-Frequency Mutations by Duplex Sequencing.” Nature Protocols 9 (11): 2586–2606. https://doi.org/10.1038/nprot.2014.170.

Khandekar, Azhar, Raviteja Vangara, Mark Barnes, Marcos Díaz-Gay, Ammal Abbasi, Erik N. Bergstrom, Christopher D. Steele, Nischalan Pillay, and Ludmil B. Alexandrov. 2023. “Visualizing and Exploring Patterns of Large Mutational Events with SigProfilerMatrixGenerator.” BMC Genomics 24 (1): 469. https://doi.org/10.1186/s12864-023-09584-y.

LeBlanc DPM, Meier M, Lo FY, Schmidt E, Valentine C 3rd, Williams A, Salk JJ, Yauk CL, Marchetti F. Duplex sequencing identifies genomic features that determine susceptibility to benzo(a)pyrene-induced in vivo mutations. BMC Genomics. 2022 Jul 28;23(1):542. doi: 10.1186/s12864-022-08752-w. PMID: 35902794; PMCID: PMC9331077.

Marchetti F, Cardoso R, Chen CL, Douglas GR, Elloway J, Escobar PA, Harper T Jr, Heflich RH, Kidd D, Lynch AM, Myers MB, Parsons BL, Salk JJ, Settivari RS, Smith-Roe SL, Witt KL, Yauk C, Young RR, Zhang S, Minocherhomji S. Error-corrected next-generation sequencing to advance nonclinical genotoxicity and carcinogenicity testing. Nat Rev Drug Discov. 2023 Mar;22(3):165-166. doi: 10.1038/d41573-023-00014-y. PMID: 36646809.

Marchetti F, Cardoso R, Chen CL, Douglas GR, Elloway J, Escobar PA, Harper T Jr, Heflich RH, Kidd D, Lynch AM, Myers MB, Parsons BL, Salk JJ, Settivari RS, Smith-Roe SL, Witt KL, Yauk CL, Young R, Zhang S, Minocherhomji S. Error-corrected next generation sequencing - Promises and challenges for genotoxicity and cancer risk assessment. Mutat Res Rev Mutat Res. 2023 Jul-Dec;792:108466. doi: 10.1016/j.mrrev.2023.108466. Epub 2023 Aug 27. PMID: 37643677.

Menon, Vijay, and Douglas E. Brash. 2023. “Next-Generation Sequencing Methodologies to Detect Low-Frequency Mutations: ‘Catch Me If You Can.’” Mutation Research/Reviews in Mutation Research 792 (July):108471. https://doi.org/10.1016/j.mrrev.2023.108471.

Piegorsch WW, Bailer AJ. Statistical approaches for analyzing mutational spectra: some recommendations for categorical data. Genetics. 1994 Jan;136(1):403-16. doi: 10.1093/genetics/136.1.403. PMID: 8138174; PMCID: PMC1205789.

Wheeler, Matthew W., Todd Blessinger, Kan Shao, Bruce C. Allen, Louis Olszyk, J. Allen Davis, and Jeffrey S. Gift. 2020. “Quantitative Risk Assessment: Developing a Bayesian Approach to Dichotomous Dose-Response Uncertainty.” Risk Analysis: An Official Publication of the Society for Risk Analysis 40 (9): 1706–22. https://doi.org/10.1111/risa.13537.

Wheeler, Matthew W., Jose Cortinas, Marc Aerts, Jeffery S. Gift, and J. Allen Davis. 2022. “Continuous Model Averaging for Benchmark Dose Analysis: Averaging Over Distributional Forms.” Environmetrics 33 (5): e2728. https://doi.org/10.1002/env.2728.

White PA, Long AS, Johnson GE. Quantitative Interpretation of Genetic Toxicity Dose-Response Data for Risk Assessment and Regulatory Decision-Making: Current Status and Emerging Priorities. Environ Mol Mutagen. 2020 Jan;61(1):66-83. doi: 10.1002/em.22351. Epub 2019 Dec 19. PMID: 31794061.






