---
editor_options: 
  markdown: 
    wrap: 72
---

# Duplex Sequencing Analysis

## Pre-abmle

To do: write a bit about the project, the scope, why it's necessary.

Provide a little background on duplex sequencing, describe the
technology briefly, and give some context for the type of data we are
meant to be processing with this package.

## Installation

Install from github with:

```{r}
# install.packages("devtools")
devtools::install_github("EHSRB-BSRSE-Bioinformatics/duplex-sequencing", auth_token = "your personal_access_token from github")
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
load in the metadata for TwinStrand's Mutagenesis Panel^(TM)^ using the
regions parameter; "mouse" or "human".

```{r}
library(DupSeqR)
# mut_data <- "file path to .mut file"
# sample_data <- "file path to sample meta data"
mutation_data <-
  import_mut_data(mut_file = mut_data,
                  sample_data_file = sample_data,
                  vaf_cutoff = 0.1,
                  regions = "mouse" # Mouse Mutagenesis Panel
                  )
```

If you are not using one of TwinStrand's Mutagenesis Panels^(TM)^, then
you will add your target regions' metadata using a custom_regions_file.
Use parameters to indicate your file's file path, delimiter, and whether
the region ranges are 0-based or 1-based.

```{r}
# mut_data <- "file path to .mut file"
# sample_data <- "file path to sample meta data"
mutation_data <-
  import_mut_data(mut_file = mut_data,
                  sample_data_file = sample_data,
                  vaf_cutoff = 0.1,
                  regions = "custom",
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

| Column      | Definition                               | Synonyms            |
|-------------|------------------------------------------|---------------------|
| Name        |                                          |                     |
| contig      | The reference sequence name.              | chr                 |
|             |                                          | chromosome          |
|             |                                          |  seqnames           |
| start       | The 0-based start position of the feature | position            |
| end         | The half-open end position of the feature  |                     |
| sample      | The sample name.                          | sample_name         |
|             |                                          |                     |
| ref         | The reference allele at this position.    |                     |
| alt         | The left-aligned, alternate allele at     | alt.value           |
|             | this position.                            |                     |
| alt_depth   | The read depth supporting the alternate   | var_depth           |
|             | allele.                                  |                     |
| depth_col   | The total read depth at this position.    | informat            |
|             | This column can be total_depth (excluding | ive_somatic_depth   |
|             | N-calls) or depth (including N-calls; if  | = total_depth       |
|             | total_depth is not available).            |                     |
| variation_type | The category to which this variant is   | type                |
|             | assigned.                                |                     |
| context     | The local reference trinucleotide        | sequence_context   |
|             | context at this position (e.g. ATC - not  | flanking_sequence  |
|             | necessarily the transcript codon).         |                     |


### Importing .vcf files

Required fields for your `.vcf` file are listed in the table below.

+------------+-----------+---------------------------------------------+
|            | **Field   | **Definition**                              |
|            | Name**    |                                             |
+------------+-----------+---------------------------------------------+
| **FIXED    | `contig`  | The reference sequence name.                |
| FIELDS**   |           |                                             |
+------------+-----------+---------------------------------------------+
|            | `start`   | The 0-based start position of the feature   |
|            |           | in contig.                                  |
+------------+-----------+---------------------------------------------+
|            | REF       | The reference allele at this position.      |
+------------+-----------+---------------------------------------------+
|            | `ALT`     | The left-aligned, alternate allele at this  |
|            |           | position.                                   |
+------------+-----------+---------------------------------------------+
| **FORMAT   | `AD`      | The allelic depths for the reference and    |
| FIELDS**   |           | alternate alleles in the order listed.      |
+------------+-----------+---------------------------------------------+
|            | `DP`      | The total read depth at this position       |
|            |           | (excluding N-calls). Equivalent to `depth`. |
+------------+-----------+---------------------------------------------+
|            | `VD`      | Variant Depth. Equivalent to `alt_depth`.   |
+------------+-----------+---------------------------------------------+
| **INFO     | `TYPE`    | The category to which this variant is       |
| FIELDS**   |           | assigned. Equivalent to `variation_type`.   |
+------------+-----------+---------------------------------------------+
|            | `END`     | The half-open end position of the feature   |
|            |           | in contig.                                  |
+------------+-----------+---------------------------------------------+
| *          | sample    | An identifying field for your samples;      |
| *SUGGESTED |           | either in the INFO field or as the header   |
| INFO       |           | to the FORMAT field.                        |
| FIELDS**   |           |                                             |
+------------+-----------+---------------------------------------------+
|            | `SVTYPE`  | Structural variant types; INV DUP DEL INS   |
|            |           | FUS.                                        |
+------------+-----------+---------------------------------------------+
|            | SVLEN\`   | Length of the structural variant in base    |
|            |           | pairs                                       |
+------------+-----------+---------------------------------------------+

The column variation_type/TYPE may contain these values:

+----------------+-----------------------------------------------------+
| `v             | Definition                                          |
| ariation_type` |                                                     |
+================+=====================================================+
| no_variant     | No variation, the null-case.                        |
+----------------+-----------------------------------------------------+
| snv            | Single nucleotide variant.                          |
+----------------+-----------------------------------------------------+
| mnv            | Multiple nucleotide variant.                        |
+----------------+-----------------------------------------------------+
| insertion      | Insertion, length of REF = 1bp.                     |
+----------------+-----------------------------------------------------+
| deletion       | Deletion, length of ALT = 1bp.                      |
+----------------+-----------------------------------------------------+
| complex        | Length of REF and ALT differ and are both \> than 1 |
|                | bp.                                                 |
+----------------+-----------------------------------------------------+
| symbolic       | Structural variant or IUPAC ambiguity code.         |
+----------------+-----------------------------------------------------+

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

+-------------------+--------------------------------------------------+
| Column Name       | Definition                                       |
+===================+==================================================+
| `nchar_ref`       | The length (in bp) of the reference allele.      |
+-------------------+--------------------------------------------------+
| `nchar_alt`       | The length (in bp) of the alternate allele.      |
+-------------------+--------------------------------------------------+
| `varlen`          | The length (in bp) of the variant.               |
+-------------------+--------------------------------------------------+
| `total_depth`     | The total read depth at this position, excluding |
|                   | N-calls.                                         |
+-------------------+--------------------------------------------------+
| `vaf`             | The variant allele fraction. Calculated as       |
|                   | `alt_depth`/`depth_col` where `depth_col` can be |
|                   | `total_depth` or `depth`.                        |
+-------------------+--------------------------------------------------+
| `is_germline`     | TRUE or FALSE. Flags ostensible germline         |
|                   | mutations (`vaf` \> `vaf_cutoff`).               |
+-------------------+--------------------------------------------------+
| `ref_depth`       | The total read depth at the position calling for |
|                   | the reference allele. Calculated as              |
|                   | `depth_col` - `alt_depth` where `depth_col` can  |
|                   | be `total_depth` or `depth`.                     |
+-------------------+--------------------------------------------------+
| `subtype`         | The substitution type for the snv variant        |
|                   | (12-base spectrum; e.g. A\>C)                    |
+-------------------+--------------------------------------------------+
| `short_ref`       | The reference base at this position.             |
+-------------------+--------------------------------------------------+
| `no               | The C/T-based substitution type for the snv      |
| rmalized_subtype` | variant (6-base spectrum; e.g. A\>C -\> T\>G).   |
+-------------------+--------------------------------------------------+
| `normalized_ref`  | The reference base in C/T-base notation for this |
|                   | position (e.g. A -\> T).                         |
+-------------------+--------------------------------------------------+
| `conte            | The substitution type fo the snv variant         |
| xt_with_mutation` | including the two flanking nucleotides           |
|                   | (192-trinucleotide spectrum; e.g. `T[A>C]G`).    |
+-------------------+--------------------------------------------------+
| `normalized_conte | The C/T-based substitution type for the snv      |
| xt_with_mutation` | variant including the two flanking nucleotides   |
|                   | (96-base spectrum e.g. `T[A>C]G` -\> `C[T>G]A`). |
+-------------------+--------------------------------------------------+
| `no               | The trinucleotide context in C/T base notation   |
| rmalized_context` | for this position (e.g. TAG -\> CTA).            |
+-------------------+--------------------------------------------------+
| `gc_content`      | \% GC of the trinucleotide context at this       |
|                   | position.                                        |
+-------------------+--------------------------------------------------+

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

To be filled in more...

## Step 2

## Step 3
