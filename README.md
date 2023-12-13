# Duplex Sequencing Analysis

## Pre-abmle

To do: write a bit about the project, the scope, why it's necessary.

Provide a little background on duplex sequencing, describe the technology
briefly, and give some context for the type of data we are meant to be
processing with this package.

## Data import

The first step is importing data. Your variant file can be imported as either
a `.mut` file using the function `import_mut_data` or as a `.vcf` file using the
function `read_vcf`.

Required columns for your `.mut` file are listed in the table below.
|---------------- |------------------------------------------------------------|
| **Column Name** | **Definition**                                             |
|---------------- |------------------------------------------------------------|
| `contig`        | The reference sequence name.                               |
| `start`         | The 0-based start position of the feature in contig.       |
| `end`           | The half-open end position of the feature in contig.       |
| `sample`        | The sample name.                                           |
| `ref`           | The reference allele at this position.                     |
| `alt`           | The left-aligned, alternate allele at this position.       | 
| `alt_depth`     | The read depth supporting the alternate allele.            |
| depth col       | The total read depth at this position. This column can be  |
|                 |`total_depth` (excluding N-calls) or `depth` (including     |
|                 | N-calls; if `total_depth` is not available).               |
| `variation_type`| The category to which this variant is assigned.            |   
| `context`       | The local reference trinucleotide context at this position | 
|                 | (e.g. ATC - not necessarily the transcript codon).         |
|-----------------|------------------------------------------------------------|

Required fields for your `.vcf` file are listed in the table below.
|--------------|---------------------------------------------------------------|
|**Field Name**| **Definition**                                                |
|--------------|---------------------------------------------------------------|
| **FIXED FIELDS**                                                            ||
|--------------|---------------------------------------------------------------|
| `contig`     | The reference sequence name.                                  |
| `start`      | The 0-based start position of the feature in contig.          |
| `ref`        | The reference allele at this position.                        |
| `alt`        | The left-aligned, alternate allele at this position.          | 
|------------- |---------------------------------------------------------------|
| **FORMAT FIELDS**                                                           ||
|--------------|---------------------------------------------------------------|
| `AD`         | The allelic depths for the reference and alternate alleles in |
|              | the order listed.                                             |
| `DP`         | The total read depth at this position (excluding N-calls).    |
|              | Equivalent to `depth`.                                        |
| `VD`         | Variant Depth. Equivalent to `alt_depth`.                     |
|--------------|---------------------------------------------------------------|
| **INFO FIELDS**                                                             ||
|--------------|---------------------------------------------------------------|
| `TYPE`       | The category to which this variant is assigned. Equivalent to |
|              |`variation_type`.                                              |
| `END`        | The half-open end position of the feature in contig.          |
|--------------|---------------------------------------------------------------|
| **SUGGESTED INFO FIELDS**                                                   ||
|--------------|---------------------------------------------------------------|
| sample       | An identifying field for your samples; either in the INFO     |
|              | field or as the header to the FORMAT field.                   | 
|`SVTYPE`      | Structural variant types; INV DUP DEL INS FUS                 |
| SVLEN`       | Length of the structural variant in base pairs                |
|--------------|---------------------------------------------------------------|

We recognize that column names may differ. Therefore, we have implemented some 
default column name synonyms. If your column name matches one of our listed synonyms,
it will automatically be changed to match our set values. For example, your `contig`
column may be named `chr` or `chromosome`. After your importing your data, this column
name will be changed to `contig`. A full list of column name synonyms are listed below:
(ENTER OP$COLUMN NAMES)

The column variation_type/TYPE may contain these values:
|------------------|------------------------------------------------------|
| `variation_type` | Definition                                           |
|------------------|------------------------------------------------------|
| no_variant       | No variation, the null-case                          |
| snv              | Single nucleotide variant                            |
| mnv              | Multiple nucleotide variant                          |
| insertion        | Insertion, length of REF = 1bp                       |
| deletion         | Deletion, length of ALT = 1bp                        |
| complex          | Length of REF and ALT differ and are both > than 1 bp|
| symbolic         | Structural variant or IUPAC ambiguity code           |
|------------------|------------------------------------------------------|

Since the main goal of this 
package is to generate summary statistics, visualization, exploratory analysis, 
and other post-processing tasks such as mutational signature analysis or 
generalized linear modeling, the main piece of information you want to import is
the `.mut` file, the schema for which is described here:

(example mut file with columns and their descriptions?)
It is critical that this .mut file gets annotated with information such as the
genomic region from which the mutation originates, some information about the 
sample from which the mutation orignates, and that each row (mutation) receives
a calculated group depth and group frequency (of interest for later analyses).

We do this by first importing the variant file as a data frame, and then joining 
it with that other data (e.g., sample_data_file, regions file), and creating some 
columns that will be helpful for calculating frequencies in later analyses.
Finally, our functions provide the option to convert the resulting data frame into a 
`granges` object. This facilitates use in other packages and makes doing 'genome
 math' on the ranges significantly easier. Columns that are added to the resulting
 data frame are listed below. 
 - `nchar_ref`: The length (in bp) of the reference allele.
 - `nchar_alt`: The length (in bp) of the alternate allele.
 - `varlen`: The length (in bp) of the variant.
 - `total_depth`: The total read depth at this position, excluding N-calls.
 - `vaf`: The variant allele fraction. Calculated as `alt_depth`/`depth_col` 
 where `depth_col` can be `total_depth` or `depth`.
 - `is_germline`: TRUE or FALSE. Flags ostensible germline mutations (`vaf` > `vaf_cutoff`).
 - `ref_depth`: The total read depth at the position calling for the reference allele. 
 Calculated as `depth_col` - `alt_depth` where `depth_col` can be `total_depth` or `depth`.
 - `subtype`: The substitution type for the snv variant (12-base spectrum; e.g. A>C)
 - `short_ref`: The reference base at this position.
 - `normalized_subtype`: The C/T-based substitution type for the snv variant 
 (6-base spectrum; e.g. A>C -> T>G)
 - `normalized_ref`: The reference base in C/T-base notation for this position 
 (e.g. A -> T).
 - `context_with_mutation`: The substitution type fo the snv variant including 
 the two flanking nucleotides (192-trinucleotide spectrum; e.g. `T[A>C]G`)
 - `normalized_context_with_mutation`: The C/T-based substitution type for the 
 snv variant including the two flanking nucleotides (96-base spectrum e.g. `T[A>C]G` -> `C[T>G]A`)
 - `normalized_context`: The trinucleotide context in C/T base notation for this 
 position (e.g. TAG -> CTA).
 - `gc_content`: % GC of the trinucleotide context at this position. 

### Metadata: an important consideration

The other important component of importing your data for proper use is to assign
 each mutation to a biological sample, and also make sure that some additional 
 information about each sample is present (e.g., a chemical treatment, a dose, 
 etc.). This is done by providing a sample data file (tab delimited, comma 
 delimited, etc.; the choice is up to the user, but the delimiter of the file 
 must be specified as a parameter in the function). Importantly, this is a file
 that would be analogous to "colData", or "column data", a term often used in 
 the `DESeq2` package. Hence, it must contain some information about an existing
 column in your variant file, which is typically going to be sample. So the first
 column in your sample data file should indeed be `sample`. Then, additional 
 columns such as `dose` or `tissue` or `treatment` can be added, and these 
 columns will be joined with your variant file to capture that information and
 associate it with each mutation.

### Other Notes

When first importing a `.mut` file, it is preferred to keep non-variant rows.
This allows the calculuation of mutation frequencies. The data set can be pared
down later to include only mutations of interest (SNVs, indels, SVs, or any 
combination).

To be filled in more...

## Step 2

## Step 3