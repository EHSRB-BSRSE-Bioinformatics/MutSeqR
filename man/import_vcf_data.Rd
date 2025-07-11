% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_vcf_data.R
\name{import_vcf_data}
\alias{import_vcf_data}
\title{Import a VCF file}
\usage{
import_vcf_data(
  vcf_file,
  sample_data = NULL,
  sd_sep = "\\t",
  regions = NULL,
  rg_sep = "\\t",
  is_0_based_rg = FALSE,
  padding = 0,
  genome = NULL,
  species = NULL,
  masked_BS_genome = FALSE,
  output_granges = FALSE
)
}
\arguments{
\item{vcf_file}{The path to the .vcf (.gvcf, gzip, bgzip) to be
imported. If you specify a directory, the function will
attempt to read all files in the directory and combine them into
a single table. VCF files should follow the VCF specifications,
version 4.5. Multisample VCF files are not supported; VCF files
must contain one sample each. Required fields are listed in details.}

\item{sample_data}{An optional file containing additional sample
metadata (dose, timepoint, etc.). This can be a data frame or a file path.
Metadata will be joined with the mutation data based on the sample column.
Required columns are \code{sample} and any additional columns you wish to
include.}

\item{sd_sep}{The delimiter for importing sample metadata tables.
Default is tab-delimited.}

\item{regions}{An optional file containing metadata of genomic regions.
Region metadata will be joined with mutation data and variants will be
checked for overlap with the regions. \code{regions} can be either a file path,
a data frame, or a GRanges object. File paths will be read using the rg_sep.
Users can also choose from the built-in TwinStrand's Mutagenesis Panels by
inputting "TSpanel_human",  "TSpanel_mouse", or "TSpanel_rat". Required
columns for the regions file are "contig", "start", and "end". For a GRanges
object, the required columns are "seqnames", "start", and "end". Default is
NULL.}

\item{rg_sep}{The delimiter for importing the custom_regions. The default is
tab-delimited "\\t".}

\item{is_0_based_rg}{A logical variable. Indicates whether the position
coordinates in \code{regions} are 0 based (TRUE) or 1 based (FALSE).
If TRUE, positions will be converted to 1-based (start + 1).
Need not be supplied for TSpanels. Default is TRUE.}

\item{padding}{Extend the range of your regions
in both directions by the given amount. Ex. Structural variants and
indels may start outside of the regions. Adjust the \code{padding} to
include these variants in your region's ranges.}

\item{genome}{The genome assembly version of the reference genome. This is
required if your data does not include a context column. The
function will install a BS genome for the given species/genome/masked
arguments to populate the context column.
Ex.Human GRCh38 = hg38 | Human GRCh37 = hg19 | Mouse GRCm38 = mm10 |
Mouse GRCm39 = mm39 | Rat RGSC 6.0 = rn6 | Rat mRatBN7.2 = rn7}

\item{species}{The species. Required if your data does not include a
context column. The function will install a BS genome for the given
species/genome/masked to populate the context column. The species can
be the common name of the species or the scientific name.
Ex. "human" or "Homo sapiens".}

\item{masked_BS_genome}{A logical value. Required when using a BS genome
to poulate the context column. Whether to use the masked version of the
BS genome (TRUE) or not (FALSE). Default is FALSE.}

\item{output_granges}{\code{TRUE} or \code{FALSE}; whether you want the mutation
data to output as a GRanges object. Default output is as a dataframe.}
}
\value{
A table where each row is a mutation, and columns indicate the
location, type, and other data. If \code{output_granges} is set to TRUE, the
mutation data will be returned as a GRanges object, otherwise mutation
data is returned as a dataframe.

Output Column Definitions:
\itemize{
\item \code{short_ref}: The reference base at the start position.
\item \code{normalized_ref}: The short_ref in C/T-base notation for
this position (e.g. A -> T, G -> C).
\item \code{context} The trinucleotide context at this position. Consists
of the reference base and the two flanking bases (e.g. TAC).
\item \code{normalized_context}: The trinucleotide context in C/T base
notation for this position (e.g. TAG -> CTA).
\item \code{variation_type} The type of variant (snv, mnv, insertion,
deletion, complex, sv, no_variant, ambiguous, uncategorized).
\item \code{subtype} The substitution type for the snv variant (12-base spectrum;
e.g. A>C).
\item \code{normalized_subtype} The C/T-based substitution type for the snv
variant (6-base spectrum; e.g. A>C -> T>G).
\item \code{context_with_mutation}: The substitution type for the snv variant
including the two flanking nucleotides (192-trinucleotide spectrum;
e.g. \verb{T[A>C]G})
\item \code{normalized_context_with_mutation}: The C/T-based substitution
type for the snv variant including the two flanking nucleotides
(96-base spectrum e.g. \verb{T[A>C]G} -> \verb{C[T>G]A}).
\item \code{nchar_ref}: The length (in bp) of the reference allele.
\item \code{nchar_alt}: The length (in bp) of the alternate allele.
\item \code{varlen}: The length (in bp) of the variant.
\item \code{ref_depth}: The depth of the reference allele. Calculated as
\code{total_depth} - \code{alt_depth}, if applicable.
\item \code{vaf} : The variant allele fraction. Calculated as
\code{alt_depth}/\code{total_depth}.
\item \code{gc_content}: \% GC of the trinucleotide context at this position.
\item \code{is_known}: TRUE or FALSE. Flags known variants (ID != ".").
\item \code{row_has_duplicate}: TRUE or FALSE. Flags rows whose position is
the same as that of at least one other row for the same sample.
\item \code{filter_mut} : A logical value, initially set to FALSE that indicates
to calculte_mf() if the variant should be excluded from mutation counts.
See the filter_mut function for more detail.
}
}
\description{
The function reads VCF file(s) and extracts the
data into a dataframe.
}
\details{
The required fields are:

\strong{FIXED FIELDS}
\itemize{
\item \code{CHROM}: The name of the reference sequence. Equivalent to \code{contig}.
\item \code{POS}: The 1-based start position of the feature. Equivalent to  \code{start}.
\item \code{REF}: The reference allele at this position.
\item \code{ALT}: The left-aligned, normalized, alternate allele at this position.
Multiple alt alleles called for a single position should be represented as
separate rows in the table.
}

\strong{INFO FIELDS}
\itemize{
\item \code{END}: The half-open end position of the feature.
\item \code{sample}: An identifying field for your samples; either in the INFO
field or as the header to the FORMAT field.
}

\strong{SUGGESTED FIELDS}

The following \strong{FORMAT} fields are not required, but are recommended for
full package functionality:
\itemize{
\item \code{AD}: The allelic depths for the reference and alternate allele
in the order listed. The sum of AD is equivalent to the \code{total_depth}
(read depth at this position excluding N-calls).
\item \code{DP}: The read depth at this position (including N-calls).
Equivalent to \code{depth}. Note that in many VCF files, the DP field
is defined as \code{total_depth}. However, in most cases, the DP field
includes N-calls.
\item \code{VD}: The read depth supporting the alternate allele. If
not included, the function will add this column, assuming a value of 1.
Equivalent to \code{alt_depth}.
}
We recommend that files include a record for every sequenced
position, regardless of whether a variant was called, along with the
\code{AD} for each record. This enables site-specific depth calculations
required for some downstream analyses. AD is used to calculate the
\code{total_depth} (the read depth excluding No-calls). If AD is not available,
the \code{DP} field will be used as the \code{total_depth}.
}
\examples{
# Example: Import a single bg-zipped vcf file. This library was sequenced
# with Duplex Sequencing using the TwinStrand Mouse Mutagenesis Panel which
# consists of 20 2.4kb targets = 48kb of sequence.
example_file <- system.file("extdata", "Example_files",
                            "example_import_vcf_data_cleaned.vcf.bgz",
                            package = "MutSeqR")
# We will create an example metadata table for this data.
sample_meta <- data.frame(sample = "dna00996.1",
                          dose = "50",
                          dose_group = "High")
# Import the data
imported_example_data <- import_vcf_data(vcf_file = example_file,
                                         sample_data = sample_meta,
                                         regions = "TSpanel_mouse",
                                         genome = "mm10",
                                         species = "mouse",
                                         masked_BS_genome = FALSE)
}
