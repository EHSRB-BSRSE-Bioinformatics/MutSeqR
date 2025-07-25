% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_mut_data.R
\name{import_mut_data}
\alias{import_mut_data}
\title{Import tabular mutation data}
\usage{
import_mut_data(
  mut_file,
  mut_sep = "\\t",
  is_0_based_mut = TRUE,
  sample_data = NULL,
  sd_sep = "\\t",
  regions = NULL,
  rg_sep = "\\t",
  is_0_based_rg = TRUE,
  padding = 0,
  genome = NULL,
  species = NULL,
  masked_BS_genome = FALSE,
  custom_column_names = NULL,
  output_granges = FALSE
)
}
\arguments{
\item{mut_file}{The mutation data file(s) to be imported.
This can be either a data frame object or a filepath
to a file or directory. If you specify a directory, the function will
attempt to read all files in the directory and combine them into
a single data frame. Mutation data should consist of a row for each
variant. Required columns are listed in details.}

\item{mut_sep}{The delimiter for importing the mutation file.
Default is tab-delimited.}

\item{is_0_based_mut}{A logical variable. Indicates whether the
position coordinates in the mutation data are 0 based (TRUE) or
1 based (FALSE). If TRUE, positions will be converted to 1-based.}

\item{sample_data}{An optional file containing additional sample
metadata (dose, timepoint, etc.). This can be a data frame or a file path.
Metadata will be joined with the mutation data based on the sample column.
Required columns are \code{sample} and any additional columns you wish to
include.}

\item{sd_sep}{The delimiter for importing sample data.
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

\item{padding}{An integer >= 0. Extend the range of your regions
in both directions by the given amount. Ex. Structural variants and
indels may start outside of the regions. Adjust the padding to
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

\item{custom_column_names}{A list of names to specify the meaning of column
headers. Since column names can vary with data, this might be necessary to
digest the mutation data properly. Typical defaults are set, but can
be substituted in the form of \code{list(my_custom_contig_name = "contig", my_custom_sample_column_name = "sample")}. You can change one or more of
these. Set column synonyms are defined in MutSeqR::op$column and will
automatically be changed to their default value.}

\item{output_granges}{A logical variable; whether you want the mutation
data to output as a GRanges object. Default output (FALSE) is as a dataframe.}
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
Imports tabular mutation file into the local R environment.
}
\details{
Required columns for mut files are:
\itemize{
\item \code{contig}: The name of the reference sequence.
\item \code{start}: The start position of the feature.
\item \code{end}: The half-open end position of the feature.
\item \code{sample}: The sample name.
\item \code{ref}: The reference allele at this position
\item \code{alt}: The left-aligned, normalized, alternate allele at this
position. Multiple alt alleles called for a single position should be
represented as separate rows in the table.
}
The following columns are not required, but are recommended for full
package functionality:
\itemize{
\item \code{alt_depth}: The read depth supporting the alternate allele. If
not included, the function will add this column, assuming a value of 1.
\item \code{total_depth}: The total read depth at this position, excluding
no-calls (N calls). If not present, the function will attempt to calculate
the \code{total_depth} as \code{depth} - \code{no_calls}. If no_calls is not present, the
function will use \code{depth} as the \code{total_depth.}
\item \code{depth}: The total read depth at this position, including no-calls.
\item \code{no_calls}: The number of no-calls (N-calls) at this position.
}
We recommend that files include a record for every sequenced
position, regardless of whether a variant was called, along with the
\code{total_depth} for each record. This enables site-specific depth calculations
required for some downstream analyses.
}
\examples{
# Example: Import a single mutation file. This library was sequenced with
# Duplex Sequencing using the TwinStrand Mouse Mutagenesis Panel which
# consists of 20 2.4kb targets = 48kb of sequence.
example_file <- system.file("extdata", "Example_files",
                            "example_import_mut_data.rds",
                            package = "MutSeqR")
example_data <- readRDS(example_file)
# We will create an example metadata table for this data.
sample_meta <- data.frame(sample = "dna00996.1",
                          dose = "50",
                          dose_group = "High")
# Import the data
imported_example_data <- import_mut_data(mut_file = example_data,
                                         sample_data = sample_meta,
                                         regions = "TSpanel_mouse",
                                         genome = "mm10",
                                         species = "mouse",
                                         masked_BS_genome = FALSE)
}
