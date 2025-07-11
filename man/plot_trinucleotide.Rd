% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_trinucleotide.R
\name{plot_trinucleotide}
\alias{plot_trinucleotide}
\title{Plot the trinucleotide spectrum}
\usage{
plot_trinucleotide(
  mf_96,
  response = "proportion",
  mf_type = "min",
  group_col = "dose",
  indiv_y = FALSE,
  sum_totals = TRUE,
  output_path = NULL,
  output_type = "svg"
)
}
\arguments{
\item{mf_96}{A data frame containing the mutation frequency data at the
96-base resolution. This should be obtained using the 'calculate_mf' with
subtype_resolution set to 'base_96'. Generally, cols_to_group should be the
same as 'group_col'.}

\item{response}{A character string specifying the type of response to plot.
Must be one of 'frequency', 'proportion', or 'sum'.}

\item{mf_type}{A character string specifying the mutation count method to
plot. Must be one of 'min' or 'max'. Default is 'min'.}

\item{group_col}{A character string specifying the column(s) in 'mf_96'
to group the data by. Default is 'sample'. The sum, proportion, or frequency
will be plotted for all unique levels of this
group. You can specify more than one column to group by. Generally the same
as the 'cols_to_group' parameter in 'calculate_mf' when generating mf_96.}

\item{indiv_y}{A logical value specifying whether the the max response value
for the y-axis should be scaled independently for each group (TRUE) or scaled
the same for all groups (FALSE). Default is FALSE.}

\item{sum_totals}{A logical value specifying whether to display the total
sum of mutations in the mutation labels. Default is TRUE.}

\item{output_path}{An optional file path to an output directory. If provided,
the plots will be automatically exported using the graphics device
specified in output_type. The function will create the output directory if it
doesn't already exist. If NULL, plots will not be exported. Default is NULL.}

\item{output_type}{A character string specifying the type of output file.
Options are  'eps', 'ps', 'tex', 'pdf', or 'jpeg', 'tiff', 'png', 'bmp',
'svg', or 'wmf' (windows only). Default is 'svg'.}
}
\value{
A named list containing ggplots.
}
\description{
Creates barplots of the trinucleotide spectrum for all levels of
a given group.
}
\details{
The function plots the trinucleotide spectrum for all levels of a
given group from the provided mf_96 data; the output of calculate_mf with
subtype_resolution = "base_96".
}
\examples{
# Load example data
example_file <- system.file(
 "extdata", "Example_files",
 "example_mutation_data_filtered.rds",
 package = "MutSeqR"
)
example_data <- readRDS(example_file)

# Calculate the mutation frequency data at the 96-base resolution
mf_96 <- calculate_mf(
 mutation_data = example_data,
 cols_to_group = "dose_group",
 subtype_resolution = "base_96",
 variant_types = "snv"
)
# Plot the trinucleotide proportions for each dose group
# Scale y-axis the same for all groups
plots <- plot_trinucleotide(
 mf_96 = mf_96,
 response = "proportion",
 mf_type = "min",
 group_col = "dose_group",
 indiv_y = FALSE,
 output_path = NULL
)
}
