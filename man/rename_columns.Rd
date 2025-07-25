% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename_columns.R
\name{rename_columns}
\alias{rename_columns}
\title{Map column names of mutation data to default column names.
A utility function that renames columns of mutation data to default columns names.}
\usage{
rename_columns(data, column_map = op$column)
}
\arguments{
\item{data}{mutation data}

\item{column_map}{a list that maps synonymous column names to their default.}
}
\value{
the mutation data with column names changed to match default.
}
\description{
Map column names of mutation data to default column names.
A utility function that renames columns of mutation data to default columns names.
}
