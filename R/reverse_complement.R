#' Get the reverse complement of a DNA or RNA sequence.
#' 
#' @param x A character vector of DNA or RNA sequences.
#' @param content c("dna", "rna") The type of sequence to be reversed.
#' @param case c("lower", "upper", "as is") The case of the output sequence.
#' @details This file is part of the source code for
#' SPGS: an R package for identifying statistical patterns in genomic sequences.
#' Copyright (C) 2015  Universidad de Chile and INRIA-Chile
#
#' A copy of Version 2 of the GNU Public License is available in the
#' share/licenses/gpl-2 file in the R installation directory or from
#' http://www.R-project.org/Licenses/GPL-2.
#' reverseComplement.R
#' @export
reverseComplement <- function(x, 
                              content=c("dna", "rna"), 
                              case=c("lower", "upper", "as is")) {
 
  # reverse character vector
  strreverse <- function(x) {
  if (!is.character(x))
    stop("x must be a character vector")
  sapply(strsplit(x, ""), function(y) paste(rev(y), collapse = ""))
} 
 #Check arguments
if (!is.character(x)) x <- as.character(x) #coerse x to a character vector
  content <- match.arg(content)
  case <- match.arg(case)
  if (length(x)==0 || (length(x)==1 && nchar(x)==0))
    return(x) #bail if input is empty
  if (case=="lower") x <- tolower(x)
  if (case=="upper") x <- toupper(x)
  if (content=="dna")
  {
    src  <- "acgturykmswbdhvnxACGTURYKMSWBDHVNX-"
    dest <- "tgcaayrmkswvhdbnxTGCAAyRMKSWVHDBNX-"
  }
  else
  {
    src  <- "acgturykmswbdhvnxACGTURYKMSWBDHVNX-"
    dest <- "ugcaayrmkswvhdbnxUGCAAyRMKSWVHDBNX-"
  } #if
  if (max(nchar(x))>1) #is x a vector of one or more non-single character strings?
    return(chartr(src, dest, strreverse(x)))
  # x is not a single string, so process it as a vector
  chartr(src, dest, rev(x))
} #function
