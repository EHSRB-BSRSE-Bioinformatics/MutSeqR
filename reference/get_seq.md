# Get sequence of genomic target regions

Create a GRanges object from the genomic target ranges and import raw
nucleotide sequences.

## Usage

``` r
get_seq(
  regions,
  rg_sep = "\t",
  is_0_based_rg = TRUE,
  padding = 0,
  BS_genome = NULL,
  ucsc = FALSE,
  species = NULL,
  genome = NULL
)
```

## Arguments

- regions:

  The regions metadata file to import. Can be either a file path, a data
  frame, or a GRanges object. File paths will be read using the rg_sep.
  Users can also choose from the built-in TwinStrand's Mutagenesis
  Panels by inputting "TSpanel_human", "TSpanel_mouse", or
  "TSpanel_rat". Required columns for the regions file are "contig",
  "start", and "end". In a GRanges object, the required columns are
  "seqnames", "start", and "end".

- rg_sep:

  The delimiter for importing the regions file. The default is
  tab-delimited ("\t").

- is_0_based_rg:

  A logical variable. Indicates whether the position coordinates in
  `regions` are 0 based (TRUE) or 1 based (FALSE). If TRUE, positions
  will be converted to 1-based (start + 1). Need not be supplied for
  TSpanels. Default is TRUE.

- padding:

  An integer value by which the function will extend the range of the
  target sequence on both sides. Start and end coordinates will be
  adjusted accordingly. Default is 0.

- BS_genome:

  The name of the appropriate BSgenome package to use for sequence
  retrieval. Ex. "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Mmusculus.UCSC.mm10", "BSgenome.Rnorvegicus.UCSC.rn6". Use
  the function find_BS_genome() to help identify the appropriate
  BSgenome package if needed. Need not be supplied for TSpanels.
  BS_genome must be installed if using this method.

- ucsc:

  A logical value. If TRUE, the function will retrieve the sequences
  from the UCSC genome browser using an API. If FALSE, the function will
  retrieve sequences using the appropriate BSgenome package, which will
  be installed as needed. Default is FALSE.

- species:

  The species for which to retrieve the sequences. Only required if
  using the UCSC method. Species may be given as the scientific name or
  the common name. Ex. "Human", "Homo sapien". Used to choose the
  appropriate BS genome. Need not be supplied for TSpanels.

- genome:

  The genome assembly version for which to retrieve the sequences. Only
  required if using the UCSC method. Ex. hg38, hg19, mm10, mm39, rn6,
  rn7. Need not be supplied for TSpanels.

## Value

a GRanges object with sequences of targeted regions.

## Details

Consult
`available.genomes(splitNameParts=FALSE, type=getOption("pkgType"))` for
a full list of the available BS genomes and their associated
species/genome/masked values. The BSgenome package will be installed if
not already available. If using the UCSC API, the function will retrieve
the sequences from the UCSC genome browser using the DAS API. See the
UCSC website for available genomes: <https://genome.ucsc.edu>.

## Examples

``` r
#  Retrieve the sequences for custom regions
# We will load the TSpanel_human regions file as an example
# and supply it to the function as a GRanges object.
human <- load_regions_file("TSpanel_human")
regions_seq <- get_seq(
  regions = human,
  is_0_based_rg = FALSE,
  BS_genome = "BSgenome.Hsapiens.UCSC.hg38",
  padding = 0
)
#> Loading reference genome: BSgenome.Hsapiens.UCSC.hg38.
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> 
#> Attaching package: ‘Biostrings’
#> The following object is masked from ‘package:base’:
#> 
#>     strsplit
```
