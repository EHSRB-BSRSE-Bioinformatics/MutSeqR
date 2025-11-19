# Imports the regions file

A helper function to import the regions metadata file and return a
GRanges object.

## Usage

``` r
load_regions_file(regions, rg_sep = "\t", is_0_based_rg = TRUE)
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

  The delimiter for importing the custom_regions. The default is
  tab-delimited "\t".

- is_0_based_rg:

  A logical variable. Indicates whether the position coordinates in
  `regions` are 0 based (TRUE) or 1 based (FALSE). If TRUE, positions
  will be converted to 1-based (start + 1). Need not be supplied for
  TSpanels. Default is TRUE.

## Value

a GRanges object of the imported regions metadata file.

## Examples

``` r
#' # Example 1: Load built-in TwinStrand's Human Mutagenesis
human_rg <- load_regions_file(regions = "TSpanel_human")
human_rg
#> GRanges object with 20 ranges and 6 metadata columns:
#>        seqnames              ranges strand | target_size description
#>           <Rle>           <IRanges>  <Rle> |   <integer> <character>
#>    [1]    chr11 108510788-108513187      * |        2400 region_1111
#>    [2]    chr13   75803913-75806312      * |        2400 region_1501
#>    [3]    chr14   74661756-74664155      * |        2400 region_1725
#>    [4]    chr18     5749265-5751664      * |        2400 region_2457
#>    [5]     chr2   40162768-40165167      * |        2400 region_2896
#>    ...      ...                 ...    ... .         ...         ...
#>   [16]    chr15   46089738-46092137      * |        2400 region_1904
#>   [17]    chr17   70672727-70675126      * |        2400 region_2378
#>   [18]    chr21   23665977-23668376      * |        2400 region_3515
#>   [19]    chr22   48262371-48264770      * |        2400 region_3703
#>   [20]    chr10 128969038-128971437      * |        2400  region_784
#>        genic_context        gene      genome       label
#>          <character> <character> <character> <character>
#>    [1]         genic       EXPH5        hg38       chr11
#>    [2]         genic        LMO7        hg38       chr13
#>    [3]         genic       AREL1        hg38       chr14
#>    [4]         genic   MIR3976HG        hg38       chr18
#>    [5]         genic  SLC8A1-AS1        hg38        chr2
#>    ...           ...         ...         ...         ...
#>   [16]    intergenic        <NA>        hg38       chr15
#>   [17]    intergenic        <NA>        hg38       chr17
#>   [18]    intergenic        <NA>        hg38       chr21
#>   [19]    intergenic        <NA>        hg38       chr22
#>   [20]    intergenic        <NA>        hg38       chr10
#>   -------
#>   seqinfo: 20 sequences from an unspecified genome; no seqlengths
# Load a custom regions file from an interval list
# We will use the human TSpanel system file for this example,
# but any file can be imported.
file <- system.file("extdata",
  "inputs",
  "metadata",
  "human_mutagenesis_panel_hg38.txt",
  package = "MutSeqR"
)
custom_rg <- load_regions_file(regions = file, rg_sep = "\t", is_0_based_rg = TRUE)
custom_rg
#> GRanges object with 20 ranges and 6 metadata columns:
#>        seqnames              ranges strand | target_size description
#>           <Rle>           <IRanges>  <Rle> |   <integer> <character>
#>    [1]    chr11 108510788-108513187      * |        2400 region_1111
#>    [2]    chr13   75803913-75806312      * |        2400 region_1501
#>    [3]    chr14   74661756-74664155      * |        2400 region_1725
#>    [4]    chr18     5749265-5751664      * |        2400 region_2457
#>    [5]     chr2   40162768-40165167      * |        2400 region_2896
#>    ...      ...                 ...    ... .         ...         ...
#>   [16]    chr15   46089738-46092137      * |        2400 region_1904
#>   [17]    chr17   70672727-70675126      * |        2400 region_2378
#>   [18]    chr21   23665977-23668376      * |        2400 region_3515
#>   [19]    chr22   48262371-48264770      * |        2400 region_3703
#>   [20]    chr10 128969038-128971437      * |        2400  region_784
#>        genic_context        gene      genome       label
#>          <character> <character> <character> <character>
#>    [1]         genic       EXPH5        hg38       chr11
#>    [2]         genic        LMO7        hg38       chr13
#>    [3]         genic       AREL1        hg38       chr14
#>    [4]         genic   MIR3976HG        hg38       chr18
#>    [5]         genic  SLC8A1-AS1        hg38        chr2
#>    ...           ...         ...         ...         ...
#>   [16]    intergenic        <NA>        hg38       chr15
#>   [17]    intergenic        <NA>        hg38       chr17
#>   [18]    intergenic        <NA>        hg38       chr21
#>   [19]    intergenic        <NA>        hg38       chr22
#>   [20]    intergenic        <NA>        hg38       chr10
#>   -------
#>   seqinfo: 20 sequences from an unspecified genome; no seqlengths
```
