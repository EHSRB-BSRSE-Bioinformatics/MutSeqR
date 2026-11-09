#' Import a VCF file
#'
#' @description The function reads VCF file(s) and extracts the
#' data into a dataframe.
#' @param vcf_file The path to the .vcf (.gvcf, gzip, bgzip) to be
#' imported. If you specify a directory, the function will
#' attempt to read all files in the directory and combine them into
#' a single table. VCF files should follow the VCF specifications,
#' version 4.5. Multisample VCF files are not supported; VCF files
#' must contain one sample each. Required fields are listed in details.
#' @param sample_data An optional file containing additional sample
#' metadata (dose, timepoint, etc.). This can be a data frame or a file path.
#' Metadata will be joined with the mutation data based on the sample column.
#' Required columns are `sample` and any additional columns you wish to
#' include.
#' @param sd_sep The delimiter for importing sample metadata tables.
#' Default is tab-delimited.
#' @param regions An optional file containing metadata of genomic regions.
#' Region metadata will be joined with mutation data and variants will be
#' checked for overlap with the regions. `regions` can be either a file path,
#' a data frame, or a GRanges object. File paths will be read using the rg_sep.
#' Users can also choose from the built-in TwinStrand's Mutagenesis Panels by
#' inputting "TSpanel_human",  "TSpanel_mouse", or "TSpanel_rat". Required
#' columns for the regions file are "contig", "start", and "end". For a GRanges
#' object, the required columns are "seqnames", "start", and "end". Default is
#' NULL.
#' @param rg_sep The delimiter for importing the custom_regions. The default is
#' tab-delimited "\\t".
#' @param is_0_based_rg A logical variable. Indicates whether the position
#' coordinates in `regions` are 0 based (TRUE) or 1 based (FALSE).
#' If TRUE, positions will be converted to 1-based (start + 1).
#' Need not be supplied for TSpanels. Default is TRUE.
#' @param padding Extend the range of your regions
#' in both directions by the given amount. Ex. Structural variants and
#' indels may start outside of the regions. Adjust the `padding` to
#' include these variants in your region's ranges.
#' @param BS_genome The pkgname of a BS genome. A BS genome must be installed
#' prior to import to populate the context column (trinucleotide context for each position).
#' Only required if data does not already include a context column, such as the
#' INFO field `CONTEXT_3mer` that vartovcf writes with `--context-3mer`. Please install the
#' appropriate BS genome using BiocManager::install("pkgname") where pkgname is the
#' name of the BSgenome package. The pkgname can be found using the find_BS_genome()
#' function, which requires the species and assembly version.
#' Ex. "BSgenome.Hsapiens.UCSC.hg38" | "BSgenome.Hsapiens.UCSC.hg19" |
#' "BSgenome.Mmusculus.UCSC.mm10" | "BSgenome.Mmusculus.UCSC.mm39" |
#' "BSgenome.Rnorvegicus.UCSC.rn6"
#' @param output_granges `TRUE` or `FALSE`; whether you want the mutation
#' data to output as a GRanges object. Default output is as a dataframe.
#' @param remove_sample_suffix An optional character string representing a regular
#' expression to remove unwanted suffixes from VCF sample names prior to joining
#' with metadata. For example, if your VCF sample is "Sample1.cons.filtered" but
#' your metadata sheet just says "Sample1", you can use `remove_sample_suffix = "\\.cons\\.filtered$"`.
#' Default is NULL.
#' @param add_chr A logical variable. If `TRUE`, prepends "chr" to contig names
#' missing it (e.g., "1" becomes "chr1") and changes "MT" to "chrM" to ensure
#' compatibility with BSgenome packages. Default is `FALSE`.
#' @details The required fields are:
#'
#' **FIXED FIELDS**
#' \itemize{
#' \item `CHROM`: The name of the reference sequence. Equivalent to `contig`.
#' \item `POS`: The 1-based start position of the feature. Equivalent to 
#' `start`.
#' \item `REF`: The reference allele at this position.
#' \item `ALT`: The left-aligned, normalized, alternate allele at this position.
#' Multiple alt alleles called for a single position should be represented as
#' separate rows in the table.
#' }
#'
#' **INFO FIELDS**
#' \itemize{
#' \item `sample`: An identifying field for your samples; either in the INFO
#' field or as the header to the FORMAT field.
#' }
#'
#' **SUGGESTED FIELDS**
#'
#' The following **FORMAT** fields are not required, but are recommended for
#' full package functionality:
#' \itemize{
#' \item `AD`: The allelic depths for the reference and alternate allele
#' in the order listed. The sum of AD is equivalent to the `total_depth`
#' (read depth at this position excluding N-calls).
#'  \item `DP`: The read depth at this position (including N-calls).
#' Equivalent to `depth`. Note that in many VCF files, the DP field
#' is defined as `total_depth`. However, in most cases, the DP field
#' includes N-calls.
#'  \item `VD`: The read depth supporting the alternate allele. If
#' not included, the function will take it from the alternate allele depths
#' of `AD` (0 for records without an alternate allele), or assume a value of
#' 1 without `AD`. Equivalent to `alt_depth`.
#' }
#' We recommend that files include a record for every sequenced
#' position, regardless of whether a variant was called, along with the
#' `AD` for each record. This enables site-specific depth calculations
#' required for some downstream analyses. AD is used to calculate the
#' `total_depth` (the read depth excluding No-calls). If AD is not available,
#' or a record's AD is missing its reference depth or the alternate depth of a
#' variant, as callers such as VarDict write for MNV and complex calls, the
#' `DP` field will be used as the `total_depth`.
#'
#' The FORMAT field `SITE_DP`, which vartovcf writes as the read depth at the
#' record's position, is read as the `total_depth`. Records whose
#' `total_depth` is missing take it from AD or DP as above.
#' @returns A table where each row is a mutation, and columns indicate the
#' location, type, and other data. If `output_granges` is set to TRUE, the
#' mutation data will be returned as a GRanges object, otherwise mutation
#' data is returned as a dataframe.
#'
#' Output Column Definitions:
#' \itemize{
#' \item 'end': The half-open end position of the feature. This is calculated
#' as the 1-based start position + the length of the reference allele - 1. For
#' structural variants, if the INFO field contains SVLEN, end is calculated as
#' the 1-based start position + SVLEN - 1. End is only calculated if
#' not already present in the VCF INFO fields.
#' \item `short_ref`: The reference base at the start position.
#' \item `normalized_ref`: The short_ref in C/T-base notation for
#' this position (e.g. A -> T, G -> C).
#' \item `context` The trinucleotide context at this position. Consists
#' of the reference base and the two flanking bases (e.g. TAC).
#' \item `normalized_context`: The trinucleotide context in C/T base
#' notation for this position (e.g. TAG -> CTA).
#'  \item `variation_type` The type of variant (snv, mnv, insertion,
#' deletion, complex, sv, no_variant, ambiguous, uncategorized).
#' \item `subtype` The substitution type for the snv variant (12-base spectrum;
#' e.g. A>C).
#' \item `normalized_subtype` The C/T-based substitution type for the snv
#' variant (6-base spectrum; e.g. A>C -> T>G).
#' \item `context_with_mutation`: The substitution type for the snv variant
#' including the two flanking nucleotides (192-trinucleotide spectrum;
#' e.g. `T[A>C]G`)
#' \item `normalized_context_with_mutation`: The C/T-based substitution
#' type for the snv variant including the two flanking nucleotides
#' (96-base spectrum e.g. `T[A>C]G` -> `C[T>G]A`).
#' \item `nchar_ref`: The length (in bp) of the reference allele.
#' \item `nchar_alt`: The length (in bp) of the alternate allele.
#' \item `varlen`: The length (in bp) of the variant.
#' \item `ref_depth`: The depth of the reference allele. Calculated as
#' `total_depth` - `alt_depth`, if applicable.
#' \item `vaf` : The variant allele fraction. Calculated as
#' `alt_depth`/`total_depth`.
#' \item `gc_content`: % GC of the trinucleotide context at this position.
#' \item `filter`: The VCF FILTER value: "PASS", the semicolon-separated
#' filters the record failed, or "." when no filters were applied. See the
#' `custom_filter_col` parameter of `filter_mut()` to act on it.
#' \item `is_known`: TRUE or FALSE. Flags known variants (ID != ".").
#' \item `row_has_duplicate`: TRUE or FALSE. Flags rows whose position is
#' the same as that of at least one other row for the same sample.
#' \item `filter_mut` : A logical value, initially set to FALSE that indicates
#' to calculte_mf() if the variant should be excluded from mutation counts.
#' See the filter_mut function for more detail.
#' }
#' @examples
#' # Mutation data is just for example purposes. It does not reflect real data
#' file <- system.file("extdata", "Example_files",
#'                    "simple_vcf_data.vcf", package = "MutSeqR")
#' # Import the data
#' imported_example_data <- import_vcf_data(
#'  vcf_file = file,
#' BS_genome = find_BS_genome("mouse", "mm10"))
#' @importFrom  VariantAnnotation alt info geno readVcf ref rbind
#' @importFrom dplyr bind_rows filter group_by left_join mutate rename select summarize ungroup
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_sub str_count
#' @importFrom SummarizedExperiment colData
#' @importFrom plyranges join_overlap_left
#' @importFrom Biostrings getSeq
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom BiocGenerics strand start end
#' @importFrom Seqinfo seqnames
#' @importFrom BSgenome getBSgenome installed.genomes
#' @export
import_vcf_data <- function(
  vcf_file,
  sample_data = NULL,
  sd_sep = "\t",
  regions = NULL,
  rg_sep = "\t",
  is_0_based_rg = FALSE,
  padding = 0,
  BS_genome = NULL,
  output_granges = FALSE,
  remove_sample_suffix = NULL,
  add_chr = FALSE
) {
  stopifnot(
    !missing(vcf_file) && is.character(vcf_file),
    is.null(sample_data) ||
      is.character(sample_data) ||
      is.data.frame(sample_data),
    is.character(sd_sep),
    is.null(remove_sample_suffix) || is.character(remove_sample_suffix),
    is.null(regions) ||
      is.character(regions) ||
      is.data.frame(regions) ||
      methods::is(regions, "GRanges"),
    is.character(rg_sep),
    is.logical(is_0_based_rg),
    is.numeric(padding) && padding >= 0,
    is.logical(output_granges),
    is.logical(add_chr)
  )

  # Load and validate sample metadata before heavy lifting
  sample_df <- NULL
  if (!is.null(sample_data)) {
    sample_df <- import_sample_data(sample_data, sd_sep)
  }

  vcf_file <- file.path(vcf_file)
  if (!file.exists(vcf_file)) {
    stop("The file path you've specified is invalid")
  }

  # Read each file to its own table so sample and fields stay file-local
  if (file.info(vcf_file)$isdir) {
    vcf_files <- list.files(
      vcf_file,
      pattern = "\\.g?vcf(\\.bgz|\\.gz)?$",
      full.names = TRUE
    )
    if (length(vcf_files) == 0) {
      stop("No VCF files found in directory: ", vcf_file)
    }
    tables <- lapply(vcf_files, read_vcf_table, add_chr = add_chr)
    names(tables) <- basename(vcf_files)
    tables <- Filter(Negate(is.null), tables)
    if (length(tables) == 0) {
      stop("No variant records found in the VCF files of: ", vcf_file)
    }
    dat <- dplyr::bind_rows(tables)
  } else {
    dat <- read_vcf_table(vcf_file, add_chr)
    if (is.null(dat)) {
      stop("VCF file contains no variant records: ", vcf_file)
    }
  }
  dat <- rename_columns(dat)

  if (!is.null(sample_df)) {
    dat <- join_sample_metadata(
      dat, sample_df, remove_sample_suffix, label = "VCF data"
    )
  }

  # The alt column can be NA in VCF-derived data
  dat <- annotate_mutations(
    dat, regions, rg_sep, is_0_based_rg, padding, BS_genome,
    required_columns = setdiff(op$base_required_mut_cols, "alt")
  )

  # Depth
  # Add alt_depth column, if it doesn't exist: from the ALT depths of AD,
  # or assumed to be 1 without them.
  if (!"alt_depth" %in% colnames(dat)) {
    alt_depth <- alt_depth_from_ad(
      dat, grep("^AD_", colnames(dat), value = TRUE)
    )
    dat$alt_depth <- if (is.null(alt_depth)) 1 else ifelse(is.na(alt_depth), 1, alt_depth)
  }
  # total_depth: from the data (e.g. vartovcf's SITE_DP), with any gaps
  # filled from depth - no_calls, AD or DP.
  total_depth <- estimate_total_depth(dat)
  if (is.null(total_depth)) {
    warning(
      "Could not find an appropriate depth column. Some package functionality may be limited.\n"
    )
  } else {
    dat$total_depth <- total_depth
  }

  finalize_mutations(dat, output_granges)
}
