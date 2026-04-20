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
#' Only required if data does not already include a context column. Please install the
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
#' @details The required fields are:
#'
#' **FIXED FIELDS**
#' \itemize{
#' \item `CHROM`: The name of the reference sequence. Equivalent to `contig`.
#' \item `POS`: The 1-based start position of the feature. Equivalent to  `start`.
#' \item `REF`: The reference allele at this position.
#' \item `ALT`: The left-aligned, normalized, alternate allele at this position.
#' Multiple alt alleles called for a single position should be represented as
#' separate rows in the table.
#' }
#'
#' **INFO FIELDS**
#' \itemize{
#' \item `END`: The half-open end position of the feature.
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
#' not included, the function will add this column, assuming a value of 1.
#' Equivalent to `alt_depth`.
#' }
#' We recommend that files include a record for every sequenced
#' position, regardless of whether a variant was called, along with the
#' `AD` for each record. This enables site-specific depth calculations
#' required for some downstream analyses. AD is used to calculate the
#' `total_depth` (the read depth excluding No-calls). If AD is not available,
#' the `DP` field will be used as the `total_depth`.
#' @returns A table where each row is a mutation, and columns indicate the
#' location, type, and other data. If `output_granges` is set to TRUE, the
#' mutation data will be returned as a GRanges object, otherwise mutation
#' data is returned as a dataframe.
#'
#' Output Column Definitions:
#' \itemize{
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
#' @importFrom dplyr filter group_by left_join mutate rename select summarize ungroup
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
  remove_sample_suffix = NULL
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
    is.logical(output_granges)
  )
  BS_genome <- match.arg(
    BS_genome,
    choices = c(
      NULL,
      BSgenome::available.genomes(splitNameParts = TRUE)$pkgname
    )
  )

  # Load and validate sample metadata before heavy lifting
  sample_df <- NULL
  if (!is.null(sample_data)) {
    sample_df <- import_sample_data(sample_data, sd_sep)
  }

  vcf_file <- file.path(vcf_file)

  # Read and bind vcfs from folder
  if (file.info(vcf_file)$isdir == TRUE) {
    vcf_files <- list.files(
      vcf_file,
      pattern = "\\.g?vcf(\\.bgz|\\.gz)?$",
      full.names = TRUE
    )
    if (length(vcf_files) == 0) {
      stop("No VCF files found in directory: ", vcf_file)
    }
    # Read and combine VCF files
    vcf_list <- lapply(vcf_files, function(file) {
      vcf <- VariantAnnotation::readVcf(file)
      vcf <- vcf_sample_fix(vcf) # fix sample column
      # Ensure consistent colData rownames so rbind doesn't complain

      # Coerce ALT to CharacterList for every VCF.
      # This prevents rbind crashes when mixing strict DNA alleles (e.g. A, C)
      # with gVCF/SV alleles (e.g. <NON_REF>, .)
      VariantAnnotation::alt(
        vcf
      ) <- IRanges::CharacterList(VariantAnnotation::alt(vcf))

      rownames(SummarizedExperiment::colData(vcf)) <- "sample_info"
      return(vcf)
    })
    vcf <- do.call(VariantAnnotation::rbind, vcf_list)
  } else {
    # Read a single vcf file
    vcf <- VariantAnnotation::readVcf(vcf_file)
    # Rename or create the "sample" column in the INFO field
    vcf <- vcf_sample_fix(vcf)
    VariantAnnotation::alt(
      vcf
    ) <- IRanges::CharacterList(VariantAnnotation::alt(vcf))
  }
  # Extract mutation data into a dataframe
  ## To Do: May want to use the expand function to unlist ALT column of a CollapsedVCF object to one row per ALT value.
  dat <- data.frame(
    contig = SummarizedExperiment::seqnames(vcf),
    start = SummarizedExperiment::start(vcf),
    end = SummarizedExperiment::end(vcf),
    ref = VariantAnnotation::ref(vcf),
    alt = VariantAnnotation::alt(vcf)
  )

  # Extract INFO fields
  info <- as.data.frame(VariantAnnotation::info(vcf))

  # Extract GENO fields depending on the type of data
  geno <- VariantAnnotation::geno(vcf)
  geno_df <- data.frame(row.names = seq_len(nrow(geno[[1]])))
  for (field_name in names(geno)) {
    field <- geno[[field_name]]
    if (is.list(field)) {
      # Ex. AD
      max_length <- max(vapply(field, length, integer(1)))
      expanded_field <- do.call(
        rbind,
        lapply(field, function(x) {
          c(x, rep(NA, max_length - length(x)))
        })
      )
      colnames(expanded_field) <- paste(
        field_name,
        seq_len(max_length),
        sep = "_"
      )
      geno_df <- cbind(geno_df, expanded_field)
    } else if (is.matrix(field)) {
      # Ex. GT, DP, VD
      geno_df[[field_name]] <- as.vector(field)
    } else if (is.array(field) && length(dim(field)) == 3) {
      # Ex. RD, ALD
      # Collapse the array over the 2nd and 3rd dimensions
      collapsed_field <- apply(field, c(1), function(x) as.vector(x))
      collapsed_field <- as.data.frame(t(collapsed_field))
      ncols <- ncol(collapsed_field)
      colnames(collapsed_field) <- paste(field_name, seq_len(ncols), sep = "_")
      geno_df <- cbind(geno_df, collapsed_field)
    } else {
      geno_df[[field_name]] <- field
    }
  }
  # Ensure info does not overwrite columns already in dat or geno_df (like END or DP)
  # We do a case-insensitive match so "END" in info is safely dropped in favor of "end" in dat
  common_cols_idx <- tolower(colnames(info)) %in%
    tolower(c(colnames(dat), colnames(geno_df)))
  info <- info[, !common_cols_idx, drop = FALSE]

  # Combine data frames
  dat <- cbind(dat, geno_df, info)
  row.names(dat) <- NULL

  prepared <- prepare_imported_mutation_data(
    dat = dat,
    sample_df = sample_df,
    source_label = "VCF data",
    mismatch_hint = "Please check for suffixes (e.g. '.cons.filtered') in your VCF files or typos in your metadata file.",
    remove_sample_suffix = remove_sample_suffix,
    allow_na_columns = "alt",
    BS_genome = BS_genome
  )
  dat <- prepared$dat
  context_exists <- prepared$context_exists

  mut_ranges <- build_imported_mutation_ranges(
    dat = dat,
    context_exists = context_exists,
    regions = regions,
    rg_sep = rg_sep,
    is_0_based_rg = is_0_based_rg,
    padding = padding,
    BS_genome = BS_genome
  )

  resolve_vcf_depth <- function(dat) {
    total_depth_exists <- "total_depth" %in% colnames(dat)
    depth_exists <- "depth" %in% colnames(dat)
    no_calls_exists <- "no_calls" %in% colnames(dat)
    ad_columns <- grep("^AD_", colnames(dat), value = TRUE)

    if (!total_depth_exists) {
      if (no_calls_exists && depth_exists) {
        dat <- dat %>%
          dplyr::mutate(total_depth = .data$depth - .data$no_calls)
      } else if (length(ad_columns) > 0) {
        dat$total_depth <- rowSums(dat[, ad_columns], na.rm = TRUE)
      } else if (depth_exists) {
        dat <- dat %>%
          dplyr::mutate(total_depth = .data$depth)
        warning(
          "Could not find total_depth column and cannot calculate. The 'total_depth' will be set to DP. You can review the diffference in the README"
        )
      } else {
        warning(
          "Could not find an appropriate depth column. Some package functionality may be limited.\n"
        )
      }
    }

    dat
  }

  finalize_imported_mutation_data(
    mut_ranges = mut_ranges,
    depth_resolver = resolve_vcf_depth,
    output_granges = output_granges
  )
}
