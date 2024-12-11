# Requirements for vcf files:
# Follow specifications for GATK GVCF files at BP_RESOLUTION.
# GVCF files at BP_RESOLUTION: a GVCF with an individual record at every site: either a variant record, or a non-variant record.
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
# <NON_REF> is indicated for all fields.
# Requires DP & VD or AD
# END tag in INFO field.
# Defining variants in VCF files:
# TS has alt be "." if there is no variant,
# but some examples show alt = ref for no variant
# TS shows CT > C as a deletion, but other examples
# show T > - as a deletion.
# REQUIREMENTS FOR VCF: Normalized left-aligned and parsimonious.
# https://genome.sph.umich.edu/wiki/Variant_Normalization
# Suggest: to improve efficiency, index your file, ensure that the header is complete. bgzip

#' import genomic vcf files
#' @description import large genomic vcf files into the R env. as mutation data
#' @param vcf_file The path to the genomic .vcf or .vcf.gz file(s)  to be imported. If
#' you specify a folder, the function will attempt to read all files in the
#' folder and combine them into on dataset. Multisample vcf files are not
#' supported; vcf files must contain one sample each. Required fields are
#' listed below
#'  - FIXED FIELDS:
#'  - `CHROM`: The reference sequence name.Equivalent to `contig`
#'  - `POS`: 0-based start position of the feature in contig.
#'  - `REF`: The reference allele at this position
#'  - `ALT`: The left-aligned, normalized, alternate allele at this position.
#' - FORMAT FIELDS:
#'  - `DP`: The total read depth at this position (including N-calls).
#'  - `VD`: Variant Depth. Equivalent to `alt_depth`.
#' - INFO FIELDS
#'  - `END`: The half-open end position of the feature in contig.
#'  - `sample`: An identifying field for your samples; either in the INFO
#' field or as the header to the FORMAT field.
#' - SUGGESTED INFO FIELDS:
#'  - `SVTYPE`: Structural variant types; INV DUP DEL INS FUS.
#'  - `SVLEN`: Length of the structural variant in base pairs.
#' @param sample_data_file An optional file containing additional sample
#' metadata (dose, timepoint, etc.). This can be a data frame or a file path.
#' @param sd_sep The delimiter for importing sample metadata tables.
#' Default is tab-delimited.
#' @param species The species of the reference genome. The value can be the
#' common name of the species or the scientific name. Ex. "human" or
#' "Homo sapiens".
#' @param genome The genome assembly of the reference genome. For a
#' complete list, refer to BSgenome::available.genomes(splitNameParts = TRUE)
#' Ex.Human GRCh38 = hg38 | Human GRCh37 = hg19 | Mouse GRCm38 = mm10 |
#' Mouse GRCm39 = mm39 | Rat RGSC 6.0 = rn6 | Rat mRatBN7.2 = rn7.
#' @param masked_BS_genome A logical value indicating whether to use
#' the masked version of the BS genome. Default is FALSE.
#' should
#' @import vcfppR
#' @importFrom dplyr left_join mutate case_when
#' @importFrom tidyr separate_rows separate pivot_wider
#' @return a data frame of your imported mutation data.
#' @export
import_genomic_vcf <- function(vcf_file,
                               sample_data_file,
                               sd_sep,
                               species,
                               genome,
                               masked_BS_genome = FALSE,
                               summarize = FALSE) {
  ## Use tryCatch to catch any errors that might pop up
  message("Installing the reference genome.
          You may be prompted to update dependencies.\n")
  ref_genome <- install_ref_genome(organism = species,
                                   genome = genome,
                                   masked = masked_BS_genome)
  # Initialize empty lists to store variant information
  chr_list <- c()
  start_list <- c()
  ref_list <- c()
  alt_list <- c()
  id_list <- c()
  qual_list <- c()
  filter_list <- c()
  variation_type_list <- c()
  sample_list <- list()
  info_list <- list()
  dp_list <- list()
  vd_list <- list()

  # read the vcf file
  vcf <- vcfppR::vcfreader$new(vcf_file)

  message("\nStarting vcf import. This may take a while...\n")
  total_variants <- sum(vcf$variant())
  progress_bar <- txtProgressBar(min = 0, max = total_variants, style = 3)
  # Loop through each variant
  while (vcf$variant()) {
    # Extract information for the current variant
    chr_list <- c(chr_list, vcf$chr())
    start_list <- c(start_list, vcf$pos())
    ref_list <- c(ref_list, vcf$ref())
    alt_list <- c(alt_list, vcf$alt())
    id_list <- c(id_list, vcf$id())
    qual_list <- c(qual_list, vcf$qual())
    filter_list <- c(filter_list, vcf$filter())
    variation_type_list <- c(variation_type_list,
                             classify_variation(vcf$ref(), vcf$alt()))
    sample_list <- c(sample_list, vcf$samples())
    info_list <- c(info_list, vcf$info())
    dp_list <- c(dp_list, vcf$formatInt("DP"))
    vd_list <- c(vd_list, vcf$formatInt("VD"))

    setTxtProgressBar(progress_bar, length(chr_list))
  }

  close(progress_bar)

  dat <- data.frame(
    contig = chr_list,
    start = start_list,
    ref = ref_list,
    alt = alt_list,
    id = id_list,
    qual = qual_list,
    filter = filter_list,
    variation_type = variation_type_list,
    sample = I(sample_list),
    info = I(info_list), # need to split up
    alt_depth = I(vd_list),
    total_depth = I(dp_list)
  )
  message("vcf imported successfully!\n")
  # Split the INFO column
  dat <- dat %>%
    tidyr::separate_rows(info, sep = ";") %>%
    tidyr::separate(info, into = c("key", "value"), sep = "=") %>%
    tidyr::pivot_wider(names_from = key, values_from = value)

  dat$sample <- as.character(dat$sample)
  dat$alt_depth <- as.integer(dat$alt_depth)
  dat$total_depth <- as.integer(dat$total_depth)

  # Ensure contig col contains "chr" prefix
  # Assuming dat is your data frame or tibble
  dat$contig <- ifelse(startsWith(dat$contig, "chr"),
                       dat$contig,
                       paste0("chr", dat$contig))

  # Remove <NON_REF> from alt column (typical for GATK gvcfs)
  dat$alt <- gsub(",\\s*<NON_REF>|<NON_REF>\\s*,|\\s*,\\s*<NON_REF>", "",
                  dat$alt)

  # Bind sample data
  # Do not throw errors because loading vcfs takes priority.
  if (!is.null(sample_data_file)) {
    message("Importing the sample data...\n")
    if (is.data.frame(sample_data_file)) {
      sampledata <- sample_data_file
      if (nrow(sampledata == 0)) {
        warning("Could not join Sample data to the Mutation data because
        the sample data frame you've provided is empty.")
      } else {
        common_columns <- intersect(names(dat), names(sampledata))
        if (length(common_columns) == 0) {
          warning("Could not join Sample data to the Mutation data because
          the sample data frame does not share any column names with Mutation
          data")
        } else {
          dat <- dplyr::left_join(dat, sampledata)
        }
      }
    } else if (is.character(sample_data_file)) {
      sample_file <- file.path(sample_data_file)
      if (!file.exists(sample_file)) {
        warning("Could not join Sample data to the Mutation data because
        the sample data file path that you provided is invalid")
      } else {
        if (file.info(sample_file)$size == 0) {
          warning("Could not join Sample data to the Mutation data because
          the sample data file is empty.")
        } else {
          sampledata <- read.table(sample_file,
                                   sep = sd_sep,
                                   header = TRUE)
          if (ncol(sampledata) <= 1) {
            warning("Could not join Sample data to the Mutation data because
            the imported sample data only has one column. Check that sd_sep is
            set to properly reflect the delimiter used for the data you are
            importing")
          } else {
            common_columns <- intersect(names(dat), names(sampledata))
            if (length(common_columns) == 0) {
              warning("Could not join Sample data to the Mutation data because
              the sample data frame does not share any column names with Mutation
              data")
            } else {
              dat <- dplyr::left_join(dat, sampledata)
            }
          }
        }
      }
    } else {
      warning("Could not join Sample data to Mutation data because the
      sample_data_file must be either a data frame or a file path.")
    }
  }
  dat <- rename_columns(dat)
  dat$end <- as.integer(dat$end)
  #dat <- check_required_columns(dat, op$base_required_mut_cols)
  message("Creating additional subtype columns...")
  # retrieve context
  dat <- makeGRangesFromDataFrame(
    df = as.data.frame(dat),
    keep.extra.columns = TRUE,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
  )
  extract_context <- function(mutations,
                              bsgenome,
                              upstream = 1,
                              downstream = 1) {
    # Resize the dat to include the context
    expanded_ranges <- GenomicRanges::resize(x = dat,
                                             width = upstream + downstream + 1,
                                             fix = "center")
    # Extract the sequences from the BSgenome
    sequences <- Biostrings::getSeq(bsgenome, expanded_ranges)
    # Return the sequences
    return(sequences)
  }
  context <- extract_context(dat, ref_genome)
  dat$context <- context
  dat <- as.data.frame(dat)

  # Create subtype columns
  sub_dict <- c(
    "G>T" = "C>A", "G>A" = "C>T", "G>C" = "C>G",
    "A>G" = "T>C", "A>C" = "T>G", "A>T" = "T>A"
  )
  dat <- dat %>%
    dplyr::mutate(
                  nchar_ref = nchar(ref),
                  nchar_alt =
                    ifelse(!.data$variation_type %in% c("sv", "ambiguous", "no_variant"),
                           nchar(alt), NA),
                  varlen =
                    ifelse(.data$variation_type %in% c("insertion", "deletion", "complex"),
                           .data$nchar_alt - .data$nchar_ref,
                           ifelse(.data$variation_type %in% c("snv", "mnv"),
                                  .data$nchar_ref, NA)),
                  subtype = ifelse(.data$variation_type == "snv",
                                   paste0(.data$ref, ">", .data$alt), .data$variation_type),
                  normalized_subtype = ifelse(.data$subtype %in% names(sub_dict),
                                              sub_dict[subtype], .data$subtype),
                  short_ref = substr(dat$ref, 1, 1),
                  normalized_ref = dplyr::case_when(
                                                    substr(ref, 1, 1) == "A" ~ "T",
                                                    substr(ref, 1, 1) == "G" ~ "C",
                                                    substr(ref, 1, 1) == "C" ~ "C",
                                                    substr(ref, 1, 1) == "T" ~ "T"),
                  context_with_mutation =
                    ifelse(.data$variation_type == "snv",
                           paste0(stringr::str_sub(.data$context, 1, 1),
                                  "[", .data$subtype, "]",
                                  stringr::str_sub(.data$context, 3, 3)),
                           .data$variation_type),
                  normalized_context =
                    ifelse(stringr::str_sub(dat$context, 2, 2) %in% c("G", "A", "g", "a"),
                          mapply(function(x) MutSeqR::reverseComplement(x, case = "upper"), dat$context),
                          dat$context),
                  normalized_context_with_mutation =
                    ifelse(.data$variation_type == "snv",
                           paste0(stringr::str_sub(.data$normalized_context, 1, 1),
                                  "[", .data$normalized_subtype, "]",
                                  stringr::str_sub(.data$normalized_context, 3, 3)),
                           .data$variation_type),
                  gc_content = (stringr::str_count(string = dat$context, pattern = "G") +
                                stringr::str_count(string = dat$context, pattern = "C"))  / stringr::str_count(dat$context))
  if (!summarize) {
    message("Done!")
    return(dat)
  } else {
    # Calculate Mut freqs
    subtype_resolutions <- c("none", "base_6", "base_12", "base_96", "base_192")
    summary <- list()
    results <- list()  # Initialize an empty list to store results
    for (subtype in subtype_resolutions) {
      result <- calculate_mut_freq(dat,
                                     subtype_resolution = subtype,
                                     filter_germ = FALSE)  # Run the function
      results[[subtype]] <- result
    }
    return(summary)
  }
}