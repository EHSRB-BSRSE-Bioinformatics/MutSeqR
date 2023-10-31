#' Import a .mut file
#'
#' Imports a .mut file into the local R environment.
#' @param mut_file The .mut file containing mutation data to be imported. If you
#' specify a folder, function will attempt to read all files in the folder and
#' combine them into a single data frame.
#' Columns required are: depth col = (depth & no_calls or total_depth), alt_depth, context, ref, variation_type, contig, start (Synonymous names are accepted)
#' @param rsids TRUE or FALSE; whether or not the .mut file contains rsID information (existing SNPs)
#' @param sample_data_file An optional file containing additional sample metadata (dose, timepoint, etc.)
#' @param sd_sep The delimiter for importing sample metadata tables
#' @param mut_sep The delimiter for importing the .mut file
#' @param regions_file "human", "mouse", or "custom". The argument refers to the TS Mutagenesis panel of the specified species, or to a custom panel. If custom, provide file path in custom_regions_file. TO DO: add rat.
#' @param custom_regions_file "filepath". If regions_file is set to custom, provide the file path for the tab-delimited file containing regions metadata. Required columns are "contig", "start", and "end".
#' @param rg_sep The delimiter for importing the custom_regions_file. Default is tab-delimited.
#' @param vaf_cutoff Add a column to identify ostensibly germline variants using a cutoff for variant allele fraction (VAF). There is no default value provided, but generally a value of 0.1 (i.e., 10%) is a good starting point. Setting this will remove variants that are present at a frequency greater than this value at a given site.
#' @returns A table where each row is a mutation, and columns indicate the location, type, and other data.
#' @importFrom dplyr bind_rows mutate left_join case_when
#' @importFrom magrittr %>%
#' @importFrom stringr str_sub str_count
#' @importFrom plyranges join_overlap_left
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom utils read.delim read.table
#' @importFrom rlang .data
#' @export

# To delete later:
# C:/Users/ADODGE/OneDrive - HC-SC PHAC-ASPC/Documents/DupSeq R Package Building/Test Data/mut files
# inst/extdata/genic_regions_mm10.txt

import_mut_data <- function(mut_file = "../../data/Jonatan_Mutations_in_blood_and_sperm_samples_221021_MM.txt",
                            rsids = F,
                            sample_data_file = NULL,
                            sd_sep = "\t",
                            mut_sep = "\t",
                            vaf_cutoff,
                            regions_file = c("human", "mouse", "custom"),
                            custom_regions_file = NULL,
                            rg_sep = "\t") {
  # col name synonyms
  column_name_mapping <- c(
    "chromosome" = "contig",
    "chr" = "contig",
    "position" = "start",
    "pos" = "start",
    "sample_id" = "sample",
    "Sample" = "sample",
    "variant_type" = "variation_type",
    "mutation_type" = "variation_type",
    "reference" = "ref",
    "ref_allele" = "ref",
    "alternate" = "alt",
    "alt_allele" = "alt",
    "alt_read_depth" = "alt_depth",
    "coverage" = "depth",
    "read_depth" = "depth",
    "no_depth" = "no_calls",
    "n_calls" = "no_calls",
    "mutation_subtype" = "subtype",
    "sequence_context" = "context",
    "flanking_sequence" = "context"
  )


  mut_file <- file.path(mut_file)
  if (file.info(mut_file)$isdir == T) {
    mut_files <- list.files(path = mut_file, full.names = T)
    # Read in the files and bind them together
    dat <- lapply(mut_files, function(file) {
      read.table(file,
        header = TRUE, sep = mut_sep,
        fileEncoding = "UTF-8-BOM"
      )
    }) %>% dplyr::bind_rows()
  } else {
    dat <- read.table(mut_file,
      header = T, sep = mut_sep,
      fileEncoding = "UTF-8-BOM"
    )
  }
  if (ncol(dat) <= 1) {
    stop("Your imported data only has one column.
                           You may want to set mut_sep to properly reflect
                           the delimiter used for the data you are importing.")
  }
  if (rsids == T) {
    if (!"id" %in% colnames(dat)) {
      stop("Error: you have set rsids to TRUE,
      but there is no id column in the mut file!")
    }
    # If we have rs IDs, add a column indicating whether the mutation is a known SNP

    dat <- dat %>% dplyr::mutate(is_known = ifelse(!id == ".", "Y", "N"))

  }

  # Read in sample data if it's provided
  if (!is.null(sample_data_file)) {

    sampledata <- read.delim(file.path(sample_data_file), sep = sd_sep,
                             header = T)
    dat <- dplyr::left_join(dat, sampledata, suffix = c("", ".sampledata"))

  }

  # Change column names to default
  for (synonym in names(column_name_mapping)) {
    matching_col <- which(colnames(dat) %in% synonym)
    if (length(matching_col) > 0) {
      colnames(dat)[matching_col] <- column_name_mapping[synonym]
    }
  }

  ################
  # Clean up data:
  # Get reverse complement of sequence context where mutation is listed on purine context
  # Change all purine substitutions to pyrimidine substitutions
  # Make new column with COSMIC-style 96 base context
  # TODO - describe better

  # Define substitution dictionary to normalize to pyrimidine context
  sub_dict <- c(
    "G>T" = "C>A", "G>A" = "C>T", "G>C" = "C>G",
    "A>G" = "T>C", "A>C" = "T>G", "A>T" = "T>A"
  )

  # The column that represents depth might vary
  depth_col <- ifelse("total_depth" %in% colnames(dat),
    "total_depth",
    ifelse("depth" %in% colnames(dat),
      "depth", stop("Error: I'm not sure which column
                             specifies depth.")
    )
  )

  dat <- dat %>%

    dplyr::mutate(
      ref_depth = .data[[depth_col]] - .data$alt_depth,
      context_with_mutation =
        ifelse(.data$subtype != ".",
          paste0(
            stringr::str_sub(.data$context, 1, 1),
            "[", .data$subtype, "]",
            stringr::str_sub(.data$context, 3, 3)
          ),
          .data$variation_type
        ),
      normalized_context = ifelse(
        stringr::str_sub(.data$context, 2, 2) %in% c("G", "A", "g", "a"),
        mapply(function(x) reverseComplement(x, case = "upper"), .data$context),
        .data$context
      ),
      normalized_subtype = ifelse(
        .data$subtype %in% names(sub_dict),
        sub_dict[.data$subtype],
        .data$subtype
      ),
      short_ref = substr(.data$ref, 1, 1),
      normalized_ref = dplyr::case_when(
        substr(.data$ref, 1, 1) == "A" ~ "T",
        substr(.data$ref, 1, 1) == "G" ~ "C",
        substr(.data$ref, 1, 1) == "C" ~ "C",
        substr(.data$ref, 1, 1) == "T" ~ "T"
      )

    ) %>%
    dplyr::mutate(
      normalized_context_with_mutation =
        ifelse(.data$subtype != ".",
          paste0(
            stringr::str_sub(.data$normalized_context, 1, 1),
            "[", .data$normalized_subtype, "]",
            stringr::str_sub(.data$normalized_context, 3, 3)
          ),
          .data$variation_type
        ),
      gc_content = (stringr::str_count(string = .data$context, pattern = "G") +
        stringr::str_count(string = .data$context, pattern = "C"))
      / stringr::str_count(.data$context)
    ) %>%
    dplyr::mutate(
      normalized_subtype = ifelse(
        .data$normalized_subtype == ".",
        .data$variation_type,
        .data$normalized_subtype
      ),
      subtype = ifelse(
        .data$subtype == ".",
        .data$variation_type,
        .data$subtype
      )
    ) %>%

    {
      if ("depth" %in% names(.)) {
        dplyr::mutate(., total_depth = .data$depth - .data$no_calls)
      } else {
        .
      }
    }

  dat <- dat %>% dplyr::mutate(VAF = .data$alt_depth / .data$total_depth) %>%
    dplyr::mutate(is_germline = ifelse(.data$VAF < vaf_cutoff, F, T))

  mut_ranges <- makeGRangesFromDataFrame(
    df = as.data.frame(dat),
    keep.extra.columns = T,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = TRUE
  )

  genic_regions <- load_regions_file(regions_file, custom_regions_file, rg_sep)

  region_ranges <- makeGRangesFromDataFrame(
    df = genic_regions,
    keep.extra.columns = T,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = TRUE
  )

  ranges_joined <- plyranges::join_overlap_left(mut_ranges, region_ranges)
  return(ranges_joined)
}
