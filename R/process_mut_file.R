#' Import a .mut file
#'
#' Imports a .mut file into the local R environment.
#' @param mut_file The .mut file containing mutation data to be imported. If you
#' specify a folder, function will attempt to read all files in the folder and
#' combine them into a single data frame.
#' Columns required are: depth col = (depth & no_calls or total_depth), 
#' alt_depth, context, ref, variation_type, contig, start 
#' (Synonymous names are accepted)
#' @param rsids TRUE or FALSE; whether or not the .mut file contains rsID information 
#' (existing SNPs)
#' @param sample_data_file An optional file containing additional sample metadata 
#' (dose, timepoint, etc.)
#' @param sd_sep The delimiter for importing sample metadata tables
#' @param mut_sep The delimiter for importing the .mut file
#' @param regions "human", "mouse", or "custom". The argument refers to the 
#' TS Mutagenesis panel of the specified species, or to a custom panel. 
#' If custom, provide file path in custom_regions_file. TO DO: add rat.
#' @param custom_regions_file "filepath". If regions is set to custom, 
#' provide the file path for the tab-delimited file containing regions metadata. 
#' Required columns are "contig", "start", and "end".
#' @param rg_sep The delimiter for importing the custom_regions_file. 
#' Default is tab-delimited.
#' @param vaf_cutoff Add a column to identify ostensibly germline variants using 
#' a cutoff for variant allele fraction (VAF). There is no default value provided, 
#' but generally a value of 0.1 (i.e., 10%) is a good starting point. Setting this 
#' will remove variants that are present at a frequency greater than this value 
#' at a given site.
#' @param depth_calc In the instance when there are two or more calls at the 
#' same location within a sample, and the depths differ, this parameter chooses 
#' the method of calculation for the total_depth. take_mean calculates the 
#' total_depth by taking the mean reference depth and then adding all the alt depths. 
#' take_del calculates the total_depth by choosing only the reference depth of 
#' the deletion in the group, or if no deletion is present, the complex variant,
#'  then adding all alt depths, if there is no deletion or complex variant, 
#'  then it takes the mean of the reference depths. Default is "take_del".
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
# sample_dat <- "C:/Users/ADODGE/OneDrive - HC-SC PHAC-ASPC/Documents/DupSeq R Package Building/Test Data/PRC_BM_sample_data.txt"
# mut_file <- "C:/Users/ADODGE/OneDrive - HC-SC PHAC-ASPC/Documents/DupSeq R Package Building/Test Data/prj00125_PRC_BM_variany-calls.genome.mut"
import_mut_data <- function(mut_file = "C:/Users/ADODGE/OneDrive - HC-SC PHAC-ASPC/Documents/DupSeq R Package Building/Test Data/mut files",
                            rsids = F,
                            sample_data_file = NULL,
                            sd_sep = "\t",
                            mut_sep = "\t",
                            vaf_cutoff,
                            regions = c("human", "mouse", "custom"),
                            custom_regions_file = NULL,
                            rg_sep = "\t",
                            depth_calc = "take_del") {

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

# Rename columns to default and check that all required columns are present.  
  dat <- rename_columns(dat)
  dat <- check_required_columns(dat, op$base_required_mut_cols)

  # Clean up data:
  # Get reverse complement of sequence context where mutation is listed on purine context
  # Change all purine substitutions to pyrimidine substitutions
  # Make new column with COSMIC-style 96 base context
  # TO DO - describe better

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
      nchar_ref = nchar(ref),
      nchar_alt = ifelse(variation_type != "symbolic" | variation_type != "sv", nchar(alt), NA),
      variation_type = tolower(dat$variation_type),
      variation_type = 
        ifelse(.data$variation_type == "sv" , "symbolic",
          ifelse(.data$variation_type == "indel" & .data$nchar_ref > .data$nchar_alt & .data$nchar_alt == 1, "deletion",
            ifelse(.data$variation_type == "indel" & .data$nchar_ref < .data$nchar_alt & .data$nchar_ref == 1, "insertion",
                   ifelse(.data$variation_type == "indel" & .data$nchar_ref != .data$nchar_alt & .data$nchar_alt > 1 & .data$nchar_ref > 1, "complex",
                 .data$variation_type)))),
      VARLEN = 
        ifelse(.data$variation_type %in% c("insertion", "deletion", "complex"), .data$nchar_alt - .data$nchar_ref,
               ifelse(.data$variation_type %in% c("snv", "mnv"), .data$nchar_ref,
                      NA)),
       ref_depth = .data[[depth_col]] - .data$alt_depth,
      subtype = 
        ifelse(.data$variation_type == "snv",
               paste0(.data$ref, ">", .data$alt),
               "."),
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

  # Calculate total_depth for the duplicated rows
  if (depth_calc == "take_del") {
    # Group by sample, contig, and start
  df_grouped <- dat %>%
    dplyr::group_by(.data$sample, .data$contig, .data$start)

 # Define a function to prioritize variation types
  prioritize_depth <- function(variation_type, total_depth) {
    if ("deletion" %in% variation_type) {
      total_depth[variation_type == "deletion"][1]  # Choose the total_depth of the first deletion
    } else if ("complex" %in% variation_type) {
      total_depth[variation_type == "complex"][1]  # Choose the total_depth of the first complex
    } else {
      total_depth[1]  # Choose the total_depth of the first row if no deletion or complex
    }
  }

  # Apply the function to each group
  dat <- df_grouped %>%
    dplyr::mutate(new_total_depth = prioritize_depth(variation_type, total_depth)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$total_depth) %>%
    dplyr::rename(total_depth = .data$new_total_depth)

   } else if (depth_calc == "take_mean") {
     df_grouped <- dat %>%
       dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
       dplyr::filter(n() > 1) %>%
       dplyr::summarize(
         new_total_depth = round(mean(.data$total_depth, na.rm = TRUE))
       ) %>%
       dplyr::ungroup()

     dat <- dat %>%
     dplyr::left_join(df_grouped, by = c("sample", "contig", "start")) 
     
     dat <- dat %>%
      dplyr::mutate(final_total_depth = ifelse(is.na(new_total_depth), total_depth, new_total_depth)) %>%
       dplyr::select(-total_depth, -new_total_depth) %>%
       dplyr::rename(total_depth = final_total_depth)

   } else {
     stop("Invalid depth_calc input. Please choose 'take_mean' or 'take_del'.")
   }

  # Create VAF column and is_germline column
  dat <- dat %>% dplyr::mutate(VAF = .data$alt_depth / .data$total_depth) %>%
    dplyr::mutate(is_germline = ifelse(.data$VAF < vaf_cutoff, F, T))

  # Turn into GRanges
  mut_ranges <- makeGRangesFromDataFrame(
    df = as.data.frame(dat),
    keep.extra.columns = T,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = TRUE
  )

  regions_df <- load_regions_file(regions, custom_regions_file, rg_sep)

  region_ranges <- makeGRangesFromDataFrame(
    df = regions_df,
    keep.extra.columns = T,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = TRUE
  )

  ranges_joined <- plyranges::join_overlap_left(mut_ranges, region_ranges)
  return(ranges_joined)
}
