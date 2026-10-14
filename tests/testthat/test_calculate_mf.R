test_that("calculate_mf works in general", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  mf_data <- calculate_mf(
    mutation_data = mutation_data,
    cols_to_group = "sample",
    subtype_resolution = "none",
    calculate_depth = TRUE,
    correct_depth = TRUE,
    correct_depth_by_indel_priority = TRUE
  )
  expect_setequal(
    colnames(mf_data),
    c("sample", "sum_min", "sum_max", "group_depth", "mf_min", "mf_max")
  )
  expect_equal(nrow(mf_data), 3)
  expect_equal(mf_data$sum_min, c(7, 10, 7))
  expect_equal(mf_data$sum_max, c(230, 447, 140))
  expect_equal(mf_data$group_depth, c(1441, 1920, 1390)) # correct depth works
  expect_equal(mf_data$mf_min, mf_data$sum_min / mf_data$group_depth)
  expect_equal(mf_data$mf_max, mf_data$sum_max / mf_data$group_depth)
})
test_that("calculate_mf works with filtering", {
  filtered_mutation_data <- readRDS("./testdata/filtered_simple_mutation_data.rds")
  mf_data <- calculate_mf(
    mutation_data = filtered_mutation_data,
    cols_to_group = "sample",
    subtype_resolution = "none",
    calculate_depth = TRUE,
    correct_depth = TRUE,
    correct_depth_by_indel_priority = TRUE
  )
  expect_equal(mf_data$sum_min, c(2, 3, 3))
  expect_equal(mf_data$sum_max, c(3, 4, 3))
  expect_equal(mf_data$group_depth, c(1440, 1909, 1390))
  expect_equal(mf_data$mf_min, mf_data$sum_min / mf_data$group_depth)
  expect_equal(mf_data$mf_max, mf_data$sum_max / mf_data$group_depth)
})

test_that("calculate_mf works with subtypes", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  mf_data <- calculate_mf(
    mutation_data = mutation_data,
    cols_to_group = "sample",
    subtype_resolution = "base_6",
    calculate_depth = TRUE,
    correct_depth = TRUE,
    correct_depth_by_indel_priority = TRUE
  )
  expect_setequal(
    colnames(mf_data),
    c("sample", "normalized_subtype", "normalized_ref",
      "sum_min", "sum_max", "group_depth", "subtype_depth",
      "mf_min", "mf_max", "proportion_min", "proportion_max")
  )
  expect_setequal(
    unique(mf_data$normalized_subtype),
    c(setdiff(MutSeqR::subtype_list$type, c("no_variant", "snv")),
      MutSeqR::subtype_list$base_6)
  )
  expect_setequal(
    unique(mf_data$normalized_ref),
    c("N", "C", "T")
  )
  expect_subtype_depth <- function(sample, subtypes, expected_depth) {
    rows <- mf_data$sample == sample & mf_data$normalized_subtype %in% subtypes
    expect_true(all(mf_data$subtype_depth[rows] == expected_depth))
  }
  expect_subtype_depth("sample1", c(MutSeqR::subtype_list$type), 1441)
  expect_subtype_depth("sample1", c("C>A", "C>G", "C>T"), 691)
  expect_subtype_depth("sample1", c("T>A", "T>C", "T>C"), 750)
  expect_subtype_depth("sample2", c(MutSeqR::subtype_list$type), 1920)
  expect_subtype_depth("sample2", c("C>A", "C>G", "C>T"), 820)
  expect_subtype_depth("sample2", c("T>A", "T>C", "T>C"), 1100)
  expect_subtype_depth("sample3", c(MutSeqR::subtype_list$type), 1390)
  expect_subtype_depth("sample3", c("C>A", "C>G", "C>T"), 670)
  expect_subtype_depth("sample3", c("T>A", "T>C", "T>C"), 720)

  expect_equal(unique(mf_data$group_depth[mf_data$sample == "sample1"]), 1441)
  expect_equal(unique(mf_data$group_depth[mf_data$sample == "sample2"]), 1920)
  expect_equal(unique(mf_data$group_depth[mf_data$sample == "sample3"]), 1390)

  expect_equal(
    mf_data$sum_min,
    c(0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 3, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0,
      0, 1, 0, 1, 1, 2, 0, 1, 0, 2, 1, 0, 0, 0, 0)
  )
  expect_equal(
    mf_data$sum_max,
    c(0, 0, 1, 0, 3, 0, 3, 0, 4, 2, 2, 2, 0, 339, 1, 1, 0, 2, 0, 0, 0, 1, 1, 0,
      0, 10, 0, 220, 90, 130, 0, 1, 0, 3, 1, 0, 0, 0, 0)
  )
  expect_equal(mf_data$mf_min, mf_data$sum_min / mf_data$subtype_depth)
  expect_equal(mf_data$mf_max, mf_data$sum_max / mf_data$subtype_depth)
  expected_prop_min <- vapply(
    seq_len(nrow(mf_data)),
    function(i) {
      sample <- mf_data$sample[i]
      mf_data$mf_min[i] / sum(mf_data$mf_min[mf_data$sample == sample])
    },
    numeric(1)
  )
  expect_equal(mf_data$proportion_min, expected_prop_min)
  expected_prop_max <- vapply(
    seq_len(nrow(mf_data)),
    function(i) {
      sample <- mf_data$sample[i]
      mf_data$mf_max[i] / sum(mf_data$mf_max[mf_data$sample == sample])
    },
    numeric(1)
  )
  expect_equal(mf_data$proportion_max, expected_prop_max)
})

test_that("calculate_mf works with precalculated depth", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  depth <- data.frame(
    sample = c("sample1", "sample2", "sample3"),
    group_depth = c(1500, 1600, 1400)
  )
  mf_data <- calculate_mf(
    mutation_data = mutation_data,
    cols_to_group = "sample",
    subtype_resolution = "none",
    calculate_depth = FALSE,
    precalc_depth_data = depth
  )
  expect_equal(mf_data$sum_min, c(7, 10, 7))
  expect_equal(mf_data$sum_max, c(230, 447, 140))
  expect_equal(mf_data$group_depth, c(1500, 1600, 1400))
  expect_equal(mf_data$mf_min, mf_data$sum_min / mf_data$group_depth)
  expect_equal(mf_data$mf_max, mf_data$sum_max / mf_data$group_depth)

  depth <- data.frame(
    sample = c(rep("sample1", 2), rep("sample2", 2), rep("sample3", 2)),
    normalized_ref = rep(c("C", "T"), 3),
    group_depth = c(rep(1300, 2), rep(700, 2), rep(900, 2)),
    subtype_depth = c(600, 700, 400, 300, 500, 400)
  )
  mf_data <- calculate_mf(
    mutation_data = mutation_data,
    cols_to_group = "sample",
    subtype_resolution = "base_6",
    calculate_depth = FALSE,
    precalc_depth_data = depth
  )
  expect_subtype_depth <- function(sample, subtypes, expected_depth) {
    rows <- mf_data$sample == sample & mf_data$normalized_subtype %in% subtypes
    expect_true(all(mf_data$subtype_depth[rows] == expected_depth))
  }
  expect_subtype_depth("sample1", c(MutSeqR::subtype_list$type), 1300)
  expect_subtype_depth("sample1", c("C>A", "C>G", "C>T"), 600)
  expect_subtype_depth("sample1", c("T>A", "T>C", "T>C"), 700)
  expect_subtype_depth("sample2", c(MutSeqR::subtype_list$type), 700)
  expect_subtype_depth("sample2", c("C>A", "C>G", "C>T"), 400)
  expect_subtype_depth("sample2", c("T>A", "T>C", "T>C"), 300)
  expect_subtype_depth("sample3", c(MutSeqR::subtype_list$type), 900)
  expect_subtype_depth("sample3", c("C>A", "C>G", "C>T"), 500)
  expect_subtype_depth("sample3", c("T>A", "T>C", "T>C"), 400)

  expect_equal(mf_data$mf_min, mf_data$sum_min / mf_data$subtype_depth)
  expect_equal(mf_data$mf_max, mf_data$sum_max / mf_data$subtype_depth)
  expected_prop_min <- vapply(
    seq_len(nrow(mf_data)),
    function(i) {
      sample <- mf_data$sample[i]
      mf_data$mf_min[i] / sum(mf_data$mf_min[mf_data$sample == sample])
    },
    numeric(1)
  )
  expect_equal(mf_data$proportion_min, expected_prop_min)
  expected_prop_max <- vapply(
    seq_len(nrow(mf_data)),
    function(i) {
      sample <- mf_data$sample[i]
      mf_data$mf_max[i] / sum(mf_data$mf_max[mf_data$sample == sample])
    },
    numeric(1)
    )
  expect_equal(mf_data$proportion_max, expected_prop_max)
})

test_that("calculate_mf derives requested resolution from finer precalc data", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  samples <- unique(mutation_data$sample)
  base96 <- expand.grid(
    normalized_context = MutSeqR::context_list$base_96,
    sample = samples,
    stringsAsFactors = FALSE
  )
  base96$subtype_depth <- seq_len(nrow(base96))
  derived <- suppressWarnings(
    MutSeqR::write_depth_data(base96, d_sep = "\t")
  )

  global_mf <- calculate_mf(
    mutation_data, subtype_resolution = "none", calculate_depth = FALSE,
    precalc_depth_data = base96
  )
  base6_mf <- calculate_mf(
    mutation_data, subtype_resolution = "base_6", calculate_depth = FALSE,
    precalc_depth_data = base96
  )
  base96_mf <- calculate_mf(
    mutation_data, subtype_resolution = "base_96", calculate_depth = FALSE,
    variant_types = "snv", precalc_depth_data = base96
  )
  base96_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(base96_file), add = TRUE)
  utils::write.table(base96, base96_file, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  global_file_mf <- calculate_mf(
    mutation_data, subtype_resolution = "none", calculate_depth = FALSE,
    precalc_depth_data = base96_file
  )
  base192 <- expand.grid(
    context = MutSeqR::context_list$base_192,
    sample = samples,
    stringsAsFactors = FALSE
  )
  base192$subtype_depth <- seq_len(nrow(base192))
  base12_from_base192 <- calculate_mf(
    mutation_data, subtype_resolution = "base_12", calculate_depth = FALSE,
    variant_types = "snv", precalc_depth_data = base192
  )
  mutation_data$dose_group <- ifelse(
    mutation_data$sample == samples[[1]], "control", "treated"
  )
  grouped_depth <- expand.grid(
    dose_group = c("control", "treated"),
    normalized_context = MutSeqR::context_list$base_96,
    stringsAsFactors = FALSE
  )
  grouped_depth$subtype_depth <- seq_len(nrow(grouped_depth))
  grouped_global_mf <- calculate_mf(
    mutation_data, cols_to_group = "dose_group", subtype_resolution = "none",
    calculate_depth = FALSE, precalc_depth_data = grouped_depth
  )
  grouped_base6_mf <- calculate_mf(
    mutation_data, cols_to_group = "dose_group", subtype_resolution = "base_6",
    calculate_depth = FALSE, precalc_depth_data = grouped_depth
  )
  grouped_base96_mf <- calculate_mf(
    mutation_data, cols_to_group = "dose_group", subtype_resolution = "base_96",
    variant_types = "snv", calculate_depth = FALSE,
    precalc_depth_data = grouped_depth
  )
  grouped_mf <- calculate_mf(
    mutation_data, cols_to_group = "dose_group", subtype_resolution = "none",
    calculate_depth = FALSE, precalc_depth_data = base96
  )
  expected_group_depth <- derived$global %>%
    dplyr::left_join(
      dplyr::distinct(mutation_data, sample, dose_group),
      by = "sample"
    ) %>%
    dplyr::group_by(dose_group) %>%
    dplyr::summarise(group_depth = sum(group_depth), .groups = "drop") %>%
    dplyr::arrange(dose_group)

  expect_equal(global_mf$group_depth, derived$global$group_depth)
  expect_equal(global_file_mf$group_depth, derived$global$group_depth)
  expect_true(all(c("mf_min", "mf_max") %in% names(base6_mf)))
  expect_true(all(c("mf_min", "mf_max") %in% names(base96_mf)))
  expect_true(all(c("mf_min", "mf_max") %in%
                    names(base12_from_base192)))
  expect_true(all(c("mf_min", "mf_max") %in%
                    names(grouped_global_mf)))
  expect_true(all(c("mf_min", "mf_max") %in%
                    names(grouped_base6_mf)))
  expect_true(all(c("mf_min", "mf_max") %in%
                    names(grouped_base96_mf)))
  expect_equal(
    as.data.frame(
      dplyr::distinct(grouped_mf, dose_group, group_depth) %>%
        dplyr::arrange(dose_group)
    ),
    as.data.frame(expected_group_depth)
  )
})

test_that("calculate_mf accepts precalculated depth at all grouping columns", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  mutation_data$contig <- rep(c("contig_1", "contig_2"),
                               length.out = nrow(mutation_data))
  groups <- dplyr::distinct(mutation_data, sample, contig)
  contig_depth <- merge(
    groups,
    data.frame(
      normalized_context = MutSeqR::context_list$base_96,
      stringsAsFactors = FALSE
    ),
    by = NULL
  )
  contig_depth$subtype_depth <- seq_len(nrow(contig_depth))

  mf_data <- calculate_mf(
    mutation_data,
    cols_to_group = c("sample", "contig"),
    subtype_resolution = "base_6",
    calculate_depth = FALSE,
    precalc_depth_data = contig_depth
  )

  expect_true(all(c("sample", "contig", "mf_min", "mf_max") %in%
                    names(mf_data)))
  expect_false(anyNA(mf_data$subtype_depth))
})

test_that("calculate_mf preserves absent contexts as zero depth and NA MF", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  present_contexts <- unique(mutation_data$normalized_context[
    mutation_data$normalized_context %in% MutSeqR::context_list$base_96
  ])
  missing_context <- present_contexts[[1]]
  mutation_data <- dplyr::filter(
    mutation_data,
    .data$normalized_context != missing_context
  )

  mf_data <- calculate_mf(
    mutation_data,
    subtype_resolution = "base_96",
    variant_types = "snv"
  )
  absent_context_rows <- mf_data$normalized_context == missing_context

  expect_true(all(mf_data$subtype_depth[absent_context_rows] == 0))
  expect_true(all(mf_data$sum_min[absent_context_rows] == 0))
  expect_true(all(mf_data$sum_max[absent_context_rows] == 0))
  expect_true(all(is.na(mf_data$mf_min[absent_context_rows])))
  expect_true(all(is.na(mf_data$mf_max[absent_context_rows])))
  expect_true(all(is.na(mf_data$proportion_min[absent_context_rows])))
  expect_true(all(is.na(mf_data$proportion_max[absent_context_rows])))
})

test_that("calculate_mf rejects precalc data too coarse for requested resolution", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  base96 <- expand.grid(
    normalized_context = MutSeqR::context_list$base_96,
    sample = unique(mutation_data$sample),
    stringsAsFactors = FALSE
  )
  base96$subtype_depth <- seq_len(nrow(base96))

  expect_error(
    calculate_mf(
      mutation_data, subtype_resolution = "base_12", calculate_depth = FALSE,
      precalc_depth_data = base96
    ),
    "at 'base_96' resolution cannot generate depths.*'base_12'"
  )
})
test_that("calculate_mf selects variation types", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  mf_data <- calculate_mf(
    mutation_data = mutation_data,
    cols_to_group = "sample",
    subtype_resolution = "none",
    variant_types = c("insertion", "deletion", "complex"),
    calculate_depth = TRUE,
    correct_depth_by_indel = TRUE
  )
  expect_equal(mf_data$sum_min, c(2, 2, 2))
  expect_equal(mf_data$sum_max, c(5, 5, 6))
  expect_equal(mf_data$group_depth, c(1441, 1920, 1390))

  mf_data <- calculate_mf(
    mutation_data = mutation_data,
    cols_to_group = "sample",
    subtype_resolution = "type",
    variant_types = c("-ambiguous", "-uncategorized", "-complex"),
    calculate_depth = TRUE,
    correct_depth_by_indel = TRUE
  )
  expect_setequal(
    unique(mf_data$variation_type),
    c("deletion", "insertion", "mnv", "snv", "sv")
  )
  expect_equal(unique(mf_data$group_depth[mf_data$sample == "sample1"]), 1441)
  expect_equal(unique(mf_data$group_depth[mf_data$sample == "sample2"]), 1920)
  expect_equal(unique(mf_data$group_depth[mf_data$sample == "sample3"]), 1390)
  expect_equal(mf_data$sum_min, c(1, 0, 1, 1, 1, 1, 0, 3, 1, 1, 0, 1, 4, 5, 2))
  expect_equal(
    mf_data$sum_max,
    c(3, 0, 4, 2, 2, 2, 0, 339, 1, 1, 0, 2, 224, 103, 130)
  )
})

test_that("calculate_mf returns mutation-level frequencies without summarizing", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  summary_data <- calculate_mf(
    mutation_data,
    subtype_resolution = "none",
    correct_depth = FALSE
  )
  mf_data <- calculate_mf(
    mutation_data,
    subtype_resolution = "none",
    correct_depth = FALSE,
    summary = FALSE
  )

  expect_equal(nrow(mf_data), nrow(mutation_data))
  expect_true(all(names(mutation_data) %in% names(mf_data)))
  expect_equal(mf_data[names(mutation_data)], tibble::as_tibble(mutation_data))
  group_rows <- match(mf_data$sample, summary_data$sample)
  for (column in c("sum_min", "sum_max", "group_depth", "mf_min", "mf_max")) {
    expect_equal(mf_data[[column]], summary_data[[column]][group_rows])
  }
})

test_that("mutation-level results skip missing summary metadata warnings", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  present_contexts <- unique(mutation_data$normalized_context[
    mutation_data$normalized_context %in% MutSeqR::context_list$base_96
  ])
  mutation_data <- dplyr::filter(
    mutation_data,
    .data$normalized_context != present_contexts[[1]]
  )
  mutation_data$metadata_label <- "retained"

  expect_warning(
    calculate_mf(
      mutation_data,
      subtype_resolution = "base_96",
      variant_types = "snv",
      correct_depth = FALSE,
      retain_metadata_cols = "metadata_label"
    ),
    "NA values detected in retained metadata columns"
  )
  expect_warning(
    mf_data <- calculate_mf(
      mutation_data,
      subtype_resolution = "base_96",
      variant_types = "snv",
      correct_depth = FALSE,
      summary = FALSE,
      retain_metadata_cols = "metadata_label"
    ),
    NA
  )
  expect_equal(nrow(mf_data), nrow(mutation_data))
  expect_equal(mf_data$metadata_label, mutation_data$metadata_label)
  expect_true(all(c("mf_min", "mf_max") %in% names(mf_data)))
  expect_false(any(c("proportion_min", "proportion_max") %in% names(mf_data)))
})


test_that("summary context mapping preserves all subtype resolutions", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  for (resolution in c("base_6", "base_12", "base_96", "base_192")) {
    mf_data <- calculate_mf(
      mutation_data,
      subtype_resolution = resolution,
      correct_depth = FALSE
    )
    subtype_column <- MutSeqR::subtype_dict[[resolution]]
    context_column <- MutSeqR::denominator_dict[[resolution]]
    # Compare with the former scalar, per-row context calculation.
    expected_contexts <- vapply(
      as.character(mf_data[[subtype_column]]),
      MutSeqR:::get_ref_of_mut,
      character(1),
      USE.NAMES = FALSE
    )
    expect_equal(
      mf_data[[context_column]],
      tidyr::replace_na(expected_contexts, "N")
    )
    expect_false(anyNA(mf_data[[context_column]]))
    expect_false(".missing_context_depth" %in% names(mf_data))
  }
})


test_that("calculate_mf rejects reserved missing-context column names", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  for (column in c("missing_context_depth", ".missing_context_depth")) {
    input <- mutation_data
    input[[column]] <- FALSE
    for (return_summary in c(TRUE, FALSE)) {
      expect_error(
        calculate_mf(input, summary = return_summary),
        "reserved column name.*Rename these columns before calling calculate_mf"
      )
    }
  }
})

# Use unequal sample depths to distinguish pooled MF from mean sample MF.
pooled_depth_fixture <- function() {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  samples <- unique(mutation_data$sample)
  mutation_data$dose <- ifelse(mutation_data$sample == samples[1], 0, 10)
  mutation_data$visit <- ifelse(mutation_data$dose == 0, "early", "late")
  depth <- data.frame(sample = samples, group_depth = c(1000, 2000, 4000))
  list(mutations = mutation_data, depth = depth)
}

test_that("global sample depths pool to arbitrary metadata combinations", {
  fixture <- pooled_depth_fixture()
  mutations <- fixture$mutations
  depth <- fixture$depth
  depth_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(depth_file), add = TRUE)
  utils::write.table(depth, depth_file, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  mapping <- dplyr::distinct(mutations, sample, dose, visit)
  supplied <- dplyr::left_join(depth, mapping, by = "sample")
  mutations$exposure <- mutations$dose

  for (groups in list("dose", "exposure", c("dose", "visit"))) {
    grouped_depth <- dplyr::left_join(
      depth, dplyr::distinct(mutations[, c("sample", groups)]), by = "sample"
    ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(groups))) %>%
      dplyr::summarise(group_depth = sum(group_depth), .groups = "drop")
    for (resolution in c("none", "type")) {
      for (summarize in c(TRUE, FALSE)) {
        expected <- calculate_mf(
          mutations, cols_to_group = groups, subtype_resolution = resolution,
          calculate_depth = FALSE, precalc_depth_data = grouped_depth,
          summary = summarize
        )
        # The direct result provides the existing schema and row-order contract.
        for (input in list(depth, depth_file)) {
          actual <- calculate_mf(
            mutations, cols_to_group = groups, subtype_resolution = resolution,
            calculate_depth = FALSE, precalc_depth_data = input, depth_group_cols = "sample",
            summary = summarize
          )
          expect_equal(actual, expected)
          if (!summarize) expect_equal(nrow(actual), nrow(mutations))
        }
      }
    }
  }
  pooled <- calculate_mf(mutations, cols_to_group = "dose",
                         calculate_depth = FALSE, precalc_depth_data = depth, depth_group_cols = "sample")
  expect_equal(pooled$group_depth[pooled$dose == 10], 6000)
  expect_equal(pooled$mf_max, pooled$sum_max / pooled$group_depth)
  expect_equal(pooled$mf_min, pooled$sum_min / pooled$group_depth)
  with_metadata <- calculate_mf(
    mutations, cols_to_group = c("dose", "visit"), calculate_depth = FALSE,
    precalc_depth_data = supplied, depth_group_cols = "sample"
  )
  without_metadata <- calculate_mf(
    mutations, cols_to_group = c("dose", "visit"), calculate_depth = FALSE,
    precalc_depth_data = depth, depth_group_cols = "sample"
  )
  expect_equal(with_metadata, without_metadata)
})

test_that("global pooling rejects ambiguous and inconsistent inputs", {
  fixture <- pooled_depth_fixture()
  mutations <- fixture$mutations
  depth <- fixture$depth
  run <- function(input, data = mutations, groups = "dose") {
    calculate_mf(data, cols_to_group = groups, calculate_depth = FALSE,
                 precalc_depth_data = input, depth_group_cols = "sample")
  }
  expect_error(run(rbind(depth, depth[1, ])), "one depth row per sample.*sample1")
  unknown <- depth
  unknown$sample[1] <- "absent_sample"
  expect_error(run(unknown), "absent from mutation_data: absent_sample")
  for (value in c(-1, NA_real_, Inf, NaN)) {
    invalid <- depth
    invalid$group_depth[1] <- value
    expect_error(run(invalid), "finite, non-negative")
  }
  invalid <- depth
  invalid$group_depth <- as.character(invalid$group_depth)
  expect_error(run(invalid), "finite, non-negative")
  conflicting <- dplyr::left_join(
    depth, dplyr::distinct(mutations, sample, dose), by = "sample"
  )
  conflicting$dose[1] <- 10
  expect_error(run(conflicting), "disagree with mutation_data.*sample1")
  for (column in c("locus", "timepoint")) {
    split <- mutations
    split[[column]] <- rep(c("a", "b"), length.out = nrow(split))
    expect_error(run(depth, split, c("dose", column)),
                 paste0("cols_to_group.*", column, ".*Ambiguous sample"))
  }
  missing_sample <- mutations
  missing_sample$sample <- NULL
  expect_error(run(depth, missing_sample), "mutation_data is missing: sample")
})

test_that("direct global inputs and missing depths retain their behavior", {
  fixture <- pooled_depth_fixture()
  mutations <- fixture$mutations
  direct <- data.frame(dose = c(0, 10), group_depth = c(1000, 6000))
  baseline <- calculate_mf(mutations, cols_to_group = "dose",
                           calculate_depth = FALSE, precalc_depth_data = direct)
  # An incidental sample column does not override an already grouped table.
  with_sample <- direct
  with_sample$sample <- c("incidental_a", "incidental_b")
  actual <- calculate_mf(mutations, cols_to_group = "dose",
                         calculate_depth = FALSE, precalc_depth_data = with_sample)
  expect_equal(actual[names(baseline)], baseline)
  expect_equal(actual$sample, with_sample$sample)
  expect_error(calculate_mf(
    mutations, cols_to_group = "dose", calculate_depth = FALSE,
    precalc_depth_data = rbind(direct, direct[1, ])
  ), "finer grouping than required")
  # Missing-depth and zero-depth handling is identical to the direct path.
  for (depth in list(fixture$depth[-1, ],
                    transform(fixture$depth, group_depth = 0))) {
    prepared <- dplyr::left_join(
      depth, dplyr::distinct(mutations, sample, dose), by = "sample"
    ) %>%
      dplyr::group_by(dose) %>%
      dplyr::summarise(group_depth = sum(group_depth), .groups = "drop")
    if (nrow(depth) < nrow(fixture$depth)) {
      expect_warning(
        actual <- calculate_mf(
          mutations, cols_to_group = "dose", calculate_depth = FALSE,
          precalc_depth_data = depth, depth_group_cols = "sample"
        ), "No precalculated depth.*sample1.*Affected groups"
      )
    } else {
      actual <- calculate_mf(
        mutations, cols_to_group = "dose", calculate_depth = FALSE,
        precalc_depth_data = depth, depth_group_cols = "sample"
      )
    }
    expect_equal(actual, calculate_mf(
      mutations, cols_to_group = "dose", calculate_depth = FALSE,
      precalc_depth_data = prepared
    ))
  }
  # Existing depth-column suffix precedence is unchanged during pooling.
  mutations$group_depth <- -1
  result <- calculate_mf(mutations, cols_to_group = "dose", summary = FALSE,
                         calculate_depth = FALSE, precalc_depth_data = fixture$depth, depth_group_cols = "sample")
  expect_true(all(result$group_depth.x == -1))
  expect_equal(result$group_depth, ifelse(mutations$dose == 0, 1000, 6000))
})

test_that("context pooling preserves separate global and context totals", {
  fixture <- pooled_depth_fixture()
  mutations <- fixture$mutations
  contexts <- expand.grid(sample = fixture$depth$sample,
                          normalized_ref = c("C", "T"),
                          stringsAsFactors = FALSE)
  contexts$subtype_depth <- ifelse(contexts$normalized_ref == "C", 100, 300)
  direct <- data.frame(dose = rep(c(0, 10), each = 2),
                       normalized_ref = rep(c("C", "T"), 2),
                       subtype_depth = c(100, 300, 200, 600))
  for (resolution in c("none", "type", "base_6")) {
    for (summarize in c(TRUE, FALSE)) {
      expect_equal(
        calculate_mf(mutations, cols_to_group = "dose",
                     subtype_resolution = resolution, summary = summarize,
                     calculate_depth = FALSE, precalc_depth_data = contexts),
        calculate_mf(mutations, cols_to_group = "dose",
                     subtype_resolution = resolution, summary = summarize,
                     calculate_depth = FALSE, precalc_depth_data = direct)
      )
    }
  }
  # Finer within-sample keys continue to work when supplied at requested groups.
  mutations$locus <- rep(c("a", "b"), length.out = nrow(mutations))
  locus_depth <- dplyr::distinct(mutations, sample, locus)
  locus_depth$group_depth <- 500
  result <- calculate_mf(mutations, cols_to_group = c("sample", "locus"),
                         calculate_depth = FALSE, precalc_depth_data = locus_depth)
  expect_true(all(result$group_depth == 500))
  expect_equal(nrow(result), nrow(locus_depth))
})

test_that("global pooling checks partial metadata and factor labels", {
  fixture <- pooled_depth_fixture()
  mutations <- fixture$mutations
  mutations$dose <- factor(mutations$dose)
  supplied <- fixture$depth
  supplied$dose <- factor(c("0", "10", "10"), levels = c("10", "0", "unused"))
  expected <- calculate_mf(
    mutations, cols_to_group = c("dose", "visit"), calculate_depth = FALSE,
    precalc_depth_data = fixture$depth, depth_group_cols = "sample"
  )
  expect_equal(calculate_mf(
    mutations, cols_to_group = c("dose", "visit"), calculate_depth = FALSE,
    precalc_depth_data = supplied, depth_group_cols = "sample"
  ), expected)
  supplied$dose[1] <- "10"
  expect_error(calculate_mf(
    mutations, cols_to_group = c("dose", "visit"), calculate_depth = FALSE,
    precalc_depth_data = supplied, depth_group_cols = "sample"
  ), "disagree with mutation_data.*sample1")
  missing_depth <- fixture$depth
  missing_depth$group_depth <- NULL
  expect_error(calculate_mf(
    mutations, cols_to_group = "dose", calculate_depth = FALSE,
    precalc_depth_data = missing_depth, depth_group_cols = "sample"
  ), "group_depth.*finite, non-negative")
  overflowing <- fixture$depth
  overflowing$group_depth <- .Machine$double.xmax
  expect_error(calculate_mf(
    mutations, cols_to_group = "dose", calculate_depth = FALSE,
    precalc_depth_data = overflowing, depth_group_cols = "sample"
  ), "pooled group_depth.*not finite")
})

test_that("explicit sample context keys resolve the dose metadata ambiguity", {
  fixture <- pooled_depth_fixture()
  mutations <- fixture$mutations
  context_depth <- expand.grid(sample = fixture$depth$sample,
                               normalized_ref = c("C", "T"),
                               stringsAsFactors = FALSE)
  context_depth$subtype_depth <- seq_len(nrow(context_depth)) * 100
  context_depth <- dplyr::left_join(
    context_depth, dplyr::distinct(mutations, sample, dose), by = "sample"
  )
  context_depth$group_depth <- ave(context_depth$subtype_depth,
                                  context_depth$sample, FUN = sum)
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  utils::write.table(context_depth, path, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  grouped <- context_depth %>%
    dplyr::group_by(dose, normalized_ref) %>%
    dplyr::summarise(subtype_depth = sum(subtype_depth), .groups = "drop")
  expect_error(calculate_mf(
    mutations, cols_to_group = "dose", subtype_resolution = "base_6",
    calculate_depth = FALSE, precalc_depth_data = context_depth
  ), "finer grouping.*depth_group_cols = 'sample'")
  for (resolution in c("none", "type", "base_6")) {
    for (summarize in c(TRUE, FALSE)) {
      expected <- calculate_mf(
        mutations, cols_to_group = "dose", subtype_resolution = resolution,
        calculate_depth = FALSE, precalc_depth_data = grouped, summary = summarize,
        depth_group_cols = "dose"
      )
      for (input in list(context_depth, path)) {
        expect_equal(calculate_mf(
          mutations, cols_to_group = "dose", subtype_resolution = resolution,
          calculate_depth = FALSE, precalc_depth_data = input, summary = summarize,
          depth_group_cols = "sample"
        ), expected)
      }
    }
  }
  conflict <- context_depth
  conflict$dose[1] <- 99
  expect_error(calculate_mf(
    mutations, cols_to_group = "dose", subtype_resolution = "base_6",
    calculate_depth = FALSE, precalc_depth_data = conflict,
    depth_group_cols = "sample"
  ), "disagree with mutation_data.*sample1")
  expect_error(calculate_mf(
    mutations, cols_to_group = "dose", subtype_resolution = "base_6",
    calculate_depth = FALSE,
    precalc_depth_data = rbind(context_depth, context_depth[1, ]),
    depth_group_cols = "sample"
  ), "duplicate combinations")
})

test_that("explicit depth declaration preserves sample replicates and default calls", {
  fixture <- pooled_depth_fixture()
  mutations <- fixture$mutations
  depth <- fixture$depth
  # Omitted and NULL declarations retain the default grouped interpretation.
  baseline <- calculate_mf(
    mutations, calculate_depth = FALSE, precalc_depth_data = depth,
    retain_metadata_cols = "dose"
  )
  expect_equal(calculate_mf(
    mutations, calculate_depth = FALSE, precalc_depth_data = depth,
    retain_metadata_cols = "dose", depth_group_cols = NULL
  ), baseline)
  expect_equal(calculate_mf(
    mutations, calculate_depth = FALSE, precalc_depth_data = depth,
    retain_metadata_cols = "dose", depth_group_cols = "sample"
  ), baseline)
  expect_equal(nrow(baseline), nrow(depth))
  expect_equal(baseline$group_depth, depth$group_depth)
  expect_equal(baseline$dose, c(0, 10, 10))
  expect_error(calculate_mf(
    mutations, cols_to_group = "dose", calculate_depth = FALSE,
    precalc_depth_data = depth
  ), "Missing columns.*dose")
  supplied <- dplyr::left_join(depth, dplyr::distinct(mutations, sample, dose),
                              by = "sample")
  expect_error(calculate_mf(
    mutations, cols_to_group = "dose", calculate_depth = FALSE,
    precalc_depth_data = supplied
  ), "finer grouping.*depth_group_cols = 'sample'")
  expect_identical(tail(names(formals(calculate_mf)), 1L), "depth_group_cols")
})

test_that("depth declarations validate keys and unsupported pooling", {
  fixture <- pooled_depth_fixture()
  run <- function(keys, depth = fixture$depth, groups = "dose") {
    calculate_mf(fixture$mutations, cols_to_group = groups,
                 calculate_depth = FALSE, precalc_depth_data = depth,
                 depth_group_cols = keys)
  }
  for (keys in list(character(), NA_character_, "", c("sample", "sample"), 1)) {
    expect_error(run(keys), "depth_group_cols must be NULL")
  }
  expect_error(run("absent"), "missing depth_group_cols: absent")
  fine <- fixture$depth
  fine$locus <- "a"
  expect_error(run(c("sample", "locus"), fine), "Summarize the depth table")
  expect_error(calculate_mf(fixture$mutations, depth_group_cols = "sample"),
               "requires calculate_depth = FALSE")
  expect_error(calculate_mf(fixture$mutations, calculate_depth = FALSE,
                           depth_group_cols = "sample"),
               "requires calculate_depth = FALSE and precalc_depth_data")
  data <- fixture$mutations
  data$locus <- rep(c("a", "b"), length.out = nrow(data))
  exact <- dplyr::distinct(data, sample, locus)
  exact$group_depth <- 100
  expect_equal(
    calculate_mf(data, cols_to_group = c("sample", "locus"),
                 calculate_depth = FALSE, precalc_depth_data = exact,
                 depth_group_cols = c("locus", "sample")),
    calculate_mf(data, cols_to_group = c("sample", "locus"),
                 calculate_depth = FALSE, precalc_depth_data = exact)
  )
})

test_that("explicit contexts warn for missing samples without changing counts", {
  fixture <- pooled_depth_fixture()
  contexts <- expand.grid(sample = fixture$depth$sample[-2],
                          normalized_ref = c("C", "T"),
                          stringsAsFactors = FALSE)
  contexts$subtype_depth <- 100
  expect_warning(
    actual <- calculate_mf(
      fixture$mutations, cols_to_group = "dose", subtype_resolution = "base_6",
      calculate_depth = FALSE, precalc_depth_data = contexts,
      depth_group_cols = "sample"
    ), "No precalculated depth.*sample2.*Affected groups.*10"
  )
  # The existing implicit context path remains free of this new warning.
  expect_warning(
    legacy <- calculate_mf(
      fixture$mutations, cols_to_group = "dose", subtype_resolution = "base_6",
      calculate_depth = FALSE, precalc_depth_data = contexts
    ), NA
  )
  expect_equal(actual, legacy)
})
