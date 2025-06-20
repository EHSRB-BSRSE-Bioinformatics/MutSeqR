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
  expected_prop_min <- sapply(seq_len(nrow(mf_data)), function(i) {
    sample <- mf_data$sample[i]
    mf_data$mf_min[i] / sum(mf_data$mf_min[mf_data$sample == sample])
  })
  expect_equal(mf_data$proportion_min, expected_prop_min)
  expected_prop_max <- sapply(seq_len(nrow(mf_data)), function(i) {
    sample <- mf_data$sample[i]
    mf_data$mf_max[i] / sum(mf_data$mf_max[mf_data$sample == sample])
  })
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
  expected_prop_min <- sapply(seq_len(nrow(mf_data)), function(i) {
    sample <- mf_data$sample[i]
    mf_data$mf_min[i] / sum(mf_data$mf_min[mf_data$sample == sample])
  })
  expect_equal(mf_data$proportion_min, expected_prop_min)
  expected_prop_max <- sapply(seq_len(nrow(mf_data)), function(i) {
    sample <- mf_data$sample[i]
    mf_data$mf_max[i] / sum(mf_data$mf_max[mf_data$sample == sample])
  })
  expect_equal(mf_data$proportion_max, expected_prop_max)
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