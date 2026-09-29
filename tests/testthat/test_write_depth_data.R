test_that("base6 depth data derives global report depth", {
  base6 <- data.frame(
    sample = c("sample1", "sample1", "sample2", "sample2"),
    normalized_ref = rep(c("C", "T"), 2),
    subtype_depth = c(0, 234, 45, 67)
  )

  prepared <- write_depth_data(base6, d_sep = "\t")

  expect_equal(prepared$global$sample, c("sample1", "sample2"))
  expect_equal(prepared$global$group_depth, c(234, 112))
  expect_equal(prepared$base6$sample, base6$sample)
  expect_equal(prepared$base6$normalized_ref, base6$normalized_ref)
  expect_equal(prepared$base6$subtype_depth, base6$subtype_depth)
  expect_equal(prepared$base6$group_depth, c(234, 234, 112, 112))
  expect_null(prepared$base12)
  expect_null(prepared$base96)
})

test_that("write_depth_data accepts arbitrary groups without sample", {
  grouped_depth <- data.frame(
    dose_group = c("control", "control", "treated", "treated"),
    normalized_ref = rep(MutSeqR::context_list$base_6, 2),
    subtype_depth = c(10, 20, 30, 40)
  )

  prepared <- write_depth_data(
    grouped_depth,
    d_sep = "\t",
    group_cols = "dose_group"
  )

  expect_equal(prepared$global$dose_group, c("control", "treated"))
  expect_equal(prepared$global$group_depth, c(30, 70))
  expect_equal(prepared$base6$dose_group, grouped_depth$dose_group)
  expect_equal(prepared$base6$normalized_ref, grouped_depth$normalized_ref)
  expect_equal(prepared$base6$subtype_depth, grouped_depth$subtype_depth)
  expect_equal(prepared$base6$group_depth, c(30, 30, 70, 70))
  expect_false("sample" %in% names(prepared$global))
})

test_that("base12 depth data derives global and base6 report depths", {
  base12 <- data.frame(
    sample = rep("sample1", 4),
    short_ref = MutSeqR::context_list$base_12,
    subtype_depth = c(10, 25, 30, 40)
  )

  prepared <- write_depth_data(base12, d_sep = "\t")

  expect_equal(prepared$base6$normalized_ref, c("C", "T"))
  expect_equal(prepared$base6$subtype_depth, c(55, 50))
  expect_equal(prepared$global$group_depth, 105)
})
test_that("base96 depth data derives global and base6 report depths", {
  base96 <- data.frame(
    sample = "sample1",
    normalized_context = MutSeqR::context_list$base_96,
    subtype_depth = seq_along(MutSeqR::context_list$base_96)
  )
  prepared <- suppressWarnings(
    MutSeqR::write_depth_data(base96, d_sep = "\t")
  )

  expect_equal(prepared$global$group_depth, sum(base96$subtype_depth))
  expect_equal(prepared$base96$group_depth,
               rep(sum(base96$subtype_depth), nrow(base96)))
  expect_equal(prepared$base6$normalized_ref, c("C", "T"))
  expect_equal(prepared$base6$subtype_depth, c(232, 296))
  expect_equal(prepared$base6$group_depth,
               rep(sum(base96$subtype_depth), 2))
})
test_that("base192 depth data derives every lower resolution", {
  contexts <- MutSeqR::context_list$base_192
  base192 <- data.frame(
    sample = rep("sample1", length(contexts)),
    context = contexts,
    subtype_depth = seq_along(contexts)
  )

  prepared <- MutSeqR::write_depth_data(base192, d_sep = "\t")

  expect_named(prepared, c("base192", "base96", "base12", "base6", "global"))
  expect_equal(nrow(prepared$base192), 64)
  expect_equal(nrow(prepared$base96), 32)
  expect_equal(nrow(prepared$base12), 4)
  expect_equal(nrow(prepared$base6), 2)
  expect_equal(prepared$global$group_depth, sum(base192$subtype_depth))
  expect_equal(sum(prepared$base96$subtype_depth), sum(base192$subtype_depth))
  expect_equal(as.numeric(prepared$base96[1, 3]), 65) # ACA + TGT
  expect_equal(prepared$base12$subtype_depth, c(424, 488, 552, 616))
  expect_equal(prepared$base6$subtype_depth, c(1040, 1040))
})
test_that("write_depth_data reads and validates file paths", {
  base6 <- data.frame(
    sample = c("sample1", "sample1"),
    normalized_ref = MutSeqR::context_list$base_6,
    subtype_depth = c(10, 20)
  )
  depth_file <- tempfile(fileext = ".tsv")
  empty_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(depth_file, empty_file)), add = TRUE)
  utils::write.table(base6, depth_file, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  file.create(empty_file)

  from_file <- MutSeqR::write_depth_data(depth_file, d_sep = "\t")
  from_data_frame <- MutSeqR::write_depth_data(base6, d_sep = "\t")
  expect_equal(from_file, from_data_frame)
  expect_error(
    MutSeqR::write_depth_data(paste0(depth_file, "_missing"), d_sep = "\t"),
    "file does not exist"
  )
  expect_error(
    MutSeqR::write_depth_data(empty_file, d_sep = "\t"),
    "file is empty"
  )
})
test_that("write_depth_data validates columns, values, and contexts", {
  valid <- data.frame(
    sample = rep("sample1", 2),
    normalized_ref = MutSeqR::context_list$base_6,
    subtype_depth = c(1, 2)
  )

  missing_column <- valid[-1]
  expect_error(
    MutSeqR::write_depth_data(missing_column, d_sep = "\t"),
    "missing required columns"
  )
  missing_value <- valid
  missing_value$subtype_depth[1] <- NA_real_
  expect_error(
    MutSeqR::write_depth_data(missing_value, d_sep = "\t"),
    "must not contain any NA values"
  )
  overflow <- valid
  overflow$subtype_depth <- c(1e308, 1e308)
  expect_error(
    MutSeqR::write_depth_data(overflow, d_sep = "\t"),
    "overflowed for group.*sample1"
  )
  invalid_context <- valid
  invalid_context$normalized_ref[1] <- "A"
  expect_error(
    MutSeqR::write_depth_data(invalid_context, d_sep = "\t"),
    "unknown normalized_ref"
  )
  duplicate <- rbind(valid, valid[1, ])
  expect_error(
    MutSeqR::write_depth_data(duplicate, d_sep = "\t"),
    "duplicate combinations"
  )
  incomplete <- valid[-1, ]
  expect_error(
    MutSeqR::write_depth_data(incomplete, d_sep = "\t"),
    "sample1: C"
  )
  multiple_resolutions <- transform(valid, short_ref = normalized_ref)
  expect_error(
    MutSeqR::write_depth_data(multiple_resolutions, d_sep = "\t"),
    "exactly one resolution context column"
  )
})
test_that("derived report depths work with calculate_mf at each resolution", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  samples <- unique(mutation_data$sample)
  base192 <- expand.grid(
    sample = samples,
    context = MutSeqR::context_list$base_192,
    stringsAsFactors = FALSE
  )
  base192$subtype_depth <- seq_len(nrow(base192))
  prepared <- MutSeqR::write_depth_data(base192, d_sep = "\t")

  global_mf <- calculate_mf(
    mutation_data, cols_to_group = "sample", subtype_resolution = "none",
    calculate_depth = FALSE, precalc_depth_data = prepared$global
  )
  base6_mf <- calculate_mf(
    mutation_data, cols_to_group = "sample", subtype_resolution = "base_6",
    calculate_depth = FALSE, precalc_depth_data = prepared$base6
  )
  base12_mf <- calculate_mf(
    mutation_data, cols_to_group = "sample", subtype_resolution = "base_12",
    calculate_depth = FALSE, precalc_depth_data = prepared$base12
  )
  base96_mf <- calculate_mf(
    mutation_data, cols_to_group = "sample", subtype_resolution = "base_96",
    variant_types = "snv", calculate_depth = FALSE,
    precalc_depth_data = prepared$base96
  )
  base192_mf <- calculate_mf(
    mutation_data, cols_to_group = "sample", subtype_resolution = "base_192",
    variant_types = "snv", calculate_depth = FALSE,
    precalc_depth_data = prepared$base192
  )

  expect_true(all(c("mf_min", "mf_max") %in% names(global_mf)))
  expect_true(all(c("mf_min", "mf_max") %in% names(base6_mf)))
  expect_true(all(c("mf_min", "mf_max") %in% names(base12_mf)))
  expect_true(all(c("mf_min", "mf_max") %in% names(base96_mf)))
  expect_true(all(c("mf_min", "mf_max") %in% names(base192_mf)))
})
