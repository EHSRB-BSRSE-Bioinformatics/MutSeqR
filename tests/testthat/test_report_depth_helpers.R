test_that("base96 depth data derives global and base6 report depths", {
  base96 <- data.frame(
    sample = "sample1",
    normalized_context = MutSeqR::context_list$base_96,
    subtype_depth = seq_along(MutSeqR::context_list$base_96)
  )

  prepared <- MutSeqR:::prepare_report_depth_data(base96)

  expect_equal(prepared$global$group_depth, sum(base96$subtype_depth))
  expect_equal(prepared$base96$group_depth,
               rep(sum(base96$subtype_depth), nrow(base96)))
  expect_equal(prepared$base6$normalized_ref, c("C", "T"))
  expect_equal(prepared$base6$subtype_depth, c(232, 296))
  expect_equal(prepared$base6$group_depth,
               rep(sum(base96$subtype_depth), 2))
})

test_that("explicit report depths take precedence independently", {
  base96 <- data.frame(
    sample = "sample1",
    normalized_context = MutSeqR::context_list$base_96,
    subtype_depth = rep(1, 32),
    group_depth = rep(700, 32)
  )
  global <- data.frame(sample = "sample1", group_depth = 900)
  base6 <- data.frame(
    sample = c("sample1", "sample1"),
    normalized_ref = c("C", "T"),
    subtype_depth = c(123, 234),
    group_depth = c(901, 901)
  )

  prepared <- MutSeqR:::prepare_report_depth_data(base96, global, base6)
  expect_equal(prepared$global$group_depth, 900)
  expect_equal(prepared$base6$subtype_depth, c(123, 234))
  expect_equal(prepared$base96$group_depth, rep(700, 32))

  base6_no_group <- base6[c("sample", "normalized_ref", "subtype_depth")]
  prepared <- MutSeqR:::prepare_report_depth_data(base96, global, base6_no_group)
  expect_equal(prepared$base6$subtype_depth, c(123, 234))
  expect_equal(prepared$base6$group_depth, c(900, 900))
})

test_that("base96 depth data rejects malformed context tables", {
  valid <- data.frame(
    sample = "sample1",
    normalized_context = MutSeqR::context_list$base_96,
    subtype_depth = rep(1, 32)
  )

  expect_error(
    MutSeqR:::prepare_report_depth_data(valid[-3]),
    "missing required columns"
  )
  expect_error(
    MutSeqR:::prepare_report_depth_data(valid[0, ]),
    "at least one sample"
  )
  invalid_depth <- valid
  invalid_depth$subtype_depth[1] <- Inf
  expect_error(MutSeqR:::prepare_report_depth_data(invalid_depth),
               "finite, non-negative")
  invalid_depth$subtype_depth[1] <- -1
  expect_error(MutSeqR:::prepare_report_depth_data(invalid_depth),
               "finite, non-negative")
  unknown <- valid
  unknown$normalized_context[1] <- "NNN"
  expect_error(MutSeqR:::prepare_report_depth_data(unknown),
               "unknown normalized_context")
  expect_error(MutSeqR:::prepare_report_depth_data(valid[-1, ]),
               "all 32 normalized_context")
  duplicate <- rbind(valid, valid[1, ])
  expect_error(MutSeqR:::prepare_report_depth_data(duplicate),
               "duplicate combinations")
  missing_context <- valid[-1, ]
  expect_error(MutSeqR:::prepare_report_depth_data(missing_context),
               "all 32 normalized_context")
})

test_that("derived report depths work with calculate_mf at each resolution", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  samples <- unique(mutation_data$sample)
  base96 <- expand.grid(
    sample = samples,
    normalized_context = MutSeqR::context_list$base_96,
    stringsAsFactors = FALSE
  )
  base96$subtype_depth <- rep(seq_len(32), each = length(samples))
  prepared <- MutSeqR:::prepare_report_depth_data(base96)

  global_mf <- calculate_mf(
    mutation_data, cols_to_group = "sample", subtype_resolution = "none",
    calculate_depth = FALSE, precalc_depth_data = prepared$global
  )
  base6_mf <- calculate_mf(
    mutation_data, cols_to_group = "sample", subtype_resolution = "base_6",
    calculate_depth = FALSE, precalc_depth_data = prepared$base6
  )
  base96_mf <- calculate_mf(
    mutation_data, cols_to_group = "sample", subtype_resolution = "base_96",
    variant_types = "snv", calculate_depth = FALSE,
    precalc_depth_data = prepared$base96
  )

  expect_true(all(c("mf_min", "mf_max") %in% names(global_mf)))
  expect_true(all(c("mf_min", "mf_max") %in% names(base6_mf)))
  expect_true(all(c("mf_min", "mf_max") %in% names(base96_mf)))
})
