test_that("filter simple mutation data", {
  mutation_data <- readRDS("./testdata/simple_mutation_data.rds")
  # VAF Cutoff
  result <- filter_mut(
    mutation_data = mutation_data,
    vaf_cutoff = 0.01
  )
  idx_t <- which(mutation_data$vaf > 0.01)
  idx_f <- which(mutation_data$vaf <= 0.01)
  expect_true(all(result$filter_mut[idx_t] == TRUE))
  expect_true(all(result$filter_reason[idx_t] == "germline"))
  expect_true(all(result$is_germline[idx_t] == TRUE))
  expect_true(all(result$filter_mut[idx_f] == FALSE))
  expect_true(all(result$filter_reason[idx_f] == ""))
  expect_true(all(result$is_germline[idx_f] == FALSE))
  expect_equal(sum(result$filter_mut), 15)
  expect_equal(sum(result$is_germline), 15)
  # SNV in germ MNV
  #  rm_filtered_mut_from_depth
  result <- filter_mut(
    mutation_data = mutation_data,
    vaf_cutoff = 0.01,
    snv_in_germ_mnv = TRUE,
    rm_filtered_mut_from_depth = TRUE
  )
  expect_equal(sum(result$filter_mut), 16)
  expect_true(result$filter_reason[17] == "snv_in_germ_mnv")
  expect_true(result$filter_reason[19] == "germline|snv_in_germ_mnv")
  expect_true(result$total_depth[17] == 339)
  expect_true(result$total_depth[19] == 180)
  # rm_abnormal_vaf
  result <- filter_mut(
    mutation_data = mutation_data,
    rm_abnormal_vaf = TRUE,
    return_filtered_rows = TRUE
  )
  expect_length(result, 2)
  expect_equal(nrow(result$mutation_data), 31)
  expect_equal(nrow(result$filtered_rows), 5)
  expect_true(unique(result$filtered_rows$filter_reason) == "abnormal_vaf")
  # Custom Filter
  result <- filter_mut(
    mutation_data = mutation_data,
    custom_filter_col = "sd_metadata",
    custom_filter_val = "C",
    custom_filter_rm = FALSE
  )
  idx <- which(mutation_data$sd_metadata == "C")
  expect_true(all(result$filter_mut[idx] == TRUE))
  expect_true(all(result$filter_reason[idx] == "C"))
  
  result <- filter_mut(
    mutation_data = mutation_data,
    custom_filter_col = "sd_metadata",
    custom_filter_val = c("A", "B"),
    custom_filter_rm = TRUE,
    return_filtered_rows = TRUE
  )
  expect_length(result, 2)
  expect_true(unique(result$mutation_data$sd_metadata) == "C")
  expect_setequal(unique(result$filtered_rows$sd_metadata), c("A", "B"))
  # Regions
  regions <- data.frame(contig = c("chr1", "chr2"),
                      start = c(101, 201),
                      end = c(110, 210),
                      rg_metadata = c("R1", "R2"))
  result <- filter_mut(
    mutation_data = mutation_data,
    regions = regions,
    regions_filter = "keep_within",
    allow_half_overlap = FALSE
  )
  expect_equal(nrow(result), 35)
  result <- filter_mut(
    mutation_data = mutation_data,
    regions = regions,
    regions_filter = "keep_within",
    allow_half_overlap = TRUE
  )
  expect_equal(nrow(result), 36)
  result <- filter_mut(
    mutation_data = mutation_data,
    regions = regions,
    regions_filter = "remove_within",
    allow_half_overlap = FALSE
  )
  expect_equal(nrow(result), 1)
})