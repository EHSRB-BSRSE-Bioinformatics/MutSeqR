test_that("signature_fitting dry run formats data without Python", {
  mutation_data <- readRDS(system.file(
    "extdata", "Example_files", "filtered_simple_mutation_data.rds",
    package = "MutSeqR"
  ))

  result <- signature_fitting(
    mutation_data = mutation_data,
    project_name = "Example",
    project_genome = "mm10",
    group = "sample",
    dry_run = TRUE
  )

  expect_s3_class(result, "data.frame")
  expect_named(result, c(
    "Project", "Samples", "ID", "Genome", "mut_type", "chrom",
    "pos_start", "pos_end", "ref", "alt", "Type"
  ))
  expect_true(nrow(result) > 0L)
  expect_true(all(result$mut_type == "SNP"))
  expect_true(all(result$Type == "SOMATIC"))
})

test_that("signature_fitting validates dry-run arguments", {
  expect_error(
    signature_fitting(data.frame(), dry_run = NA),
    "dry_run must be TRUE or FALSE"
  )
  expect_error(
    signature_fitting(data.frame()),
    "python_version is required"
  )
})
