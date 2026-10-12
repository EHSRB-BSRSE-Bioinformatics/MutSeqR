test_that("plot_radar works with facet = TRUE and facet = FALSE", {
  skip_if_not_installed("fmsb")

  test_df <- data.frame(
    dose_group = rep(c("Control", "Low", "High"), each = 4),
    label = rep(c("C>A", "C>G", "C>T", "T>A"), times = 3),
    mean_mf = c(1e-5, 2e-5, 3e-5, 4e-5,
                2e-5, 3e-5, 4e-5, 5e-5,
                5e-5, 6e-5, 7e-5, 8e-5)
  )

  # Test facet = TRUE
  tmp1 <- tempfile(fileext = ".png")
  png(tmp1)
  expect_silent(
    plot_radar(
      mf_data = test_df,
      response_col = "mean_mf",
      label_col = "label",
      group_col = "dose_group",
      facet = TRUE,
      indiv_y = FALSE
    )
  )
  dev.off()
  unlink(tmp1)

  # Test facet = FALSE with custom_palette and label_colors
  tmp2 <- tempfile(fileext = ".png")
  png(tmp2)
  expect_silent(
    plot_radar(
      mf_data = test_df,
      response_col = "mean_mf",
      label_col = "label",
      group_col = "dose_group",
      facet = FALSE,
      custom_palette = c("Control" = "#45567f", "Low" = "#6a74a7", "High" = "#e2aab9"),
      label_colors = c("C>A" = "red", "C>T" = "blue"),
      plot_legend = TRUE
    )
  )
  dev.off()
  unlink(tmp2)

  # Test backward compatibility with facet_col
  tmp3 <- tempfile(fileext = ".png")
  png(tmp3)
  expect_silent(
    plot_radar(
      mf_data = test_df,
      response_col = "mean_mf",
      label_col = "label",
      facet_col = "dose_group"
    )
  )
  dev.off()
  unlink(tmp3)

  # Test validation
  expect_error(
    plot_radar(
      mf_data = test_df,
      response_col = "nonexistent_col",
      label_col = "label",
      group_col = "dose_group"
    ),
    "Columns not found"
  )
})
