calculate_total_depths <-
  function (mut_df,
            unique_group_by_columns,
            summarize_by_columns) {
    mut_df %>% group_by_at(unique_group_by_columns) %>% summarise(total_bases = first(.data[[op$column.total_depth]]),
                                                                  .groups = "keep") %>% ungroup() %>% group_by_at(summarize_by_columns) %>%
      summarise(total_bases = sum(.data[["total_bases"]])) %>%
      ungroup()
  }