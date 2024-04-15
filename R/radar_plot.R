#' Create a radar plot

# radar chart
radar_plot <- function(data,
                      group){

plot_data <- MF_per_target %>%
  group_by(label, treatment) %>%
  summarise(MF = mean(sample_label_MF_unique), Std.Err = sd(sample_label_MF_unique)/sqrt(n()))

radar_data <- plot_data %>%
  select(label, treatment, MF) %>%
  pivot_wider(names_from = label, values_from = MF)
View(radar_data)
treatments <- radar_data$treatment  
View(radar_data)
plots <- list()
for (i in seq_along(treatments)) {
 treat = treatments[i]
 df_i <- radar_data %>%
  dplyr::filter(treatment == treat) %>%
  select(-treatment)

count <- ncol(df_i)
# Add rows for max and min values
max <- max(df_i, na.rm = TRUE) * 1.1
max_row <- rep(max, count)
min_row <- rep(0, count)
df <- rbind(max_row, min_row, df_i)
axis_labels <- seq(from = 0, to = max, length.out = 5)
axis_labels <- sprintf("%.1e", axis_labels)
title <- paste(treat)
plot <- radarchart(df = df,
            axistype = 1, # 1: center axis, 2: circle axis, 3 both center and around
            caxislabels = axis_labels,
            sep = 5, # number of axis ticks
            caxiscol = "grey",
            vlabels = NULL, # variable labels
            axislabcol = "grey",
            vlcex = 0.7,  # variable label font size
            title = title,
            pcol = "black", # color of the polygon
            pfcol = NULL, # color of the polygon fill
            plwd = 2, # width of the polygon line
            plty = 1, # type of the polygon line
            cglcol = "grey",
            cglty = 1,
            cglwd = 0.7,
            )
plots[[i]] <- plot
}
}