## ----setup, include=FALSE-----------------------------------------------------
BiocStyle::markdown()
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(DT)
library(bs4Dash)
library(htmltools)

## ----install-git, eval = FALSE------------------------------------------------
# install.packages("devtools")
# 
# devtools::install_github(
#   "EHSRB-BSRSE-Bioinformatics/MutSeqR",
#   auth_token = "your personal_access_token from GitHub"
# )

## ----install-bioc, eval=FALSE-------------------------------------------------
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("MutSeqR")

## ----load-lib-----------------------------------------------------------------
library(MutSeqR)

## ----import-vcf---------------------------------------------------------------
example_file <- system.file(
  "extdata",
  "Example_files",
  "example_import_vcf_data_cleaned.vcf.bgz",
  package = "MutSeqR"
)
sample_metadata <- data.frame(
  sample = "dna00996.1",
  dose = "50",
  dose_group = "High"
)
# Import the data
imported_example_data <- import_vcf_data(
  vcf_file = example_file,
  sample_data = sample_metadata,
  genome = "mm10",
  species = "mouse",
  masked_BS_genome = FALSE
)

## ----import-vcf-table, echo=FALSE---------------------------------------------
DT::datatable(head(imported_example_data))

## ----import-mut---------------------------------------------------------------
example_file <- system.file(
  "extdata",
  "Example_files",
  "example_import_mut_data.rds",
  package = "MutSeqR"
)
example_data <- readRDS(example_file)

sample_metadata <- data.frame(
  sample = "dna00996.1",
  dose = "50",
  dose_group = "High"
)
# Import the data
imported_example_data <- import_mut_data(
  mut_file = example_data,
  sample_data = sample_metadata,
  genome = "mm10",
  species = "mouse",
  masked_BS_genome = FALSE,
  is_0_based_mut = TRUE # indicates that the genomic coordinates are 0-based.
# Coordinates will be changed to 1-based upon import.
)

## ----import-mut-table, echo=FALSE---------------------------------------------
DT::datatable(head(imported_example_data))

## ----import-mut-rg------------------------------------------------------------
imported_example_data <- import_mut_data(
  mut_file = example_data,
  sample_data = sample_metadata,
  genome = "mm10",
  species = "mouse",
  masked_BS_genome = FALSE,
  is_0_based_mut = TRUE,
  regions = "TSpanel_mouse"
)

## ----import-mut-table-rg, echo=FALSE------------------------------------------
DT::datatable(head(imported_example_data))

## ----load-rg------------------------------------------------------------------
region_example <- load_regions_file("TSpanel_mouse")
region_example

## ----import-mut-custom-cols, message=TRUE-------------------------------------
mut_data <- system.file(
  "extdata", "Example_files",
  "example_import_mut_data_custom_col_names.txt",
  package = "MutSeqR"
)
imported_example_data_custom <- import_mut_data(
  mut_file = mut_data,
  custom_column_names = list(my_contig_name = "contig",
                             my_sample_name = "sample")
)

## ----import-mut-cust-table, echo=FALSE----------------------------------------
DT::datatable(head(imported_example_data_custom))

## ----filter-mut, message=TRUE-------------------------------------------------
# load the example data
example_file <- system.file(
  "extdata", "Example_files",
  "example_mutation_data.rds",
  package = "MutSeqR"
)
example_data <- readRDS(example_file)

# Filter
filtered_example_mutation_data <- filter_mut(
  mutation_data = example_data,
  correct_depth = TRUE,
  vaf_cutoff = 0.01,
  regions = "TSpanel_mouse",
  regions_filter = "keep_within",
  custom_filter_col = "filter",
  custom_filter_val = "EndRepairFillInArtifact",
  custom_filter_rm = FALSE,
  snv_in_germ_mnv = TRUE,
  rm_filtered_mut_from_depth = TRUE,
  return_filtered_rows = FALSE
)

## ----filter-mut-table, echo=FALSE---------------------------------------------
filter_table <- filtered_example_mutation_data %>%
  dplyr::filter(filter_mut == TRUE)
DT::datatable(head(filter_table))

## ----mf-global----------------------------------------------------------------
# load example data:
example_file <- system.file(
  "extdata", "Example_files",
  "example_mutation_data_filtered.rds",
  package = "MutSeqR"
)
example_data <- readRDS(example_file)

mf_data_global <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "none",
  retain_metadata_cols = c("dose_group", "dose")
)

## ----mf-global-table, echo=FALSE----------------------------------------------
DT::datatable(mf_data_global,
  options = list(
    columnDefs = list(
      list(targets = c("mf_min", "mf_max"),
           render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }"))
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
)

## ----mf-rg--------------------------------------------------------------------
mf_data_rg <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = c("sample", "label"),
  subtype_resolution = "none",
  retain_metadata_cols = "dose_group"
)

## ----mf-rg-table, echo=FALSE--------------------------------------------------
DT::datatable(mf_data_rg,
  options = list(
    columnDefs = list(
      list(targets = c("mf_min", "mf_max"),
           render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }"))
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
)

## ----mean-mf------------------------------------------------------------------
mf_data_global <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "none",
  retain_metadata_cols = c("dose_group", "dose")
)
mean_mf <- mf_data_global %>%
  dplyr::group_by(dose_group) %>%
  dplyr::summarise(mean_mf_min = mean(mf_min),
                   SE = sd(mf_min) / sqrt(dplyr::n()))

## ----mf-mean-table, echo=FALSE------------------------------------------------
DT::datatable(mean_mf,
  options = list(
    columnDefs = list(
      list(targets = c("mean_mf_min", "SE"),
           render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }"))
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
)

## ----mf-6---------------------------------------------------------------------
mf_data_6 <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "base_6"
)

## ----mf-6-table, echo=FALSE---------------------------------------------------
DT::datatable(mf_data_6,
  options = list(
    columnDefs = list(
      list(targets = c("mf_min", "mf_max"),
           render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }"))
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
)

## ----mf-indels----------------------------------------------------------------
mf_data_global_indels <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "none",
  variant_types = c("insertion", "deletion")
)

## ----mf-indels-table, echo=FALSE----------------------------------------------
DT::datatable(mf_data_global_indels,
  options = list(
    columnDefs = list(
      list(targets = c("mf_min", "mf_max"),
           render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }"))
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
)

## ----mf-types-----------------------------------------------------------------
mf_data_types <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "type",
  variant_types = c("-ambiguous", "-uncategorized")
)

## ----mf-types-table, echo=FALSE-----------------------------------------------
DT::datatable(mf_data_types,
  options = list(
    columnDefs = list(
      list(targets = c("mf_min", "mf_max"),
           render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }"))
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
)

## ----mf-96--------------------------------------------------------------------
mf_data_96 <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "base_96",
  variant_types = "snv"
)

## ----mf-96_tables, echo=FALSE-------------------------------------------------
DT::datatable(mf_data_types,
  options = list(
    columnDefs = list(
      list(targets = c("mf_min", "mf_max"),
           render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }"))
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
)

## ----precalc1-depth-----------------------------------------------------------
sample_depth <- data.frame(
  sample = unique(example_data$sample),
  group_depth = c(565395266, 755574283, 639909215, 675090988, 598104021,
                  611295330, 648531765, 713240735, 669734626, 684951248,
                  716913381, 692323218, 297661400, 172863681, 672259724,
                  740901132, 558051386, 733727643, 703349287, 884821671,
                  743311822, 799605045, 677693752, 701163532)
)

DT::datatable(sample_depth)

## ----mf-precalc-global--------------------------------------------------------
mf_data_global_precalc <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "none",
  calculate_depth = FALSE,
  precalc_depth_data = sample_depth
)

## ----mf-precalc1-table, echo=FALSE--------------------------------------------
DT::datatable(mf_data_global_precalc,
  options = list(
    columnDefs = list(
      list(targets = c("mf_min", "mf_max"),
           render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }"))
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
)

## ----mf-precalc-6-------------------------------------------------------------
mf_data_6_precalc <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "base_6",
  calculate_depth = FALSE,
  precalc_depth = system.file("extdata", "Example_files",
                              "precalc_depth_base_6_example.txt",
                              package = "MutSeqR")
)

## ----precalc2-depth, echo=FALSE-----------------------------------------------
depth <- read.table(
  file = system.file(
    "extdata", "Example_files",
    "precalc_depth_base_6_example.txt",
    package = "MutSeqR"
  ),
  header = TRUE
)

DT::datatable(depth)

## ----mf-precalc2-table, echo=FALSE--------------------------------------------
DT::datatable(mf_data_6_precalc,
  options = list(
    columnDefs = list(
      list(targets = c("mf_min", "mf_max"),
           render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }"))
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
)

## ----plot-mf-caption, include=FALSE-------------------------------------------
caption <- paste(
  "Mutation Frequency (MF) Minimum and Maximum (mutations/bp) per Sample.",
  "Light colored bars represent MFmin and dark coloured bars represent MFmax.",
  "Bars are coloured and grouped by dose.",
  "Data labels are the number of mutations per sample."
)

## ----plot-mf, fig.width=12, fig.height=6, fig.cap=caption---------------------
# Define the order for dose groups
mf_data_global$dose_group <- factor(
  mf_data_global$dose_group,
  levels = c("Control", "Low", "Medium", "High")
)
plot <- plot_mf(
  mf_data = mf_data_global,
  group_col = "sample",
  plot_type = "bar",
  mf_type = "both",
  fill_col = "dose_group",
  group_order = "arranged",
  group_order_input = "dose_group"
)
plot

## ----plot-mean-mf-caption, include=FALSE--------------------------------------
caption <- paste(
  "Mean Mutation Frequency (MF) Minimum per Dose. Lines are mean ± S.E.M.",
  "Points are individual samples, coloured by dose."
)

## ----plot-mean-mf, fig.cap=caption--------------------------------------------
plot_mean <- plot_mean_mf(
  mf_data = mf_data_global,
  group_col = "dose_group",
  mf_type = "min",
  fill_col = "dose_group",
  add_labels = "none",
  group_order = "arranged",
  group_order_input = "dose_group",
  plot_legend = FALSE,
  x_lab = "Dose Group"
)
plot_mean

## ----model-mf-dose, fig.show='hide'-------------------------------------------
# Create a contrasts table for pairwise comparisons
contrasts_table <- data.frame(
  col1 = c("Low", "Medium", "High"),
  col2 = c("Control", "Control", "Control")
)
# Run the model
model_by_dose <- model_mf(
  mf_data = mf_data_global,
  fixed_effects = "dose_group",
  muts = "sum_min",
  total_count = "group_depth",
  reference_level = "Control",
  contrasts = contrasts_table
)

## ----hist-caption1, include=FALSE---------------------------------------------
caption <- paste(
  "GLM residuals of MFmin modelled as an effect of Dose.",
  "x is pearson's residuals, y is frequency.",
  "Plotted to validate model assumptions. n = 24."
)

## ----model-dose-hist, echo=FALSE, fig.small=TRUE, fig.cap=caption-------------
model_data <- model_by_dose$model_data
hist_data <- hist(model_data$residuals, plot = FALSE)
ylim_max <- max(hist_data$counts) + 1
plot <- hist(model_data$residuals,
  main = "Pearson Residuals",
  col = "yellow",
  ylim = c(0, ylim_max)
)

## ----qq-caption1, include=FALSE-----------------------------------------------
caption <- paste(
  "GLM residuals of MFmin modelled as an effect of Dose",
  "expressed as a quantile-quantile plot.",
  "Y is the pearson's residuals of the model in ascending order",
  "x is the quantiles of standard normal distribution for n of 24.",
  "Plotted to validate model assumptions."
)

## ----model-dose-qq, echo=FALSE, fig.small=TRUE, fig.cap=caption---------------
qqplot <- stats::qqnorm(model_data$residuals, main = "QQ Plot of Residuals")
stats::qqline(model_data$residuals, col = "red")

## ----model-dose-summary-------------------------------------------------------
model_by_dose$summary

## ----model-dose-est, echo=FALSE-----------------------------------------------
DT::datatable(
  model_by_dose$point_estimates,
  options = list(
    columnDefs = list(
      list(targets = 1:4,
           render = DT::JS("function(data, type, row, meta) {
            return data.toExponential(2); 
          }"))
    ),
    rowCallback = DT::JS("function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }")
  )
)

## ----model-dose-comp, echo=FALSE----------------------------------------------
DT::datatable(
  model_by_dose$pairwise_comparisons,
  options = list(
    columnDefs = list(
      list(targets = c("p.value", "adj_p.value"),
        render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }")
      )
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
      $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
) %>%
  DT::formatRound(columns = c(1:3, 6:7), digits = 2)

