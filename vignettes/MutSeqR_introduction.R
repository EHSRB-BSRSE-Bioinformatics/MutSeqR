## ----setup, include=FALSE-----------------------------------------------------
BiocStyle::markdown()
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(DT)
library(htmltools)


## ----install-git, eval = FALSE------------------------------------------------
## install.packages("devtools")
## 
## devtools::install_github(
##   "EHSRB-BSRSE-Bioinformatics/MutSeqR",
##   auth_token = "your personal_access_token from GitHub"
## )


## ----install-bioc, eval=FALSE-------------------------------------------------
## if (!require("BiocManager", quietly = TRUE)) {
##   install.packages("BiocManager")
## }
## BiocManager::install("MutSeqR")


## ----load-lib-----------------------------------------------------------------
library(MutSeqR)


## -----------------------------------------------------------------------------
#| label: tbl-required-columns
#| tbl-cap: "Required columns for mutation data import."

cols <- data.frame(
  Column = c("contig", "start", "end", "ref", "alt", "sample", "*SUGGESTED FIELDS*", "alt_depth", "total_depth or depth"),
  `VCF Specification` = c(
    "CHROM", "POS", "INFO(END)", "REF", "ALT", "INFO(sample) or Genotype Header", "", "VD", "AD or DP"
  ),
  Definition = c(
    "The name of the reference sequence.",
    "The start position of the feature. 0-based coordinates are accepted but will be changed to 1-based during import.",
    "The half-open end position of the feature in contig.",
    "The reference allele at this position.",
    "The left-aligned, normalized, alternate allele at this position.",
    "A unique identifier for the sample library. For VCF files, this field may be provided in either the INFO field, or as the header to the GENOTYPE field.",
    "",
    "The read depth supporting the alternate allele. If not included, the function will assume an alt_depth of 1 at variant sites.",
    "The total read depth at this position. This column can be “total_depth” which excludes N-calls, or “depth”, which includes N-calls, if “total_depth” is not available. For VCF files, the total_depth is calculated as the sum of AD. DP is equivalent to \"depth\"."
  )
)

knitr::kable(cols, format = "html", escape = FALSE)


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
DT::datatable(head(imported_example_data), options = list(scrollX = TRUE))


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
DT::datatable(head(imported_example_data), options = list(scrollX = TRUE))


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
DT::datatable(head(imported_example_data), options = list(scrollX = TRUE))


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
DT::datatable(head(imported_example_data_custom), options = list(scrollX = TRUE))


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
DT::datatable(head(filter_table), options = list(scrollX = TRUE))


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
#| tbl-cap: "datatable {#tbl-datatable}"
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
    scrollX = TRUE,
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
    scrollX = TRUE,
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
    scrollX = TRUE,
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
    scrollX = TRUE,
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
    scrollX = TRUE,
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


## ----model-mf-rg, fig.show='hide'---------------------------------------------
# Create a contrasts table for the pairwise comparisons.
combinations <- expand.grid(dose_group = unique(mf_data_rg$dose_group),
                            label = unique(mf_data_rg$label))
combinations <- combinations[combinations$dose_group != "Control", ]
combinations$col1 <- with(combinations, paste(dose_group, label, sep = ":"))
combinations$col2 <- with(combinations, paste("Control", label, sep = ":"))
contrasts2 <- combinations[, c("col1", "col2")]

# Run the model
# To improve convergence of the model, we will supply the control argument
# to the  glmer function
model_by_target <- model_mf(mf_data = mf_data_rg,
  fixed_effects = c("dose_group", "label"),
  test_interaction = TRUE,
  random_effects = "sample",
  muts = "sum_min",
  total_count = "group_depth",
  contrasts = contrasts2,
  reference_level = c("Control", "chr1"),
  control = lme4::glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 2e5))
)


## ----hist-caption2, include=FALSE---------------------------------------------
caption <- paste(
  "GLMM residuals of MFmin modelled as an effect of Dose and Genomic Target.",
  "x is pearson's residuals, y is frequency.",
  "Plotted to validate model assumptions. n = 24."
)


## ----model-rg-hist, echo=FALSE, fig.small=TRUE, fig.cap=caption---------------
model_data <- model_by_target$model_data
hist_data <- hist(model_data$residuals, plot = FALSE)
ylim_max <- max(hist_data$counts) + 1
plot <- hist(model_data$residuals,
  main = "Pearson Residuals",
  col = "yellow",
  ylim = c(0, ylim_max)
)


## ----qq-caption2, include=FALSE-----------------------------------------------
caption <- paste(
  "GLMM residuals of MFmin modelled as an effect of Dose and Genomic Target",
  "expressed as a quantile-quantile plot.",
  "Y is the pearson's residuals of the model in ascending order",
  "x is the quantiles of standard normal distribution for n of 24.",
  "Plotted to validate model assumptions."
)


## ----model-rg-qq, echo=FALSE, fig.small=TRUE, fig.cap=caption-----------------
qqplot <- stats::qqnorm(model_data$residuals, main = "QQ Plot of Residuals")
stats::qqline(model_data$residuals, col = "red")


## ----model-rg-summary---------------------------------------------------------
model_by_target$summary


## ----model-rg-anova-----------------------------------------------------------
model_by_target$anova


## ----model-rg-est, echo=FALSE-------------------------------------------------
DT::datatable(
  model_by_target$point_estimates,
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


## ----model-rg-comp, echo=FALSE------------------------------------------------
DT::datatable(
  model_by_target$pairwise_comparisons,
  options = list(
    scrollX = TRUE,
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


## ----plot-model-cap, include=FALSE--------------------------------------------
caption <- paste(
  "Mean Mutation Frequency Minimum (mutations/bp) per Dose",
  "estimated using a generalized linear model. Error bars are the S.E.M.",
  "Symbols indicate significance differences (p < 0.05)."
)


## ----plot-model-dose, fig.cap=caption-----------------------------------------
plot <- plot_model_mf(
  model_by_dose,
  plot_type = "bar",
  x_effect = "dose",
  plot_error_bars = TRUE,
  plot_signif = TRUE,
  x_order = c("Control", "Low", "Medium", "High"),
  x_label = "Dose Group",
  y_label = "Estimated Mean Mutation Frequency (mutations/bp)"
)
plot


## ----plot-model-cap2, include=FALSE-------------------------------------------
caption <- paste(
  "Mean Mutation Frequency Minimum (mutations/bp) per Genomic Target and Dose",
  "estimated using a generalized linear mixed model. Error bars are the SEM.",
  "Symbols indicate significance differences (p < 0.05) between dose levels",
  "for individual genomic regions."
)


## ----plot-model-rg, fig.wide=TRUE, fig.cap=caption----------------------------
# Define the order of the genomic targets for the x-axis:
# We will order them from lowest to highest MF at the High dose.
label_order <- model_by_target$point_estimates %>%
  dplyr::filter(dose_group == "High") %>%
  dplyr::arrange(Estimate) %>%
  dplyr::pull(label)

# Define the order of the doses for the fill
dose_order <- c("Control", "Low", "Medium", "High")

plot <- plot_model_mf(
  model = model_by_target,
  plot_type = "bar",
  x_effect = "label",
  plot_error_bars = TRUE,
  plot_signif = TRUE,
  ref_effect = "dose_group",
  x_order = label_order,
  fill_order = dose_order,
  x_label = "Target",
  y_label = "Mutation Frequency (mutations/bp)",
  fill_label = "Dose",
  plot_title = "",
  custom_palette = c("#ef476f",
                     "#ffd166",
                     "#06d6a0",
                     "#118ab2")
)
# Rotate the x-axis labels for clarity using ggplot2 functions.
plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
plot


## ----proast, results='hide'---------------------------------------------------
proast_results <- bmd_proast(
  mf_data = mf_data_global,
  dose_col = "dose",
  response_col = c("mf_min", "mf_max"),
  bmr = 0.5,
  model_averaging = TRUE,
  num_bootstraps = 10, # recommended value 200
  plot_results = FALSE
)


## ----proast-results, echo=FALSE-----------------------------------------------
DT::datatable(proast_results) %>%
  DT::formatRound(columns = 3:8, digits = 2)


## ----install-toxicr, eval=FALSE-----------------------------------------------
## devtools::install_github("NIEHS/ToxicR")


## ----toxicr, eval=FALSE-------------------------------------------------------
## toxicr_results <- bmd_toxicr(
##   mf_data = mf_data_global,
##   dose_col = "dose",
##   response_col = c("mf_min", "mf_max"),
##   bmr_type = "rel",
##   bmr = 0.5,
##   model_averaging = TRUE,
##   ma_summary = TRUE,
##   plot_results = FALSE
## )
## toxicr_results


## ----ci-plot-caption, include=FALSE-------------------------------------------
caption <- paste(
  "Benchmark dose with 90% confidence intervals representing the dose at,",
  "which a 50% increase in mutation frequency occurs from reference level.",
  "Calculated using ROAST software. Black points represent the BMD,",
  "red points the BMDL, and blue points, the BMDU"
)


## ----plot-ci, fig.cap=caption-------------------------------------------------
plot_results <- data.frame(
  Response = c("PROAST", "ToxicR"),
  BMD = c(9.111, 9.641894),
  BMDL = c(7.38, 8.032936),
  BMDU = c(10.9, 10.97636)
)
plot <- plot_ci(
  data = plot_results,
  order = "asc",
  x_lab = "Dose (mg/kg-bw/d)",
  y_lab = "BMD Method"
)
plot


## ----spectra-comparison-------------------------------------------------------
# Calculate the MF per dose at the base_6 resolution.
mf_data_6_dose <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "dose_group",
  subtype_resolution = "base_6"
)

# Create the contrast table
contrasts_table <- data.frame(
  col1 = c("Low", "Medium", "High"),
  col2 = c("Control", "Control", "Control")
)
# Run the analysis
ex_spectra_comp <- spectra_comparison(
  mf_data = mf_data_6_dose,
  exp_variable = "dose_group",
  contrasts = contrasts_table
)


## ----spectra-comp-table, echo=FALSE-------------------------------------------
DT::datatable(
  ex_spectra_comp,
  options = list(
    columnDefs = list(
      list(targets = c("p.value", "adj_p.value"),
           render = DT::JS("function(data, type, row, meta) { return data.toExponential(2); }"))
    ),
    rowCallback = DT::JS(
      "function(row, data, dataIndex) {
        $('td:eq(2)', row).css('text-align', 'right'); // Align the content of column 2 to the right
      }"
    )
  )
) %>%
  DT::formatRound(columns = "G2", digits = 2)


## ----install-sigprof-depends, eval=FALSE--------------------------------------
## # Install reticulate
## install.packages("reticulate")
## 
## # Install python
## reticulate::install_python()
## 
## # Install SigProfilerMatrixGeneratorR from github using devtools.
## devtools::install_github("AlexandrovLab/SigProfilerMatrixGeneratorR")
## 


## ----sig-fit, eval=FALSE------------------------------------------------------
## # Run Analysis
## signature_fitting(
##   mutation_data = example_data, # filtered mutation data
##   project_name = "Example",
##   project_genome = "mm10",
##   env_name = "MutSeqR",
##   group = "dose_group",
##   python_version = "3.11",
##   output_path = NULL
## )


## ----mut-calling-file, eval=FALSE---------------------------------------------
## write_mutation_calling_file(
##   mutation_data = example_data,
##   project_name = "Example",
##   project_genome = "mm10",
##   output_path = NULL
## )


## ----write-mut-matrix, eval=FALSE---------------------------------------------
## write_mutational_matrix(
##   mutation_data = example_data,
##   group = "dose_group",
##   subtype_resolution = "base_96",
##   mf_type = "min",
##   output_path = NULL
## )


## ----plot-spectra-cap1,  include=FALSE----------------------------------------
caption <- paste(
  "Mutation spectrum (minimum) per Dose. Subtypes include single-nucleotide,",
  "variants at 6-base resolution, complex variants, deletions, insertions,",
  "multi-nucleotide variants (mnv) and structural variants (sv). Subtypes",
  "are represented by colour. Data is the proportion normalized to",
  "sequencing depth."
)


## ----plot-spectra-dose-p, fig.cap=caption-------------------------------------
# Calculate the mf data at the 6-base resolution for each dose
# We will exclude ambiguous or uncategorized variants since we don't
# have any in this data.
mf_data_6_dose <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "dose_group",
  subtype_resolution = "base_6",
  variant_types = c("-ambiguous", "-uncategorized")
)
# Set the desired order for the dose group:
mf_data_6_dose$dose_group <- factor(
  mf_data_6_dose$dose_group,
  levels = c("Control", "Low", "Medium", "High")
)
# Plot
plot <- plot_spectra(
  mf_data = mf_data_6_dose,
  group_col = "dose_group",
  subtype_resolution = "base_6",
  response = "proportion",
  group_order = "arranged",
  group_order_input = "dose_group",
  x_lab = "Dose Group",
  y_lab = "Subtype Proportion"
)
plot


## ----plot-spectra-cap2, include=FALSE-----------------------------------------
caption <- paste(
  "Mutation spectrum (minimum) per Dose. Subtypes include single-nucleotide,",
  "variants at 6-base resolution, complex variants, deletions, insertions,",
  "multi-nucleotide variants (mnv) and structural variants (sv). Subtypes",
  "are represented by colour. Data is subtype frequency (mutations/bp)."
)


## ----plot-spectra-dose-f, fig.cap=caption-------------------------------------
# Plot
plot <- plot_spectra(
  mf_data = mf_data_6_dose,
  group_col = "dose_group",
  subtype_resolution = "base_6",
  response = "mf",
  group_order = "arranged",
  group_order_input = "dose_group",
  x_lab = "Dose Group",
  y_lab = "Subtype Frequency (mutations/bp)"
)
plot


## ----plot-spectra-cap3,  include=FALSE----------------------------------------
caption <- paste(
  "Mutation spectrum (minimum) per Dose. Subtypes include single-nucleotide,",
  "variants at 6-base resolution, complex variants, deletions, insertions,",
  "multi-nucleotide variants (mnv) and structural variants (sv). Subtypes",
  "are represented by colour. Data is the mutation count."
)


## ----plot-spectra-dose-s------------------------------------------------------
# Plot
plot <- plot_spectra(
  mf_data = mf_data_6_dose,
  group_col = "dose_group",
  subtype_resolution = "base_6",
  response = "sum",
  group_order = "arranged",
  group_order_input = "dose_group",
  x_lab = "Dose Group",
  y_lab = "Subtype Mutation Count"
)
plot


## ----plot-spectra-cap4,  include=FALSE----------------------------------------
caption <- paste(
  "Mutation spectrum (minimum) per Sample. Subtypes include",
  "single-nucleotide variants at 6-base resolution, complex variants,",
  "deletions, insertions, multi-nucleotide variants (mnv) and structural",
  "variants (sv). Subtypes are represented by colour. Data is the proportion",
  "normalized to sequencing depth. Samples are clustered based on the",
  "Euclidean distance between their subtype proportions."
)


## ----plot-spectra-clustered, fig.wide=TRUE, fig.cap=caption-------------------
# Calculate the mf data at the 6-base resolution for each sample
mf_data_6 <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "new_sample_id",
  subtype_resolution = "base_6",
  variant_types = c("-ambiguous", "-uncategorized")
)
# Plot
plot <- plot_spectra(
  mf_data = mf_data_6,
  group_col = "new_sample_id",
  subtype_resolution = "base_6",
  response = "proportion",
  group_order = "clustered",
  x_lab = "Sample",
  y_lab = "Subtype Proportion"
)
plot


## ----trinuc-plot-caption, include=FALSE---------------------------------------
caption <- paste(
  "96-base trinucleotide spectra (minimum) per Dose.",
  "Bars are the proportion of SNV subtypes within their trinucleotide context",
  "normalized to the sequencing depth. ",
  "Bars are coloured based on SNV subtype. Data labels indicate the total",
  "number of mutations for each SNV subtype."
)


## ----plot-trinucleotide, fig.wide=TRUE, fig.show='hide'-----------------------
# Calculate the mf data at the 96-base resolution for each dose
mf_data_96_dose <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "dose_group",
  subtype_resolution = "base_96",
  variant_types = "snv"
)
# Plot
plots <- plot_trinucleotide(
  mf_96 = mf_data_96_dose,
  group_col = "dose_group",
  response = "proportion",
  mf_type = "min",
  output_path = NULL
)


## ----trinuc-C, fig.wide=TRUE, echo=FALSE, fig.cap=caption---------------------
plots[1]


## ----trinuc-L, fig.wide=TRUE, echo=FALSE, fig.cap=caption---------------------
plots[2]


## ----trinuc-M, fig.wide=TRUE, echo=FALSE, fig.cap=caption---------------------
plots[3]


## ----trinuc-H, fig.wide=TRUE, echo=FALSE, fig.cap=caption---------------------
plots[4]


## ----heatmap-cap, include=FALSE-----------------------------------------------
caption <- paste(
  "96-base trinucleotide spectra (minimum) per Sample, facetted by Dose.",
  "Colour represents the proportion of SNV subtypes within their",
  "trinucleotide context normalized to the sequencing depth."
)


## ----plot-heatmap, fig.wide=TRUE, fig.cap=caption-----------------------------
# Calculate the mf data at the 96-base resolution for each sample
mf_data_96 <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = "sample",
  subtype_resolution = "base_96",
  variant_types = "snv",
  retain_metadata_cols = "dose_group"
)
mf_data_96$dose_group <- factor(
  mf_data_96$dose_group,
  levels = c("Control", "Low", "Medium", "High")
)
# Plot
plot <- plot_trinucleotide_heatmap(
  mf_data = mf_data_96,
  group_col = "sample",
  facet_col = "dose_group"
)
plot


## ----bubbles-cap, include=FALSE-----------------------------------------------
caption <- paste(
  "Multiplet mutations plotted per Dose. Each circle represents a mutation,",
  "coloured by mutation subtype. The size of the circle is scaled by the",
  "mutation's alternative depth."
)


## ----bubble-plot, fig.wide=TRUE, fig.cap=caption------------------------------
plot <- plot_bubbles(
  mutation_data = example_data,
  size_by = "alt_depth",
  facet_col = "dose_group",
  color_by = "normalized_subtype"
)
plot


## ----get-seq------------------------------------------------------------------
regions_seq <- get_seq(regions = "TSpanel_mouse")
regions_seq


## ----get-seq-custom-----------------------------------------------------------
# We will load the TSpanel_human regions file as an example
human <- load_regions_file("TSpanel_human")
regions_seq <- get_seq(regions = human,
                       is_0_based_rg = FALSE,
                       species = "human",
                       genome = "hg38",
                       masked = FALSE,
                       padding = 0)
regions_seq


## ----write-ref-fasta, eval=FALSE----------------------------------------------
## write_reference_fasta(regions_seq, output_path = NULL)


## ----write-excel, eval=FALSE--------------------------------------------------
## # save a single data frame to an Excel file
## write_excel(mf_data_global, workbook_name = "example_mf_data")
## 
## # Write multiple data frames to a list to export all at once.
## list <- list(mf_data_global, mf_data_rg, mf_data_6)
## names(list) <- c("mf_per_sample", "mf_per_region", "pf_6spectra")
## 
## #save a list of data frames to an Excel file
## write_excel(list, workbook_name = "example_mf_data_list")
## 


## ----write-excel-model, eval=FALSE--------------------------------------------
## write_excel(
##   model_by_dose,
##   workbook_name = "Example_model",
##   model_results = TRUE
## )


## ----write-vcf, eval=FALSE----------------------------------------------------
## write_vcf_from_mut(example_data)


## ----download-yaml, eval=FALSE------------------------------------------------
## config <- system.file("extdata", "inputs", "summary_config.yaml", package = "MutSeqR")
## file.copy(from = config, to = "path/to/summary_config.yaml")


## ----render-report, eval=FALSE------------------------------------------------
## render_report(config_file_path = "path/to/summary_config.yaml",
##               output_file = "path/to/output_file.html",
##               output_format = "html_document")


## ----session-info, echo=FALSE-------------------------------------------------
sessionInfo()

