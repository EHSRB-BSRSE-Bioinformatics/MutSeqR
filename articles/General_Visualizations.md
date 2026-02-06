# MutSeqR: General Visualization

## Load MutSeqR and Example Data

``` r
library(MutSeqR)
```

The example data consists of 24 mouse bone marrow DNA samples imported
using import_mut_data() and filtered with filter_mut. Data was sequenced
on the TS Mouse Mutagenesis Panel. Example data is retrieved from
MutSeqRData, an ExperimentHub data package.

``` r
library(ExperimentHub)
# load the index
eh <- ExperimentHub()
```

## Bubble Plots

`plot_bubbles` is used to visually represent the distribution and
density of mutations that are observed in multiple reads (multiplets).
Each mutation in a given group is represented by a bubble whose size is
scaled on either the `alt_depth` or the `vaf`. Thus mutations that occur
in a large number of reads are represented by a large bubble. These
plots make it easy to determine if MFmax is driven by a few highly
expanded mutations versus serveral moderately expanded mutations.

Plots can be facetted by user-defined groups, and bubbles can be
coloured by any variable of interest to help discern patterns in
mutation recurrence.

*Example 1. Plot mutations per dose group, bubbles coloured by base-6
subtype*

``` r
# load example data:
example_data <- eh[["EH9861"]]

plot <- plot_bubbles(
  mutation_data = example_data,
  size_by = "alt_depth",
  facet_col = "dose_group",
  color_by = "normalized_subtype"
)
plot
```

![Multiplet mutations plotted per Dose. Each circle represents a
mutation, coloured by mutation subtype. The size of the circle is scaled
by the mutation's alternative
depth.](General_Visualizations_files/figure-html/bubble-plot-1.png)

Multiplet mutations plotted per Dose. Each circle represents a mutation,
coloured by mutation subtype. The size of the circle is scaled by the
mutation’s alternative depth.

## Radar Plots

[`plot_radar()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/plot_radar.md)
is used to visualize mutation frequencies across specified groups as a
radar/spider plot. Plots may also be facetted across a second group.

*Example 2. Plot the mutation frequency for each of the 20 genomic
targets of the TwinStrand Mutagenesis Panel. Facet panels by dose
group.*

First we will calculate the average MFmin for each genomic target across
dose groups. We will also define the order in which the genomic targets
should appear on our plot. We will arrange genomic targets based on
their genic context so that we can visualize differences in mutation
frequency driven be genic context.

``` r
mf <- calculate_mf(
    mutation_data = example_data,
    cols_to_group = c("sample", "label"),
    retain_metadata_cols = c("dose_group", "genic_context")
)

label_order <- mf %>%
    dplyr::arrange(genic_context) %>%
    dplyr::pull(label) %>%
    unique()

avg <- mf %>%
    dplyr::group_by(dose_group, label) %>%
    dplyr::summarise(mean_mf = mean(mf_min))

avg$label <- factor(avg$label, levels = label_order)
avg$dose_group <- factor(avg$dose_group,
    levels = c("Control", "Low", "Medium", "High")
)
```

``` r
plot <- plot_radar(
    mf_data = avg,
    response_col = "mean_mf",
    label_col = "label",
    facet_col = "dose_group",
    indiv_y = FALSE
)
```

![Average Minimum Mutation Frequency (mutation/bp) per genomic target.
Plots are facetted by dose
group.](General_Visualizations_files/figure-html/radar-plot-1.png)

Average Minimum Mutation Frequency (mutation/bp) per genomic target.
Plots are facetted by dose group.

``` r
plot
```

    ## NULL

This plot shows us that average mutation frequency increases with dose.
Second, intergenic targets (targets on the left of the plot) show higher
mutation frequencies compared to genic targets (targets on the right
side of the plot).

## Lollipop Plots

[`plot_lollipop()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/plot_lollipop.md)is
used to visualize recurrent mutations. Recurrent mutations may be
analysed across the entire data sets or across user-specified
variable(s) (ex. sample, contig, dose). Mutations may also be analysed
at the various subtype resolutions. Mutations i a group that occur above
a specified recurrence threshold are plotted by genomic position. If
specified, mutations are colored by their mutation subtype. For each
level in the specified group, a lollipop plot will be generated and
stored in a list.

*Example 3. We will visualize mutations that occur a minimum of 2 times
in a group* *Mutations will be grouped by genomic target (label) and by
dose group* *Recurring mutations will be resolved at the base_6 level*

``` r
plot_list <- plot_lollipop(
    mutation_data = example_data,
    min_recurrence = 2,
    group_col = c("label", "dose_group"),
    subtype_resolution = "base_6"
)

plot_list[50]
```

    ## $`chr2 | High`

![Recurrent mutations in target chr 2 found in the High dose group
plotted by genomic position and coloured by normalized subtype. Plots
are separated by genomic
target.](General_Visualizations_files/figure-html/plot-lollipop-1.png)

Recurrent mutations in target chr 2 found in the High dose group plotted
by genomic position and coloured by normalized subtype. Plots are
separated by genomic target.

## References

## Appendix

### Session Info

    ## R Under development (unstable) (2026-02-04 r89376)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ExperimentHub_3.1.0 AnnotationHub_4.1.0 BiocFileCache_3.1.0
    ## [4] dbplyr_2.5.1        BiocGenerics_0.57.0 generics_0.1.4     
    ## [7] MutSeqR_0.99.9      htmltools_0.5.9     DT_0.34.0          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] DBI_1.2.3                   bitops_1.0-9               
    ##   [3] httr2_1.2.2                 rlang_1.1.7                
    ##   [5] magrittr_2.0.4              otel_0.2.0                 
    ##   [7] matrixStats_1.5.0           compiler_4.6.0             
    ##   [9] RSQLite_2.4.6               GenomicFeatures_1.63.1     
    ##  [11] png_0.1-8                   systemfonts_1.3.1          
    ##  [13] vctrs_0.7.1                 stringr_1.6.0              
    ##  [15] pkgconfig_2.0.3             crayon_1.5.3               
    ##  [17] fastmap_1.2.0               backports_1.5.0            
    ##  [19] XVector_0.51.0              labeling_0.4.3             
    ##  [21] Rsamtools_2.27.0            rmarkdown_2.30             
    ##  [23] ragg_1.5.0                  purrr_1.2.1                
    ##  [25] bit_4.6.0                   xfun_0.56                  
    ##  [27] cachem_1.1.0                cigarillo_1.1.0            
    ##  [29] jsonlite_2.0.0              blob_1.3.0                 
    ##  [31] DelayedArray_0.37.0         BiocParallel_1.45.0        
    ##  [33] parallel_4.6.0              R6_2.6.1                   
    ##  [35] plyranges_1.31.1            VariantAnnotation_1.57.1   
    ##  [37] bslib_0.10.0                stringi_1.8.7              
    ##  [39] RColorBrewer_1.1-3          rtracklayer_1.71.3         
    ##  [41] GenomicRanges_1.63.1        jquerylib_0.1.4            
    ##  [43] Rcpp_1.1.1                  Seqinfo_1.1.0              
    ##  [45] SummarizedExperiment_1.41.0 knitr_1.51                 
    ##  [47] IRanges_2.45.0              Matrix_1.7-4               
    ##  [49] tidyselect_1.2.1            dichromat_2.0-0.1          
    ##  [51] abind_1.4-8                 yaml_2.3.12                
    ##  [53] codetools_0.2-20            curl_7.0.0                 
    ##  [55] lattice_0.22-7              tibble_3.3.1               
    ##  [57] withr_3.0.2                 Biobase_2.71.0             
    ##  [59] KEGGREST_1.51.1             S7_0.2.1                   
    ##  [61] evaluate_1.0.5              desc_1.4.3                 
    ##  [63] Biostrings_2.79.4           pillar_1.11.1              
    ##  [65] BiocManager_1.30.27         filelock_1.0.3             
    ##  [67] MatrixGenerics_1.23.0       checkmate_2.3.4            
    ##  [69] stats4_4.6.0                rprojroot_2.1.1            
    ##  [71] RCurl_1.98-1.17             BiocVersion_3.23.1         
    ##  [73] S4Vectors_0.49.0            ggplot2_4.0.2              
    ##  [75] scales_1.4.0                glue_1.8.0                 
    ##  [77] tools_4.6.0                 BiocIO_1.21.0              
    ##  [79] data.table_1.18.2.1         BSgenome_1.79.1            
    ##  [81] GenomicAlignments_1.47.0    fmsb_0.7.6                 
    ##  [83] fs_1.6.6                    XML_3.99-0.20              
    ##  [85] grid_4.6.0                  tidyr_1.3.2                
    ##  [87] AnnotationDbi_1.73.0        restfulr_0.0.16            
    ##  [89] cli_3.6.5                   rappdirs_0.3.4             
    ##  [91] textshaping_1.0.4           S4Arrays_1.11.1            
    ##  [93] ggdendro_0.2.0              dplyr_1.2.0                
    ##  [95] gtable_0.3.6                sass_0.4.10                
    ##  [97] digest_0.6.39               SparseArray_1.11.10        
    ##  [99] rjson_0.2.23                htmlwidgets_1.6.4          
    ## [101] farver_2.1.2                memoise_2.0.1              
    ## [103] pkgdown_2.2.0               lifecycle_1.0.5            
    ## [105] httr_1.4.7                  here_1.0.2                 
    ## [107] packcircles_0.3.7           bit64_4.6.0-1              
    ## [109] MASS_7.3-65
