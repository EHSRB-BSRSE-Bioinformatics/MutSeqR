# MutSeqR: Analysing Genomic Targets

In this vignette, we will provide examples on how to perform an analysis
of mutation frequencies across a panel of genomic targets. We are using
TwinStrand’s [Mouse Mutagenesis
Panel](https://github.com/twinstrandbio/twinstrandbio-reference-data/)
as an example which consists of 20 2.4kb targets; one on each autosome,
with the exeption of chromosome 1 which has 2 targets (chr1 & chr1.2).

The example data is taken from @leblanc-2022. This data consists of 24
mouse bone marrow samples sequenced with Duplex Sequencing. Mice were
exposed to three doses of benzo\[a\]pyrene (BaP) alongside vehicle
controls, n = 6. The goal of this project was to quntify the mutagenic
effects of BaP in mouse bone marrow. Example data is retrieved from
MutSeqRData, an ExperimentHub data package available via Bioconductor.
Example data is retrieved from the ExperimentHub index (eh) through
specific accessors.

## Load MutSeqR and Example Data

``` r

library(MutSeqR)

library(ExperimentHub)
# load the index
eh <- ExperimentHub()
```

## Generalized Linear Modeling by Target

### Model by Target

We will measure the effect of genomic target and BaP dose on mutation
frequency (MF) using a generalized linear mixed model using
[`model_mf()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/model_mf.md).
We will determine whether the MF of each BaP dose group is significantly
different from the Control individually for all 20 genomic targets. In
this model, dose group and target label will be our fixed effects. We
include the interaction between the two fixed effects. Sample will be
used as a random effect to control for repeated measurements.

First, load the example data. This mutation data has already undergone
necessary filtering with `filter_mut` (demonstrated in the vignette;
MutSeqR Introduction). We will then calculate MF for each sample and
genomic locus, while retaining the dose column. The “label” column
refers to the label of the genomic target.

``` r

example_data <- eh[["EH9861"]]

mf_data_rg <- calculate_mf(
  mutation_data = example_data,
  cols_to_group = c("sample", "label"),
  subtype_resolution = "none",
  retain_metadata_cols = "dose_group"
)
```

Next, create a contrasts table that compares each dose group back to the
control for each genomic locus:

``` r

combinations <- expand.grid(dose_group = unique(mf_data_rg$dose_group),
                            label = unique(mf_data_rg$label))
combinations <- combinations[combinations$dose_group != "Control", ]
combinations$col1 <- with(combinations, paste(dose_group, label, sep = ":"))
combinations$col2 <- with(combinations, paste("Control", label, sep = ":"))
contrasts2 <- combinations[, c("col1", "col2")]
```

Finally, run the model. As this model is more complex, we will improve
covergence by supplying the control argument, which is passed directly
to [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html).

``` r

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
```

#### Residuals Histogram

![GLMM residuals of MFmin modelled as an effect of Dose and Genomic
Target. x is pearson's residuals, y is frequency. Plotted to validate
model assumptions. n =
24.](Target_Analysis_files/figure-html/model-rg-hist-1.png)

GLMM residuals of MFmin modelled as an effect of Dose and Genomic
Target. x is pearson’s residuals, y is frequency. Plotted to validate
model assumptions. n = 24.

#### Residuals QQ-plot

![GLMM residuals of MFmin modelled as an effect of Dose and Genomic
Target expressed as a quantile-quantile plot. Y is the pearson's
residuals of the model in ascending order x is the quantiles of standard
normal distribution for n of 24. Plotted to validate model
assumptions.](Target_Analysis_files/figure-html/model-rg-qq-1.png)

GLMM residuals of MFmin modelled as an effect of Dose and Genomic Target
expressed as a quantile-quantile plot. Y is the pearson’s residuals of
the model in ascending order x is the quantiles of standard normal
distribution for n of 24. Plotted to validate model assumptions.

#### Model Summary

``` r

model_by_target$summary
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: cbind(sum_min, group_depth) ~ dose_group * label + (1 | sample)
    ##    Data: mf_data
    ## Control: ..1
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##    2825.5    3163.6   -1331.7    2663.5       399 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7513 -0.7140 -0.0326  0.6138  5.4812 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  sample (Intercept) 0.01088  0.1043  
    ## Number of obs: 480, groups:  sample, 24
    ## 
    ## Fixed effects:
    ##                                Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                  -1.596e+01  2.225e-01 -71.723  < 2e-16 ***
    ## dose_groupHigh                1.903e+00  2.398e-01   7.938 2.05e-15 ***
    ## dose_groupLow                 6.208e-01  2.741e-01   2.265 0.023521 *  
    ## dose_groupMedium              1.293e+00  2.614e-01   4.948 7.51e-07 ***
    ## labelchr1.2                   1.382e-01  2.828e-01   0.489 0.625060    
    ## labelchr10                    4.517e-01  2.685e-01   1.682 0.092520 .  
    ## labelchr11                    8.366e-01  2.644e-01   3.164 0.001557 ** 
    ## labelchr12                    5.183e-01  2.664e-01   1.946 0.051662 .  
    ## labelchr13                    1.984e-02  2.911e-01   0.068 0.945662    
    ## labelchr14                    6.560e-01  2.674e-01   2.453 0.014166 *  
    ## labelchr15                    5.084e-01  2.664e-01   1.909 0.056323 .  
    ## labelchr16                    7.757e-01  2.696e-01   2.877 0.004016 ** 
    ## labelchr17                    4.593e-01  2.888e-01   1.590 0.111755    
    ## labelchr18                   -7.519e-03  2.962e-01  -0.025 0.979746    
    ## labelchr19                   -2.378e-01  3.088e-01  -0.770 0.441157    
    ## labelchr2                     6.258e-01  2.664e-01   2.349 0.018810 *  
    ## labelchr3                     2.568e-01  2.867e-01   0.896 0.370334    
    ## labelchr4                     5.106e-01  2.777e-01   1.839 0.065945 .  
    ## labelchr5                     4.209e-01  2.828e-01   1.488 0.136683    
    ## labelchr6                     2.937e-01  2.810e-01   1.045 0.295958    
    ## labelchr7                    -3.983e-04  2.962e-01  -0.001 0.998927    
    ## labelchr8                     4.383e-01  2.911e-01   1.506 0.132186    
    ## labelchr9                     7.000e-01  2.674e-01   2.618 0.008852 ** 
    ## dose_groupHigh:labelchr1.2   -2.483e-01  3.026e-01  -0.821 0.411880    
    ## dose_groupLow:labelchr1.2    -2.540e-03  3.463e-01  -0.007 0.994147    
    ## dose_groupMedium:labelchr1.2 -2.494e-02  3.293e-01  -0.076 0.939613    
    ## dose_groupHigh:labelchr10    -4.750e-01  2.885e-01  -1.646 0.099669 .  
    ## dose_groupLow:labelchr10     -3.356e-01  3.364e-01  -0.998 0.318378    
    ## dose_groupMedium:labelchr10  -3.318e-01  3.171e-01  -1.046 0.295381    
    ## dose_groupHigh:labelchr11    -5.554e-02  2.812e-01  -0.197 0.843447    
    ## dose_groupLow:labelchr11      4.252e-01  3.178e-01   1.338 0.180873    
    ## dose_groupMedium:labelchr11   4.359e-01  3.032e-01   1.438 0.150523    
    ## dose_groupHigh:labelchr12    -1.676e-01  2.839e-01  -0.591 0.554833    
    ## dose_groupLow:labelchr12     -1.544e-01  3.291e-01  -0.469 0.638845    
    ## dose_groupMedium:labelchr12  -1.756e-01  3.119e-01  -0.563 0.573379    
    ## dose_groupHigh:labelchr13    -3.978e-01  3.128e-01  -1.272 0.203483    
    ## dose_groupLow:labelchr13      2.170e-03  3.565e-01   0.006 0.995143    
    ## dose_groupMedium:labelchr13  -2.386e-01  3.429e-01  -0.696 0.486507    
    ## dose_groupHigh:labelchr14     1.207e-01  2.836e-01   0.426 0.670457    
    ## dose_groupLow:labelchr14      2.822e-01  3.231e-01   0.873 0.382465    
    ## dose_groupMedium:labelchr14   4.877e-01  3.058e-01   1.595 0.110778    
    ## dose_groupHigh:labelchr15    -1.145e+00  2.938e-01  -3.898 9.68e-05 ***
    ## dose_groupLow:labelchr15     -4.122e-01  3.351e-01  -1.230 0.218786    
    ## dose_groupMedium:labelchr15  -6.712e-01  3.218e-01  -2.086 0.036991 *  
    ## dose_groupHigh:labelchr16    -2.607e-01  2.881e-01  -0.905 0.365530    
    ## dose_groupLow:labelchr16      1.133e-01  3.283e-01   0.345 0.729968    
    ## dose_groupMedium:labelchr16   6.542e-02  3.128e-01   0.209 0.834334    
    ## dose_groupHigh:labelchr17     1.326e-01  3.062e-01   0.433 0.664948    
    ## dose_groupLow:labelchr17      3.705e-01  3.461e-01   1.071 0.284375    
    ## dose_groupMedium:labelchr17   3.523e-01  3.311e-01   1.064 0.287312    
    ## dose_groupHigh:labelchr18     3.822e-01  3.121e-01   1.225 0.220757    
    ## dose_groupLow:labelchr18      4.491e-01  3.529e-01   1.273 0.203175    
    ## dose_groupMedium:labelchr18   5.986e-01  3.350e-01   1.787 0.073985 .  
    ## dose_groupHigh:labelchr19     6.478e-02  3.277e-01   0.198 0.843289    
    ## dose_groupLow:labelchr19      3.157e-01  3.693e-01   0.855 0.392672    
    ## dose_groupMedium:labelchr19   3.342e-01  3.526e-01   0.948 0.343213    
    ## dose_groupHigh:labelchr2     -6.041e-01  2.872e-01  -2.104 0.035416 *  
    ## dose_groupLow:labelchr2      -2.459e-01  3.316e-01  -0.742 0.458322    
    ## dose_groupMedium:labelchr2   -2.993e-01  3.141e-01  -0.953 0.340667    
    ## dose_groupHigh:labelchr3     -1.072e+00  3.181e-01  -3.371 0.000749 ***
    ## dose_groupLow:labelchr3      -5.817e-01  3.705e-01  -1.570 0.116395    
    ## dose_groupMedium:labelchr3   -6.911e-01  3.507e-01  -1.971 0.048744 *  
    ## dose_groupHigh:labelchr4     -2.788e-01  2.971e-01  -0.939 0.347874    
    ## dose_groupLow:labelchr4      -1.492e-01  3.438e-01  -0.434 0.664228    
    ## dose_groupMedium:labelchr4   -5.001e-02  3.238e-01  -0.154 0.877254    
    ## dose_groupHigh:labelchr5     -2.603e-01  3.026e-01  -0.860 0.389725    
    ## dose_groupLow:labelchr5      -4.738e-02  3.475e-01  -0.136 0.891550    
    ## dose_groupMedium:labelchr5   -1.236e-01  3.315e-01  -0.373 0.709370    
    ## dose_groupHigh:labelchr6     -2.353e-01  3.002e-01  -0.784 0.433083    
    ## dose_groupLow:labelchr6      -4.394e-02  3.452e-01  -0.127 0.898716    
    ## dose_groupMedium:labelchr6   -2.395e-01  3.311e-01  -0.724 0.469372    
    ## dose_groupHigh:labelchr7      1.376e-01  3.135e-01   0.439 0.660753    
    ## dose_groupLow:labelchr7       4.021e-01  3.534e-01   1.138 0.255143    
    ## dose_groupMedium:labelchr7    5.092e-01  3.359e-01   1.516 0.129484    
    ## dose_groupHigh:labelchr8      2.822e-01  3.078e-01   0.917 0.359175    
    ## dose_groupLow:labelchr8       5.457e-01  3.458e-01   1.578 0.114497    
    ## dose_groupMedium:labelchr8    6.388e-01  3.302e-01   1.934 0.053087 .  
    ## dose_groupHigh:labelchr9     -4.547e-01  2.870e-01  -1.584 0.113151    
    ## dose_groupLow:labelchr9      -1.701e-01  3.308e-01  -0.514 0.607104    
    ## dose_groupMedium:labelchr9   -4.265e-01  3.175e-01  -1.343 0.179168    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### ANOVA

``` r

model_by_target$anova
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: cbind(sum_min, group_depth)
    ##                    Chisq Df Pr(>Chisq)    
    ## dose_group        602.15  3  < 2.2e-16 ***
    ## label            1238.55 19  < 2.2e-16 ***
    ## dose_group:label  129.18 57  1.597e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### Model Estimates

**Table 20.** model_by_target\$point_estimates: Estimated Mean Mutation
Frequency per Dose and Genomic Target.

#### Pairwise Comparisons

**Table 21.** model_by_target\$pairwise_comparisons.

#### Plot Model by Dose and Target

``` r

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
```

![Mean Mutation Frequency Minimum (mutations/bp) per Genomic Target and
Dose estimated using a generalized linear mixed model. Error bars are
the SEM. Symbols indicate significance differences (p \< 0.05) between
dose levels for individual genomic
regions.](Target_Analysis_files/figure-html/plot-model-rg-1.png)

Mean Mutation Frequency Minimum (mutations/bp) per Genomic Target and
Dose estimated using a generalized linear mixed model. Error bars are
the SEM. Symbols indicate significance differences (p \< 0.05) between
dose levels for individual genomic regions.

## Visualizae Mutations along Target Regions

[`plot_lollipop()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/plot_lollipop.md)

## Retrieve Sequences of Target regions

[`get_seq()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/get_seq.md)
will retrive raw nucleotide sequences for specified genomic intervals.
The function can retrieve sequences in one of two ways. Either specify
an installed BSgenome from which to retrieve sequences or retrieve
sequences through API request to the UCSC database based on supplied
species and genome assembly.

Supply regions with a filepath, data frame or GRanges object containing
the specified genomic intervals. TwinStrand’s Mutagenesis Panels are
stored in package files and can easily be retrieved.

Supply the function with the appropriate BS_genome from which to
retrieve sequences. You can use the find_BS_genome() function to
identify a BS genome to use based on species and genome assembly. Be
sure to install the BS genome prior to running get_seq() with
BiocManager().

To retrieve sequences from the UCSC data base, please input the species
and assembly genome to the species and genome parameters.

\**BS_genome, species, or genome* *do not need to be specified for the
Mutagenesis Panels as this information* *is stored*.

Sequences are returned within a *GRanges* object.

*Example 2.1. Retrieve the sequences for our example’s target panel,*
*TwinStrand’s Mouse Mutagenesis Panel*

``` r

regions_seq <- get_seq(regions = "TSpanel_mouse")
regions_seq
```

    ## GRanges object with 20 ranges and 6 metadata columns:
    ##        seqnames              ranges strand | target_size       label
    ##           <Rle>           <IRanges>  <Rle> |   <integer> <character>
    ##    [1]     chr1   69304218-69306617      * |        2400        chr1
    ##    [2]     chr1 155235939-155238338      * |        2400      chr1.2
    ##    [3]     chr2   50833176-50835575      * |        2400        chr2
    ##    [4]     chr3 109633161-109635560      * |        2400        chr3
    ##    [5]     chr4   96825281-96827680      * |        2400        chr4
    ##    ...      ...                 ...    ... .         ...         ...
    ##   [16]    chr15   66779763-66782162      * |        2400       chr15
    ##   [17]    chr16   72381581-72383980      * |        2400       chr16
    ##   [18]    chr17   94009029-94011428      * |        2400       chr17
    ##   [19]    chr18   81262079-81264478      * |        2400       chr18
    ##   [20]    chr19     4618814-4621213      * |        2400       chr19
    ##        genic_context region_GC_content      genome                sequence
    ##          <character>         <numeric> <character>          <DNAStringSet>
    ##    [1]    intergenic              37.3        mm10 CAATCTTTCT...CAAAATGCAA
    ##    [2]         genic              54.0        mm10 AATCTCCAGG...CAAGCACTGG
    ##    [3]    intergenic              45.3        mm10 TGTGCCCCAT...TGCTTGCCAC
    ##    [4]         genic              39.2        mm10 AACGATGAAT...GCACTCAAGA
    ##    [5]    intergenic              39.4        mm10 ATTGTTTGAA...CTCAGGGCCT
    ##    ...           ...               ...         ...                     ...
    ##   [16]         genic              44.0        mm10 GTGTCATTTT...CAGGTAGAGG
    ##   [17]    intergenic              38.3        mm10 TCTGTAGCAA...ATAAAAACTC
    ##   [18]    intergenic              35.2        mm10 TAAGGAAACT...TATCAAAGAT
    ##   [19]    intergenic              47.3        mm10 AGCCATCTCC...GGGACTCAGA
    ##   [20]         genic              56.1        mm10 TCCAGGCTGT...GCCCATGGAG
    ##   -------
    ##   seqinfo: 19 sequences from an unspecified genome; no seqlengths

*Example 2.2. Retrieve sequences for a custom interval of regions. We
will use* *the Human Mutagenesis Panel as an example.*

``` r

# We will load the TSpanel_human regions file to obtain
# an example list of regions.
human <- load_regions_file("TSpanel_human")
regions_seq <- get_seq(regions = human,
                       is_0_based_rg = FALSE,
                       BS_genome = find_BS_genome("human", "hg38"),
                       padding = 0)
regions_seq
```

    ## GRanges object with 20 ranges and 7 metadata columns:
    ##        seqnames              ranges strand | target_size description
    ##           <Rle>           <IRanges>  <Rle> |   <integer> <character>
    ##    [1]    chr11 108510788-108513187      * |        2400 region_1111
    ##    [2]    chr13   75803913-75806312      * |        2400 region_1501
    ##    [3]    chr14   74661756-74664155      * |        2400 region_1725
    ##    [4]    chr18     5749265-5751664      * |        2400 region_2457
    ##    [5]     chr2   40162768-40165167      * |        2400 region_2896
    ##    ...      ...                 ...    ... .         ...         ...
    ##   [16]    chr15   46089738-46092137      * |        2400 region_1904
    ##   [17]    chr17   70672727-70675126      * |        2400 region_2378
    ##   [18]    chr21   23665977-23668376      * |        2400 region_3515
    ##   [19]    chr22   48262371-48264770      * |        2400 region_3703
    ##   [20]    chr10 128969038-128971437      * |        2400  region_784
    ##        genic_context        gene      genome       label
    ##          <character> <character> <character> <character>
    ##    [1]         genic       EXPH5        hg38       chr11
    ##    [2]         genic        LMO7        hg38       chr13
    ##    [3]         genic       AREL1        hg38       chr14
    ##    [4]         genic   MIR3976HG        hg38       chr18
    ##    [5]         genic  SLC8A1-AS1        hg38        chr2
    ##    ...           ...         ...         ...         ...
    ##   [16]    intergenic        <NA>        hg38       chr15
    ##   [17]    intergenic        <NA>        hg38       chr17
    ##   [18]    intergenic        <NA>        hg38       chr21
    ##   [19]    intergenic        <NA>        hg38       chr22
    ##   [20]    intergenic        <NA>        hg38       chr10
    ##                       sequence
    ##                 <DNAStringSet>
    ##    [1] GTTTCCTTCA...CTTTCCTGGA
    ##    [2] AGAATTATTT...TCAGACAACC
    ##    [3] TTCCCTGGTT...AAGATACACT
    ##    [4] TGGCAACTTG...AATGAAAACA
    ##    [5] CTAGATTTTC...AGCATATCAC
    ##    ...                     ...
    ##   [16] TGAACAGACA...ATAAATTGCT
    ##   [17] GTGGTGATCA...TAAAGATTCT
    ##   [18] TGTAATAATG...TCCAGTCATT
    ##   [19] GTAAAGGCAG...TCCACAGCAG
    ##   [20] GTGGACTGAT...TTCTCACACT
    ##   -------
    ##   seqinfo: 20 sequences from an unspecified genome; no seqlengths

*Example 2.3. Retrieve sequences using UCSC database.*

``` r

regions_seq <- get_seq(regions = human,
                       is_0_based_rg = FALSE,
                       padding = 0,
                       species = "human",
                       genome = "hg38",
                       ucsc = TRUE)
```

Sequences can be exported as FASTA files with
[`write_reference_fasta()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/write_reference_fasta.md).
Supply this function with the GRanges object containing the sequences of
the regions. Each one will be written to a single FASTA file. The FASTA
file will be saved to the `output_path`. If NULL, the file will be saved
to the working directory.

``` r

write_reference_fasta(regions_seq, output_path = NULL)
```

## References

## Appendix

### Session Info

    ## R Under development (unstable) (2026-10-08 r90650)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.5  GenomeInfoDb_1.48.0               
    ##  [3] BSgenome.Mmusculus.UCSC.mm10_1.4.3 BSgenome_1.80.0                   
    ##  [5] rtracklayer_1.72.0                 BiocIO_1.22.0                     
    ##  [7] Biostrings_2.80.2                  XVector_0.52.0                    
    ##  [9] GenomicRanges_1.64.0               Seqinfo_1.2.0                     
    ## [11] IRanges_2.46.0                     S4Vectors_0.50.3                  
    ## [13] ExperimentHub_3.2.2                AnnotationHub_4.2.2               
    ## [15] BiocFileCache_3.2.0                dbplyr_2.6.0                      
    ## [17] BiocGenerics_0.58.1                generics_0.1.4                    
    ## [19] MutSeqR_1.1.1                      htmltools_0.5.9                   
    ## [21] DT_0.34.0                          BiocStyle_2.40.0                  
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3          ggdendro_0.2.0             
    ##   [3] jsonlite_2.0.0              magrittr_2.0.5             
    ##   [5] GenomicFeatures_1.64.0      nloptr_2.2.1               
    ##   [7] farver_2.1.2                rmarkdown_2.32             
    ##   [9] fs_2.1.0                    ragg_1.5.2                 
    ##  [11] vctrs_0.7.3                 minqa_1.2.8                
    ##  [13] memoise_2.0.1               Rsamtools_2.28.0           
    ##  [15] RCurl_1.98-1.20             S4Arrays_1.12.1            
    ##  [17] curl_8.0.0                  broom_1.0.13               
    ##  [19] SparseArray_1.12.3          Formula_1.2-6              
    ##  [21] sass_0.4.10                 bslib_0.12.0               
    ##  [23] htmlwidgets_1.6.4           desc_1.4.3                 
    ##  [25] httr2_1.3.0                 zoo_1.9-1                  
    ##  [27] cachem_1.1.0                GenomicAlignments_1.48.0   
    ##  [29] lifecycle_1.0.5             pkgconfig_2.0.3            
    ##  [31] Matrix_1.7-6                R6_2.6.1                   
    ##  [33] fastmap_1.2.0               rbibutils_2.4.1            
    ##  [35] MatrixGenerics_1.24.0       digest_0.6.39              
    ##  [37] colorspace_2.1-4            AnnotationDbi_1.74.0       
    ##  [39] rprojroot_2.1.1             crosstalk_1.2.2            
    ##  [41] textshaping_1.0.5           RSQLite_3.53.3             
    ##  [43] labeling_0.4.3              filelock_1.0.3             
    ##  [45] httr_1.4.9                  abind_1.4-8                
    ##  [47] compiler_4.7.0              here_1.0.2                 
    ##  [49] bit64_4.8.6                 withr_3.0.3                
    ##  [51] S7_0.2.2                    backports_1.5.1            
    ##  [53] BiocParallel_1.46.0         carData_3.0-6              
    ##  [55] DBI_1.3.0                   MASS_7.3-66                
    ##  [57] rappdirs_0.3.4              DelayedArray_0.38.2        
    ##  [59] rjson_0.2.23                tools_4.7.0                
    ##  [61] otel_0.2.0                  glue_1.8.1                 
    ##  [63] restfulr_0.0.17             nlme_3.1-171               
    ##  [65] grid_4.7.0                  gtable_0.3.6               
    ##  [67] tidyr_1.3.2                 data.table_1.18.6.1        
    ##  [69] doBy_4.7.2                  xml2_1.6.0                 
    ##  [71] car_3.1-5                   Deriv_4.3.5                
    ##  [73] BiocVersion_3.23.1          pillar_1.11.1              
    ##  [75] stringr_1.6.0               splines_4.7.0              
    ##  [77] dplyr_1.2.1                 lattice_0.23-1             
    ##  [79] bit_4.6.0                   tidyselect_1.2.1           
    ##  [81] knitr_1.52                  reformulas_0.4.4           
    ##  [83] bookdown_0.48               urca_1.3-4                 
    ##  [85] SummarizedExperiment_1.42.0 forecast_9.0.2             
    ##  [87] xfun_0.61                   Biobase_2.72.0             
    ##  [89] timeDate_4052.112           matrixStats_1.5.0          
    ##  [91] UCSC.utils_1.8.0            stringi_1.8.9              
    ##  [93] yaml_2.3.12                 boot_1.3-32                
    ##  [95] evaluate_1.0.5              codetools_0.2-20           
    ##  [97] cigarillo_1.2.1             tibble_3.3.1               
    ##  [99] BiocManager_1.30.27         cli_3.6.6                  
    ## [101] Rdpack_2.6.6                systemfonts_1.3.2          
    ## [103] jquerylib_0.1.4             dichromat_2.0-1            
    ## [105] modelr_0.1.11               Rcpp_1.1.2                 
    ## [107] png_0.1-9                   XML_3.99-0.25              
    ## [109] parallel_4.7.0              pkgdown_2.2.1              
    ## [111] fracdiff_1.5-4              ggplot2_4.0.3              
    ## [113] blob_1.3.0                  plyranges_1.32.0           
    ## [115] bitops_1.1-0                lme4_2.0-6                 
    ## [117] VariantAnnotation_1.58.0    scales_1.4.0               
    ## [119] purrr_1.2.2                 crayon_1.5.3               
    ## [121] rlang_1.3.0                 cowplot_1.2.0              
    ## [123] KEGGREST_1.52.2
